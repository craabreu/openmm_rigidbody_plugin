/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "OpenCLRigidBodyKernels.h"
#include "OpenCLRigidBodyKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/opencl/OpenCLBondedUtilities.h"
#include "openmm/opencl/OpenCLForceInfo.h"
#include "openmm/opencl/OpenCLIntegrationUtilities.h"

using namespace RigidBodyPlugin;
using namespace OpenMM;
using namespace std;

class OpenCLRigidBodyForceInfo : public OpenCLForceInfo {
public:
    OpenCLRigidBodyForceInfo(const RigidBodyForce& force) : OpenCLForceInfo(0), force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumBonds();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2;
        double length, k;
        force.getBondParameters(index, particle1, particle2, length, k);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2;
        double length1, length2, k1, k2;
        force.getBondParameters(group1, particle1, particle2, length1, k1);
        force.getBondParameters(group2, particle1, particle2, length2, k2);
        return (length1 == length2 && k1 == k2);
    }
private:
    const RigidBodyForce& force;
};

static void setPosqCorrectionArg(OpenCLContext& cl, cl::Kernel& kernel, int index) {
    if (cl.getUseMixedPrecision())
        kernel.setArg<cl::Buffer>(index, cl.getPosqCorrection().getDeviceBuffer());
    else
        kernel.setArg<void*>(index, NULL);
}

OpenCLCalcRigidBodyForceKernel::~OpenCLCalcRigidBodyForceKernel() {
    if (params != NULL)
        delete params;
}

void OpenCLCalcRigidBodyForceKernel::initialize(const System& system, const RigidBodyForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    if (numBonds == 0)
        return;
    vector<vector<int> > atoms(numBonds, vector<int>(2));
    params = OpenCLArray::create<mm_float2>(cl, numBonds, "bondParams");
    vector<mm_float2> paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        double length, k;
        force.getBondParameters(startIndex+i, atoms[i][0], atoms[i][1], length, k);
        paramVector[i] = mm_float2((cl_float) length, (cl_float) k);
    }
    params->upload(paramVector);
    map<string, string> replacements;
    replacements["PARAMS"] = cl.getBondedUtilities().addArgument(params->getDeviceBuffer(), "float2");
    cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(OpenCLRigidBodyKernelSources::rigidbodyForce, replacements), force.getForceGroup());
    cl.addForce(new OpenCLRigidBodyForceInfo(force));
}

double OpenCLCalcRigidBodyForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    return 0.0;
}

void OpenCLCalcRigidBodyForceKernel::copyParametersToContext(ContextImpl& context, const RigidBodyForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumBonds()/numContexts;
    if (numBonds != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0)
        return;
    
    // Record the per-bond parameters.
    
    vector<mm_float2> paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        int atom1, atom2;
        double length, k;
        force.getBondParameters(startIndex+i, atom1, atom2, length, k);
        paramVector[i] = mm_float2((cl_float) length, (cl_float) k);
    }
    params->upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cl.invalidateMolecules();
}

OpenCLIntegrateRigidBodyStepKernel::~OpenCLIntegrateRigidBodyStepKernel() {
}

void OpenCLIntegrateRigidBodyStepKernel::initialize(const System& system, const RigidBodyIntegrator& integrator) {
    cl.getPlatformData().initializeContexts(system);
    cl::Program program = cl.createProgram(OpenCLRigidBodyKernelSources::rigidbodyintegrator, "");
    kernel1 = cl::Kernel(program, "integrateRigidBodyPart1");
    kernel2 = cl::Kernel(program, "integrateRigidBodyPart2");
}

void OpenCLIntegrateRigidBodyStepKernel::execute(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    OpenCLIntegrationUtilities& integration = cl.getIntegrationUtilities();
    int numAtoms = cl.getNumAtoms();
    double dt = integrator.getStepSize();
    if (!hasInitializedKernels) {
        hasInitializedKernels = true;
        kernel1.setArg<cl_int>(0, numAtoms);
        kernel1.setArg<cl::Buffer>(1, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
        setPosqCorrectionArg(cl, kernel1, 3);
        kernel1.setArg<cl::Buffer>(4, cl.getVelm().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(5, cl.getForce().getDeviceBuffer());
        kernel1.setArg<cl::Buffer>(6, integration.getPosDelta().getDeviceBuffer());
        kernel2.setArg<cl_int>(0, numAtoms);
        kernel2.setArg<cl::Buffer>(1, cl.getIntegrationUtilities().getStepSize().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(2, cl.getPosq().getDeviceBuffer());
        setPosqCorrectionArg(cl, kernel2, 3);
        kernel2.setArg<cl::Buffer>(4, cl.getVelm().getDeviceBuffer());
        kernel2.setArg<cl::Buffer>(5, integration.getPosDelta().getDeviceBuffer());
    }
    cl.getIntegrationUtilities().setNextStepSize(dt);

    // Call the first integration kernel.

    cl.executeKernel(kernel1, numAtoms);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    cl.executeKernel(kernel2, numAtoms);
    integration.computeVirtualSites();

    // Update the time and step count.

    cl.setTime(cl.getTime()+dt);
    cl.setStepCount(cl.getStepCount()+1);
    cl.reorderAtoms();
    
    // Reduce UI lag.
    
#ifdef WIN32
    cl.getQueue().flush();
#endif
}

double OpenCLIntegrateRigidBodyStepKernel::computeKineticEnergy(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    return cl.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}
