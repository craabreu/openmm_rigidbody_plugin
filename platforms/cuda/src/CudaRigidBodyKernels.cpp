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

#include "CudaRigidBodyKernels.h"
#include "CudaRigidBodyKernelSources.h"
#include "RigidBodyIntegrator.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/cuda/CudaPlatform.h"
#include "openmm/cuda/CudaBondedUtilities.h"
#include "openmm/cuda/CudaForceInfo.h"
#include "openmm/cuda/CudaIntegrationUtilities.h"

using namespace RigidBodyPlugin;
using namespace OpenMM;
using namespace std;

int CudaIntegrateRigidBodyStepKernel::getBodyDataSize(CUmodule& module) {
    CUfunction kernel = cu.getKernel(module, "getBodyDataSize");
    CUdeviceptr pointer;
    CUresult result;
    result = cuMemAlloc(&pointer, sizeof(int));
    if (result != CUDA_SUCCESS)
        throw OpenMMException("Error creating variable for rigid body data size");
    void* arg[] = {&pointer};
    cu.executeKernel(kernel, arg, 1, 128);
    int bodyDataSize;
    result = cuMemcpyDtoH(&bodyDataSize, pointer, sizeof(int));
    if (result != CUDA_SUCCESS)
        throw OpenMMException("Error retrieving rigid body data size from device memory");
    return bodyDataSize;
}

void CudaIntegrateRigidBodyStepKernel::initialize(const System& system, const RigidBodyIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(system);
    cu.setAsCurrent();
    map<string, string> defines;
    CUmodule module = cu.createModule(CudaRigidBodyKernelSources::rigidbodyintegrator, defines, "");
    kernel1 = cu.getKernel(module, "integrateRigidBodyPart1");
    kernel2 = cu.getKernel(module, "integrateRigidBodyPart2");

//    numAtoms = integrator.getNumActualAtoms();
//    numBodies = integrator.getNumBodies();
//    numFree = integrator.getNumFreeAtoms();

    // Allocate array of atom indices:
    int TileSize = cu.TileSize;
    paddedNumAtoms = TileSize*((numAtoms + TileSize - 1)/TileSize);
    atomIndex.initialize<int>(cu, paddedNumAtoms, "atomIndex");

    // Allocate array of body data:
    if (numBodies) {
        int bodyDataSize = getBodyDataSize(module);
        paddedNumBodies = cu.TileSize*((numBodies + cu.TileSize - 1)/cu.TileSize);
        bodyData.initialize(cu, paddedNumBodies, bodyDataSize, "bodyData");
    }
}

void CudaIntegrateRigidBodyStepKernel::execute(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    double dt = integrator.getStepSize();
    cu.getIntegrationUtilities().setNextStepSize(dt);

    // Call the first integration kernel.

    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    void* args1[] = {&numAtoms, &paddedNumAtoms, &cu.getIntegrationUtilities().getStepSize().getDevicePointer(), &cu.getPosq().getDevicePointer(), &posCorrection,
            &cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel1, args1, numAtoms, 128);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    void* args2[] = {&numAtoms, &cu.getIntegrationUtilities().getStepSize().getDevicePointer(), &cu.getPosq().getDevicePointer(), &posCorrection,
            &cu.getVelm().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numAtoms, 128);
    integration.computeVirtualSites();

    // Update the time and step count.

    cu.setTime(cu.getTime()+dt);
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
}

double CudaIntegrateRigidBodyStepKernel::computeKineticEnergy(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    return cu.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}
