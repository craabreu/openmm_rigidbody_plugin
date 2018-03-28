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

#include "ReferenceRigidBodyKernels.h"
#include "RigidBodyForce.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferenceConstraints.h"
#include "openmm/reference/ReferencePlatform.h"

using namespace RigidBodyPlugin;
using namespace OpenMM;
using namespace std;

static vector<RealVec>& extractPositions(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->positions);
}

static vector<RealVec>& extractVelocities(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->velocities);
}

static vector<RealVec>& extractForces(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *((vector<RealVec>*) data->forces);
}

static ReferenceConstraints& extractConstraints(ContextImpl& context) {
    ReferencePlatform::PlatformData* data = reinterpret_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    return *(ReferenceConstraints*) data->constraints;
}

/**
 * Compute the kinetic energy of the system, possibly shifting the velocities in time to account
 * for a leapfrog integrator.
 */
static double computeShiftedKineticEnergy(ContextImpl& context, vector<double>& masses, double timeShift) {
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3>& forceData = extractForces(context);
    int numParticles = context.getSystem().getNumParticles();
    
    // Compute the shifted velocities.
    
    vector<Vec3> shiftedVel(numParticles);
    for (int i = 0; i < numParticles; ++i) {
        if (masses[i] > 0)
            shiftedVel[i] = velData[i]+forceData[i]*(timeShift/masses[i]);
        else
            shiftedVel[i] = velData[i];
    }
    
    // Apply constraints to them.
    
    vector<double> inverseMasses(numParticles);
    for (int i = 0; i < numParticles; i++)
        inverseMasses[i] = (masses[i] == 0 ? 0 : 1/masses[i]);
    extractConstraints(context).applyToVelocities(posData, shiftedVel, inverseMasses, 1e-4);
    
    // Compute the kinetic energy.
    
    double energy = 0.0;
    for (int i = 0; i < numParticles; ++i)
        if (masses[i] > 0)
            energy += masses[i]*(shiftedVel[i].dot(shiftedVel[i]));
    return 0.5*energy;
}

void ReferenceCalcRigidBodyForceKernel::initialize(const System& system, const RigidBodyForce& force) {
    // Initialize bond parameters.
    
    int numBonds = force.getNumBonds();
    particle1.resize(numBonds);
    particle2.resize(numBonds);
    length.resize(numBonds);
    k.resize(numBonds);
    for (int i = 0; i < numBonds; i++)
        force.getBondParameters(i, particle1[i], particle2[i], length[i], k[i]);
}

double ReferenceCalcRigidBodyForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    vector<RealVec>& pos = extractPositions(context);
    vector<RealVec>& force = extractForces(context);
    int numBonds = particle1.size();
    double energy = 0;
    
    // Compute the interactions.
    
    for (int i = 0; i < numBonds; i++) {
        int p1 = particle1[i];
        int p2 = particle2[i];
        RealVec delta = pos[p1]-pos[p2];
        RealOpenMM r2 = delta.dot(delta);
        RealOpenMM r = sqrt(r2);
        RealOpenMM dr = (r-length[i]);
        RealOpenMM dr2 = dr*dr;
        energy += k[i]*dr2*dr2;
        RealOpenMM dEdR = 4*k[i]*dr2*dr;
        dEdR = (r > 0) ? (dEdR/r) : 0;
        force[p1] -= delta*dEdR;
        force[p2] += delta*dEdR;
    }
    return energy;
}

void ReferenceCalcRigidBodyForceKernel::copyParametersToContext(ContextImpl& context, const RigidBodyForce& force) {
    if (force.getNumBonds() != particle1.size())
        throw OpenMMException("updateParametersInContext: The number of RigidBody bonds has changed");
    for (int i = 0; i < force.getNumBonds(); i++) {
        int p1, p2;
        force.getBondParameters(i, p1, p2, length[i], k[i]);
        if (p1 != particle1[i] || p2 != particle2[i])
            throw OpenMMException("updateParametersInContext: A particle index has changed");
    }
}

ReferenceIntegrateRigidBodyStepKernel::~ReferenceIntegrateRigidBodyStepKernel() {
    if (dynamics)
        delete dynamics;
}

void ReferenceIntegrateRigidBodyStepKernel::initialize(const System& system, const RigidBodyIntegrator& integrator) {
    int numParticles = system.getNumParticles();
    masses.resize(numParticles);
    for (int i = 0; i < numParticles; ++i)
        masses[i] = system.getParticleMass(i);
}

void ReferenceIntegrateRigidBodyStepKernel::execute(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    double stepSize = integrator.getStepSize();
    vector<Vec3>& posData = extractPositions(context);
    vector<Vec3>& velData = extractVelocities(context);
    vector<Vec3>& forceData = extractForces(context);
    if (dynamics == 0 || stepSize != prevStepSize) {
        // Recreate the computation objects with the new parameters.
        
        if (dynamics)
            delete dynamics;
        dynamics = new ReferenceRigidBodyDynamics(context.getSystem().getNumParticles(), stepSize);
        dynamics->setReferenceConstraintAlgorithm(&extractConstraints(context));
        prevStepSize = stepSize;
    }
    dynamics->update(context.getSystem(), posData, velData, forceData, masses, integrator.getConstraintTolerance());
    data.time += stepSize;
    data.stepCount++;
}

double ReferenceIntegrateRigidBodyStepKernel::computeKineticEnergy(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    return computeShiftedKineticEnergy(context, masses, 0.5*integrator.getStepSize());
}
