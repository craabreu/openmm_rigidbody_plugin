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
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/reference/RealVec.h"
#include "openmm/reference/ReferenceConstraints.h"
#include "openmm/reference/ReferencePlatform.h"
#include "openmm/reference/ReferenceVirtualSites.h"

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

ReferenceIntegrateRigidBodyStepKernel::~ReferenceIntegrateRigidBodyStepKernel() {
}

void ReferenceIntegrateRigidBodyStepKernel::initialize(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    const System& system = context.getSystem();
    int numAtoms = system.getNumParticles();
    hasMass.resize(numAtoms);
    invMass.resize(numAtoms);
    for (int i = 0; i < numAtoms; ++i) {
        double mass = system.getParticleMass(i);
        hasMass[i] = mass != 0.0;
        invMass[i] = hasMass[i] ? 1.0/mass : 0.0;
    }
    posDelta.resize(numAtoms);
}

void ReferenceIntegrateRigidBodyStepKernel::uploadBodySystem(RigidBodySystem& bodySystem) {
    this->bodySystem.copy(bodySystem);
}

void ReferenceIntegrateRigidBodyStepKernel::execute(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    double dt = integrator.getStepSize();
    double halfDt = 0.5*dt;
    double tol = integrator.getConstraintTolerance();
    vector<Vec3> savedPos;

    vector<Vec3>& R = extractPositions(context);
    vector<Vec3>& V = extractVelocities(context);
    vector<Vec3>& F = extractForces(context);

    int numFree = bodySystem.getNumFree();
    if (numFree != 0) {
        savedPos.resize(numFree);
        for (int k = 0; k < numFree; k++) {
            int i = bodySystem.getAtomIndex(k);
            if (hasMass[i]) {
                V[i] += F[i]*invMass[i]*halfDt;
                posDelta[i] = V[i]*dt;
                savedPos[k] = R[i] + posDelta[i];
            }
        }
        context.applyConstraints(tol);
        for (int k = 0; k < numFree; k++) {
            int i = bodySystem.getAtomIndex(k);
            if (hasMass[i])
                R[i] += posDelta[i];
        }
    }

    bodySystem.integratePart1(dt, R);
    ReferenceVirtualSites::computePositions(context.getSystem(), R);
    context.calcForcesAndEnergy(true, false);
    bodySystem.integratePart2(dt, F, V);

    if (numFree != 0) {
        double invDt = 1.0/dt;
        for (int k = 0; k < numFree; k++) {
            int i = bodySystem.getAtomIndex(k);
            if (hasMass[i])
                V[i] += F[i]*invMass[i]*halfDt + (R[i] - savedPos[k])*invDt;
        }
        context.applyVelocityConstraints(tol);
    }

    data.time += dt;
    data.stepCount++;
}

double ReferenceIntegrateRigidBodyStepKernel::computeKineticEnergy(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    vector<Vec3>& V = extractVelocities(context);
    double KE = 0.0;
    for (int k = 0; k < bodySystem.getNumFree(); k++) {
        int i = bodySystem.getAtomIndex(k);
        KE += V[i].dot(V[i])/invMass[i];
    }
    return 0.5*KE + bodySystem.getKineticEnergy();
}
