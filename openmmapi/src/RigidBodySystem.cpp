/* ---------------------------------------------------------------------------------------------- *
 *                                    OpenMM Rigid Body Plugin                                    *
 * ---------------------------------------------------------------------------------------------- */

#include "RigidBodySystem.h"
#include "internal/MatVec.h"
#include "openmm/Vec3.h"
#include "openmm/System.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "internal/eigenDecomposition.h"
#include <vector>
#include <cmath>
#include <algorithm>

#include <iostream> // TEMPORARY

using namespace RigidBodyPlugin;
using namespace OpenMM;
using std::vector;

/*--------------------------------------------------------------------------------------------------
  Create a clean list of rigid body indices. The rule is: if bodyIndices[i] <= 0 or is unique, then
  index[i] = 0. Otherwise, 1 <= index[i] <= number of distinct positive entries which are
  non-unique.
--------------------------------------------------------------------------------------------------*/

vector<int> cleanBodyIndices(const vector<int>& bodyIndices) {
    int listSize = bodyIndices.size();
    int maxIndex = *std::max_element(std::begin(bodyIndices), std::end(bodyIndices));
    vector<int> head(maxIndex, -1), next(listSize, -1), index(listSize, 0);
    for (int i = 0; i < listSize; i++)
        if (bodyIndices[i] > 0) {
            next[i] = head[bodyIndices[i]-1];
            head[bodyIndices[i]-1] = i;
        }
    int body = 0;
    for (int i = 0; i < maxIndex; i++) {
        int j = head[i];
        if (j != -1) {
            body++;
            while (j != -1) {
                index[j] = body;
                j = next[j];
            }
        }
    }
    return index;
}

/*--------------------------------------------------------------------------------------------------
  Create a data structure for the system of rigid bodies and free atoms.
--------------------------------------------------------------------------------------------------*/

void RigidBodySystem::initialize(ContextImpl& context, const vector<int>& bodyIndices, int rotationMode) {
    bodyIndex = cleanBodyIndices(bodyIndices);
    this->rotationMode = rotationMode;

    numBodies = 0;
    for (auto index : bodyIndex)
        numBodies = std::max(numBodies, index);

    const System& system = context.getSystem();
    int numAtoms = system.getNumParticles();
    numActualAtoms = numAtoms;
    for (int i = 0; i < numAtoms; i++)
        if (system.isVirtualSite(i))
            numActualAtoms--;
    atomIndex.resize(numActualAtoms);

    numFree = 0;
    body.resize(numBodies);
    for (int i = 0; i < numAtoms; i++)
        if (!(system.isVirtualSite(i) || system.getParticleMass(i) == 0.0)) {
            int ibody = bodyIndex[i];
            if (ibody == 0)
                atomIndex[numFree++] = i;
            else
                body[ibody-1].N++;
        }
    numBodyAtoms = numActualAtoms - numFree;
    bodyFixedPositions.resize(numBodyAtoms);

    freeAtom.resize(numFree);
    int loc = 0;
    for (auto& a : freeAtom) {
        a.index = atomIndex[loc++];
        a.invMass = 1.0/system.getParticleMass(a.index);
    }

    loc = 0;
    for (auto& b : body) {
        b.loc = loc;
        b.atom = &atomIndex[numFree+loc];
        b.d = &bodyFixedPositions[loc];
        loc += b.N;
    }

    vector<int> iatom(numBodies, 0);
    for (int i = 0; i < numAtoms; i++) {
        int ibody = bodyIndex[i];
        if (ibody > 0) {
            body[ibody-1].atom[iatom[ibody-1]++] = i;
        }
    }

    for (int i = 0; i < system.getNumConstraints(); i++) {
        int atom1, atom2;
        double distance;
        system.getConstraintParameters(i, atom1, atom2, distance);
        if (bodyIndex[atom1] != 0 || bodyIndex[atom2] != 0)
            throw OpenMMException("Constraints involving rigid-body atoms are not allowed");
    }
}

/*--------------------------------------------------------------------------------------------------
  Update the kinematic properties of all rigid bodies.
--------------------------------------------------------------------------------------------------*/

void RigidBodySystem::update(ContextImpl& context, bool geometry, bool velocities) {
    const System& system = context.getSystem();
    int N = system.getNumParticles();
    vector<double> M(N);
    for (int i = 0; i < N; i++)
        M[i] = system.getParticleMass(i);
    if (geometry) {
        vector<Vec3> R(N), F(N);
        context.getPositions(R);
        context.getForces(F);
        numDOF = numFree - system.getNumConstraints();
        for (auto& b : body) {
            b.buildGeometry(R, F, M);
            numDOF += b.dof;
        }
    }
    if (velocities) {
        vector<Vec3> V(N);
        context.getVelocities(V);
        for (auto& b : body)
            b.buildDynamics(V, M);
    }
}

/*--------------------------------------------------------------------------------------------------
  Copy structure from another rigid body system.
--------------------------------------------------------------------------------------------------*/

void RigidBodySystem::copy(const RigidBodySystem& bodySystem) {
    rotationMode = bodySystem.rotationMode;
    numFree = bodySystem.numFree;
    numBodies = bodySystem.numBodies;
    numActualAtoms = bodySystem.numActualAtoms;
    numBodyAtoms = bodySystem.numBodyAtoms;
    numDOF = bodySystem.numDOF;
    bodyIndex = bodySystem.bodyIndex;
    atomIndex = bodySystem.atomIndex;
    bodyFixedPositions = bodySystem.bodyFixedPositions;
    body = bodySystem.body;
    for (auto& b : body) {
        b.atom = &atomIndex[numFree+b.loc];
        b.d = &bodyFixedPositions[b.loc];
    }
}

/*--------------------------------------------------------------------------------------------------
  First part of rigid body NVE integration step.
--------------------------------------------------------------------------------------------------*/

void RigidBodySystem::integratePart1(double dt, const vector<Vec3>& F, vector<Vec3>& V, vector<Vec3>& R) {
    double halfDt = 0.5*dt;
    for (auto& a : freeAtom) {
        V[a.index] += F[a.index]*a.invMass*halfDt;
        R[a.index] += V[a.index]*dt;
        a.savedPos = R[a.index];
    }
    for (auto& b : body) {
        b.pcm += b.force*halfDt;
        b.pi += b.torque*dt;
        b.rcm += b.pcm*(b.invMass*dt);
        if (rotationMode == 0)
            b.exactRotation(dt);
        else
            b.noSquishRotation(dt, rotationMode);
        b.updateAtomicPositions(R);
    }
}

/*--------------------------------------------------------------------------------------------------
  Second part of rigid body NVE integration step.
--------------------------------------------------------------------------------------------------*/

void RigidBodySystem::integratePart2(double dt, const vector<Vec3>& R, const vector<Vec3>& F, vector<Vec3>& V) {
    double halfDt = 0.5*dt;
    double invDt = 1.0/dt;
    for (auto& a : freeAtom)
        V[a.index] += F[a.index]*a.invMass*halfDt + (R[a.index] - a.savedPos)*invDt;
    for (auto& b : body) {
        b.forceAndTorque(F);
        b.pcm += b.force*halfDt;
        b.pi += b.torque*dt;
        b.updateAtomicVelocities(V);
    }
}

/*--------------------------------------------------------------------------------------------------
  Compute kinetic energy.
--------------------------------------------------------------------------------------------------*/

void RigidBodySystem::computeKineticEnergies(const vector<Vec3>& V) {
    transKE = rotKE = 0.0;
    for (auto& a : freeAtom)
        transKE += V[a.index].dot(V[a.index])/a.invMass;
    for (auto& b : body) {
        transKE += b.twoKt;
        rotKE += b.twoKr;
    }
    transKE *= 0.5;
    rotKE *= 0.5;
}
