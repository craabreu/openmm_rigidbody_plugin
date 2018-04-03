/* -------------------------------------------------------------------------- *
 *                          OpenMM Rigid Body Plugin                          *
 * -------------------------------------------------------------------------- */

#include "RigidBodySystem.h"
#include "internal/Vec4.h"
#include "openmm/Vec3.h"
#include "openmm/System.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "internal/diagonalization.h"
#include <vector>

#include <iostream> // TEMPORARY

using namespace RigidBodyPlugin;
using namespace OpenMM;
using std::vector;

void RigidBody::update(vector<Vec3>& R, vector<Vec3>& V, vector<double>& M) {

    // Compute total mass and center-of-mass position
    mass = 0.0;
    rcm = Vec3();
    for (int j = 0; j < N; j++) {
        int i = atom[j];
        mass += M[i];
        rcm += R[i]*M[i];
    }
    rcm /= mass;

    // Compute center-of-mass displacements and upper-triangular inertia tensor
    vector<Vec3> delta(N);
    double inertia[3][3] = {0.0};
    for (int j = 0; j < N; j++) {
        int i = atom[j];
        Vec3 disp = R[i] - rcm;
        delta[j] = disp;
        inertia[0][0] += M[i]*(disp[1]*disp[1] + disp[2]*disp[2]);
        inertia[1][1] += M[i]*(disp[0]*disp[0] + disp[2]*disp[2]);
        inertia[2][2] += M[i]*(disp[0]*disp[0] + disp[1]*disp[1]);
        inertia[0][1] -= M[i]*disp[0]*disp[1];
        inertia[0][2] -= M[i]*disp[0]*disp[2];
        inertia[1][2] -= M[i]*disp[1]*disp[2];
    }

    // Compute rotation matrix and moments of inertia
    double A[3][3], I[3];
    int result = dsyevh3(inertia, A, I);
    if (result != 0)
        throw OpenMMException("Diagonalization of rigid body inertia tensor failed");
    MoI = Vec3(I[0], I[1], I[2]);

    cout<<MoI<<"\n";
    // TODO: Compute the quaternions from the rotation matrices:
    // TODO: Compute the body-fixed atom positions
}

/*
 * Update linear and angular velocity based on individual atomic velocities. If necessary,
 * also update these atomic velocities so as to eliminate central components.
 */

void RigidBody::updateVelocities(vector<Vec3>& V, vector<double>& M) {
    vcm = Vec3();
    for (int j = 0; j < N; j++) {
        int i = atom[j];
        vcm += V[i]*M[i];
    }
    vcm /= mass;

    // TODO: Update angular velocities
    // TODO: Eliminate central components
}

/*
 * Create a clean list of rigid body indices. The rule is: if bodyIndices[i] <= 0 or
 * bodyIndices[i] is unique, then index[i] = 0. Otherwise, 1 <= index[i] <= number of
 * distinct positive entries which are non-unique.
 */

vector<int> cleanBodyIndices(const vector<int>& bodyIndices) {
    int size = bodyIndices.size();
    vector<int> index(size, 0), saved(size, 0), amount(size, 0), first(size, 0);
    int nsaved = 0;
    for (int i = 0; i < size; i++) {
        int ibody = bodyIndices[i];
        if (ibody > 0) {
            int j;
            int found = 0;
            for (j = 0; j < nsaved; j++) {
                found = saved[j] == ibody;
                if (found) break;
            }
            if (found) {
                amount[j]++;
                index[i] = j+1;
                index[first[j]] = j+1;
            }
            else {
                amount[nsaved] = 1;
                saved[nsaved] = ibody;
                first[nsaved] = i;
                index[i] = 0;
                nsaved++;
            }
        }
        else
            index[i] = 0;
    }
    int n = 0;
    for (int i = 0; i < nsaved; i++)
        saved[i] = amount[i] > 1 ? ++n : 0;
    for (int i = 0; i < size; i++)
        if (index[i] > 0)
            index[i] = saved[index[i]-1];
    return index;
}

/*
 * Create a data structure for the system of rigid bodies and free atoms
 */

RigidBodySystem::RigidBodySystem(ContextImpl& contextRef, const vector<int>& bodyIndices) {
    context = &contextRef;
    bodyIndex = cleanBodyIndices(bodyIndices);

    numBodies = 0;
    for (auto index : bodyIndex)
        numBodies = std::max(numBodies, index);

    const System *system = &context->getSystem();
    int numAtoms = system->getNumParticles();
    numActualAtoms = numAtoms;
    for (int i = 0; i < numAtoms; i++)
        if (system->isVirtualSite(i))
            numActualAtoms--;
    atomIndex.resize(numActualAtoms);

    int numFree = 0;
    body.resize(numBodies);
    for (int i = 0; i < numAtoms; i++)
        if (!system->isVirtualSite(i)) {
            int ibody = bodyIndex[i];
            if (ibody == 0)
                atomIndex[numFree++] = i;
            else
                body[ibody-1].N++;
        }
    numBodyAtoms = numActualAtoms - numFree;
    bodyFixedPositions.resize(numBodyAtoms);

    int loc = 0;
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

    cout<<"Number of bodies = "<<numBodies<<"\n"
        <<"Number of actual atoms = "<<numActualAtoms<<"\n"
        <<"Number of free atoms = "<<numFree<<"\n";
}

/*
 * Update the kinematic properties of all rigid bodies
 */

void RigidBodySystem::update() {
    const System* system = &context->getSystem();
    int N = system->getNumParticles();
    vector<Vec3> positions(N), velocities(N);
    vector<double> masses(N);
    context->getPositions(positions);
    context->getVelocities(velocities);
    for (int i = 0; i < N; i++)
      masses[i] = system->getParticleMass(i);
    for (auto&b : body)
        b.update(positions, velocities, masses);
}


/*
 * Update the linear and angular velocities of all rigid bodies
 */

void RigidBodySystem::updateVelocities() {
//    const System* system = &context->getSystem();
//    int N = system->getNumParticles();
//    vector<Vec3> positions(N), velocities(N);
//    vector<double> masses(N);
//    context->getPositions(positions);
//    context->getVelocities(velocities);
//    for (int i = 0; i < N; i++)
//      masses[i] = system->getParticleMass(i);
//    for (auto&b : body)
//        b.update(positions, velocities, masses, atomIndex);
}
