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
#include <math.h>

#include <iostream> // TEMPORARY

using namespace RigidBodyPlugin;
using namespace OpenMM;
using std::vector;

/*--------------------------------------------------------------------------------------------------
  Update the geometric and kinematic properties of a rigid body based on the positions and
  velocities of individual atoms.
--------------------------------------------------------------------------------------------------*/

void RigidBody::update(vector<Vec3>& R, vector<Vec3>& V, vector<double>& M) {

    // Total mass and center-of-mass position
    mass = 0.0;
    rcm = Vec3();
    for (int j = 0; j < N; j++) {
        int i = atom[j];
        mass += M[i];
        rcm += R[i]*M[i];
    }
    rcm /= mass;
    invMass = 1.0/mass;

    // Center-of-mass displacements
    vector<Vec3> delta(N);
    for (int j = 0; j < N; j++)
        delta[j] = R[atom[j]] - rcm;

    // Inertia tensor
    Mat3 inertia;
    for (int j = 0; j < N; j++)
        inertia += Projection(delta[j])*M[atom[j]];

    // Principal moments of inertia, rotation matrix, and orientation quaternion
    Mat3 A;
    eigenDecomposition(inertia, A, MoI);
    A = A.t();
    q = Quat(A);
    invMoI = Vec3(1.0/MoI[0], 1.0/MoI[1], 1.0/MoI[2]);

    // Atom positions in the body-fixed frame of reference
    for (int j = 0; j < N; j++)
        d[j] = A*delta[j];
}

/*--------------------------------------------------------------------------------------------------
  Update linear and angular velocity based on individual atomic velocities. If necessary, also
  update these atomic velocities so as to eliminate central components.
--------------------------------------------------------------------------------------------------*/

void RigidBody::updateVelocities(vector<Vec3>& V, vector<double>& M) {

    // Total kinetic energy and center-of-mass velocity
    double K = 0.0;
    Vec3 pcm;
    for (int j = 0; j < N; j++) {
        int i = atom[j];
        Vec3 p = V[i]*M[i];
        K += 0.5*p.dot(V[i]);
        pcm += p;
    }
    vcm = pcm/mass;
    double Kt = 0.5*pcm.dot(vcm);

    // Body-fixed angular velocity
    Vec3 L;
    for (int j = 0; j < N; j++) {
        int i = atom[j];
        L += d[j].cross(q.A(V[i] - vcm)*M[i]);
    }
    Vec3 omega = Diag3(invMoI)*L;
    double Kr = 0.5*L.dot(omega);
}

/*--------------------------------------------------------------------------------------------------
  Create a clean list of rigid body indices. The rule is: if bodyIndices[i] <= 0 or is unique, then
  index[i] = 0. Otherwise, 1 <= index[i] <= number of distinct positive entries which are
  non-unique.
--------------------------------------------------------------------------------------------------*/

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

/*--------------------------------------------------------------------------------------------------
  Create a data structure for the system of rigid bodies and free atoms.
--------------------------------------------------------------------------------------------------*/

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

/*--------------------------------------------------------------------------------------------------
  Update the kinematic properties of all rigid bodies.
--------------------------------------------------------------------------------------------------*/

void RigidBodySystem::update() {
    const System* system = &context->getSystem();
    int N = system->getNumParticles();
    vector<Vec3> positions(N), velocities(N);
    vector<double> masses(N);
    context->getPositions(positions);
    context->getVelocities(velocities);
    for (int i = 0; i < N; i++)
      masses[i] = system->getParticleMass(i);
    for (auto&b : body) {
        b.update(positions, velocities, masses);
//        b.updateVelocities(positions, velocities, masses);
    }
}


/*--------------------------------------------------------------------------------------------------
  Update the linear and angular velocities of all rigid bodies.
--------------------------------------------------------------------------------------------------*/

void RigidBodySystem::updateVelocities() {
    const System* system = &context->getSystem();
    int N = system->getNumParticles();
    vector<Vec3> positions(N), velocities(N);
    vector<double> masses(N);
    context->getPositions(positions);
    context->getVelocities(velocities);
    for (int i = 0; i < N; i++)
      masses[i] = system->getParticleMass(i);
    for (auto&b : body) {
//        b.update(positions, velocities, masses);
        b.updateVelocities(velocities, masses);
    }
}
