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

#include <iostream> // TEMPORARY

using namespace RigidBodyPlugin;
using namespace OpenMM;
using std::vector;

/*--------------------------------------------------------------------------------------------------
  Check whether all 3D vectors of a given set lie in the same straight line. If so, return a unit
  vector denoting the direction of such line.
--------------------------------------------------------------------------------------------------*/

bool collinear(vector<Vec3>& delta, vector<double>& d2, Vec3& u) {
    const double TOL = 1.0E-5;
    double d0d0 = d2[0];
    int jmax = 0;
    double d2max = d0d0;
    int N = delta.size();
    for (int j = 1; j < N; j++)
        if (d2[j] > d0d0) {
            jmax = j;
            d2max = d2[j];
        }
    u = delta[jmax]/sqrt(d2max);
    bool isCollinear = true;
    for (int j = 0; isCollinear && j < N; j++) {
        double djdj = d2[j];
        double udj = u.dot(delta[j]);
        isCollinear = isCollinear && (djdj < TOL*d2max || abs(udj*udj/djdj - 1.0) < TOL);
    }
    return isCollinear;
}

/*--------------------------------------------------------------------------------------------------
  Return a unit vector orthogonal to a given vector u
--------------------------------------------------------------------------------------------------*/

Vec3 orthonormal(Vec3& u) {
    int imin = u[0] < u[1] ? 0 : 1;
    if (u[2] < u[imin]) imin = 2;
    Vec3 v;
    v[imin] = 1.0;
    v = Projection(u)*v;
    return v/sqrt(v.dot(v));
}

/*--------------------------------------------------------------------------------------------------
  Update the geometric and kinematic properties of a rigid body based on the positions and
  velocities of individual atoms.
--------------------------------------------------------------------------------------------------*/

void RigidBody::updateGeometry(vector<Vec3>& R, vector<double>& M) {

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
    vector<double> d2(N);
    for (int j = 0; j < N; j++) {
        delta[j] = R[atom[j]] - rcm;
        d2[j] = delta[j].dot(delta[j]);
    }

    // Principal moments of inertia and inverse rotation matrix
    Vec3 u;
    Mat3 At;
    if (collinear(delta, d2, u)) {
        double I = 0.0;
        for (int j = 0; j < N; j++)
            I += M[atom[j]]*d2[j];
        MoI = Vec3(I, I, 0.0);
        invMoI = Vec3(1.0/I, 1.0/I, 0.0);
        Vec3 v = orthonormal(u);
        At = Mat3(v, v.cross(u), u);
        dof = 5;
    }
    else {
        Mat3 inertia;
        for (int j = 0; j < N; j++)
            inertia += Projection(delta[j])*M[atom[j]];
        eigenDecomposition(inertia, At, MoI);
        invMoI = Vec3(1.0/MoI[0], 1.0/MoI[1], 1.0/MoI[2]);
        dof = 6;
    }

    // Rotation matrix
    Mat3 A = At.t();
    q = Quat(A);

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
    for (int j = 0; j < N; j++) {
        int i = atom[j];
        Vec3 p = V[i]*M[i];
        pcm += p;
    }
    Vec3 vcm = pcm/mass;

    // Quaternion-conjugated momentum
    Vec3 L;
    for (int j = 0; j < N; j++) {
        int i = atom[j];
        L += d[j].cross(q.A(V[i] - vcm)*M[i]);
    }
    pi = q.B(L)*2.0;
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

void RigidBodySystem::update(bool geometry, bool velocities) {
    const System* system = &context->getSystem();
    int N = system->getNumParticles();
    vector<double> M(N);
    for (int i = 0; i < N; i++)
        M[i] = system->getParticleMass(i);
    if (geometry) {
        vector<Vec3> R(N);
        context->getPositions(R);
        for (auto& b : body)
            b.updateGeometry(R, M);
    }
    if (velocities) {
        vector<Vec3> V(N);
        context->getVelocities(V);
        for (auto& b : body)
            b.updateVelocities(V, M);
    }
}
