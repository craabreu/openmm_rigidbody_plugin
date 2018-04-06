/* -------------------------------------------------------------------------- *
 *                          OpenMM Rigid Body Plugin                          *
 * -------------------------------------------------------------------------- */

#include "RigidBodySystem.h"
#include "internal/MatVec.h"
#include "openmm/Vec3.h"
#include "openmm/System.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "internal/eigendecomposition.h"
#include <vector>
#include <math.h>

#include <iostream> // TEMPORARY

using namespace RigidBodyPlugin;
using namespace OpenMM;
using std::vector;

/*--------------------------------------------------------------------------------------------------
  Compute quaternion components from a given rotation matrix (Ref: S. W. Shepperd, Journal of
  Guidance and Control 1, p. 223, 1978).
--------------------------------------------------------------------------------------------------*/

Vec4 quaternion( const Mat3& A ) {
    double a11 = A[0][0];
    double a22 = A[1][1];
    double a33 = A[2][2];
    Vec4 Q2(1.0+a11+a22+a33, 1.0+a11-a22-a33, 1.0-a11+a22-a33, 1.0-a11-a22+a33);
    int imax = Q2.maxloc();
    double Q2max = Q2[imax];
    double factor = 0.5/sqrt(Q2max);
    if (imax == 0)
        return Vec4(Q2max, A[1][2]-A[2][1], A[2][0]-A[0][2], A[0][1]-A[1][0])*factor;
    else if (imax == 1)
        return Vec4(A[1][2]-A[2][1], Q2max, A[0][1]+A[1][0], A[0][2]+A[2][0])*factor;
    else if (imax == 2)
        return Vec4(A[2][0]-A[0][2], A[0][1]+A[1][0], Q2max, A[1][2]+A[2][1])*factor;
    else
        return Vec4(A[0][1]-A[1][0], A[0][2]+A[2][0], A[1][2]+A[2][1], Q2max)*factor;
}

/*--------------------------------------------------------------------------------------------------
  Premultiply matrix B(q) by a vector x
--------------------------------------------------------------------------------------------------*/

Vec4 multiplyB(Vec4 q, Vec3 x) {
    return Vec4(-q[1]*x[0] - q[2]*x[1] - q[3]*x[2],
                 q[0]*x[0] - q[3]*x[1] + q[2]*x[2],
                 q[3]*x[0] + q[0]*x[1] - q[1]*x[2],
                -q[2]*x[0] + q[1]*x[1] + q[0]*x[2]);
}

/*--------------------------------------------------------------------------------------------------
  Premultiply matrix C(q) by a vector x
--------------------------------------------------------------------------------------------------*/

Vec4 multiplyC(Vec4 q, Vec3 x) {
    return Vec4(-q[1]*x[0] - q[2]*x[1] - q[3]*x[2],
                 q[0]*x[0] + q[3]*x[1] - q[2]*x[2],
                -q[3]*x[0] + q[0]*x[1] + q[1]*x[2],
                 q[2]*x[0] - q[1]*x[1] + q[0]*x[2]);
}

/*--------------------------------------------------------------------------------------------------
  Premultiply matrix B^t(q) by a quaternion y
--------------------------------------------------------------------------------------------------*/

Vec3 multiplyBt(Vec4 q, Vec4 y) {
    return Vec3(-q[1]*y[0] + q[0]*y[1] + q[3]*y[2] - q[2]*y[3],
                -q[2]*y[0] - q[3]*y[1] + q[0]*y[2] + q[1]*y[3],
                -q[3]*y[0] + q[2]*y[1] - q[1]*y[2] + q[0]*y[3]);
}

/*--------------------------------------------------------------------------------------------------
  Premultiply matrix C^t(q) by a quaternion y
--------------------------------------------------------------------------------------------------*/

Vec3 multiplyCt(Vec4 q, Vec4 y) {
    return Vec3(-q[1]*y[0] + q[0]*y[1] - q[3]*y[2] + q[2]*y[3],
                -q[2]*y[0] + q[3]*y[1] + q[0]*y[2] - q[1]*y[3],
                -q[3]*y[0] - q[2]*y[1] + q[1]*y[2] + q[0]*y[3]);
}

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

    // Center-of-mass displacements and upper-triangular inertia tensor
    vector<Vec3> delta(N);
    Mat3 inertia;
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

    // Principal moments of inertia, rotation matrix, and orientational quaternion
    Mat3 A;
    eigendecomposition(inertia, A, MoI);
    A = A.t();
    q = quaternion(A);

    // Atom positions in the body-fixed frame of reference
    for (int i = 0; i < N; i++)
        d[i] = A*delta[i];
}

/*--------------------------------------------------------------------------------------------------
  Update linear and angular velocity based on individual atomic velocities. If necessary, also
  update these atomic velocities so as to eliminate central components.
--------------------------------------------------------------------------------------------------*/

void RigidBody::updateVelocities(vector<Vec3>& R, vector<Vec3>& V, vector<double>& M) {

    cout<<"velocities\n";
    // Compute total kinetic energy and center-of-mass velocity
    double K = 0.0;
    vcm = Vec3();
    for (int j = 0; j < N; j++) {
        int i = atom[j];
        double m = M[i];
        Vec3 p = V[i]*m;
        K += V[i].dot(V[i])*(0.5*m);
        vcm += p;
    }
    vcm /= mass;

    // Compute angular velocity via least square regression of V_i = vcm + omega x (R_i - rcm)
    /*                      |S(d_1)|                              |vcm-V_1|
       -|S(d_1) ... S(d_N)| | ...  | omega = -|S(d_1) ... S(d_N)| |  ...  |
                            |S(d_N)|                              |vcm-V_N|
        -sum_i{S^2(d_i)} omega = sum_i{S(d_i)(V_i-vcm) = sum_i{d_i x (V_i-vcm)}
        A omega = sum_i{d_i x (V_i-vcm)}, where A = -sum_i{S^2(d_i)} is symmetric */
    Mat3 A;
    Vec3 b;
    for (int j = 0; j < N; j++) {
        int i = atom[j];
        Vec3 d = R[i] - rcm;
        for (int k = 0; k < 3; k++) {
            A[k][k] += d.dot(d);
            for (int l = k; l < 3; l++)
                A[k][l] -= d[k]*d[l];
        }
        b += d.cross(V[i] - vcm);
        Mat3 Q;
        Vec3 lambda;
        eigendecomposition(A, Q, lambda);
//        omega
//        omega = Q.t()*b/lambda;
//        double Q[3][3], lambda[3];
//        int result = dsyevh3(A, Q, lambda); // eigendecomposition
//        omega = matVec(Q, vecDiv(transMatVec(Q, b), lambda));
    }
    // TODO: Solve the symmetric linear system A*omega = b, noting that only the upper triangular
    // part of A has been computed above


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
        b.updateVelocities(positions, velocities, masses);
    }
}


/*--------------------------------------------------------------------------------------------------
  Update the linear and angular velocities of all rigid bodies.
--------------------------------------------------------------------------------------------------*/

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
