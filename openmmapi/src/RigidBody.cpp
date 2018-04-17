/* ---------------------------------------------------------------------------------------------- *
 *                                    OpenMM Rigid Body Plugin                                    *
 * ---------------------------------------------------------------------------------------------- */

#include "RigidBody.h"
#include "internal/MatVec.h"
#include "openmm/Vec3.h"
#include "openmm/OpenMMException.h"
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

void RigidBody::updateGeometry(vector<Vec3>& R, vector<Vec3>& F, vector<double>& M) {

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

    // Principal moments of inertia, rotation matrix, and quaternion
    Vec3 u;
    Mat3 A;
    if (collinear(delta, d2, u)) {
        double I = 0.0;
        for (int j = 0; j < N; j++)
            I += M[atom[j]]*d2[j];
        MoI = Vec3(I, I, 0.0);
        invMoI = Vec3(1.0/I, 1.0/I, 0.0);
        Vec3 v = orthonormal(u);
        A = Mat3(v, u.cross(v), u).t();
        dof = 5;
    }
    else {
        Mat3 inertia;
        for (int j = 0; j < N; j++)
            inertia += Projection(delta[j])*M[atom[j]];
        MoI = eigenvalues(inertia);
        invMoI = Vec3(1.0/MoI[0], 1.0/MoI[1], 1.0/MoI[2]);
        A = eigenvectors(inertia, MoI);
        dof = 6;
    }
    q = Quat(A);

    // Atom positions in the body-fixed frame of reference
    for (int j = 0; j < N; j++)
        d[j] = A*delta[j];

    // Resultant force and quaternion-frame resultant torque
    force = Vec3();
    Vec3 tau;
    for (int j = 0; j < N; j++) {
        int i = atom[j];
        force += F[i];
        tau += delta[j].cross(F[i]);
    }
    torque = q.C(tau);
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
