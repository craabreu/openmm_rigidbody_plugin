/* ---------------------------------------------------------------------------------------------- *
 *                                    OpenMM Rigid Body Plugin                                    *
 * ---------------------------------------------------------------------------------------------- */

#include "RigidBody.h"
#include "internal/MatVec.h"
#include "openmm/Vec3.h"
#include "openmm/OpenMMException.h"
#include "internal/eigenDecomposition.h"
#include "internal/ellipticFunctions.h"
#include <vector>
#include <cmath>
#include <limits>

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

void RigidBody::buildGeometry(vector<Vec3>& R, vector<Vec3>& F, vector<double>& M) {

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
    delta.resize(N);
    vector<double> d2(N);
    for (int j = 0; j < N; j++) {
        delta[j] = R[atom[j]] - rcm;
        d2[j] = delta[j].dot(delta[j]);
    }

    // Principal moments of inertia, rotation matrix, and quaternion
    Vec3 u;
    Mat3 A;
    if (collinear(delta, d2, u)) {
        double MoI = 0.0;
        for (int j = 0; j < N; j++)
            MoI += M[atom[j]]*d2[j];
        I = Vec3(MoI, MoI, 0.0);
        invI = Vec3(1.0/MoI, 1.0/MoI, 0.0);
        Vec3 v = orthonormal(u);
        A = Mat3(v, u.cross(v), u).t();
        dof = 5;
    }
    else {
        Mat3 inertia;
        for (int j = 0; j < N; j++)
            inertia += Projection(delta[j])*M[atom[j]];
        I = eigenvalues(inertia);
        invI = Vec3(1.0/I[0], 1.0/I[1], 1.0/I[2]);
        A = eigenvectors(inertia, I);
        dof = 6;
    }
    q = Quat(A);

    // Atom positions in the body-fixed frame of reference
    for (int j = 0; j < N; j++)
        d[j] = A*delta[j];

    // Resultant force and quaternion-frame resultant torque
    forceAndTorque(F);
}

/*--------------------------------------------------------------------------------------------------
  Update linear and angular velocity based on individual atomic velocities. If necessary, also
  update these atomic velocities so as to eliminate central components.
--------------------------------------------------------------------------------------------------*/

void RigidBody::buildDynamics(vector<Vec3>& V, vector<double>& M) {

    // Center-of-mass velocity
    for (int j = 0; j < N; j++) {
        int i = atom[j];
        Vec3 p = V[i]*M[i];
        pcm += p;
    }
    Vec3 vcm = pcm/mass;
    twoKt = pcm.dot(vcm);

    // Quaternion-conjugated momentum
    Vec3 L;
    for (int j = 0; j < N; j++) {
        int i = atom[j];
        L += d[j].cross(q.A(V[i] - vcm)*M[i]);
    }
    pi = q.B(L)*2.0;
    twoKr = L.dot(Diag3(invI)*L);
}

/*--------------------------------------------------------------------------------------------------
  Update the positions atoms based on current rigid-body state.
--------------------------------------------------------------------------------------------------*/

void RigidBody::updateAtomicPositions(vector<Vec3>& R) {
    for (int j = 0; j < N; j++) {
        delta[j] = q.At(d[j]);
        R[atom[j]] = rcm + delta[j];
    }
}

/*--------------------------------------------------------------------------------------------------
  Update the velocities atoms based on current rigid-body state.
--------------------------------------------------------------------------------------------------*/

void RigidBody::updateAtomicVelocities(vector<Vec3>& V) {
    Vec3 L = q.Bt(pi)*0.5;
    Vec3 omega = Diag3(invI)*L;
    Vec3 spaceFixedOmega = q.At(omega);
    Vec3 vcm = pcm*invMass;
    for (int j = 0; j < N; j++)
        V[atom[j]] = vcm + spaceFixedOmega.cross(delta[j]);
    twoKt = pcm.dot(vcm);
    twoKr = L.dot(omega);
}

/*--------------------------------------------------------------------------------------------------
  Update resultant force and torque exerted on the body.
--------------------------------------------------------------------------------------------------*/

void RigidBody::forceAndTorque(const vector<Vec3>& F) {
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
  Perform rotation around a principal axis
--------------------------------------------------------------------------------------------------*/

void RigidBody::uniaxialRotationAxis1(double dt) {
    Quat Bq = q.B1();
    double omegaDtBy2 = 0.25*pi.dot(Bq)*dt*invI[0];
    double vsin = sin(omegaDtBy2);
    double vcos = cos(omegaDtBy2);
    q = q*vcos + Bq*vsin;
    pi = pi*vcos + pi.B1()*vsin;
}

void RigidBody::uniaxialRotationAxis2(double dt) {
    Quat Bq = q.B2();
    double omegaDtBy2 = 0.25*pi.dot(Bq)*dt*invI[1];
    double vsin = sin(omegaDtBy2);
    double vcos = cos(omegaDtBy2);
    q = q*vcos + Bq*vsin;
    pi = pi*vcos + pi.B2()*vsin;
}

void RigidBody::uniaxialRotationAxis3(double dt) {
    Quat Bq = q.B3();
    double omegaDtBy2 = 0.25*pi.dot(Bq)*dt*invI[2];
    double vsin = sin(omegaDtBy2);
    double vcos = cos(omegaDtBy2);
    q = q*vcos + Bq*vsin;
    pi = pi*vcos + pi.B3()*vsin;
}

/*--------------------------------------------------------------------------------------------------
  Perform a torque-free rotation using the NO_SQUISH method.
--------------------------------------------------------------------------------------------------*/

void RigidBody::noSquishRotation(double dt, int n) {
    double dtByN = dt/n;
    double halfDtByN = 0.5*dtByN;
    bool axis3 = dof == 6;
    for (int i = 0; i < n; i++) {
        if (axis3) uniaxialRotationAxis3(halfDtByN);
        uniaxialRotationAxis2(halfDtByN);
        uniaxialRotationAxis1(dtByN);
        uniaxialRotationAxis2(halfDtByN);
        if (axis3) uniaxialRotationAxis3(halfDtByN);
    }
}

/*--------------------------------------------------------------------------------------------------
  Perform a torque-free rotation using exact solution.
--------------------------------------------------------------------------------------------------*/

#define SIGN(x)      ((x) >= 0 ? 1 : -1)
#define stairCase(x) ((x) > 0 ? (int)ceil((x) - 0.5) : (int)floor((x) + 0.5))
#define PI           3.14159265358979323846264338328

void RigidBody::exactRotation(double dt) {
    const double EPSILON = numeric_limits<double>::epsilon();
    Vec3 Iw = q.Bt(pi)*0.5;
    Vec3 w0 = Diag3(invI)*Iw;
    double Lsq = Iw[1]*Iw[1] + Iw[2]*Iw[2];
    if (Lsq < EPSILON)
        return uniaxialRotationAxis1(dt);
    Lsq += Iw[0]*Iw[0];
    double L = sqrt(Lsq);
    double twoKr = Iw.dot(w0);
    Quat z0(Iw[2], Iw[1], L - Iw[0], 0.0);
    double r1 = Lsq - twoKr*I[2];
    double r3 = twoKr*I[0] - Lsq;
    double l1 = r1*invI[1]/(I[1] - I[2]);
    double l3 = r3*invI[1]/(I[0] - I[1]);
    double lmin = min(l1, l3);
    double c13 = 1.0/(I[0] - I[2]);
    Vec3 a(SIGN(w0[0])*sqrt(r1*invI[0]*c13), sqrt(lmin), SIGN(w0[2])*sqrt(r3*invI[2]*c13));
    double m = lmin/max(l1, l3);
    double K = carlsonRF(0.0, 1.0 - m, 1.0);
    double inv2K = 0.5/K;
    double s0 = w0[1]/a[1];
    double c0, u0;
    int i0;
    if (fabs(s0) < 1.0) {
        c0 = l1 < l3 ? w0[0]/a[0] : w0[2]/a[2];
        u0 = s0*carlsonRF(1.0 - s0*s0, 1.0 - m*s0*s0, 1.0);
        i0 = stairCase(u0*inv2K);
    }
    else {
        a[1] = fabs(w0[1]);
        s0 = SIGN(s0);
        c0 = 0.0;
        u0 = s0*K;
        i0 = 0;
    }
    double wp = -invI[1]*a[0]*a[2]/(a[1]*c13);
    double u = wp*dt + u0;
    int jump = stairCase(u*inv2K) - i0;
    double sn, cn, dn, deltaF;
    jacobi(u, m, sn, cn, dn);
    double alpha = I[0]*a[0]/L;
    double eta = alpha*alpha;
    eta /= 1.0 - eta;
    Diag3 Ia = Diag3(I)*Diag3(a);
    if (l1 < l3) {
        double C = sqrt(m + eta);
        deltaF = u - u0 + SIGN(cn)*Omega(sn, eta, m) - SIGN(c0)*Omega(s0, eta, m)
                         + (alpha/C)*(atan(C*sn/dn) - atan(C*s0*a[2]/w0[2]));
        if (jump != 0) deltaF += jump*2.0*Omega(1.0, eta, m);
        Iw = Ia*Vec3(cn, sn, dn);
    }
    else {
        double k2eta = m*eta;
        double C = sqrt(1.0 + k2eta);
        deltaF = u - u0 + SIGN(cn)*Omega(sn, k2eta, m) - SIGN(c0)*Omega(s0, k2eta, m)
                        + (alpha/C)*(atan(C*sn/cn) - atan(C*s0/c0));
        if (jump != 0) deltaF += jump*(2.0*Omega(1.0, k2eta, m) + (alpha/C)*PI);
        Iw = Ia*Vec3(dn, sn, cn);
    }
    deltaF *= 1.0 + eta;
    double theta = (Lsq*(u - u0) + r3*deltaF)/(2.0*L*I[0]*wp);
    Quat z = Quat( Iw[2], Iw[1], L - Iw[0], 0.0)*cos(theta) +
             Quat(-Iw[1], Iw[2], 0.0, L - Iw[0])*sin(theta);
    q = z*z0.dot(q) + z.C(z0.Ct(q));
    q /= q.norm();
    pi = q.B(Iw*2.0);
}
