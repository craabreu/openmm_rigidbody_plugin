#define zero    ((mixed)0.00)
#define quarter ((mixed)0.25)
#define half    ((mixed)0.50)
#define one     ((mixed)1.00)
#define two     ((mixed)2.00)

/*--------------------------------------------------------------------------------------------------
  Define struct for rigid body data.
--------------------------------------------------------------------------------------------------*/

extern "C" typedef struct {
    int    N;     // number of atoms
    int    loc;   // pointer to set of atoms
    mixed  invm;  // mass
    mixed3 invI;  // principal moments of inertia
    mixed3 r;     // center-of-mass position
    mixed3 v;     // center-of-mass velocity
    mixed3 F;     // resultant force
    mixed3 rdot;  // time-derivative of center-of-mass position
    mixed4 q;     // orientation quaternion
    mixed4 pi;    // quaternion-conjugated momentum
    mixed4 Ctau;  // quaternion-frame resultant torque
    mixed4 qdot;  // time-derivative of orientation quaternion
} BodyData;

/*--------------------------------------------------------------------------------------------------
  Load the position of a particle.
--------------------------------------------------------------------------------------------------*/

inline __device__ mixed3 loadPos(const real4* __restrict__ posq,
                                 const real4* __restrict__ posqCorrection,
                                 int index) {
    const real4& R = posq[index];
#ifdef USE_MIXED_PRECISION
    const real4& corr = posqCorrection[index];
    return make_mixed3(R.x + (mixed)corr.x, R.y + (mixed)corr.y, R.z + (mixed)corr.z);
#else
    return make_mixed3(R.x, R.y, R.z);
#endif
}

/*--------------------------------------------------------------------------------------------------
  Store the position of a particle.
--------------------------------------------------------------------------------------------------*/

inline __device__ void storePos(real4* __restrict__ posq,
                                real4* __restrict__ posqCorrection,
                                int index,
                                mixed3 pos) {
    real4& R = posq[index];
#ifdef USE_MIXED_PRECISION
    real4& corr = posqCorrection[index];
    R.x = (real)pos.x;
    R.y = (real)pos.y;
    R.z = (real)pos.z;
    corr.x = pos.x - (real)pos.x;
    corr.y = pos.y - (real)pos.y;
    corr.z = pos.z - (real)pos.z;
#else
    R.x = pos.x;
    R.y = pos.y;
    R.z = pos.z;
#endif
}

/*--------------------------------------------------------------------------------------------------
  Pre-multiplication with rotation matrix factors.
--------------------------------------------------------------------------------------------------*/

inline __device__ mixed4 B(mixed4 q, mixed3 v) {
    return make_mixed4(-q.y*v.x - q.z*v.y - q.w*v.z,
                        q.x*v.x - q.w*v.y + q.z*v.z,
                        q.w*v.x + q.x*v.y - q.y*v.z,
                       -q.z*v.x + q.y*v.y + q.x*v.z);
}

inline __device__ mixed4 C(mixed4 q, mixed3 v) {
    return make_mixed4(-q.y*v.x - q.z*v.y - q.w*v.z,
                        q.x*v.x + q.w*v.y - q.z*v.z,
                       -q.w*v.x + q.x*v.y + q.y*v.z,
                        q.z*v.x - q.y*v.y + q.x*v.z);
}

inline __device__ mixed3 Bt(mixed4 q, mixed4 v) {
    return make_mixed3(-q.y*v.x + q.x*v.y + q.w*v.z - q.z*v.w,
                       -q.z*v.x - q.w*v.y + q.x*v.z + q.y*v.w,
                       -q.w*v.x + q.z*v.y - q.y*v.z + q.x*v.w);
}

inline __device__ mixed3 Ct(mixed4 q, mixed4 v) {
    return make_mixed3(-q.y*v.x + q.x*v.y - q.w*v.z + q.z*v.w,
                       -q.z*v.x + q.w*v.y + q.x*v.z - q.y*v.w,
                       -q.w*v.x - q.z*v.y + q.y*v.z + q.x*v.w);
}

/*--------------------------------------------------------------------------------------------------
  Pre-multiplication with permutation matrices.
--------------------------------------------------------------------------------------------------*/

inline __device__ mixed4 B1(mixed4 q) {
    return make_mixed4(-q.y,  q.x,  q.w, -q.z);
}

inline __device__ mixed4 B2(mixed4 q) {
    return make_mixed4(-q.z, -q.w,  q.x,  q.y);
}

inline __device__ mixed4 B3(mixed4 q) {
    return make_mixed4(-q.w,  q.z, -q.y,  q.x);
}

/*--------------------------------------------------------------------------------------------------
  Uniaxial rotations
--------------------------------------------------------------------------------------------------*/

inline __device__ void uniaxialRotationAxis1(mixed4& q, mixed4& pi, mixed3 invI, mixed dt) {
    mixed4 Bq = B1(q);
    mixed omegaDtBy2 = quarter*dot(pi, Bq)*dt*invI.x;
    mixed vsin = sin(omegaDtBy2);
    mixed vcos = cos(omegaDtBy2);
    q = q*vcos + Bq*vsin;
    pi = pi*vcos + B1(pi)*vsin;
}

inline __device__ void uniaxialRotationAxis2(mixed4& q, mixed4& pi, mixed3 invI, mixed dt) {
    mixed4 Bq = B2(q);
    mixed omegaDtBy2 = quarter*dot(pi, Bq)*dt*invI.y;
    mixed vsin = sin(omegaDtBy2);
    mixed vcos = cos(omegaDtBy2);
    q = q*vcos + Bq*vsin;
    pi = pi*vcos + B2(pi)*vsin;
}

inline __device__ void uniaxialRotationAxis3(mixed4& q, mixed4& pi, mixed3 invI, mixed dt) {
    mixed4 Bq = B3(q);
    mixed omegaDtBy2 = quarter*dot(pi, Bq)*dt*invI.z;
    mixed vsin = sin(omegaDtBy2);
    mixed vcos = cos(omegaDtBy2);
    q = q*vcos + Bq*vsin;
    pi = pi*vcos + B3(pi)*vsin;
}

/*--------------------------------------------------------------------------------------------------
  NO_SQUISH rotation
--------------------------------------------------------------------------------------------------*/

inline __device__ void noSquishRotation(mixed4& q, mixed4& pi, mixed3 invI, mixed dt) {
    mixed dtByN = dt/NSPLIT;
    mixed halfDtByN = 0.5*dtByN;
    bool axis3 = invI.z != zero;
    for (int i = 0; i < NSPLIT; i++) {
        if (axis3) uniaxialRotationAxis3(q, pi, invI, halfDtByN);
        uniaxialRotationAxis2(q, pi, invI, halfDtByN);
        uniaxialRotationAxis1(q, pi, invI, dtByN);
        uniaxialRotationAxis2(q, pi, invI, halfDtByN);
        if (axis3) uniaxialRotationAxis3(q, pi, invI, halfDtByN);
    }
}

/*--------------------------------------------------------------------------------------------------
  Exact rotation
--------------------------------------------------------------------------------------------------*/

#define SIGN(x)      ((x) >= 0 ? 1 : -1)
#define stairCase(x) ((x) > 0 ? (int)ceil((x) - 0.5) : (int)floor((x) + 0.5))
#define PI           3.14159265358979323846264338328

inline __device__ void exactRotation(mixed4& q, mixed4& pi, mixed3 invI, mixed dt) {
    mixed3 I = make_mixed3(one/invI.x, one/invI.y, invI.z != zero ? one/invI.z : zero);
    mixed3 Iw = Bt(q, pi)*half;
    mixed3 w0 = invI*Iw;
    mixed Lsq = Iw.y*Iw.y + Iw.z*Iw.z;
    if (Lsq < EPSILON) return uniaxialRotationAxis1(q, pi, invI, dt);
    Lsq += Iw.x*Iw.x;
    mixed L = sqrt(Lsq);
    mixed twoKr = dot(Iw, w0);
    mixed4 z0 = make_mixed4(Iw.z, Iw.y, L - Iw.x, zero);
    mixed r1 = Lsq - twoKr*I.z;
    mixed r3 = twoKr*I.x - Lsq;
    mixed l1 = r1*invI.y/(I.y - I.z);
    mixed l3 = r3*invI.y/(I.x - I.y);
    mixed lmin = min(l1, l3);
    mixed c13 = one/(I.x - I.z);
    mixed3 a = make_mixed3(SIGN(w0.x)*sqrt(r1*invI.x*c13), sqrt(lmin), SIGN(w0.z)*sqrt(r3*invI.z*c13));
    mixed m = lmin/max(l1, l3);
    mixed K = carlsonRF(zero, one - m, one);
    mixed inv2K = half/K;
    mixed s0 = w0.y/a.y;
    mixed c0, u0;
    int i0;
    if (fabs(s0) < one) {
        c0 = l1 < l3 ? w0.x/a.x : w0.z/a.z;
        u0 = s0*carlsonRF(one - s0*s0, one - m*s0*s0, one);
        i0 = stairCase(u0*inv2K);
    }
    else {
        a.y = fabs(w0.y);
        s0 = SIGN(s0);
        c0 = zero;
        u0 = s0*K;
        i0 = 0;
    }
    mixed wp = -invI.y*a.x*a.z/(a.y*c13);
    mixed u = wp*dt + u0;
    int jump = stairCase(u*inv2K) - i0;
    mixed sn, cn, dn, deltaF;
    jacobi(u, m, sn, cn, dn);
    mixed alpha = I.x*a.x/L;
    mixed eta = alpha*alpha;
    eta /= one - eta;
    if (l1 < l3) {
        mixed C = sqrt(m + eta);
        deltaF = u - u0 + SIGN(cn)*Omega(sn, eta, m) - SIGN(c0)*Omega(s0, eta, m)
                         + (alpha/C)*(atan(C*sn/dn) - atan(C*s0*a.z/w0.z));
        if (jump != 0) deltaF += jump*two*Omega(one, eta, m);
        Iw = I*a*make_mixed3(cn, sn, dn);
    }
    else {
        mixed k2eta = m*eta;
        mixed C = sqrt(one + k2eta);
        deltaF = u - u0 + SIGN(cn)*Omega(sn, k2eta, m) - SIGN(c0)*Omega(s0, k2eta, m)
                        + (alpha/C)*(atan(C*sn/cn) - atan(C*s0/c0));
        if (jump != 0) deltaF += jump*(two*Omega(one, k2eta, m) + (alpha/C)*PI);
        Iw = I*a*make_mixed3(dn, sn, cn);
    }
    deltaF *= one + eta;
    mixed theta = (Lsq*(u - u0) + r3*deltaF)/(two*L*I.x*wp);
    mixed4 z = make_mixed4( Iw.z, Iw.y, L - Iw.x, zero)*cos(theta) +
               make_mixed4(-Iw.y, Iw.z, zero, L - Iw.x)*sin(theta);
    q = normalize(z*dot(z0, q) + C(z, Ct(z0, q)));
    pi = B(q, Iw*two);
}

/*--------------------------------------------------------------------------------------------------
  Virtual rotation.
--------------------------------------------------------------------------------------------------*/

inline __device__ mixed4 virtualRotation(BodyData& body, mixed dt) {
    mixed4 q = body.q;
    mixed4 pi = body.pi;
    body.pi += body.Ctau*dt;
    ROTATION(body, dt);
    mixed4 qnew = body.q;
    body.q = q;
    body.pi = pi;
    return qnew;
}

/*--------------------------------------------------------------------------------------------------
  Perform the initial step of Rigid Body integration.
--------------------------------------------------------------------------------------------------*/

extern "C" __global__ void integrateRigidBodyPart1(int numAtoms,
                                                   int stride,
                                                   int numFree,
                                                   int numBodies,
                                                   const mixed dt,
                                                   real4* __restrict__ posq,
                                                   real4* __restrict__ posqCorrection,
                                                   mixed4* __restrict__ velm,
                                                   const long long* __restrict__ force,
                                                   mixed4* __restrict__ posDelta,
                                                   BodyData* __restrict__ bodyData,
                                                   const int* __restrict__ atomLocation,
                                                   const mixed3* __restrict__ bodyFixedPos,
                                                   mixed3* __restrict__ savedPos) {

    const mixed scale = one/(mixed)0x100000000;
    const mixed halfDt = half*dt;
    for (int j = blockIdx.x*blockDim.x+threadIdx.x; j < numFree; j += blockDim.x*gridDim.x) {
        int i = atomLocation[j];
        mixed4& velocity = velm[i];
        mixed3 f = make_mixed3(force[i], force[i+stride], force[i+stride*2])*scale;
        mixed3 v = trim(velocity) + f*(velocity.w*halfDt);
        mixed3 delta = v*dt;
        velocity = fuse(v, velocity.w);
        posDelta[i] = fuse(delta, zero);
        savedPos[j] = loadPos(posq, posqCorrection, i) + delta;
    }
}

/*--------------------------------------------------------------------------------------------------
  Perform the final step of Rigid Body integration.
--------------------------------------------------------------------------------------------------*/

extern "C" __global__ void integrateRigidBodyPart2(int numAtoms,
                                                   int stride,
                                                   int numFree,
                                                   int numBodies,
                                                   const mixed dt,
                                                   real4* __restrict__ posq,
                                                   real4* __restrict__ posqCorrection,
                                                   mixed4* __restrict__ velm,
                                                   const long long* __restrict__ force,
                                                   mixed4* __restrict__ posDelta,
                                                   BodyData* __restrict__ bodyData,
                                                   const int* __restrict__ atomLocation,
                                                   const mixed3* __restrict__ bodyFixedPos,
                                                   mixed3* __restrict__ savedPos) {

    for (int j = blockIdx.x*blockDim.x+threadIdx.x; j < numFree; j += blockDim.x*gridDim.x) {
        int i = atomLocation[j];
        mixed3 pos = loadPos(posq, posqCorrection, i) + trim(posDelta[i]);
        storePos(posq, posqCorrection, i, pos);
    }

    const mixed halfDt = half*dt;
    for (int k = blockIdx.x*blockDim.x+threadIdx.x; k < numBodies; k += blockDim.x*gridDim.x) {
        BodyData &body = bodyData[k];

#if COMPMOD == 1
        body.rdot = body.r*(-2.5) + (body.F*(body.invm*halfDt) - body.v)*halfDt;
        body.qdot = virtualRotation( body, -dt )*0.5 - body.q*(1.0/3.0);
#endif

        // Half-step integration of velocities
        body.v += body.F*(body.invm*halfDt);
        body.pi += body.Ctau*dt;

        // Full-step translation and rotation
        body.r += body.v*dt;
        ROTATION(body.q, body.pi, body.invI, dt);

        // Update of atomic positions and their displacements from the center of mass
        int loc = body.loc;
        for (int j = 0; j < body.N; j++) {
            int i = atomLocation[numFree + loc];
            mixed3 delta = Ct(body.q, B(body.q, bodyFixedPos[loc]));
            storePos(posq, posqCorrection, i, body.r + delta);
            posDelta[i] = fuse(delta, zero);
            loc++;
        }
    }
}

/*--------------------------------------------------------------------------------------------------
  Perform the final step of Rigid Body integration.
--------------------------------------------------------------------------------------------------*/

extern "C" __global__ void integrateRigidBodyPart3(int numAtoms,
                                                   int stride,
                                                   int numFree,
                                                   int numBodies,
                                                   const mixed dt,
                                                   real4* __restrict__ posq,
                                                   real4* __restrict__ posqCorrection,
                                                   mixed4* __restrict__ velm,
                                                   const long long* __restrict__ force,
                                                   mixed4* __restrict__ posDelta,
                                                   BodyData* __restrict__ bodyData,
                                                   const int* __restrict__ atomLocation,
                                                   const mixed3* __restrict__ bodyFixedPos,
                                                   mixed3* __restrict__ savedPos) {

    const mixed scale = one/(mixed)0x100000000;
    const mixed halfDt = half*dt;
    const mixed invDt = one/dt;

    for (int j = blockIdx.x*blockDim.x+threadIdx.x; j < numFree; j += blockDim.x*gridDim.x) {
        int i = atomLocation[j];
        mixed4& velocity = velm[i];
        mixed3 f = make_mixed3(force[i], force[i+stride], force[i+stride*2])*scale;
        mixed3 r = loadPos(posq, posqCorrection, i);
        mixed3 v = trim(velocity) + f*(velocity.w*halfDt) + (r - savedPos[j])*invDt;
        velocity = fuse(v, velocity.w);
    }

    for (int k = blockIdx.x*blockDim.x+threadIdx.x; k < numBodies; k += blockDim.x*gridDim.x) {
        BodyData &body = bodyData[k];

        // Computation of resultant force and torque
        body.F = make_mixed3(zero, zero, zero);
        mixed3 tau = make_mixed3(zero, zero, zero);
        int loc = numFree + body.loc;
        for (int j = 0; j < body.N; j++) {
            int i = atomLocation[loc++];
            mixed3 delta = trim(posDelta[i]);
            mixed3 f = make_mixed3(force[i], force[i+stride], force[i+stride*2])*scale;
            body.F += f;
            tau += cross(delta, f);
        }
        body.Ctau = C(body.q, tau);

        // Half-step integration of velocities
        body.v += body.F*(body.invm*halfDt);
        body.pi += body.Ctau*dt;

#if COMPMOD == 1
        body.rdot += body.r*2.5 + (body.v + body.F*(body.invm*halfDt))*halfDt;
        body.qdot += body.q*1.5 + virtualRotation( body, dt );
        body.qdot -= body.q*dot(body.qdot, body.q);
#endif

        // Update of atomic velocities
        mixed3 omega = (body.invI*Bt(body.q, body.pi))*half;
        mixed3 spaceFixedOmega = Ct(body.q, B(body.q, omega));
        loc = numFree + body.loc;
        for (int j = 0; j < body.N; j++) {
            int i = atomLocation[loc++];
            mixed3 delta = trim(posDelta[i]);
            velm[i] = fuse(body.v + cross(spaceFixedOmega, delta), velm[i].w);
        }
    }
}

/*--------------------------------------------------------------------------------------------------
  Computation of kinetic energy.
--------------------------------------------------------------------------------------------------*/

extern "C" __global__ void computeKineticEnergies(int numFree,
                                                  int numBodies,
                                                  mixed4* __restrict__ velm,
                                                  BodyData* __restrict__ bodyData,
                                                  const int* __restrict__ atomLocation,
                                                  mixed* __restrict__ atomKE,
                                                  mixed2* __restrict__ bodyKE) {

    for (int j = blockIdx.x*blockDim.x+threadIdx.x; j < numFree; j += blockDim.x*gridDim.x) {
        mixed4& v = velm[atomLocation[j]];
        atomKE[j] = half*(v.x*v.x + v.y*v.y + v.z*v.z)/v.w;
    }

    for (int k = blockIdx.x*blockDim.x+threadIdx.x; k < numBodies; k += blockDim.x*gridDim.x) {
        BodyData &body = bodyData[k];
        mixed3 L = Bt(body.q, body.pi)*half;
        bodyKE[k].x = dot(body.v/body.invm, body.v)*half;
        bodyKE[k].y = dot(L, L*body.invI)*half;
    }
}
