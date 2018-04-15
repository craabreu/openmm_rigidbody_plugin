#define zero ((mixed)0.0)
#define half ((mixed)0.5)
#define one  ((mixed)1.0)

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
    mixed4 q;     // orientation quaternion
    mixed4 pi;    // quaternion-conjugated momentum
    mixed4 Ctau;  // quaternion-frame resultant torque
} BodyData;

/*--------------------------------------------------------------------------------------------------
  Retrive the size in bytes of a BodyData structure.
--------------------------------------------------------------------------------------------------*/

extern "C" __global__ void getBodyDataSize(size_t* __restrict__ size) {
    size[0] = sizeof(BodyData);
}

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
  Multiplications with rotation matrix factors
--------------------------------------------------------------------------------------------------*/

inline __device__ mixed4 multiplyB(mixed4 q, mixed3 v) {
    return make_mixed4(-q.y*v.x - q.z*v.y - q.w*v.z,
                        q.x*v.x - q.w*v.y + q.z*v.z,
                        q.w*v.x + q.x*v.y - q.y*v.z,
                       -q.z*v.x + q.y*v.y + q.x*v.z);
}

inline __device__ mixed4 multiplyC(mixed4 q, mixed3 v) {
    return make_mixed4(-q.y*v.x - q.z*v.y - q.w*v.z,
                        q.x*v.x + q.w*v.y - q.z*v.z,
                       -q.w*v.x + q.x*v.y + q.y*v.z,
                        q.z*v.x - q.y*v.y + q.x*v.z);
}

inline __device__ mixed3 multiplyBt(mixed4 q, mixed4 v) {
    return make_mixed3(-q.y*v.x + q.x*v.y + q.w*v.z - q.z*v.w,
                       -q.z*v.x - q.w*v.y + q.x*v.z + q.y*v.w,
                       -q.w*v.x + q.z*v.y - q.y*v.z + q.x*v.w);
}

inline __device__ mixed3 multiplyCt(mixed4 q, mixed4 v) {
    return make_mixed3(-q.y*v.x + q.x*v.y - q.w*v.z + q.z*v.w,
                       -q.z*v.x + q.w*v.y + q.x*v.z - q.y*v.w,
                       -q.w*v.x - q.z*v.y + q.y*v.z + q.x*v.w);
}

/*--------------------------------------------------------------------------------------------------
  Uniaxial rotation
--------------------------------------------------------------------------------------------------*/

inline __device__ void uniaxialRotation(BodyData& body, mixed dtBy4, int axis) {
    mixed4& q = body.q;
    mixed4& p = body.pi;
    mixed4 Bq, Bp;
    mixed invMoI;
    if (axis == 0) {
        Bq = make_mixed4(-q.y,  q.x,  q.w, -q.z);
        Bp = make_mixed4(-p.y,  p.x,  p.w, -p.z);
        invMoI = body.invI.x;
    }
    else if (axis == 1) {
        Bq = make_mixed4(-q.z, -q.w,  q.x,  q.y);
        Bp = make_mixed4(-p.z, -p.w,  p.x,  p.y);
        invMoI = body.invI.y;
    }
    else {
        Bq = make_mixed4(-q.w,  q.z, -q.y,  q.x);
        Bp = make_mixed4(-p.w,  p.z, -p.y,  p.x);
        invMoI = body.invI.z;
    }
    mixed omegaDtBy2 = dot(p, Bq)*dtBy4*invMoI;
    mixed vsin = sin(omegaDtBy2);
    mixed vcos = cos(omegaDtBy2);
    q = q*vcos + Bq*vsin;
    p = p*vcos + Bp*vsin;
}

/*--------------------------------------------------------------------------------------------------
  NO_SQUISH rotation
--------------------------------------------------------------------------------------------------*/

inline __device__ void noSquishRotation(BodyData& body, mixed dt, int n) {
    mixed dtBy8n = dt/(8.0*n);
    mixed dtBy4n = 2.0*dtBy8n;
    bool axis3 = body.invI.z != zero;
    for (int i = 0; i < n; i++) {
        if (axis3) uniaxialRotation(body, dtBy8n, 3);
        uniaxialRotation(body, dtBy8n, 2);
        uniaxialRotation(body, dtBy4n, 1);
        uniaxialRotation(body, dtBy8n, 2);
        if (axis3) uniaxialRotation(body, dtBy8n, 3);
    }
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
        if (velocity.w != zero) {
            mixed3 f = make_mixed3(force[i], force[i+stride], force[i+stride*2])*scale;
            mixed3 v = trimTo3(velocity) + f*(velocity.w*halfDt);
            mixed3 delta = v*dt;
            velocity = growTo4(v, velocity.w);
            posDelta[i] = growTo4(delta, 0.0);
            savedPos[j] = loadPos(posq, posqCorrection, i) + delta;
        }
    }
}

/*--------------------------------------------------------------------------------------------------
  Perform the final step of Rigid Body integration.
--------------------------------------------------------------------------------------------------*/

extern "C" __global__ void integrateRigidBodyPart2(int numAtoms,
                                                   int paddedNumAtoms,
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
        if (velm[i].w != zero) {
            mixed3 pos = loadPos(posq, posqCorrection, i) + trimTo3(posDelta[i]);
            storePos(posq, posqCorrection, i, pos);
        }
    }

    const mixed halfDt = half*dt;
    for (int k = blockIdx.x*blockDim.x+threadIdx.x; k < numBodies; k += blockDim.x*gridDim.x) {
        BodyData &body = bodyData[k];

        // Half-step integration of velocities
        body.v += body.F*(body.invm*halfDt);
        body.pi += body.Ctau*dt;

        // Full-step translation and rotation
        body.r += body.v*dt;
        noSquishRotation(body, dt, 1);

        // Update of atomic positions and deltas from center of mass
        int loc = body.loc;
        for (int j = 0; j < body.N; j++) {
            int i = atomLocation[numFree + loc];
            mixed3 delta = multiplyCt(body.q, multiplyB(body.q, bodyFixedPos[loc]));
            storePos(posq, posqCorrection, i, body.r + delta);
            posDelta[i] = growTo4(delta, 0.0);
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
        if (velocity.w != zero) {
            mixed3 f = make_mixed3(force[i], force[i+stride], force[i+stride*2])*scale;
            mixed3 v = trimTo3(velocity) + f*(velocity.w*halfDt);
            velocity = growTo4(v, velocity.w);
        }
    }

    for (int k = blockIdx.x*blockDim.x+threadIdx.x; k < numBodies; k += blockDim.x*gridDim.x) {
        BodyData &body = bodyData[k];

        // Computation of resultant force and torque
        body.F = make_mixed3(0.0);
        mixed3 tau = make_mixed3(0.0);
        int loc = numFree + body.loc;
        for (int j = 0; j < body.N; j++) {
            int i = atomLocation[loc++];
            mixed3 delta = trimTo3(posDelta[i]);
            mixed3 f = make_mixed3(force[i], force[i+stride], force[i+stride*2])*scale;
            body.F += f;
            tau += cross(delta, f);
        }
        body.Ctau = multiplyC(body.q, tau);

        // Half-step integration of velocities
        body.v += body.F*(body.invm*halfDt);
        body.pi += body.Ctau*dt;

        // Update of atomic velocities
        mixed3 omega = (body.invI*multiplyBt(body.q, body.pi))*half;
        mixed3 spaceFixedOmega = multiplyCt(body.q, multiplyB(body.q, omega));
        loc = numFree + body.loc;
        for (int j = 0; j < body.N; j++) {
            int i = atomLocation[loc++];
            mixed3 delta = trimTo3(posDelta[i]);
            velm[i] = growTo4(body.v + cross(spaceFixedOmega, delta), velm[i].w);
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
                                                  mixed* __restrict__ freeAtomKE,
                                                  mixed2* __restrict__ bodyKE) {

    for (int j = blockIdx.x*blockDim.x+threadIdx.x; j < numFree; j += blockDim.x*gridDim.x) {
        mixed4& v = velm[atomLocation[j]];
        freeAtomKE[j] = half*(v.x*v.x + v.y*v.y + v.z*v.z)/v.w;
    }

    for (int k = blockIdx.x*blockDim.x+threadIdx.x; k < numBodies; k += blockDim.x*gridDim.x) {
        BodyData &body = bodyData[k];
        mixed3& v = body.v;
        bodyKE[k].x = half*(v.x*v.x + v.y*v.y + v.z*v.z)/body.invm;
        mixed3 L = multiplyBt(body.q, body.pi)*half;
        mixed3 omega = body.invI*L;
        bodyKE[k].y = half*(L.x*omega.x + L.y*omega.y + L.z*omega.z);
    }
}
