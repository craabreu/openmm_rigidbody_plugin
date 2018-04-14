/*--------------------------------------------------------------------------------------------------
  Define struct for rigid body data.
--------------------------------------------------------------------------------------------------*/

extern "C" typedef struct {
    int    N;     // number of atoms
    int    loc;   // pointer to set of atoms
    mixed3 I;     // principal moments of inertia
    mixed  m;     // mass
    mixed3 r;     // center-of-mass position
    mixed3 F;     // resultant force
    mixed3 p;     // center-of-mass momentum
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

    const mixed scale = (mixed)1.0/(mixed)0x100000000;
    const mixed dtVel = (mixed)0.5*dt;
    const mixed dtVelScaled = scale*dtVel;
    for (int j = blockIdx.x*blockDim.x+threadIdx.x; j < numFree; j += blockDim.x*gridDim.x) {
        int i = atomLocation[j];
        mixed4& velocity = velm[i];
        if (velocity.w != (mixed)0.0) {
            velocity.x += dtVelScaled*force[i]*velocity.w;
            velocity.y += dtVelScaled*force[i+stride]*velocity.w;
            velocity.z += dtVelScaled*force[i+stride*2]*velocity.w;
            mixed4& delta = posDelta[i];
            delta.x = velocity.x*dt;
            delta.y = velocity.y*dt;
            delta.z = velocity.z*dt;
            mixed3 pos = loadPos(posq, posqCorrection, i);
            savedPos[j] = pos + trimTo3(delta);
        }
    }

    for (int j = blockIdx.x*blockDim.x+threadIdx.x; j < numBodies; j += blockDim.x*gridDim.x) {
        BodyData &body = bodyData[j];
        body.p += body.F*dtVel/body.m;
        body.r += body.p*dt;
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
        int index = atomLocation[j];
        if (velm[index].w != (mixed)0.0) {
            mixed3 pos = loadPos(posq, posqCorrection, index);
            mixed4& delta = posDelta[index];
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
            storePos(posq, posqCorrection, index, pos);
        }
    }

    for (int k = blockIdx.x*blockDim.x+threadIdx.x; k < numBodies; k += blockDim.x*gridDim.x) {
        BodyData &body = bodyData[k];
        int loc = body.loc;
        for (int j = 0; j < body.N; j++) {
            mixed3 pos = body.r + multiplyCt(body.q, multiplyB(body.q, bodyFixedPos[loc]));
//            storePos(posq, posqCorrection, atomLocation[numFree + loc], pos);
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

    const mixed scale = (mixed)1.0/(mixed)0x100000000;
    const mixed dtVel = (mixed)0.5*dt;
    const mixed dtVelScaled = scale*dtVel;
    const mixed invDt = (mixed)1.0/dt;

    for (int j = blockIdx.x*blockDim.x+threadIdx.x; j < numFree; j += blockDim.x*gridDim.x) {
        int i = atomLocation[j];
        mixed4& velocity = velm[i];
        if (velocity.w != (mixed)0.0) {
            mixed3 pos = loadPos(posq, posqCorrection, i);
            mixed3& saved = savedPos[j];
            velocity.x += dtVelScaled*force[i]*velocity.w + (pos.x - saved.x)*invDt;
            velocity.y += dtVelScaled*force[i+stride]*velocity.w + (pos.y - saved.y)*invDt;
            velocity.z += dtVelScaled*force[i+stride*2]*velocity.w + (pos.z - saved.z)*invDt;
        }
    }

    for (int k = blockIdx.x*blockDim.x+threadIdx.x; k < numBodies; k += blockDim.x*gridDim.x) {
        BodyData &body = bodyData[k];
        body.F = make_mixed3(0.0);
        mixed3 tau = make_mixed3(0.0);
        int loc = numFree + body.loc;
        for (int j = 0; j < body.N; j++) {
            int i = atomLocation[loc++];
            const mixed3 delta = loadPos(posq, posqCorrection, i) - body.r;
            const mixed3 f = make_mixed3(force[i], force[i+stride], force[i+stride*2])*scale;
            body.F += f;
            tau += cross(delta, f);
        }
        body.Ctau = multiplyC(body.q, tau*2.0);
        body.p += body.F*(dtVel/body.m);
    }
}
