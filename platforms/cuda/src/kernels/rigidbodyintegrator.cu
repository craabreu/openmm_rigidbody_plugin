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
    mixed4 Ctau2; // quaternion-frame resultant torque
} bodyData;

/*--------------------------------------------------------------------------------------------------
  Retrive the size in bytes of a bodyData structure.
--------------------------------------------------------------------------------------------------*/

extern "C" __global__ void getBodyDataSize(size_t* __restrict__ size) {
    size[0] = sizeof(bodyData);
}

/*--------------------------------------------------------------------------------------------------
  Load the position of a particle.
--------------------------------------------------------------------------------------------------*/

inline __device__ mixed4 loadPos(const real4* __restrict__ posq,
                                 const real4* __restrict__ posqCorrection,
                                 int index) {
#ifdef USE_MIXED_PRECISION
    real4 pos1 = posq[index];
    real4 pos2 = posqCorrection[index];
    return make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
    return posq[index];
#endif
}

/*--------------------------------------------------------------------------------------------------
  Store the position of a particle.
--------------------------------------------------------------------------------------------------*/

inline __device__ void storePos(real4* __restrict__ posq,
                                real4* __restrict__ posqCorrection,
                                int index,
                                mixed4 pos) {
#ifdef USE_MIXED_PRECISION
    posq[index] = make_real4((real)pos.x, (real)pos.y, (real)pos.z, (real)pos.w);
    posqCorrection[index] = make_real4(pos.x-(real)pos.x, pos.y-(real)pos.y, pos.z-(real)pos.z, 0);
#else
    posq[index] = pos;
#endif
}

/*--------------------------------------------------------------------------------------------------
  Premultiply matrix C
--------------------------------------------------------------------------------------------------*/

inline __device__ mixed4 multiplyC(mixed4 q, mixed3 v) {
    return make_mixed4(-q.y*v.x - q.z*v.y - q.w*v.z,
                        q.x*v.x + q.w*v.y - q.z*v.z,
                       -q.w*v.x + q.x*v.y + q.y*v.z,
                        q.z*v.x - q.y*v.y + q.x*v.z);
}

/*--------------------------------------------------------------------------------------------------
  Compute resultant force and resultant torque.
--------------------------------------------------------------------------------------------------*/

inline __device__ void forceAndTorque(int paddedNumAtoms,
                                      int numFree, const real4* __restrict__ posq,
                                      const real4* __restrict__ posqCorrection,
                                      const long long* __restrict__ force,
                                      const int* __restrict__ atomIndex,
                                      const mixed scale,
                                      bodyData& b) {
    mixed3 tau;
    b.F = tau = make_mixed3((mixed)0.0, (mixed)0.0, (mixed)0.0);
    int loc = numFree + b.loc;
    for (int j = 0; j < b.N; j++) {
        int i = atomIndex[loc++];
        mixed4 pos = loadPos(posq, posqCorrection, i);
        mixed3 f = make_mixed3(force[i], force[i+paddedNumAtoms], force[i+paddedNumAtoms*2]);
        mixed3 delta = make_mixed3(pos.x - b.r.x, pos.y - b.r.y, pos.z - b.r.z);
        b.F.x += f.x;
        b.F.y += f.y;
        b.F.z += f.z;
        tau.x += delta.y*f.z - delta.z*f.y;
        tau.y += delta.z*f.x - delta.x*f.z;
        tau.z += delta.x*f.y - delta.y*f.x;
    }
    b.F.x *= scale;
    b.F.y *= scale;
    b.F.z *= scale;
    mixed scale2 = scale*(mixed)2.0;
    tau.x *= scale2;
    tau.y *= scale2;
    tau.z *= scale2;
    b.Ctau2 = multiplyC(b.q, tau);
}

/*--------------------------------------------------------------------------------------------------
  Perform the initial step of Rigid Body integration.
--------------------------------------------------------------------------------------------------*/

extern "C" __global__ void integrateRigidBodyPart1(int numAtoms,
                                                   int paddedNumAtoms,
                                                   const mixed2* __restrict__ dt,
                                                   const real4* __restrict__ posq,
                                                   const real4* __restrict__ posqCorrection,
                                                   mixed4* __restrict__ velm,
                                                   const long long* __restrict__ force,
                                                   mixed4* __restrict__ posDelta,
                                                   int numBodies,
                                                   int numFree,
                                                   bodyData* __restrict__ body,
                                                   const int* __restrict__ atomIndex,
                                                   const mixed3* __restrict__ bodyFixedPos) {
    const mixed2 stepSize = dt[0];
    const mixed dtPos = stepSize.y;
    const mixed dtVel = 0.5f*(stepSize.x+stepSize.y);
    const mixed scale = (mixed) 1.0/(mixed) 0x100000000;
    const mixed dtVelScaled = scale*dtVel;

    // Half-step integration of free-atom velocities, which are then synchronized with positions
    // for constraint application
    for (int j = blockIdx.x*blockDim.x+threadIdx.x; j < numFree; j += blockDim.x*gridDim.x) {
        int index = atomIndex[j];
        mixed4& velocity = velm[index];
        velocity.x += dtVelScaled*force[index]*velocity.w;
        velocity.y += dtVelScaled*force[index+paddedNumAtoms]*velocity.w;
        velocity.z += dtVelScaled*force[index+paddedNumAtoms*2]*velocity.w;
        mixed4& delta = posDelta[index];
        delta.x = velocity.x*dtPos;
        delta.y = velocity.y*dtPos;
        delta.z = velocity.z*dtPos;
    }

    // Full-step integration of rigid body velocities, positions, and orientations
    for (int k = blockIdx.x*blockDim.x+threadIdx.x; k < numBodies; k += blockDim.x*gridDim.x) {
        bodyData& b = body[k];
        forceAndTorque(paddedNumAtoms, numFree, posq, posqCorrection, force, atomIndex, scale, b);
    }
}

/*--------------------------------------------------------------------------------------------------
  Perform the final step of Rigid Body integration.
--------------------------------------------------------------------------------------------------*/

extern "C" __global__ void integrateRigidBodyPart2(int numAtoms,
                                                   mixed2* __restrict__ dt,
                                                   real4* __restrict__ posq,
                                                   real4* __restrict__ posqCorrection,
                                                   mixed4* __restrict__ velm,
                                                   const mixed4* __restrict__ posDelta,
                                                   int numBodies,
                                                   int numFree,
                                                   bodyData* __restrict__ body,
                                                   const int* __restrict__ atomIndex,
                                                   const mixed3* __restrict__ bodyFixedPos) {
    mixed2 stepSize = dt[0];
#if __CUDA_ARCH__ >= 130
    double oneOverDt = 1.0/stepSize.y;
#else
    float oneOverDt = 1.0f/stepSize.y;
    float correction = (1.0f-oneOverDt*stepSize.y)/stepSize.y;
#endif
    int j = blockIdx.x*blockDim.x+threadIdx.x;
    if (j == 0)
        dt[0].x = stepSize.y;
    for (; j < numFree; j += blockDim.x*gridDim.x) {
        int index = atomIndex[j];
        mixed4 pos = loadPos(posq, posqCorrection, index);;
        const mixed4& delta = posDelta[index];
        mixed4& velocity = velm[index];
        pos.x += delta.x;
        pos.y += delta.y;
        pos.z += delta.z;
#if __CUDA_ARCH__ >= 130
        velocity = make_mixed4((mixed) (delta.x*oneOverDt),
                               (mixed) (delta.y*oneOverDt),
                               (mixed) (delta.z*oneOverDt), velocity.w);
#else
        velocity = make_mixed4((mixed) (delta.x*oneOverDt+delta.x*correction),
                               (mixed) (delta.y*oneOverDt+delta.y*correction),
                               (mixed) (delta.z*oneOverDt+delta.z*correction), velocity.w);
#endif
        storePos(posq, posqCorrection, index, pos);
    }
}
