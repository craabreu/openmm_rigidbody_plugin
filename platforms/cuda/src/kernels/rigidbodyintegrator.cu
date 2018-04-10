/**
 * Define struct for rigid body data.
 */

extern "C" typedef struct {
    int    N;    // number of atoms
    mixed  m;    // mass
    mixed3 r;    // center-of-mass position
    mixed3 p;    // center-of-mass momentum
    mixed3 I;    // principal moments of inertia
    mixed4 q;    // orientation quaternion
    mixed4 pi;   // quaternion-conjugated momentum
    int    loc;  // pointer to set of atoms
} bodyData;


/**
 * Retrive the size in bytes of a bodyData structure.
 */
extern "C" __global__ void getBodyDataSize(size_t* __restrict__ size) {
    size[0] = sizeof(bodyData);
}

/**
 * Load the position of a particle.
 */
inline __device__ mixed4 loadPos(const real4* __restrict__ posq, const real4* __restrict__ posqCorrection, int index) {
#ifdef USE_MIXED_PRECISION
    real4 pos1 = posq[index];
    real4 pos2 = posqCorrection[index];
    return make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
    return posq[index];
#endif
}

/**
 * Store the position of a particle.
 */
inline __device__ void storePos(real4* __restrict__ posq, real4* __restrict__ posqCorrection, int index, mixed4 pos) {
#ifdef USE_MIXED_PRECISION
    posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
    posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
    posq[index] = pos;
#endif
}

/**
 * Compute resultant force and resultant torque.
 */
inline __device__ void forceAndTorque(int paddedNumAtoms, int numFree, const real4* __restrict__ posq,
                const real4* __restrict__ posqCorrection, const long long* __restrict__ force,
                const int* __restrict__ atomIndex, bodyData& b, mixed3& F, mixed3& tau) {
    F = tau = make_mixed3((mixed)0.0, (mixed)0.0, (mixed)0.0);
    int loc = numFree + b.loc;
    for (int j = 0; j < b.N; j++) {
        int i = atomIndex[loc++];
        mixed4 pos = loadPos(posq, posqCorrection, i);
        mixed3 f = make_mixed3(force[i], force[i+paddedNumAtoms], force[i+paddedNumAtoms*2]);
        mixed3 delta = make_mixed3(pos.x - b.r.x, pos.y - b.r.y, pos.z - b.r.z);
        F.x += f.x;
        F.y += f.y;
        F.z += f.z;
        tau.x += delta.y*f.z - delta.z*f.y;
        tau.y += delta.z*f.x - delta.x*f.z;
        tau.z += delta.x*f.y - delta.y*f.x;
    }
}

/**
 * Perform the first step of RigidBody integration.
 */

extern "C" __global__ void integrateRigidBodyPart1(int numAtoms, int paddedNumAtoms, const mixed2* __restrict__ dt, const real4* __restrict__ posq,
        const real4* __restrict__ posqCorrection, mixed4* __restrict__ velm, const long long* __restrict__ force, mixed4* __restrict__ posDelta,
        int numBodies, int numFree, bodyData* __restrict__ body, const int* __restrict__ atomIndex, const mixed3* __restrict__ bodyFixedPos) {
    const mixed2 stepSize = dt[0];
    const mixed dtPos = stepSize.y;
    const mixed dtVel = 0.5f*(stepSize.x+stepSize.y);
    const mixed scale = dtVel/(mixed) 0x100000000;

    // Half-step integration of free-atom velocities, which are then synchronized with positions
    // for constraint application
    for (int j = blockIdx.x*blockDim.x+threadIdx.x; j < numFree; j += blockDim.x*gridDim.x) {
        int index = atomIndex[j];
        mixed4& velocity = velm[index];
        velocity.x += scale*force[index]*velocity.w;
        velocity.y += scale*force[index+paddedNumAtoms]*velocity.w;
        velocity.z += scale*force[index+paddedNumAtoms*2]*velocity.w;
        mixed4& delta = posDelta[index];
        delta.x = velocity.x*dtPos;
        delta.y = velocity.y*dtPos;
        delta.z = velocity.z*dtPos;
    }

    // Full-step integration of rigid body velocities, positions, and orientations
    mixed3 F, tau;
    for (int k = blockIdx.x*blockDim.x+threadIdx.x; k < numBodies; k += blockDim.x*gridDim.x) {
        bodyData& b = body[k];
        forceAndTorque(paddedNumAtoms, numFree, posq, posqCorrection, force, atomIndex, b, F, tau);
    }



/*
    for (int index = blockIdx.x*blockDim.x+threadIdx.x; index < numAtoms; index += blockDim.x*gridDim.x) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            velocity.x += scale*force[index]*velocity.w;
            velocity.y += scale*force[index+paddedNumAtoms]*velocity.w;
            velocity.z += scale*force[index+paddedNumAtoms*2]*velocity.w;
            pos.x = velocity.x*dtPos;
            pos.y = velocity.y*dtPos;
            pos.z = velocity.z*dtPos;
            posDelta[index] = pos;
            velm[index] = velocity;
        }
    }
*/
}

/**
 * Perform the second step of RigidBody integration.
 */

extern "C" __global__ void integrateRigidBodyPart2(int numAtoms, mixed2* __restrict__ dt, real4* __restrict__ posq,
        real4* __restrict__ posqCorrection, mixed4* __restrict__ velm, const mixed4* __restrict__ posDelta,
        int numBodies, int numFree, bodyData* __restrict__ body, const int* __restrict__ atomIndex, const mixed3* __restrict__ bodyFixedPos) {
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
        velocity = make_mixed4((mixed) (delta.x*oneOverDt), (mixed) (delta.y*oneOverDt), (mixed) (delta.z*oneOverDt), velocity.w);
#else
        velocity = make_mixed4((mixed) (delta.x*oneOverDt+delta.x*correction), (mixed) (delta.y*oneOverDt+delta.y*correction), (mixed) (delta.z*oneOverDt+delta.z*correction), velocity.w);
#endif
        storePos(posq, posqCorrection, index, pos);
    }
}
