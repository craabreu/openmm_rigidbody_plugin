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


extern "C" __global__ void getBodyDataSize(size_t* __restrict__ size) {
    size[0] = sizeof(bodyData);
}


/**
 * Perform the first step of RigidBody integration.
 */

extern "C" __global__ void integrateRigidBodyPart1(int numAtoms, int paddedNumAtoms, const mixed2* __restrict__ dt, const real4* __restrict__ posq,
        const real4* __restrict__ posqCorrection, mixed4* __restrict__ velm, const long long* __restrict__ force, mixed4* __restrict__ posDelta) {
    const mixed2 stepSize = dt[0];
    const mixed dtPos = stepSize.y;
    const mixed dtVel = 0.5f*(stepSize.x+stepSize.y);
    const mixed scale = dtVel/(mixed) 0x100000000;
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
}

/**
 * Perform the second step of RigidBody integration.
 */

extern "C" __global__ void integrateRigidBodyPart2(int numAtoms, mixed2* __restrict__ dt, real4* __restrict__ posq,
        real4* __restrict__ posqCorrection, mixed4* __restrict__ velm, const mixed4* __restrict__ posDelta) {
    mixed2 stepSize = dt[0];
#if __CUDA_ARCH__ >= 130
    double oneOverDt = 1.0/stepSize.y;
#else
    float oneOverDt = 1.0f/stepSize.y;
    float correction = (1.0f-oneOverDt*stepSize.y)/stepSize.y;
#endif
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    if (index == 0)
        dt[0].x = stepSize.y;
    for (; index < numAtoms; index += blockDim.x*gridDim.x) {
        mixed4 velocity = velm[index];
        if (velocity.w != 0.0) {
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[index];
            real4 pos2 = posqCorrection[index];
            mixed4 pos = make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[index];
#endif
            mixed4 delta = posDelta[index];
            pos.x += delta.x;
            pos.y += delta.y;
            pos.z += delta.z;
#if __CUDA_ARCH__ >= 130
            velocity = make_mixed4((mixed) (delta.x*oneOverDt), (mixed) (delta.y*oneOverDt), (mixed) (delta.z*oneOverDt), velocity.w);
#else
            velocity = make_mixed4((mixed) (delta.x*oneOverDt+delta.x*correction), (mixed) (delta.y*oneOverDt+delta.y*correction), (mixed) (delta.z*oneOverDt+delta.z*correction), velocity.w);
#endif
#ifdef USE_MIXED_PRECISION
            posq[index] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
            posqCorrection[index] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[index] = pos;
#endif
            velm[index] = velocity;
        }
    }
}
