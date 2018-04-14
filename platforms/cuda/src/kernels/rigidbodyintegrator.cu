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
  Determine positions and velocities of body atoms
--------------------------------------------------------------------------------------------------*/

inline __device__ void positionsAndVelocities(int numFree,
                                              real4* __restrict__ posq,
                                              real4* __restrict__ posqCorrection,
                                              const int* __restrict__ atomLocation,
                                              const mixed3* __restrict__ bodyFixedPos,
                                              BodyData& body) {
    int loc = body.loc;
    for (int j = 0; j < body.N; j++) {
        int index = atomLocation[numFree + loc];
        mixed4 pos = loadPos(posq, posqCorrection, index);
        const mixed3 newPos = body.r + multiplyCt(body.q, multiplyB(body.q, bodyFixedPos[loc]));
        
    }
}


/*--------------------------------------------------------------------------------------------------
  Perform the initial step of Rigid Body integration.
--------------------------------------------------------------------------------------------------*/

extern "C" __global__ void integrateRigidBodyPart1(int numAtoms,
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
                                                   mixed4* __restrict__ savedPos) {

    const mixed scale = (mixed)1.0/(mixed)0x100000000;
    const mixed dtVel = (mixed)0.5*dt;
    const mixed dtVelScaled = scale*dtVel;

    for (int j = blockIdx.x*blockDim.x+threadIdx.x; j < numFree; j += blockDim.x*gridDim.x) {
        int index = atomLocation[j];
        mixed4& velocity = velm[index];
        if (velocity.w != (mixed)0.0) {
            velocity.x += dtVelScaled*force[index]*velocity.w;
            velocity.y += dtVelScaled*force[index+paddedNumAtoms]*velocity.w;
            velocity.z += dtVelScaled*force[index+paddedNumAtoms*2]*velocity.w;
            mixed4& delta = posDelta[index];
            delta.x = velocity.x*dt;
            delta.y = velocity.y*dt;
            delta.z = velocity.z*dt;
            mixed4 pos = loadPos(posq, posqCorrection, index);
            mixed4& saved = savedPos[j];
            saved.x = pos.x + delta.x;
            saved.y = pos.y + delta.y;
            saved.z = pos.z + delta.z;
        }
    }

//    for (int j = blockIdx.x*blockDim.x+threadIdx.x; j < numBodies; j += blockDim.x*gridDim.x) {
//        BodyData &body = bodyData[j];
//        body.p += body.F*dtVel/body.m;
//        body.r += body.p*dt;
//    }
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
                                                   mixed4* __restrict__ savedPos) {

    // Apply deltas to positions of free atoms:
    for (int j = blockIdx.x*blockDim.x+threadIdx.x; j < numFree; j += blockDim.x*gridDim.x) {
        int index = atomLocation[j];
        if (velm[index].w != (mixed)0.0) {
            mixed4 pos = loadPos(posq, posqCorrection, index);
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
//        printf("q[%d] = %f\n",k,body.q.x*body.q.x+body.q.y*body.q.y+body.q.z*body.q.z+body.q.w*body.q.w);
        for (int j = 0; j < body.N; j++) {
//            mixed3 pos = body.r + multiplyCt(body.q, multiplyB(body.q, bodyFixedPos[loc]));
            mixed3 pos = bodyFixedPos[loc];
//            printf("pos1[%d] = %f %f %f  loc = %d\n",j,pos.x,pos.y,pos.z,loc);
            int index = atomLocation[numFree + loc];
            mixed4 pos2 = loadPos(posq, posqCorrection, index);
//            printf("pos2[%d] = %f %f %f\n",j,pos2.x,pos2.y,pos2.z);
            loc++;
        }
    }
}

/*--------------------------------------------------------------------------------------------------
  Perform the final step of Rigid Body integration.
--------------------------------------------------------------------------------------------------*/

extern "C" __global__ void integrateRigidBodyPart3(int numAtoms,
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
                                                   mixed4* __restrict__ savedPos) {

    const mixed scale = (mixed)1.0/(mixed)0x100000000;
    const mixed dtVel = (mixed)0.5*dt;
    const mixed dtVelScaled = scale*dtVel;
    const mixed oneOverDt = (mixed)1.0/dt;

    for (int j = blockIdx.x*blockDim.x+threadIdx.x; j < numFree; j += blockDim.x*gridDim.x) {
        int index = atomLocation[j];
        mixed4& velocity = velm[index];
        if (velocity.w != (mixed)0.0) {
            mixed4 pos = loadPos(posq, posqCorrection, index);
            mixed4& saved = savedPos[j];
            velocity.x += dtVelScaled*force[index]*velocity.w + (pos.x - saved.x)*oneOverDt;
            velocity.y += dtVelScaled*force[index+paddedNumAtoms]*velocity.w + (pos.y - saved.y)*oneOverDt;
            velocity.z += dtVelScaled*force[index+paddedNumAtoms*2]*velocity.w + (pos.z - saved.z)*oneOverDt;
        }
    }

//    for (int k = blockIdx.x*blockDim.x+threadIdx.x; k < numBodies; k += blockDim.x*gridDim.x) {
//        BodyData &body = bodyData[k];
//        body.F = make_mixed3(0.0);
//        mixed3 tau = make_mixed3(0.0);
//        int loc = numFree + body.loc;
//        for (int j = 0; j < body.N; j++) {
//            int i = atomLocation[loc++];
//            const mixed3 delta = trimTo3(loadPos(posq, posqCorrection, i)) - body.r;
//            const mixed3 f = make_mixed3(scale*force[i],
//                                         scale*force[i+paddedNumAtoms],
//                                         scale*force[i+paddedNumAtoms*2]);
//            body.F += f;
//            tau += cross(delta, f);
//        }
//        body.Ctau = multiplyC(body.q, tau*2.0);
//        body.p += body.F*dtVel/body.m;
//    }
}
