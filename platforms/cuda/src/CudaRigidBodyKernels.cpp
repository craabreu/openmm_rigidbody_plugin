/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "CudaRigidBodyKernels.h"
#include "CudaRigidBodyKernelSources.h"
#include "RigidBodyIntegrator.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/cuda/CudaPlatform.h"
#include "openmm/cuda/CudaBondedUtilities.h"
#include "openmm/cuda/CudaForceInfo.h"
#include "openmm/cuda/CudaIntegrationUtilities.h"

#include <iostream> // TEMPORARY

using namespace RigidBodyPlugin;
using namespace OpenMM;
using namespace std;

typedef struct {
    int    N;     // number of atoms
    int    loc;   // pointer to set of atoms
    float  m;     // mass
    float3 I;     // principal moments of inertia
    float3 r;     // center-of-mass position
    float3 p;     // center-of-mass momentum
    float3 F;     // resultant force
    float4 q;     // orientation quaternion
    float4 pi;    // quaternion-conjugated momentum
    float4 Ctau2; // quaternion-frame resultant torque
} bodyDataFloat;

typedef struct {
    int     N;     // number of atoms
    int     loc;   // pointer to set of atoms
    double  m;     // mass
    double3 I;     // principal moments of inertia
    double3 r;     // center-of-mass position
    double3 p;     // center-of-mass momentum
    double3 F;     // resultant force
    double4 q;     // orientation quaternion
    double4 pi;    // quaternion-conjugated momentum
    double4 Ctau2; // quaternion-frame resultant torque
} bodyDataDouble;

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

size_t CudaIntegrateRigidBodyStepKernel::getBodyDataSize(CUmodule& module) {
    CUfunction kernel = cu.getKernel(module, "getBodyDataSize");
    CUdeviceptr pointer;
    void* arg[] = {&pointer};
    CUresult result = cuMemAlloc(&pointer, sizeof(size_t));
    if (result != CUDA_SUCCESS)
        throw OpenMMException("Error creating variable for retrieving rigid body data size");
    cu.executeKernel(kernel, arg, 1, 128);
    size_t bodyDataSize;
    result = cuMemcpyDtoH(&bodyDataSize, pointer, sizeof(size_t));
    if (result != CUDA_SUCCESS)
        throw OpenMMException("Error retrieving rigid body data size from device memory");
    return bodyDataSize;
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

class CudaIntegrateRigidBodyStepKernel::ReorderListener : public CudaContext::ReorderListener {
public:
    ReorderListener(CudaContext& cu, RigidBodySystem& bodySystem, CudaArray& atomLocation) :
        cu(cu), bodySystem(bodySystem), atomLocation(atomLocation) {
    }
    void execute() {
        int numActualAtoms = bodySystem.getNumActualAtoms();
        int paddedNumAtoms = cu.getPaddedNumAtoms();
        vector<int> location(atomLocation.getSize());
        atomLocation.download(location);

        long long* force = (long long*) cu.getPinnedBuffer();
        cu.getForce().download(force);
        vector<long long> F(3*numActualAtoms);
        int k = 0;
        for (int i = 0; i < numActualAtoms; i++) {
            int loc = location[i];
            F[k++] = force[loc];
            F[k++] = force[loc+paddedNumAtoms];
            F[k++] = force[loc+paddedNumAtoms*2];
        }

        const vector<int>& order = cu.getAtomIndex();
        vector<int> invOrder(order.size());
        for (int i = 0; i < order.size(); i++)
            invOrder[order[i]] = i;
        for (int i = 0; i < bodySystem.getNumActualAtoms(); i++)
            location[i] = invOrder[bodySystem.getAtomIndex(i)];

        k = 0;
        for (int i = 0; i < numActualAtoms; i++) {
            int loc = location[i];
            force[loc] = F[k++];
            force[loc+paddedNumAtoms] = F[k++];
            force[loc+paddedNumAtoms*2] = F[k++];
        }

        atomLocation.upload(location);
        cu.getForce().upload(force);
    }
private:
    CudaContext& cu;
    RigidBodySystem& bodySystem;
    CudaArray& atomLocation;
};

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

void CudaIntegrateRigidBodyStepKernel::initialize(const System& system,
                                                  const RigidBodyIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(system);
    cu.setAsCurrent();
    map<string, string> defines;
    CUmodule module = cu.createModule(CudaRigidBodyKernelSources::vectorOps +
                                      CudaRigidBodyKernelSources::rigidbodyintegrator,
                                      defines, "");
    kernel1 = cu.getKernel(module, "integrateRigidBodyPart1");
    kernel2 = cu.getKernel(module, "integrateRigidBodyPart2");
    kernel3 = cu.getKernel(module, "integrateRigidBodyPart3");

    size_t bodyDataSize = getBodyDataSize(module);
    bool mixedOrDouble = cu.getUseMixedPrecision() || cu.getUseDoublePrecision();
    if (mixedOrDouble && bodyDataSize != sizeof(bodyDataDouble))
        throw OpenMMException("Error in size of memory chunk allocated for body data");
    else if (!mixedOrDouble && bodyDataSize != sizeof(bodyDataFloat))
        throw OpenMMException("Error in size of memory chunk allocated for body data");

    RigidBodySystem* bodySystem = integrator.getRigidBodySystem();
    numBodies = bodySystem->getNumBodies();
    numFree = bodySystem->getNumFree();
    int numActualAtoms = bodySystem->getNumActualAtoms();
    int numBodyAtoms = bodySystem->getNumBodyAtoms();

    // Compute padded numbers:
    int TileSize = cu.TileSize;
    paddedNumActualAtoms = TileSize*((numActualAtoms + TileSize - 1)/TileSize);
    paddedNumBodies = TileSize*((numBodies + TileSize - 1)/TileSize);
    paddedNumBodyAtoms = TileSize*((numBodyAtoms + TileSize - 1)/TileSize);

    // Create pinned buffer for fast memory transfer:
    size_t bodyFixedPosSize = mixedOrDouble ? sizeof(double3) : sizeof(float3);
    size_t pinnedBufferSize = max(paddedNumBodies*bodyDataSize,
                              max(paddedNumActualAtoms*sizeof(int),
                                  paddedNumBodyAtoms*bodyFixedPosSize));
    CUresult result = cuMemHostAlloc(&pinnedBuffer, pinnedBufferSize, 0);
    if (result != CUDA_SUCCESS)
        throw OpenMMException("Error creating pinned buffer for rigid body integrator");

    // Allocate array of atom locations:
    atomLocation.initialize<int>(cu, paddedNumActualAtoms, "atomLocation");
    const vector<int>& order = cu.getAtomIndex();
    vector<int> invOrder(order.size());
    for (int i = 0; i < order.size(); i++)
        invOrder[order[i]] = i;
    int* location = (int*) pinnedBuffer;
    for (int i = 0; i < numActualAtoms; i++)
        location[i] = invOrder[bodySystem->getAtomIndex(i)];
    atomLocation.upload(location);

    reorderListener = new ReorderListener(cu, *bodySystem, atomLocation);
    cu.addReorderListener(reorderListener);

    // Allocate arrays of body data and body-fixed atom positions, if necessary:
    if (numBodies != 0) {
        bodyData.initialize(cu, paddedNumBodies, bodyDataSize, "bodyData");
        if (mixedOrDouble)
            bodyFixedPos.initialize<double3>(cu, paddedNumBodyAtoms, "bodyFixedPos");
        else
            bodyFixedPos.initialize<float3>(cu, paddedNumBodyAtoms, "bodyFixedPos");
    }

    // Allocate array for saving positions:
    int paddedNumFree = TileSize*((numFree + TileSize - 1)/TileSize);
    if (mixedOrDouble)
      savedPos.initialize<double4>(cu, paddedNumFree, "savedPos");
    else
      savedPos.initialize<float4>(cu, paddedNumFree, "savedPos");
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

void CudaIntegrateRigidBodyStepKernel::uploadBodySystem(RigidBodySystem& bodySystem) {
    cu.setAsCurrent();

    int numBodies = bodySystem.getNumBodies();
    if (numBodies != 0) {
        int numBodyAtoms = bodySystem.getNumBodyAtoms();
        if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
            bodyDataDouble* data = (bodyDataDouble*) pinnedBuffer;
            for (int i = 0; i < numBodies; i++) {
                RigidBody& b = *bodySystem.getRigidBody(i);
                bodyDataDouble& body = data[i];
                body.N = b.N;
                body.loc = b.loc;
                body.m = b.mass;
                body.I = make_double3(b.MoI[0], b.MoI[1], b.MoI[2]);
                body.r = make_double3(b.rcm[0], b.rcm[1], b.rcm[2]);
                body.p = make_double3(b.pcm[0], b.pcm[1], b.pcm[2]);
                body.F = make_double3(b.force[0], b.force[1], b.force[2]);
                body.q = make_double4(b.q[0], b.q[1], b.q[2], b.q[3]);
                body.pi = make_double4(b.pi[0], b.pi[1], b.pi[2], b.pi[3]);
                body.Ctau2 = make_double4(b.torque[0], b.torque[1], b.torque[2], b.torque[3]);
            }
            bodyData.upload(data);

            double3* d = (double3*) pinnedBuffer;
            for (int i = 0; i < numBodyAtoms; i++) {
                Vec3 x = bodySystem.getBodyFixedPosition(i);
                d[i] = make_double3(x[0], x[1], x[2]);
            }
            bodyFixedPos.upload(d);
        }
        else {
            bodyDataFloat* data = (bodyDataFloat*) pinnedBuffer;
            for (int i = 0; i < numBodies; i++) {
                RigidBody& b = *bodySystem.getRigidBody(i);
                bodyDataFloat& body = data[i];
                body.N = b.N;
                body.loc = b.loc;
                body.m = (float)b.mass;
                body.I = make_float3((float)b.MoI[0], (float)b.MoI[1], (float)b.MoI[2]);
                body.r = make_float3((float)b.rcm[0], (float)b.rcm[1], (float)b.rcm[2]);
                body.p = make_float3((float)b.pcm[0], (float)b.pcm[1], (float)b.pcm[2]);
                body.F = make_float3((float)b.force[0], (float)b.force[1], (float)b.force[2]);
                body.q = make_float4((float)b.q[0], (float)b.q[1], (float)b.q[2], (float)b.q[3]);
                body.pi = make_float4((float)b.pi[0], (float)b.pi[1], (float)b.pi[2], (float)b.pi[3]);
                body.Ctau2 = make_float4((float)b.torque[0], (float)b.torque[1], (float)b.torque[2], (float)b.torque[3]);
            }
            bodyData.upload(data);

            float3* d = (float3*) pinnedBuffer;
            for (int i = 0; i < numBodyAtoms; i++) {
                Vec3 x = bodySystem.getBodyFixedPosition(i);
                d[i] = make_float3((float)x[0], (float)x[1], (float)x[2]);
            }
            bodyFixedPos.upload(d);
        }
    }
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

void CudaIntegrateRigidBodyStepKernel::execute(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    context.updateContextState();
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    double timeStepDouble = integrator.getStepSize();
    float timeStepFloat = (float)timeStepDouble;
    bool useDouble = cu.getUseDoublePrecision() || cu.getUseMixedPrecision();

    // Call the first integration kernel.

    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    void* args[] = {&numAtoms, &paddedNumAtoms, &numFree, &numBodies,
                    useDouble ? (void*) &timeStepDouble : (void*) &timeStepFloat,
                    &cu.getPosq().getDevicePointer(),
                    &posCorrection,
                    &cu.getVelm().getDevicePointer(),
                    &cu.getForce().getDevicePointer(),
                    &integration.getPosDelta().getDevicePointer(),
                    &bodyData.getDevicePointer(),
                    &atomLocation.getDevicePointer(),
                    &bodyFixedPos.getDevicePointer(),
                    &savedPos.getDevicePointer()};

    cu.executeKernel(kernel1, args, numAtoms, 128);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());
    cu.executeKernel(kernel2, args, numAtoms, 128);
    integration.computeVirtualSites();

    context.calcForcesAndEnergy(true, false);

    cu.executeKernel(kernel3, args, numAtoms, 128);

    integration.applyVelocityConstraints(integrator.getConstraintTolerance());

    // Update the time and step count.

    cu.setTime(cu.getTime()+integrator.getStepSize());
    cu.setStepCount(cu.getStepCount()+1);
//    cu.reorderAtoms();
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

double CudaIntegrateRigidBodyStepKernel::computeKineticEnergy(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    return cu.getIntegrationUtilities().computeKineticEnergy(0.0);
}
