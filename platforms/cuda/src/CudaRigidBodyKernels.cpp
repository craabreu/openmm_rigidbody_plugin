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
#include <vector>
#include <algorithm>

#include <iostream> // TEMPORARY

using namespace RigidBodyPlugin;
using namespace OpenMM;
using namespace std;

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

class CudaIntegrateRigidBodyStepKernel::ReorderListener : public CudaContext::ReorderListener {
public:
    ReorderListener(CudaContext& cu, const RigidBodySystem& bodySystem, CudaArray& atomLocation) :
        cu(cu), bodySystem(bodySystem), atomLocation(atomLocation) {
    }
    void execute() {
        int numActualAtoms = bodySystem.getNumActualAtoms();
        int stride = cu.getPaddedNumAtoms();
        vector<int> location(atomLocation.getSize());
        atomLocation.download(location);

        long long* force = (long long*) cu.getPinnedBuffer();
        cu.getForce().download(force);
        vector<long long> F(3*numActualAtoms);
        int k = 0;
        for (int i = 0; i < numActualAtoms; i++) {
            int loc = location[i];
            F[k++] = force[loc];
            F[k++] = force[loc+stride];
            F[k++] = force[loc+stride*2];
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
            force[loc+stride] = F[k++];
            force[loc+stride*2] = F[k++];
        }

        atomLocation.upload(location);
        cu.getForce().upload(force);
    }
private:
    CudaContext& cu;
    const RigidBodySystem& bodySystem;
    CudaArray& atomLocation;
};

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

template <class real, class real2>
vector<double> CudaIntegrateRigidBodyStepKernel::kineticEnergy(ContextImpl& context,
                                                               const RigidBodyIntegrator& integrator) {
    cu.setAsCurrent();
    int numAtoms = cu.getNumAtoms();

    // Call the first integration kernel.

    void* args[] = {&numFree, &numBodies,
                    &cu.getVelm().getDevicePointer(),
                    &bodyData.getDevicePointer(),
                    &atomLocation.getDevicePointer(),
                    &atomKE.getDevicePointer(),
                    &bodyKE.getDevicePointer()};

    cu.executeKernel(kineticEnergyKernel, args, numAtoms, 128);

    vector<real> KE(2, 0.0);
    if (numFree != 0) {
        real* aKE = (real*)pinnedBuffer;
        atomKE.download(aKE);
        for (int i = 0; i < numFree; i++)
            KE[0] += aKE[i];
    }
    if (numBodies != 0) {
        real2* bKE = (real2*)pinnedBuffer;
        bodyKE.download(bKE);
        for (int i = 0; i < numBodies; i++) {
            KE[0] += bKE[i].x;
            KE[1] += bKE[i].y;
        }
    }

    vector<double> dKE(KE.begin(), KE.end());
    return dKE;
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

template <class real, class real2, class real3>
real CudaIntegrateRigidBodyStepKernel::allocateArrays(size_t bodyDataSize) {
    if (numBodies != 0) {
        bodyData.initialize(cu, paddedNumBodies, bodyDataSize, "bodyData");
        bodyFixedPos.initialize<real3>(cu, paddedNumBodyAtoms, "bodyFixedPos");
        bodyKE.initialize<real2>(cu, paddedNumBodies, "bodyKE");
    }
    if (numFree != 0) {
        savedPos.initialize<real3>(cu, paddedNumFree, "savedPos");
        atomKE.initialize<real>(cu, paddedNumFree, "atomKE");
    }
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

void CudaIntegrateRigidBodyStepKernel::initialize(ContextImpl& context,
                                                  const RigidBodyIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(context.getSystem());
    cu.setAsCurrent();
    map<string, string> defines;

    int rotationMode = integrator.getRotationMode();
    defines["ROTATION"] = rotationMode == 0 ? "exactRotation" : "noSquishRotation";
    defines["NSPLIT"] = cu.intToString(rotationMode);
    CUmodule module = cu.createModule(CudaRigidBodyKernelSources::vectorOps +
                                      CudaRigidBodyKernelSources::elliptic +
                                      CudaRigidBodyKernelSources::rigidbodyintegrator,
                                      defines, "");
    kernel1 = cu.getKernel(module, "integrateRigidBodyPart1");
    kernel2 = cu.getKernel(module, "integrateRigidBodyPart2");
    kernel3 = cu.getKernel(module, "integrateRigidBodyPart3");
    kineticEnergyKernel = cu.getKernel(module, "computeKineticEnergies");

    bool mixedOrDouble = cu.getUseMixedPrecision() || cu.getUseDoublePrecision();
    size_t bodyDataSize = mixedOrDouble ? sizeof(bodyType<double,double3,double4>) : sizeof(bodyType<float,float3,float4>);

    const RigidBodySystem& bodySystem = integrator.getRigidBodySystem();
    numBodies = bodySystem.getNumBodies();
    numFree = bodySystem.getNumFree();
    int numActualAtoms = bodySystem.getNumActualAtoms();
    int numBodyAtoms = bodySystem.getNumBodyAtoms();

    // Compute padded numbers:
    int TileSize = cu.TileSize;
    paddedNumActualAtoms = TileSize*((numActualAtoms + TileSize - 1)/TileSize);
    paddedNumBodies = TileSize*((numBodies + TileSize - 1)/TileSize);
    paddedNumBodyAtoms = TileSize*((numBodyAtoms + TileSize - 1)/TileSize);
    paddedNumFree = TileSize*((numFree + TileSize - 1)/TileSize);

    // Create pinned buffer for fast memory transfer:
    size_t CoordinateSize = mixedOrDouble ? sizeof(double3) : sizeof(float3);
    size_t pinnedBufferSize = max(paddedNumBodies*bodyDataSize,
                              max(paddedNumActualAtoms*sizeof(int),
                              max(paddedNumBodyAtoms*CoordinateSize,
                                  paddedNumFree*CoordinateSize)));
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
        location[i] = invOrder[bodySystem.getAtomIndex(i)];
    atomLocation.upload(location);

    reorderListener = new ReorderListener(cu, bodySystem, atomLocation);
    cu.addReorderListener(reorderListener);

    // Allocate arrays
    if (mixedOrDouble)
        allocateArrays<double,double2,double3>(bodyDataSize);
    else
        allocateArrays<float,float2,float3>(bodyDataSize);
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

void CudaIntegrateRigidBodyStepKernel::uploadBodySystem(RigidBodySystem& bodySystem) {
    cu.setAsCurrent();

    int numBodies = bodySystem.getNumBodies();
    if (numBodies != 0) {
        int numBodyAtoms = bodySystem.getNumBodyAtoms();
        if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision()) {
            using bodyDouble = bodyType<double,double3,double4>;
            bodyDouble* data = (bodyDouble*) pinnedBuffer;
            for (int i = 0; i < numBodies; i++) {
                RigidBody b = bodySystem.getRigidBody(i);
                bodyDouble& body = data[i];
                body.N = b.N;
                body.loc = b.loc;
                body.invm = b.invMass;
                body.invI = make_double3(b.invI[0], b.invI[1], b.invI[2]);
                body.r = make_double3(b.rcm[0], b.rcm[1], b.rcm[2]);
                body.v = make_double3(b.pcm[0]/b.mass, b.pcm[1]/b.mass, b.pcm[2]/b.mass);
                body.F = make_double3(b.force[0], b.force[1], b.force[2]);
                body.q = make_double4(b.q[0], b.q[1], b.q[2], b.q[3]);
                body.pi = make_double4(b.pi[0], b.pi[1], b.pi[2], b.pi[3]);
                body.Ctau = make_double4(b.torque[0], b.torque[1], b.torque[2], b.torque[3]);
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
            using bodyFloat = bodyType<float,float3,float4>;
            bodyFloat* data = (bodyFloat*) pinnedBuffer;
            for (int i = 0; i < numBodies; i++) {
                RigidBody b = bodySystem.getRigidBody(i);
                bodyFloat& body = data[i];
                body.N = b.N;
                body.loc = b.loc;
                body.invm = b.invMass;
                body.invI = make_float3(b.invI[0], b.invI[1], b.invI[2]);
                body.r = make_float3(b.rcm[0], b.rcm[1], b.rcm[2]);
                body.v = make_float3(b.pcm[0]/b.mass, b.pcm[1]/b.mass, b.pcm[2]/b.mass);
                body.F = make_float3(b.force[0], b.force[1], b.force[2]);
                body.q = make_float4(b.q[0], b.q[1], b.q[2], b.q[3]);
                body.pi = make_float4(b.pi[0], b.pi[1], b.pi[2], b.pi[3]);
                body.Ctau = make_float4(b.torque[0], b.torque[1], b.torque[2], b.torque[3]);
            }
            bodyData.upload(data);

            float3* d = (float3*) pinnedBuffer;
            for (int i = 0; i < numBodyAtoms; i++) {
                Vec3 x = bodySystem.getBodyFixedPosition(i);
                d[i] = make_float3(x[0], x[1], x[2]);
            }
            bodyFixedPos.upload(d);
        }
    }
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

void CudaIntegrateRigidBodyStepKernel::execute(ContextImpl& context,
                                               const RigidBodyIntegrator& integrator) {
    context.updateContextState();
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    double timeStepDouble = integrator.getStepSize();
    float timeStepFloat = (float)timeStepDouble;
    bool useDouble = cu.getUseDoublePrecision() || cu.getUseMixedPrecision();
    int numThreads = max(numFree, numBodies);

    void* args[] = {&numAtoms, &paddedNumAtoms, &numFree, &numBodies,
                    useDouble ? (void*) &timeStepDouble : (void*) &timeStepFloat,
                    &cu.getPosq().getDevicePointer(),
                    &cu.getPosqCorrection().getDevicePointer(),
                    &cu.getVelm().getDevicePointer(),
                    &cu.getForce().getDevicePointer(),
                    &integration.getPosDelta().getDevicePointer(),
                    &bodyData.getDevicePointer(),
                    &atomLocation.getDevicePointer(),
                    &bodyFixedPos.getDevicePointer(),
                    &savedPos.getDevicePointer()};

    if (numFree != 0) {
        cu.executeKernel(kernel1, args, numFree, 128);
        integration.applyConstraints(integrator.getConstraintTolerance());
    }

    cu.executeKernel(kernel2, args, numThreads, 128);
    integration.computeVirtualSites();

    context.calcForcesAndEnergy(true, false);

    cu.executeKernel(kernel3, args, numThreads, 128);

    if (numFree != 0)
        integration.applyVelocityConstraints(integrator.getConstraintTolerance());

    cu.setTime(cu.getTime()+integrator.getStepSize());
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

double CudaIntegrateRigidBodyStepKernel::computeKineticEnergy(ContextImpl& context, const RigidBodyIntegrator& integrator)
{
    vector<double> KE;
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision())
        KE = kineticEnergy<double,double2>(context, integrator);
    else
        KE = kineticEnergy<float,float2>(context, integrator);
    return std::accumulate(KE.begin(), KE.end(), 0.0);
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

vector<double> CudaIntegrateRigidBodyStepKernel::getKineticEnergies(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision())
        return kineticEnergy<double,double2>(context, integrator);
    else
        return kineticEnergy<float,float2>(context, integrator);
}
