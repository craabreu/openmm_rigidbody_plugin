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

CudaIntegrateRigidBodyStepKernel::~CudaIntegrateRigidBodyStepKernel() {
    if (pinnedBuffer != NULL) {
        cuMemFreeHost(pinnedBuffer);
        delete atomLocation;
        delete bodyData;
        delete bodyFixedPos;
        delete savedPos;
        delete atomE;
        delete bodyE1;
        delete bodyE2;
        delete rdot;
    }
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

class CudaIntegrateRigidBodyStepKernel::ReorderListener : public CudaContext::ReorderListener {
public:
    ReorderListener(CudaContext& cu, const RigidBodySystem& bodySystem, CudaArray* atomLocation) :
        cu(cu), bodySystem(bodySystem), atomLocation(atomLocation) {
    }
    void execute() {
        int numActualAtoms = bodySystem.getNumActualAtoms();
        int stride = cu.getPaddedNumAtoms();
        vector<int> location(atomLocation->getSize());
        atomLocation->download(location);

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

        atomLocation->upload(location);
        cu.getForce().upload(force);
    }
private:
    CudaContext& cu;
    const RigidBodySystem& bodySystem;
    CudaArray* atomLocation;
};

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

template <class real, class real2>
vector<double> CudaIntegrateRigidBodyStepKernel::kineticEnergy(ContextImpl& context,
                                                               const RigidBodyIntegrator& integrator,
                                                               bool refined) {
    cu.setAsCurrent();
    int numAtoms = cu.getNumAtoms();

    // Call the first integration kernel.

    CUdeviceptr bodyDataPtr = numBodies != 0 ? bodyData->getDevicePointer() : 0;
    CUdeviceptr atomEPtr = numFree != 0 ? atomE->getDevicePointer() : 0;
    CUdeviceptr bodyE1Ptr = numBodies != 0 ? bodyE1->getDevicePointer() : 0;
    CUdeviceptr bodyE2Ptr = numBodies != 0 ? bodyE2->getDevicePointer() : 0;
    void* args[] = {&numFree, &numBodies,
                    &cu.getVelm().getDevicePointer(),
                    &bodyDataPtr, &atomLocation->getDevicePointer(),
                    &rdot->getDevicePointer(),
                    &atomEPtr, &bodyE1Ptr, &bodyE2Ptr};

    if (refined)
        cu.executeKernel(refinedKineticEnergyKernel, args, numAtoms, 128);
    else
        cu.executeKernel(kineticEnergyKernel, args, numAtoms, 128);

    real KEt = (real) 0.0;
    real KEr = (real) 0.0;
    real* E = (real*) pinnedBuffer;
    if (numFree != 0) {
        atomE->download(E);
        for (int i = 0; i < numFree; i++)
            KEt += E[i];
    }
    if (numBodies != 0) {
        bodyE1->download(E);
        for (int i = 0; i < numBodies; i++)
            KEt += E[i];
        bodyE2->download(E);
        for (int i = 0; i < numBodies; i++)
            KEr += E[i];
    }

    vector<double> KE(2);
    double factor = refined ? 1.0/(6.0*integrator.getStepSize()) : 1.0;
    KE[0] = KEt*factor;
    KE[1] = KEr*factor;
    return KE;
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

template <class real>
real CudaIntegrateRigidBodyStepKernel::potentialEnergyRefinement() {
    cu.setAsCurrent();
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    CUdeviceptr bodyDataPtr = numBodies != 0 ? bodyData->getDevicePointer() : 0;
    CUdeviceptr atomEPtr = numFree != 0 ? atomE->getDevicePointer() : 0;
    CUdeviceptr bodyE1Ptr = numBodies != 0 ? bodyE1->getDevicePointer() : 0;
    void* args[] = {&paddedNumAtoms, &numFree, &numBodies,
                    &cu.getVelm().getDevicePointer(),
                    &cu.getForce().getDevicePointer(),
                    &bodyDataPtr, &atomLocation->getDevicePointer(),
                    &atomEPtr, &bodyE1Ptr};
    cu.executeKernel(potentialEnergyRefinementKernel, args, max(numFree, numBodies), 128);
    real U = (real) 0.0;
    real* E = (real*) pinnedBuffer;
    if (numFree != 0) {
        atomE->download(E);
        for (int i = 0; i < numFree; i++)
            U += E[i];
    }
    if (numBodies != 0) {
        bodyE1->download(E);
        for (int i = 0; i < numBodies; i++)
            U += E[i];
    }
    return U;
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

template <class real, class real2, class real3, class real4>
size_t CudaIntegrateRigidBodyStepKernel::allocateArrays() {

    size_t sizeOne = sizeof(int);
    size_t maxSize = paddedNumActualAtoms*sizeOne;
    atomLocation = new CudaArray(cu, paddedNumActualAtoms, sizeOne, "atomLocation");

    if (numBodies != 0) {
        bodyE1 = new CudaArray(cu, paddedNumBodies, sizeof(real2), "bodyE1");
        bodyE2 = new CudaArray(cu, paddedNumBodies, sizeof(real2), "bodyE2");

        sizeOne = sizeof(bodyType<real,real3,real4>);
        maxSize = max(maxSize, paddedNumBodies*sizeOne);
        bodyData = new CudaArray(cu, paddedNumBodies, sizeOne, "bodyData");

        sizeOne = sizeof(real3);
        maxSize = max(maxSize, paddedNumBodyAtoms*sizeOne);
        bodyFixedPos = new CudaArray(cu, paddedNumBodyAtoms, sizeOne, "bodyFixedPos");
    }
    if (numFree != 0) {
        atomE = new CudaArray(cu, paddedNumFree, sizeof(real), "atomE");

        sizeOne = sizeof(real3);
        maxSize = max(maxSize, paddedNumFree*sizeOne);
        savedPos = new CudaArray(cu, paddedNumFree, sizeOne, "savedPos");
        rdot = new CudaArray(cu, paddedNumFree, sizeOne, "rdot");
    }
    return maxSize;
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

void CudaIntegrateRigidBodyStepKernel::initialize(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(context.getSystem());
    cu.setAsCurrent();
    map<string, string> defines;

    int rotationMode = integrator.getRotationMode();
    defines["ROTATION"] = rotationMode == 0 ? "exactRotation" : "noSquishRotation";
    defines["NSPLIT"] = cu.intToString(rotationMode);
    defines["COMPMOD"] = cu.intToString(integrator.getComputeRefinedEnergies() ? 1 : 0);
    CUmodule module = cu.createModule(CudaRigidBodyKernelSources::vectorOps +
                                      CudaRigidBodyKernelSources::elliptic +
                                      CudaRigidBodyKernelSources::rigidbodyintegrator,
                                      defines, "");
    freeAtomsDelta = cu.getKernel(module, "freeAtomsDelta");
    freeAtomsDot = cu.getKernel(module, "freeAtomsDot");
    rigidBodiesPart1 = cu.getKernel(module, "integrateRigidBodyPart1");
    rigidBodiesPart2 = cu.getKernel(module, "integrateRigidBodyPart2");
    kineticEnergyKernel = cu.getKernel(module, "kineticEnergies");
    refinedKineticEnergyKernel = cu.getKernel(module, "refinedKineticEnergies");
    potentialEnergyRefinementKernel = cu.getKernel(module, "potentialEnergyRefinement");

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

    // Allocate array and create pinned buffer for fast memory transfer
    size_t pinnedBufferSize;
    if (cu.getUseMixedPrecision() || cu.getUseDoublePrecision())
        pinnedBufferSize = allocateArrays<double,double2,double3,double4>();
    else
        pinnedBufferSize = allocateArrays<float,float2,float3,float4>();
    CUresult result = cuMemHostAlloc(&pinnedBuffer, pinnedBufferSize, 0);
    if (result != CUDA_SUCCESS)
        throw OpenMMException("Error creating pinned buffer for rigid body integrator");

    // Fill array of atom locations
    const vector<int>& order = cu.getAtomIndex();
    vector<int> invOrder(order.size());
    for (int i = 0; i < order.size(); i++)
        invOrder[order[i]] = i;
    int* location = (int*) pinnedBuffer;
    for (int i = 0; i < numActualAtoms; i++)
        location[i] = invOrder[bodySystem.getAtomIndex(i)];
    atomLocation->upload(location);

    reorderListener = new ReorderListener(cu, bodySystem, atomLocation);
    cu.addReorderListener(reorderListener);
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

void CudaIntegrateRigidBodyStepKernel::uploadBodySystem(RigidBodySystem& bodySystem) {
    cu.setAsCurrent();
    bool useDouble = cu.getUseDoublePrecision() || cu.getUseMixedPrecision();

    int numFree = bodySystem.getNumFree();
    if (numFree != 0) {
        if (useDouble) {
            double3* data = (double3*) pinnedBuffer;
            for (int i = 0; i < numFree; i++)
            data[i] = make_double3(0.0, 0.0, 0.0);
            rdot->upload(data);
        }
        else {
            float3* data = (float3*) pinnedBuffer;
            for (int i = 0; i < numFree; i++)
            data[i] = make_float3(0.0f, 0.0f, 0.0f);
            rdot->upload(data);
        }
    }

    int numBodies = bodySystem.getNumBodies();
    if (numBodies != 0) {
        int numBodyAtoms = bodySystem.getNumBodyAtoms();
        if (useDouble) {
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
                body.rdot = make_double3(0.0, 0.0, 0.0);
                body.q = make_double4(b.q[0], b.q[1], b.q[2], b.q[3]);
                body.pi = make_double4(b.pi[0], b.pi[1], b.pi[2], b.pi[3]);
                body.Ctau = make_double4(b.torque[0], b.torque[1], b.torque[2], b.torque[3]);
                body.qdot = make_double4(0.0, 0.0, 0.0, 0.0);
            }
            bodyData->upload(data);
            double3* d = (double3*) pinnedBuffer;
            for (int i = 0; i < numBodyAtoms; i++) {
                Vec3 x = bodySystem.getBodyFixedPosition(i);
                d[i] = make_double3(x[0], x[1], x[2]);
            }
            bodyFixedPos->upload(d);
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
                body.rdot = make_float3(0.0f, 0.0f, 0.0f);
                body.q = make_float4(b.q[0], b.q[1], b.q[2], b.q[3]);
                body.pi = make_float4(b.pi[0], b.pi[1], b.pi[2], b.pi[3]);
                body.Ctau = make_float4(b.torque[0], b.torque[1], b.torque[2], b.torque[3]);
                body.qdot = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
            }
            bodyData->upload(data);

            float3* d = (float3*) pinnedBuffer;
            for (int i = 0; i < numBodyAtoms; i++) {
                Vec3 x = bodySystem.getBodyFixedPosition(i);
                d[i] = make_float3(x[0], x[1], x[2]);
            }
            bodyFixedPos->upload(d);
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
    int numThreads = max(numFree, numBodies);
    int rdotRezero, rdotFactor;

    CUdeviceptr posCorrection = cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0;
    CUdeviceptr bodyDataPtr = numBodies != 0 ? bodyData->getDevicePointer() : 0;
    CUdeviceptr bodyFixedPosPtr = numBodies != 0 ? bodyFixedPos->getDevicePointer() : 0;
    CUdeviceptr savedPosPtr = numFree != 0 ? savedPos->getDevicePointer() : 0;

    void* args[] = {&paddedNumAtoms, &numFree, &numBodies,
                    useDouble ? (void*) &timeStepDouble : (void*) &timeStepFloat,
                    &cu.getPosq().getDevicePointer(),
                    &posCorrection, &cu.getVelm().getDevicePointer(),
                    &cu.getForce().getDevicePointer(),
                    &integration.getPosDelta().getDevicePointer(),
                    &bodyDataPtr, &atomLocation->getDevicePointer(),
                    &bodyFixedPosPtr, &savedPosPtr, &rdotRezero, &rdotFactor,
                    &rdot->getDevicePointer()};

    if (numFree != 0) {
        if (integrator.getComputeRefinedEnergies()) {
            timeStepDouble = -timeStepDouble;
            timeStepFloat = -timeStepFloat;
            integration.applyConstraints(integrator.getConstraintTolerance());
            cu.executeKernel(freeAtomsDelta, args, numFree, 128);
            rdotRezero = 1;
            rdotFactor = -1;
            cu.executeKernel(freeAtomsDot, args, numFree, 128);
            timeStepDouble = -timeStepDouble;
            timeStepFloat = -timeStepFloat;
        }
        cu.executeKernel(freeAtomsDelta, args, numFree, 128);
        integration.applyConstraints(integrator.getConstraintTolerance());
    }

    cu.executeKernel(rigidBodiesPart1, args, numThreads, 128);
    integration.computeVirtualSites();

    context.calcForcesAndEnergy(true, false);

    cu.executeKernel(rigidBodiesPart2, args, numThreads, 128);

    if (numFree != 0)  {
        integration.applyVelocityConstraints(integrator.getConstraintTolerance());
        if (integrator.getComputeRefinedEnergies()) {
            rdotRezero = 0;
            rdotFactor = 5;
            cu.executeKernel(freeAtomsDot, args, numFree, 128);
            cu.executeKernel(freeAtomsDelta, args, numFree, 128);
            integration.applyConstraints(integrator.getConstraintTolerance());
            rdotFactor = 2;
            cu.executeKernel(freeAtomsDot, args, numFree, 128);
        }
    }

    cu.setTime(cu.getTime()+integrator.getStepSize());
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

double CudaIntegrateRigidBodyStepKernel::computeKineticEnergy(ContextImpl& context, const RigidBodyIntegrator& integrator)
{
    cu.setAsCurrent();
    vector<double> KE;
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision())
        KE = kineticEnergy<double,double2>(context, integrator, false);
    else
        KE = kineticEnergy<float,float2>(context, integrator, false);
    return std::accumulate(KE.begin(), KE.end(), 0.0);
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

vector<double> CudaIntegrateRigidBodyStepKernel::getKineticEnergies(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    cu.setAsCurrent();
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision())
        return kineticEnergy<double,double2>(context, integrator, false);
    else
        return kineticEnergy<float,float2>(context, integrator, false);
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

vector<double> CudaIntegrateRigidBodyStepKernel::getRefinedKineticEnergies(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    cu.setAsCurrent();
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision())
        return kineticEnergy<double,double2>(context, integrator, true);
    else
        return kineticEnergy<float,float2>(context, integrator, true);
}

/*--------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------*/

double CudaIntegrateRigidBodyStepKernel::getPotentialEnergyRefinement(const RigidBodyIntegrator& integrator) {
    cu.setAsCurrent();
    double Uref;
    if (cu.getUseDoublePrecision() || cu.getUseMixedPrecision())
        Uref = potentialEnergyRefinement<double>();
    else
        Uref = potentialEnergyRefinement<float>();
    double dt = integrator.getStepSize();
    return -Uref*dt*dt/24.0;
}
