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

void CudaIntegrateRigidBodyStepKernel::initialize(const System& system, const RigidBodyIntegrator& integrator) {
    cu.getPlatformData().initializeContexts(system);
    cu.setAsCurrent();
    map<string, string> defines;
    CUmodule module = cu.createModule(CudaRigidBodyKernelSources::rigidbodyintegrator, defines, "");
    kernel1 = cu.getKernel(module, "integrateRigidBodyPart1");
    kernel2 = cu.getKernel(module, "integrateRigidBodyPart2");

    RigidBodySystem* bodySystem = integrator.getRigidBodySystem();
    int numActualAtoms = bodySystem->getNumActualAtoms();
    int numBodies = bodySystem->getNumBodies();
    int numBodyAtoms = bodySystem->getNumBodyAtoms();

    // Compute padded numbers:
    int TileSize = cu.TileSize;
    paddedNumActualAtoms = TileSize*((numActualAtoms + TileSize - 1)/TileSize);
    paddedNumBodies = TileSize*((numBodies + TileSize - 1)/TileSize);
    paddedNumBodyAtoms = TileSize*((numBodyAtoms + TileSize - 1)/TileSize);

    // Create pinned buffer for fast memory transfer:
    size_t bodyDataSize = getBodyDataSize(module);
    size_t bodyFixedPosSize = cu.getUseDoublePrecision() ? sizeof(double3) : sizeof(float3);
    size_t pinnedBufferSize = max(paddedNumBodies*bodyDataSize, max(
                                  paddedNumActualAtoms*sizeof(int),
                                  paddedNumBodyAtoms*bodyFixedPosSize));
    CUresult result = cuMemHostAlloc(&pinnedBuffer, pinnedBufferSize, 0);
    if (result != CUDA_SUCCESS)
        throw OpenMMException("Error creating pinned buffer for rigid body integration");

    // Allocate and upload array of atom indices:
    atomIndex.initialize<int>(cu, paddedNumActualAtoms, "atomIndex");
    int* index = (int*) pinnedBuffer;
    for (int i = 0; i < numActualAtoms; ++i)
        index[i] = bodySystem->getAtomIndex(i);
    atomIndex.upload(index);

    // Allocate arrays of body data and body-fixed atom positions, if necessary:
    if (numBodies != 0) {
        bodyData.initialize(cu, paddedNumBodies, bodyDataSize, "bodyData");
        if (cu.getUseDoublePrecision())
            bodyFixedPos.initialize<double3>(cu, paddedNumBodyAtoms, "bodyFixedPos");
        else
            bodyFixedPos.initialize<float3>(cu, paddedNumBodyAtoms, "bodyFixedPos");
    }
}

void CudaIntegrateRigidBodyStepKernel::uploadBodySystem(const RigidBodyIntegrator& integrator) {
    cu.setAsCurrent();
    RigidBodySystem* bodySystem = integrator.getRigidBodySystem();
    int numBodies = bodySystem->getNumBodies();
    if (numBodies != 0) {
        int numBodyAtoms = bodySystem->getNumBodyAtoms();
        if (cu.getUseDoublePrecision()) {
            double3* d = (double3*) pinnedBuffer;
            for (int i = 0; i < numBodyAtoms; ++i) {
                Vec3 x = bodySystem->getBodyFixedPosition(i);
                d[i] = make_double3(x[0], x[1], x[2]);
            }
            bodyFixedPos.upload(d);
        }
        else {
            float3* d = (float3*) pinnedBuffer;
            for (int i = 0; i < numBodyAtoms; ++i) {
                Vec3 x = bodySystem->getBodyFixedPosition(i);
                d[i] = make_float3((float)x[0], (float)x[1], (float)x[2]);
            }
            bodyFixedPos.upload(d);
        }
    }
    cout<<"TEST\n";
}

void CudaIntegrateRigidBodyStepKernel::execute(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    cu.setAsCurrent();
    CudaIntegrationUtilities& integration = cu.getIntegrationUtilities();
    int numAtoms = cu.getNumAtoms();
    int paddedNumAtoms = cu.getPaddedNumAtoms();
    double dt = integrator.getStepSize();
    cu.getIntegrationUtilities().setNextStepSize(dt);

    // Call the first integration kernel.

    CUdeviceptr posCorrection = (cu.getUseMixedPrecision() ? cu.getPosqCorrection().getDevicePointer() : 0);
    void* args1[] = {&numAtoms, &paddedNumAtoms, &cu.getIntegrationUtilities().getStepSize().getDevicePointer(), &cu.getPosq().getDevicePointer(), &posCorrection,
            &cu.getVelm().getDevicePointer(), &cu.getForce().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel1, args1, numAtoms, 128);

    // Apply constraints.

    integration.applyConstraints(integrator.getConstraintTolerance());

    // Call the second integration kernel.

    void* args2[] = {&numAtoms, &cu.getIntegrationUtilities().getStepSize().getDevicePointer(), &cu.getPosq().getDevicePointer(), &posCorrection,
            &cu.getVelm().getDevicePointer(), &integration.getPosDelta().getDevicePointer()};
    cu.executeKernel(kernel2, args2, numAtoms, 128);
    integration.computeVirtualSites();

    // Update the time and step count.

    cu.setTime(cu.getTime()+dt);
    cu.setStepCount(cu.getStepCount()+1);
    cu.reorderAtoms();
}

double CudaIntegrateRigidBodyStepKernel::computeKineticEnergy(ContextImpl& context, const RigidBodyIntegrator& integrator) {
    return cu.getIntegrationUtilities().computeKineticEnergy(0.5*integrator.getStepSize());
}
