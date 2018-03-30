#ifndef OPENCL_RIGIDBODY_KERNELS_H_
#define OPENCL_RIGIDBODY_KERNELS_H_

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

#include "RigidBodyKernels.h"
#include "openmm/opencl/OpenCLContext.h"
#include "openmm/opencl/OpenCLArray.h"

namespace RigidBodyPlugin {

/**
 * This kernel is invoked by RigidBodyIntegrator to take one time step.
 */
class OpenCLIntegrateRigidBodyStepKernel : public IntegrateRigidBodyStepKernel {
public:
    OpenCLIntegrateRigidBodyStepKernel(std::string name, const OpenMM::Platform& platform, OpenMM::OpenCLContext& cl) : IntegrateRigidBodyStepKernel(name, platform), cl(cl),
            hasInitializedKernels(false) {
    }
    ~OpenCLIntegrateRigidBodyStepKernel();
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param integrator the RigidBodyIntegrator this kernel will be used for
     */
    void initialize(const OpenMM::System& system, const RigidBodyIntegrator& integrator);
    /**
     * Execute the kernel.
     *
     * @param context    the context in which to execute this kernel
     * @param integrator the RigidBodyIntegrator this kernel is being used for
     */
    void execute(OpenMM::ContextImpl& context, const RigidBodyIntegrator& integrator);
    /**
     * Compute the kinetic energy.
     * 
     * @param context    the context in which to execute this kernel
     * @param integrator the RigidBodyIntegrator this kernel is being used for
     */
    double computeKineticEnergy(OpenMM::ContextImpl& context, const RigidBodyIntegrator& integrator);
private:
    OpenMM::OpenCLContext& cl;
    bool hasInitializedKernels;
    cl::Kernel kernel1, kernel2;
};

} // namespace RigidBodyPlugin

#endif /*OPENCL_RIGIDBODY_KERNELS_H_*/
