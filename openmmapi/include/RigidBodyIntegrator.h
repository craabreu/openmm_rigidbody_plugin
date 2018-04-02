#ifndef OPENMM_RIGIDBODYINTEGRATOR_H_
#define OPENMM_RIGIDBODYINTEGRATOR_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-2012 Stanford University and the Authors.      *
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

#include "RigidBodySystem.h"
#include "openmm/Context.h"
#include "openmm/Integrator.h"
#include "openmm/Kernel.h"
#include <vector>
#include "internal/windowsExportRigidBody.h"

namespace RigidBodyPlugin {

/**
 * This is an Integrator which simulates a System using the leap-frog RigidBody algorithm.
 */

class OPENMM_EXPORT_RIGIDBODY RigidBodyIntegrator : public OpenMM::Integrator {
public:
    /**
     * Create a RigidBodyIntegrator.
     * 
     * @param stepSize the step size with which to integrate the system (in picoseconds)
     */
    explicit RigidBodyIntegrator(double stepSize, const std::vector<int>& bodyIndices);
    /**
     * Advance a simulation through time by taking a series of time steps.
     * 
     * @param steps   the number of time steps to take
     */
    void step(int steps);
    /**
     * Retrieve a vector of integers containing the rigid-body indices.
     */
    std::vector<int> getBodyIndices() const;
protected:
    /**
     * This will be called by the Context when it is created.  It informs the Integrator
     * of what context it will be integrating, and gives it a chance to do any necessary initialization.
     * It will also get called again if the application calls reinitialize() on the Context.
     */
    void initialize(OpenMM::ContextImpl& context);
    /**
     * This will be called by the Context when it is destroyed to let the Integrator do any necessary
     * cleanup.  It will also get called again if the application calls reinitialize() on the Context.
     */
    void cleanup();
    /**
     * Get the names of all Kernels used by this Integrator.
     */
    std::vector<std::string> getKernelNames();
    /**
     * This will be called by the Context when the user modifies aspects of the context state, such
     * as positions, velocities, or parameters.
     * 
     * @param changed     this specifies what aspect of the Context was changed
     */
    void stateChanged(State::DataType changed);
    /**
     * Compute the kinetic energy of the system at the current time.
     */
    double computeKineticEnergy();
private:
    std::vector<int> bodyIndices;
    RigidBodySystem* bodySystem;
    OpenMM::Kernel kernel;
};

} // namespace OpenMM

#endif /*OPENMM_RIGIDBODYINTEGRATOR_H_*/
