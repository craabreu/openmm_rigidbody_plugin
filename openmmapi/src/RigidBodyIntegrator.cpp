/* -------------------------------------------------------------------------- *
 *                          OpenMM Rigid Body Plugin                          *
 * -------------------------------------------------------------------------- */

#include "RigidBodyIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include "RigidBodyKernels.h"

using namespace RigidBodyPlugin;
using namespace OpenMM;
using std::string;
using std::vector;

RigidBodyIntegrator::RigidBodyIntegrator(double stepSize, const vector<int>& bodyIndices) {
    setStepSize(stepSize);
    setConstraintTolerance(1e-5);
    this->bodyIndices = bodyIndices;
}

void RigidBodyIntegrator::setRotationMode(int mode) {
    if (mode < 0)
        throw OpenMMException("Rotation mode cannot be negative");
    if (owner != NULL)
        throw OpenMMException("Cannot set rotation mode if integrator is bound to a context");
    rotationMode = mode;
}

void RigidBodyIntegrator::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    context = &contextRef;
    owner = &contextRef.getOwner();
    const System *system = &contextRef.getSystem();
    if (system->getNumParticles() != bodyIndices.size())
        throw OpenMMException("Number of body indices differs from that of atoms in Context");
    bodySystem = new RigidBodySystem(*context, bodyIndices);
    kernel = context->getPlatform().createKernel(IntegrateRigidBodyStepKernel::Name(), contextRef);
    kernel.getAs<IntegrateRigidBodyStepKernel>().initialize(*system, *this);
}

void RigidBodyIntegrator::cleanup() {
    kernel = Kernel();
}

vector<string> RigidBodyIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateRigidBodyStepKernel::Name());
    return names;
}

void RigidBodyIntegrator::stateChanged(State::DataType changed) {
    if (changed == State::Positions || changed == State::Velocities) {
        if (changed == State::Positions) {
            context->updateContextState();
            context->calcForcesAndEnergy(true, false);
            bodySystem->update(true, true);
        }
        else
            bodySystem->update(false, true);
        kernel.getAs<IntegrateRigidBodyStepKernel>().uploadBodySystem(*this->bodySystem);
    }
}

double RigidBodyIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateRigidBodyStepKernel>().computeKineticEnergy(*context, *this);
}

std::vector<int> RigidBodyIntegrator::getBodyIndices() const {
    return bodyIndices;
}

void RigidBodyIntegrator::step(int steps) {
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");
    for (int i = 0; i < steps; ++i)
        kernel.getAs<IntegrateRigidBodyStepKernel>().execute(*context, *this);
}
