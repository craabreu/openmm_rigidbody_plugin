#ifndef OPENMM_RIGIDBODY_INTEGRATOR_PROXY_H_
#define OPENMM_RIGIDBODY_INTEGRATOR_PROXY_H_

#include "internal/windowsExportRigidBody.h"
#include "openmm/serialization/SerializationProxy.h"
//#include "openmm/serialization/XmlSerializer.h"

namespace RigidBodyPlugin {

class RigidBodyIntegratorProxy : public OpenMM::SerializationProxy {
public:
    RigidBodyIntegratorProxy();
    void serialize(const void* object, OpenMM::SerializationNode& node) const;
    void* deserialize(const OpenMM::SerializationNode& node) const;
};

}

#endif /*OPENMM_RIGIDBODY_INTEGRATOR_PROXY_H_*/
