/* -------------------------------------------------------------------------- *
 *                          OpenMM Rigid Body Plugin                          *
 * -------------------------------------------------------------------------- */

#include "RigidBody.h"
#include "internal/MatVec.h"
#include "openmm/Vec3.h"
#include "openmm/internal/ContextImpl.h"
#include <vector>

using namespace OpenMM;
using namespace std;

namespace RigidBodyPlugin {

/**
 * This class represents a system of rigid bodies and free atoms
 */

class RigidBodySystem {
public:
    RigidBodySystem(ContextImpl& contextRef, const vector<int>& bodyIndices);

    void update(bool geometry, bool velocities);
    int getNumFree() const { return numFree; }
    int getNumBodies() const { return numBodies; }
    int getNumActualAtoms() const { return numActualAtoms; }
    int getNumBodyAtoms() const { return numBodyAtoms; }
    int getAtomIndex(int i) const { return atomIndex[i]; }
    Vec3 getBodyFixedPosition(int i) const { return bodyFixedPositions[i]; }
    RigidBody getRigidBody(int i) const { return body[i]; }
private:
    ContextImpl* context;
    std::vector<int> bodyIndex;
    std::vector<int> atomIndex;
    std::vector<RigidBody> body;
    std::vector<OpenMM::Vec3> bodyFixedPositions;
    int numFree;
    int numBodies;
    int numActualAtoms;
    int numBodyAtoms;
};

}
