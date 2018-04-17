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
    int getNumFree() { return numFree; }
    int getNumBodies() { return numBodies; }
    int getNumActualAtoms() { return numActualAtoms; }
    int getNumBodyAtoms() { return numBodyAtoms; }
    int getAtomIndex(int i) { return atomIndex[i]; }
    Vec3 getBodyFixedPosition(int i) { return bodyFixedPositions[i]; }
    RigidBody* getRigidBody(int i) { return &body[i]; }
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
