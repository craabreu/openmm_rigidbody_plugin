/* -------------------------------------------------------------------------- *
 *                          OpenMM Rigid Body Plugin                          *
 * -------------------------------------------------------------------------- */

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

class RigidBody {
public:
    void updateGeometry(vector<Vec3>& R, vector<double>& M);
    void updateVelocities(vector<Vec3>& V, vector<double>& M);

    int    N = 0;         // number of atoms
    int    loc;           // location of first atom index
    int*   atom;          // pointer to the index of the first atom
    Vec3*  d;             // pointer to the body-fixed position of the first atom
private:
    int    dof;           // Number of degrees of freedom
    double mass, invMass; // total body mass and its inverse
    Vec3   MoI, invMoI;   // Principal moments of inertia and their inverses
    Vec3   rcm;           // Center-of-mass position
    Vec3   vcm;           // Center-of-mass velocity
    Quat   q;             // Unit orientation quaternion
    Quat   pi;            // Quaternion-conjugated momentum
};


class RigidBodySystem {
public:
    RigidBodySystem(ContextImpl& contextRef, const vector<int>& bodyIndices);

    void update(bool geometry, bool velocities);
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
    bool upToDate = false;
};

}
