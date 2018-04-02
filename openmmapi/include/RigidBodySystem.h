/* -------------------------------------------------------------------------- *
 *                          OpenMM Rigid Body Plugin                          *
 * -------------------------------------------------------------------------- */

#include "Vec4.h"
#include "openmm/Vec3.h"
#include "openmm/Context.h"
#include <vector>

using namespace OpenMM;
using namespace std;

namespace RigidBodyPlugin {

/**
 * This class represents a system of rigid bodies and free atoms
 */

class RigidBody {
public:
    void update(vector<int> atomIndices, ContextImpl& context);

    int    N = 0;       // number of atoms
    int    loc;         // location of first atom index
    int*   atom;        // pointer to the location of first atom index
    double mass = 0.0;  // total body mass
private:
    int    dof = 6; // Number of degrees of freedom

    OpenMM::Vec3   MoI;     // Principal moments of inertia
    OpenMM::Vec3   rcm;     // Center-of-mass position
    OpenMM::Vec3   pcm;     // Center-of-mass momentum vector
    OpenMM::Vec3   f;       // Resultant force
    OpenMM::Vec3   omega;   // Angular velocities
    OpenMM::Vec3   tau;     // Resultant torque

    Vec4   q;       // Unit quaternion of orientation
    Vec4   pi;      // Quaternion momentum

    vector<double*> M;
    vector<Vec3*>   R;
    vector<Vec3*>   V;
    vector<Vec3*>   F;
};


class RigidBodySystem {
public:
    RigidBodySystem(const System* system, const vector<int>& bodyIndices);
private:
    std::vector<int> atomIndex;
    std::vector<int> bodyIndex;
    std::vector<RigidBody> body;
    int numFree;
    int numBodies;
};

}
