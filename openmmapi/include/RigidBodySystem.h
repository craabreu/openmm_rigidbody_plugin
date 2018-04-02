/* -------------------------------------------------------------------------- *
 *                          OpenMM Rigid Body Plugin                          *
 * -------------------------------------------------------------------------- */

#include "Vec4.h"
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
    void update(vector<Vec3>& R, vector<Vec3>& V, vector<double>& M);
    void updateVelocities(vector<Vec3>& V, vector<double>& M);

    int    N = 0;       // number of atoms
    int    loc;         // location of first atom index
    int*   atom;        // pointer to the index of the first atom
    Vec3*  d;           // pointer to the body-fixed position of the first atom
private:
    int    dof = 6; // Number of degrees of freedom
    double mass = 0.0;  // total body mass
    Vec3   MoI;     // Principal moments of inertia
    Vec3   rcm;     // Center-of-mass position
    Vec3   vcm;     // Center-of-mass velocity
    Vec3   f;       // Resultant force
    Vec3   omega;   // Angular velocities
    Vec3   tau;     // Resultant torque

    Vec4   q;       // Unit quaternion of orientation
    Vec4   pi;      // Quaternion momentum

    vector<double*> M;
    vector<Vec3*>   R;
    vector<Vec3*>   V;
    vector<Vec3*>   F;
};


class RigidBodySystem {
public:
    RigidBodySystem(ContextImpl& contextRef, const vector<int>& bodyIndices);

    void update();
    void updateVelocities();
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
