/* -------------------------------------------------------------------------- *
 *                          OpenMM Rigid Body Plugin                          *
 * -------------------------------------------------------------------------- */

#include "internal/MatVec.h"
#include "openmm/Vec3.h"
#include <vector>

using namespace OpenMM;
using namespace std;

namespace RigidBodyPlugin {

/**
 * This class represents a rigid body
 */

class RigidBody {
public:
    void updateGeometry(vector<Vec3>& R, vector<Vec3>& F, vector<double>& M);
    void updateVelocities(vector<Vec3>& V, vector<double>& M);

    int    N = 0;         // number of atoms
    int    dof;           // Number of degrees of freedom
    double mass;          // total body mass and its inverse
    Vec3   MoI;           // Principal moments of inertia and their inverses
    Vec3   rcm;           // Center-of-mass position
    Vec3   pcm;           // Center-of-mass momentum
    Quat   q;             // Unit orientation quaternion
    Quat   pi;            // Quaternion-conjugated momentum

    Vec3   force;
    Quat   torque;

    int    loc;           // location of first atom index
    int*   atom;          // pointer to the index of the first atom
    Vec3*  d;             // pointer to the body-fixed position of the first atom

    double invMass;       // total body mass and its inverse
    Vec3   invMoI;        // Principal moments of inertia and their inverses
};

}
