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
    void buildGeometry(vector<Vec3>& R, vector<Vec3>& F, vector<double>& M);
    void buildDynamics(vector<Vec3>& V, vector<double>& M);
    void updateAtomicPositions(vector<Vec3>& R);
    void updateAtomicVelocities(vector<Vec3>& V);
    void forceAndTorque(const vector<Vec3>& F);
    void uniaxialRotationAxis1(double dt);
    void uniaxialRotationAxis2(double dt);
    void uniaxialRotationAxis3(double dt);
    void noSquishRotation(double dt, int n);
    void exactRotation(double dt);

    int    N;             // number of atoms
    int    dof;           // Number of degrees of freedom
    double mass;          // total body mass and its inverse
    Vec3   I;             // Principal moments of inertia and their inverses
    Vec3   rcm;           // Center-of-mass position
    Vec3   pcm;           // Center-of-mass momentum
    Quat   q;             // Unit orientation quaternion
    Quat   pi;            // Quaternion-conjugated momentum

    Vec3   force;         // Resultant force
    Quat   torque;        // Resultant quaternion-fixed torque

    int    loc;           // location of first atom index
    int*   atom;          // pointer to the index of the first atom
    Vec3*  d;             // pointer to the body-fixed position of the first atom

    double invMass;       // total body mass and its inverse
    Vec3   invI;          // Principal moments of inertia and their inverses

    double twoKt;         // Twice the translational kinetic energy
    double twoKr;         // Twice the rotational kinetic energy
    vector<Vec3> delta;   // Space-fixed displacements from the center of mass
};

template<class real, class real3, class real4>
class bodyType {
public:
    int   N;     // number of atoms
    int   loc;   // pointer to set of atoms
    real  invm;  // inverse mass
    real3 invI;  // inverse principal moments of inertia
    real3 r;     // center-of-mass position
    real3 v;     // center-of-mass velocity
    real3 F;     // resultant force
    real3 rdot;  // time-derivative of center-of-mass position
    real4 q;     // orientation quaternion
    real4 pi;    // quaternion-conjugated momentum
    real4 Ctau;  // quaternion-frame resultant torque
    real4 qdot;  // time-derivative of orientation quaternion
};

}
