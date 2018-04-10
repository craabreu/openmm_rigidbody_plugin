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
    void updateGeometry(vector<Vec3>& R, vector<Vec3>& F, vector<double>& M);
    void updateVelocities(vector<Vec3>& V, vector<double>& M);

    int    N = 0;         // number of atoms
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
private:
    int    dof;           // Number of degrees of freedom
    double invMass;       // inverse of total body mass
    Vec3   invMoI;        // inverses of principal moments of inertia
};


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
