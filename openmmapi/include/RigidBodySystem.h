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
 
 class FreeAtom {
 public:
    int index;
    double invMass;
    OpenMM::Vec3 savedPos;
 };

class RigidBodySystem {
public:
    RigidBodySystem() {}
    void initialize(ContextImpl& contextRef, const vector<int>& bodyIndices);
    void update(ContextImpl& contextRef, bool geometry, bool velocities);
    void copy(const RigidBodySystem& bodySystem);
    void integratePart1(double dt, const vector<Vec3>& F, vector<Vec3>& V, vector<Vec3>& R);
    void integratePart2(double dt, const vector<Vec3>& R, const vector<Vec3>& F, vector<Vec3>& V);
    void computeKineticEnergies(const vector<Vec3>& V);
    int getNumDOF() const {return numDOF; }
    int getNumFree() const { return numFree; }
    int getNumBodies() const { return numBodies; }
    int getNumActualAtoms() const { return numActualAtoms; }
    int getNumBodyAtoms() const { return numBodyAtoms; }
    int getAtomIndex(int i) const { return atomIndex[i]; }
    Vec3 getBodyFixedPosition(int i) const { return bodyFixedPositions[i]; }
    RigidBody getRigidBody(int i) const { return body[i]; }
    double getTranslationalEnergy() const { return transKE; }
    double getRotationalEnergy() const { return rotKE; }
    double getKineticEnergy() const { return transKE + rotKE; }
private:
    std::vector<int> bodyIndex;
    std::vector<int> atomIndex;
    std::vector<RigidBody> body;
    std::vector<OpenMM::Vec3> bodyFixedPositions;
    std::vector<FreeAtom> freeAtom;
    int numFree;
    int numBodies;
    int numActualAtoms;
    int numBodyAtoms;
    int numDOF;
    double transKE;
    double rotKE;
};

}
