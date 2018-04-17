/* ---------------------------------------------------------------------------------------------- *
 *                                    OpenMM Rigid Body Plugin                                    *
 * ---------------------------------------------------------------------------------------------- */

#include "RigidBodySystem.h"
#include "internal/MatVec.h"
#include "openmm/Vec3.h"
#include "openmm/System.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "internal/eigenDecomposition.h"
#include <vector>
#include <cmath>

#include <iostream> // TEMPORARY

using namespace RigidBodyPlugin;
using namespace OpenMM;
using std::vector;

/*--------------------------------------------------------------------------------------------------
  Create a clean list of rigid body indices. The rule is: if bodyIndices[i] <= 0 or is unique, then
  index[i] = 0. Otherwise, 1 <= index[i] <= number of distinct positive entries which are
  non-unique.
--------------------------------------------------------------------------------------------------*/

vector<int> cleanBodyIndices(const vector<int>& bodyIndices) {
    int size = bodyIndices.size();
    vector<int> index(size, 0), saved(size, 0), amount(size, 0), first(size, 0);
    int nsaved = 0;
    for (int i = 0; i < size; i++) {
        int ibody = bodyIndices[i];
        if (ibody > 0) {
            int j;
            int found = 0;
            for (j = 0; j < nsaved; j++) {
                found = saved[j] == ibody;
                if (found) break;
            }
            if (found) {
                amount[j]++;
                index[i] = j+1;
                index[first[j]] = j+1;
            }
            else {
                amount[nsaved] = 1;
                saved[nsaved] = ibody;
                first[nsaved] = i;
                index[i] = 0;
                nsaved++;
            }
        }
        else
            index[i] = 0;
    }
    int n = 0;
    for (int i = 0; i < nsaved; i++)
        saved[i] = amount[i] > 1 ? ++n : 0;
    for (int i = 0; i < size; i++)
        if (index[i] > 0)
            index[i] = saved[index[i]-1];
    return index;
}

/*--------------------------------------------------------------------------------------------------
  Create a data structure for the system of rigid bodies and free atoms.
--------------------------------------------------------------------------------------------------*/

RigidBodySystem::RigidBodySystem(ContextImpl& contextRef, const vector<int>& bodyIndices) {
    context = &contextRef;
    bodyIndex = cleanBodyIndices(bodyIndices);

    numBodies = 0;
    for (auto index : bodyIndex)
        numBodies = std::max(numBodies, index);

    const System *system = &context->getSystem();
    int numAtoms = system->getNumParticles();
    numActualAtoms = numAtoms;
    for (int i = 0; i < numAtoms; i++)
        if (system->isVirtualSite(i))
            numActualAtoms--;
    atomIndex.resize(numActualAtoms);

    numFree = 0;
    body.resize(numBodies);
    for (int i = 0; i < numAtoms; i++)
        if (!system->isVirtualSite(i)) {
            int ibody = bodyIndex[i];
            if (ibody == 0)
                atomIndex[numFree++] = i;
            else
                body[ibody-1].N++;
        }
    numBodyAtoms = numActualAtoms - numFree;
    bodyFixedPositions.resize(numBodyAtoms);

    int loc = 0;
    for (auto& b : body) {
        b.loc = loc;
        b.atom = &atomIndex[numFree+loc];
        b.d = &bodyFixedPositions[loc];
        loc += b.N;
    }

    vector<int> iatom(numBodies, 0);
    for (int i = 0; i < numAtoms; i++) {
        int ibody = bodyIndex[i];
        if (ibody > 0) {
            body[ibody-1].atom[iatom[ibody-1]++] = i;
        }
    }

    cout<<"Number of bodies = "<<numBodies<<"\n"               // TEMPORARY
        <<"Number of actual atoms = "<<numActualAtoms<<"\n"    // TEMPORARY
        <<"Number of free atoms = "<<numFree<<"\n";            // TEMPORARY

    for (int i = 0; i < system->getNumConstraints(); i++) {
        int atom1, atom2;
        double distance;
        system->getConstraintParameters(i, atom1, atom2, distance);
        if (bodyIndex[atom1] != 0 || bodyIndex[atom2] != 0)
            throw OpenMMException("Constraints involving rigid-body atoms are not allowed");
    }
}

/*--------------------------------------------------------------------------------------------------
  Update the kinematic properties of all rigid bodies.
--------------------------------------------------------------------------------------------------*/

void RigidBodySystem::update(bool geometry, bool velocities) {
    const System* system = &context->getSystem();
    int N = system->getNumParticles();
    vector<double> M(N);
    for (int i = 0; i < N; i++)
        M[i] = system->getParticleMass(i);
    if (geometry) {
        vector<Vec3> R(N), F(N);
        context->getPositions(R);
        context->getForces(F);
        for (auto& b : body)
            b.updateGeometry(R, F, M);
    }
    if (velocities) {
        vector<Vec3> V(N);
        context->getVelocities(V);
        for (auto& b : body)
            b.updateVelocities(V, M);
    }
}
