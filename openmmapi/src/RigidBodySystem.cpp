/* -------------------------------------------------------------------------- *
 *                          OpenMM Rigid Body Plugin                          *
 * -------------------------------------------------------------------------- */

#include "RigidBodySystem.h"
#include "Vec4.h"
#include "openmm/Vec3.h"
#include "openmm/Context.h"
#include "openmm/System.h"
#include <vector>

using namespace RigidBodyPlugin;
using namespace OpenMM;
using std::vector;

void RigidBody::update(vector<int> atomIndices, ContextImpl& context) {
}

/*
 * Creates a clean list of rigid body indices. The rule is: if bodyIndices[i] <= 0 or
 * bodyIndices[i] is unique, then index[i] = 0. Otherwise, 1 <= index[i] <= number of
 * distinct positive entries which are non-unique.
 */

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

/*
 * Creates a data structure for the system of rigid bodies and free atoms
 */
RigidBodySystem::RigidBodySystem(const System* system, const vector<int>& bodyIndices) {

    bodyIndex = cleanBodyIndices(bodyIndices);

    numBodies = 0;
    for (auto index : bodyIndex)
        numBodies = std::max(numBodies, index);

    int numAtoms = system->getNumParticles();
    int numActualAtoms = numAtoms;
    for (int i = 0; i < numAtoms; i++)
        if (system->isVirtualSite(i))
            numActualAtoms--;
    atomIndex.resize(numActualAtoms);

    int numFree = 0;
    body.resize(numBodies);
    for (int i = 0; i < numAtoms; i++)
        if (!system->isVirtualSite(i)) {
            int ibody = bodyIndex[i];
            if (ibody == 0)
                atomIndex[numFree++] = i;
            else
                body[ibody-1].N++;
        }

    int loc = numFree;
    for (auto& b : body) {
        b.loc = loc;
        b.atom = &atomIndex[loc];
        loc += b.N;
    }

//    cout<<"Number of bodies = "<<numBodies<<"\n"
//        <<"Number of actual atoms = "<<numActualAtoms<<"\n"
//        <<"Number of free atoms = "<<numFree<<"\n";

    vector<int> iatom(numBodies, 0);
    for (int i = 0; i < numAtoms; i++) {
        int ibody = bodyIndex[i];
        if (ibody > 0) {
            body[ibody-1].atom[iatom[ibody-1]++] = i;
        }
    }
}
