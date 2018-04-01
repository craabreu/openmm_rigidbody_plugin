/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2015 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Yutong Zhao                                        *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/internal/AssertionUtilities.h"
#include "RigidBodyIntegrator.h"
#include "openmm/serialization/XmlSerializer.h"
#include <iostream>
#include <sstream>
#include <stdlib.h>

#include <fstream>
#include <vector>

using namespace RigidBodyPlugin;
using namespace OpenMM;
using namespace std;

extern "C" void registerRigidBodySerializationProxies();

void testSerializeRigidBodyIntegrator() {
    vector<int> bodyIndices;
    int indices[] = {5, 4, 3, 2, 1};
    bodyIndices.assign(indices, indices+5);
    RigidBodyIntegrator *intg = new RigidBodyIntegrator(0.00342, bodyIndices);
    stringstream ss;
    XmlSerializer::serialize<Integrator>(intg, "RigidBodyIntegrator", ss);
    RigidBodyIntegrator *intg2 = dynamic_cast<RigidBodyIntegrator*>(XmlSerializer::deserialize<Integrator>(ss));
    ASSERT_EQUAL(intg->getConstraintTolerance(), intg2->getConstraintTolerance());
    ASSERT_EQUAL(intg->getStepSize(), intg2->getStepSize());
    int i = 0;
    for (auto index : intg2->getBodyIndices())
        ASSERT_EQUAL(index, indices[i++]);
    delete intg;
    delete intg2;
}

int main() {
    try {
        testSerializeRigidBodyIntegrator();
    }
    catch(const exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
