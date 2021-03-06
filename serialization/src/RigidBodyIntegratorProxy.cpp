/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

#include "RigidBodyIntegratorProxy.h"
#include "RigidBodyIntegrator.h"
#include <OpenMM.h>

using namespace std;
using namespace RigidBodyPlugin;
using namespace OpenMM;

RigidBodyIntegratorProxy::RigidBodyIntegratorProxy() : SerializationProxy("RigidBodyIntegrator") {
}

void RigidBodyIntegratorProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const RigidBodyIntegrator& integrator = *reinterpret_cast<const RigidBodyIntegrator*>(object);
    node.setDoubleProperty("stepSize", integrator.getStepSize());
    node.setDoubleProperty("constraintTolerance", integrator.getConstraintTolerance());
    vector<int> bodyIndices;
    bodyIndices = integrator.getBodyIndices();
    SerializationNode& bodyIndicesNode = node.createChildNode("bodyIndices");
    for (auto index : bodyIndices)
        bodyIndicesNode.createChildNode("bodyIndex").setIntProperty("index", index);
}

void* RigidBodyIntegratorProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    double stepSize = node.getDoubleProperty("stepSize");
    double constraintTolerance = node.getDoubleProperty("constraintTolerance");
    const SerializationNode& bodyIndicesNode = node.getChildNode("bodyIndices");
    vector<int> bodyIndices;
    for (auto& child : bodyIndicesNode.getChildren())
        bodyIndices.push_back(child.getIntProperty("index"));
    RigidBodyIntegrator *integrator = new RigidBodyIntegrator(stepSize, bodyIndices);
    integrator->setConstraintTolerance(constraintTolerance);
    return integrator;
}
