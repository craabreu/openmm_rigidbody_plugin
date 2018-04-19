%module rigidbodyplugin

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

%include "std_vector.i"
namespace std {
  %template(vectori) vector<int>;
};

%{
#include "RigidBodyIntegrator.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%include "imports.py"
%include "forcefield.py"
%include "statedatareporter.py"

namespace RigidBodyPlugin {

class RigidBodySystem {
public:
    explicit RigidBodySystem(ContextImpl& contextRef, const vector<int>& bodyIndices);
    int getNumDOF() const {return numDOF; }
};

class RigidBodyIntegrator : public OpenMM::Integrator {
public:
    explicit RigidBodyIntegrator(double stepSize, const std::vector<int>& bodyIndices);
    void setRotationMode(int mode);
    int getRotationMode();
    RigidBodySystem getRigidBodySystem();
    void step(int steps);
};

}
