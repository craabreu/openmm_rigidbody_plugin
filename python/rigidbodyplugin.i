%module rigidbodyplugin

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

%include "std_vector.i"
namespace std {
  %template(vectori) vector<int>;
  %template(vectord) std::vector<double>;
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

%pythonappend RigidBodyPlugin::RigidBodySystem::getTranslationalEnergy() const %{
    val = unit.Quantity(val, unit.kilojoule_per_mole)
%}

%pythonappend RigidBodyPlugin::RigidBodySystem::getRotationalEnergy() const %{
    val = unit.Quantity(val, unit.kilojoule_per_mole)
%}

%pythonappend RigidBodyPlugin::RigidBodySystem::getKineticEnergy() const %{
    val = unit.Quantity(val, unit.kilojoule_per_mole)
%}

%pythonappend RigidBodyPlugin::RigidBodyIntegrator::getKineticEnergies() %{
    val = tuple(unit.Quantity(v, unit.kilojoule_per_mole) for v in val)
%}

%pythonappend RigidBodyPlugin::RigidBodyIntegrator::getRefinedKineticEnergies() %{
    val = tuple(unit.Quantity(v, unit.kilojoule_per_mole) for v in val)
%}

%pythonappend RigidBodyPlugin::RigidBodyIntegrator::getPotentialEnergyRefinement() %{
    val = unit.Quantity(val, unit.kilojoule_per_mole)
%}

namespace RigidBodyPlugin {

class RigidBodySystem {
public:
    int getNumDOF() const;
    double getTranslationalEnergy() const;
    double getRotationalEnergy() const;
    double getKineticEnergy() const;
};

class RigidBodyIntegrator : public OpenMM::Integrator {
public:
    explicit RigidBodyIntegrator(double stepSize, const std::vector<int>& bodyIndices);
    void setRotationMode(int mode);
    int getRotationMode();
    void setComputeRefinedEnergies(bool compute);
    bool getComputeRefinedEnergies();
    RigidBodySystem getRigidBodySystem();
    std::vector<double> getKineticEnergies();
    std::vector<double> getRefinedKineticEnergies();
    double getPotentialEnergyRefinement();
    void step(int steps);
};

}
