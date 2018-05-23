OpenMM Rigid Body Plugin
========================

[![Build Status](https://travis-ci.org/craabreu/openmm_rigidbody_plugin.svg?branch=master)](https://travis-ci.org/craabreu/openmm_rigidbody_plugin)
[![codecov.io](http://codecov.io/github/craabreu/openmm_rigidbody_plugin/coverage.svg?branch=master)](http://codecov.io/github/craabreu/openmm_rigidbody_plugin?branch=master)

In this project, the dynamics of rigid bodies in implemented for [OpenMM](http://openmm.org), which
is provided with the following new features:

1. A `RigidBodyIntegrator` class for performing a symplectic, time-reversible, and volume-preserving
integration of the Hamiltonian equations of motion for a system composed of rigid bodies as well as
free atoms.

2. An extension of the `ForceField` class which enables the definition of rigid-body templates and
the creation of systems containing rigid bodies.

3. An extension of the `StateDataReporter` class which is able to display the correct temperature
(considering the reduced number of degrees of freedom), as well as the translational and rotational
parts of the kinetic energy.

This plugin was built upon the official [OpenMM Example Plugin](https://github.com/peastman/openmmexampleplugin).

Formulation
===========

If you use this plugin, please cite the paper which contains the employed formulation:

A. J. Silveira and C. R. A. Abreu, Molecular dynamics with rigid bodies: Alternative formulation
and assessment of its limitations when employed to simulate liquid water, Journal of Chemical
Physics 2017, 147, 124104, doi: [10.1063/1.5003636](https://doi.org/10.1063%2F1.5003636)

Installation
============

In the following instructions, it is assumed that [OpenMM](http://openmm.org) is already installed
in your system.

Downloading and compiling the source code
-----------------------------------------

Clone the git repository locally:

    git clone https://github.com/craabreu/openmm_rigidbody_plugin.git

You will need [CMake](http://www.cmake.org) to build the plugin. Once it is installed, you can make:

    cd openmm_rigidbody_plugin
    mkdir build && cd build
    cmake ..
    make
    make install
    make PythonInstall

The commands `make install` and `make PythonInstall` might require administrator privileges. The
former will install the plugin libraries and the latter will install the python module. In order
to test the installation, please run:

    make test

Python API Extension
====================

```python
class ForceField(simtk.openmm.app.ForceField):

    def registerBodyTemplate(self, name, residue, pattern=None):
    """
        Register a new rigid body template.

        Parameters
        ----------
        name : string
            A name for the rigid body template being registered.
        residue : string
            The name of the residue template to which this body belongs.
        pattern : string, optional, default=None
            A string containing a regular expression for selecting the atoms in residue which will
            form the body. A list of atom names can be passed as '(name1|name2|...)', for instance.
            Only actual atoms will be considered (i.e. virtual sites will be ignored). If patttern
            is None, then the whole residue will form the rigid body.
      """

      def createSystem(self, topology, **kwargs):
      """
          Construct an OpenMM System representing a Topology with this force field. This extends
          the original method by creating rigid bodies from previously registered body templates.

          Parameters
          ----------
          topology : Topology
              The Topology for which to create a System.
          kwargs : multiple types
              All keyword arguments of the original `ForceField.createSystem` method, plus:
          mergeList : list(list(int)), optional, default=None
              A nested list of indices of rigid bodies to be merged together after the templates
              are applied. If this is None, then no merging will occur.
          removeConstraints : bool, optional, default=False
              If true, every constraint involving at least one atom that belongs to a rigid body
              will be silently removed from the system.
          removeForces : bool, optional, default=False
              If true, all interactions whose all involved atoms belong to the same rigid body
              will be silently removed from the system.

          Returns
          -------
          system : System
              The newly created System.
          bodyIndices : list(int)
              A list containing the index of the rigid body to which every atom belongs (note:
              index=0 for free atoms).
      """
```

```python
class RigidBodyIntegrator(simtk.openmm.Integrator):

    def __init__(self, stepSize, bodyIndices):
    """
        Create a Rigid Body Integrator.

        Parameters
        ----------
        stepSize : double
            The step size with which to integrate the system (in picoseconds).
        bodyIndices : list(int)
            A list containing the index of the rigid body to which every atom belongs (note:
            index=0 for free atoms).
    """

    def setRotationMode(mode):
    """
        Defines the mode of rotational move integration.

        Parameters
        ----------
        mode : int
            The mode of rotation, where mode=0 consisting in using an exact solution for rotations
            spanning the whole time step, while mode=n consists in using a NO-SQUISH solution with
            the time step split into n small steps.
    """

    def step(self, steps):
    """
       	Advance a simulation through time by taking a series of time steps.

        Parameters
        ----------
        steps : int
            The number of time steps to be taken.
    """
```

```python
class StateDataReporter(simtk.openmm.app.StateDataReporter):

    def __init__(self, *args, **kwargs):
    """
        Create a State Data Reporter.

        Parameters
        ----------
        args : multiple types
            All positional arguments of the original `StateDataReporter` constructor.
        kwargs : multiple types
            All keyword arguments of the original `StateDataReporter` constructor, plus:
        translationalEnergy : bool, optional, default=False
            Whether to write the translational part of the kinetic energy to the file.
        rotationalEnergy : bool, optional, default=False
            Whether to write the rotational part of the kinetic energy to the file.
    """
```

Python Code Example
===================

```python
import simtk.openmm as mm
from simtk.openmm import app
import rigidbodyplugin as rbmm
pdb = app.PDBFile('config.pdb')
forcefield = rbmm.ForceField('model.xml')
forcefield.registerBodyTemplate('water', 'HOH')
(system, bodies) = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME,
                                           removeConstraints=True, removeForces=True)
integrator = rbmm.RigidBodyIntegrator(1.0*unit.femtoseconds, bodies)
integrator.setRotationMode(0)
reporter = rbmm.StateDataReporter(stdout, 1, step=True, temperature=True,
                                  translationalEnergy=True, rotationalEnergy=True)
simulation = app.Simulation(pdb.topology, system, integrator,
                            platform=mm.Platform.getPlatformByName('CUDA'),
                            properties={'Precision': 'mixed'})
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
simulation.reporters.append(reporter)
simulation.step(1000)
```

CPU and GPU Computing
=====================

In its current version, this plugin can be executed using either the Reference or the CUDA platform
of OpenMM. The OpenCL platform is not supported.

License
=======

MIT License

Copyright (c) 2018, Charlles R. A. Abreu

Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
