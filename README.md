OpenMM Rigid Body Plugin
========================

[![Build Status](https://travis-ci.org/craabreu/openmm_rigidbody_plugin.svg?branch=master)](https://travis-ci.org/craabreu/openmm_rigidbody_plugin)
[![codecov.io](http://codecov.io/github/craabreu/openmm_rigidbody_plugin/coverage.svg?branch=master)](http://codecov.io/github/craabreu/openmm_rigidbody_plugin?branch=master)

This project is an implementation of rigid-body dynamics for [OpenMM](http://openmm.org). In its
current version, it provides OpenMM with the following new features:

1. A class for performing symplectic, time-reversible, volume-preserving integration of Hamiltonian
equations of motion for a system composed of rigid bodies and free atoms.

2. An extension of the `ForceField` class which allows the definition of rigid-body templates and
creation of systems containing rigid bodies.

3. An extension of the `StateDataReporter` class which is able to compute the correct temperature
(considering the reduced number of degrees of freedom) and report the translational and rotational
parts of the kinetic energy.

If you use this plugin, please cite the paper which contains the employed formulation:

A. J. Silveira and C. R. A. Abreu, Molecular dynamics with rigid bodies: Alternative formulation
and assessment of its limitations when employed to simulate liquid water, Journal of Chemical
Physics 2017, 147, 124104, doi: [10.1063/1.5003636](https://doi.org/10.1063%2F1.5003636)

This plugin was built upon the official [OpenMM Example Plugin](https://github.com/peastman/openmmexampleplugin).


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


CPU and GPU Computing
=====================

In its current version, this plugin can be executed using either the Reference or the CUDA platform
of OpenMM. The OpenCL platform is not supported.


Python API Extension
====================

--> Content will be added soon.

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
