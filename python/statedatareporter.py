%pythoncode %{

def _isRigid(simulation):
    return isinstance(simulation.integrator, RigidBodyIntegrator)

class StateDataReporter(app.StateDataReporter):

    def __init__(self, *args, **kwargs):
        self._translationalEnergy = kwargs.pop('translationalEnergy', False)
        self._rotationalEnergy = kwargs.pop('rotationalEnergy', False)
        self._refinedPotentialEnergy = kwargs.pop('refinedPotentialEnergy', False)
        self._refinedKineticEnergy = kwargs.pop('refinedKineticEnergy', False)
        self._refinedTotalEnergy = kwargs.pop('refinedTotalEnergy', False)
        self._refinedTemperature = kwargs.pop('refinedTemperature', False)
        self._refinedTranslationalEnergy = kwargs.pop('refinedTranslationalEnergy', False)
        self._refinedRotationalEnergy = kwargs.pop('refinedRotationalEnergy', False)
        super(StateDataReporter, self).__init__(*args, **kwargs)

    def _initializeConstants(self, simulation):
        super(StateDataReporter, self)._initializeConstants(simulation)
        if (self._temperature or self._refinedTemperature) and _isRigid(simulation):
            self._dof = simulation.integrator.getRigidBodySystem().getNumDOF()
            system = simulation.system
            forces = [system.getForce(i) for i in range(system.getNumForces())]
            if any(isinstance(f, mm.CMMotionRemover) for f in forces):
                self._dof -= 3
        elif self._refinedTemperature:
            self._dof = 1

    def _constructHeaders(self):
        headers = super(StateDataReporter, self)._constructHeaders()
        if self._translationalEnergy:
            headers.append('Translational Energy (kJ/mole)')
        if self._rotationalEnergy:
            headers.append('Rotational Energy (kJ/mole)')
        if self._refinedPotentialEnergy:
            headers.append('Refined Potential Energy (kJ/mole)')
        if self._refinedKineticEnergy:
            headers.append('Refined Kinetic Energy (kJ/mole)')
        if self._refinedTotalEnergy:
            headers.append('Refined Total Energy (kJ/mole)')
        if self._refinedTemperature:
            headers.append('Refined Temperature (K)')
        if self._refinedTranslationalEnergy:
            headers.append('Refined Translational Energy (kJ/mole)')
        if self._refinedRotationalEnergy:
            headers.append('Refined Rotational Energy (kJ/mole)')
        return headers

    def _constructReportValues(self, simulation, state):
        values = super(StateDataReporter, self)._constructReportValues(simulation, state)

        if (self._translationalEnergy or self._rotationalEnergy):
            if _isRigid(simulation):
                KE = simulation.integrator.getKineticEnergies()
            else:
                KE = [state.getKineticEnergy(), 0.0*unit.kilojoules_per_mole]
            if self._translationalEnergy:
                values.append(KE[0].value_in_unit(unit.kilojoules_per_mole))
            if self._rotationalEnergy:
                values.append(KE[1].value_in_unit(unit.kilojoules_per_mole))

        if (self._refinedPotentialEnergy or self._refinedTotalEnergy):
            U = state.getPotentialEnergy()
            if _isRigid(simulation):
                U += simulation.integrator.getPotentialEnergyRefinement()
            if self._refinedPotentialEnergy:
                values.append(U.value_in_unit(unit.kilojoules_per_mole))

        if (self._refinedKineticEnergy or self._refinedTotalEnergy or self._refinedTemperature or
            self._refinedTranslationalEnergy or self._refinedRotationalEnergy):
            if _isRigid(simulation):
                KE = simulation.integrator.getRefinedKineticEnergies()
            else:
                KE = [state.getKineticEnergy(), 0.0*unit.kilojoules_per_mole]
            if self._refinedKineticEnergy:
                values.append((KE[0] + KE[1]).value_in_unit(unit.kilojoules_per_mole))
            if self._refinedTotalEnergy:
                values.append((KE[0] + KE[1] + U).value_in_unit(unit.kilojoules_per_mole))
            if self._refinedTemperature:
                T = 2*(KE[0] + KE[1])/(self._dof*unit.MOLAR_GAS_CONSTANT_R)
                values.append(T.value_in_unit(unit.kelvin))
            if self._refinedTranslationalEnergy:
                values.append(KE[0].value_in_unit(unit.kilojoules_per_mole))
            if self._refinedRotationalEnergy:
                values.append(KE[1].value_in_unit(unit.kilojoules_per_mole))

        return values

%}
