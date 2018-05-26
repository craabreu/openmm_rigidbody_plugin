%pythoncode %{

def _isRigid(simulation):
    return isinstance(simulation.integrator, RigidBodyIntegrator)

class StateDataReporter(app.StateDataReporter):

    def __init__(self, *args, **kwargs):
        self._translationalEnergy = kwargs.pop('translationalEnergy', False)
        self._rotationalEnergy = kwargs.pop('rotationalEnergy', False)
        self._refinedKineticEnergy = kwargs.pop('refinedKineticEnergy', False)
#         self._refinedPotentialEnergy = kwargs.pop('refinedPotentialEnergy', False)
        self._refinedTranslationalEnergy = kwargs.pop('refinedTranslationalEnergy', False)
        self._refinedRotationalEnergy = kwargs.pop('refinedRotationalEnergy', False)
#         self._refinedTotalEnergy = kwargs.pop('_refinedTotalEnergy', False)
        super(StateDataReporter, self).__init__(*args, **kwargs)

    def _initializeConstants(self, simulation):
        super(StateDataReporter, self)._initializeConstants(simulation)
        if self._temperature and _isRigid(simulation):
            self.dof = simulation.integrator.getRigidBodySystem().getNumDOF()
            system = simulation.system
            forces = [system.getForce(i) for i in range(system.getNumForces())]
            if any(isinstance(f, mm.CMMotionRemover) for f in forces):
                dof -= 3

    def _constructHeaders(self):
        headers = super(StateDataReporter, self)._constructHeaders()
        if self._translationalEnergy:
            headers.append('Translational Energy (kJ/mole)')
        if self._rotationalEnergy:
            headers.append('Rotational Energy (kJ/mole)')
        if self._refinedKineticEnergy:
            headers.append('Refined Kinetic Energy (kJ/mole)')
        if self._refinedTranslationalEnergy:
            headers.append('Refined Translational Energy (kJ/mole)')
        if self._refinedRotationalEnergy:
            headers.append('Refined Rotational Energy (kJ/mole)')
        return headers

    def _constructReportValues(self, simulation, state):
        values = super(StateDataReporter, self)._constructReportValues(simulation, state)
        if _isRigid(simulation):
            if (self._translationalEnergy or self._rotationalEnergy):
                KE = simulation.integrator.getKineticEnergies()
            if (self._refinedKineticEnergy or self._refinedTranslationalEnergy or self._refinedRotationalEnergy):
                modKE = simulation.integrator.getRefinedKineticEnergies()
        elif (self._translationalEnergy or self._refinedTranslationalEnergy or self._refinedKineticEnergy):
            KE = modKE = [state.getKineticEnergy(), 0.0]
        if self._translationalEnergy:
            value = KE[0]
            values.append(value.value_in_unit(unit.kilojoules_per_mole))
        if self._rotationalEnergy:
            value = KE[1]
            values.append(value.value_in_unit(unit.kilojoules_per_mole))
        if self._refinedKineticEnergy:
            value = modKE[0] + modKE[1]
            values.append(value.value_in_unit(unit.kilojoules_per_mole))
        if self._refinedTranslationalEnergy:
            value = modKE[0]
            values.append(value.value_in_unit(unit.kilojoules_per_mole))
        if self._refinedRotationalEnergy:
            value = modKE[1]
            values.append(value.value_in_unit(unit.kilojoules_per_mole))
        return values

%}
