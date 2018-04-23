%pythoncode %{

def _isRigid(simulation):
    return isinstance(simulation.integrator, RigidBodyIntegrator)

class StateDataReporter(app.StateDataReporter):

    def __init__(self, *args, **kwargs):
        self._translationalEnergy = kwargs.pop('translationalEnergy', False)
        self._rotationalEnergy = kwargs.pop('rotationalEnergy', False)
        self._modifiedKineticEnergy = kwargs.pop('modifiedKineticEnergy', False)
        self._modifiedPotentialEnergy = kwargs.pop('modifiedPotentialEnergy', False)
        self._modifiedTranslationalEnergy = kwargs.pop('modifiedTranslationalEnergy', False)
        self._modifiedRotationalEnergy = kwargs.pop('modifiedRotationalEnergy', False)
        self._modifiedTotalEnergy = kwargs.pop('_modifiedTotalEnergy', False)
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
        return headers

    def _constructReportValues(self, simulation, state):
        values = super(StateDataReporter, self)._constructReportValues(simulation, state)
        if _isRigid(simulation) and (self._translationalEnergy or self._rotationalEnergy):
            KE = simulation.integrator.getKineticEnergies()
        if self._translationalEnergy:
            value = state.getKineticEnergy() if not _isRigid(simulation) else KE[0]
            values.append(value.value_in_unit(unit.kilojoules_per_mole))
        if self._rotationalEnergy:
            value = 0.0 if not _isRigid(simulation) else KE[1]
            values.append(value.value_in_unit(unit.kilojoules_per_mole))
        return values

%}
