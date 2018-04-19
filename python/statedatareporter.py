%pythoncode %{

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
        if self._temperature and isinstance(simulation.integrator, RigidBodyIntegrator):
            self.dof = simulation.integrator.getRigidBodySystem().getNumDOF()
            system = simulation.system
            forces = [system.getForce(i) for i in range(system.getNumForces())]
            if any(isinstance(f, mm.CMMotionRemover) for f in forces):
                dof -= 3

%}
