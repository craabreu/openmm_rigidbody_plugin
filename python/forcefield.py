%pythoncode %{

class ForceField(app.ForceField):
    def __init__(self, *files):
        super(ForceField, self).__init__(*files)
        self._bodyTemplates = {}

    class _bodyTemplateData(object):
        """Inner class used to encapsulate data about a rigid body template definition."""
        def __init__(self, residue, atoms):
            self.residue = residue
            self.atoms = set(atoms)

        def __str__(self):
            return "%s: [%s]" % (self.residue, " ".join(self.atoms))

    def getBodyTemplate(self, name):
        return self._bodyTemplates[name]

    def registerBodyTemplate(self, name, residue, pattern=None):
        """Register a new rigid body template.

        Parameters
        ----------
        name : string
            A name for the rigid body template being registered.
        residue : string
            The name of the residue template to which this body belongs.
        pattern : string=None
            A string containing a regular expression for selecting the atoms in residue which will
            form the body. A list of atom names can be passed as '(name1|name2|...)', for instance.
            Only actual atoms will be considered (i.e. virtual sites will be ignored). If patttern
            is None, then the whole residue will form the rigid body.

        """
        if name in self._bodyTemplates:
            raise ValueError('A rigid body template named %s has already been registered.' % name)
        if residue not in self._templates:
            raise ValueError('Unknown residue %s in rigid body registration.' % residue)
        template = self._templates[residue]
        allAtoms = [a.name for a in template.atoms]
        virtualSites = set(allAtoms[vs.index] for vs in template.virtualSites)
        actualAtoms = list(set(allAtoms) - virtualSites)
        bodyAtoms = [a for a in actualAtoms if re.match(pattern, a)] if pattern else actualAtoms
        if not bodyAtoms:
            raise ValueError('No actual atom in %s matches the provided pattern.' % residue)
        self._bodyTemplates[name] = self._bodyTemplateData(residue, bodyAtoms)

    def _disjointSets(self, sets):
        new = []
        mixed = set()
        for s in sets:
            if s.isdisjoint(mixed):
                new.append(s)
            else:
                i = 0
                while s.isdisjoint(new[i]):
                    i += 1
                new[i] = new[i].union(s)
            mixed = mixed.union(s)
        return new if sum([len(n) for n in new]) == len(mixed) else self._disjointSets(new)


    def resolveBodies(self, topology, merge=None):
        """For each atom in topology, determine the index of the rigid body to which it belongs
        (an index 0 means that the atom does not belong to any body).

        Parameters
        ----------
        topology : Topology
            The Topology for which to resolve the rigid body indices.
        merge : list of list(int)=None
            A nested list of indices of rigid bodies to be merged together after the templates are
            applied. If this is None, then no merging will happen.

        Returns
        -------
        list(int)
            the newly created list of rigid body indices.

        """
        index = [0 for i in range(topology.getNumAtoms())]
        bodies = self._bodyTemplates.values()
        n = 0
        for res in topology.residues():
            for body in filter(lambda x: x.residue == res.name, bodies):
                n += 1
                for atom in filter(lambda x: x.name in body.atoms, res.atoms()):
                    index[atom.index] = n
        if merge is not None:
            if not hasattr(merge, '__iter__'):
                raise ValueError('merge parameter is not a sequence')
            if all(hasattr(a, '__iter__') for a in merge):
                mergeSets = self._disjointSets([set(a) for a in merge])
            else:
                mergeSets = [set(merge)]
            newIndex = [min(a) for a in mergeSets]
            for (i, body) in enumerate(index):
                for (j, s) in enumerate(mergeSets):
                    if body in s:
                        index[i] = newIndex[j]
        return index

    def _intraBodyPairs(self, bodyIndices):
        head = [-1]*max(bodyIndices)
        next = [-1]*len(bodyIndices)
        for i in range(len(bodyIndices)):
            if bodyIndices[i] != 0:
                next[i] = head[bodyIndices[i]-1]
                head[bodyIndices[i]-1] = i
        pairs = list()
        for k in filter(lambda x: x != -1, head):
            i = k
            while (i != -1):
                j = next[i]
                while (j != -1):
                    pairs.append((i, j))
                    j = next[j]
                i = next[i]
        return pairs

    def removeConstraints(self, system, bodyIndices):
        for i in reversed(range(system.getNumConstraints())):
            (atom1, atom2, distance) = system.getConstraintParameters(i)
            if bodyIndices[atom1] != 0 or bodyIndices[atom2] != 0:
                system.removeConstraint(i)

    def removeForces(self, system, bodyIndices):
        def isNonbonded(force):
            return isinstance(force, (mm.NonbondedForce, mm.CustomNonbondedForce))
        forces = [system.getForce(i) for i in range(system.getNumForces())]
        nonbondedForces = [f for f in forces if isNonbonded(f)]
        for (i, j) in self._intraBodyPairs(bodyIndices):
            for force in nonbondedForces:
                force.addException(i, j, 0, 1, 0, replace=True)

    def createSystem(self, topology, **kwargs):
        """Construct an OpenMM System representing a Topology with this force field.
        This extends the original method by creating rigid bodies from previously
        registered body templates.

        Parameters
        ----------
        topology : Topology
            The Topology for which to create a System
        kwargs : Multiple types
            All keyword arguments of the original OpenMM's `ForceField.createSystem` method.
        mergeList : list(list(int))=None
            A nested list of indices of rigid bodies to be merged together after the templates
            are applied. If this is None, then no merging will occur.
        removeConstraints : bool=False
            If true, every constraint involving at least one atom that belongs to a rigid body
            will be silently removed from the system.
        removeForces : bool=False
            If true, all interactions whose all involved atoms belong to the same rigid body
            will be silently removed from the system.

        Returns
        -------
        system : System
            the newly created System
        bodyIndices : list(int)
            a list with the index of each atom's rigid body (index=0 for free atoms)

        """
        mergeList = kwargs.pop('mergeList', None)
        removeForces = kwargs.pop('removeForces', False)
        removeConstraints = kwargs.pop('removeConstraints', False)
        system = super(ForceField, self).createSystem(topology, **kwargs)
        bodyIndices = self.resolveBodies(topology, merge=mergeList)
        if removeConstraints:
            self.removeConstraints(system, bodyIndices)
        if removeForces:
            self.removeForces(system, bodyIndices)
        return (system, bodyIndices)

%}
