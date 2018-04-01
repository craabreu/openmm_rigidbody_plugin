%pythoncode %{
import simtk.openmm.app as app
import re


def _coalescedSets(sets):
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
    return new if sum([len(n) for n in new]) == len(mixed) else _coalescedSets(new)


class ForceField(app.ForceField):
    def __init__(self,*args):
        super(ForceField, self).__init__(*args)
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
                mergeSets = _coalescedSets([set(a) for a in merge])
            else:
                mergeSets = set(set(merge))
            newIndex = [min(a) for a in mergeSets]
            for (i, body) in enumerate(index):
                for (j, s) in enumerate(mergeSets):
                    if body in s:
                        index[i] = newIndex[j]
        return index
%}
