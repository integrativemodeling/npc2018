# Modify the model read from an old RMF file so that it looks like
# one read from a modern RMF file. (Modern RMF files contain provenance
# information that is needed to generate an mmCIF file.)

from __future__ import print_function
import IMP.pmi.topology
import IMP.mmcif
import os
import RMF

class RMFFrame(IMP.mmcif.RMFFrame):
    def __init__(self, *args, **keys):
        super(RMFFrame, self).__init__(*args, **keys)
        self._seqs = _SequenceMapper()
        self._producer = None

    def create(self, model):
        if not self._producer:
            # Get producer info from file
            rmf = RMF.open_rmf_file_read_only(self.filename)
            prod = rmf.get_producer().split()
            assert(len(prod) == 2)
            assert(prod[0] == 'IMP')
            self._producer = prod[1]
        hiers, restraints = super(RMFFrame, self).create(model)
        self._add_chains_sequences(hiers)
        return hiers, restraints

    def _add_chains_sequences(self, hiers):
        """Add Chain objects with primary sequence to the model"""
        assert(len(hiers) == 1)
        mols = [IMP.atom.Molecule(m)
                for m in IMP.atom.get_by_type(hiers[0], IMP.atom.MOLECULE_TYPE)]
        for m in mols:
            # set up Molecule as a Chain as well
            chain_h = IMP.atom.Chain.setup_particle(m, 'X')
            seq = self._seqs[m.get_name()]
            chain_h.set_sequence(seq)
            chain_h.set_chain_type(IMP.atom.Protein)
            #print("Molecule %s, added chain with sequence %s"
            #      % (m.get_name(), seq))


class _SequenceMapper(object):
    """Map molecule names to primary sequence"""

    # All directories containing FASTA files
    paths = ["../../data_nup82", "../../data_nic96",
             "../../data_npc"]

    def __init__(self):
        self._seqmap = {}

    def __getitem__(self, molname):
        seqname = molname.split('@')[0].split('.')[0]
        if seqname not in self._seqmap:
            self._seqmap[seqname] = self._read_fasta(seqname)
        return self._seqmap[seqname]

    def _read_fasta(self, seqname):
        for p in self.paths:
            fname = os.path.join(p, "protein_fasta.%s.txt" % seqname)
            if os.path.exists(fname):
                s = IMP.pmi.topology.Sequences(fname)
                assert(len(s) == 1)
                return s[0]
        raise ValueError("Could not find sequence for %s" % seqname)
