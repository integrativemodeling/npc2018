from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb

import fnmatch
import os

log.verbose()
env = environ()

#It begins from here
aln = alignment(env)
aln.append(file='Gle1N.ali')
aln.write(file='Gle1N_f.ali', alignment_format='PIR')
aln.write(file='Gle1N_f.pap', alignment_format='PAP')
aln.id_table(matrix_file='Gle1N_f.mat')
aln.check()

class MyModel(automodel):
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['A'],
                             renumber_residues=[121])

    def special_restraints(self, aln):
        rsr = self.restraints
        at = self.atoms
        
        rsr.add(secondary_structure.alpha(self.residue_range('121:A', '140:A')))


a = MyModel(env,
            alnfile='Gle1N_f.ali',
            knowns=('4wijA'),  
            sequence='Gle1N',
            assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 1

a.make()
