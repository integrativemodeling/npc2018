from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb

import fnmatch
import os

log.verbose()
env = environ()

#It begins from here
aln = alignment(env)
aln.append(file='Nic96N_align.ali')
aln.write(file='Nic96N_f.ali', alignment_format='PIR')
aln.write(file='Nic96N_f.pap', alignment_format='PAP')
aln.id_table(matrix_file='Nic96N_f.mat')
aln.check()

class MyModel(automodel):
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['F'],
                             renumber_residues=[20])

a = MyModel(env,
            alnfile='Nic96N_f.ali',
            knowns=('5cws_F'),  
            sequence='Nic96_N',
            assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 1

a.make()
