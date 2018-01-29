from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb

import fnmatch
import os

log.verbose()
env = environ()

#It begins from here
aln = alignment(env)
aln.append(file='Nup188C.ali')
aln.write(file='Nup188C_f.ali', alignment_format='PIR')
aln.write(file='Nup188C_f.pap', alignment_format='PAP')
aln.id_table(matrix_file='Nup188C_f.mat')
aln.check()

class MyModel(automodel):
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['A'],
                             renumber_residues=[1273])

a = MyModel(env,
            alnfile='Nup188C_f.ali',
            knowns=('4kf8_A'),  
            sequence='Nup188C',
            assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 1

a.make()
