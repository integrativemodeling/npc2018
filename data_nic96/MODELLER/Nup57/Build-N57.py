from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb

import fnmatch
import os

log.verbose()
env = environ()

#It begins from here
aln = alignment(env)
aln.append(file='Nup57.ali')
aln.write(file='Nup57_f.ali', alignment_format='PIR')
aln.write(file='Nup57_f.pap', alignment_format='PAP')
aln.id_table(matrix_file='Nup57_f.mat')
aln.check()

class MyModel(automodel):
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['E'],
                             renumber_residues=[287])

a = MyModel(env,
            alnfile='Nup57_f.ali',
            knowns=('5cws_E'),  
            sequence='N57',
            assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 1

a.make()
