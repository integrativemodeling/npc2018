from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb

import fnmatch
import os

log.verbose()
env = environ()

#It begins from here
aln = alignment(env)
aln.append(file='Nup59.ali')
aln.write(file='Nup59_f.ali', alignment_format='PIR')
aln.write(file='Nup59_f.pap', alignment_format='PAP')
aln.id_table(matrix_file='Nup59_f.mat')
aln.check()

class MyModel(automodel):
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['A', 'B'],
                             renumber_residues=[266, 266])


a = MyModel(env,
            alnfile='Nup59_f.ali',
            knowns=('Nup53'),  
            sequence='Nup59',
            assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 1

a.make()
