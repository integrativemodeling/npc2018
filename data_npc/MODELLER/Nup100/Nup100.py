from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb

import fnmatch
import os

log.verbose()
env = environ()

#It begins from here
aln = alignment(env)
aln.append(file='Nup100.ali')
aln.write(file='Nup100_f.ali', alignment_format='PIR')
aln.write(file='Nup100_f.pap', alignment_format='PAP')
aln.id_table(matrix_file='Nup100_f.mat')
aln.check()

class MyModel(automodel):
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['A', 'B'],
                             renumber_residues=[803, 803])


a = MyModel(env,
            alnfile='Nup100_f.ali',
            knowns=('3nf5A'),  
            sequence='Nup100',
            assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 1

a.make()
