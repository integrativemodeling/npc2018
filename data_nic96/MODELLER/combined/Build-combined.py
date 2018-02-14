from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb

log.verbose()
env = environ()
a = alignment(env, file='alignment.ali')
a.write(file='alignment.pap', alignment_format='PAP')

class MyModel(automodel):
    def special_patches(self, aln):
        # Chain A: Nic96; Chain B: Nsp1 (3 segments)
        # Chain C: Nup49; Chain D: Nup57
        self.rename_segments(segment_ids=['A', 'B', 'B', 'B', 'C', 'D'],
                             renumber_residues=[20, 637, 742, 788, 270, 287])

a = MyModel(env,
            alnfile='alignment.ali',
            knowns=('5cwsF', '5cwsC_1', '5cwsC_2', '5cwsC_3', '5cwsD',
                    '5cwsE'),
            sequence='Nic96_complex',
            assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 1

a.make()
