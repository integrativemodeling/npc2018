import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.display
import IMP.isd
import IMP.isd.gmm_tools
import IMP.pmi1
import IMP.pmi1.tools
import IMP.pmi1.output
import IMP.pmi1.topology
from math import pi, sqrt
import os


###################################
# GMM to MRC
###################################
m = IMP.Model()
ps = []
filename = "./Nup192"
IMP.isd.gmm_tools.decorate_gmm_from_text("%s.txt"%filename, ps, m)
print (len(ps))

#IMP.isd.gmm_tools.write_gmm_to_text(ps, "avg_monomer_final_sj2_ring_6rot.456.txt")
IMP.isd.gmm_tools.write_gmm_to_map(ps, "%s.mrc"%filename, voxel_size=1.0, fast=True)

###################################
# UCSF Chimera command to resample
# vop resample #1 onGrid #2
###################################
