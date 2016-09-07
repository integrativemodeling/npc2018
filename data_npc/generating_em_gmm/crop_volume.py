#!/usr/bin/end python
import IMP
import IMP.em
import os, sys
import numpy as np

#dmap = IMP.em.read_map(sys.argv[1])
dmap = IMP.em.read_map("./SJ_ynpc_eman_06_01_sym_cleaned_rotated_adjusted90.mrc")
step_size = 1.0
voxel_size = 5.3
threshold = 0.0035

dmap2 = dmap.get_cropped(0.000001)
IMP.em.write_map(dmap2, "./SJ_cropped_ynpc_eman_06_01_sym_cleaned_rotated_adjusted90.mrc")
