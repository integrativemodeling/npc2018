#!/usr/bin/end python
import IMP
import IMP.em
import os, sys
import numpy as np


###################################
# UCSF Chimera command to resample
# vop resample #1 onGrid #0
###################################

#dmap = IMP.em.read_map(sys.argv[1])
dmap = IMP.em.read_map("./SJ_ynpc_eman_06_01_sym_cleaned_rotated.mrc")
step_size = 1.0
voxel_size = 5.3
threshold = 0.0035

for i in np.arange(0, dmap.get_header().get_nz(), step_size):
    z = (i - 149.5) * voxel_size
    for j in np.arange(0, dmap.get_header().get_ny(), step_size):
        y = (j - 149.5) * voxel_size
        for k in np.arange(0, dmap.get_header().get_nx(), step_size):
            x = (k - 149.5) * voxel_size
            if ( (x < 0) or (x*x + y*y < 36*36*voxel_size*voxel_size) ):
                dmap.set_value(x,y,z, 0.0)
            """
            elif ( abs(x) < abs(y) ):
                dmap.set_value(k,j,i, 0.0)
            elif(dmap.get_value(dmap.xyz_ind2voxel(k,j,i)) >= threshold):
                dmap.set_value(k,j,i, threshold)
            else:
                dmap.set_value(k,j,i, 0.0)
            """
#dmap2 = dmap.get_cropped(0.00001)
#IMP.em.write_map(dmap2, "SJ_Modified_Map.mrc")
IMP.em.write_map(dmap, "./SJ_ynpc_eman_06_01_sym_cleaned_rotated_half.mrc")
