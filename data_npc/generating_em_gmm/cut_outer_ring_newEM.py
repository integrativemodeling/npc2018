#!/usr/bin/end python
import IMP
import IMP.em
import os, sys
import numpy as np

#dmap = IMP.em.read_map(sys.argv[1])
dmap = IMP.em.read_map("./SJ_cropped_sym8_avg_monomer_final_rotated_adjusted_inner_ring.mrc")
dmap2 = IMP.em.read_map("./SJ_cropped_ynpc_eman_06_01_sym_cleaned_rotated_adjusted90.mrc")
step_size = 1.0
voxel_size = 5.3
threshold = 0.0035

for i in np.arange(0, dmap.get_header().get_nz(), step_size):
    z = (i - 29.5) * voxel_size
    for j in np.arange(0, dmap.get_header().get_ny(), step_size):
        #y = j - 149.5
        y = (j - 58.5) * voxel_size
        for k in np.arange(0, dmap.get_header().get_nx(), step_size):
            #x = k - 149.5
            x = (k + 28.5) * voxel_size
            """
            if ( dmap.get_value(dmap.xyz_ind2voxel(x,y,z)) > 0.0001 ):
                dmap2.set_value(x,y,z, 0.0)
            elif ( abs(x) < abs(y) ):
                dmap.set_value(k,j,i, 0.0)
            elif(dmap.get_value(dmap.xyz_ind2voxel(k,j,i)) >= threshold):
                dmap.set_value(k,j,i, threshold)
            else:
                dmap.set_value(k,j,i, 0.0)
            """

for i in np.arange(0, dmap2.get_header().get_nz(), step_size):
    z = (i - 55.5) * voxel_size
    for j in np.arange(0, dmap2.get_header().get_ny(), step_size):
        #y = j - 149.5
        y = (j - 64.5) * voxel_size
        for k in np.arange(0, dmap2.get_header().get_nx(), step_size):
            #x = k - 149.5
            x = (k + 26.5) * voxel_size

            if ( z > -70 and z < 50 ):
                dmap2.set_value(x,y,z, 0.0)

            if ( (z > 0) and (z < 105) and (x*x + y*y < 465*465) ):
                dmap2.set_value(x,y,z, 0.0)
            if ( (z > 0) and (z < 130) and (x*x + y*y < 340*340) ):
                dmap2.set_value(x,y,z, 0.0)

            if ( (z < 0) and (z > -100) and (x*x + y*y < 440*440) ):
                dmap2.set_value(x,y,z, 0.0)
            if ( (z < 0) and (z > -130) and (x*x + y*y < 330*330) ):
                dmap2.set_value(x,y,z, 0.0)

            if ( (y < -160) and (z > 0) ):
                dmap2.set_value(x,y,z, 0.0)
            if ( (y > 90) and (z < 0) ):
                dmap2.set_value(x,y,z, 0.0)

            """
            #if ( (y > 10) and (y < 150) and (x*x + y*y < 370*370) and (z > 0) ):
            #    dmap2.set_value(x,y,z, 0.0)

            if ( (y < -30) and (x*x + y*y < 360*360) and (z < 0) ):
                dmap2.set_value(x,y,z, 0.0)

            if ( (z < 0) and (z > -130) ):
                dmap2.set_value(x,y,z, 0.0)
            """

dmap3 = dmap2.get_cropped(0.00001)
#IMP.em.write_map(dmap2, "SJ_Modified_Map.mrc")
IMP.em.write_map(dmap3, "SJ_outer_ring_newEM_raw.mrc")
