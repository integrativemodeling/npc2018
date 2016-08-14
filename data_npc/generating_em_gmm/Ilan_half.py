#!/usr/bin/end python
import IMP
import IMP.em
import os, sys
import numpy as np

#dmap = IMP.em.read_map(sys.argv[1])
dmap = IMP.em.read_map("./SJ_cropped_sym8_avg_monomer_final_rotated_adjusted_inner_ring.mrc")
step_size = 1.0
voxel_size = 5.3
threshold = 0.0057

for i in np.arange(0, dmap.get_header().get_nz(), step_size):
    for j in np.arange(0, dmap.get_header().get_ny(), step_size):
        #y = j - 149.5
        y = (j - 58.5) * voxel_size
        for k in np.arange(0, dmap.get_header().get_nx(), step_size):
            #x = k - 149.5
            x = (k + 28.5) * voxel_size
            z = (i - 29.5) * voxel_size
            if ( (y > -30) and (x*x + y*y < 245*245) ):
                dmap.set_value(x,y,z, 0.0)

            if ( (y < -80) and (x*x + y*y < 245*245) ):
                dmap.set_value(x,y,z, 0.0)

            if ( (z < 0) or (abs(y/x) > 0.65) ):
                dmap.set_value(x,y,z, 0.0)

            if ( (x < 280) and (y > 90) ):
                dmap.set_value(x,y,z, 0.0)

            if ( (x < 280) and (y < -130) ):
                dmap.set_value(x,y,z, 0.0)

            if ( (x > 280) and (y < -90) ):
                dmap.set_value(x,y,z, 0.0)

            if ( (x > 390) ):
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
IMP.em.write_map(dmap, "SJ_SamplingBoundary.mrc")
