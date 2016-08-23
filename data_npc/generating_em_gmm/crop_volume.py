#!/usr/bin/end python
import IMP
import IMP.em
import os, sys
import numpy as np

#dmap = IMP.em.read_map(sys.argv[1])
dmap = IMP.em.read_map("./SJ_outer_ring.mrc")
step_size = 1.0
threshold = 0.0013

dmap2 = dmap.get_cropped(0.000001)
IMP.em.write_map(dmap2, "./aSJ_outer_ring.mrc")
