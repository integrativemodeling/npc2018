import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.display
import IMP.isd
import IMP.isd.gmm_tools
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.output
import IMP.pmi.topology
from math import pi, sqrt
import os

###################################
# Distance_to_point restraint
###################################
if (False):
    m = IMP.Model()
    p = IMP.Particle(m)
    d = IMP.core.XYZR.setup_particle(p)
    d.set_coordinates((417.0, 195.0, 150.0))        # -> 435.099, 150.343, 150  (6.0)   # -> 433.77, 154.135, 150   (5.5)
    #d.set_coordinates((338.0, -170.0, 170.0))      # -> 318.379, -204.399, 170 (6.0)   # -> 320.15, -201.613, 170  (5.5)
    #d.set_coordinates((556.0, -160.0, 110.0))      # -> 536.23, -217.241, 110  (6.0)   # -> 538.105, -212.554, 110 (5.5)

    rotation_angle = -pi * (5.5 / 180.0)
    rot = IMP.algebra.get_rotation_about_axis(IMP.algebra.Vector3D(0., 0., 1.), rotation_angle)
    tr = IMP.algebra.Transformation3D(rot, IMP.algebra.Vector3D(0., 0., 0.))

    print("before : ", p.get_name(), IMP.core.XYZ(p).get_coordinates())
    d = IMP.core.XYZ.setup_particle(p, IMP.core.XYZ(p).get_coordinates())
    IMP.core.transform(d, tr)
    print("after  : ", p.get_name(), IMP.core.XYZ(p).get_coordinates())
    exit(0)

"""
for p in IMP.atom.get_leaves(simo.prot):
    if IMP.core.RigidBodyMember.get_is_setup(p):
        rb = IMP.core.RigidBodyMember(p).get_rigid_body()
        rigid_bodies.add(rb)
    elif IMP.core.XYZR.get_is_setup(p):
        XYZRs.add(p)

for rb in list(rigid_bodies):
    IMP.core.transform(rb, transformation)

for p in list(XYZRs):
    IMP.core.transform(IMP.core.XYZ(p), transformation)
"""

###################################
# GMM rotation
###################################
m = IMP.Model()
ps = []

IMP.isd.gmm_tools.decorate_gmm_from_text("../data_npc/em_gmm_model/avg_monomer_final_sj2_ring.456.txt", ps, m)
print (len(ps))

rotation_angle = -pi * (5.5 / 180.0)
rot = IMP.algebra.get_rotation_about_axis(IMP.algebra.Vector3D(0., 0., 1.), rotation_angle)
tr = IMP.algebra.Transformation3D(rot, IMP.algebra.Vector3D(0., 0., 0.))

for p in ps:
    print("before : ", p.get_name(), IMP.core.XYZ(p).get_coordinates())
    d = IMP.core.XYZ.setup_particle(p, IMP.core.XYZ(p).get_coordinates())
    IMP.core.transform(d, tr)
    print("after  : ", p.get_name(), IMP.core.XYZ(p).get_coordinates())

IMP.isd.gmm_tools.write_gmm_to_text(ps, "avg_monomer_final_sj2_ring_6rot.456.txt")
IMP.isd.gmm_tools.write_gmm_to_map(ps, "avg_monomer_final_sj2_ring_6rot.456.mrc", voxel_size=5.6, fast=True)

###################################
# UCSF Chimera command to resample
# vop resample #0 onGrid #1
###################################
