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
# Distance_to_point restraint
###################################
if (False):
    m = IMP.Model()
    p = IMP.Particle(m)
    d = IMP.core.XYZR.setup_particle(p)
    #d.set_coordinates((487.7892857, 99.75357143, 159.1892857))  # 21.5 degree rotation  -> 417.288, 271.588, 159.189
    #d.set_coordinates((271.4357143, -165.0571429, 218.1517857))  # 21.5 degree rotation  -> 313.042, -54.0905, 218.152
    d.set_coordinates((462.9928571, -237.8375, 191.4625))       # 21.5 degree rotation  -> 517.944, -51.6007, 191.462

    #rotation_angle = -pi * (5.5 / 180.0)
    rotation_angle = pi * (21.5 / 180.0)
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

#IMP.isd.gmm_tools.decorate_gmm_from_text("../data_npc/em_gmm_model/avg_monomer_final_sj2_ring.456.txt", ps, m)
IMP.isd.gmm_tools.decorate_gmm_from_text("./SJ_cropped_sym8_avg_monomer_final_rotated_adjusted90.gmm.100.txt", ps, m)
print (len(ps))

#rotation_angle = -pi * (5.5 / 180.0)
#rot = IMP.algebra.get_rotation_about_axis(IMP.algebra.Vector3D(0., 0., 1.), rotation_angle)
tr = IMP.algebra.Transformation3D(IMP.algebra.Vector3D(0., 0., 23.83))

for p in ps:
    print("before : ", p.get_name(), IMP.core.XYZ(p).get_coordinates())
    d = IMP.core.XYZ.setup_particle(p, IMP.core.XYZ(p).get_coordinates())
    IMP.core.transform(d, tr)
    print("after  : ", p.get_name(), IMP.core.XYZ(p).get_coordinates())

IMP.isd.gmm_tools.write_gmm_to_text(ps, "2SJ_cropped_sym8_avg_monomer_final_rotated_adjusted90.gmm.100.txt")
IMP.isd.gmm_tools.write_gmm_to_map(ps, "2SJ_cropped_sym8_avg_monomer_final_rotated_adjusted90.gmm.100.mrc", voxel_size=6.0, fast=True)


dmap = IMP.em.read_map("./2SJ_cropped_sym8_avg_monomer_final_rotated_adjusted90.gmm.100.mrc")
dmap2 = dmap.get_cropped(1.0e-9)
IMP.em.write_map(dmap2, "./2SJ_cropped_sym8_avg_monomer_final_rotated_adjusted90.gmm.100.mrc")


###################################
# UCSF Chimera command to resample
# vop resample #0 onGrid #1
###################################
