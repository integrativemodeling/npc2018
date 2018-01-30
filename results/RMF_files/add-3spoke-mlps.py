#!/usr/bin/env python

# Simple script to splice Mlp1/Mlp2 beads from final refinement into the
# existing 3-spoke output

import IMP.rmf
import IMP.atom
import RMF
import os
import math

def find_beads(mh, comp):
    for mol in mh.get_children():
        if mol.get_name() == comp:
            for c in mol.get_children():
                if c.get_name() == 'Beads':
                    return c
    raise ValueError("Could not find %s" % comp)

m = IMP.Model()
r = RMF.open_rmf_file_read_only('cluster0_47-35_3spokes.rmf3')
scaffold = IMP.rmf.create_hierarchies(r, m)
producer = r.get_producer()
assert(len(scaffold) == 1)
scaffold = scaffold[0]

r = RMF.open_rmf_file_read_only('Mlps.rmf3')
mlps = IMP.rmf.create_hierarchies(r, m)
assert(len(mlps) == 1)
mlps = mlps[0]

def rotate(hier, angle):
    rot = IMP.algebra.get_rotation_about_axis(IMP.algebra.Vector3D(0., 0., 1.),
                                              angle)
    tr = IMP.algebra.Transformation3D(rot)
    # Find all rigid bodies
    rbs = {}
    for nup in hier.get_children():
        beads, = nup.get_children()
        for bead in beads.get_children():
            assert(IMP.core.RigidBodyMember.get_is_setup(bead))
            rb = IMP.core.RigidBodyMember(bead).get_rigid_body()
            rbs[rb] = None
    # Rotate each rigid body
    for rb in rbs:
        IMP.core.transform(rb, tr)

# Make symmetry units
r = RMF.open_rmf_file_read_only('Mlps.rmf3')
mlps2, = IMP.rmf.create_hierarchies(r, m)
rotate(mlps2, 0.25 * math.pi)

r = RMF.open_rmf_file_read_only('Mlps.rmf3')
mlps3, = IMP.rmf.create_hierarchies(r, m)
rotate(mlps3, -0.25 * math.pi)

for comp, in_mlp in (('Mlp1', mlps), ('Mlp2', mlps),
                     ('Mlp1@2', mlps2), ('Mlp2@2', mlps2),
                     ('Mlp1@3', mlps3), ('Mlp2@3', mlps3)):
    orig_comp = comp.split('@')[0]
    scaffold_beads = find_beads(scaffold, comp)
    while scaffold_beads.get_number_of_children() > 0:
        scaffold_beads.remove_child(0)
    extra_mlp_beads = find_beads(in_mlp, orig_comp)
    for bead in extra_mlp_beads.get_children():
        if orig_comp != comp:
            bead.set_name(bead.get_name().replace(orig_comp, comp))
        scaffold_beads.add_child(bead)

r = RMF.create_rmf_file('out.rmf3')
IMP.rmf.add_hierarchies(r, [scaffold])
IMP.rmf.save_frame(r)
r.set_producer(producer)
del r

os.rename('out.rmf3', 'cluster0_47-35_3spokes.rmf3')
