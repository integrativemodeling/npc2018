#!/usr/bin/env python

# Simple script to splice Mlp1/Mlp2 beads from final refinement into the
# existing 1-spoke output

import IMP.rmf
import RMF
import os

def find_beads(mh, comp):
    for mol in mh.get_children():
        if mol.get_name() == comp:
            for c in mol.get_children():
                if c.get_name() == 'Beads':
                    return c
    raise ValueError("Could not find %s" % comp)

m = IMP.Model()
r = RMF.open_rmf_file_read_only('cluster0_47-35_1spoke.rmf3')
scaffold = IMP.rmf.create_hierarchies(r, m)
producer = r.get_producer()
assert(len(scaffold) == 1)
scaffold = scaffold[0]

r = RMF.open_rmf_file_read_only('Mlps.rmf3')
mlps = IMP.rmf.create_hierarchies(r, m)
assert(len(mlps) == 1)
mlps = mlps[0]

for comp in ('Mlp1', 'Mlp2'):
    scaffold_beads = find_beads(scaffold, comp)
    while scaffold_beads.get_number_of_children() > 0:
        scaffold_beads.remove_child(0)
    extra_mlp_beads = find_beads(mlps, comp)
    for bead in extra_mlp_beads.get_children():
        scaffold_beads.add_child(bead)

r = RMF.create_rmf_file('out.rmf3')
IMP.rmf.add_hierarchies(r, [scaffold])
IMP.rmf.save_frame(r)
r.set_producer(producer)

os.rename('out.rmf3', 'cluster0_47-35_1spoke.rmf3')
