#!/usr/bin/env python
#####################################################
# Last Update: June 30th, 2016 by Seung Joong Kim
# Andrej Sali group, University of California San Francisco (UCSF)
#####################################################
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi1.restraints.crosslinking
import IMP.pmi1.restraints.stereochemistry
import IMP.pmi1.restraints.em
import IMP.pmi1.restraints.basic
import IMP.pmi1.restraints.proteomics
import IMP.pmi1.representation
import IMP.pmi1.macros
import IMP.pmi1.restraints
import IMP.pmi1.tools
import IMP.pmi1.output
import IMP.pmi1.samplers
#import IMP.pmi1.topology
#import IMP.pmi1.dof
import IMP.npc
import IMP.pmi1.restraints.npc
import random
import os
import math

#####################################################
# Parsing parameter inputs
#####################################################
import argparse

parser = argparse.ArgumentParser(description='Performing the INITIAL/REFINEMENT Monte Carlo job, with crosslinks and selected/ALL domain mapping data. Example of usage: setup_environment.sh python ./sj_SEA_XLDM.py -f models_1877.rmf -n 0')
parser.add_argument('-copy', action="store", dest="ncopy", help="copy numbers (stoichiometry) for SEA4 and Seh1" )
parser.add_argument('-sym', action="store", dest="symmetry", help="symmetry option for SEA4 and Seh1" )
parser.add_argument('-rmf', action="store", dest="rmf_input", help="rmf file name to continue" )
parser.add_argument('-rmf_n', action="store", dest="rmf_frame_number", help="rmf frame number to continue" )
parser.add_argument('-em2d', action="store", dest="em2d_input", help="em2d image file name to read" )
parser.add_argument('-r', action="store", dest="nrepeats", help="number of Monte Carlo cycles" )
parser.add_argument('-x', action="store", dest="XL_input", help="Cross-links file name to read" )
parser.add_argument('-out', action="store", dest="folder_output", help="folder name for output" )
parser.add_argument('-o', action="store", dest="rmf_output", help="rmf file name for output" )
parser.add_argument('-s', action="store", dest="stat_output", help="stat file name for output" )
parser.add_argument('-REFINE', action="store", dest="refinement", help="refinement True or False" )
parser.add_argument('-weight', action="store", dest="weight", help="weight for the EM 2D restraint" )
parser.add_argument('-res_cry', action="store", dest="res_cry", help="resolution of the crystal structures" )
parser.add_argument('-res_hom', action="store", dest="res_hom", help="resolution of the comparative (homology) models" )
parser.add_argument('-res_ev', action="store", dest="res_ev", help="resolution of the excluded volume restraints" )
parser.add_argument('-res_compo', action="store", dest="res_compo", help="resolution of the composite restraints" )
parser.add_argument('-draw_hierarchy', action="store", dest="draw_hierarchy", help="draw hierarchy" )
inputs = parser.parse_args()

#####################################################
# Setting up the input parameters
#####################################################
if (inputs.ncopy is None) :
    inputs.ncopy = "2"
if (inputs.symmetry == "True") or (inputs.symmetry == "true") or (inputs.symmetry == "Yes") or (inputs.symmetry == "yes") :
    inputs.symmetry = True
else:
    inputs.symmetry = False

if (inputs.rmf_input is not None) :
    f=open(inputs.rmf_input,"r")
    f.close()
if (inputs.rmf_frame_number is None) :
    inputs.rmf_frame_number = 0
if (inputs.em2d_input is not None) :
    f=open(inputs.em2d_input,"r")
    f.close()

if (inputs.XL_input is None) :
    inputs.XL_input = "./p85_xl1.txt"
else :
    f=open(inputs.XL_input,"r")
    f.close()
if (inputs.nrepeats is None) :
    inputs.nrepeats = 1000
if (inputs.folder_output is None) :
    inputs.folder_output = "output"
if (inputs.rmf_output is None) :
    inputs.rmf_output = "models.rmf"
if (inputs.stat_output is None) :
    inputs.stat_output = "stat.dat"
if (inputs.refinement == "True") or (inputs.refinement == "true") or (inputs.refinement == "Yes") or (inputs.refinement == "yes") :
    inputs.refinement = True
else:
    inputs.refinement = False
if (inputs.weight is None) :
    inputs.weight = 10000.0

if (inputs.res_cry is None) :
    inputs.res_cry = 1.0
if (inputs.res_hom is None) :
    inputs.res_hom = 5.0
if (inputs.res_ev is None) :
    inputs.res_ev = 10.0
if (inputs.res_compo is None) :
    inputs.res_compo = 100.0
if (inputs.draw_hierarchy == "True") or (inputs.draw_hierarchy == "true") or (inputs.draw_hierarchy == "Yes") or (inputs.draw_hierarchy == "yes") :
    inputs.draw_hierarchy = True
else:
    inputs.draw_hierarchy = False
print inputs


#####################################################
# setting up topology and parameters
#####################################################
import RMF
#rmf_file_name = "../data_npc/WholeNPC_rmfs/47-35_1spoke.rmf3"
rmf_file_name = "../data_npc/WholeNPC_rmfs/47-35_truncated.rmf3"
#rmf_file_name = "../data_npc/WholeNPC_rmfs/70-193_truncated.rmf3"

hier_raw = { "Nup84complex" : ['Nup84', 'Nup85', 'Nup120', 'Nup133', 'Nup145c', 'Seh1', 'Sec13'],
             "Nup82complex" : ['Dyn2', 'Nup82', 'Nup159', 'Nsp1.1', 'Nsp1.2'],
             "Nic96complex" : ['Nic96.1', 'Nsp1.3', 'Nup49.1', 'Nup57.1', 'Nic96.2', 'Nsp1.4', 'Nup49.2', 'Nup57.2'],
             "Inner_ring"   : ['Nup157', 'Nup170', 'Nup188', 'Nup192'],
             "Membrane_ring": ['Nup53', 'Nup59', 'Ndc1', 'Pom34', 'Pom152'],
             "Cytoplasm"    : ['Nup116', 'Nup100', 'Nup42', 'Gle1'],
             "Nucleoplasm"  : ['Nup145.1', 'Nup145.2', 'Nup60', 'Nup1', 'Nup1@2', 'Nup1@3'],
             "Mlps"         : ['Mlp1', 'Mlp2'],
            }
for hier_name, hier_dict in hier_raw.items():
    print (hier_name, hier_dict)

    m = IMP.Model()
    rh = RMF.open_rmf_file_read_only(rmf_file_name)
    try:
        h = IMP.rmf.create_hierarchies(rh, m)
    except:
        continue
    ps = IMP.atom.get_leaves(h[0]); print (len(ps))
    for component in h[0].get_children():
        flag = 0
        for hier in hier_dict:
            if (hier == "Nup1" and hier != component.get_name()):
                continue
            if (hier in component.get_name()):
                flag = 1
                print ("component.get_name() = ", component.get_name())

        if (flag == 0):
            IMP.atom.destroy(component)

    ps = IMP.atom.get_leaves(h[0]); print (len(ps))
    rh = RMF.create_rmf_file(hier_name + ".rmf3")
    IMP.rmf.add_hierarchies(rh, h)
    IMP.rmf.save_frame(rh)
    del rh

    m = IMP.Model()
    prots = IMP.pmi1.analysis.get_hiers_from_rmf(m, 0, hier_name + ".rmf3")
    if not prots:
        raise ValueError("cannot read hiearchy from rmf")
    prot=prots[0]

    o = IMP.pmi1.output.Output()
    o.init_pdb(hier_name + ".pdb", prot)
    o.write_pdb(hier_name + ".pdb")
    del o
    #IMP.atom.write_pdb(prot, "bb.pdb")

    ps = IMP.atom.get_leaves(prot)
    image_resolution = 30.0
    map = IMP.em.SampledDensityMap(ps, image_resolution, 10.0)
    IMP.em.write_map(map, hier_name + ".mrc")
