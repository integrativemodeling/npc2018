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

import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.proteomics
import IMP.pmi.representation
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.tools
import IMP.pmi.output
import IMP.pmi.samplers
#import IMP.pmi.topology
#import IMP.pmi.dof
import IMP.npc
import IMP.npc.npc_restraints
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
m = IMP.Model()
#s = IMP.pmi.topology.System(m)
#st = s.create_state()
simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)

try:
    from mpi4py import MPI
except ImportError:
    MPI = None

if MPI:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
else:
    rank = 0
print "rank = ", rank

# rigid body movement params
rbmaxtrans = 3.00
rbmaxrot = 0.03

# flexible bead movement
fbmaxtrans = 3.00
outputobjects = []
sampleobjects = []

res_cry = int(inputs.res_cry)
res_ev = int(inputs.res_ev)
res_conn = int(inputs.res_cry)
beadsize = 10
beadsize25 = 25
beadsize50 = 50
beadsize100 = 100
em2d_weight = float(inputs.weight)

n82   = "../data_nup82/"
f_n82 = "../data_nup82/protein_fasta."
n96   = "../data_nic96/"
f_n96 = "../data_nic96/protein_fasta."
npc   = "../data_npc/"
f_npc = "../data_npc/protein_fasta."
gmm_f = npc + "em_gmm_model/"

n84_fastafile  = '../data_nup84_2016/protein_fasta.Nup84.txt'
n85_fastafile  = '../data_nup84_2016/protein_fasta.Nup85.txt'
n120_fastafile = '../data_nup84_2016/protein_fasta.Nup120.txt'
n133_fastafile = '../data_nup84_2016/protein_fasta.Nup133.txt'
n145c_fastafile= '../data_nup84_2016/protein_fasta.Nup145c.txt'
seh1_fastafile = '../data_nup84_2016/protein_fasta.Seh1.txt'
sec13_fastafile= '../data_nup84_2016/protein_fasta.Sec13.txt'

n84_pdb        = n82 + 'rmfs/B_8_1-95_noBEA.pdb'    #n84_pdb        = n82 + 'rmfs/B_8_1-95r.pdb'    # for GMM generation
n82_pdb        = n82 + 'rmfs/B_8_1-95_noBEA.pdb'    #n82_pdb        = n82 + 'rmfs/B_8_1-95r.pdb'    # for GMM generation
n96_pdb        = n96 + 'rmfs/Cluster2_noBEA.pdb'    #n96_pdb        = n96 + 'rmfs/cluster2r.pdb'    # for GMM generation
n961_pdb       = npc + 'Inner_ring_rmfs/Nic96.1.pdb'    #n961_pdb       = npc + 'Inner_ring_rmfs/Nic96.1.rebuilt.pdb'   # for GMM generation
n962_pdb       = npc + 'Inner_ring_rmfs/Nic96.2.pdb'    #n962_pdb       = npc + 'Inner_ring_rmfs/Nic96.2.rebuilt.pdb'   # for GMM generation

n157N_pdbfile  = npc + "Nup157_0_noBEA.pdb"         #n157N_pdbfile  = npc + "Nup157_4mhc_A_88_892ca.pdb"
n157C_pdbfile  = npc + "Nup157_0_noBEA.pdb"         #n157C_pdbfile  = npc + "Nup157_3i5p_A_900_1390ca.pdb"
n170N_pdbfile  = npc + "Nup170_0_noBEA.pdb"         #n170N_pdbfile  = npc + "Nup170_4mhc_A_98_992ca.pdb"
n170C_pdbfile  = npc + "Nup170_0_noBEA.pdb"         #n170C_pdbfile  = npc + "Nup170_3i5p_A_1000_1502ca.pdb"
n188_pdbfile   = npc + "Nup188_2_noBEA.pdb"         #n188_pdbfile   = npc + "Nup188_12_1652ca.pdb"
n192_pdbfile   = npc + "Nup192_4_noBEA.pdb"         #n192_pdbfile   = npc + "Nup192_1_1683_v2_new_ca.pdb"

n53_pdbfile    = npc + "Nup53_9_noBEA.pdb"          #n53_pdbfile    = npc + "Nup53_dimer_ca.pdb"
n59_pdbfile    = npc + "Nup59_0_noBEA.pdb"          #n59_pdbfile    = npc + "Nup59_dimer_ca.pdb"
pom152_pdb     = npc + "Pom152_noBEA.pdb"           #pom152_pdb     = npc + "Pom152r.pdb"   # for GMM generation
n100_pdbfile   = npc + "Nup100_3nf5_AB_816_958ca.pdb"
n116_pdbfile   = npc + "Outer_ring_rmfs/Nup116ca.pdb"
n145N_pdbfile  = npc + "Nup145N_3kep_AB_459_605r_ca.pdb"

#Dbp5_pdbfile   = npc + "Dbp5_3rrm_A_91_482.pdb"    # Not using it for now
Gle1N_pdbfile  = npc + "Gle1_4wij_A_121_239ca.pdb"
Gle1C_pdbfile  = npc + "Gle1_3rrm_B_244_538ca.pdb"
Gle2_pdbfile   = npc + "Gle2_3mmy_A_4_362ca.pdb"

#####################################################
# Parameters for Debugging
#####################################################
is_n84 = True
is_n82 = True
is_nic96 = True
is_inner_ring = True
is_membrane = True
is_cytoplasm = True
is_nucleoplasm = True
is_basket = True
is_FG = False

use_neighboring_spokes = True
#Stopwatch_None_delta_seconds  ~22   (1 spoke for OR / IR + 3 spokes for others, 3.0G memory) with XL
#Stopwatch_None_delta_seconds  ~25   (1 spoke for OR / IR + 3 spokes for others, 3.0G memory) with XL + EM
#Stopwatch_None_delta_seconds  ~65   (1 spoke for OR / IR + 3 spokes for others, 5.0G memory) with XL + EM + EV
#Stopwatch_None_delta_seconds  ~150   (3 spokes for all, ~8.0G memory) with XL + EM + EV
use_shuffle = False
use_ExcludedVolume = True
use_Immuno_EM = False
use_FG_anchor = False
use_sampling_boundary = True
use_MembraneExclusion = True
use_XL = True
use_EM3D = True

#####################################################
# REPRESENTATION
#####################################################
# comp_name, hier_name, color, fasta_file, fasta_id, pdb_name, chain_id, res_range, read_em_files, bead_size, rb, super_rb, em_num_components, em_txt_file_name, em_mrc_file_name, chain_of_super_rb, keep_gaussian_on_flexible_beads
domains = []
if (use_EM3D):  gmm = True
else:           gmm = None
if (use_neighboring_spokes):
    clones_range_A = range(2,4)+range(11,14)
    clones_range_B = range(2,4)
else:
    clones_range_A = range(11,12)
    clones_range_B = []
##########################
# Nup84 complex
##########################
if (is_n84):
    #n84_rb = None;    n133_rb = None;  n120_rb = None
    n84_rb = 84;    n133_rb = 133;  n120_rb = 120;  n120c_rb = 84;    n85_rb = 84
    domains.append(('Nup84',  "Nup84",    0.0,  n84_fastafile,   "Nup84",   n84_pdb, "K",  (  1, 726,0),  gmm,  beadsize,  n84_rb, [n84_rb], 3,  " ",   " ",  None, False))
    domains.append(('Nup85',  "Nup85_1",  0.2,  n85_fastafile,   "Nup85",   n84_pdb, "L",  (  1, 492,0),  gmm,  beadsize,  n85_rb, [n84_rb], 3,  " ",   " ",  None, False))
    domains.append(('Nup85',  "Nup85_2",  0.25, n85_fastafile,   "Nup85",   n84_pdb, "L",  (493, 744,0),  gmm,  beadsize,  n84_rb, [n84_rb], 2,  " ",   " ",  None, False))
    domains.append(('Nup120', "Nup120_1", 0.35, n120_fastafile,  "Nup120",  n84_pdb, "M",  (  1, 714,0),  gmm,  beadsize,  n120_rb,[n84_rb], 3,  " ",   " ",  None, False))
    domains.append(('Nup120', "Nup120_2", 0.4,  n120_fastafile,  "Nup120",  n84_pdb, "M",  (715,1037,0),  gmm,  beadsize,  n120c_rb,[n84_rb],2,  " ",   " ",  None, False))
    domains.append(('Nup133', "Nup133_1", 0.5,  n133_fastafile,  "Nup133",  n84_pdb, "N",  (  1, 489,0),  gmm,  beadsize,  n133_rb,[n84_rb], 2,  " ",   " ",  None, False))
    domains.append(('Nup133', "Nup133_2", 0.55, n133_fastafile,  "Nup133",  n84_pdb, "N",  (490,1157,0),  gmm,  beadsize,  n84_rb, [n84_rb], 3,  " ",   " ",  None, False))
    domains.append(('Nup145c',"Nup145c_1",0.65, n145c_fastafile, "Nup145c", n84_pdb, "O",  (  1, 125,0),  gmm,  beadsize,  n84_rb, [n84_rb], 1,  " ",   " ",  None, False))
    domains.append(('Nup145c',"Nup145c_2",0.7,  n145c_fastafile, "Nup145c", n84_pdb, "O",  (126, 712,0),  gmm,  beadsize,  n84_rb, [n84_rb], 3,  " ",   " ",  None, False))
    domains.append(('Seh1',   "Seh1",     0.8,  seh1_fastafile,  "Seh1",    n84_pdb, "P",  (  1, 349,0),  gmm,  beadsize,  n85_rb, [n84_rb], 2,  " ",   " ",  None, False))
    domains.append(('Sec13',  "Sec13",    0.95, sec13_fastafile, "Sec13",   n84_pdb, "Q",  (  1, 297,0),  gmm,  beadsize,  n84_rb, [n84_rb], 2,  " ",   " ",  None, False))
    for i in clones_range_A:
        if (i==11): gmm_c = gmm
        else:       gmm_c = None
        domains.append(('Nup84@%d'%i,  "Nup84@%d"%i,    0.0, n84_fastafile,   "Nup84",   n84_pdb, "K",  (  1, 726,0),  gmm_c,  beadsize,   None,   None,   3,  gmm_f+"Nup84.txt",     gmm_f+"Nup84.mrc",     None, False))
        domains.append(('Nup85@%d'%i,  "Nup85_1@%d"%i,  0.2, n85_fastafile,   "Nup85",   n84_pdb, "L",  (  1, 492,0),  gmm_c,  beadsize,   None,   None,   3,  gmm_f+"Nup85_1.txt",   gmm_f+"Nup85_1.mrc",   None, False))
        domains.append(('Nup85@%d'%i,  "Nup85_2@%d"%i,  0.25,n85_fastafile,   "Nup85",   n84_pdb, "L",  (493, 744,0),  gmm_c,  beadsize,   None,   None,   2,  gmm_f+"Nup85_2.txt",   gmm_f+"Nup85_2.mrc",   None, False))
        domains.append(('Nup120@%d'%i, "Nup120_1@%d"%i, 0.35,n120_fastafile,  "Nup120",  n84_pdb, "M",  (  1, 714,0),  gmm_c,  beadsize,   None,   None,   3,  gmm_f+"Nup120_1.txt",  gmm_f+"Nup120_1.mrc",  None, False))
        domains.append(('Nup120@%d'%i, "Nup120_2@%d"%i, 0.4, n120_fastafile,  "Nup120",  n84_pdb, "M",  (715,1037,0),  gmm_c,  beadsize,   None,   None,   2,  gmm_f+"Nup120_2.txt",  gmm_f+"Nup120_2.mrc",  None, False))
        domains.append(('Nup133@%d'%i, "Nup133_1@%d"%i, 0.5, n133_fastafile,  "Nup133",  n84_pdb, "N",  (  1, 489,0),  gmm_c,  beadsize,   None,   None,   2,  gmm_f+"Nup133_1.txt",  gmm_f+"Nup133_1.mrc",  None, False))
        domains.append(('Nup133@%d'%i, "Nup133_2@%d"%i, 0.55,n133_fastafile,  "Nup133",  n84_pdb, "N",  (490,1157,0),  gmm_c,  beadsize,   None,   None,   3,  gmm_f+"Nup133_2.txt",  gmm_f+"Nup133_2.mrc",  None, False))
        domains.append(('Nup145c@%d'%i,"Nup145c_1@%d"%i,0.65,n145c_fastafile, "Nup145c", n84_pdb, "O",  (  1, 125,0),  gmm_c,  beadsize,   None,   None,   1,  gmm_f+"Nup145c_1.txt", gmm_f+"Nup145c_1.mrc", None, False))
        domains.append(('Nup145c@%d'%i,"Nup145c_2@%d"%i,0.7, n145c_fastafile, "Nup145c", n84_pdb, "O",  (126, 712,0),  gmm_c,  beadsize,   None,   None,   3,  gmm_f+"Nup145c_2.txt", gmm_f+"Nup145c_2.mrc", None, False))
        domains.append(('Seh1@%d'%i,   "Seh1@%d"%i,     0.8, seh1_fastafile,  "Seh1",    n84_pdb, "P",  (  1, 349,0),  gmm_c,  beadsize,   None,   None,   2,  gmm_f+"Seh1.txt",      gmm_f+"Seh1.mrc",      None, False))
        domains.append(('Sec13@%d'%i,  "Sec13@%d"%i,    0.95,sec13_fastafile, "Sec13",   n84_pdb, "Q",  (  1, 297,0),  gmm_c,  beadsize,   None,   None,   2,  gmm_f+"Sec13.txt",     gmm_f+"Sec13.mrc",     None, False))

##########################
# Nup82 complex
##########################
if (is_n82):
    n82_rb = 82
    domains.append(("Dyn2.1",  "Dyn2.1",      0.48,  f_n82+"Dyn2.txt",   "Dyn2",   n82_pdb,  "A", (   1,  92,0),  gmm,   beadsize,   n82_rb, [n82_rb], 1,  " ",   " ",  None, False))
    domains.append(("Dyn2.2",  "Dyn2.2",      0.65,  f_n82+"Dyn2.txt",   "Dyn2",   n82_pdb,  "B", (   1,  92,0),  gmm,   beadsize,   n82_rb, [n82_rb], 1,  " ",   " ",  None, False))
    domains.append(("Nup82.1", "Nup82.1_1",   0.0,   f_n82+"Nup82.txt",  "Nup82",  n82_pdb,  "C", (   1, 452,0),  gmm,   beadsize,   n82_rb, [n82_rb], 2,  " ",   " ",  None, False))
    domains.append(("Nup82.1", "Nup82.1_2",   0.0,   f_n82+"Nup82.txt",  "Nup82",  n82_pdb,  "C", ( 453, 713,0),  gmm,   beadsize,   n82_rb, [n82_rb], 4,  " ",   " ",  None, False))
    domains.append(("Nup82.2", "Nup82.2_1",   0.15,  f_n82+"Nup82.txt",  "Nup82",  n82_pdb,  "D", (   1, 452,0),  gmm,   beadsize,   n82_rb, [n82_rb], 2,  " ",   " ",  None, False))
    domains.append(("Nup82.2", "Nup82.2_2",   0.15,  f_n82+"Nup82.txt",  "Nup82",  n82_pdb,  "D", ( 453, 713,0),  gmm,   beadsize,   n82_rb, [n82_rb], 4,  " ",   " ",  None, False))
    if (is_FG):
        domains.append(("Nup159.1","Nup159.1_1",  1.0,   f_n82+"Nup159.txt", "Nup159", n82_pdb,  "E", (   1, 381,0),  None,  beadsize,     1159, [n82_rb], 0,  None,  None, None))
        domains.append(("Nup159.1","Nup159.1_10", 1.0,   f_n82+"Nup159.txt", "Nup159", "BEADS",  " ", ( 382,1116,0),  None,  beadsize100,  1159, [n82_rb], 0,  None,  None, None))
        domains.append(("Nup159.2","Nup159.2_1",  0.9,   f_n82+"Nup159.txt", "Nup159", n82_pdb,  "F", (   1, 381,0),  None,  beadsize,     2159, [n82_rb], 0,  None,  None, None))
        domains.append(("Nup159.2","Nup159.2_10", 0.9,   f_n82+"Nup159.txt", "Nup159", "BEADS",  " ", ( 382,1116,0),  None,  beadsize100,  2159, [n82_rb], 0,  None,  None, None))
        domains.append(("Nsp1.1",  "Nsp1.1_10",   0.3,   f_n82+"Nsp1.txt",   "Nsp1",   "BEADS",  " ", (   1, 636,0),  None,  beadsize100,  1062, [n82_rb], 0,  None,  None, None))
        domains.append(("Nsp1.2",  "Nsp1.2_10",   0.38,  f_n82+"Nsp1.txt",   "Nsp1",   "BEADS",  " ", (   1, 636,0),  None,  beadsize100,  2062, [n82_rb], 0,  None,  None, None))
        domains.append(("Nup116.1","Nup116.1_10", 0.75,  f_n82+"Nup116.txt", "Nup116", "BEADS",  " ", (   1, 750,0),  None,  beadsize100,  1116, [n82_rb], 0,  None,  None, None))
        domains.append(("Nup116.2","Nup116.2_10", 0.8,   f_n82+"Nup116.txt", "Nup116", "BEADS",  " ", (   1, 750,0),  None,  beadsize100,  2116, [n82_rb], 0,  None,  None, None))
    else:
        domains.append(("Nup159.1","Nup159.1_10", 1.0,   f_n82+"Nup159.txt", "Nup159", "BEADS",  " ", (1082,1116,0),  None,  beadsize100, n82_rb,[n82_rb], 0,  None,  None, None))
        domains.append(("Nup159.2","Nup159.2_10", 0.9,   f_n82+"Nup159.txt", "Nup159", "BEADS",  " ", (1082,1116,0),  None,  beadsize100, n82_rb,[n82_rb], 0,  None,  None, None))
        domains.append(("Nsp1.1",  "Nsp1.1_10",   0.3,   f_n82+"Nsp1.txt",   "Nsp1",   "BEADS",  " ", ( 601, 636,0),  None,  beadsize100, n82_rb,[n82_rb], 0,  None,  None, None))
        domains.append(("Nsp1.2",  "Nsp1.2_10",   0.38,  f_n82+"Nsp1.txt",   "Nsp1",   "BEADS",  " ", ( 601, 636,0),  None,  beadsize100, n82_rb,[n82_rb], 0,  None,  None, None))
    domains.append(("Nup159.1","Nup159.1",    1.0,   f_n82+"Nup159.txt", "Nup159", n82_pdb,  "E", (1117,1460,0),  gmm,   beadsize,   n82_rb, [n82_rb], 6,  " ",   " ",  None, False))
    domains.append(("Nup159.2","Nup159.2",    0.9,   f_n82+"Nup159.txt", "Nup159", n82_pdb,  "F", (1117,1460,0),  gmm,   beadsize,   n82_rb, [n82_rb], 6,  " ",   " ",  None, False))
    domains.append(("Nsp1.1",  "Nsp1.1",      0.3,   f_n82+"Nsp1.txt",   "Nsp1",   n82_pdb,  "G", ( 637, 823,0),  gmm,   beadsize,   n82_rb, [n82_rb], 4,  " ",   " ",  None, False))
    domains.append(("Nsp1.2",  "Nsp1.2",      0.38,  f_n82+"Nsp1.txt",   "Nsp1",   n82_pdb,  "H", ( 637, 823,0),  gmm,   beadsize,   n82_rb, [n82_rb], 4,  " ",   " ",  None, False))
    domains.append(("Nup116.1", "Nup116.1",   0.75,  f_n82+"Nup116.txt", "Nup116", n116_pdbfile,"I",(751,1113,0), None,  beadsize25, n82_rb, [n82_rb], 1,  " ",   " ",  None, False))
    domains.append(("Nup116.2", "Nup116.2",   0.8,   f_n82+"Nup116.txt", "Nup116", n116_pdbfile,"J",(751,1113,0), None,  beadsize25, n82_rb, [n82_rb], 1,  " ",   " ",  None, False))
    for i in clones_range_B:
        domains.append(("Dyn2.1@%d"%i,  "Dyn2.1@%d"%i,      0.48,  f_n82+"Dyn2.txt",   "Dyn2",   n82_pdb,  "A", (   1,  92,0),  None, beadsize,    None,  None,  1,  None,  None,  None))
        domains.append(("Dyn2.2@%d"%i,  "Dyn2.2@%d"%i,      0.65,  f_n82+"Dyn2.txt",   "Dyn2",   n82_pdb,  "B", (   1,  92,0),  None, beadsize,    None,  None,  1,  None,  None,  None))
        domains.append(("Nup82.1@%d"%i, "Nup82.1_1@%d"%i,   0.0,   f_n82+"Nup82.txt",  "Nup82",  n82_pdb,  "C", (   1, 452,0),  None, beadsize,    None,  None,  2,  None,  None,  None))
        domains.append(("Nup82.1@%d"%i, "Nup82.1_2@%d"%i,   0.0,   f_n82+"Nup82.txt",  "Nup82",  n82_pdb,  "C", ( 453, 713,0),  None, beadsize,    None,  None,  4,  None,  None,  None))
        domains.append(("Nup82.2@%d"%i, "Nup82.2_1@%d"%i,   0.15,  f_n82+"Nup82.txt",  "Nup82",  n82_pdb,  "D", (   1, 452,0),  None, beadsize,    None,  None,  2,  None,  None,  None))
        domains.append(("Nup82.2@%d"%i, "Nup82.2_2@%d"%i,   0.15,  f_n82+"Nup82.txt",  "Nup82",  n82_pdb,  "D", ( 453, 713,0),  None, beadsize,    None,  None,  4,  None,  None,  None))
        if (is_FG):
            domains.append(("Nup159.1@%d"%i,"Nup159.1_1@%d"%i,  1.0,   f_n82+"Nup159.txt", "Nup159", n82_pdb,  "E", (   1, 381,0),  None, beadsize,    None,  None,  0,  None,  None,  None))
            domains.append(("Nup159.1@%d"%i,"Nup159.1_10@%d"%i, 1.0,   f_n82+"Nup159.txt", "Nup159", "BEADS",  " ", ( 382,1116,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
            domains.append(("Nup159.2@%d"%i,"Nup159.2_1@%d"%i,  0.9,   f_n82+"Nup159.txt", "Nup159", n82_pdb,  "F", (   1, 381,0),  None, beadsize,    None,  None,  0,  None,  None,  None))
            domains.append(("Nup159.2@%d"%i,"Nup159.2_10@%d"%i, 0.9,   f_n82+"Nup159.txt", "Nup159", "BEADS",  " ", ( 382,1116,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
            domains.append(("Nsp1.1@%d"%i,  "Nsp1.1_10@%d"%i,   0.3,   f_n82+"Nsp1.txt",   "Nsp1",   "BEADS",  " ", (   1, 636,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
            domains.append(("Nsp1.2@%d"%i,  "Nsp1.2_10@%d"%i,   0.38,  f_n82+"Nsp1.txt",   "Nsp1",   "BEADS",  " ", (   1, 636,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
            domains.append(("Nup116.1@%d"%i,"Nup116.1_10@%d"%i, 0.75,  f_n82+"Nup116.txt", "Nup116", "BEADS",  " ", (   1, 750,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
            domains.append(("Nup116.2@%d"%i,"Nup116.2_10@%d"%i, 0.8,   f_n82+"Nup116.txt", "Nup116", "BEADS",  " ", (   1, 750,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
        else:
            domains.append(("Nup159.1@%d"%i,"Nup159.1_10@%d"%i, 1.0,   f_n82+"Nup159.txt", "Nup159", "BEADS",  " ", (1082,1116,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
            domains.append(("Nup159.2@%d"%i,"Nup159.2_10@%d"%i, 0.9,   f_n82+"Nup159.txt", "Nup159", "BEADS",  " ", (1082,1116,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
            domains.append(("Nsp1.1@%d"%i,  "Nsp1.1_10@%d"%i,   0.3,   f_n82+"Nsp1.txt",   "Nsp1",   "BEADS",  " ", ( 601, 636,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
            domains.append(("Nsp1.2@%d"%i,  "Nsp1.2_10@%d"%i,   0.38,  f_n82+"Nsp1.txt",   "Nsp1",   "BEADS",  " ", ( 601, 636,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
        domains.append(("Nup159.1@%d"%i,"Nup159.1@%d"%i,    1.0,   f_n82+"Nup159.txt", "Nup159", n82_pdb,  "E", (1117,1460,0),  None, beadsize,    None,  None,  6,  None,  None,  None))
        domains.append(("Nup159.2@%d"%i,"Nup159.2@%d"%i,    0.9,   f_n82+"Nup159.txt", "Nup159", n82_pdb,  "F", (1117,1460,0),  None, beadsize,    None,  None,  6,  None,  None,  None))
        domains.append(("Nsp1.1@%d"%i,  "Nsp1.1@%d"%i,      0.3,   f_n82+"Nsp1.txt",   "Nsp1",   n82_pdb,  "G", ( 637, 823,0),  None, beadsize,    None,  None,  4,  None,  None,  None))
        domains.append(("Nsp1.2@%d"%i,  "Nsp1.2@%d"%i,      0.38,  f_n82+"Nsp1.txt",   "Nsp1",   n82_pdb,  "H", ( 637, 823,0),  None, beadsize,    None,  None,  4,  None,  None,  None))
        domains.append(("Nup116.1@%d"%i,"Nup116.1@%d"%i,    0.75,  f_n82+"Nup116.txt", "Nup116", n116_pdbfile,"I",( 751,1113,0),None, beadsize25,  None,  None,  1,  None,  None,  None))
        domains.append(("Nup116.2@%d"%i,"Nup116.2@%d"%i,    0.8,   f_n82+"Nup116.txt", "Nup116", n116_pdbfile,"J",( 751,1113,0),None, beadsize25,  None,  None,  1,  None,  None,  None))

##########################
# Nic96 complex
##########################
if (is_nic96):
    n96_rb = 196
    domains.append(("Nic96.1",  "Nic96.1_1",  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (  1,106,0),  gmm,  beadsize,    n96_rb, [n96_rb], 2,  " ",   " ",  None, False))
    domains.append(("Nic96.1",  "Nic96.1_2",  0.25, f_n96+"Nic96.txt", "YFR002W", "BEADS",  " ",  (107,204,0),  None, beadsize25,  1096,   [n96_rb], 0,  None,  None, None))
    #domains.append(("Nic96.1",  "Nic96.1_2",  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (107,204,0),  gmm,  beadsize,    1097,   [n96_rb], 1,  " ",   " ",  None, False))
    domains.append(("Nic96.1",  "Nic96.1",    0.25, f_n96+"Nic96.txt", "YFR002W", n961_pdb, "A",  (205,839,0),  gmm,  beadsize,    1096,   [n96_rb], 3,  " ",   " ",  None, False))
    if (is_FG):
        domains.append(("Nsp1.3",   "Nsp1.3_10",  0.50, f_n96+"Nsp1.txt",  "YJL041W", "BEADS",  " ",  (  1,636,0),  None, beadsize100, n96_rb, [n96_rb], 0,  None,  None, None))
        domains.append(("Nup49.1",  "Nup49.1_10", 0.75, f_n96+"Nup49.txt", "YGL172W", "BEADS",  " ",  (  1,269,0),  None, beadsize100, n96_rb, [n96_rb], 0,  None,  None, None))
        domains.append(("Nup57.1",  "Nup57.1_10", 1.0,  f_n96+"Nup57.txt", "YGR119C", "BEADS",  " ",  (  1,286,0),  None, beadsize100, n96_rb, [n96_rb], 0,  None,  None, None))
    else:
        domains.append(("Nsp1.3",   "Nsp1.3_10",  0.50, f_n96+"Nsp1.txt",  "YJL041W", "BEADS",  " ",  (601,636,0),  None, beadsize100, n96_rb, [n96_rb], 0,  None,  None, None))
        domains.append(("Nup49.1",  "Nup49.1_10", 0.75, f_n96+"Nup49.txt", "YGL172W", "BEADS",  " ",  (201,269,0),  None, beadsize100, n96_rb, [n96_rb], 0,  None,  None, None))
        domains.append(("Nup57.1",  "Nup57.1_10", 1.0,  f_n96+"Nup57.txt", "YGR119C", "BEADS",  " ",  (201,286,0),  None, beadsize100, n96_rb, [n96_rb], 0,  None,  None, None))
    domains.append(("Nsp1.3",   "Nsp1.3",     0.50, f_n96+"Nsp1.txt",  "YJL041W", n96_pdb,  "B",  (637,823,0),  gmm,  beadsize,    n96_rb, [n96_rb], 4,  " ",   " ",  None, False))
    domains.append(("Nup49.1",  "Nup49.1",    0.75, f_n96+"Nup49.txt", "YGL172W", n96_pdb,  "C",  (270,472,0),  gmm,  beadsize,    n96_rb, [n96_rb], 4,  " ",   " ",  None, False))
    domains.append(("Nup57.1",  "Nup57.1",    1.0,  f_n96+"Nup57.txt", "YGR119C", n96_pdb,  "D",  (287,541,0),  gmm,  beadsize,    n96_rb, [n96_rb], 5,  " ",   " ",  None, False))

    n96_rb = 296
    domains.append(("Nic96.2",  "Nic96.2_1",  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (  1,106,0),  gmm,  beadsize,    n96_rb, [n96_rb], 2,  gmm_f+"Nic96.1_1.txt", gmm_f+"Nic96.1_1.mrc", None, False))
    domains.append(("Nic96.2",  "Nic96.2_2",  0.25, f_n96+"Nic96.txt", "YFR002W", "BEADS",  " ",  (107,204,0),  None, beadsize25,  2096,   [n96_rb], 0,  None,  None, None))
    #domains.append(("Nic96.2",  "Nic96.2_2",  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (107,204,0),  gmm,  beadsize,    2097,   [n96_rb], 1,  gmm_f+"Nic96.1_2.txt", gmm_f+"Nic96.1_2.mrc", None, False))
    domains.append(("Nic96.2",  "Nic96.2",    0.25, f_n96+"Nic96.txt", "YFR002W", n962_pdb, "B",  (205,839,0),  gmm,  beadsize,    2096,   [n96_rb], 3,  " ",   " ",  None, False))
    if (is_FG):
        domains.append(("Nsp1.4",   "Nsp1.4_10",  0.50, f_n96+"Nsp1.txt",  "YJL041W", "BEADS",  " ",  (  1,636,0),  None, beadsize100, n96_rb, [n96_rb], 0,  None,  None, None))
        domains.append(("Nup49.2",  "Nup49.2_10", 0.75, f_n96+"Nup49.txt", "YGL172W", "BEADS",  " ",  (  1,269,0),  None, beadsize100, n96_rb, [n96_rb], 0,  None,  None, None))
        domains.append(("Nup57.2",  "Nup57.2_10", 1.0,  f_n96+"Nup57.txt", "YGR119C", "BEADS",  " ",  (  1,286,0),  None, beadsize100, n96_rb, [n96_rb], 0,  None,  None, None))
    else:
        domains.append(("Nsp1.4",   "Nsp1.4_10",  0.50, f_n96+"Nsp1.txt",  "YJL041W", "BEADS",  " ",  (601,636,0),  None, beadsize100, n96_rb, [n96_rb], 0,  None,  None, None))
        domains.append(("Nup49.2",  "Nup49.2_10", 0.75, f_n96+"Nup49.txt", "YGL172W", "BEADS",  " ",  (201,269,0),  None, beadsize100, n96_rb, [n96_rb], 0,  None,  None, None))
        domains.append(("Nup57.2",  "Nup57.2_10", 1.0,  f_n96+"Nup57.txt", "YGR119C", "BEADS",  " ",  (201,286,0),  None, beadsize100, n96_rb, [n96_rb], 0,  None,  None, None))
    domains.append(("Nsp1.4",   "Nsp1.4",     0.50, f_n96+"Nsp1.txt",  "YJL041W", n96_pdb,  "B",  (637,823,0),  gmm,  beadsize,    n96_rb, [n96_rb], 4,  gmm_f+"Nsp1.3.txt",    gmm_f+"Nsp1.3.mrc",    None, False))
    domains.append(("Nup49.2",  "Nup49.2",    0.75, f_n96+"Nup49.txt", "YGL172W", n96_pdb,  "C",  (270,472,0),  gmm,  beadsize,    n96_rb, [n96_rb], 4,  gmm_f+"Nup49.1.txt",   gmm_f+"Nup49.1.mrc",   None, False))
    domains.append(("Nup57.2",  "Nup57.2",    1.0,  f_n96+"Nup57.txt", "YGR119C", n96_pdb,  "D",  (287,541,0),  gmm,  beadsize,    n96_rb, [n96_rb], 5,  gmm_f+"Nup57.1.txt",   gmm_f+"Nup57.1.mrc",   None, False))

    for i in clones_range_A:
        if (i==11): gmm_c = gmm
        else:       gmm_c = None
        domains.append(("Nic96.1@%d"%i,  "Nic96.1_1@%d"%i,  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (  1,106,0),  gmm_c,  beadsize,    None, None, 2,  gmm_f+"Nic96.1_1.txt", gmm_f+"Nic96.1_1.mrc", None, False))
        domains.append(("Nic96.1@%d"%i,  "Nic96.1_2@%d"%i,  0.25, f_n96+"Nic96.txt", "YFR002W", "BEADS",  " ",  (107,204,0),  None,   beadsize25,  None, None, 0,  None,  None, None))
        #domains.append(("Nic96.1@%d"%i,  "Nic96.1_2@%d"%i,  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (107,204,0),  gmm_c,  beadsize,    None, None, 1,  gmm_f+"Nic96.1_2.txt", gmm_f+"Nic96.1_2.mrc", None, False))
        domains.append(("Nic96.1@%d"%i,  "Nic96.1@%d"%i,    0.25, f_n96+"Nic96.txt", "YFR002W", n961_pdb, "A",  (205,839,0),  gmm_c,  beadsize,    None, None, 3,  gmm_f+"Nic96.1.txt",   gmm_f+"Nic96.1.mrc",   None, False))
        if (is_FG):
            domains.append(("Nsp1.3@%d"%i,   "Nsp1.3_10@%d"%i,  0.50, f_n96+"Nsp1.txt",  "YJL041W", "BEADS",  " ",  (  1,636,0),  None,   beadsize100, None, None, 0,  None,                  None,                  None))
            domains.append(("Nup49.1@%d"%i,  "Nup49.1_10@%d"%i, 0.75, f_n96+"Nup49.txt", "YGL172W", "BEADS",  " ",  (  1,269,0),  None,   beadsize100, None, None, 0,  None,                  None,                  None))
            domains.append(("Nup57.1@%d"%i,  "Nup57.1_10@%d"%i, 1.0,  f_n96+"Nup57.txt", "YGR119C", "BEADS",  " ",  (  1,286,0),  None,   beadsize100, None, None, 0,  None,                  None,                  None))
        else:
            domains.append(("Nsp1.3@%d"%i,   "Nsp1.3_10@%d"%i,  0.50, f_n96+"Nsp1.txt",  "YJL041W", "BEADS",  " ",  (601,636,0),  None,   beadsize100, None, None, 0,  None,                  None,                  None))
            domains.append(("Nup49.1@%d"%i,  "Nup49.1_10@%d"%i, 0.75, f_n96+"Nup49.txt", "YGL172W", "BEADS",  " ",  (201,269,0),  None,   beadsize100, None, None, 0,  None,                  None,                  None))
            domains.append(("Nup57.1@%d"%i,  "Nup57.1_10@%d"%i, 1.0,  f_n96+"Nup57.txt", "YGR119C", "BEADS",  " ",  (201,286,0),  None,   beadsize100, None, None, 0,  None,                  None,                  None))
        domains.append(("Nsp1.3@%d"%i,   "Nsp1.3@%d"%i,     0.50, f_n96+"Nsp1.txt",  "YJL041W", n96_pdb,  "B",  (637,823,0),  gmm_c,  beadsize,    None, None, 4,  gmm_f+"Nsp1.3.txt",    gmm_f+"Nsp1.3.mrc",    None, False))
        domains.append(("Nup49.1@%d"%i,  "Nup49.1@%d"%i,    0.75, f_n96+"Nup49.txt", "YGL172W", n96_pdb,  "C",  (270,472,0),  gmm_c,  beadsize,    None, None, 4,  gmm_f+"Nup49.1.txt",   gmm_f+"Nup49.1.mrc",   None, False))
        domains.append(("Nup57.1@%d"%i,  "Nup57.1@%d"%i,    1.0,  f_n96+"Nup57.txt", "YGR119C", n96_pdb,  "D",  (287,541,0),  gmm_c,  beadsize,    None, None, 5,  gmm_f+"Nup57.1.txt",   gmm_f+"Nup57.1.mrc",   None, False))

        domains.append(("Nic96.2@%d"%i,  "Nic96.2_1@%d"%i,  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (  1,106,0),  gmm_c,  beadsize,    None, None, 2,  gmm_f+"Nic96.1_1.txt", gmm_f+"Nic96.1_1.mrc", None, False))
        domains.append(("Nic96.2@%d"%i,  "Nic96.2_2@%d"%i,  0.25, f_n96+"Nic96.txt", "YFR002W", "BEADS",  " ",  (107,204,0),  None,   beadsize25,  None, None, 0,  None,  None, None))
        #domains.append(("Nic96.2@%d"%i,  "Nic96.2_2@%d"%i,  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (107,204,0),  gmm_c,  beadsize,    None, None, 1,  gmm_f+"Nic96.1_2.txt", gmm_f+"Nic96.1_2.mrc", None, False))
        domains.append(("Nic96.2@%d"%i,  "Nic96.2@%d"%i,    0.25, f_n96+"Nic96.txt", "YFR002W", n962_pdb, "B",  (205,839,0),  gmm_c,  beadsize,    None, None, 3,  gmm_f+"Nic96.2.txt",   gmm_f+"Nic96.2.mrc",   None, False))
        if (is_FG):
            domains.append(("Nsp1.4@%d"%i,   "Nsp1.4_10@%d"%i,  0.50, f_n96+"Nsp1.txt",  "YJL041W", "BEADS",  " ",  (  1,636,0),  None,   beadsize100, None, None, 0,  None,                  None,                  None))
            domains.append(("Nup49.2@%d"%i,  "Nup49.2_10@%d"%i, 0.75, f_n96+"Nup49.txt", "YGL172W", "BEADS",  " ",  (  1,269,0),  None,   beadsize100, None, None, 0,  None,                  None,                  None))
            domains.append(("Nup57.2@%d"%i,  "Nup57.2_10@%d"%i, 1.0,  f_n96+"Nup57.txt", "YGR119C", "BEADS",  " ",  (  1,286,0),  None,   beadsize100, None, None, 0,  None,                  None,                  None))
        else:
            domains.append(("Nsp1.4@%d"%i,   "Nsp1.4_10@%d"%i,  0.50, f_n96+"Nsp1.txt",  "YJL041W", "BEADS",  " ",  (601,636,0),  None,   beadsize100, None, None, 0,  None,                  None,                  None))
            domains.append(("Nup49.2@%d"%i,  "Nup49.2_10@%d"%i, 0.75, f_n96+"Nup49.txt", "YGL172W", "BEADS",  " ",  (201,269,0),  None,   beadsize100, None, None, 0,  None,                  None,                  None))
            domains.append(("Nup57.2@%d"%i,  "Nup57.2_10@%d"%i, 1.0,  f_n96+"Nup57.txt", "YGR119C", "BEADS",  " ",  (201,286,0),  None,   beadsize100, None, None, 0,  None,                  None,                  None))
        domains.append(("Nsp1.4@%d"%i,   "Nsp1.4@%d"%i,     0.50, f_n96+"Nsp1.txt",  "YJL041W", n96_pdb,  "B",  (637,823,0),  gmm_c,  beadsize,    None, None, 4,  gmm_f+"Nsp1.3.txt",    gmm_f+"Nsp1.3.mrc",    None, False))
        domains.append(("Nup49.2@%d"%i,  "Nup49.2@%d"%i,    0.75, f_n96+"Nup49.txt", "YGL172W", n96_pdb,  "C",  (270,472,0),  gmm_c,  beadsize,    None, None, 4,  gmm_f+"Nup49.1.txt",   gmm_f+"Nup49.1.mrc",   None, False))
        domains.append(("Nup57.2@%d"%i,  "Nup57.2@%d"%i,    1.0,  f_n96+"Nup57.txt", "YGR119C", n96_pdb,  "D",  (287,541,0),  gmm_c,  beadsize,    None, None, 5,  gmm_f+"Nup57.1.txt",   gmm_f+"Nup57.1.mrc",   None, False))

##########################
# Inner Ring scaffold
##########################
if (is_inner_ring):
    domains.append(("Nup157",    "Nup157n",       0.0,  f_npc+"Nup157.txt", "YER105C", n157N_pdbfile,  "A",  (  1, 892,0),  gmm,  beadsize25,  1571, [157], 3, " ",   " ",  None, False))
    domains.append(("Nup157",    "Nup157c",       0.0,  f_npc+"Nup157.txt", "YER105C", n157C_pdbfile,  "A",  (893,1391,0),  gmm,  beadsize25,  1572, [157], 3, " ",   " ",  None, False))
    domains.append(("Nup170",    "Nup170n",       0.1,  f_npc+"Nup170.txt", "Nup170",  n170N_pdbfile,  "A",  (  1, 992,0),  gmm,  beadsize25,  1701, [170], 3, " ",   " ",  None, False))
    domains.append(("Nup170",    "Nup170c",       0.1,  f_npc+"Nup170.txt", "Nup170",  n170C_pdbfile,  "A",  (993,1502,0),  gmm,  beadsize25,  1702, [170], 3, " ",   " ",  None, False))
    domains.append(("Nup188",    "Nup188",        0.85, f_npc+"Nup188.txt", "YML103C", n188_pdbfile,   "A",  (  1,1655,0),  gmm,  beadsize25,  188,  None,  6, " ",   " ",  None, False))
    domains.append(("Nup192",    "Nup192",        0.75, f_npc+"Nup192.txt", "YJL039C", n192_pdbfile,   "A",  (  1,1683,0),  gmm,  beadsize25,  192,  None,  6, " ",   " ",  None, False))
    for i in clones_range_A:
        if (i==11): gmm_c = gmm
        else:       gmm_c = None
        domains.append(("Nup157@%d"%i,  "Nup157n@%d"%i,  0.0,  f_npc+"Nup157.txt", "YER105C", n157N_pdbfile,  "A",  (  1, 892,0),  gmm_c,  beadsize25,  None, None, 3, gmm_f+"Nup157n.txt", gmm_f+"Nup157n.mrc", None, False))
        domains.append(("Nup157@%d"%i,  "Nup157c@%d"%i,  0.0,  f_npc+"Nup157.txt", "YER105C", n157C_pdbfile,  "A",  (893,1391,0),  gmm_c,  beadsize25,  None, None, 3, gmm_f+"Nup157c.txt", gmm_f+"Nup157c.mrc", None, False))
        domains.append(("Nup170@%d"%i,  "Nup170n@%d"%i,  0.1,  f_npc+"Nup170.txt", "Nup170",  n170N_pdbfile,  "A",  (  1, 992,0),  gmm_c,  beadsize25,  None, None, 3, gmm_f+"Nup170n.txt", gmm_f+"Nup170n.mrc", None, False))
        domains.append(("Nup170@%d"%i,  "Nup170c@%d"%i,  0.1,  f_npc+"Nup170.txt", "Nup170",  n170C_pdbfile,  "A",  (993,1502,0),  gmm_c,  beadsize25,  None, None, 3, gmm_f+"Nup170c.txt", gmm_f+"Nup170c.mrc", None, False))
        domains.append(("Nup188@%d"%i,  "Nup188@%d"%i,   0.85, f_npc+"Nup188.txt", "YML103C", n188_pdbfile,   "A",  (  1,1655,0),  gmm_c,  beadsize25,  None, None, 6, gmm_f+"Nup188.txt",  gmm_f+"Nup188.mrc",  None, False))
        domains.append(("Nup192@%d"%i,  "Nup192@%d"%i,   0.75, f_npc+"Nup192.txt", "YJL039C", n192_pdbfile,   "A",  (  1,1683,0),  gmm_c,  beadsize25,  None, None, 6, gmm_f+"Nup192.txt",  gmm_f+"Nup192.mrc",  None, False))

##########################
# Membrane nups
##########################
if (is_membrane):
    domains.append(("Nup53",     "Nup53",         0.0,  f_npc+"Nup53.txt",  "YMR153W", n53_pdbfile,  "A", (   1, 475,0),  gmm,  beadsize50,  53,   [53],  2,  " ",   " ",  None, False))
    domains.append(("Nup59",     "Nup59",         0.66, f_npc+"Nup59.txt",  "YDL088C", n59_pdbfile,  "A", (   1, 528,0),  gmm,  beadsize50,  59,   [53],  2,  " ",   " ",  None, False))
    domains.append(("Ndc1",      "Ndc1",          0.8,  f_npc+"Ndc1.txt",   "YML031W", "BEADS",      " ", (   1, 655,0),  gmm,  beadsize100, 101,  [101], 0,  None,  None, None, False))
    domains.append(("Pom34",     "Pom34",         0.9,  f_npc+"Pom34.txt",  "YLR018C", "BEADS",      " ", (   1, 299,0),  gmm,  beadsize50,  34,   [34],  0,  None,  None, None, False))
    domains.append(("Pom152",    "Pom152_1" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", (   1, 496,0),  gmm,  beadsize50,  None, None,  2,  " ",   " ",  None, False))
    domains.append(("Pom152",    "Pom152_2" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", ( 497, 613,0),  gmm,  beadsize100, None, None,  2,  " ",   " ",  None, False))
    domains.append(("Pom152",    "Pom152_3" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", ( 614, 718,0),  gmm,  beadsize100, None, None,  2,  " ",   " ",  None, False))
    domains.append(("Pom152",    "Pom152_4" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", ( 719, 821,0),  gmm,  beadsize100, None, None,  2,  " ",   " ",  None, False))
    domains.append(("Pom152",    "Pom152_5" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", ( 822, 924,0),  gmm,  beadsize100, None, None,  2,  " ",   " ",  None, False))
    domains.append(("Pom152",    "Pom152_6" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", ( 925,1031,0),  gmm,  beadsize100, None, None,  2,  " ",   " ",  None, False))
    domains.append(("Pom152",    "Pom152_7" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", (1032,1145,0),  gmm,  beadsize100, None, None,  2,  " ",   " ",  None, False))
    domains.append(("Pom152",    "Pom152_8" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", (1146,1236,0),  gmm,  beadsize100, None, None,  2,  " ",   " ",  None, False))
    domains.append(("Pom152",    "Pom152_9" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", (1237,1337,0),  gmm,  beadsize100, None, None,  2,  " ",   " ",  None, False))
    for i in clones_range_A:
        if (i==11): gmm_c = gmm
        else:       gmm_c = None
        domains.append(("Nup53@%d"%i,  "Nup53@%d"%i,    0.0, f_npc+"Nup53.txt",  "YMR153W", n53_pdbfile, "A", (   1, 475,0), gmm_c,  beadsize50,  None, None, 2,  gmm_f+"Nup53.txt",    gmm_f+"Nup53.mrc",    None, False))
        domains.append(("Nup59@%d"%i,  "Nup59@%d"%i,    0.66,f_npc+"Nup59.txt",  "YDL088C", n59_pdbfile, "A", (   1, 528,0), gmm_c,  beadsize50,  None, None, 2,  gmm_f+"Nup59.txt",    gmm_f+"Nup59.mrc",    None, False))
        domains.append(("Ndc1@%d"%i,   "Ndc1@%d"%i,     0.8, f_npc+"Ndc1.txt",   "YML031W", "BEADS",     " ", (   1, 655,0), gmm_c,  beadsize100, None, None, 0,  None,                 None,                 None, False))
        domains.append(("Pom34@%d"%i,  "Pom34@%d"%i,    0.9, f_npc+"Pom34.txt",  "YLR018C", "BEADS",     " ", (   1, 299,0), gmm_c,  beadsize50,  None, None, 0,  None,                 None,                 None, False))
        domains.append(("Pom152@%d"%i, "Pom152_1@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", (   1, 496,0), gmm_c,  beadsize50,  None, None, 2,  gmm_f+"Pom152_1.txt", gmm_f+"Pom152_1.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_2@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", ( 497, 613,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_2.txt", gmm_f+"Pom152_2.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_3@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", ( 614, 718,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_3.txt", gmm_f+"Pom152_3.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_4@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", ( 719, 821,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_4.txt", gmm_f+"Pom152_4.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_5@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", ( 822, 924,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_5.txt", gmm_f+"Pom152_5.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_6@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", ( 925,1031,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_6.txt", gmm_f+"Pom152_6.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_7@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", (1032,1145,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_7.txt", gmm_f+"Pom152_7.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_8@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", (1146,1236,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_8.txt", gmm_f+"Pom152_8.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_9@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", (1237,1337,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_9.txt", gmm_f+"Pom152_9.mrc", None, False))

##########################
# Cytoplasm only - TODO: Dbp5 / Gle2.1 / Gle2.2
##########################
if (is_cytoplasm):
    if (is_FG):
        domains.append(("Nup100.1", "Nup100.1_10", 0.2, f_npc+"Nup100.txt", "YKL068W", "BEADS",      " ", (  1, 550,0), None, beadsize100, 1100, [1100], 0,  None,  None, None, False))
        domains.append(("Nup100.2", "Nup100.2_10", 0.4, f_npc+"Nup100.txt", "YKL068W", "BEADS",      " ", (  1, 550,0), None, beadsize100, 2100, [2100], 0,  None,  None, None, False))
        #domains.append(("Nup116.1", "Nup116.1_10", 0.75,f_n82+"Nup116.txt", "Nup116",  "BEADS",      " ", (  1, 750,0), None, beadsize100, 1116, [1116], 0,  None,  None, None, False))
        #domains.append(("Nup116.2", "Nup116.2_10", 0.8, f_n82+"Nup116.txt", "Nup116",  "BEADS",      " ", (  1, 750,0), None, beadsize100, 2116, [2116], 0,  None,  None, None, False))
        domains.append(("Nup42",    "Nup42_10",    0.6, f_npc+"Nup42.txt",  "YDR192C", "BEADS",      " ", (  1, 363,0), None, beadsize100, 42,   [611],  0,  None,  None, None, False))
    domains.append(("Nup100.1", "Nup100.1",    0.2, f_npc+"Nup100.txt", "YKL068W", n100_pdbfile, "A", (551, 959,0), None, beadsize25,  1100, [1100],2,  " ",   " ",  None, False))
    domains.append(("Nup100.2", "Nup100.2",    0.4, f_npc+"Nup100.txt", "YKL068W", n100_pdbfile, "B", (551, 959,0), None, beadsize25,  2100, [2100],2,  " ",   " ",  None, False))
    #domains.append(("Nup116.1", "Nup116.1",    0.75,f_n82+"Nup116.txt", "Nup116",  n116_pdbfile, "I", (751,1113,0), None, beadsize25,  n82_rb, [n82_rb], 1,  " ",   " ",  None, False))
    #domains.append(("Nup116.2", "Nup116.2",    0.8, f_n82+"Nup116.txt", "Nup116",  n116_pdbfile, "J", (751,1113,0), None, beadsize25,  n82_rb, [n82_rb], 1,  " ",   " ",  None, False))
    domains.append(("Nup42",    "Nup42",       0.6, f_npc+"Nup42.txt",  "YDR192C", "BEADS",      " ", (364, 430,0), None, beadsize50,  42,   [611], 0,  None,  None, None, False))
    domains.append(("Gle1",     "Gle1_10",     0.8, f_npc+"Gle1.txt",   "YDL207W", Gle1N_pdbfile,"A", (  1, 239,0), None, beadsize50,  610,  [611], 2,  " ",   " ",  None, False))
    domains.append(("Gle1",     "Gle1",        0.8, f_npc+"Gle1.txt",   "YDL207W", Gle1C_pdbfile,"B", (240, 538,0), None, beadsize25,  611,  [611], 2,  " ",   " ",  None, False))
    #domains.append(("Gle2.1",   "Gle2.1",      0.9, f_npc+"Gle2.txt",   "YER107C", Gle2_pdbfile, "A", (  1, 365,0), None, beadsize25,  1612, [1612],2,  " ",                 " ",                None, False))
    #domains.append(("Gle2.2",   "Gle2.2",      1.0, f_npc+"Gle2.txt",   "YER107C", Gle2_pdbfile, "A", (  1, 365,0), None, beadsize25,  2612, [2612],2,  gmm_f+"Gle2.1.txt",  gmm_f+"Gle2.1.mrc", None, False))
    for i in clones_range_B:
        if (is_FG):
            domains.append(("Nup100.1@%d"%i, "Nup100.1_10@%d"%i, 0.2, f_npc+"Nup100.txt", "YKL068W", "BEADS",      " ", (  1, 550,0), None, beadsize100, None, None, 0,  None,  None, None))
            domains.append(("Nup100.2@%d"%i, "Nup100.2_10@%d"%i, 0.4, f_npc+"Nup100.txt", "YKL068W", "BEADS",      " ", (  1, 550,0), None, beadsize100, None, None, 0,  None,  None, None))
            #domains.append(("Nup116.1@%d"%i, "Nup116.1_10@%d"%i, 0.75,f_n82+"Nup116.txt", "Nup116",  "BEADS",      " ", (  1, 750,0), None, beadsize100, None, None, 0,  None,  None, None))
            #domains.append(("Nup116.2@%d"%i, "Nup116.2_10@%d"%i, 0.8, f_n82+"Nup116.txt", "Nup116",  "BEADS",      " ", (  1, 750,0), None, beadsize100, None, None, 0,  None,  None, None))
            domains.append(("Nup42@%d"%i,    "Nup42_10@%d"%i,    0.6, f_npc+"Nup42.txt",  "YDR192C", "BEADS",      " ", (  1, 363,0), None, beadsize100, None, None, 0,  None,  None, None))
        domains.append(("Nup100.1@%d"%i, "Nup100.1@%d"%i,    0.2, f_npc+"Nup100.txt", "YKL068W", n100_pdbfile, "A", (551, 959,0), None, beadsize25,  None, None, 2,  None,  None, None))
        domains.append(("Nup100.2@%d"%i, "Nup100.2@%d"%i,    0.4, f_npc+"Nup100.txt", "YKL068W", n100_pdbfile, "B", (551, 959,0), None, beadsize25,  None, None, 2,  None,  None, None))
        #domains.append(("Nup116.1@%d"%i, "Nup116.1@%d"%i,    0.75,f_n82+"Nup116.txt", "Nup116",  n116_pdbfile, "I", (751,1113,0), None, beadsize25,  None, None, 1,  None,  None, None))
        #domains.append(("Nup116.2@%d"%i, "Nup116.2@%d"%i,    0.8, f_n82+"Nup116.txt", "Nup116",  n116_pdbfile, "J", (751,1113,0), None, beadsize25,  None, None, 1,  None,  None, None))
        domains.append(("Nup42@%d"%i,    "Nup42@%d"%i,       0.6, f_npc+"Nup42.txt",  "YDR192C", "BEADS",      " ", (364, 430,0), None, beadsize50,  None, None, 0,  None,  None, None))
        domains.append(("Gle1@%d"%i,     "Gle1_10@%d"%i,     0.8, f_npc+"Gle1.txt",   "YDL207W", Gle1N_pdbfile,"A", (  1, 239,0), None, beadsize50,  None, None, 2,  None,  None, None))
        domains.append(("Gle1@%d"%i,     "Gle1@%d"%i,        0.8, f_npc+"Gle1.txt",   "YDL207W", Gle1C_pdbfile,"B", (240, 538,0), None, beadsize25,  None, None, 2,  None,  None, None))
        #domains.append(("Gle2.1@%d"%i,   "Gle2.1@%d"%i,      0.9, f_npc+"Gle2.txt",   "YER107C", Gle2_pdbfile, "A", (  1, 365,0), None, beadsize25,  None, None, 2,  None,  None, None))
        #domains.append(("Gle2.2@%d"%i,   "Gle2.2@%d"%i,      1.0, f_npc+"Gle2.txt",   "YER107C", Gle2_pdbfile, "A", (  1, 365,0), None, beadsize25,  None, None, 2,  None,  None, None))

##########################
# Nucleoplasm only
##########################
if (is_nucleoplasm):
    if (is_FG):
        domains.append(("Nup145.1", "Nup145.1_10", 0.7, f_npc+"Nup145.txt", "YGL092W", "BEADS",      " ", (  1, 200,0), None, beadsize100, 1145, [1145], 0,  None,  None, None, False))
        domains.append(("Nup145.2", "Nup145.2_10", 0.75,f_npc+"Nup145.txt", "YGL092W", "BEADS",      " ", (  1, 200,0), None, beadsize100, 2145, [2145], 0,  None,  None, None, False))
        domains.append(("Nup1",     "Nup1_10",     0.8, f_npc+"Nup1.txt",   "YOR098C", "BEADS",      " ", (352,1076,0), None, beadsize100, 1,    [1],    0,  None,  None, None, False))
        domains.append(("Nup60.1",  "Nup60.1_10",  0.86,f_npc+"Nup60.txt",  "YAR002W", "BEADS",      " ", (399, 539,0), None, beadsize100, 60,   [60],   0,  None,  None, None, False))
        domains.append(("Nup60.2",  "Nup60.2_10",  0.93,f_npc+"Nup60.txt",  "YAR002W", "BEADS",      " ", (399, 539,0), None, beadsize100, 60,   [60],   0,  None,  None, None, False))
    domains.append(("Nup145.1", "Nup145.1",    0.7, f_npc+"Nup145.txt", "YGL092W", n145N_pdbfile,"A", (201, 605,0), None, beadsize25,  1145, [1145], 2,  " ",   " ",  None, False))
    domains.append(("Nup145.2", "Nup145.2",    0.75,f_npc+"Nup145.txt", "YGL092W", n145N_pdbfile,"B", (201, 605,0), None, beadsize25,  2145, [2145], 2,  " ",   " ",  None, False))
    domains.append(("Nup1",     "Nup1",        0.8, f_npc+"Nup1.txt",   "YOR098C", "BEADS",      " ", (  1, 351,0), None, beadsize50,  1,    [1],    0,  None,  None, None, False))
    domains.append(("Nup60.1",  "Nup60.1",     0.86,f_npc+"Nup60.txt",  "YAR002W", "BEADS",      " ", (  1, 398,0), None, beadsize50,  60,   [60],   0,  None,  None, None, False))
    domains.append(("Nup60.2",  "Nup60.2",     0.93,f_npc+"Nup60.txt",  "YAR002W", "BEADS",      " ", (  1, 398,0), None, beadsize50,  60,   [60],   0,  None,  None, None, False))
    for i in clones_range_B:
        if (is_FG):
            domains.append(("Nup145.1@%d"%i, "Nup145.1_10@%d"%i, 0.7, f_npc+"Nup145.txt", "YGL092W", "BEADS",      " ", (  1, 200,0), None, beadsize100, None, None, 0,  None,  None, None))
            domains.append(("Nup145.2@%d"%i, "Nup145.2_10@%d"%i, 0.75,f_npc+"Nup145.txt", "YGL092W", "BEADS",      " ", (  1, 200,0), None, beadsize100, None, None, 0,  None,  None, None))
            domains.append(("Nup1@%d"%i,     "Nup1_10@%d"%i,     0.8, f_npc+"Nup1.txt",   "YOR098C", "BEADS",      " ", (352,1076,0), None, beadsize100, None, None, 0,  None,  None, None))
            domains.append(("Nup60.1@%d"%i,  "Nup60.1_10@%d"%i,  0.86,f_npc+"Nup60.txt",  "YAR002W", "BEADS",      " ", (399, 539,0), None, beadsize100, None, None, 0,  None,  None, None))
            domains.append(("Nup60.2@%d"%i,  "Nup60.2_10@%d"%i,  0.93,f_npc+"Nup60.txt",  "YAR002W", "BEADS",      " ", (399, 539,0), None, beadsize100, None, None, 0,  None,  None, None))
        domains.append(("Nup145.1@%d"%i, "Nup145.1@%d"%i,    0.7, f_npc+"Nup145.txt", "YGL092W", n145N_pdbfile,"A", (201, 605,0), None, beadsize25,  None, None, 2,  None,  None, None))
        domains.append(("Nup145.2@%d"%i, "Nup145.2@%d"%i,    0.75,f_npc+"Nup145.txt", "YGL092W", n145N_pdbfile,"B", (201, 605,0), None, beadsize25,  None, None, 2,  None,  None, None))
        domains.append(("Nup1@%d"%i,     "Nup1@%d"%i,        0.8, f_npc+"Nup1.txt",   "YOR098C", "BEADS",      " ", (  1, 351,0), None, beadsize50,  None, None, 0,  None,  None, None))
        domains.append(("Nup60.1@%d"%i,  "Nup60.1@%d"%i,     0.86,f_npc+"Nup60.txt",  "YAR002W", "BEADS",      " ", (  1, 398,0), None, beadsize50,  None, None, 0,  None,  None, None))
        domains.append(("Nup60.2@%d"%i,  "Nup60.2@%d"%i,     0.93,f_npc+"Nup60.txt",  "YAR002W", "BEADS",      " ", (  1, 398,0), None, beadsize50,  None, None, 0,  None,  None, None))

##########################
# Basket proteins - Mlp1, Mlp2
##########################
if (is_basket):
    #domains.append((    "Mlp1",      "Mlp1n",     0.0, f_npc+"Mlp1.txt",   "YKR095W", "BEADS",   " ", (  1, 237,0), gmm,  beadsize100, 9191,  [9191],  0,  None,  None, None))
    domains.append((    "Mlp1",      "Mlp1",      0.0, f_npc+"Mlp1.txt",   "YKR095W", "BEADS",   " ", (238, 716,0), gmm,  beadsize50,  9191,  [9191], 0,  None,  None, None))
    domains.append((    "Mlp1",      "Mlp1c",     0.0, f_npc+"Mlp1.txt",   "YKR095W", "BEADS",   " ", (717,1875,0), gmm,  beadsize100, 9192,  [9191], 0,  None,  None, None, False))
    #domains.append((    "Mlp2",      "Mlp2n",     0.2, f_npc+"Mlp2.txt",   "YIL149C", "BEADS",   " ", (  1, 214,0), gmm,  beadsize100, 9193,  [9191],  0,  None,  None, None))
    domains.append((    "Mlp2",      "Mlp2",      0.2, f_npc+"Mlp2.txt",   "YIL149C", "BEADS",   " ", (215, 690,0), gmm,  beadsize50,  9191,  [9191], 0,  None,  None, None))
    domains.append((    "Mlp2",      "Mlp2c",     0.2, f_npc+"Mlp2.txt",   "YIL149C", "BEADS",   " ", (691,1679,0), gmm,  beadsize100, 9192,  [9191], 0,  None,  None, None, False))
    for i in clones_range_B:
        #domains.append(("Mlp1@%d"%i, "Mlp1n@%d"%i,0.0, f_npc+"Mlp1.txt",   "YKR095W", "BEADS",   " ", (  1, 237,0), None, beadsize100, None,  None,   0,  None,  None, None))
        domains.append(("Mlp1@%d"%i, "Mlp1@%d"%i, 0.0, f_npc+"Mlp1.txt",   "YKR095W", "BEADS",   " ", (238, 716,0), None, beadsize50,  None,  None,   0,  None,  None, None))
        domains.append(("Mlp1@%d"%i, "Mlp1c@%d"%i,0.0, f_npc+"Mlp1.txt",   "YKR095W", "BEADS",   " ", (717,1875,0), None, beadsize100, None,  None,   0,  None,  None, None))
        #domains.append(("Mlp2@%d"%i, "Mlp2n@%d"%i,0.2, f_npc+"Mlp2.txt",   "YIL149C", "BEADS",   " ", (  1, 214,0), None, beadsize100, None,  None,   0,  None,  None, None))
        domains.append(("Mlp2@%d"%i, "Mlp2@%d"%i, 0.2, f_npc+"Mlp2.txt",   "YIL149C", "BEADS",   " ", (215, 690,0), None, beadsize50,  None,  None,   0,  None,  None, None))
        domains.append(("Mlp2@%d"%i, "Mlp2c@%d"%i,0.2, f_npc+"Mlp2.txt",   "YIL149C", "BEADS",   " ", (691,1679,0), None, beadsize100, None,  None,   0,  None,  None, None))

#####################################################
# Model Building
#####################################################
bm1 = IMP.pmi.macros.BuildModel1(simo)
bm1.set_gmm_models_directory(gmm_f)

if (True):
    if (is_n84):
        n84=['Nup84', 'Nup85', 'Nup120', 'Nup133', 'Nup145c', 'Seh1', 'Sec13']
        for d in list(n84):
            bm1.set_rmf_file(d, inputs.rmf_input, 0)
            #bm1.set_rmf_file(d, "../data_npc/Outer_ring_rmfs/OR_84844_newEM_NE.rmf3", 0)

    if (is_n82):
        n82=['Dyn2.1', 'Dyn2.2', 'Nup82.1', 'Nup82.2', 'Nup159.1', 'Nup159.2', 'Nsp1.1', 'Nsp1.2', 'Nup116.1', 'Nup116.2']
        for d in list(n82):
            bm1.set_rmf_file(d, inputs.rmf_input, 0)
            #bm1.set_rmf_file(d, "../data_npc/Outer_ring_rmfs/OR_846_0_best_newEM_cytoplasm.rmf3", 0)

    if (is_nic96):
        nic96=['Nic96.1', 'Nsp1.3', 'Nup49.1', 'Nup57.1', 'Nic96.2', 'Nsp1.4', 'Nup49.2', 'Nup57.2']
        for d in list(nic96):
            bm1.set_rmf_file(d, inputs.rmf_input, 0)
            #bm1.set_rmf_file(d, "../data_npc/Inner_ring_rmfs/Nic96complex_initial.rmf3", 0)

    if (is_inner_ring):
        inner_ring=['Nup157', 'Nup170', 'Nup188', 'Nup192']
        for d in list(inner_ring):
            bm1.set_rmf_file(d, inputs.rmf_input, 0)
            #bm1.set_rmf_file(d, "../data_npc/Inner_ring_rmfs/IR_865_0_final.rmf3", 0)

    if (is_membrane):
        membrane=['Nup53', 'Nup59', 'Ndc1', 'Pom34']
        #membrane=['Nup53', 'Nup59', 'Ndc1', 'Pom34', 'Pom152']
        for d in list(membrane):
            bm1.set_rmf_file(d, inputs.rmf_input, 0)
            #bm1.set_rmf_file(d, "../data_npc/Inner_ring_rmfs/IR_865_0_final.rmf3", 0)
        bm1.set_rmf_file('Pom152', "../data_npc/Pom152_rmfs/Pom152_new_final.rmf3", 0)

    if (is_cytoplasm):
        cytoplasm=['Nup100.1', 'Nup100.2', 'Nup42', 'Gle1']
        for d in list(cytoplasm):
            bm1.set_rmf_file(d, inputs.rmf_input, 0)
            #if (is_FG): bm1.set_rmf_file(d, "../data_nup82/rmfs/B_8_1-95.rmf3", 0)
            #else:       bm1.set_rmf_file(d, "../data_nup82/rmfs/B_8_1-95_FGtruncated.rmf3", 0)

    if (is_nucleoplasm):
        nucleoplasm=['Nup145.1', 'Nup145.2', 'Nup60.1', 'Nup60.2', 'Nup1']
        for d in list(nucleoplasm):
            bm1.set_rmf_file(d, inputs.rmf_input, 0)

    if (is_basket):
        Mlps=['Mlp1', 'Mlp2']
        for d in list(Mlps):
            bm1.set_rmf_file(d, inputs.rmf_input.replace("NPC", "MLPs"), 0)
            #bm1.set_rmf_file(d, inputs.rmf_input, 0)
            #bm1.set_rmf_file(d, "../data_npc/Mlps_1.rmf3", 0)
else:
    main_spoke = [entry[0] for entry in domains if not '@' in entry[0]]
    main_spoke_unique = sorted(list(set(main_spoke)))
    print ("main_spoke_unique = ", main_spoke_unique)
    for entry in main_spoke_unique:
        bm1.set_rmf_file(entry, inputs.rmf_input, 0)

# remove connectivity for clones
clone_list = [entry[0] for entry in domains if '@' in entry[0]]
clone_list_unique = sorted(list(set(clone_list)))   # Make a unique list
print ("clone_list_unique = ", clone_list_unique)

bm1.build_model(data_structure = domains, sequence_connectivity_scale=3.0, sequence_connectivity_resolution=1.0,
                skip_connectivity_these_domains=clone_list_unique, skip_gaussian_in_rmf=False, skip_gaussian_in_representation=False)
                #skip_connectivity_these_domains=clone_list_unique, skip_gaussian_in_rmf=True, skip_gaussian_in_representation=use_EM3D)
#exit(0)
bm1.scale_bead_radii(100, 0.8)


#####################################################
# apply the rotational symmetry
#####################################################
if (use_neighboring_spokes):
    if (is_n84):
        for protein in ['Nup84', 'Nup85', 'Nup120', 'Nup133', 'Nup145c', 'Seh1', 'Sec13']:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
            simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
            simo.create_rotational_symmetry(protein+'@11', [protein+'@%d'%i for i in range(12,14)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
    if (is_n82):
        for protein in ['Dyn2.1', 'Dyn2.2', 'Nup82.1', 'Nup82.2', 'Nup159.1', 'Nup159.2', 'Nsp1.1', 'Nsp1.2', 'Nup116.1', 'Nup116.2']:
            simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
    if (is_nic96):
        for protein in ['Nic96.1', 'Nic96.2', 'Nsp1.3', 'Nsp1.4', 'Nup49.1', 'Nup49.2', 'Nup57.1', 'Nup57.2']:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
            simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
            simo.create_rotational_symmetry(protein+'@11', [protein+'@%d'%i for i in range(12,14)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
    if (is_inner_ring):
        for protein in ['Nup157', 'Nup170', 'Nup188', 'Nup192']:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
            simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
            simo.create_rotational_symmetry(protein+'@11', [protein+'@%d'%i for i in range(12,14)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
    if (is_membrane):
        for protein in ['Nup53', 'Nup59', 'Ndc1', 'Pom34', 'Pom152']:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
            simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
            simo.create_rotational_symmetry(protein+'@11', [protein+'@%d'%i for i in range(12,14)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
    if (is_cytoplasm):
        for protein in ['Nup100.1', 'Nup100.2', 'Nup42', 'Gle1']:
            simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
    if (is_nucleoplasm):
        for protein in ['Nup145.1', 'Nup145.2', 'Nup60.1', 'Nup60.2', 'Nup1']:
            simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
    if (is_basket):
        for protein in ['Mlp1', 'Mlp2']:
            simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
else:
    if (is_n84):
        for protein in ['Nup84', 'Nup85', 'Nup120', 'Nup133', 'Nup145c', 'Seh1', 'Sec13']:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
    if (is_nic96):
        for protein in ['Nic96.1', 'Nic96.2', 'Nsp1.3', 'Nsp1.4', 'Nup49.1', 'Nup49.2', 'Nup57.1', 'Nup57.2']:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
    if (is_inner_ring):
        for protein in ['Nup157', 'Nup170', 'Nup188', 'Nup192']:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
    if (is_membrane):
        for protein in ['Nup53', 'Nup59', 'Ndc1', 'Pom34', 'Pom152']:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))


#####################################################
# rigidify floppy bodies
#####################################################
rigid_tuples = []
#for protein in ['Nup84', 'Nup85', 'Nup120', 'Nup133', 'Nup145c', 'Seh1', 'Sec13']:
#for protein in ['Nup84', 'Nup85', (715,1037,'Nup120'), (490,1157,'Nup133'), 'Nup145c', 'Seh1', 'Sec13']:
for protein in ['Nup84', (715,1037,'Nup120'), (490,1157,'Nup133')]:
    rigid_tuples.append(protein)
for protein in ['Dyn2.1', 'Dyn2.2', 'Nup82.1', 'Nup82.2']:
    rigid_tuples.append(protein)
#for protein in ['Nup159.1', 'Nup159.2', 'Nsp1.1', 'Nsp1.2', (966,1113,'Nup116.1'),(966,1113,'Nup116.2')]:
#for protein in [(1117,1460,'Nup159.1'),(1117,1460,'Nup159.2'), (637,823,'Nsp1.1'),(637,823,'Nsp1.2'), (966,1113,'Nup116.1'),(966,1113,'Nup116.2')]:
for protein in [(1211,1460,'Nup159.1'),(1211,1460,'Nup159.2'), (637,823,'Nsp1.1'),(637,823,'Nsp1.2'), (966,1113,'Nup116.1'),(966,1113,'Nup116.2')]:
    rigid_tuples.append(protein)
#for protein in [(1,56,'Nic96.1'),(205,839,'Nic96.1'), (637,823,'Nsp1.3'), (270,472,'Nup49.1'), (287,541,'Nup57.1')]:
for protein in [(637,823,'Nsp1.3'), (270,472,'Nup49.1'), (287,541,'Nup57.1')]:
    rigid_tuples.append(protein)
#for protein in [(1,56,'Nic96.2'),(205,839,'Nic96.2'), (637,823,'Nsp1.4'), (270,472,'Nup49.2'), (287,541,'Nup57.2')]:
for protein in [(637,823,'Nsp1.4'), (270,472,'Nup49.2'), (287,541,'Nup57.2')]:
    rigid_tuples.append(protein)
#for protein in ['Nup157', 'Nup170', 'Nup188', 'Nup192']:
for protein in [(88,892,'Nup157'),(900,1391,'Nup157'), (98,992,'Nup170'),(1000,1502,'Nup170'), 'Nup188', 'Nup192']:
    rigid_tuples.append(protein)
#for protein in [(248,475,'Nup53'), (266,528,'Nup59'), 'Ndc1', 'Pom34', (379,1337,'Pom152'), 'Mlp1', 'Mlp2']:
for protein in [(248,360,'Nup53'), (266,402,'Nup59'), (379,1337,'Pom152'), (238,716,'Mlp1'), (717,1875,'Mlp1'), (215,690,'Mlp2'), (691,1679,'Mlp2')]:
    rigid_tuples.append(protein)
# Remove flexible movers for all clones
for protein in clone_list_unique:
    rigid_tuples.append(protein)

print ("rigid_tuples = ", rigid_tuples)
for rt in rigid_tuples:
    hs = IMP.pmi.tools.select_by_tuple(simo,rt)
    simo.remove_floppy_bodies(hs)


#####################################################
# randomize the initial configuration
#####################################################
if (use_shuffle) :
    #simo.shuffle_configuration(max_translation=1, avoidcollision=False, ignore_initial_coordinates=True)
    #simo.shuffle_configuration(bounding_box=((350, -150, 50), (1050, 200, 550)), ignore_initial_coordinates=True, cutoff=1.0, niterations=1000)
    simo.shuffle_configuration(bounding_box=((150, -150, 0), (750, 150, 350)), ignore_initial_coordinates=True, cutoff=1.0, niterations=1000)


#####################################################
# defines the movers
#####################################################
simo.set_rigid_bodies_max_rot(rbmaxrot)
simo.set_floppy_bodies_max_trans(fbmaxtrans)
simo.set_rigid_bodies_max_trans(rbmaxtrans)
### These two below are already executed in BuildModel1
#simo.set_floppy_bodies()
#simo.setup_bonds()

outputobjects.append(simo)
sampleobjects.append(simo)


#####################################################
# Restraints setup - Immuno-EM
# Supplementary Table 7. Upper and lower bounds on R-radial restraints of C-terminal bead of nups
# NupType : (min R value, max R value) (in Angstrom)
# Supplementary Table 7. Upper and lower bounds on Z-axial restraints of C-terminal bead of nups
# NupType : (min Z value, max Z value) (in Angstrom)
#####################################################
radial_weight = zaxial_weight = yaxial_weight = 10.0
if (use_Immuno_EM):
    nup_list = [entry[0] for entry in domains if not '@' in entry[0]]
    nup_list_unique = sorted(list(set(nup_list)))   # Make a unique list

    RADIAL = {
        "Ndc1"     : [106, 448],
        "Nic96.1"  : [  0, 633],
        "Nic96.2"  : [  0, 633],
        "Nup120"   : [190, 550],
        "Nup133"   : [190, 550],
        "Nup145c"  : [190, 550],
        "Nup157"   : [106, 448],
        "Nup159.1" : [ 43, 637],
        "Nup159.2" : [ 43, 637],
        "Nup170"   : [106, 448],
        "Nup188"   : [152, 368],
        "Nup192"   : [152, 368],
        "Nup49.1"  : [  0, 633],
        "Nup49.2"  : [  0, 633],
        "Nup53"    : [106, 448],
        "Nup57.1"  : [  0, 633],
        "Nup57.2"  : [  0, 633],
        "Nup59"    : [106, 448],
        "Nup82.1"  : [ 43, 637],
        "Nup82.2"  : [ 43, 637],
        "Nup84"    : [190, 550],
        "Nup85"    : [190, 550],
        "Pom34"    : [106, 448],
        "Pom152"   : [266, 734]
    }
    print "\nXYRadialPositionRestraint !!"
    for protein, r in RADIAL.iteritems():
        if (protein not in nup_list_unique):
            continue
        xyr = IMP.npc.npc_restraints.XYRadialPositionRestraint(simo, protein, lower_bound=r[0], upper_bound=r[1], consider_radius=False, sigma=1.0)
        xyr.set_label('Lower_%d_Upper_%d_%s' % (r[0], r[1], protein))
        xyr.set_weight(radial_weight)
        xyr.add_to_model()
        outputobjects.append(xyr)
        print (xyr.get_output())

    ZAXIAL = {
        "Ndc1"     : [-70, 152],
        "Nic96.1"  : [-60, 235],
        "Nic96.2"  : [-60, 235],
        "Nup120"   : [ 38, 240],
        "Nup133"   : [ 38, 240],
        "Nup145c"  : [ 38, 240],
        "Nup157"   : [-70, 152],
        "Nup159.1" : [ 72, 355],
        "Nup159.2" : [ 72, 355],
        "Nup170"   : [-70, 152],
        "Nup188"   : [ 16, 124],
        "Nup192"   : [-12, 132],
        "Nup49.1"  : [-60, 235],
        "Nup49.2"  : [-60, 235],
        "Nup53"    : [-70, 152],
        "Nup57.1"  : [-60, 235],
        "Nup57.2"  : [-60, 235],
        "Nup59"    : [-70, 152],
        "Nup82.1"  : [ 72, 355],
        "Nup82.2"  : [ 72, 355],
        "Nup84"    : [ 38, 240],
        "Nup85"    : [ 38, 240],
        "Pom34"    : [-70, 152],
        "Pom152"   : [-68, 148]
    }
    print "\nZAxialPositionRestraint !!"
    for protein, z in ZAXIAL.iteritems():
        if (protein not in nup_list_unique):
            continue
        zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, protein, lower_bound=z[0], upper_bound=z[1], consider_radius=False, sigma=1.0)
        zax.set_label('Lower_%d_Upper_%d_%s' % (z[0], z[1], protein))
        zax.set_weight(zaxial_weight)
        zax.add_to_model()
        outputobjects.append(zax)
        print (zax.get_output())


#####################################################
# Restraints setup - Cytoplasm / Nucleoplasm / Basket
#####################################################
if (is_cytoplasm or is_nucleoplasm or is_basket):
    nup_list = [entry[0] for entry in domains if not '@' in entry[0]]
    nup_list_unique = sorted(list(set(nup_list)))   # Make a unique list

    ZAXIAL = {
        "Nup100.1" : [   0, 400],
        "Nup100.2" : [   0, 400],
        "Nup116.1" : [   0, 400],
        "Nup116.2" : [   0, 400],
        "Nup42"    : [   0, 400],
        "Gle1"     : [   0, 400],
        "Gle2.1"   : [   0, 400],
        "Gle2.2"   : [   0, 400],
        "Nup145.1" : [-400,   0],
        "Nup145.2" : [-400,   0],
        "Nup1"     : [-400,   0],
        "Nup60.1"  : [-400,   0],
        "Nup60.2"  : [-400,   0],
        "Mlp1"     : [-700,-200],
        "Mlp2"     : [-700,-200]
    }
    print "\nZAxialPositionRestraints for Cytoplasm / Nucleoplasm / Basket !!"
    for protein, z in ZAXIAL.iteritems():
        if (protein not in nup_list_unique):
            continue
        # applying Z-position restraints on the C-termini of each protein
        zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, protein, lower_bound=z[0], upper_bound=z[1], consider_radius=False, sigma=1.0)
        zax.set_label('Lower_%d_Upper_%d_%s' % (z[0], z[1], protein))
        zax.set_weight(zaxial_weight)
        zax.add_to_model()
        outputobjects.append(zax)
        print (zax.get_output())

    # localizing the Mlp1-Mlp2 hetero-dimer near the outer-ring on the nucleoplasm
    dist_min = -400.0
    dist_max = -200.0
    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (716,716,"Mlp1"), lower_bound=dist_min, upper_bound=dist_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (dist_min, dist_max, "Mlp1_716"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (690,690,"Mlp2"), lower_bound=dist_min, upper_bound=dist_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (dist_min, dist_max, "Mlp2_690"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())


#####################################################
# Restraints setup - FG anchor restraints
#####################################################
if (use_FG_anchor and is_nic96 and not is_FG):
    nup_list = [entry[0] for entry in domains if not '@' in entry[0]]
    nup_list_unique = sorted(list(set(nup_list)))   # Make a unique list

    RADIAL = {
        "Nsp1.3"   : [180, 230],
        "Nup49.1"  : [180, 230],
        "Nup57.1"  : [180, 230],
        "Nsp1.4"   : [180, 230],
        "Nup49.2"  : [180, 230],
        "Nup57.2"  : [180, 230],
        #"Nup1"     : [180, 350],
        #"Nup60.1"  : [180, 350],
        #"Nup60.2"  : [180, 350]
    }
    print "\nFG_XYRadialPositionRestraint !!"
    radial_weight = 10.0
    for protein, r in RADIAL.iteritems():
        if (protein not in nup_list_unique):
            continue
        if (protein in ['Nup1', 'Nup60.1', 'Nup60.2']):
            xyr = IMP.npc.npc_restraints.XYRadialPositionRestraint(simo, protein, lower_bound=r[0], upper_bound=r[1], consider_radius=False, sigma=1.0, term='C')
        else:
            xyr = IMP.npc.npc_restraints.XYRadialPositionRestraint(simo, protein, lower_bound=r[0], upper_bound=r[1], consider_radius=False, sigma=1.0, term='N')
        xyr.set_label('Lower_%d_Upper_%d_%s' % (r[0], r[1], protein))
        xyr.set_weight(radial_weight)
        xyr.add_to_model()
        outputobjects.append(xyr)
        print (xyr.get_output())


#####################################################
# Restraints setup
# Distance restraints
#####################################################
dist_min = 3.0
dr_weight = 10.0

# Nup145n - Nup145c
if (is_n84 and is_nucleoplasm):
    dist_max = 10.0
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(605,605,"Nup145.1"), (1,1,"Nup145c@11"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Nup145n-Nup145c@11")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

# Nup60 homo-dimer cross-link
if (is_nucleoplasm):
    dist_max = 25.0
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(151,151,"Nup60.1"), (151,151,"Nup60.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Nup60.1_151-Nup60.2_151")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

"""
# Nup120 - Nup133 to form the outer ring  (Seo et al, PNAS 2009) ; Not sure if it is real
if (is_n84 and use_neighboring_spokes):
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(11,11,"Nup133"), (641,641,"Nup120@2"), distancemin=3.0, distancemax=35.0, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Nup133-Nup120@2")
    dr.set_weight(10.0)
    outputobjects.append(dr)
    print(dr.get_output())
"""

if (is_inner_ring):
    dist_max = 15.0
    # connection between NTD and CTD
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(892,892,"Nup157"), (900,900,"Nup157"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Nup157N-C")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    # connection between NTD and CTD
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(992,992,"Nup170"), (1000,1000,"Nup170"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Nup170N-C")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

if (is_basket):
    # end-to-end distance of Mlps (230-350 angstrom)
    dist_min = 330.0    #dist_min = 230.0
    dist_max = 450.0    #dist_max = 350.0
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(500,500,"Mlp1"), (1875,1875,"Mlp1"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Mlp1_end-to-end")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(500,500,"Mlp2"), (1679,1679,"Mlp2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Mlp2_end-to-end")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    # Distal Ring
    dist_min = 1.35*130.0
    dist_max = 1.35*170.0
    xyr = IMP.npc.npc_restraints.XYRadialPositionRestraint(simo, (1875,1875,"Mlp1"), lower_bound=dist_min, upper_bound=dist_max, consider_radius=False, sigma=1.0, term='M')
    xyr.set_label('Lower_%d_Upper_%d_%s' % (dist_min, dist_max, "Mlp1_1875"))
    xyr.set_weight(radial_weight)
    xyr.add_to_model()
    outputobjects.append(xyr)
    print (xyr.get_output())

    xyr = IMP.npc.npc_restraints.XYRadialPositionRestraint(simo, (1679,1679,"Mlp2"), lower_bound=dist_min, upper_bound=dist_max, consider_radius=False, sigma=1.0, term='M')
    xyr.set_label('Lower_%d_Upper_%d_%s' % (dist_min, dist_max, "Mlp2_1679"))
    xyr.set_weight(radial_weight)
    xyr.add_to_model()
    outputobjects.append(xyr)
    print (xyr.get_output())

    # localizing the Mlp1-Mlp2 hetero-dimer near the x-axis
    dist_min = 0.0
    dist_max = 40.0
    yax = IMP.npc.npc_restraints.YAxialPositionRestraint(simo, (1875,1875,"Mlp1"), lower_bound=dist_min, upper_bound=dist_max, consider_radius=False, sigma=1.0, term='M')
    yax.set_label('Lower_%d_Upper_%d_%s' % (dist_min, dist_max, "Mlp1_1875"))
    yax.set_weight(yaxial_weight)
    yax.add_to_model()
    outputobjects.append(yax)
    print (yax.get_output())

    yax = IMP.npc.npc_restraints.YAxialPositionRestraint(simo, (1679,1679,"Mlp2"), lower_bound=dist_min, upper_bound=dist_max, consider_radius=False, sigma=1.0, term='M')
    yax.set_label('Lower_%d_Upper_%d_%s' % (dist_min, dist_max, "Mlp2_1679"))
    yax.set_weight(yaxial_weight)
    yax.add_to_model()
    outputobjects.append(yax)
    print (yax.get_output())

    """
    # Mlp1-Mlp2 hetero-dimer
    dist_min = 3.0
    dist_max = 15.0
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(1875,1875,"Mlp1"), (1679,1679,"Mlp2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Mlp1-Mlp2")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())
    """


#####################################################
# Restraints setup - Membrane Localization + ALPS Motif
# Campelo et al, PLOS CompBio, 2014 (PMC3983069)
#####################################################
tor_th      = 45.0
tor_th_ALPS = 12.0
tor_R       = 390.0 + 150.0
tor_r       = 150.0 - tor_th/2.0
tor_r_ALPS  = 150.0 - tor_th_ALPS/2.0
msl_sigma   = 1.0
msl_weight  = 1000.0

# Transmembrane domains
if (is_membrane):
    print "\nMembraneSurfaceLocationRestraint !!"
    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (111,194,'Pom152'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Pom152_111_194')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    # Pom152 201-1337 should be outside of the lipid bilayer
    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (201,1337,'Pom152'), tor_R=tor_R, tor_r=52.5, tor_th=105.0, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Pom152_201_1337')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (29,247,'Ndc1'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Ndc1')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (64,150,'Pom34'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Pom34')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (475,475,'Nup53'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup53')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (528,528,'Nup59'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup59')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

# The Pom152 ring
if (is_membrane):
    dist_min = 3.0
    dist_max = 25.0

    # same residue cross-link of Pom152 62-62
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(62,62,"Pom152"), (62,62,"Pom152@11"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Pom152-Pom152@11")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    """
    if (use_neighboring_spokes):
        # TODO: Pom152 orientation?  (clockwise or counter-clockwise?)
        #dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(351,351,"Pom152"), (1337,1337,"Pom152@12"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(351,351,"Pom152"), (1337,1337,"Pom152@13"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr.add_to_model()
        #dr.set_label("Pom152-Pom152@12")
        dr.set_label("Pom152-Pom152@13")
        dr.set_weight(dr_weight)
        outputobjects.append(dr)
        print(dr.get_output())

    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(351,351,"Pom152"), (351,351,"Pom152@11"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Pom152-Pom152@11")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    xyr = IMP.npc.npc_restraints.XYRadialPositionRestraint(simo, (859,859,"Pom152"), lower_bound=515, upper_bound=550, consider_radius=False, sigma=1.0, term='M')
    xyr.set_label('Lower_%d_Upper_%d_%s' % (515, 550, "Pom152_859"))
    xyr.set_weight(radial_weight)
    xyr.add_to_model()
    outputobjects.append(xyr)
    print (xyr.get_output())

    # TODO: Pom152 orientation?  (clockwise or counter-clockwise?)
    #yax = IMP.npc.npc_restraints.YAxialPositionRestraint(simo, (859,859,"Pom152"), lower_bound=180, upper_bound=205, consider_radius=False, sigma=1.0, term='M')
    yax = IMP.npc.npc_restraints.YAxialPositionRestraint(simo, (859,859,"Pom152"), lower_bound=-205, upper_bound=-180, consider_radius=False, sigma=1.0, term='M')
    #yax.set_label('Lower_%d_Upper_%d_%s' % (180, 205, "Pom152_859"))
    yax.set_label('Lower_%d_Upper_%d_%s' % (-205, -180, "Pom152_859"))
    yax.set_weight(yaxial_weight)
    yax.add_to_model()
    outputobjects.append(yax)
    print (yax.get_output())

    pom152_min = 5;     pom152_max = 55;
    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (379,379,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_379"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (520,520,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_520"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (616,616,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_616"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (722,722,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_722"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (824,824,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_824"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (931,931,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_931"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (1036,1036,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_1036"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (1150,1150,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_1150"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (1244,1244,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_1244"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())
    """

# ALPS Motifs
if (is_nucleoplasm):
    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (1,32,'Nup1'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup1')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (27,47,'Nup60.1'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup60.1')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (27,47,'Nup60.2'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup60.2')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

# ALPS Motifs
if (is_inner_ring):
    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (310,334,'Nup157'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup157')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    yax = IMP.npc.npc_restraints.YAxialPositionRestraint(simo, (310,334,'Nup157'), lower_bound=-60, upper_bound=-10, consider_radius=False, sigma=1.0, term='M')
    yax.set_label('Lower_%d_Upper_%d_%s' % (-60, -10, "Nup157_ALPS"))
    yax.set_weight(yaxial_weight)
    yax.add_to_model()
    outputobjects.append(yax)
    print (yax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (310,334,'Nup157'), lower_bound=-15, upper_bound=35, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (-15, 35, "Nup157_ALPS"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())


    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (320,344,'Nup170'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup170')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    yax = IMP.npc.npc_restraints.YAxialPositionRestraint(simo, (320,344,'Nup170'), lower_bound=-170, upper_bound=-120, consider_radius=False, sigma=1.0, term='M')
    yax.set_label('Lower_%d_Upper_%d_%s' % (-170, -120, "Nup170_ALPS"))
    yax.set_weight(yaxial_weight)
    yax.add_to_model()
    outputobjects.append(yax)
    print (yax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (320,344,'Nup170'), lower_bound=-60, upper_bound=-10, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (-60, -10, "Nup170_ALPS"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

# ALPS Motifs
if (is_n84):
    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (252,270,'Nup133'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup133')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationConditionalRestraint(simo, protein1=(135,152,'Nup120'), protein2=(197,216,'Nup120'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup120_135-152_197-216')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())


#####################################################
# Restraints setup - Membrane Exclusion
#####################################################
tor_th      = 150.0 - tor_th_ALPS
tor_R       = 390.0 + 150.0
tor_r       = tor_th/2.0
mex_sigma   = 0.2
mex_weight  = 1000.0

if (use_MembraneExclusion):
    nup_list = [entry[0] for entry in domains if not '@' in entry[0]]
    nup_list_unique = sorted(list(set(nup_list)))   # Make a unique list

    MEX_LIST = [
        [0, 0, "Nup84"],
        #[0, 0, "Nup85"],
        [0, 0, "Nup120"],
        [0, 0, "Nup133"],
        [0, 0, "Nup145c"],
        #[0, 0, "Seh1"],
        #[0, 0, "Sec13"],

        #[0, 0, "Nup82.1"],
        #[0, 0, "Nup82.2"],
        #[0, 0, "Nsp1.1"],
        #[0, 0, "Nsp1.2"],
        #[0, 0, "Nup159.1"],
        #[0, 0, "Nup159.2"],
        [1, 965, "Nup116.1"],
        [1, 965, "Nup116.2"],
        #[0, 0, "Gle1"],
        #[0, 0, "Nup42"],

        #[0, 0, "Nic96.1"],
        #[0, 0, "Nic96.2"],
        #[0, 0, "Nsp1.3"],
        #[0, 0, "Nsp1.4"],
        #[0, 0, "Nup57.1"],
        #[0, 0, "Nup57.2"],
        #[0, 0, "Nup49.1"],
        #[0, 0, "Nup49.2"],

        [1, 110, "Pom152"],
        [248, 655, "Ndc1"],
        [151, 299, "Pom34"],
        [0, 0, "Nup157"],
        [0, 0, "Nup170"],
        [0, 0, "Nup53"],
        [0, 0, "Nup59"],
        #[0, 0, "Nup188"],
        #[0, 0, "Nup192"],

        [0, 0, "Nup1"],
        [0, 0, "Nup60.1"],
        [0, 0, "Nup60.2"],
        [0, 0, "Nup145.1"],
        [0, 0, "Nup145.2"],
        [0, 0, "Nup100.1"],
        [0, 0, "Nup100.2"],
        [0, 0, "Mlp1"],
        [0, 0, "Mlp2"],
    ]
    print "\nMembraneExclusionRestraint !!"
    for z in MEX_LIST:
        if (z[2] not in nup_list_unique):
            continue
        if (z[0] > 0):
            mex = IMP.npc.npc_restraints.MembraneExclusionRestraint(simo, (z[0], z[1], z[2]), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=mex_sigma, resolution = res_ev)
        else:
            mex = IMP.npc.npc_restraints.MembraneExclusionRestraint(simo, z[2], tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=mex_sigma, resolution = res_ev)
        mex.set_label('%s_mex_%d_%d' % (z[2], z[0], z[1]))
        mex.set_weight(mex_weight)
        mex.add_to_model()
        outputobjects.append(mex)
        print (mex.get_output())


#####################################################
# Restraints setup
# Cross-link restraints using the whole NPC DSS XL data
#####################################################
if (use_XL):
    columnmap = {}
    columnmap["Protein1"] = "Protein 1"
    columnmap["Protein2"] = "Protein 2"
    columnmap["Residue1"] = "Residue 1"
    columnmap["Residue2"] = "Residue 2"
    columnmap["IDScore"] = "p value"
    columnmap["XLUniqueID"] = "XLUniqueID"
    ids_map = IMP.pmi.tools.map()
    ids_map.set_map_element(1.0, 1.0)

    xl1 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(simo,
                                                        '../data_npc/XL_optimized_ambiguity.csv',
                                                        length = 26.0,
                                                        slope = 0.00,
                                                        columnmapping = columnmap,
                                                        ids_map = ids_map,
                                                        resolution = 1.0,
                                                        #inner_slope = 0.02,
                                                        inner_slope = 0.01,
                                                        filelabel = "wtDSS",
                                                        label = "wtDSS",
                                                        attributes_for_label = ["XLUniqueID"],
                                                        csvfile = True)
    xl1.add_to_model()
    xl1.set_weight(10.0)        # play with the weight
    sampleobjects.append(xl1)
    outputobjects.append(xl1)
    xl1.set_psi_is_sampled(False)
    psi2 = xl1.get_psi(1.0)[0]
    psi2.set_scale(0.05)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
    print "\nEVAL 1 : ", sf.evaluate(False), " (after applying the XL restraint) - ", rank
    XL_restraints = [xl1]
else:
    XL_restraints = None


#####################################################
# Restraints setup
# Sampling Boundary Restraint for the outer ring
#####################################################
if (use_sampling_boundary):
    main_spoke = [];  other_spokes = [];    main_spoke_hier_name = []
    for entry in domains:
        if '@' in entry[0]:
            other_spokes.append(entry[0])
        elif ('Nup84' in entry[0]) or ('Nup85' in entry[0]) or ('Nup120' in entry[0]) or ('Nup133' in entry[0]):
            main_spoke.append(entry[0])
            main_spoke_hier_name.append(entry[1])
        elif ('Nup145c' in entry[0]) or ('Seh1' in entry[0]) or ('Sec13' in entry[0]):
            main_spoke.append(entry[0])
            main_spoke_hier_name.append(entry[1])
        elif ('Nup82' in entry[0]) or ('Nup159' in entry[0]) or ('Nsp1.1' in entry[0]) or ('Nsp1.2' in entry[0]):
            main_spoke.append(entry[0])
            main_spoke_hier_name.append(entry[1])
        else:
            other_spokes.append(entry[0])
    main_spoke_unique = sorted(list(set(main_spoke)))
    main_spoke_hier_name = sorted(main_spoke_hier_name)
    other_spokes_unique = sorted(list(set(other_spokes)))
    print ("main_spoke_hier_name = ", main_spoke_hier_name)
    print ("main_spoke_unique = ", main_spoke_unique)
    print ("other_spokes_unique = ", other_spokes_unique)

    #resdensities = bm1.get_density_hierarchies([t[1] for t in domains])
    resdensities = bm1.get_density_hierarchies(main_spoke_hier_name)
    print ("resdensities=", resdensities)       ####  TODO: make sure resdensities are correct

    mass = sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))
    #mass *= 1.2 * 2.0           # 1.2 for adjustment of the GMM (after removing flexible GMMs) and 2.0 for approximation of the NPC spoke mass
    mass *= 1.2           # 1.2 for adjustment of the GMM (after removing flexible GMMs)
    print ("Total mass for the Sampling Boundary EM restraint for the outer-ring = ", mass)
    sbr = IMP.pmi.restraints.em.GaussianEMRestraint(resdensities,
                                                    '../data_npc/em_gmm_model/SJ_outer_ring_newnewEM.gmm.200.txt',
                                                    target_mass_scale=mass,
                                                    #slope=0.01,
                                                    slope=0.0000001,
                                                    target_radii_scale=3.0)
    sbr.add_to_model()
    sbr.set_weight(10000.0)        # play with the weight
    sbr.set_label("Sampling_Boundary_outer-ring")
    #sbr.center_model_on_target_density(simo)
    outputobjects.append(sbr)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
    print "\nEVAL 2 : ", sf.evaluate(False), " (after applying the Sampling Boundary EM restraint) - ", rank


"""
#####################################################
# 1st Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "\nEVAL 3 : ", sf.evaluate(False), " (initial) - ", rank

simo.optimize_floppy_bodies(300)
print "\nEVAL 4 : ", sf.evaluate(False), " (after calling optimize_floppy_bodies(300)) - ", rank

mc1 = IMP.pmi.macros.ReplicaExchange0(m,
                                    simo,
                                    monte_carlo_sample_objects = sampleobjects,
                                    output_objects = outputobjects,
                                    crosslink_restraints = XL_restraints,
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 5.0,
                                    number_of_best_scoring_models = 0,
                                    monte_carlo_steps = 10,
                                    number_of_frames = 500,
                                    write_initial_rmf = True,
                                    initial_rmf_name_suffix = "initial",
                                    stat_file_name_suffix = "stat",
                                    best_pdb_name_suffix = "model",
                                    do_clean_first = True,
                                    do_create_directories = True,
                                    global_output_directory = "1_pre_output",
                                    rmf_dir = "rmfs/",
                                    best_pdb_dir = "pdbs/",
                                    replica_stat_file_suffix = "stat_replica")
mc1.execute_macro()
rex1 = mc1.get_replica_exchange_object()
print "\nEVAL 5 : ", sf.evaluate(False), " (after performing the pre_sampling) - ", rank
#exit(0)
"""

#####################################################
# Restraints setup
# EM 3D restraint using GMM
#####################################################
if (use_EM3D):
    main_spoke = [];  other_spokes = [];    main_spoke_hier_name = []
    for entry in domains:
        # ignore GMMs of the Nup84 complex in nucleoplasm
        if ( ('Nup84' in entry[0]) and ('@11' in entry[0]) ):
            other_spokes.append(entry[0])
        elif ( ('Nup85' in entry[0]) and ('@11' in entry[0]) ):
            other_spokes.append(entry[0])
        elif ( ('Nup120' in entry[0]) and ('@11' in entry[0]) ):
            other_spokes.append(entry[0])
        elif ( ('Nup133' in entry[0]) and ('@11' in entry[0]) ):
            other_spokes.append(entry[0])
        elif ( ('Nup145c' in entry[0]) and ('@11' in entry[0]) ):
            other_spokes.append(entry[0])
        elif ( ('Seh1' in entry[0]) and ('@11' in entry[0]) ):
            other_spokes.append(entry[0])
        elif ( ('Sec13' in entry[0]) and ('@11' in entry[0]) ):
            other_spokes.append(entry[0])
        elif 'Dyn2' in entry[0]:
            other_spokes.append(entry[0])
        elif '@11' in entry[0]:
            main_spoke.append(entry[0])
            main_spoke_hier_name.append(entry[1])
        elif '@' in entry[0]:
            other_spokes.append(entry[0])
        else:
            main_spoke.append(entry[0])
            main_spoke_hier_name.append(entry[1])
    main_spoke_unique = sorted(list(set(main_spoke)))
    main_spoke_hier_name = sorted(main_spoke_hier_name)
    other_spokes_unique = sorted(list(set(other_spokes)))
    print ("main_spoke_hier_name = ", main_spoke_hier_name)
    print ("main_spoke_unique = ", main_spoke_unique)
    print ("other_spokes_unique = ", other_spokes_unique)

    #resdensities = bm1.get_density_hierarchies([t[1] for t in domains])
    resdensities = bm1.get_density_hierarchies(main_spoke_hier_name)
    print ("resdensities=", resdensities)       ####  TODO: make sure resdensities are correct

    mass = sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))
    mass *= 1.2 * 2.0           # 1.2 for adjustment of the GMM (after removing flexible GMMs) and 2.0 for approximation of the NPC spoke mass
    print ("Total mass for the EM restraint = ", mass)
    gem = IMP.pmi.restraints.em.GaussianEMRestraint(resdensities,
                                                    '../data_npc/em_gmm_model/SJ_cropped_r_09_02_resfilt143_rotated_adjusted90.gmm.1750.txt',
                                                    target_mass_scale=mass,
                                                    #slope=0.0000005,
                                                    slope=0.0000001,
                                                    target_radii_scale=3.0)
    gem.add_to_model()
    gem.set_weight(3571.4)        # play with the weight
    #gem.center_model_on_target_density(simo)
    outputobjects.append(gem)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
    print "\nEVAL 6 : ", sf.evaluate(False), " (after applying the EM 3D restraint) - ", rank

"""
#####################################################
# 2nd Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
mc2 = IMP.pmi.macros.ReplicaExchange0(m,
                                    simo,
                                    monte_carlo_sample_objects = sampleobjects,
                                    output_objects = outputobjects,
                                    crosslink_restraints = XL_restraints,
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 5.0,
                                    number_of_best_scoring_models = 0,
                                    monte_carlo_steps = 10,
                                    number_of_frames = 1000,
                                    write_initial_rmf = True,
                                    initial_rmf_name_suffix = "initial",
                                    stat_file_name_suffix = "stat",
                                    best_pdb_name_suffix = "model",
                                    do_clean_first = True,
                                    do_create_directories = True,
                                    global_output_directory = "2_XL_EM_output",
                                    rmf_dir = "rmfs/",
                                    best_pdb_dir = "pdbs/",
                                    #replica_stat_file_suffix = "stat_replica")
                                    replica_stat_file_suffix = "stat_replica",
                                    replica_exchange_object = rex1)
mc2.execute_macro()
rex2 = mc2.get_replica_exchange_object()
print "\nEVAL 7 : ", sf.evaluate(False), " (after performing the XL_EM_sampling) - ", rank
#exit(0)
"""

#####################################################
# Restraints setup
# Excluded Volume restraint for components in the main spoke
#####################################################
if (use_ExcludedVolume):
    main_spoke = [];  other_spokes = [];    main_spoke_hier_name = []
    for entry in domains:
        if '@' in entry[0]:
            other_spokes.append(entry[0])
        else:
            main_spoke.append(entry[0])
            main_spoke_hier_name.append(entry[1])
    main_spoke_unique = sorted(list(set(main_spoke)))
    main_spoke_hier_name = sorted(main_spoke_hier_name)
    other_spokes_unique = sorted(list(set(other_spokes)))
    print ("main_spoke_hier_name = ", main_spoke_hier_name)
    print ("main_spoke_unique = ", main_spoke_unique)
    print ("other_spokes_unique = ", other_spokes_unique)

    included_objects = [];  other_objects = []
    for entry in main_spoke_unique:
        obj = simo.hier_dict[entry]
        included_objects.append(obj)
        #other_objects.append(obj)
    for entry in other_spokes_unique:
        other_objects.append(simo.hier_dict[entry])
    print ("EV included_objects in the main spoke = ", included_objects)
    print ("EV other_objects = ", other_objects)
    print ("resolution for EV = ", res_ev)

    ev1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo,
                                                                 included_objects = included_objects,
                                                                 #other_objects = other_objects,
                                                                 resolution = res_ev)
    ev1.add_to_model()
    ev1.set_label('main_spoke')
    ev1.set_weight(1.0)
    outputobjects.append(ev1)
    print(ev1.get_output())
    print "ExcludedVolumeSphere1 for the main spoke !!\n"

    if (use_neighboring_spokes):
        ev2 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo,
                                                                     included_objects = included_objects,
                                                                     other_objects = other_objects,
                                                                     resolution = res_ev)
        ev2.add_to_model()
        ev2.set_label('bipartite')
        ev2.set_weight(10.0)
        outputobjects.append(ev2)
        print(ev2.get_output())
        print "ExcludedVolumeSphere2 between the main spoke and the neighboring spokes !!\n"

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
    print "\nEVAL 8 : ", sf.evaluate(False), " (after applying the Excluded Volume restraint) - ", rank

"""
#####################################################
# 3rd Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
mc3 = IMP.pmi.macros.ReplicaExchange0(m,
                                    simo,
                                    monte_carlo_sample_objects = sampleobjects,
                                    output_objects = outputobjects,
                                    crosslink_restraints = XL_restraints,
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 5.0,
                                    number_of_best_scoring_models = 0,
                                    monte_carlo_steps = 10,
                                    number_of_frames = 2000,
                                    write_initial_rmf = True,
                                    initial_rmf_name_suffix = "initial",
                                    stat_file_name_suffix = "stat",
                                    best_pdb_name_suffix = "model",
                                    do_clean_first = True,
                                    do_create_directories = True,
                                    global_output_directory = "3_EV_output",
                                    rmf_dir = "rmfs/",
                                    best_pdb_dir = "pdbs/",
                                    replica_stat_file_suffix = "stat_replica",
                                    replica_exchange_object = rex2)
mc3.execute_macro()
rex3 = mc3.get_replica_exchange_object()
print "\nEVAL 9 : ", sf.evaluate(False), " (after performing the EV_sampling) - ", rank
"""


#####################################################
# 4th Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
mc4 = IMP.pmi.macros.ReplicaExchange0(m,
                                    simo,
                                    monte_carlo_sample_objects = sampleobjects,
                                    output_objects = outputobjects,
                                    crosslink_restraints = XL_restraints,
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 5.0,
                                    number_of_best_scoring_models = 0,
                                    monte_carlo_steps = 10,
                                    number_of_frames = int(inputs.nrepeats),
                                    write_initial_rmf = True,
                                    initial_rmf_name_suffix = "initial",
                                    stat_file_name_suffix = "stat",
                                    best_pdb_name_suffix = "model",
                                    do_clean_first = True,
                                    do_create_directories = True,
                                    global_output_directory = inputs.folder_output,
                                    rmf_dir = "rmfs/",
                                    best_pdb_dir = "pdbs/",
                                    replica_stat_file_suffix = "stat_replica")
                                    #replica_stat_file_suffix = "stat_replica",
                                    #replica_exchange_object = rex2)
                                    #replica_exchange_object = rex3)
mc4.execute_macro()
print "\nEVAL 10 : ", sf.evaluate(False), " (final evaluation) - ", rank
#exit(0)
