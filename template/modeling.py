#!/usr/bin/env python
#####################################################
# Last Update: September 10th, 2015
# by Seung Joong Kim and Riccardo Pellarin
# at Andrej Sali group, University of California San Francisco (UCSF)
#####################################################
import IMP
import IMP.core
#import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

#import crosslinking_nup82
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

#----------------------------------------------------
# Setting up the input parameters
#----------------------------------------------------
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
# setting up topology
#####################################################
m = IMP.Model()
#s = IMP.pmi.topology.System(m)
#st = s.create_state()
#simo = IMP.npc.npc_restraints.Representation(m,upperharmonic=True,disorderedlength=False)
simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)


#####################################################
# setting up parameters
#####################################################
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
rbmaxtrans = 4.00
rbmaxrot = 0.04

# flexible bead movement
fbmaxtrans = 4.00
outputobjects = []
sampleobjects = []

res_cry = int(inputs.res_cry)
res_ev = int(inputs.res_ev)
res_conn = int(inputs.res_cry)
beadsize = 10
beadsize20 = 20
beadsize25 = 25
beadsize100 = 100
em2d_weight = float(inputs.weight)
initial_nframes = 100

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

#n84_pdb        = n82 + 'rmfs/B_8_1-95r.pdb'         # for GMM generation
#n82_pdb        = n82 + 'rmfs/B_8_1-95r.pdb'         # for GMM generation
#n96_pdb        = n96 + 'rmfs/cluster2r.pdb'         # for GMM generation
n84_pdb        = n82 + 'rmfs/B_8_1-95_noBEA.pdb'
n82_pdb        = n82 + 'rmfs/B_8_1-95_noBEA.pdb'
n96_pdb        = n96 + 'rmfs/Cluster2_noBEA.pdb'

n157N_pdbfile  = npc + "Nup157_4mhc_A_88_892.pdb"
n157C_pdbfile  = npc + "Nup157_3i5p_A_900_1390.pdb"
n170N_pdbfile  = npc + "Nup170_4mhc_A_98_992.pdb"
n170C_pdbfile  = npc + "Nup170_3i5p_A_1000_1502.pdb"
n188_pdbfile   = npc + "Nup188_12_1652.pdb"
n192_pdbfile   = npc + "Nup192_1_1683_v2_new.pdb"

n53_pdbfile    = npc + "Nup53_dimer.pdb"
n59_pdbfile    = npc + "Nup59_dimer.pdb"
#pom152_pdb     = npc + "Pom152r.pdb"               # for GMM generation
pom152_pdb     = npc + "Pom152_noBEA.pdb"
n100_pdbfile   = npc + "Nup100_3nf5_AB_816_958.pdb"
n145N_pdbfile  = npc + "Nup145N_3kep_AB_459_605r.pdb"

#Dbp5_pdbfile   = npc + "Dbp5_3rrm_A_91_482.pdb"    # Not using it for now
Gle1N_pdbfile  = npc + "Gle1_4wij_A_121_239.pdb"
Gle1C_pdbfile  = npc + "Gle1_3rrm_B_244_538.pdb"
Gle2_pdbfile   = npc + "Gle2_3mmy_A_4_362.pdb"

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
is_basket = False
is_FG = False

use_neighboring_spokes = True
#Stopwatch_None_delta_seconds  32~38 (1 spoke) / 100~110 sec (3 spokes)
#Stopwatch_None_delta_seconds  34~40 (1 spoke) / 105~115 sec (3 spokes) with XL
#Stopwatch_None_delta_seconds  51~57 (1 spoke) / 120~130 sec (3 spokes) with XL + EM
use_shuffle = True
use_Distance_to_Point = True
use_Immuno_EM = True
use_Composite = False
use_XL = True
use_EM3D = True

#####################################################
# REPRESENTATION
#####################################################
# comp_name, hier_name, color, fasta_file, fasta_id, pdb_name, chain_id, res_range, read_em_files, bead_size, rb, super_rb, em_num_components, em_txt_file_name, em_mrc_file_name, chain_of_super_rb, keep_gaussian_flexible_beads
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
    n84_rb = 84
    domains.append(('Nup84',  "Nup84",    0.0,  n84_fastafile,   "Nup84",   n84_pdb, "K",  (  1, 726,0),  gmm,  beadsize,  n84_rb, [n84_rb], 3,  " ",   " ",  None, False))
    domains.append(('Nup85',  "Nup85_1",  0.2,  n85_fastafile,   "Nup85",   n84_pdb, "L",  (  1, 492,0),  gmm,  beadsize,  n84_rb, [n84_rb], 3,  " ",   " ",  None, False))
    domains.append(('Nup85',  "Nup85_2",  0.25, n85_fastafile,   "Nup85",   n84_pdb, "L",  (493, 744,0),  gmm,  beadsize,  n84_rb, [n84_rb], 2,  " ",   " ",  None, False))
    domains.append(('Nup120', "Nup120_1", 0.35, n120_fastafile,  "Nup120",  n84_pdb, "M",  (  1, 711,0),  gmm,  beadsize,  n84_rb, [n84_rb], 3,  " ",   " ",  None, False))
    domains.append(('Nup120', "Nup120_2", 0.4,  n120_fastafile,  "Nup120",  n84_pdb, "M",  (712,1037,0),  gmm,  beadsize,  n84_rb, [n84_rb], 2,  " ",   " ",  None, False))
    domains.append(('Nup133', "Nup133_1", 0.5,  n133_fastafile,  "Nup133",  n84_pdb, "N",  (  1, 480,0),  gmm,  beadsize,  n84_rb, [n84_rb], 2,  " ",   " ",  None, False))
    domains.append(('Nup133', "Nup133_2", 0.55, n133_fastafile,  "Nup133",  n84_pdb, "N",  (481,1157,0),  gmm,  beadsize,  n84_rb, [n84_rb], 3,  " ",   " ",  None, False))
    domains.append(('Nup145c',"Nup145c_1",0.65, n145c_fastafile, "Nup145c", n84_pdb, "O",  (  1, 125,0),  gmm,  beadsize,  n84_rb, [n84_rb], 1,  " ",   " ",  None, False))
    domains.append(('Nup145c',"Nup145c_2",0.7,  n145c_fastafile, "Nup145c", n84_pdb, "O",  (126, 712,0),  gmm,  beadsize,  n84_rb, [n84_rb], 3,  " ",   " ",  None, False))
    domains.append(('Seh1',   "Seh1",     0.8,  seh1_fastafile,  "Seh1",    n84_pdb, "P",  (  1, 349,0),  gmm,  beadsize,  n84_rb, [n84_rb], 2,  " ",   " ",  None, False))
    domains.append(('Sec13',  "Sec13",    0.95, sec13_fastafile, "Sec13",   n84_pdb, "Q",  (  1, 297,0),  gmm,  beadsize,  n84_rb, [n84_rb], 2,  " ",   " ",  None, False))
    for i in clones_range_A:
        if (i==11): gmm_c = gmm
        else:       gmm_c = None
        domains.append(('Nup84@%d'%i,  "Nup84@%d"%i,    0.0, n84_fastafile,   "Nup84",   n84_pdb, "K",  (  1, 726,0),  gmm_c,  beadsize,   None,   None,   3,  gmm_f+"Nup84.txt",     gmm_f+"Nup84.mrc",     None, False))
        domains.append(('Nup85@%d'%i,  "Nup85_1@%d"%i,  0.2, n85_fastafile,   "Nup85",   n84_pdb, "L",  (  1, 492,0),  gmm_c,  beadsize,   None,   None,   3,  gmm_f+"Nup85_1.txt",   gmm_f+"Nup85_1.mrc",   None, False))
        domains.append(('Nup85@%d'%i,  "Nup85_2@%d"%i,  0.25,n85_fastafile,   "Nup85",   n84_pdb, "L",  (493, 744,0),  gmm_c,  beadsize,   None,   None,   2,  gmm_f+"Nup85_2.txt",   gmm_f+"Nup85_2.mrc",   None, False))
        domains.append(('Nup120@%d'%i, "Nup120_1@%d"%i, 0.35,n120_fastafile,  "Nup120",  n84_pdb, "M",  (  1, 711,0),  gmm_c,  beadsize,   None,   None,   3,  gmm_f+"Nup120_1.txt",  gmm_f+"Nup120_1.mrc",  None, False))
        domains.append(('Nup120@%d'%i, "Nup120_2@%d"%i, 0.4, n120_fastafile,  "Nup120",  n84_pdb, "M",  (712,1037,0),  gmm_c,  beadsize,   None,   None,   2,  gmm_f+"Nup120_2.txt",  gmm_f+"Nup120_2.mrc",  None, False))
        domains.append(('Nup133@%d'%i, "Nup133_1@%d"%i, 0.5, n133_fastafile,  "Nup133",  n84_pdb, "N",  (  1, 480,0),  gmm_c,  beadsize,   None,   None,   2,  gmm_f+"Nup133_1.txt",  gmm_f+"Nup133_1.mrc",  None, False))
        domains.append(('Nup133@%d'%i, "Nup133_2@%d"%i, 0.55,n133_fastafile,  "Nup133",  n84_pdb, "N",  (481,1157,0),  gmm_c,  beadsize,   None,   None,   3,  gmm_f+"Nup133_2.txt",  gmm_f+"Nup133_2.mrc",  None, False))
        domains.append(('Nup145c@%d'%i,"Nup145c_1@%d"%i,0.65,n145c_fastafile, "Nup145c", n84_pdb, "O",  (  1, 125,0),  gmm_c,  beadsize,   None,   None,   1,  gmm_f+"Nup145c_1.txt", gmm_f+"Nup145c_1.mrc", None, False))
        domains.append(('Nup145c@%d'%i,"Nup145c_2@%d"%i,0.7, n145c_fastafile, "Nup145c", n84_pdb, "O",  (126, 712,0),  gmm_c,  beadsize,   None,   None,   3,  gmm_f+"Nup145c_2.txt", gmm_f+"Nup145c_2.mrc", None, False))
        domains.append(('Seh1@%d'%i,   "Seh1@%d"%i,     0.8, seh1_fastafile,  "Seh1",    n84_pdb, "P",  (  1, 349,0),  gmm_c,  beadsize,   None,   None,   2,  gmm_f+"Seh1.txt",      gmm_f+"Seh1.mrc",      None, False))
        domains.append(('Sec13@%d'%i,  "Sec13@%d"%i,    0.95,sec13_fastafile, "Sec13",   n84_pdb, "Q",  (  1, 297,0),  gmm_c,  beadsize,   None,   None,   2,  gmm_f+"Sec13.txt",     gmm_f+"Sec13.mrc",     None, False))

##########################
# Nup82 complex
##########################
if (is_n82):
    n82_rb = 84
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
    domains.append(("Nup116.1","Nup116.1",    0.75,  f_n82+"Nup116.txt", "Nup116", n82_pdb,  "I", ( 751,1113,0),  gmm,   beadsize25, n82_rb, [n82_rb], 1,  " ",   " ",  None, False))
    domains.append(("Nup116.2","Nup116.2",    0.8,   f_n82+"Nup116.txt", "Nup116", n82_pdb,  "J", ( 751,1113,0),  gmm,   beadsize25, n82_rb, [n82_rb], 1,  " ",   " ",  None, False))
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
        domains.append(("Nup116.1@%d"%i,"Nup116.1@%d"%i,    0.75,  f_n82+"Nup116.txt", "Nup116", n82_pdb,  "I", ( 751,1113,0),  None, beadsize25,  None,  None,  1,  None,  None,  None))
        domains.append(("Nup116.2@%d"%i,"Nup116.2@%d"%i,    0.8,   f_n82+"Nup116.txt", "Nup116", n82_pdb,  "J", ( 751,1113,0),  None, beadsize25,  None,  None,  1,  None,  None,  None))

##########################
# Nic96 complex
##########################
if (is_nic96):
    n96_rb = 196
    domains.append(("Nic96.1",  "Nic96.1_1",  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (  1,106,0),  gmm,  beadsize,    n96_rb, [n96_rb], 2,  " ",   " ",  None))
    domains.append(("Nic96.1",  "Nic96.1_2",  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (107,204,0),  gmm,  beadsize,    1097,   [n96_rb], 1,  " ",   " ",  None))
    domains.append(("Nic96.1",  "Nic96.1",    0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (205,839,0),  gmm,  beadsize,    1096,   [n96_rb], 3,  " ",   " ",  None, False))
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
    domains.append(("Nic96.2",  "Nic96.2_1",  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (  1,106,0),  gmm,  beadsize,    n96_rb, [n96_rb], 2,  gmm_f+"Nic96.1_1.txt", gmm_f+"Nic96.1_1.mrc", None))
    domains.append(("Nic96.2",  "Nic96.2_2",  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (107,204,0),  gmm,  beadsize,    2097,   [n96_rb], 1,  gmm_f+"Nic96.1_2.txt", gmm_f+"Nic96.1_2.mrc", None))
    domains.append(("Nic96.2",  "Nic96.2",    0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (205,839,0),  gmm,  beadsize,    2096,   [n96_rb], 3,  gmm_f+"Nic96.1.txt",   gmm_f+"Nic96.1.mrc",   None, False))
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
        domains.append(("Nic96.1@%d"%i,  "Nic96.1_1@%d"%i,  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (  1,106,0),  gmm_c,  beadsize,    None, None, 2,  gmm_f+"Nic96.1_1.txt", gmm_f+"Nic96.1_1.mrc", None))
        domains.append(("Nic96.1@%d"%i,  "Nic96.1_2@%d"%i,  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (107,204,0),  gmm_c,  beadsize,    None, None, 1,  gmm_f+"Nic96.1_2.txt", gmm_f+"Nic96.1_2.mrc", None))
        domains.append(("Nic96.1@%d"%i,  "Nic96.1@%d"%i,    0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (205,839,0),  gmm_c,  beadsize,    None, None, 3,  gmm_f+"Nic96.1.txt",   gmm_f+"Nic96.1.mrc",   None, False))
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

        domains.append(("Nic96.2@%d"%i,  "Nic96.2_1@%d"%i,  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (  1,106,0),  gmm_c,  beadsize,    None, None, 2,  gmm_f+"Nic96.1_1.txt", gmm_f+"Nic96.1_1.mrc", None))
        domains.append(("Nic96.2@%d"%i,  "Nic96.2_2@%d"%i,  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (107,204,0),  gmm_c,  beadsize,    None, None, 1,  gmm_f+"Nic96.1_2.txt", gmm_f+"Nic96.1_2.mrc", None))
        domains.append(("Nic96.2@%d"%i,  "Nic96.2@%d"%i,    0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (205,839,0),  gmm_c,  beadsize,    None, None, 3,  gmm_f+"Nic96.1.txt",   gmm_f+"Nic96.1.mrc",   None, False))
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
    domains.append(("Nup157",    "Nup157n",       0.0,  f_npc+"Nup157.txt", "YER105C", n157N_pdbfile,  "A",  (  1, 896,0),  gmm,  beadsize25,  1571, [157], 3, " ",   " ",  None, False))
    domains.append(("Nup157",    "Nup157c",       0.0,  f_npc+"Nup157.txt", "YER105C", n157C_pdbfile,  "A",  (897,1391,0),  gmm,  beadsize25,  1572, [157], 3, " ",   " ",  None, False))
    domains.append(("Nup170",    "Nup170n",       0.25, f_npc+"Nup170.txt", "Nup170",  n170N_pdbfile,  "A",  (  1, 996,0),  gmm,  beadsize25,  1701, [170], 3, " ",   " ",  None, False))
    domains.append(("Nup170",    "Nup170c",       0.25, f_npc+"Nup170.txt", "Nup170",  n170C_pdbfile,  "A",  (997,1502,0),  gmm,  beadsize25,  1702, [170], 4, " ",   " ",  None, False))
    domains.append(("Nup188",    "Nup188",        0.5,  f_npc+"Nup188.txt", "YML103C", n188_pdbfile,   "A",  (  1,1655,0),  gmm,  beadsize25,  188,  [196], 6, " ",   " ",  None, False))
    domains.append(("Nup192",    "Nup192",        0.75, f_npc+"Nup192.txt", "YJL039C", n192_pdbfile,   "A",  (  1,1683,0),  gmm,  beadsize25,  192,  [296], 5, " ",   " ",  None, False))
    for i in clones_range_A:
        if (i==11): gmm_c = gmm
        else:       gmm_c = None
        domains.append(("Nup157@%d"%i,  "Nup157n@%d"%i,  0.0,  f_npc+"Nup157.txt", "YER105C", n157N_pdbfile,  "A",  (  1, 896,0),  gmm_c,  beadsize25,  None, None, 3, gmm_f+"Nup157n.txt", gmm_f+"Nup157n.mrc", None, False))
        domains.append(("Nup157@%d"%i,  "Nup157c@%d"%i,  0.0,  f_npc+"Nup157.txt", "YER105C", n157C_pdbfile,  "A",  (897,1391,0),  gmm_c,  beadsize25,  None, None, 3, gmm_f+"Nup157c.txt", gmm_f+"Nup157c.mrc", None, False))
        domains.append(("Nup170@%d"%i,  "Nup170n@%d"%i,  0.25, f_npc+"Nup170.txt", "Nup170",  n170N_pdbfile,  "A",  (  1, 996,0),  gmm_c,  beadsize25,  None, None, 3, gmm_f+"Nup170n.txt", gmm_f+"Nup170n.mrc", None, False))
        domains.append(("Nup170@%d"%i,  "Nup170c@%d"%i,  0.25, f_npc+"Nup170.txt", "Nup170",  n170C_pdbfile,  "A",  (997,1502,0),  gmm_c,  beadsize25,  None, None, 4, gmm_f+"Nup170c.txt", gmm_f+"Nup170c.mrc", None, False))
        domains.append(("Nup188@%d"%i,  "Nup188@%d"%i,   0.5,  f_npc+"Nup188.txt", "YML103C", n188_pdbfile,   "A",  (  1,1655,0),  gmm_c,  beadsize25,  None, None, 6, gmm_f+"Nup188.txt",  gmm_f+"Nup188.mrc",  None, False))
        domains.append(("Nup192@%d"%i,  "Nup192@%d"%i,   0.75, f_npc+"Nup192.txt", "YJL039C", n192_pdbfile,   "A",  (  1,1683,0),  gmm_c,  beadsize25,  None, None, 5, gmm_f+"Nup192.txt",  gmm_f+"Nup192.mrc",  None, False))

##########################
# Membrane nups
##########################
if (is_membrane):
    domains.append(("Nup53",     "Nup53",         0.0,  f_npc+"Nup53.txt",  "YMR153W", n53_pdbfile,  "A", (   1, 475,0),  gmm,  beadsize100, 53,   [53], 2,  " ",   " ",  None))
    domains.append(("Nup59",     "Nup59",         0.66, f_npc+"Nup59.txt",  "YDL088C", n59_pdbfile,  "A", (   1, 528,0),  gmm,  beadsize100, 59,   [59], 2,  " ",   " ",  None))
    domains.append(("Ndc1",      "Ndc1",          0.8,  f_npc+"Ndc1.txt",   "YML031W", "BEADS",      " ", (   1, 655,0),  gmm,  beadsize100, 101,  [101],0,  None,  None, None))
    domains.append(("Pom34",     "Pom34",         0.9,  f_npc+"Pom34.txt",  "YLR018C", "BEADS",      " ", (   1, 299,0),  gmm,  beadsize100, 34,   [34], 0,  None,  None, None))
    domains.append(("Pom152",    "Pom152_1" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", (   1, 496,0),  gmm,  beadsize100, 1521, [152],2,  " ",   " ",  None))
    domains.append(("Pom152",    "Pom152_2" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", ( 497, 613,0),  gmm,  beadsize100, 1522, [152],2,  " ",   " ",  None, False))
    domains.append(("Pom152",    "Pom152_3" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", ( 614, 718,0),  gmm,  beadsize100, 1523, [152],2,  " ",   " ",  None, False))
    domains.append(("Pom152",    "Pom152_4" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", ( 719, 821,0),  gmm,  beadsize100, 1524, [152],2,  " ",   " ",  None, False))
    domains.append(("Pom152",    "Pom152_5" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", ( 822, 924,0),  gmm,  beadsize100, 1525, [152],2,  " ",   " ",  None, False))
    domains.append(("Pom152",    "Pom152_6" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", ( 925,1031,0),  gmm,  beadsize100, 1526, [152],2,  " ",   " ",  None, False))
    domains.append(("Pom152",    "Pom152_7" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", (1032,1145,0),  gmm,  beadsize100, 1527, [152],2,  " ",   " ",  None, False))
    domains.append(("Pom152",    "Pom152_8" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", (1146,1236,0),  gmm,  beadsize100, 1528, [152],2,  " ",   " ",  None, False))
    domains.append(("Pom152",    "Pom152_9" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", (1237,1337,0),  gmm,  beadsize100, 1529, [152],2,  " ",   " ",  None, False))
    for i in clones_range_A:
        if (i==11): gmm_c = gmm
        else:       gmm_c = None
        domains.append(("Nup53@%d"%i,  "Nup53@%d"%i,    0.0, f_npc+"Nup53.txt",  "YMR153W", n53_pdbfile, "A", (   1, 475,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Nup53.txt",    gmm_f+"Nup53.mrc",    None))
        domains.append(("Nup59@%d"%i,  "Nup59@%d"%i,    0.66,f_npc+"Nup59.txt",  "YDL088C", n59_pdbfile, "A", (   1, 528,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Nup59.txt",    gmm_f+"Nup59.mrc",    None))
        domains.append(("Ndc1@%d"%i,   "Ndc1@%d"%i,     0.8, f_npc+"Ndc1.txt",   "YML031W", "BEADS",     " ", (   1, 655,0), gmm_c,  beadsize100, None, None, 0,  None,                 None,                 None))
        domains.append(("Pom34@%d"%i,  "Pom34@%d"%i,    0.9, f_npc+"Pom34.txt",  "YLR018C", "BEADS",     " ", (   1, 299,0), gmm_c,  beadsize100, None, None, 0,  None,                 None,                 None))
        domains.append(("Pom152@%d"%i, "Pom152_1@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", (   1, 496,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_1.txt", gmm_f+"Pom152_1.mrc", None))
        domains.append(("Pom152@%d"%i, "Pom152_2@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", ( 497, 613,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_2.txt", gmm_f+"Pom152_2.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_3@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", ( 614, 718,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_3.txt", gmm_f+"Pom152_3.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_4@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", ( 719, 821,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_4.txt", gmm_f+"Pom152_4.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_5@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", ( 822, 924,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_5.txt", gmm_f+"Pom152_5.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_6@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", ( 925,1031,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_6.txt", gmm_f+"Pom152_6.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_7@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", (1032,1145,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_7.txt", gmm_f+"Pom152_7.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_8@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", (1146,1236,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_8.txt", gmm_f+"Pom152_8.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_9@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", (1237,1337,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_9.txt", gmm_f+"Pom152_9.mrc", None, False))

##########################
# Cytoplasm only - TODO: Dbp5
##########################
if (is_cytoplasm):
    if (is_FG):
        domains.append(("Nup100.1", "Nup100.1_10", 0.2, f_npc+"Nup100.txt", "YKL068W", "BEADS",      " ", (  1, 550,0), None, beadsize100, 1100, None, 0,  None,  None, None))
        domains.append(("Nup100.2", "Nup100.2_10", 0.4, f_npc+"Nup100.txt", "YKL068W", "BEADS",      " ", (  1, 550,0), None, beadsize100, 2100, None, 0,  None,  None, None))
        domains.append(("Nup42",    "Nup42_10",    0.6, f_npc+"Nup42.txt",  "YDR192C", "BEADS",      " ", (  1, 363,0), None, beadsize100, 42,   None, 0,  None,  None, None))
    domains.append(("Nup100.1", "Nup100.1",    0.2, f_npc+"Nup100.txt", "YKL068W", n100_pdbfile, "A", (551, 959,0), gmm,  beadsize25,  1100, [1100],2,  " ",   " ",  None, False))
    domains.append(("Nup100.2", "Nup100.2",    0.4, f_npc+"Nup100.txt", "YKL068W", n100_pdbfile, "B", (551, 959,0), gmm,  beadsize25,  2100, [2100],2,  " ",   " ",  None, False))
    domains.append(("Nup42",    "Nup42",       0.6, f_npc+"Nup42.txt",  "YDR192C", "BEADS",      " ", (364, 430,0), gmm,  beadsize100, 42,   [42],  0,  None,  None, None))
    domains.append(("Gle1",     "Gle1_10",     0.8, f_npc+"Gle1.txt",   "YDL207W", Gle1N_pdbfile,"A", (  1, 239,0), gmm,  beadsize25,  610,  [611], 2,  " ",   " ",  None, False))
    domains.append(("Gle1",     "Gle1",        0.8, f_npc+"Gle1.txt",   "YDL207W", Gle1C_pdbfile,"B", (240, 538,0), gmm,  beadsize25,  611,  [611], 2,  " ",   " ",  None, False))
    domains.append(("Gle2.1",   "Gle2.1",      0.9, f_npc+"Gle2.txt",   "YER107C", Gle2_pdbfile, "A", (  1, 365,0), gmm,  beadsize25,  1612, [1612],2,  " ",                 " ",                None, False))
    domains.append(("Gle2.2",   "Gle2.2",      1.0, f_npc+"Gle2.txt",   "YER107C", Gle2_pdbfile, "A", (  1, 365,0), gmm,  beadsize25,  2612, [2612],2,  gmm_f+"Gle2.1.txt",  gmm_f+"Gle2.1.mrc", None, False))
    for i in clones_range_B:
        if (is_FG):
            domains.append(("Nup100.1@%d"%i, "Nup100.1_10@%d"%i, 0.2, f_npc+"Nup100.txt", "YKL068W", "BEADS",      " ", (  1, 550,0), None, beadsize100, None, None, 0,  None,  None, None))
            domains.append(("Nup100.2@%d"%i, "Nup100.2_10@%d"%i, 0.4, f_npc+"Nup100.txt", "YKL068W", "BEADS",      " ", (  1, 550,0), None, beadsize100, None, None, 0,  None,  None, None))
            domains.append(("Nup42@%d"%i,    "Nup42_10@%d"%i,    0.6, f_npc+"Nup42.txt",  "YDR192C", "BEADS",      " ", (  1, 363,0), None, beadsize100, None, None, 0,  None,  None, None))
        domains.append(("Nup100.1@%d"%i, "Nup100.1@%d"%i,    0.2, f_npc+"Nup100.txt", "YKL068W", n100_pdbfile, "A", (551, 959,0), None, beadsize25,  None, None, 2,  None,  None, None))
        domains.append(("Nup100.2@%d"%i, "Nup100.2@%d"%i,    0.4, f_npc+"Nup100.txt", "YKL068W", n100_pdbfile, "B", (551, 959,0), None, beadsize25,  None, None, 2,  None,  None, None))
        domains.append(("Nup42@%d"%i,    "Nup42@%d"%i,       0.6, f_npc+"Nup42.txt",  "YDR192C", "BEADS",      " ", (364, 430,0), None, beadsize100, None, None, 0,  None,  None, None))
        domains.append(("Gle1@%d"%i,     "Gle1_10@%d"%i,     0.8, f_npc+"Gle1.txt",   "YDL207W", Gle1N_pdbfile,"A", (  1, 239,0), None, beadsize25,  None, None, 2,  None,  None, None))
        domains.append(("Gle1@%d"%i,     "Gle1@%d"%i,        0.8, f_npc+"Gle1.txt",   "YDL207W", Gle1C_pdbfile,"B", (240, 538,0), None, beadsize25,  None, None, 2,  None,  None, None))
        domains.append(("Gle2.1@%d"%i,   "Gle2.1@%d"%i,      0.9, f_npc+"Gle2.txt",   "YER107C", Gle2_pdbfile, "A", (  1, 365,0), None, beadsize25,  None, None, 2,  None,  None, None))
        domains.append(("Gle2.2@%d"%i,   "Gle2.2@%d"%i,      1.0, f_npc+"Gle2.txt",   "YER107C", Gle2_pdbfile, "A", (  1, 365,0), None, beadsize25,  None, None, 2,  None,  None, None))

##########################
# Nucleoplasm only - TODO: Mlp1, Mlp2
##########################
if (is_nucleoplasm):
    if (is_FG):
        domains.append(("Nup145.1", "Nup145.1_10", 0.7, f_npc+"Nup145.txt", "YGL092W", "BEADS",      " ", (  1, 200,0), None, beadsize100, 1145, None, 0,  None,  None, None))
        domains.append(("Nup145.2", "Nup145.2_10", 0.75,f_npc+"Nup145.txt", "YGL092W", "BEADS",      " ", (  1, 200,0), None, beadsize100, 2145, None, 0,  None,  None, None))
        domains.append(("Nup1",     "Nup1_10",     0.8, f_npc+"Nup1.txt",   "YOR098C", "BEADS",      " ", (352,1076,0), None, beadsize100, 1,    None, 0,  None,  None, None))
        domains.append(("Nup60",    "Nup60_10",    0.9, f_npc+"Nup60.txt",  "YAR002W", "BEADS",      " ", (399, 539,0), None, beadsize100, 60,   None, 0,  None,  None, None))
    domains.append(("Nup145.1", "Nup145.1",    0.7, f_npc+"Nup145.txt", "YGL092W", n145N_pdbfile,"A", (201, 605,0), gmm,  beadsize25,  1145, None, 2,  " ",   " ",  None, False))
    domains.append(("Nup145.2", "Nup145.2",    0.75,f_npc+"Nup145.txt", "YGL092W", n145N_pdbfile,"B", (201, 605,0), gmm,  beadsize25,  2145, None, 2,  " ",   " ",  None, False))
    domains.append(("Nup1",     "Nup1",        0.8, f_npc+"Nup1.txt",   "YOR098C", "BEADS",      " ", (  1, 351,0), gmm,  beadsize100, 1,    None, 0,  None,  None, None))
    domains.append(("Nup60",    "Nup60",       0.9, f_npc+"Nup60.txt",  "YAR002W", "BEADS",      " ", (  1, 398,0), gmm,  beadsize100, 60,   None, 0,  None,  None, None))

    for i in clones_range_B:
        if (is_FG):
            domains.append(("Nup145.1@%d"%i, "Nup145.1_10@%d"%i, 0.7, f_npc+"Nup145.txt", "YGL092W", "BEADS",      " ", (  1, 200,0), None, beadsize100, None, None, 0,  None,  None, None))
            domains.append(("Nup145.2@%d"%i, "Nup145.2_10@%d"%i, 0.75,f_npc+"Nup145.txt", "YGL092W", "BEADS",      " ", (  1, 200,0), None, beadsize100, None, None, 0,  None,  None, None))
            domains.append(("Nup1@%d"%i,     "Nup1_10@%d"%i,     0.8, f_npc+"Nup1.txt",   "YOR098C", "BEADS",      " ", (352,1076,0), None, beadsize100, None, None, 0,  None,  None, None))
            domains.append(("Nup60@%d"%i,    "Nup60_10@%d"%i,    0.9, f_npc+"Nup60.txt",  "YAR002W", "BEADS",      " ", (399, 539,0), None, beadsize100, None, None, 0,  None,  None, None))
        domains.append(("Nup145.1@%d"%i, "Nup145.1@%d"%i,    0.7, f_npc+"Nup145.txt", "YGL092W", n145N_pdbfile,"A", (201, 605,0), None, beadsize25,  None, None, 2,  None,  None, None))
        domains.append(("Nup145.2@%d"%i, "Nup145.2@%d"%i,    0.75,f_npc+"Nup145.txt", "YGL092W", n145N_pdbfile,"B", (201, 605,0), None, beadsize25,  None, None, 2,  None,  None, None))
        domains.append(("Nup1@%d"%i,     "Nup1@%d"%i,        0.8, f_npc+"Nup1.txt",   "YOR098C", "BEADS",      " ", (  1, 351,0), None, beadsize100, None, None, 0,  None,  None, None))
        domains.append(("Nup60@%d"%i,    "Nup60@%d"%i,       0.9, f_npc+"Nup60.txt",  "YAR002W", "BEADS",      " ", (  1, 398,0), None, beadsize100, None, None, 0,  None,  None, None))

#####################################################
# Model Building
#####################################################
bm1 = IMP.pmi.macros.BuildModel1(simo)
bm1.set_gmm_models_directory(gmm_f)

if (is_n82):
    n82=['Dyn2.1', 'Dyn2.2', 'Nup82.1', 'Nup82.2', 'Nup159.1', 'Nup159.2', 'Nsp1.1', 'Nsp1.2', 'Nup116.1', 'Nup116.2']
    for d in list(n82):
        if (is_FG): bm1.set_rmf_file(d, "../data_nup82/rmfs/B_8_1-95.rmf3", 0)
        else:       bm1.set_rmf_file(d, "../data_nup82/rmfs/B_8_1-95_FGtruncated.rmf3", 0)

if (is_n84):
    n84=['Nup84', 'Nup85', 'Nup120', 'Nup133', 'Nup145c', 'Seh1', 'Sec13']
    for d in list(n84):
        bm1.set_rmf_file(d, "../data_nup82/rmfs/B_8_1-95.rmf3", 0)

if (is_nic96):
    nic96_dict={'Nic96.1':'Nic96', 'Nsp1.3':'Nsp1', 'Nup49.1':'Nup49', 'Nup57.1':'Nup57', 'Nic96.2':'Nic96', 'Nsp1.4':'Nsp1', 'Nup49.2':'Nup49', 'Nup57.2':'Nup57'}
    for key,val in nic96_dict.items():
        if (is_FG): bm1.set_rmf_file(key, "../data_nic96/rmfs/cluster2.rmf3", 0, rmf_component_name=val)
        else:       bm1.set_rmf_file(key, "../data_nic96/rmfs/cluster2_FGtruncated.rmf3", 0, rmf_component_name=val)
        #print("{} = {}".format(key, val))

# remove connectivity for clones
clone_list = [entry[0] for entry in domains if '@' in entry[0]]
clone_list_unique = sorted(list(set(clone_list)))   # Make a unique list
print ("clone_list_unique = ", clone_list_unique)

bm1.build_model(data_structure = domains, sequence_connectivity_scale=2.0, sequence_connectivity_resolution=1.0,
                skip_connectivity_these_domains=clone_list_unique, skip_gaussian_in_representation=use_EM3D)
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
        for protein in ['Nup100.1', 'Nup100.2', 'Nup42', 'Gle1', 'Gle2.1', 'Gle2.2']:
            simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
    if (is_nucleoplasm):
        #for protein in ['Nup145.1', 'Nup145.2', 'Nup60', 'Nup1', 'Mlp1', 'Mlp2']:
        for protein in ['Nup145.1', 'Nup145.2', 'Nup60', 'Nup1']:
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
for protein in ['Nup84', 'Nup85', 'Nup120', 'Nup133', 'Nup145c', 'Seh1', 'Sec13']:
    rigid_tuples.append(protein)
for protein in ['Dyn2.1', 'Dyn2.2', 'Nup82.1', 'Nup82.2']:
    rigid_tuples.append(protein)
for protein in [(1117,1460,'Nup159.1'),(1117,1460,'Nup159.2'), (637,823,'Nsp1.1'),(637,823,'Nsp1.2'), (966,1113,'Nup116.1'),(966,1113,'Nup116.2')]:
    rigid_tuples.append(protein)
for protein in [(1,56,'Nic96.1'),(205,839,'Nic96.1'), (637,823,'Nsp1.3'), (270,472,'Nup49.1'), (287,541,'Nup57.1')]:
    rigid_tuples.append(protein)
for protein in [(1,56,'Nic96.2'),(205,839,'Nic96.2'), (637,823,'Nsp1.4'), (270,472,'Nup49.2'), (287,541,'Nup57.2')]:
    rigid_tuples.append(protein)
"""
for protein in [(1,56,'Nic96.1'),(205,394,'Nic96.1'),(405,839,'Nic96.1'), (637,737,'Nsp1.3'),(742,823,'Nsp1.3'), (270,359,'Nup49.1'),(369,427,'Nup49.1'),(433,472,'Nup49.1'), (287,496,'Nup57.1'),(505,541,'Nup57.1')]:
    rigid_tuples.append(protein)
for protein in [(1,56,'Nic96.2'),(205,394,'Nic96.2'),(405,839,'Nic96.2'), (637,737,'Nsp1.4'),(742,823,'Nsp1.4'), (270,359,'Nup49.2'),(369,427,'Nup49.2'),(433,472,'Nup49.2'), (287,496,'Nup57.2'),(505,541,'Nup57.2')]:
    rigid_tuples.append(protein)
"""
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
    simo.shuffle_configuration(bounding_box=((250, -175, 50), (1250, 175, 750)), ignore_initial_coordinates=True, cutoff=1.0, niterations=1000)


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
# Restraints setup
# Excluded Volume restraint for components in the main spoke
#####################################################
main_spoke = [];  other_spokes = [];    main_spoke_hier_name = []
for entry in domains:
    if '@11' in entry[0]:
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
    outputobjects.append(ev2)
    print(ev2.get_output())
    print "ExcludedVolumeSphere2 between the main spoke and the neighboring spokes !!\n"


#####################################################
# Restraints setup
# External Barrier restraint
#####################################################
if (False):
    eb = IMP.pmi.restraints.basic.ExternalBarrier(simo, radius = 2000)
    eb.add_to_model()
    outputobjects.append(eb)
    print(eb.get_output())
    print "ExternalBarrier !!\n"


#####################################################
# Restraints setup
# Composite restraint (Proximity / Contact)
#####################################################
class CompositeRepresentation(object):
    """Simple representation for composites. For speed, each protein
       is represented with a single sphere that covers the 10 residue
       beads used for excluded volume. This is very approximate but
       that should be fine for composites. This should dramatically
       reduce the number of distances the composite restraint needs
       to calculate at each evaluation.

       The class acts like a dictionary where the keys are protein
       names and the values are the one-bead-per-protein particles,
       which are calculated on demand."""

    # All proteins that have multiple copies in the primary spoke
    # (proteins not explicitly listed here have only a single copy)
    COPIES = { 'Gle2': 3, 'Nic96': 3, 'Nsp1': 5, 'Nup100': 3,
               'Nup116': 3, 'Nup145': 3, 'Nup159': 3, 'Nup49': 3,
               'Nup57': 3, 'Nup82': 3, 'Dyn2': 3}

    def __init__(self, simo):
        self.simo = simo
        self.protein_beads = {}

    def __getitem__(self, protein):
        if protein not in self.protein_beads:
            # Make a single particle that covers the entire protein
            low_res = IMP.pmi.tools.select(self.simo, resolution=res_ev,
                                           name=protein)
            p = IMP.Particle(self.simo.m, "Sphere covering " + protein)
            c = IMP.core.Cover.setup_particle(p, low_res)
            self.protein_beads[protein] = c
        return self.protein_beads[protein]

    def get_all_copies(self, protein):
        """Return a list of particles for all copies of the given protein
           (e.g. r.get_all_copies('Gle2') should be equivalent to
           [r['Gle2'], r['Gle2.1'], r['Gle2.2']]"""
        copies = self.COPIES.get(protein, 1)
        return [self[protein]] + [self['%s.%d'] % i for i in range(1, copies)]

if use_Composite:
    COMPOSITE = {
        "C1"  : ("Gle2", "Nup116", 1),
        "C2"  : ("Nup82", "Nsp1", 1),
        "C3"  : ("Nsp1", "Nup57", 1),
        "C4"  : ("Nup188", "Nic96", 1),
        "C5"  : ("Nsp1", "Nic96", 1),
        "C6"  : ("Nup192", "Nic96", 1),
        "C7"  : ("Seh1", "Nup85", 1),
        "C8"  : ("Nup49", "Nup57", 1),
        "C9"  : ("Gle1", "Nup42", 1),
        "C10" : ("Nup42", "Gle1", 1),
        "C11" : ("Nsp1", "Nup49", 1),
        "C12" : ("Pom152", "Pom34", 1),
        "C13" : ("Nup53", "Nup170", 1),
        "C14" : ("Nup84", "Nup145C", 1),
        "C15" : ("Nup85", "Seh1", 1),
        "C16" : ("Nup57", "Nsp1", 1),
        "C17" : ("Nup116", "Gle2", 1),
        "C18" : ("Nup57", "Nup49", 1),
        "C19" : ("Nup82", "Nsp1", "Nup159", 2),
        "C20" : ("Nup84", "Nup145C", "Sec13", 2),
        "C21" : ("Nup192", "Nup60", "Pom152", 3),
        "C22" : ("Nup170", "Ndc1", "Pom152", 3),
        "C24" : ("Nup82", "Nsp1", "Nup159", 3),
        "C25" : ("Nup170", "Nup192", "Nup59", 2),
        "C26" : ("Nup84", "Nup145C", "Sec13", 1),
        "C27" : ("Nup82", "Nup116", "Gle2", 1),
        "C28" : ("Nsp1", "Nup57", "Nup49", "Nic96", 1),
        "C30" : ("Nup145C", "Nup84", "Sec13", "Nup133", 2),
        "C31" : ("Nup82", "Nsp1", "Nup159", 1),
        "C32" : ("Nup82", "Nsp1", "Nup116", "Nup159", 2),
        "C33" : ("Nup84", "Nup145C", "Nup120", "Nup85", 1),
        "C34" : ("Nsp1", "Nup82", "Nup159", "Nup116", "Nup82", 3),
        "C35" : ("Nup53", "Nup170", "Pom152", "Nic96", "Nup145C", 2),
        "C36" : ("Nup159", "Nup82", "Nsp1", "Gle2", "Nup116", 2),
        "C37" : ("Pom34", "Pom152", "Nup192", "Nup170", "Nup157", 2),
        "C38" : ("Nup57", "Nup49", "Nsp1", "Nic96", "Nup192", "Pom152", 1),
        "C39" : ("Nup84", "Nup145C", "Nup85", "Nup120", "Nup133", 1),
        "C40" : ("Pom152", "Nup192", "Nup170", "Ndc1", "Nup59", 3),
        "C41" : ("Pom34", "Pom152", "Nup170", "Ndc1", "Nup157", 2),
        "C42" : ("Seh1", "Nup120", "Nup145N", "Nup170", "Nup53", 3),
        "C43" : ("Nup42", "Gle1", "Nup82", "Nup100", "Nup116", "Nup159", 3),
        "C44" : ("Nup145C", "Nup84", "Sec13", "Nup85", "Seh1", "Nup120", 1),
        "C46" : ("Gle2", "Nup116", "Nup82", "Nsp1", "Nup159", "Nup84", 2),
        "C47" : ("Nup57", "Nup49", "Nsp1", "Nic96", "Nup192", "Pom152", 1),
        "C48" : ("Nup49", "Nup57", "Nsp1", "Nic96", "Nup192", "Pom152", 2),
        "C49" : ("Pom152", "Nup192", "Nup170", "Nup157", "Nup188", "Nup1", 3),
        "C50" : ("Nup57", "Nsp1", "Nup49", "Nic96", "Nup82", "Nup159", 1),
        "C51" : ("Nup145C", "Nup84", "Sec13", "Nup85", "Seh1", "Nup120", 1),
        "C52" : ("Nsp1", "Nup188", "Nup192", "Pom152", "Nup157", "Nup60", 3),
        "C53" : ("Nup120", "Seh1", "Nup85", "Nup84", "Nup145C", "Sec13", "Nup133", 1),
        "C55" : ("Nup53", "Nup170", "Nup157", "Nup192", "Nup188", "Nup145N",
                 "Pom152", 3),
        "C56" : ("Nup159", "Nup82", "Nsp1", "Nup188", "Nup192", "Nup170", "Pom152", 1),
        "C57" : ("Nup120", "Seh1", "Nup85", "Nup84", "Nup145C", "Sec13", "Nup133", 1),
        "C58" : ("Nup49", "Nsp1", "Nup57", "Nic96", "Nup82", "Nup192", "Nup159", 2),
        "C59" : ("Nup57", "Nsp1", "Nic96", "Nup192", "Pom152", "Nup82", "Nup84",
                 "Nup116", 3),
        "C60" : ("Nup145N", "Nup120", "Seh1", "Nup85", "Nup84", "Sec13", "Nup145C",
                 "Nup133", 3),
        "C61" : ("Nup159", "Nup82", "Nsp1", "Nup188", "Nup192", "Nup170", "Pom152", 1),
        "C62" : ("Nup57", "Nup49", "Nic96", "Nup170", "Nup157", "Nup59", "Pom152",
                 "Ndc1", 3),
        "C63" : ("Nup133", "Nup145C", "Nup84", "Sec13", "Nup85", "Seh1", "Nup120",
                 "Nup145N", "Nup157", 1),
        "C64" : ("Nup116", "Nup82", "Nsp1", "Nic96", "Nup192", "Pom152", "Nup170",
                 "Nup157", 3),
        "C65" : ("Nup57", "Nup49", "Nsp1", "Nic96", "Nup116", "Nup133", "Nup192",
                 "Nup170", "Nup157", 2),
        "C66" : ("Nup133", "Nup145C", "Nup84", "Sec13", "Nup85", "Seh1", "Nup120",
                 "Nup145N", "Nup157", 2),
        "C67" : ("Gle2", "Nup116", "Nup82", "Nup159", "Nsp1", "Nic96", "Nup188",
                 "Nup192", "Pom152", "Nup170", 1),
        "C68" : ("Sec13", "Nup84", "Nup145C", "Nup133", "Nup116", "Nup82", "Nup85",
                 "Seh1", "Nup120", "Nup145N", 2),
        "C69" : ("Nup53", "Nup170", "Nup192", "Nic96", "Nup188", "Nsp1", "Nup82",
                 "Nup159", "Nup116", "Nup84", 3),
        "C70" : ("Pom34", "Pom152", "Nup192", "Nup170", "Nic96", "Nup188", "Nsp1",
                 "Nup82", "Nup159", "Nup84", 2),
        "C71" : ("Nup145N", "Nup120", "Seh1", "Nup85", "Nup84", "Nup145C", "Sec13",
                 "Nup133", "Nic96", "Nup157", "Pom152", 3),
        "C72" : ("Pom152", "Nup192", "Nup170", "Ndc1", "Nup53", "Nic96", "Nup188",
                 "Nsp1", "Nup82", "Nup159", "Nup84", 2),
        "C73" : ("Gle2", "Nup116", "Nup82", "Nup159", "Nsp1", "Nic96", "Nup188",
                 "Nup192", "Nup170", "Pom152", "Nup157", 3),
        "C74" : ("Nup133", "Nup84", "Nup82", "Nup116", "Nup159", "Nsp1", "Nic96",
                 "Nup188", "Nup192", "Nup170", "Nup53", "Nup157", 2),
        "C75" : ("Nup120", "Nup145C", "Nup84", "Nup82", "Nup159", "Nsp1", "Nup188",
                 "Nup1", "Nup192", "Pom152", "Nup170", "Nup157", 3),
        "C76" : ("Pom34", "Pom152", "Nup192", "Nup170", "Nup53", "Nup59", "Nup157",
                 "Nic96", "Nup188", "Nsp1", "Nup82", "Gle2", "Nup159", 2),
        "C77" : ("Nup100", "Nup188", "Nic96", "Nup192", "Nup170", "Nup53", "Nup59",
                 "Nsp1", "Nup57", "Nup82", "Nup84", "Nup116", "Nup159", "Gle1",
                 "Nup42", 3),
        "C78" : ("Nup133", "Nup145C", "Nup120", "Nup84", "Nup82", "Nup116",
                 "Nup159", "Nsp1", "Nic96", "Nup188", "Nup100", "Nup192", "Pom152",
                 "Nup170", "Nup59", 2),
        "C79" : ("Nup120", "Seh1", "Nup85", "Nup84", "Nup145C", "Nup133", "Sec13",
                 "Nup82", "Nup159", "Nsp1", "Nic96", "Nup188", "Nup192",
                 "Nup170", "Pom152", 1),
        "C80" : ("Nup120", "Seh1", "Nup85", "Nup84", "Sec13", "Nup145C", "Nup133",
                 "Nic96", "Nsp1", "Nup57", "Nup188", "Nup192", "Pom152",
                 "Nup170", "Nup157", "Nup53", 1),
        "C81" : ("Nup42", "Gle1", "Nup82", "Nup116", "Gle2", "Nup159", "Nup84",
                 "Nup133", "Nup120", "Nup145N", "Nsp1", "Nup57", "Nic96", "Nup188",
                 "Nup192", "Nup170", "Nup157", "Pom152", "Pom34", 1),
        "C82" : ("Nup133", "Nup84", "Nup82", "Gle1", "Nup42", "Nup159", "Nup116",
                 "Gle2", "Nsp1", "Nup57", "Nic96", "Nup188", "Nup100", "Nup192",
                 "Pom152", "Pom34", "Nup170", "Nup59", "Nup157", "Nup53", 1) }
    cr = CompositeRepresentation(simo)
    # Score distance between protein pairs considering that one protein might
    # not be in the primary spoke (i.e. rotated +/- 45 degrees)
    axis = IMP.algebra.Vector3D(0,0,1)
    rot = IMP.algebra.get_rotation_about_axis(axis, math.pi / 2.)
    transforms = [IMP.algebra.Transformation3D(r, IMP.algebra.Vector3D(0,0,0))
                  for r in rot, rot.get_inverse()]
    # Penalize configurations where spheres are not in contact
    penalty = 1.
    uf = IMP.core.HarmonicUpperBound(0., penalty)
    ps = IMP.npc.MinimumSphereDistancePairScore(uf, transforms)

    for c, res in COMPOSITE.items():
        rsr = IMP.npc.CompositeRestraint(simo.m, ps)
        rsr.set_name("CompositeRestraint " + c)
        # Speed up evaluation by immediately discarding subtrees with
        # high scores
        rsr.set_maximum_score(penalty * 5.)
        for protein in res[:-1]:
            rsr.add_type(cr.get_all_copies(protein))
        # Add restraint to the PMI scoring function
        IMP.pmi.tools.add_restraint_to_model(simo.m, rsr)


#####################################################
# Restraints setup - Immuno-EM
# Supplementary Table 7. Upper and lower bounds on R-radial restraints of C-terminal bead of nups
# NupType : (min R value, max R value) (in Angstrom)
# Supplementary Table 7. Upper and lower bounds on Z-axial restraints of C-terminal bead of nups
# NupType : (min Z value, max Z value) (in Angstrom)
#####################################################
if (use_Immuno_EM):
    nup_list = [entry[0] for entry in domains if not '@' in entry[0]]
    nup_list_unique = sorted(list(set(nup_list)))   # Make a unique list

    RADIAL = {
        "Gle1" : [200, 320],
        "Gle2.1" : [150, 330],
        "Gle2.2" : [150, 330],
        "Ndc1" : [280, 400],
        "Nic96.1" : [255, 525],
        "Nic96.2" : [255, 525],
        "Nsp1.1" : [155, 425],
        "Nsp1.2" : [155, 425],
        "Nsp1.3" : [155, 425],
        "Nsp1.4" : [155, 425],
        "Nup1" : [200, 360],
        "Nup100.1" : [230, 330],
        "Nup100.2" : [230, 330],
        "Nup116.1" : [200, 350],    #"Nup116.1" : [250, 350],
        "Nup116.2" : [200, 350],    #"Nup116.2" : [250, 350],
        "Nup120" : [250, 450],      #"Nup120" : [250, 370],
        "Nup133" : [300, 520],      #"Nup133" : [300, 420],
        "Nup145c" : [270, 520],     #"Nup145c" : [270, 470],
        "Nup145.1" : [125, 395],
        "Nup145.2" : [125, 395],
        "Nup157" : [190, 350],
        "Nup159.1" : [200, 430],    #"Nup159.1" : [250, 430],
        "Nup159.2" : [200, 430],    #"Nup159.2" : [250, 430],
        "Nup170" : [170, 330],
        "Nup188" : [200, 320],
        "Nup192" : [200, 320],
        "Nup42" : [220, 400],
        "Nup49.1" : [200, 300],
        "Nup49.2" : [200, 300],
        "Nup53" : [280, 380],
        "Nup57.1" : [80, 300],
        "Nup57.2" : [80, 300],
        "Nup59" : [250, 370],
        "Nup60" : [240, 400],
        "Nup82.1" : [175, 505],
        "Nup82.2" : [175, 505],
        "Nup84" : [290, 520],       #"Nup84" : [290, 450],
        "Nup85" : [300, 520],       #"Nup85" : [300, 420],
        "Pom34" : [280, 380],
        "Seh1" : [250, 370],
        "Pom152" : [470, 630]       #"Pom152" : [370, 630]
    }
    print "\nXYRadialPositionRestraint !!"
    radial_weight = 10.0
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
        "Gle1" : [110, 190],        #"Gle1" : [120, 180],
        "Gle2.1" : [20, 120],
        "Gle2.2" : [20, 120],
        "Ndc1" : [0, 90],
        "Nic96.1" : [25, 175],
        "Nic96.2" : [25, 175],
        "Nsp1.1" : [0, 150],        #"Nsp1.1" : [0, 120],
        "Nsp1.2" : [0, 150],        #"Nsp1.2" : [0, 120],
        "Nsp1.3" : [0, 150],        #"Nsp1.2" : [0, 120],
        "Nsp1.4" : [0, 150],        #"Nsp1.2" : [0, 120],
        "Nup1" : [-220, -140],
        "Nup100.1" : [40, 120],
        "Nup100.2" : [40, 120],
        "Nup116.1" : [70, 170],     #"Nup116.1" : [70, 150],
        "Nup116.2" : [70, 170],     #"Nup116.2" : [70, 150],
        "Nup120" : [70, 150],
        "Nup133" : [100, 200],
        "Nup145c" : [70, 150],
        "Nup145.1" : [-170, -50],
        "Nup145.2" : [-170, -50],
        "Nup157" : [0, 95],
        "Nup159.1" : [100, 240],    #"Nup159.1" : [120, 240],
        "Nup159.2" : [100, 240],    #"Nup159.2" : [120, 240],
        "Nup170" : [0, 100],
        "Nup188" : [40, 100],
        "Nup192" : [20, 100],
        "Nup42" : [70, 150],
        "Nup49.1" : [40, 100],
        "Nup49.2" : [40, 100],
        "Nup53" : [20, 100],
        "Nup57.1" : [0, 80],
        "Nup57.2" : [0, 80],
        "Nup59" : [40, 120],
        "Nup60" : [-200, -100],
        "Nup82.1" : [125, 295],     #"Nup82.1" : [145, 295],
        "Nup82.2" : [125, 295],     #"Nup82.2" : [145, 295],
        "Nup84" : [150, 190],       #"Nup84" : [150, 170],
        "Nup85" : [140, 220],       #"Nup85" : [140, 200],
        "Pom34" : [0, 65],
        "Seh1" : [50, 170],
        "Pom152" : [0, 35]          #"Pom152" : [0, 95]
    }
    print "\nZAxialPositionRestraint !!"
    zaxial_weight = 10.0
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
# Restraints setup - Membrane Localization + ALPS Motif
#####################################################
tor_th = 45.0
tor_R = 390.0 + 150.0
tor_r = 150.0 - tor_th/2.0
msl_sigma = 1.0
msl_weight = 1.0

if (is_membrane):
    print "\nMembraneSurfaceLocationRestraint !!"
    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (101,200,'Pom152'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma)
    msl.set_label('Pom152_101_200')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, 'Ndc1', tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma)
    msl.set_label('Ndc1')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, 'Pom34', tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma)
    msl.set_label('Pom34')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (475,475,'Nup53'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma)
    msl.set_label('Nup53')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (528,528,'Nup59'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma)
    msl.set_label('Nup59')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())
    """
    perinuclear_restraints = []
    if PERINUCLEAR:
    # Perinuclear volume location
      for i, beads in enumerate(T.get_beads("Pom152", 7)):
        peri = IMP.npc.PerinuclearVolumeLocationRestraint(m, tor_R, tor_r, tor_th, False, 0.2)
        for b in beads:
          peri.add_particle(b)
        peri.set_name('Perinuclear_Pom152_%d' % i)
        all_restraints.add_restraint(peri)
        perinuclear_restraints.append(peri)


    pore_side_restraints = []
    if PORE_SIDE:
    # Pore-side volume location
      for i, beads in enumerate(T.get_beads("Ndc1", 8)):
        pore = IMP.npc.PoreSideVolumeLocationRestraint(m, tor_R, tor_r, tor_th, False, 0.2)
        for b in beads:
          pore.add_particle(b)
        pore.set_name('Pore_Side_Ndc1_%d' % i)
        all_restraints.add_restraint(pore)
        pore_side_restraints.append(pore)

      for i, beads in enumerate(T.get_beads("Pom152", 8)):
        pore = IMP.npc.PoreSideVolumeLocationRestraint(m, tor_R, tor_r, tor_th, False, 0.2)
        for b in beads:
          pore.add_particle(b)
        pore.set_name('Pore_Side_Pom152_%d' % i)
        all_restraints.add_restraint(pore)
        pore_side_restraints.append(pore)

      for i, beads in enumerate(T.get_beads("Pom34", 8)):
        pore = IMP.npc.PoreSideVolumeLocationRestraint(m, tor_R, tor_r, tor_th, False, 0.2)
        for b in beads:
          pore.add_particle(b)
        pore.set_name('Pore_Side_Pom34_%d' % i)
        all_restraints.add_restraint(pore)
        pore_side_restraints.append(pore)
    """

if (is_n84):
    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (252,270,'Nup133'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma)
    msl.set_label('Nup133')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (135,152,'Nup120'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma)
    msl.set_label('Nup120_1')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (197,216,'Nup120'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma)
    msl.set_label('Nup120_2')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

if (is_inner_ring):
    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (310,338,'Nup157'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma)
    msl.set_label('Nup157')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (320,352,'Nup170'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma)
    msl.set_label('Nup170')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())


#####################################################
# Restraints setup
# Distance_to_point restraints for orientation of the Nup84 complex
#####################################################
if (is_n84 and use_Distance_to_Point):
    dpr_weight = 100.0
    dpr_radius = 100.0

    dpr = IMP.pmi.restraints.basic.DistanceToPointRestraint(simo, tuple_selection=(230,230,"Nup133"), anchor_point=IMP.algebra.Vector3D(433.8, 154.1, 150.0), radius=dpr_radius, kappa=10.0)
    dpr.set_label("Nup133")
    dpr.set_weight(dpr_weight)
    dpr.add_to_model()
    outputobjects.append(dpr)
    print(dpr.get_output())

    dpr = IMP.pmi.restraints.basic.DistanceToPointRestraint(simo, tuple_selection=(324,324,"Nup85"), anchor_point=IMP.algebra.Vector3D(320.2, -201.6, 170.0), radius=dpr_radius, kappa=10.0)
    dpr.set_label("Nup85")
    dpr.set_weight(dpr_weight)
    dpr.add_to_model()
    outputobjects.append(dpr)
    print(dpr.get_output())

    dpr = IMP.pmi.restraints.basic.DistanceToPointRestraint(simo, tuple_selection=(465,465,"Nup120"), anchor_point=IMP.algebra.Vector3D(538.1, -212.6, 110.0), radius=dpr_radius, kappa=10.0)
    dpr.set_label("Nup120")
    dpr.set_weight(dpr_weight)
    dpr.add_to_model()
    outputobjects.append(dpr)
    print(dpr.get_output())
    print "DistanceToPointRestraint !!\n"


#####################################################
# Restraints setup
# Distance restraints
#####################################################
dist_min = 3.0
dr_weight = 10.0

# Nup145n - Nup145c
if (is_n84 and is_nucleoplasm):
    dist_max = 15.0
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(605,605,"Nup145.1"), (1,1,"Nup145c@11"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Nup145n-Nup145c")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

# The Pom152 ring
if (is_membrane):
    dist_max = 30.0
    if (use_neighboring_spokes):
        # TODO: Pom152 orientation?  (clockwise or counter-clockwise?)
        #dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(379,379,"Pom152"), (379,379,"Pom152@13"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(379,379,"Pom152"), (379,379,"Pom152@12"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr.add_to_model()
        dr.set_label("Pom152-Pom152@12")
        dr.set_weight(dr_weight)
        outputobjects.append(dr)
        print(dr.get_output())

    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(379,379,"Pom152"), (1337,1337,"Pom152@11"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Pom152-Pom152@11")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    xyr = IMP.npc.npc_restraints.XYRadialPositionRestraint(simo, (859,859,"Pom152"), lower_bound=530, upper_bound=580, consider_radius=False, sigma=1.0, term='M')
    xyr.set_label('Lower_%d_Upper_%d_%s' % (530, 580, "Pom152_859"))
    xyr.set_weight(radial_weight)
    xyr.add_to_model()
    outputobjects.append(xyr)
    print (xyr.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (859,859,"Pom152"), lower_bound=0, upper_bound=25, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (0, 25, "Pom152_859"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())


#####################################################
# Restraints setup
# Distance restraints for homo-dimers
#####################################################
if (False):
    """
    dist_min = 3.0
    dist_max = 30.0
    dr_weight = 100.0

    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(1417,1417,"Nup159.1"), (1417,1417,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Nup159_1417-1417")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())
    """


#####################################################
# Restraints setup
# Distance restraints for XL cliques
#####################################################
if (False):
    """
    protein1 = columnmap["Protein1"]
    protein2 = columnmap["Protein2"]
    residue1 = columnmap["Residue1"]
    residue2 = columnmap["Residue2"]
    idscore = columnmap["IDScore"]
    xluniqueid = columnmap["XLUniqueID"]

    db = IMP.pmi.tools.get_db_from_csv('../data/XL_cliques_2copies.csv')

    dist_min = 3.0
    dist_max = 35.0
    dr_weight = 10.0

    for nxl, entry in enumerate(db):
        #print nxl, entry

        mol1 = entry[protein1]
        res1 = int(entry[residue1])
        mol2 = entry[protein2]
        res2 = int(entry[residue2])
        id_score = float(entry[idscore])
        xlunique_id = int(entry[xluniqueid])

        dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(res1,res1,mol1), (res2,res2,mol2), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr.add_to_model()
        temp_label = mol1 + "_" + str(res1) + "-" + mol2 + "_" + str(res2)
        dr.set_label(temp_label)
        dr.set_weight(dr_weight)
        outputobjects.append(dr)
        print(dr.get_output())

    print "\nDistance Restraints applied for XL cliques !!"
    print "weight = ", dr_weight, "dist_min = ", dist_min, "dist_max = ", dist_max, "\n"
    """


#####################################################
# 1st Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "\nEVAL 1 : ", sf.evaluate(False), " (initial) - ", rank

simo.optimize_floppy_bodies(300)
print "\nEVAL 2 : ", sf.evaluate(False), " (after calling optimize_floppy_bodies(300)) - ", rank

XL_restraints = None
mc1 = IMP.pmi.macros.ReplicaExchange0(m,
                                    simo,
                                    monte_carlo_sample_objects = sampleobjects,
                                    output_objects = outputobjects,
                                    crosslink_restraints = XL_restraints,
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 2.5,
                                    number_of_best_scoring_models = 0,
                                    monte_carlo_steps = 10,
                                    number_of_frames = initial_nframes,
                                    write_initial_rmf = True,
                                    initial_rmf_name_suffix = "initial",
                                    stat_file_name_suffix = "stat",
                                    best_pdb_name_suffix = "model",
                                    do_clean_first = True,
                                    do_create_directories = True,
                                    global_output_directory = "pre-XL_output",
                                    rmf_dir = "rmfs/",
                                    best_pdb_dir = "pdbs/",
                                    replica_stat_file_suffix = "stat_replica")
mc1.execute_macro()
rex1 = mc1.get_replica_exchange_object()
print "\nEVAL 3 : ", sf.evaluate(False), " (after performing the pre-XL_sampling) - ", rank


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
                                                        '../data_npc/XL_Merged_wholeNPC_sorted_uniq_0927_2013_noFG.csv',
                                                        length = 21.0,
                                                        slope = 0.00,
                                                        columnmapping = columnmap,
                                                        ids_map = ids_map,
                                                        resolution = 1.0,
                                                        inner_slope = 0.03,
                                                        #inner_slope = 0.01,
                                                        filelabel = "wtDSS",
                                                        label = "wtDSS",
                                                        attributes_for_label = ["XLUniqueID"],
                                                        csvfile = True)
    xl1.add_to_model()
    sampleobjects.append(xl1)
    outputobjects.append(xl1)
    xl1.set_psi_is_sampled(False)
    psi2 = xl1.get_psi(1.0)[0]
    psi2.set_scale(0.05)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
    print "\nEVAL 4 : ", sf.evaluate(False), " (after applying the XL restraint) - ", rank
    XL_restraints = [xl1]
else:
    XL_restraints = None


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
                                    replica_exchange_maximum_temperature = 2.5,
                                    number_of_best_scoring_models = 0,
                                    monte_carlo_steps = 10,
                                    number_of_frames = initial_nframes,
                                    write_initial_rmf = True,
                                    initial_rmf_name_suffix = "initial",
                                    stat_file_name_suffix = "stat",
                                    best_pdb_name_suffix = "model",
                                    do_clean_first = True,
                                    do_create_directories = True,
                                    global_output_directory = "pre-EM3D_output",
                                    rmf_dir = "rmfs/",
                                    best_pdb_dir = "pdbs/",
                                    replica_stat_file_suffix = "stat_replica",
                                    replica_exchange_object = rex1)
mc2.execute_macro()
rex2 = mc2.get_replica_exchange_object()
print "\nEVAL 5 : ", sf.evaluate(False), " (after performing the pre-EM3D_sampling) - ", rank


#####################################################
# Restraints setup
# EM 3D restraint using GMM
#####################################################
if (use_EM3D):
    #resdensities = bm1.get_density_hierarchies([t[1] for t in domains])
    resdensities = bm1.get_density_hierarchies(main_spoke_hier_name)
    print ("resdensities=", resdensities)       ####  TODO: make sure resdensities are correct

    mass = sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))
    mass *= 1.2 * 2.0           # 1.2 for adjustment of the GMM (after removing flexible GMMs) and 2.0 for approximation of the NPC spoke mass
    print ("Total mass for the EM restraint = ", mass)
    gem = IMP.pmi.restraints.em.GaussianEMRestraint(resdensities,
                                                    '../data_npc/em_gmm_model/avg_monomer_final_sj2_ring_5p5rot.456.txt',
                                                    target_mass_scale=mass,
                                                    slope=0.0000001,
                                                    target_radii_scale=3.0)
    gem.add_to_model()
    gem.set_weight(10.0)        # play with weight
    #gem.center_model_on_target_density(simo)
    outputobjects.append(gem)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
    print "\nEVAL 6 : ", sf.evaluate(False), " (after applying the EM 3D restraint) - ", rank


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
                                    replica_exchange_maximum_temperature = 2.5,
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
                                    replica_stat_file_suffix = "stat_replica",
                                    replica_exchange_object = rex2)
mc3.execute_macro()
print "\nEVAL 7 : ", sf.evaluate(False), " (final evaluation) - ", rank
exit(0)




"""
if RES_CONTACT:
  print 'ADDING PROTEIN CONTACT RESTRAINTS'
  for res1, res2 in [x.split('-') for x in CONTACT]:
    ri = 1
    for k in [2, 4, 9]:
      beads1 = T.get_beads(res1, k)
      beads2 = T.get_beads(res2, k)
      if len(beads1) != len(beads2):
        print 'WARNING: residue %s has different number of copies than %s - skipping for now' % (res1, res2)
        continue
      for b1, b2 in zip(beads1, beads2):
        n = min(len(b1), len(b2))
        for i in xrange(n):
          c = IMP.npc.ProteinContactRestraint(m, 1.3, 0.01)
          c.set_name('Protein_Contact_%s_%s_%d' % (res1, res2, ri))
          ri += 1
          c.add_particle(b1[i])
          c.add_particle(b2[i])
          all_restraints.add_restraint(c)
          protein_contact_restraints.append(c)
  print_report()
  print '-- after adding contact restraint score is %s' % (all_restraints.evaluate(False))


if PROXIMITY:
  print 'ADDING PROXIMITY RESTRAINTS'
  for cname, compo in COMPOSITE.iteritems():
    ri = 1
    if len(compo) > 5:
      q = compo[-1]
      n = determine_nresidues_complex(compo[:-1])
      if n <= 0:
        continue
      D_max = 1.35*0.495*((q*n)**(1.0/3))
      print '-- adding proximity restraint for complex %s' % cname
      for k in [2, 4, 9]:
        for chain in 'ABCDEFGHIJKLMNOP':
          all_beads = []
          for res in compo[:-1]:
            all_beads.extend(T.get_beads_by_chain(res, chain, k))
          if all_beads:
            pr = IMP.npc.ProteinProximityRestraint(m, D_max, 0.01)
            pr.set_name('Proximity_%s_%d' % (cname, ri))
            ri += 1
            pr.add_particles(all_beads)
            all_restraints.add_restraint(pr)
            proximity_restraints.append(pr)

  print_report()
  print '-- after adding proximity restraint score is %s' % (all_restraints.evaluate(False))

if DIAMETER:
  print 'ADDING DIAMETER RESTRAINT'
  for chain in 'ABCDEFGHIJKLMNOP':
    all_beads = []
    for res in COMPOSITE['C51'][:-1]:
      all_beads.extend(T.get_beads_by_chain(res, chain, 1))
    if all_beads:
      h = IMP.core.HarmonicUpperBound(0, 1)
      r = IMP.core.DiameterRestraint(h, IMP.container.ListSingletonContainer(all_beads), 19.2)
      r.set_name('Diameter_%s' % chain)
      all_restraints.add_restraint(r)
      diameter_restraints.append(r)

  print_report()
  print '-- after adding diameter restraint score is %s' % (all_restraints.evaluate(False))
"""
