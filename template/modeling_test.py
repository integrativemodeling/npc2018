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
#import IMP.pmi.restraints.em
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
simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)


#####################################################
# setting up parameters
#####################################################
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
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

n82 = "../data_nup82/"
f_n82 = "../data_nup82/protein_fasta."

data_nic96 = "../data_nic96/"
fasta_nic96 = "../data_nic96/protein_fasta."

data_npc = "../data_npc/"
fasta_npc = "../data_npc/protein_fasta."

n84_fastafile  ='../data_nup84_2016/protein_fasta.Nup84.txt'
n85_fastafile  ='../data_nup84_2016/protein_fasta.Nup85.txt'
n120_fastafile ='../data_nup84_2016/protein_fasta.Nup120.txt'
n133_fastafile ='../data_nup84_2016/protein_fasta.Nup133.txt'
n145c_fastafile='../data_nup84_2016/protein_fasta.Nup145c.txt'
seh1_fastafile ='../data_nup84_2016/protein_fasta.Seh1.txt'
sec13_fastafile='../data_nup84_2016/protein_fasta.Sec13.txt'

# After removal of disordered regions in PDB files.
n84_pdbfile    ='../data_nup84_2016/4ycz_3ewe/ScNup84_7-488_506-726_4xmmF.pdb'
n85n_pdbfile   ='../data_nup84_2016/4ycz_3ewe/ScNup85_1-492_4yczB_3ewe.pdb'
n85c_pdbfile   ='../data_nup84_2016/4ycz_3ewe/ScNup85_496-744_4yczB.pdb'
n120_pdbfile   ='../data_nup84_2016/4ycz_3ewe/ScNup120_4yczC.pdb'
n133n_pdbfile  ='../data_nup84_2016/4ycz_3ewe/ScNup133N_56-480.pdb'
n133c_pdbfile  ='../data_nup84_2016/4ycz_3ewe/ScNup133C_490_1157.pdb'
n145c_pdbfile  ='../data_nup84_2016/4ycz_3ewe/ScNup145C_4xmmB.pdb'
n145c_pdbfile2 ='../data_nup84_2016/4ycz_3ewe/ScNup145C_92-99_4yczB.pdb'
seh1_pdbfile   ='../data_nup84_2016/4ycz_3ewe/ScSeh1_4yczB_3ewe.pdb'
sec13_pdbfile  ='../data_nup84_2016/4ycz_3ewe/ScSec13_2-296_4xmmA.pdb'

n157N_pdbfile  = data_npc + "Nup157_4mhc_A_88_892.pdb"
n157C_pdbfile  = data_npc + "Nup157_3i5p_A_900_1390.pdb"
n170N_pdbfile  = data_npc + "Nup170_4mhc_A_98_992.pdb"
n170C_pdbfile  = data_npc + "Nup170_3i5p_A_1000_1502.pdb"
n188_pdbfile   = data_npc + "Nup188_12_1652.pdb"
n192_pdbfile   = data_npc + "Nup192_1_1683_v2_new.pdb"

n53_pdbfile    = data_npc + "Nup53_dimer.pdb"
n59_pdbfile    = data_npc + "Nup59_dimer.pdb"
n100_pdbfile   = data_npc + "Nup100_3nf5_AB_803_958.pdb"
n145N_pdbfile  = data_npc + "Nup145N_3kep_AB_459_605.pdb"

Dbp5_pdbfile   = data_npc + "Dbp5_3rrm_A_91_482.pdb"
Gle1N_pdbfile  = data_npc + "Gle1_4wij_A_121_239.pdb"
Gle1C_pdbfile  = data_npc + "Gle1_3rrm_B_244_538.pdb"
Gle2_pdbfile   = data_npc + "Gle2_3mmy_A_4_362.pdb"

# for Debugging
is_n84 = True
is_n82 = True
is_clone = True

use_XL = True
use_Immuno_EM = True
use_Distance_to_Point = True

# Rigid body assignment
n84_rb = 84
n82_rb = 84
domains = []

#####################################################
# REPRESENTATION
#####################################################
# compname  hier_name      color   fastafile                fastaid    pdbname      chain  res_range   read_em_files bead_size  rigid_body super_rigid_body em_num_components em_txt_file_name em_mrc_file_name chain_of_super_rb

##########################
# Nup84 complex
##########################
if (is_n84):
    domains.append(('Nup84',  "Nup84",    0.0, n84_fastafile,   "Nup84",   n84_pdbfile,  "A",  (  1, 726,0),   None,     beadsize,   n84_rb,   None,   7,  None,   None, None))
    domains.append(('Nup85',  "Nup85_1",  0.2, n85_fastafile,   "Nup85",   n85n_pdbfile, "D",  (  1, 492,0),   None,     beadsize,   n84_rb,   None,   8,  None,   None, None))
    domains.append(('Nup85',  "Nup85_2",  0.25,n85_fastafile,   "Nup85",   n85c_pdbfile, "D",  (493, 744,0),   None,     beadsize,   n84_rb,   None,   3,  None,   None, None))
    domains.append(('Nup120', "Nup120_1", 0.35,n120_fastafile,  "Nup120",  n120_pdbfile, "E",  (  1, 711,0),   None,     beadsize,   n84_rb,   None,   7,  None,   None, None))
    domains.append(('Nup120', "Nup120_2", 0.4, n120_fastafile,  "Nup120",  n120_pdbfile, "E",  (712,1037,0),   None,     beadsize,   n84_rb,   None,   3,  None,   None, None))
    domains.append(('Nup133', "Nup133_1", 0.5, n133_fastafile,  "Nup133",  n133n_pdbfile,"D",  (  1, 480,0),   None,     beadsize,   n84_rb,   None,   5,  None,   None, None))
    domains.append(('Nup133', "Nup133_2", 0.55,n133_fastafile,  "Nup133",  n133c_pdbfile,"D",  (481,1157,0),   None,     beadsize,   n84_rb,   None,   7,  None,   None, None))
    domains.append(('Nup145c',"Nup145c_1",0.65,n145c_fastafile, "Nup145c", n145c_pdbfile2,"B", (  1, 125,0),   None,     beadsize,   n84_rb,   None,   7,  None,   None, None))
    domains.append(('Nup145c',"Nup145c_2",0.7, n145c_fastafile, "Nup145c", n145c_pdbfile,"B",  (126, 712,0),   None,     beadsize,   n84_rb,   None,   7,  None,   None, None))
    domains.append(('Seh1',   "Seh1",     0.8, seh1_fastafile,  "Seh1",    seh1_pdbfile, "C",  (  1, 349,0),   None,     beadsize,   n84_rb,   None,   4,  None,   None, None))
    domains.append(('Sec13',  "Sec13",    0.95,sec13_fastafile, "Sec13",   sec13_pdbfile,"G",  (  1, 297,0),   None,     beadsize,   n84_rb,   None,   3,  None,   None, None))

    if (is_clone):
        for i in range(2,4)+range(11,14):
            domains.append(('Nup84@%d'%i,  "Nup84@%d"%i,    0.0, n84_fastafile,   "Nup84",   n84_pdbfile,   "A",  (  1, 726,0),   None,  beadsize,   None,   None,   7,  None,   None, None))
            domains.append(('Nup85@%d'%i,  "Nup85_1@%d"%i,  0.2, n85_fastafile,   "Nup85",   n85n_pdbfile,  "D",  (  1, 492,0),   None,  beadsize,   None,   None,   8,  None,   None, None))
            domains.append(('Nup85@%d'%i,  "Nup85_2@%d"%i,  0.25,n85_fastafile,   "Nup85",   n85c_pdbfile,  "D",  (493, 744,0),   None,  beadsize,   None,   None,   3,  None,   None, None))
            domains.append(('Nup120@%d'%i, "Nup120_1@%d"%i, 0.35,n120_fastafile,  "Nup120",  n120_pdbfile,  "E",  (  1, 711,0),   None,  beadsize,   None,   None,   7,  None,   None, None))
            domains.append(('Nup120@%d'%i, "Nup120_2@%d"%i, 0.4, n120_fastafile,  "Nup120",  n120_pdbfile,  "E",  (712,1037,0),   None,  beadsize,   None,   None,   3,  None,   None, None))
            domains.append(('Nup133@%d'%i, "Nup133_1@%d"%i, 0.5, n133_fastafile,  "Nup133",  n133n_pdbfile, "D",  (  1, 480,0),   None,  beadsize,   None,   None,   5,  None,   None, None))
            domains.append(('Nup133@%d'%i, "Nup133_2@%d"%i, 0.55,n133_fastafile,  "Nup133",  n133c_pdbfile, "D",  (481,1157,0),   None,  beadsize,   None,   None,   7,  None,   None, None))
            domains.append(('Nup145c@%d'%i,"Nup145c_1@%d"%i,0.65,n145c_fastafile, "Nup145c", n145c_pdbfile2,"B",  (  1, 125,0),   None,  beadsize,   None,   None,   7,  None,   None, None))
            domains.append(('Nup145c@%d'%i,"Nup145c_2@%d"%i,0.7, n145c_fastafile, "Nup145c", n145c_pdbfile, "B",  (126, 712,0),   None,  beadsize,   None,   None,   7,  None,   None, None))
            domains.append(('Seh1@%d'%i,   "Seh1@%d"%i,     0.8, seh1_fastafile,  "Seh1",    seh1_pdbfile,  "C",  (  1, 349,0),   None,  beadsize,   None,   None,   4,  None,   None, None))
            domains.append(('Sec13@%d'%i,  "Sec13@%d"%i,    0.95,sec13_fastafile, "Sec13",   sec13_pdbfile, "G",  (  1, 297,0),   None,  beadsize,   None,   None,   3,  None,   None, None))

##########################
# Nup82 complex
##########################
if (is_n82):
    domains.append(("Dyn2.1",  "Dyn2.1",      0.48,  f_n82+"Dyn2.txt",    "Dyn2",    n82+"4DS1.pdb",     "A",  (  1, 92,0),    None, beadsize, n82_rb, None, 4,  None,   None, None))
    domains.append(("Dyn2.2",  "Dyn2.2",      0.65,  f_n82+"Dyn2.txt",    "Dyn2",    n82+"4DS1.pdb",     "C",  (  1, 92,0),    None, beadsize, n82_rb, None, 4,  None,   None, None))
    domains.append(("Nup82.1", "Nup82.1_1",   0.0,   f_n82+"Nup82.txt",   "Nup82",   n82+"3PBP1.pdb",    "G",  (  1,452,0),    None, beadsize, n82_rb, None, 18, None,   None, None))
    domains.append(("Nup82.1", "Nup82.1_11",  0.0,   f_n82+"Nup82.txt",   "Nup82",   "BEADS",            " ",  (453,521,0),    None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup82.1", "Nup82.1_2",   0.0,   f_n82+"Nup82.txt",   "Nup82",   n82+"cc_tr1_1.pdb", "E",  (522,612,0),    None, beadsize, n82_rb, None, 4,  None,   None, None))
    domains.append(("Nup82.1", "Nup82.1_21",  0.0,   f_n82+"Nup82.txt",   "Nup82",   "BEADS",            " ",  (613,624,0),    None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup82.1", "Nup82.1_3",   0.0,   f_n82+"Nup82.txt",   "Nup82",   n82+"cc_tr2_1.pdb", "E",  (625,669,0),    None, beadsize, n82_rb, None, 2,  None,   None, None))
    domains.append(("Nup82.1", "Nup82.1_31",  0.0,   f_n82+"Nup82.txt",   "Nup82",   "BEADS",            " ",  (670,677,0),    None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup82.1", "Nup82.1_4",   0.0,   f_n82+"Nup82.txt",   "Nup82",   n82+"cc_tr3_1.pdb", "E",  (678,713,0),    None, beadsize, n82_rb, None, 2,  None,   None, None))
    domains.append(("Nup82.2", "Nup82.2_1",   0.15,  f_n82+"Nup82.txt",   "Nup82",   n82+"3PBP2.pdb",    "G",  (  1,452,0),    None, beadsize, n82_rb, None, 18, None,   None, None))
    domains.append(("Nup82.2", "Nup82.2_11",  0.15,  f_n82+"Nup82.txt",   "Nup82",   "BEADS",            " ",  (453,521,0),    None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup82.2", "Nup82.2_2",   0.15,  f_n82+"Nup82.txt",   "Nup82",   n82+"cc_tr1_2.pdb", "E",  (522,612,0),    None, beadsize, n82_rb, None, 4,  None,   None, None))
    domains.append(("Nup82.2", "Nup82.2_21",  0.15,  f_n82+"Nup82.txt",   "Nup82",   "BEADS",            " ",  (613,624,0),    None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup82.2", "Nup82.2_3",   0.15,  f_n82+"Nup82.txt",   "Nup82",   n82+"cc_tr2_2.pdb", "E",  (625,669,0),    None, beadsize, n82_rb, None, 2,  None,   None, None))
    domains.append(("Nup82.2", "Nup82.2_31",  0.15,  f_n82+"Nup82.txt",   "Nup82",   "BEADS",            " ",  (670,677,0),    None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup82.2", "Nup82.2_4",   0.15,  f_n82+"Nup82.txt",   "Nup82",   n82+"cc_tr3_2.pdb", "E",  (678,713,0),    None, beadsize, n82_rb, None, 2,  None,   None, None))
    domains.append(("Nup159.1","Nup159.1_1",  1.0,   f_n82+"Nup159.txt",  "Nup159",  n82+"1XIP1.pdb",    "A",  (   1, 381,0),  None, beadsize,    821, None, 0,  None,   None, None))
    domains.append(("Nup159.1","Nup159.1_11", 1.0,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  ( 382,1116,0),  None, beadsize100, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup159.1","Nup159.1_2",  1.0,   f_n82+"Nup159.txt",  "Nup159",  n82+"4DS1.pdb",     "B",  (1117,1126,0),  None, beadsize, n82_rb, None, 1,  None,   None, None))
    domains.append(("Nup159.1","Nup159.1_21", 1.0,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1127,1210,0),  None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup159.1","Nup159.1_3",  1.0,   f_n82+"Nup159.txt",  "Nup159",  n82+"cc_tr1_1.pdb", "D",  (1211,1321,0),  None, beadsize, n82_rb, None, 4,  None,   None, None))
    domains.append(("Nup159.1","Nup159.1_31", 1.0,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1322,1331,0),  None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup159.1","Nup159.1_4",  1.0,   f_n82+"Nup159.txt",  "Nup159",  n82+"cc_tr2_1.pdb", "D",  (1332,1372,0),  None, beadsize, n82_rb, None, 2,  None,   None, None))
    domains.append(("Nup159.1","Nup159.1_41", 1.0,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1373,1381,0),  None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup159.1","Nup159.1_5",  1.0,   f_n82+"Nup159.txt",  "Nup159",  n82+"cc_tr3_1.pdb", "D",  (1382,1412,0),  None, beadsize, n82_rb, None, 2,  None,   None, None))
    domains.append(("Nup159.1","Nup159.1_51", 1.0,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1413,1428,0),  None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup159.1","Nup159.1_6",  1.0,   f_n82+"Nup159.txt",  "Nup159",  n82+"3PBP1.pdb",    "I",  (1429,1460,0),  None, beadsize, n82_rb, None, 2,  None,   None, None))
    domains.append(("Nup159.2","Nup159.2_1",  0.9,   f_n82+"Nup159.txt",  "Nup159",  n82+"1XIP2.pdb",    "A",  (   1, 381,0),  None, beadsize,    822, None, 0,  None,   None, None))
    domains.append(("Nup159.2","Nup159.2_11", 0.9,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  ( 382,1116,0),  None, beadsize100, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup159.2","Nup159.2_2",  0.9,   f_n82+"Nup159.txt",  "Nup159",  n82+"4DS1.pdb",     "D",  (1117,1126,0),  None, beadsize, n82_rb, None, 1,  None,   None, None))
    domains.append(("Nup159.2","Nup159.2_21", 0.9,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1127,1210,0),  None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup159.2","Nup159.2_3",  0.9,   f_n82+"Nup159.txt",  "Nup159",  n82+"cc_tr1_2.pdb", "D",  (1211,1321,0),  None, beadsize, n82_rb, None, 4,  None,   None, None))
    domains.append(("Nup159.2","Nup159.2_31", 0.9,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1322,1331,0),  None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup159.2","Nup159.2_4",  0.9,   f_n82+"Nup159.txt",  "Nup159",  n82+"cc_tr2_2.pdb", "D",  (1332,1372,0),  None, beadsize, n82_rb, None, 2,  None,   None, None))
    domains.append(("Nup159.2","Nup159.2_41", 0.9,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1373,1381,0),  None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup159.2","Nup159.2_5",  0.9,   f_n82+"Nup159.txt",  "Nup159",  n82+"cc_tr3_2.pdb", "D",  (1382,1412,0),  None, beadsize, n82_rb, None, 2,  None,   None, None))
    domains.append(("Nup159.2","Nup159.2_51", 0.9,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1413,1428,0),  None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup159.2","Nup159.2_6",  0.9,   f_n82+"Nup159.txt",  "Nup159",  n82+"3PBP2.pdb",    "I",  (1429,1460,0),  None, beadsize, n82_rb, None, 2,  None,   None, None))
    domains.append(("Nsp1.1",  "Nsp1.1_10",   0.3,   f_n82+"Nsp1.txt",    "Nsp1",    "BEADS",            " ",  (  1,636,0),    None, beadsize100, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nsp1.1",  "Nsp1.1_1",    0.3,   f_n82+"Nsp1.txt",    "Nsp1",    n82+"cc_tr1_1.pdb", "C",  (637,727,0),    None, beadsize, n82_rb, None, 4,  None,   None, None))
    domains.append(("Nsp1.1",  "Nsp1.1_11",   0.3,   f_n82+"Nsp1.txt",    "Nsp1",    "BEADS",            " ",  (728,741,0),    None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nsp1.1",  "Nsp1.1_2",    0.3,   f_n82+"Nsp1.txt",    "Nsp1",    n82+"cc_tr2_1.pdb", "C",  (742,778,0),    None, beadsize, n82_rb, None, 2,  None,   None, None))
    domains.append(("Nsp1.1",  "Nsp1.1_21",   0.3,   f_n82+"Nsp1.txt",    "Nsp1",    "BEADS",            " ",  (779,787,0),    None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nsp1.1",  "Nsp1.1_3",    0.3,   f_n82+"Nsp1.txt",    "Nsp1",    n82+"cc_tr3_1.pdb", "C",  (788,823,0),    None, beadsize, n82_rb, None, 2,  None,   None, None))
    domains.append(("Nsp1.2",  "Nsp1.2_10",   0.38,  f_n82+"Nsp1.txt",    "Nsp1",    "BEADS",            " ",  (  1,636,0),    None, beadsize100, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nsp1.2",  "Nsp1.2_1",    0.38,  f_n82+"Nsp1.txt",    "Nsp1",    n82+"cc_tr1_2.pdb", "C",  (637,727,0),    None, beadsize, n82_rb, None, 4,  None,   None, None))
    domains.append(("Nsp1.2",  "Nsp1.2_11",   0.38,  f_n82+"Nsp1.txt",    "Nsp1",    "BEADS",            " ",  (728,741,0),    None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nsp1.2",  "Nsp1.2_2",    0.38,  f_n82+"Nsp1.txt",    "Nsp1",    n82+"cc_tr2_2.pdb", "C",  (742,778,0),    None, beadsize, n82_rb, None, 2,  None,   None, None))
    domains.append(("Nsp1.2",  "Nsp1.2_21",   0.38,  f_n82+"Nsp1.txt",    "Nsp1",    "BEADS",            " ",  (779,787,0),    None, beadsize, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nsp1.2",  "Nsp1.2_3",    0.38,  f_n82+"Nsp1.txt",    "Nsp1",    n82+"cc_tr3_2.pdb", "C",  (788,823,0),    None, beadsize, n82_rb, None, 2,  None,   None, None))
    domains.append(("Nup116.1","Nup116.1_10", 0.75,  f_n82+"Nup116.txt",  "Nup116",  "BEADS",            " ",  (  1, 750,0),   None, beadsize100, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup116.1","Nup116.1_1",  0.75,  f_n82+"Nup116.txt",  "Nup116",  n82+"3PBP1.pdb",    "H",  (751,1113,0),   None, beadsize25,  n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup116.2","Nup116.2_10", 0.8,   f_n82+"Nup116.txt",  "Nup116",  "BEADS",            " ",  (  1, 750,0),   None, beadsize100, n82_rb, None, 0,  None,   None, None))
    domains.append(("Nup116.2","Nup116.2_1",  0.8,   f_n82+"Nup116.txt",  "Nup116",  n82+"3PBP2.pdb",    "H",  (751,1113,0),   None, beadsize25,  n82_rb, None, 0,  None,   None, None))

    if (is_clone):
        for i in range(2,4):
            domains.append(("Dyn2.1@%d"%i,  "Dyn2.1@%d"%i,      0.48,  f_n82+"Dyn2.txt",    "Dyn2",    n82+"4DS1.pdb",     "A",  (  1, 92,0),    None, beadsize, None, None, 4,  None,   None, None))
            domains.append(("Dyn2.2@%d"%i,  "Dyn2.2@%d"%i,      0.65,  f_n82+"Dyn2.txt",    "Dyn2",    n82+"4DS1.pdb",     "C",  (  1, 92,0),    None, beadsize, None, None, 4,  None,   None, None))
            domains.append(("Nup82.1@%d"%i, "Nup82.1_1@%d"%i,   0.0,   f_n82+"Nup82.txt",   "Nup82",   n82+"3PBP1.pdb",    "G",  (  1,452,0),    None, beadsize, None, None, 18, None,   None, None))
            domains.append(("Nup82.1@%d"%i, "Nup82.1_11@%d"%i,  0.0,   f_n82+"Nup82.txt",   "Nup82",   "BEADS",            " ",  (453,521,0),    None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nup82.1@%d"%i, "Nup82.1_2@%d"%i,   0.0,   f_n82+"Nup82.txt",   "Nup82",   n82+"cc_tr1_1.pdb", "E",  (522,612,0),    None, beadsize, None, None, 4,  None,   None, None))
            domains.append(("Nup82.1@%d"%i, "Nup82.1_21@%d"%i,  0.0,   f_n82+"Nup82.txt",   "Nup82",   "BEADS",            " ",  (613,624,0),    None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nup82.1@%d"%i, "Nup82.1_3@%d"%i,   0.0,   f_n82+"Nup82.txt",   "Nup82",   n82+"cc_tr2_1.pdb", "E",  (625,669,0),    None, beadsize, None, None, 2,  None,   None, None))
            domains.append(("Nup82.1@%d"%i, "Nup82.1_31@%d"%i,  0.0,   f_n82+"Nup82.txt",   "Nup82",   "BEADS",            " ",  (670,677,0),    None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nup82.1@%d"%i, "Nup82.1_4@%d"%i,   0.0,   f_n82+"Nup82.txt",   "Nup82",   n82+"cc_tr3_1.pdb", "E",  (678,713,0),    None, beadsize, None, None, 2,  None,   None, None))
            domains.append(("Nup82.2@%d"%i, "Nup82.2_1@%d"%i,   0.15,  f_n82+"Nup82.txt",   "Nup82",   n82+"3PBP2.pdb",    "G",  (  1,452,0),    None, beadsize, None, None, 18, None,   None, None))
            domains.append(("Nup82.2@%d"%i, "Nup82.2_11@%d"%i,  0.15,  f_n82+"Nup82.txt",   "Nup82",   "BEADS",            " ",  (453,521,0),    None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nup82.2@%d"%i, "Nup82.2_2@%d"%i,   0.15,  f_n82+"Nup82.txt",   "Nup82",   n82+"cc_tr1_2.pdb", "E",  (522,612,0),    None, beadsize, None, None, 4,  None,   None, None))
            domains.append(("Nup82.2@%d"%i, "Nup82.2_21@%d"%i,  0.15,  f_n82+"Nup82.txt",   "Nup82",   "BEADS",            " ",  (613,624,0),    None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nup82.2@%d"%i, "Nup82.2_3@%d"%i,   0.15,  f_n82+"Nup82.txt",   "Nup82",   n82+"cc_tr2_2.pdb", "E",  (625,669,0),    None, beadsize, None, None, 2,  None,   None, None))
            domains.append(("Nup82.2@%d"%i, "Nup82.2_31@%d"%i,  0.15,  f_n82+"Nup82.txt",   "Nup82",   "BEADS",            " ",  (670,677,0),    None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nup82.2@%d"%i, "Nup82.2_4@%d"%i,   0.15,  f_n82+"Nup82.txt",   "Nup82",   n82+"cc_tr3_2.pdb", "E",  (678,713,0),    None, beadsize, None, None, 2,  None,   None, None))
            domains.append(("Nup159.1@%d"%i,"Nup159.1_1@%d"%i,  1.0,   f_n82+"Nup159.txt",  "Nup159",  n82+"1XIP1.pdb",    "A",  (   1, 381,0),  None, beadsize,    None, None, 0,  None,   None, None))
            domains.append(("Nup159.1@%d"%i,"Nup159.1_11@%d"%i, 1.0,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  ( 382,1116,0),  None, beadsize100, None, None, 0,  None,   None, None))
            domains.append(("Nup159.1@%d"%i,"Nup159.1_2@%d"%i,  1.0,   f_n82+"Nup159.txt",  "Nup159",  n82+"4DS1.pdb",     "B",  (1117,1126,0),  None, beadsize, None, None, 1,  None,   None, None))
            domains.append(("Nup159.1@%d"%i,"Nup159.1_21@%d"%i, 1.0,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1127,1210,0),  None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nup159.1@%d"%i,"Nup159.1_3@%d"%i,  1.0,   f_n82+"Nup159.txt",  "Nup159",  n82+"cc_tr1_1.pdb", "D",  (1211,1321,0),  None, beadsize, None, None, 4,  None,   None, None))
            domains.append(("Nup159.1@%d"%i,"Nup159.1_31@%d"%i, 1.0,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1322,1331,0),  None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nup159.1@%d"%i,"Nup159.1_4@%d"%i,  1.0,   f_n82+"Nup159.txt",  "Nup159",  n82+"cc_tr2_1.pdb", "D",  (1332,1372,0),  None, beadsize, None, None, 2,  None,   None, None))
            domains.append(("Nup159.1@%d"%i,"Nup159.1_41@%d"%i, 1.0,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1373,1381,0),  None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nup159.1@%d"%i,"Nup159.1_5@%d"%i,  1.0,   f_n82+"Nup159.txt",  "Nup159",  n82+"cc_tr3_1.pdb", "D",  (1382,1412,0),  None, beadsize, None, None, 2,  None,   None, None))
            domains.append(("Nup159.1@%d"%i,"Nup159.1_51@%d"%i, 1.0,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1413,1428,0),  None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nup159.1@%d"%i,"Nup159.1_6@%d"%i,  1.0,   f_n82+"Nup159.txt",  "Nup159",  n82+"3PBP1.pdb",    "I",  (1429,1460,0),  None, beadsize, None, None, 2,  None,   None, None))
            domains.append(("Nup159.2@%d"%i,"Nup159.2_1@%d"%i,  0.9,   f_n82+"Nup159.txt",  "Nup159",  n82+"1XIP2.pdb",    "A",  (   1, 381,0),  None, beadsize,    None, None, 0,  None,   None, None))
            domains.append(("Nup159.2@%d"%i,"Nup159.2_11@%d"%i, 0.9,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  ( 382,1116,0),  None, beadsize100, None, None, 0,  None,   None, None))
            domains.append(("Nup159.2@%d"%i,"Nup159.2_2@%d"%i,  0.9,   f_n82+"Nup159.txt",  "Nup159",  n82+"4DS1.pdb",     "D",  (1117,1126,0),  None, beadsize, None, None, 1,  None,   None, None))
            domains.append(("Nup159.2@%d"%i,"Nup159.2_21@%d"%i, 0.9,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1127,1210,0),  None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nup159.2@%d"%i,"Nup159.2_3@%d"%i,  0.9,   f_n82+"Nup159.txt",  "Nup159",  n82+"cc_tr1_2.pdb", "D",  (1211,1321,0),  None, beadsize, None, None, 4,  None,   None, None))
            domains.append(("Nup159.2@%d"%i,"Nup159.2_31@%d"%i, 0.9,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1322,1331,0),  None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nup159.2@%d"%i,"Nup159.2_4@%d"%i,  0.9,   f_n82+"Nup159.txt",  "Nup159",  n82+"cc_tr2_2.pdb", "D",  (1332,1372,0),  None, beadsize, None, None, 2,  None,   None, None))
            domains.append(("Nup159.2@%d"%i,"Nup159.2_41@%d"%i, 0.9,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1373,1381,0),  None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nup159.2@%d"%i,"Nup159.2_5@%d"%i,  0.9,   f_n82+"Nup159.txt",  "Nup159",  n82+"cc_tr3_2.pdb", "D",  (1382,1412,0),  None, beadsize, None, None, 2,  None,   None, None))
            domains.append(("Nup159.2@%d"%i,"Nup159.2_51@%d"%i, 0.9,   f_n82+"Nup159.txt",  "Nup159",  "BEADS",            " ",  (1413,1428,0),  None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nup159.2@%d"%i,"Nup159.2_6@%d"%i,  0.9,   f_n82+"Nup159.txt",  "Nup159",  n82+"3PBP2.pdb",    "I",  (1429,1460,0),  None, beadsize, None, None, 2,  None,   None, None))
            domains.append(("Nsp1.1@%d"%i,  "Nsp1.1_10@%d"%i,   0.3,   f_n82+"Nsp1.txt",    "Nsp1",    "BEADS",            " ",  (  1,636,0),    None, beadsize100, None, None, 0,  None,   None, None))
            domains.append(("Nsp1.1@%d"%i,  "Nsp1.1_1@%d"%i,    0.3,   f_n82+"Nsp1.txt",    "Nsp1",    n82+"cc_tr1_1.pdb", "C",  (637,727,0),    None, beadsize, None, None, 4,  None,   None, None))
            domains.append(("Nsp1.1@%d"%i,  "Nsp1.1_11@%d"%i,   0.3,   f_n82+"Nsp1.txt",    "Nsp1",    "BEADS",            " ",  (728,741,0),    None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nsp1.1@%d"%i,  "Nsp1.1_2@%d"%i,    0.3,   f_n82+"Nsp1.txt",    "Nsp1",    n82+"cc_tr2_1.pdb", "C",  (742,778,0),    None, beadsize, None, None, 2,  None,   None, None))
            domains.append(("Nsp1.1@%d"%i,  "Nsp1.1_21@%d"%i,   0.3,   f_n82+"Nsp1.txt",    "Nsp1",    "BEADS",            " ",  (779,787,0),    None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nsp1.1@%d"%i,  "Nsp1.1_3@%d"%i,    0.3,   f_n82+"Nsp1.txt",    "Nsp1",    n82+"cc_tr3_1.pdb", "C",  (788,823,0),    None, beadsize, None, None, 2,  None,   None, None))
            domains.append(("Nsp1.2@%d"%i,  "Nsp1.2_10@%d"%i,   0.38,  f_n82+"Nsp1.txt",    "Nsp1",    "BEADS",            " ",  (  1,636,0),    None, beadsize100, None, None, 0,  None,   None, None))
            domains.append(("Nsp1.2@%d"%i,  "Nsp1.2_1@%d"%i,    0.38,  f_n82+"Nsp1.txt",    "Nsp1",    n82+"cc_tr1_2.pdb", "C",  (637,727,0),    None, beadsize, None, None, 4,  None,   None, None))
            domains.append(("Nsp1.2@%d"%i,  "Nsp1.2_11@%d"%i,   0.38,  f_n82+"Nsp1.txt",    "Nsp1",    "BEADS",            " ",  (728,741,0),    None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nsp1.2@%d"%i,  "Nsp1.2_2@%d"%i,    0.38,  f_n82+"Nsp1.txt",    "Nsp1",    n82+"cc_tr2_2.pdb", "C",  (742,778,0),    None, beadsize, None, None, 2,  None,   None, None))
            domains.append(("Nsp1.2@%d"%i,  "Nsp1.2_21@%d"%i,   0.38,  f_n82+"Nsp1.txt",    "Nsp1",    "BEADS",            " ",  (779,787,0),    None, beadsize, None, None, 0,  None,   None, None))
            domains.append(("Nsp1.2@%d"%i,  "Nsp1.2_3@%d"%i,    0.38,  f_n82+"Nsp1.txt",    "Nsp1",    n82+"cc_tr3_2.pdb", "C",  (788,823,0),    None, beadsize, None, None, 2,  None,   None, None))
            domains.append(("Nup116.1@%d"%i,"Nup116.1_10@%d"%i, 0.75,  f_n82+"Nup116.txt",  "Nup116",  "BEADS",            " ",  (  1, 750,0),   None, beadsize100, None, None, 0,  None,   None, None))
            domains.append(("Nup116.1@%d"%i,"Nup116.1_1@%d"%i,  0.75,  f_n82+"Nup116.txt",  "Nup116",  n82+"3PBP1.pdb",    "H",  (751,1113,0),   None, beadsize25,  None, None, 0,  None,   None, None))
            domains.append(("Nup116.2@%d"%i,"Nup116.2_10@%d"%i, 0.8,   f_n82+"Nup116.txt",  "Nup116",  "BEADS",            " ",  (  1, 750,0),   None, beadsize100, None, None, 0,  None,   None, None))
            domains.append(("Nup116.2@%d"%i,"Nup116.2_1@%d"%i,  0.8,   f_n82+"Nup116.txt",  "Nup116",  n82+"3PBP2.pdb",    "H",  (751,1113,0),   None, beadsize25,  None, None, 0,  None,   None, None))

# missing - Mlp1, Mlp2, Dbp5, Gle1, Gle2, Nup1, Nup2, Pom152
#("POM152",  "POM152",      1.0,   fasta_npc+"POM152.txt",    "YMR129W", "BEADS",                      " ",  (  1,1337,0),   None,         beadsize100,    1337,   [59],       0,               None,            None, [59]),

bm1 = IMP.pmi.macros.BuildModel1(simo)
#bm1.set_gmm_models_directory(datadirectory + "em_gmm_model/")

if (inputs.rmf_input is not None) :
    n82=['Dyn2.1', 'Dyn2.2', 'Nup82.1', 'Nup82.2', 'Nup159.1', 'Nup159.2', 'Nsp1.1', 'Nsp1.2', 'Nup116.1', 'Nup116.2']
    n84=['Nup84', 'Nup85', 'Nup120', 'Nup133', 'Nup145c', 'Seh1', 'Sec13']
    nic96=['Nic96', 'Nsp1', 'Nup49', 'Nup57']
    nup170=['Nup157', 'Nup170', 'Nup188', 'Nup192']
    nup42=['Nup42', 'Nup53', 'Nup59']

    #n82=set([s[0] for s in domains])
    for d in list(n82):
        bm1.set_rmf_file(d, inputs.rmf_input, int(inputs.rmf_frame_number))
        #bm1.set_rmf_file(d, "../data_nup82/rmfs/B_8_1-95.rmf3", 0)

    for d in list(n84):
        bm1.set_rmf_file(d, inputs.rmf_input, int(inputs.rmf_frame_number))
        #bm1.set_rmf_file(d, "../data_nup82/rmfs/B_8_1-95.rmf3", 0)

    for d in list(nic96):
        bm1.set_rmf_file(d, "../data_nic96/rmfs/30-162.rmf3", 0)

# remove connectivity for clones
clone_list = [entry[0] for entry in domains if '@' in entry[0]]
clone_list_unique = list(set(clone_list))   # Make a unique list
print ("clone_list_unique = ", clone_list_unique)

bm1.build_model(data_structure = domains, sequence_connectivity_scale=2.0, sequence_connectivity_resolution=1.0,
                skip_connectivity_these_domains=clone_list_unique)
#bm1.scale_bead_radii(40,0.8)
#resdensities = bm1.get_density_hierarchies([t[1] for t in domainxl_cliques_psi = 0.25s])
#print resdensities; exit()

#model_ps = []
#for h in self.densities:
#    model_ps += IMP.atom.get_leaves(h)


#####################################################
# apply the rotational symmetry
#####################################################
if (is_n84 and is_clone):
    for protein in ['Nup84','Nup85','Nup120','Nup133','Nup145c','Seh1','Sec13']:
        simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis="X")
        simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis="Z", nSymmetry=8)
        simo.create_rotational_symmetry(protein+'@11', [protein+'@%d'%i for i in range(12,14)], rotational_axis="Z", nSymmetry=8)
if (is_n82 and is_clone):
    for protein in ['Dyn2.1', 'Dyn2.2', 'Nup82.1', 'Nup82.2', 'Nup159.1', 'Nup159.2', 'Nsp1.1', 'Nsp1.2', 'Nup116.1', 'Nup116.2']:
        simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis="Z", nSymmetry=8)


#####################################################
# rigidify floppy bodies
#####################################################
rigid_tuples = []
for protein in ['Nup84','Nup85','Nup120','Nup133','Nup145c','Seh1','Sec13']:
    rigid_tuples.append(protein)

for protein in ['Dyn2.1', 'Dyn2.2', 'Nup82.1', 'Nup82.2']:
    rigid_tuples.append(protein)

#for protein in [(1,381,'Nup159.1'),(1,381,'Nup159.2'), (1117,1460,'Nup159.1'),(1117,1460,'Nup159.2'), (637,823,'Nsp1.1'),(637,823,'Nsp1.2'), (751,1113,'Nup116.1'),(751,1113,'Nup116.2')]:
for protein in [(1,381,'Nup159.2'), (1117,1460,'Nup159.1'),(1117,1460,'Nup159.2'), (637,823,'Nsp1.1'),(637,823,'Nsp1.2'), (751,1113,'Nup116.1'),(751,1113,'Nup116.2')]:
    rigid_tuples.append(protein)

# Remove flexible movers for all clones
for protein in clone_list_unique:
    rigid_tuples.append(protein)

print ("rigid_tuples = ", rigid_tuples)
for rt in rigid_tuples:
    hs = IMP.pmi.tools.select_by_tuple(simo,rt)
    simo.remove_floppy_bodies(hs)

"""
rigid_tuples = [(20,56,'Nic96'),(107,159,'Nic96'),(205,839,'Nic96'), (637,823,'Nsp1'), (270,472,'Nup49'), (287,541,'Nup57')]
for rt in rigid_tuples:
    hs = IMP.pmi.tools.select_by_tuple(simo,rt)
    simo.remove_floppy_bodies(hs)
"""


#####################################################
# randomize the initial configuration
#####################################################
if (True) :
    #simo.shuffle_configuration(max_translation=1, avoidcollision=False, ignore_initial_coordinates=True)
    simo.shuffle_configuration(bounding_box=((200, -125, 50), (700, 125, 300)), ignore_initial_coordinates=True)


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
# Excluded Volume restraint
#####################################################
included_objects = [];  other_objects = []
for entry in domains:
    if '@11' in entry[0]:
        included_objects.append(entry[0])
    elif '@' in entry[0]:
        other_objects.append(entry[0])
    else:
        included_objects.append(entry[0])
# Make unique lists
included_objects_unique = list(set(included_objects))
other_objects_unique = list(set(other_objects))

included_objects = [];  other_objects = []
for entry in included_objects_unique:
    obj = simo.hier_dict[entry]
    included_objects.append(obj)
    other_objects.append(obj)
for entry in other_objects_unique:
    other_objects.append(simo.hier_dict[entry])
print ("EV included_objects = ", included_objects)
print ("EV other_objects = ", other_objects)

ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo,
                                                             included_objects = included_objects,
                                                             other_objects = other_objects,
                                                             resolution = res_ev)
#ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo, resolution = res_ev)
ev.add_to_model()
outputobjects.append(ev)
print(ev.get_output())
print "ExcludedVolumeSphere !!\n"


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
# Restraints setup - Immuno-EM
# Supplementary Table 7. Upper and lower bounds on R-radial restraints of C-terminal bead of nups
# NupType : (min R value, max R value) (in Angstrom)
# Supplementary Table 7. Upper and lower bounds on Z-axial restraints of C-terminal bead of nups
# NupType : (min Z value, max Z value) (in Angstrom)
#####################################################
if (use_Immuno_EM):
    nup_list = [entry[0] for entry in domains if not '@' in entry[0]]
    nup_list_unique = list(set(nup_list))   # Make a unique list

    RADIAL = {
        "Gle1" : [200, 320],
        "Gle2" : [150, 330],
        "Ndc1" : [280, 1400],
        "Nic96" : [255, 525],
        "Nsp1.1" : [155, 425],
        "Nsp1.2" : [155, 425],
        "Nup1" : [200, 360],
        "Nup100" : [230, 330],
        "Nup116.1" : [200, 350],    #"Nup116.1" : [250, 350],
        "Nup116.2" : [200, 350],    #"Nup116.2" : [250, 350],
        "Nup120" : [250, 450],      #"Nup120" : [250, 380],
        "Nup133" : [300, 520],      #"Nup133" : [300, 420],
        "Nup145C" : [270, 520],     #"Nup145C" : [270, 470],
        "Nup145N" : [125, 395],
        "Nup157" : [190, 350],
        "Nup159.1" : [200, 430],    #"Nup159.1" : [250, 430],
        "Nup159.2" : [200, 430],    #"Nup159.2" : [250, 430],
        "Nup170" : [170, 330],
        "Nup188" : [200, 330],
        "Nup192" : [200, 320],
        "Nup42" : [220, 400],
        "Nup49" : [200, 300],
        "Nup53" : [280, 380],
        "Nup57" : [80, 300],
        "Nup59" : [250, 370],
        "Nup60" : [240, 400],
        "Nup82.1" : [175, 505],
        "Nup82.2" : [175, 505],
        "Nup84" : [290, 520],       #"Nup84" : [290, 450],
        "Nup85" : [300, 520],       #"Nup85" : [300, 420],
        "Pom152" : [370, 1630],
        "Pom34" : [280, 1380],
        "Seh1" : [250, 370],
        "Sec13" : [50, 550]
    }
    print "\nXYRadialPositionLowerRestraint !!"
    print "XYRadialPositionUpperRestraint !!\n"
    radial_weight = 1.0
    for protein, r in RADIAL.iteritems():
        if (protein not in nup_list_unique):
            continue

        xyr = IMP.npc.npc_restraints.XYRadialPositionLowerRestraint(simo, protein, r[0], False, 1.0)
        xyr.set_label('Radial_Position_Lower_%d_%s' % (r[0], protein))
        xyr.set_weight(radial_weight)
        xyr.add_to_model()
        outputobjects.append(xyr)
        print (xyr.get_output())

        xyr = IMP.npc.npc_restraints.XYRadialPositionUpperRestraint(simo, protein, r[1], False, 1.0)
        xyr.set_label('Radial_Position_Upper_%d_%s' % (r[1], protein))
        xyr.set_weight(radial_weight)
        xyr.add_to_model()
        outputobjects.append(xyr)
        print (xyr.get_output())

    ZAXIAL = {
        "Gle1" : [110, 170],
        "Gle2" : [20, 120],
        "Ndc1" : [0, 90],
        "Nic96" : [25, 175],
        "Nsp1.1" : [0, 150],        #"Nsp1.1" : [0, 120],
        "Nsp1.2" : [0, 150],        #"Nsp1.2" : [0, 120],
        "Nup1" : [-220, -140],
        "Nup100" : [40, 120],
        "Nup116.1" : [70, 170],     #"Nup116.1" : [70, 150],
        "Nup116.2" : [70, 170],     #"Nup116.2" : [70, 150],
        "Nup120" : [70, 150],
        "Nup133" : [100, 200],
        "Nup145C" : [70, 150],
        "Nup145N" : [-170, -50],
        "Nup157" : [0, 95],
        "Nup159.1" : [100, 240],    #"Nup159.1" : [120, 240],
        "Nup159.2" : [100, 240],    #"Nup159.2" : [120, 240],
        "Nup170" : [0, 75],
        "Nup188" : [40, 100],
        "Nup192" : [20, 100],
        "Nup42" : [70, 150],
        "Nup49" : [40, 100],
        "Nup53" : [20, 100],
        "Nup57" : [0, 75],
        "Nup59" : [40, 120],
        "Nup60" : [-200, -100],
        "Nup82.1" : [125, 295],     #"Nup82.1" : [145, 295],
        "Nup82.2" : [125, 295],     #"Nup82.2" : [145, 295],
        "Nup84" : [130, 170],       #"Nup84" : [150, 170],
        "Nup85" : [120, 200],       #"Nup85" : [140, 200],
        "Pom152" : [0, 95],
        "Pom34" : [0, 65],
        "Seh1" : [50, 170],
        "Sec13" : [0, 500]
    }
    print "\nZAxialPositionLowerRestraint !!"
    print "ZAxialPositionUpperRestraint !!\n"
    zaxial_weight = 1.0
    for protein, z in ZAXIAL.iteritems():
        if (protein not in nup_list_unique):
            continue

        zax = IMP.npc.npc_restraints.ZAxialPositionLowerRestraint(simo, protein, z[0], False, 1.0)
        zax.set_label('Zaxial_Position_Lower_%d_%s' % (z[0], protein))
        zax.set_weight(zaxial_weight)
        zax.add_to_model()
        outputobjects.append(zax)
        print (zax.get_output())

        zax = IMP.npc.npc_restraints.ZAxialPositionUpperRestraint(simo, protein, z[1], False, 1.0)
        zax.set_label('Zaxial_Position_Upper_%d_%s' % (z[1], protein))
        zax.set_weight(zaxial_weight)
        zax.add_to_model()
        outputobjects.append(zax)
        print (zax.get_output())

"""
membrane_surface_location_restraints = []
if MEMBRANE_SURFACE_LOCATION:
# Membrane-surface location
  for i, beads in enumerate(T.get_beads("Ndc1", 6)):
    msl = IMP.npc.MembraneSurfaceLocationRestraint(m, tor_R, tor_r, tor_th, 0.2)
    for b in beads:
      msl.add_particle(b)
    msl.set_name('Membrane_Surface_Location_Ndc1_%d' % i)
    all_restraints.add_restraint(msl)
    membrane_surface_location_restraints.append(msl)
  for i, beads in enumerate(T.get_beads("Pom152", 6)):
    msl = IMP.npc.MembraneSurfaceLocationRestraint(m, tor_R, tor_r, tor_th, 0.2)
    for b in beads:
      msl.add_particle(b)
    msl.set_name('Membrane_Surface_Location_Pom152_%d' % i)
    all_restraints.add_restraint(msl)
    membrane_surface_location_restraints.append(msl)
  for i, beads in enumerate(T.get_beads("Pom34", 6)):
    msl = IMP.npc.MembraneSurfaceLocationRestraint(m, tor_R, tor_r, tor_th, 0.2)
    for b in beads:
      msl.add_particle(b)
    msl.set_name('Membrane_Surface_Location_Pom34_%d' % i)
    all_restraints.add_restraint(msl)
    membrane_surface_location_restraints.append(msl)

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


#####################################################
# Restraints setup
# Distance_to_ponit restraints for orientation of the Nup84 complex
#####################################################
if (use_Distance_to_Point):
    dpr_weight = 1.0
    dpr_radius = 100.0

    dpr = IMP.pmi.restraints.basic.Distance_to_Point_Restraint(simo, tuple_selection=(230,230,"Nup133"), anchor_point=IMP.algebra.Vector3D(417.155, 195.314, 150.0), radius=dpr_radius, kappa=10.0)
    dpr.set_label("Nup133")
    dpr.set_weight(dpr_weight)
    dpr.add_to_model()
    outputobjects.append(dpr)
    print(dpr.get_output())

    dpr = IMP.pmi.restraints.basic.Distance_to_Point_Restraint(simo, tuple_selection=(324,324,"Nup85"), anchor_point=IMP.algebra.Vector3D(338.093, -170.387, 170.0), radius=dpr_radius, kappa=10.0)
    dpr.set_label("Nup85")
    dpr.set_weight(dpr_weight)
    dpr.add_to_model()
    outputobjects.append(dpr)
    print(dpr.get_output())

    dpr = IMP.pmi.restraints.basic.Distance_to_Point_Restraint(simo, tuple_selection=(465,465,"Nup120"), anchor_point=IMP.algebra.Vector3D(555.849, -159.949, 110.0), radius=dpr_radius, kappa=10.0)
    dpr.set_label("Nup120")
    dpr.set_weight(dpr_weight)
    dpr.add_to_model()
    outputobjects.append(dpr)
    print(dpr.get_output())
    print "Distance_to_Point_Restraint !!\n"


#####################################################
# Restraints setup
# Cross-link restraints using the whole NPC DSS XL data
#####################################################
sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "\nEVAL 0 : ", sf.evaluate(False), " (before XL restraints) - ", rank

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


#####################################################
# Restraints setup
# Distance restraints for homo-dimers
#####################################################
if (False):
    """
    dist_min = 3.0
    dist_max = 30.0
    dr_weight = 100.0

    dr1 = IMP.pmi.restraints.basic.DistanceRestraint(simo,(1417,1417,"Nup159.1"), (1417,1417,"Nup159.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr1.add_to_model()
    dr1.set_label("Nup159_1417-1417")
    dr1.set_weight(dr_weight)
    outputobjects.append(dr1)
    print(dr1.get_output())
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
# Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "\nEVAL 1 : ", sf.evaluate(False), " (initial) - ", rank

if (True):
    simo.optimize_floppy_bodies(150)
    print "\nEVAL 2 : ", sf.evaluate(False), " (after calling optimize_floppy_bodies(150)) - ", rank

    initial_nframes = 200
    if (use_XL):
        XL_restraints = [xl1]
    else:
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
                                        global_output_directory = "pre-EM3D_output",
                                        rmf_dir = "rmfs/",
                                        best_pdb_dir = "pdbs/",
                                        replica_stat_file_suffix = "stat_replica")
    mc1.execute_macro()
    rex1 = mc1.get_replica_exchange_object()
    print "\nEVAL 3 : ", sf.evaluate(False), " (after performing the pre-sampling) - ", rank
else:
    rex1 = None
    print "\n>> NO pre-sampling"

exit(0)

#####################################################
# Restraints setup
# EM 3D restraint using GMM
#####################################################
if (False):
    """
    # tail module em density
    mass = sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))
    gem = IMP.pmi.restraints.em.GaussianEMRestraint(resdensities,
                                                    '../data/em_density/spidervol_overlap.mrc.gmm.100.txt',
                                                    #'../data/em_density/emanvol.mrc.gmm.100.txt',
                                                    target_mass_scale=mass,
                                                    slope=0.000001,
                                                    target_radii_scale=3.0)
    gem.add_to_model()
    gem.set_weight(100.0)
    #gem.set_weight(200.0)
    #gem.center_model_on_target_density(simo)
    outputobjects.append(gem)


    # tail module em density
    mass2 = sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))
    gem2 = IMP.pmi.restraints.em.GaussianEMRestraint(resdensities,
                                                    #'../data/em_density/spidervol_overlap.mrc.gmm.100.txt',
                                                    '../data/em_density/emanvol.mrc.gmm.100.txt',
                                                    target_mass_scale=mass2,
                                                    slope=0.000001,
                                                    target_radii_scale=3.0)
    gem2.add_to_model()
    gem2.set_weight(100.0)
    #gem2.set_weight(200.0)
    #gem2.center_model_on_target_density(simo)
    outputobjects.append(gem2)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
    print "\nEVAL 4 : ", sf.evaluate(False), " (after applying the EM 3D restraint) - ", rank
    """

#####################################################
# Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
mc2 = IMP.pmi.macros.ReplicaExchange0(m,
                                    simo,
                                    monte_carlo_sample_objects = sampleobjects,
                                    output_objects = outputobjects,
                                    #crosslink_restraints = [xl1, xl2],
                                    crosslink_restraints = [xl1],
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 2.5,
                                    number_of_best_scoring_models = 100,
                                    #number_of_best_scoring_models = int(inputs.nrepeats)-200,
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
                                    replica_exchange_object = rex1)
mc2.execute_macro()
print "\nEVAL 5 : ", sf.evaluate(False), " (final evaluation) - ", rank
