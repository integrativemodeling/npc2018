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
import IMP.npc.npc_restraints
import random
import os
import math
from math import pi,log,sqrt

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

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
#s = IMP.pmi1.topology.System(m)
#st = s.create_state()
simo = IMP.pmi1.representation.Representation(m,upperharmonic=True,disorderedlength=False)

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
n145N_pdbfile  = npc + "Nup145N_3kep_AB_459_605r_ca.pdb"

#Dbp5_pdbfile   = npc + "Dbp5_3rrm_A_91_482.pdb"    # Not using it for now
Gle1N_pdbfile  = npc + "Gle1_4wij_A_121_239ca.pdb"
Gle1C_pdbfile  = npc + "Gle1_3rrm_B_244_538ca.pdb"
Gle2_pdbfile   = npc + "Gle2_3mmy_A_4_362ca.pdb"

#####################################################
# Parameters for Debugging
#####################################################
is_n84 = False
is_n82 = False
is_nic96 = False
is_inner_ring = False
is_membrane = True
is_cytoplasm = False
is_nucleoplasm = False
is_basket = False
is_FG = False

use_neighboring_spokes = True
#Stopwatch_None_delta_seconds  20~25 (1 spoke) / 60-70 sec (3 spokes)
#Stopwatch_None_delta_seconds  22~27 (1 spoke) / 65-75 sec (3 spokes) with XL
#Stopwatch_None_delta_seconds  42~47 (1 spoke) / 90~100 sec (3 spokes) with XL + EM
use_shuffle = False
use_ExcludedVolume = True
use_Immuno_EM = False
use_sampling_boundary = True
use_XL = False
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
# Membrane nups
##########################
if (is_membrane):
    #domains.append(("Pom152",    "Pom152_1" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", "BEADS",      " ", (   1, 378,0),  gmm,  beadsize50,  1521, [152],0,  None,  None, None, False))
    domains.append(("Pom152",    "Pom152_1" ,     1.0,  f_npc+"Pom152.txt", "YMR129W", pom152_pdb,   "A", (   1, 496,0),  gmm,  beadsize50,  1521, [152],2,  " ",   " ",  None, False))
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
        #domains.append(("Pom152@%d"%i, "Pom152_1@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", "BEADS",     " ", (   1, 378,0), gmm_c,  beadsize50,  None, None, 0,  None,                 None,                 None, False))
        domains.append(("Pom152@%d"%i, "Pom152_1@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", (   1, 496,0), gmm_c,  beadsize50,  None, None, 2,  gmm_f+"Pom152_1.txt", gmm_f+"Pom152_1.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_2@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", ( 497, 613,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_2.txt", gmm_f+"Pom152_2.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_3@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", ( 614, 718,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_3.txt", gmm_f+"Pom152_3.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_4@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", ( 719, 821,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_4.txt", gmm_f+"Pom152_4.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_5@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", ( 822, 924,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_5.txt", gmm_f+"Pom152_5.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_6@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", ( 925,1031,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_6.txt", gmm_f+"Pom152_6.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_7@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", (1032,1145,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_7.txt", gmm_f+"Pom152_7.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_8@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", (1146,1236,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_8.txt", gmm_f+"Pom152_8.mrc", None, False))
        domains.append(("Pom152@%d"%i, "Pom152_9@%d"%i, 1.0, f_npc+"Pom152.txt", "YMR129W", pom152_pdb,  "A", (1237,1337,0), gmm_c,  beadsize100, None, None, 2,  gmm_f+"Pom152_9.txt", gmm_f+"Pom152_9.mrc", None, False))

#####################################################
# Model Building
#####################################################
bm1 = IMP.pmi1.macros.BuildModel1(simo)
bm1.set_gmm_models_directory(gmm_f)

if (True):
    if (is_membrane):
        #bm1.set_rmf_file('Pom152', "../data_npc/Pom152_rmfs/Pom152_0_final.rmf3", 0)
        bm1.set_rmf_file('Pom152', "../data_npc/Pom152_rmfs/Pom152_new_final.rmf3", 0)

# remove connectivity for clones
clone_list = [entry[0] for entry in domains if '@' in entry[0]]
clone_list_unique = sorted(list(set(clone_list)))   # Make a unique list
print ("clone_list_unique = ", clone_list_unique)

bm1.build_model(data_structure = domains, sequence_connectivity_scale=0.1, sequence_connectivity_resolution=1.0,
                skip_connectivity_these_domains=clone_list_unique, skip_gaussian_in_rmf=False, skip_gaussian_in_representation=False)
#exit(0)
bm1.scale_bead_radii(100, 0.6)


#####################################################
# apply the rotational symmetry
#####################################################
if (use_neighboring_spokes):
    if (is_membrane):
        for protein in ['Pom152']:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
            simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
            simo.create_rotational_symmetry(protein+'@11', [protein+'@%d'%i for i in range(12,14)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
else:
    if (is_membrane):
        for protein in ['Pom152']:
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
for protein in [(248,360,'Nup53'), (266,402,'Nup59'), (379,1337,'Pom152'), 'Mlp1', 'Mlp2']:
    rigid_tuples.append(protein)
# Remove flexible movers for all clones
for protein in clone_list_unique:
    rigid_tuples.append(protein)

print ("rigid_tuples = ", rigid_tuples)
for rt in rigid_tuples:
    hs = IMP.pmi1.tools.select_by_tuple(simo,rt)
    simo.remove_floppy_bodies(hs)


#####################################################
# randomize the initial configuration
#####################################################
if (use_shuffle) :
    #simo.shuffle_configuration(max_translation=1, avoidcollision=False, ignore_initial_coordinates=True)
    simo.shuffle_configuration(bounding_box=((400, -150, 0), (650, 150, 150)), ignore_initial_coordinates=True, cutoff=1.0, niterations=1000)


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
        "Pom152"   : [470, 630]     #(Table S2 and Table S7 are different)
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
        "Pom152"   : [   5,  60]        #"Pom152"   : [   0,  95]   (Table S2 and Table S7 are different)
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
# Restraints setup
# Distance restraints
#####################################################
dist_min = 3.0
dr_weight = 10.0


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
msl_weight  = 10.0

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

# The Pom152 ring
if (is_membrane):
    dist_max = 20.0

    # same residue cross-link of Pom152 62-62
    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(62,62,"Pom152"), (62,62,"Pom152@11"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Pom152-Pom152@11")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(301,301,"Pom152"), (301,301,"Pom152@11"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Pom152-Pom152@11_301")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(351,351,"Pom152"), (351,351,"Pom152@11"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.add_to_model()
    dr.set_label("Pom152-Pom152@11_351")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    if (use_neighboring_spokes):
        dist_max = 23.0
        # TODO: Pom152 orientation?  (clockwise or counter-clockwise?)
        #dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(351,351,"Pom152"), (1337,1337,"Pom152@12"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(351,351,"Pom152"), (1337,1337,"Pom152@13"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr.add_to_model()
        #dr.set_label("Pom152-Pom152@12")
        dr.set_label("Pom152-Pom152@13")
        dr.set_weight(dr_weight)
        outputobjects.append(dr)
        print(dr.get_output())


    pom152_min = 2.5;     pom152_max = 15;
    POM_LIST = [ 379, 392, 472, 520, 611, 616, 714, 722, 818, 824, 866, 918, 931, 1010, 1026, 1036, 1064, 1141, 1150, 1182, 1229, 1244, 1282, 1337 ]
    for z in POM_LIST:
        zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (z, z, "Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
        zax.set_label('Lower_%d_Upper_%d_Pom152_%s' % (pom152_min, pom152_max, z))
        zax.set_weight(zaxial_weight)
        zax.add_to_model()
        outputobjects.append(zax)
        print (zax.get_output())

    """
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
    """


#####################################################
# Restraints setup - Membrane Exclusion
#####################################################
tor_th      = 150.0 - tor_th_ALPS
tor_R       = 390.0 + 150.0
tor_r       = tor_th/2.0
mex_sigma   = 0.2
mex_weight  = 10.0

if (is_membrane):
    nup_list = [entry[0] for entry in domains if not '@' in entry[0]]
    nup_list_unique = sorted(list(set(nup_list)))   # Make a unique list

    MEX_LIST = [
        [1, 110, "Pom152"],
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
    ids_map = IMP.pmi1.tools.map()
    ids_map.set_map_element(1.0, 1.0)

    xl1 = IMP.pmi1.restraints.crosslinking.ISDCrossLinkMS(simo,
                                                        '../data_npc/XL_optimized_ambiguity.csv',
                                                        length = 26.0,
                                                        slope = 0.00,
                                                        columnmapping = columnmap,
                                                        ids_map = ids_map,
                                                        resolution = 1.0,
                                                        inner_slope = 0.01,
                                                        filelabel = "wtDSS",
                                                        label = "wtDSS",
                                                        attributes_for_label = ["XLUniqueID"],
                                                        csvfile = True)
    xl1.add_to_model()
    xl1.set_weight(2.0)        # play with the weight
    sampleobjects.append(xl1)
    outputobjects.append(xl1)
    xl1.set_psi_is_sampled(False)
    psi2 = xl1.get_psi(1.0)[0]
    psi2.set_scale(0.05)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
    print "\nEVAL 0 : ", sf.evaluate(False), " (after applying the XL restraint) - ", rank
    XL_restraints = [xl1]
else:
    XL_restraints = None


#####################################################
# Sampling Boundary Restraint
#####################################################
if (use_sampling_boundary):
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

    #resdensities = bm1.get_density_hierarchies([t[1] for t in domains])
    resdensities = bm1.get_density_hierarchies(main_spoke_hier_name)
    print ("resdensities=", resdensities)       ####  TODO: make sure resdensities are correct

    mass = sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))
    mass *= 1.2 * 2.0           # 1.2 for adjustment of the GMM (after removing flexible GMMs)
    print ("Total mass for the Sampling Boundary EM restraint = ", mass)
    sbr = IMP.pmi1.restraints.em.GaussianEMRestraint(resdensities,
                                                    '../data_npc/em_gmm_model/SJ_Pom152.gmm.15.txt',
                                                    target_mass_scale=mass,
                                                    slope=0.01,
                                                    #slope=0.0000001,
                                                    target_radii_scale=3.0)
    sbr.add_to_model()
    sbr.set_weight(0.5)        # play with the weight
    sbr.set_label("Sampling_Boundary")
    #sbr.center_model_on_target_density(simo)
    outputobjects.append(sbr)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
    print "\nEVAL 1 : ", sf.evaluate(False), " (after applying the Sampling Boundary EM restraint) - ", rank


#####################################################
# 1st Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
print "\nEVAL 2 : ", sf.evaluate(False), " (initial) - ", rank

if (False):
    simo.optimize_floppy_bodies(300)
    print "\nEVAL 2 : ", sf.evaluate(False), " (after calling optimize_floppy_bodies(300)) - ", rank

    #XL_restraints = None
    mc1 = IMP.pmi1.macros.ReplicaExchange0(m,
                                        simo,
                                        monte_carlo_sample_objects = sampleobjects,
                                        output_objects = outputobjects,
                                        crosslink_restraints = XL_restraints,
                                        monte_carlo_temperature = 1.0,
                                        replica_exchange_minimum_temperature = 1.0,
                                        replica_exchange_maximum_temperature = 2.5,
                                        number_of_best_scoring_models = 0,
                                        monte_carlo_steps = 10,
                                        number_of_frames = 200,
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
else:
    rex1 = None
sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
print "\nEVAL 3 : ", sf.evaluate(False), " (after performing the pre_sampling) - ", rank

#####################################################
# Restraints setup
# EM 3D restraint using GMM
#####################################################
if (use_EM3D):
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

    #resdensities = bm1.get_density_hierarchies([t[1] for t in domains])
    resdensities = bm1.get_density_hierarchies(main_spoke_hier_name)
    print ("resdensities=", resdensities)       ####  TODO: make sure resdensities are correct

    mass = sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))
    mass *= 1.2 * 2.0           # 1.2 for adjustment of the GMM (after removing flexible GMMs)
    print ("Total mass for the EM 3D restraint = ", mass)
    gem = IMP.pmi1.restraints.em.GaussianEMRestraint(resdensities,
                                                    '../data_npc/em_gmm_model/SJ_Pom152.gmm.15.txt',
                                                    target_mass_scale=mass,
                                                    slope=0.0000001,
                                                    target_radii_scale=3.0)
    gem.add_to_model()
    gem.set_weight(2000.0)        # play with the weight
    #gem.center_model_on_target_density(simo)
    outputobjects.append(gem)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
    print "\nEVAL 4 : ", sf.evaluate(False), " (after applying the EM 3D restraint) - ", rank


#####################################################
# 2nd Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
if (False):
    mc2 = IMP.pmi1.macros.ReplicaExchange0(m,
                                        simo,
                                        monte_carlo_sample_objects = sampleobjects,
                                        output_objects = outputobjects,
                                        crosslink_restraints = XL_restraints,
                                        monte_carlo_temperature = 1.0,
                                        replica_exchange_minimum_temperature = 1.0,
                                        replica_exchange_maximum_temperature = 2.5,
                                        number_of_best_scoring_models = 0,
                                        monte_carlo_steps = 10,
                                        number_of_frames = 3000,
                                        write_initial_rmf = True,
                                        initial_rmf_name_suffix = "initial",
                                        stat_file_name_suffix = "stat",
                                        best_pdb_name_suffix = "model",
                                        do_clean_first = True,
                                        do_create_directories = True,
                                        global_output_directory = "2_XL_EM_output",
                                        rmf_dir = "rmfs/",
                                        best_pdb_dir = "pdbs/",
                                        replica_stat_file_suffix = "stat_replica",
                                        replica_exchange_object = rex1)
    mc2.execute_macro()
    rex2 = mc2.get_replica_exchange_object()
else:
    rex2 = None
sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
print "\nEVAL 5 : ", sf.evaluate(False), " (after performing the XL_EM_sampling) - ", rank


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

    ev1 = IMP.pmi1.restraints.stereochemistry.ExcludedVolumeSphere(simo,
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
        ev2 = IMP.pmi1.restraints.stereochemistry.ExcludedVolumeSphere(simo,
                                                                     included_objects = included_objects,
                                                                     other_objects = other_objects,
                                                                     resolution = res_ev)
        ev2.add_to_model()
        ev2.set_label('bipartite')
        ev2.set_weight(10000.0)
        outputobjects.append(ev2)
        print(ev2.get_output())
        print "ExcludedVolumeSphere2 between the main spoke and the neighboring spokes !!\n"


#####################################################
# Restraints setup
# Distance restraints between neighboring Ig domains
# Domain connectivity
#####################################################
if (True):
    dr_weight = 100.0
    IG_LIST = [ [471,520], [611,616], [713,722], [816,824], [916,933], [1025,1038], [1141,1150], [1229,1244] ]
    for z in IG_LIST:
        if (z[0] == 611):
            dist_min = 12.0
        else:
            dist_min = 9.0

        dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo, (z[0],z[0],"Pom152"), (z[1],z[1],"Pom152"), distancemin=dist_min, distancemax=dist_min+5.0, resolution=1)
        dr.add_to_model()
        dr.set_label('Pom152_%d_%d' % (z[0], z[1]))
        dr.set_weight(dr_weight)
        outputobjects.append(dr)
        print(dr.get_output())


#####################################################
"""Restrain the dihedral between planes defined by three particles.

This restraint is useful for restraining the twist of a string of
more or less identical rigid bodies, so long as the curvature is mild.
"""
#####################################################
class PlaneDihedralRestraint(object):
    def __init__(self, representation = None, particle_triplets = None, angle=0., k=1., hier = None):
        """Constructor
        @param particle_triplets List of selection_tuple triplets. Each triplet
                                 defines a plane. Dihedrals of adjacent planes
                                 in list are scored.
        @param angle Angle of plane dihedral in degrees
        @param k Strength of restraint
        @param label Label for output
        @param weight Weight of restraint
        \note Particles defining planes should be rigid and more or less
              parallel for proper behavior
        """

        # PMI1/2 selection
        if representation is None and hier is not None:
            self.m = hier.get_model()
        elif hier is None and representation is not None:
            self.m = representation.prot.get_model()
        else:
            raise Exception("PlaneDihedralRestraint: must pass hier or representation")

        #self.m = particle_triplets[0][0].get_model()
        #super(PlaneDihedralRestraint, self).__init__(m, label=label,
        #                                             weight=weight)

        self.rs = IMP.RestraintSet(self.m, 'PlaneDihedralRestraint')
        self.weight = 1.0
        self.label = "None"

        angle = pi * angle / 180.
        ds = IMP.core.Cosine(.5 * k, 1., -angle)
        for i, t1 in enumerate(particle_triplets[:-1]):
            if len(t1) != 3:
                raise ValueError("wrong length of quadruplet")
            t2 = particle_triplets[i + 1]

            ps1 = []
            for selection_tuple in t1:
                #print (selection_tuple)
                p = IMP.pmi1.tools.select_by_tuple(representation, selection_tuple, resolution=1)
                #print(IMP.atom.Residue.get_is_setup(p[0]))
                ps1.append(p[0])
            print (ps1)

            ps2 = []
            for selection_tuple in t2:
                #print (selection_tuple)
                p = IMP.pmi1.tools.select_by_tuple(representation, selection_tuple, resolution=1)
                #print(IMP.atom.Residue.get_is_setup(p[0]))
                ps2.append(p[0])
            print (ps2)

            q1 = [ps1[1], ps1[0], ps2[0], ps2[1]]
            q2 = [ps1[2], ps1[0], ps2[0], ps2[2]]
            self.rs.add_restraint(IMP.core.DihedralRestraint(self.m, ds, q1[0], q1[1], q1[2], q1[3]))
            self.rs.add_restraint(IMP.core.DihedralRestraint(self.m, ds, q2[0], q2[1], q2[2], q2[3]))

    def set_label(self, label):
        self.label = label

    def add_to_model(self):
        IMP.pmi1.tools.add_restraint_to_model(self.m, self.rs)

    def get_restraint(self):
        return self.rs

    def set_weight(self, weight):
        self.weight = weight
        self.rs.set_weight(self.weight)

    def get_output(self):
        self.m.update()
        output = {}
        score = self.weight * self.rs.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["PlaneDihedralRestraint_" + self.label] = str(score)
        return output

    def evaluate(self):
        return self.weight * self.rs.unprotected_evaluate(None)


#####################################################
# Restrain the dihedral between planes defined by three particles.
#####################################################
if (True):
    pdr_weight = 1000.0
    PDR_LIST = [ [388,398,447], [530,540,589], [624,635,689], [730,740,793], [834,844,894], \
                 [941,951,1001], [1045,1055,1116], [1155,1163,1206], [1251,1261,1319] ]
    for i, t1 in enumerate(PDR_LIST[:-1]):
        t2 = PDR_LIST[i + 1]

        ps1=[]
        for z in t1:
            ps1.append((z,z,"Pom152"))
        ps2=[]
        for z in t2:
            ps2.append((z,z,"Pom152"))
        print (ps1, ps2)

        pdr = PlaneDihedralRestraint(simo, [ps1, ps2], angle=90.0, k=1.)
        pdr.add_to_model()
        pdr.set_label('Pom152_pdr_%d_%d' % (t1[0], t2[0]))
        pdr.set_weight(pdr_weight)
        outputobjects.append(pdr)
        print(pdr.get_output())


#####################################################
# 4th Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
mc4 = IMP.pmi1.macros.ReplicaExchange0(m,
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
                                    #replica_exchange_object = rex3)
mc4.execute_macro()
sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
print "\nEVAL 8 : ", sf.evaluate(False), " (final evaluation) - ", rank
exit(0)
