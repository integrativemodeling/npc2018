#!/usr/bin/env python
#####################################################
# Last Update: Feb 8, 2015
# by Seung Joong Kim
# at Andrej Sali group, University of California San Francisco (UCSF)
#####################################################
import IMP
import IMP.core
#import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.em2d
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.proteomics
import IMP.pmi.representation
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.output
import IMP.pmi.samplers
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
parser.add_argument('-refine', action="store", dest="refinement", help="refinement True or False (XL distance restraints for dimers)" )
parser.add_argument('-weight', action="store", dest="weight", help="weight for XL distance restaints for dimer" )
parser.add_argument('-res_cry', action="store", dest="res_cry", help="resolution of the crystal structures" )
parser.add_argument('-res_hom', action="store", dest="res_hom", help="resolution of the comparative (homology) models" )
parser.add_argument('-res_ev', action="store", dest="res_ev", help="resolution of the excluded volume restraints" )
parser.add_argument('-res_compo', action="store", dest="res_compo", help="resolution of the composite restraints" )
parser.add_argument('-draw_hierarchy', action="store", dest="draw_hierarchy", help="draw hierarchy" )
inputs = parser.parse_args()

# Setting up the input parameters
if inputs.ncopy==None:
    inputs.ncopy = "2"
if (inputs.symmetry=="True") or (inputs.symmetry=="true") or (inputs.symmetry=="Yes") or (inputs.symmetry=="yes") :
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

if inputs.XL_input==None:
    inputs.XL_input = "../data/XL.csv"
else:
    f=open(inputs.XL_input,"r")
    f.close()
if inputs.nrepeats==None:
    inputs.nrepeats = 1000
if inputs.folder_output==None:
    inputs.folder_output = "output"
if inputs.rmf_output==None:
    inputs.rmf_output = "models.rmf"
if inputs.stat_output==None:
    inputs.stat_output = "stat.dat"
if (inputs.refinement=="True") or (inputs.refinement=="true") or (inputs.refinement=="Yes") or (inputs.refinement=="yes") :
    inputs.refinement = True
else:
    inputs.refinement = False
if inputs.weight==None:
    inputs.weight = 100.0

if inputs.res_cry==None:
    inputs.res_cry = 1.0
if inputs.res_hom==None:
    inputs.res_hom = 5.0
if inputs.res_ev==None:
    inputs.res_ev = 1.0
if inputs.res_compo==None:
    inputs.res_compo = 100.0
if (inputs.draw_hierarchy=="True") or (inputs.draw_hierarchy=="true") or (inputs.draw_hierarchy=="Yes") or (inputs.draw_hierarchy=="yes") :
    inputs.draw_hierarchy = True
else:
    inputs.draw_hierarchy = False
print inputs


#####################################################
# Create hierarchies and rigid bodies and flexible parts
# for bead representations
#####################################################
m = IMP.Model()
simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)
#simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=True)


#####################################################
# setting up parameters
#####################################################
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print "rank = ", rank

#rbmaxtrans = 0.5
#fbmaxtrans = 0.5
#rbmaxrot=0.02
rbmaxtrans = 2.00
fbmaxtrans = 3.00
rbmaxrot=0.04
#rbmaxtrans = 1.50
#fbmaxtrans = 1.50
#rbmaxrot = 0.025
outputobjects = []
sampleobjects = []
partialscore1 = []
partialscore2 = []

beadsize = int(inputs.res_cry)
beadsize5 = 5
beadsize10 = 10
beadsize20 = 20
beadsize50 = 50
beadsize100 = 100
beadsize_det = 1
det_ini = 9000
det_pos = det_ini + beadsize_det
res_ev = int(inputs.res_ev)
res_conn = int(inputs.res_cry)
res_str = int(inputs.res_cry)
em_weight = float(inputs.weight)
sc_scale = 0.1

datadirectory = "../data/"
fasta_files = "../data/FASTA_"
pdb_files = "../data/PDB_"


#####################################################
## REPRESENTATION
#####################################################
# compname  hier_name      color  fasta_file                  fasta_id    pdbname                       chain res_range       read_em_files  bead_size  rigid_body super_rigid_body em_num_components em_txt_file_name em_mrc_file_name chain_of_super_rb
domains = \
[
# ("pom152",  "detgnt_99",  1.0,  fasta_files+"detgnt.txt",    "detgnt",   "BEADS",                      " ",  (det_ini,det_pos,0), True,     beadsize_det,   99,     [11],          1,                None,            None, [0]),

 ("pom152",  "pom152_11a", 0.0,  fasta_files+"pom152.txt",    "pom152",   "BEADS",                      " ",  (1,100,0),      True,          beadsize20,     111,    [11],          1,                None,            None, [0]),
 ("pom152",  "pom152_11b", 1.0,  fasta_files+"pom152.txt",    "pom152",   "BEADS",                      " ",  (101,200,0),    True,          beadsize100,    112,    [11],          1,                None,            None, [0]),
 ("pom152",  "pom152_11c", 0.0,  fasta_files+"pom152.txt",    "pom152",   "BEADS",                      " ",  (201,378,0),    True,          beadsize20,     113,    [11],          1,                None,            None, [0]),

# ("pom152",  "pom152_1",   0.1,  fasta_files+"pom152.txt",    "pom152",   "BEADS",                      " ",  (379,472,0),    True,          beadsize50,     1,      [1],           1,                None,            None, [0]),
 ("pom152",  "pom152_1",   0.1,  fasta_files+"pom152.txt",    "pom152",   pdb_files+"375_482.pdb",      "A",  (379,472,0),    True,          beadsize,       1,      [1],           4,                None,            None, [0]),
# ("pom152",  "pom152_12",  0.1,  fasta_files+"pom152.txt",    "pom152",   "BEADS",                      " ",  (473,519,0),    True,          beadsize50,     12,     [1],           1,                None,            None, [0]),
 ("pom152",  "pom152_12",  0.1,  fasta_files+"pom152.txt",    "pom152",   "BEADS",                      " ",  (473,519,0),    True,          beadsize10,     12,     [1],           1,                None,            None, [0]),

 ("pom152",  "pom152_2",   0.2,  fasta_files+"pom152.txt",    "pom152",   pdb_files+"516_611.pdb",      "A",  (520,611,0),    True,          beadsize,       2,      [2],           4,                None,            None, [0]),
 ("pom152",  "pom152_22",  0.2,  fasta_files+"pom152.txt",    "pom152",   "BEADS",                      " ",  (612,615,0),    True,          beadsize10,     22,     [2],           1,                None,            None, [0]),

 ("pom152",  "pom152_3",   0.3,  fasta_files+"pom152.txt",    "pom152",   pdb_files+"27005.pdb",        "A",  (616,714,0),    True,          beadsize,       3,      [3],           4,                None,            None, [0]),
 ("pom152",  "pom152_32",  0.3,  fasta_files+"pom152.txt",    "pom152",   "BEADS",                      " ",  (715,721,0),    True,          beadsize10,     32,     [3],           1,                None,            None, [0]),

 ("pom152",  "pom152_4",   0.4,  fasta_files+"pom152.txt",    "pom152",   pdb_files+"27005.pdb",        "A",  (722,818,0),    True,          beadsize,       4,      [4],           4,                None,            None, [0]),
 ("pom152",  "pom152_42",  0.4,  fasta_files+"pom152.txt",    "pom152",   "BEADS",                      " ",  (819,823,0),    True,          beadsize10,     42,     [4],           1,                None,            None, [0]),

 ("pom152",  "pom152_5",   0.5,  fasta_files+"pom152.txt",    "pom152",   pdb_files+"26996.pdb",        "A",  (824,918,0),    True,          beadsize,       5,      [5],           4,                None,            None, [0]),
 ("pom152",  "pom152_52",  0.5,  fasta_files+"pom152.txt",    "pom152",   "BEADS",                      " ",  (919,930,0),    True,          beadsize10,     52,     [5],           1,                None,            None, [0]),

 ("pom152",  "pom152_6",   0.6,  fasta_files+"pom152.txt",    "pom152",   pdb_files+"26996.pdb",        "A",  (931,1026,0),   True,          beadsize,       6,      [6],           4,                None,            None, [0]),
 ("pom152",  "pom152_62",  0.6,  fasta_files+"pom152.txt",    "pom152",   "BEADS",                      " ",  (1027,1035,0),  True,          beadsize10,     62,     [6],           1,                None,            None, [0]),

 ("pom152",  "pom152_7",   0.7,  fasta_files+"pom152.txt",    "pom152",   pdb_files+"26996.pdb",        "A",  (1036,1141,0),  True,          beadsize,       7,      [7],           4,                None,            None, [0]),
 ("pom152",  "pom152_72",  0.7,  fasta_files+"pom152.txt",    "pom152",   "BEADS",                      " ",  (1142,1149,0),  True,          beadsize10,     72,     [7],           1,                None,            None, [0]),
                                                                                                                                                                                                                                
 ("pom152",  "pom152_8",   0.8,  fasta_files+"pom152.txt",    "pom152",   pdb_files+"1146_1237.pdb",    "A",  (1150,1229,0),  True,          beadsize,       8,      [8],           4,                None,            None, [0]),
 ("pom152",  "pom152_82",  0.8,  fasta_files+"pom152.txt",    "pom152",   "BEADS",                      " ",  (1230,1243,0),  True,          beadsize10,     82,     [8],           1,                None,            None, [0]),
                                                                                                                                                                                                                                
 ("pom152",  "pom152_9",   0.9,  fasta_files+"pom152.txt",    "pom152",   pdb_files+"1238_1337.pdb",    "A",  (1244,1337,0),  True,          beadsize,       9,      [9],           4,                None,            None, [0]),
]

bm1=IMP.pmi.macros.BuildModel1(simo)
#bm1.set_gmm_models_directory(datadirectory+"em_gmm_model/")

if (inputs.rmf_input is not None) :
    dom = set([s[0] for s in domains])
    for d in list(dom):
        bm1.set_rmf_file(d, inputs.rmf_input, int(inputs.rmf_frame_number))

bm1.build_model(data_structure = domains, sequence_connectivity_scale=sc_scale, sequence_connectivity_resolution=1.0)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "\nEVAL 1 : ", sf.evaluate(False), " (before applying the EM 2D restraint) - ", rank


ps = []
resolution = 1.0
"""
################################
# Detergent
################################
particles = IMP.pmi.tools.select(simo,
                                resolution=resolution,
                                name="detgnt")
for p in particles:
    ps.append(p)
    print p
    print IMP.core.XYZR(p).get_radius()
    print IMP.atom.Mass(p).get_mass()
"""

################################
# Pom152
################################
particles = IMP.pmi.tools.select(simo,
                                resolution=resolution,
                                name="pom152")
for p in particles:
    ps.append(p)
    print p
    print IMP.core.XYZR(p).get_radius()
    print IMP.atom.Mass(p).get_mass()
    

print "len(ps) = ", len(ps)   #len(particles_all)=907
################################

# writing the particle densities
map = IMP.em.SampledDensityMap(ps, 30, 5.0)
IMP.em.write_map(map, "em2d_particle_selected.mrc")

# writing the pdb files
output = IMP.pmi.output.Output()
output.init_pdb("test_pdb_writing.pdb", simo.prot)
output.write_pdb("test_pdb_writing.pdb")

images = [inputs.em2d_input]
pixel_size = 2.03
image_resolution = 30.0
projection_number = 500

"""
import RMF
m = IMP.Model()
inf = RMF.open_rmf_file_read_only(inputs.rmf_input)
h = IMP.rmf.create_hierarchies(inf, m)[0]
ps = IMP.core.get_leaves(h)
IMP.rmf.load_frame(inf, 0)
children = h.get_children()
"""

rs = IMP.RestraintSet(m, 'em2d')

em2d = IMP.em2d.PCAFitRestraint(
    ps, images, pixel_size, image_resolution, projection_number, True)
rs.add_restraint(em2d)

"""
sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "\nEVAL 3 : ", sf.evaluate(False), " (before applying the EM 2D restraint) - ", rank

em2d = IMP.pmi.restraints.em2d.ElectronMicroscopy2D_FFT(simo,
                                                    images,
                                                    resolution = resolution,
                                                    pixel_size = pixel_size,
                                                    image_resolution = image_resolution,
                                                    projection_number = projection_number)
em2d.add_to_model()
#outputobjects.append(em2d)

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "\nEVAL 4 : ", sf.evaluate(False), " (after applying the EM 2D restraint) - ", rank
exit(1)
"""

m.update()
#score = em2d_weight*rs.unprotected_evaluate(None)
score = em2d.unprotected_evaluate(None)
print("PCAFitRestraint score = ", score)
print("CCC = ", math.exp(-float(score)))
em2d.write_best_projections("best_projections.pgm")

sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "\nEVAL 2 : ", sf.evaluate(False), " (after applying the EM 2D restraint) - ", rank
