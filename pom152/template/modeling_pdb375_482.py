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
#import IMP.pmi.representation
import representation_pom152
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.tools
import IMP.pmi.output
import IMP.pmi.samplers
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
#simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)
simo = representation_pom152.Representation(m,upperharmonic=True,disorderedlength=False)
#simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=True)
"""
    import IMP.pmi.topology
    import sys
    
    # each list contains list of domain names (from topology) that move together
    #  flexible beads are automatically added to missing regions and sampled
    rigid_bodies = [["pom152_11"]]
    super_rigid_bodies = [["pom152_11"]]
    chain_of_super_rigid_bodies = [["pom152_11"]]
                                   
    # Create list of components from topology file
    topology = IMP.pmi.topology.TopologyReader("./topology.txt")
    domains = topology.component_list
    
    
    bm = IMP.pmi.macros.BuildModel(m,
                        component_topologies=domains,
                        list_of_rigid_bodies=rigid_bodies,
                        list_of_super_rigid_bodies=super_rigid_bodies,
                        chain_of_super_rigid_bodies=chain_of_super_rigid_bodies)
    representation = bm.get_representation()
    exit()
"""
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

# rigid body movement params
rb_max_trans = 4.00
rb_max_rot = 0.04

# flexible bead movement
bead_max_trans = 4.00

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
resdensities=bm1.get_density_hierarchies([t[1] for t in domains])
#print resdensities

#bm1.scale_bead_radii(40,0.8)
#resdensities = bm1.get_density_hierarchies([t[1] for t in domainxl_cliques_psi = 0.25s])
#print resdensities; exit()

#model_ps = []
#for h in self.densities:
#    model_ps += IMP.atom.get_leaves(h)

"""
#####################################################
# randomize the initial configuration
#####################################################
if (inputs.rmf_input is None) :
    #simo.shuffle_configuration(50)
    simo.shuffle_configuration(250)
"""
#####################################################
# defines the movers
#####################################################
# Add default mover parameters to simulation
simo.set_rigid_bodies_max_rot(rb_max_rot)
simo.set_floppy_bodies_max_trans(bead_max_trans)
simo.set_rigid_bodies_max_trans(rb_max_trans)

### These two below are already executed in BuildModel1
#simo.set_floppy_bodies()
#simo.setup_bonds()

#prot = simo.prot
outputobjects.append(simo)
sampleobjects.append(simo)


#####################################################
# Create rigid bodies and flexible parts
# for bead representations (OLDER VERSION)
#####################################################
if (False) : 
    """
    simo.set_rigid_body_from_hierarchies(PiNup53_1)
    simo.set_super_rigid_body_from_hierarchies(PiNup53_1)
    #simo.set_rigid_body_from_hierarchies(PiNup53_1 + PiNup53_2)
    #simo.set_super_rigid_body_from_hierarchies(PiNup53_1 + PiNup53_2)

    simo.set_rigid_bodies_max_rot(rbmaxrot)
    simo.set_floppy_bodies_max_trans(fbmaxtrans)
    simo.set_rigid_bodies_max_trans(rbmaxtrans)
    simo.set_floppy_bodies()
    simo.setup_bonds()

    #prot = simo.prot
    outputobjects.append(simo)
    sampleobjects.append(simo)
    """


#####################################################
# Restraints setup
# Excluded Volume restraint
#####################################################
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo, resolution = res_ev)
ev.add_to_model()
outputobjects.append(ev)
print(ev.get_output())
print "ExcludedVolumeSphere !!\n"


#####################################################
# Restraints setup
# External Barrier restraint
#####################################################
eb = IMP.pmi.restraints.basic.ExternalBarrier(simo, radius = 1000)
eb.add_to_model()
outputobjects.append(eb)
print(eb.get_output())
print "ExternalBarrier !!\n"


#####################################################
# Restraints setup
# Distance restraints for the transmembrane domain
#####################################################
if (False):
    dist_min = 1.0
    dist_max = 12.5
    dr_weight = 100

    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(124,124,"pom152"), (det_pos,det_pos,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_TM1")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(156,156,"pom152"), (det_pos,det_pos,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_TM2")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(185,185,"pom152"), (det_pos,det_pos,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_TM3")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    print "Distance Restraints between the transmembrane domain and the detergent michelle !!"
    print "dr_weight = ", dr_weight, "dist_min = ", dist_min, "dist_max = ", dist_max, "\n"

if (True):
    dist_min = 240.0
    dist_max = 420.0
    dr_weight = 100
    
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(379,379,"pom152"), (1337,1337,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_elongated")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    print "Distance Restraints to elongate POM152 !!"
    print "dr_weight = ", dr_weight, "dist_min = ", dist_min, "dist_max = ", dist_max, "\n"

if (True):
    dist_min = 1.0
    dist_max = 35.0
    dr_weight = 100
    
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(444,444,"pom152"), (322,322,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_XL_444_322")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    print "Distance Restraints for XL 444-322 !!"
    print "dr_weight = ", dr_weight, "dist_min = ", dist_min, "dist_max = ", dist_max, "\n"


#####################################################
# Metropolis Monte Carlo sampling with Replica Exchange (PRE-SAMPLING)
#####################################################
sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "\nEVAL 1 : ", sf.evaluate(False), " (initial) - ", rank

if (False):
    simo.optimize_floppy_bodies(150)
    print "\nEVAL 2 : ", sf.evaluate(False), " (after calling optimize_floppy_bodies(150)) - ", rank

    initial_nframes = 5
    mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                        simo,
                                        monte_carlo_sample_objects = sampleobjects,
                                        output_objects = outputobjects,
                                        #crosslink_restraints=[xl1,xl2],
                                        #crosslink_restraints=[xl1],
                                        monte_carlo_temperature = 1.0,
                                        replica_exchange_minimum_temperature = 1.0,
                                        replica_exchange_maximum_temperature = 2.5,
                                        number_of_best_scoring_models = initial_nframes,
                                        #number_of_best_scoring_models = int(inputs.nrepeats),                                    
                                        #number_of_best_scoring_models = int(inputs.nrepeats)-10,
                                        monte_carlo_steps=10,
                                        number_of_frames = initial_nframes,
                                        #number_of_frames=50000,
                                        write_initial_rmf = True,
                                        initial_rmf_name_suffix = "initial",
                                        stat_file_name_suffix = "stat",
                                        best_pdb_name_suffix = "model",
                                        do_clean_first = True,
                                        do_create_directories = True,
                                        global_output_directory = "pre-EM2D_output",
                                        rmf_dir = "rmfs/",
                                        best_pdb_dir = "pdbs/",
                                        replica_stat_file_suffix = "stat_replica")
    mc1.execute_macro()
    rex1 = mc1.get_replica_exchange_object()
    print "\nEVAL 3 : ", sf.evaluate(False), " (after performing the pre-sampling) - ", rank
else:
    rex1 = None
    print "\n>> NO pre-sampling"
    

#####################################################
# Restraints setup
# EM 3D restraint using GMM
#####################################################
if (True):
    #### optimize a bit before adding the EM restraint
    simo.optimize_floppy_bodies(30000)
    ####

    # tail module em density
    mass = sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))
    print mass
    gem = IMP.pmi.restraints.em.GaussianEMRestraint(resdensities,
                                                    '../data/pom152_relion_s40.gmm.50.txt',
                                                    target_mass_scale=mass,
                                                    slope=0.000001,
                                                    target_radii_scale=3.0)
    gem.add_to_model()
    gem.set_weight(em_weight)
    #gem.center_model_on_target_density(simo)
    #gem.set_label("GaussianEMRestraint_Relion")
    outputobjects.append(gem)

    """
    # Electron Microscopy Restraint
    #  The GaussianEMRestraint uses a density overlap function to compare model to data
    #   First the EM map is approximated with a Gaussian Mixture Model (done separately)
    #   Second, the componets of the model are represented with Gaussians (forming the model GMM)
    #   Other options: scale_to_target_mass ensures the total mass of model and map are identical
    #                  slope: nudge model closer to map when far away
    #                  weight: experimental, needed becaues the EM restraint is quasi-bayesian
    em_components = bm1.get_density_hierarchies([t[1] for t in domains])
    print em_components
    gemt = IMP.pmi.restraints.em.GaussianEMRestraint(em_components,
                                                     '../data/run1_class001.gmm.50.txt',
                                                     scale_target_to_mass=True,
                                                     slope=0.000001,
                                                     weight=100.0)
    gemt.add_to_model()
    #gem.set_weight(100.0)
    gem.center_model_on_target_density(simo)
    #gem.set_label("GaussianEMRestraint_Relion")
    outputobjects.append(gemt)
    """

"""
#####################################################
# Restraints setup
# EM 2D restraint for each class
#####################################################
if (inputs.em2d_input is not None):
    if (True):
        images = [inputs.em2d_input]
    else:
        images = []
        for class_num in range(0, 5):
            pgm = "../data/em2d/pom152_wt_resized_classes." + str(class_num) + ".pgm"
            images.append(pgm)
    print images

    em2d = IMP.pmi.restraints.em2d.ElectronMicroscopy2D(simo,
                                                        images,
                                                        resolution = 1.0,
                                                        pixel_size = 2.03,
                                                        image_resolution = 10.0,
                                                        projection_number = 100)
                                                        #projection_number = 400)
    em2d.add_to_model()
    em2d.set_weight(em_weight)
    outputobjects.append(em2d)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
    print "\nEVAL 4 : ", sf.evaluate(False), " (after applying the EM 2D restraint) - ", rank
"""


#####################################################
# Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
# TODO: Ask how to save pdb files in the correct sequence order
mc2=IMP.pmi.macros.ReplicaExchange0(m,
                                    simo,
                                    monte_carlo_sample_objects = sampleobjects,
                                    output_objects = outputobjects,
                                    #crosslink_restraints=[xl1,xl2],
                                    #crosslink_restraints=[xl1],
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 2.5,
                                    number_of_best_scoring_models = 100,
                                    #number_of_best_scoring_models = int(inputs.nrepeats),                                    
                                    #number_of_best_scoring_models = int(inputs.nrepeats)-10,
                                    monte_carlo_steps=10,
                                    number_of_frames = int(inputs.nrepeats),
                                    #number_of_frames=50000,
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

