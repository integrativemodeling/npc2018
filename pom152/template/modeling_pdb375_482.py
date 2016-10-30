#!/usr/bin/env python
#####################################################
# Last Update: Oct 24, 2016
# by Seung Joong Kim
# at Andrej Sali group, University of California San Francisco (UCSF)
#####################################################
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
#import IMP.pmi.restraints.em2d
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.proteomics
import IMP.pmi.representation
#import representation_pom152
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.tools
import IMP.pmi.output
import IMP.pmi.samplers
#import IMP.pmi.topology
#import IMP.pmi.dof
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
parser.add_argument('-refine', action="store", dest="refinement", help="refinement True or False (XL distance restraints for dimers)" )
parser.add_argument('-weight', action="store", dest="weight", help="weight for XL distance restaints for dimer" )
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
    inputs.XL_input = "../data/XL.csv"
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
simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)
#simo = representation_pom152.Representation(m,upperharmonic=True,disorderedlength=False)
#simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=True)

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
# comp_name, hier_name, color,  fasta_file,                   fasta_id,    pdb_name,              chain_id, res_range, read_em_files, bead_size,   rb, super_rb, em_num_components, em_txt_file_name, em_mrc_file_name, chain_of_super_rb, keep_gaussian_on_flexible_beads
domains = \
[
# ("pom152",  "detgnt_99",  1.0,  fasta_files+"detgnt.txt",   "detgnt",  "BEADS",                    " ",  (det_ini,det_pos,0), True, beadsize_det,  99,   [1],   1,   None, None, None),

 ("pom152",  "pom152_11a", 0.0,  fasta_files+"pom152.txt",   "pom152",  "BEADS",                    " ",  (  1,100,0),    True,      beadsize20,    111,  [1],   1,   None, None, None),
 ("pom152",  "pom152_11b", 0.0,  fasta_files+"pom152.txt",   "pom152",  "BEADS",                    " ",  (101,200,0),    True,      beadsize100,   112,  [1],   1,   None, None, None),
 ("pom152",  "pom152_11c", 0.0,  fasta_files+"pom152.txt",   "pom152",  "BEADS",                    " ",  (201,378,0),    True,      beadsize20,    113,  [1],   1,   None, None, None),

 ("pom152",  "pom152_1",   0.1,  fasta_files+"pom152.txt",   "pom152",  pdb_files+"375_482.pdb",    "A",  (379,472,0),    True,      beadsize,      1,    [1],   4,   None, None, None),
 ("pom152",  "pom152_12",  0.1,  fasta_files+"pom152.txt",   "pom152",  "BEADS",                    " ",  (473,519,0),    True,      beadsize10,    1,    [1],   1,   None, None, None),

 ("pom152",  "pom152_2",   0.2,  fasta_files+"pom152.txt",   "pom152",  pdb_files+"516_611.pdb",    "A",  (520,611,0),    True,      beadsize,      2,    [1],   4,   None, None, None),
 ("pom152",  "pom152_22",  0.2,  fasta_files+"pom152.txt",   "pom152",  "BEADS",                    " ",  (612,615,0),    True,      beadsize10,    2,    [1],   1,   None, None, None),

 ("pom152",  "pom152_3",   0.3,  fasta_files+"pom152.txt",   "pom152",  pdb_files+"27005.pdb",      "A",  (616,714,0),    True,      beadsize,      3,    [1],   4,   None, None, None),
 ("pom152",  "pom152_32",  0.3,  fasta_files+"pom152.txt",   "pom152",  "BEADS",                    " ",  (715,721,0),    True,      beadsize10,    3,    [1],   1,   None, None, None),

 ("pom152",  "pom152_4",   0.4,  fasta_files+"pom152.txt",   "pom152",  pdb_files+"27005.pdb",      "A",  (722,818,0),    True,      beadsize,      4,    [1],   4,   None, None, None),
 ("pom152",  "pom152_42",  0.4,  fasta_files+"pom152.txt",   "pom152",  "BEADS",                    " ",  (819,823,0),    True,      beadsize10,    4,    [1],   1,   None, None, None),

 ("pom152",  "pom152_5",   0.5,  fasta_files+"pom152.txt",   "pom152",  pdb_files+"26996.pdb",      "A",  (824,918,0),    True,      beadsize,      5,    [1],   4,   None, None, None),
 ("pom152",  "pom152_52",  0.5,  fasta_files+"pom152.txt",   "pom152",  "BEADS",                    " ",  (919,930,0),    True,      beadsize10,    5,    [1],   1,   None, None, None),

 ("pom152",  "pom152_6",   0.6,  fasta_files+"pom152.txt",   "pom152",  pdb_files+"26996.pdb",      "A",  (931,1026,0),   True,      beadsize,      6,    [1],   4,   None, None, None),
 ("pom152",  "pom152_62",  0.6,  fasta_files+"pom152.txt",   "pom152",  "BEADS",                    " ",  (1027,1035,0),  True,      beadsize10,    6,    [1],   1,   None, None, None),

 ("pom152",  "pom152_7",   0.7,  fasta_files+"pom152.txt",   "pom152",  pdb_files+"26996.pdb",      "A",  (1036,1141,0),  True,      beadsize,      7,    [1],   4,   None, None, None),
 ("pom152",  "pom152_72",  0.7,  fasta_files+"pom152.txt",   "pom152",  "BEADS",                    " ",  (1142,1149,0),  True,      beadsize10,    7,    [1],   1,   None, None, None),

 ("pom152",  "pom152_8",   0.8,  fasta_files+"pom152.txt",   "pom152",  pdb_files+"1146_1237.pdb",  "A",  (1150,1229,0),  True,      beadsize,      8,    [1],   4,   None, None, None),
 ("pom152",  "pom152_82",  0.8,  fasta_files+"pom152.txt",   "pom152",  "BEADS",                    " ",  (1230,1243,0),  True,      beadsize10,    8,    [1],   1,   None, None, None),

 ("pom152",  "pom152_9",   0.9,  fasta_files+"pom152.txt",   "pom152",  pdb_files+"1238_1337.pdb",  "A",  (1244,1337,0),  True,      beadsize,      9,    [1],   4,   None, None, None),
]

bm1=IMP.pmi.macros.BuildModel1(simo)
#bm1.set_gmm_models_directory(datadirectory+"em_gmm_model/")

if (inputs.rmf_input is not None) :
    dom = set([s[0] for s in domains])
    for entry in list(dom):
        bm1.set_rmf_file(entry, inputs.rmf_input, int(inputs.rmf_frame_number))

bm1.build_model(data_structure = domains, sequence_connectivity_scale=sc_scale, sequence_connectivity_resolution=1.0)
#bm1.scale_bead_radii(40,0.8)

#####################################################
# randomize the initial configuration
#####################################################
if (inputs.rmf_input is None) :
    #simo.shuffle_configuration(50)
    simo.shuffle_configuration(250)


#####################################################
# defines the movers
#####################################################
# Add default mover parameters to simulation
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
if (True):
    ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo, resolution = res_ev)
    ev.add_to_model()
    outputobjects.append(ev)
    print(ev.get_output())
    print "ExcludedVolumeSphere !!\n"


#####################################################
# Restraints setup
# External Barrier restraint
#####################################################
if (False):
    eb = IMP.pmi.restraints.basic.ExternalBarrier(simo, radius = 1000)
    eb.add_to_model()
    outputobjects.append(eb)
    print(eb.get_output())
    print "ExternalBarrier !!\n"


#####################################################
# Restraints setup
# End-to-End Distance restraint
#####################################################
if (True):
    dist_min = 350.0
    dist_max = 400.0
    dr_weight = 100.0

    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(379,379,"pom152"), (1337,1337,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_EndToEnd")
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    print "End-to-End Distance Restraints to elongate pom152 !!"
    print "dr_weight = ", dr_weight, "dist_min = ", dist_min, "dist_max = ", dist_max, "\n"


#####################################################
# Restraints setup
# Distance restraints between neighboring Ig domains
# Domain connectivity
#####################################################
if (True):
    dist_min = 7.5
    dr_weight = 100.0

    IG_LIST = [ [472,520], [611,616], [714,722], [818,824], [918,931], [1026,1036], [1141,1150], [1229,1244] ]
    for z in IG_LIST:
        if (z[0] == 472):
            dist_max = 17.5
        else:
            dist_max = 12.5
            
        dr = IMP.pmi.restraints.basic.DistanceRestraint(simo, (z[0],z[0],"pom152"), (z[1],z[1],"pom152"), distancemin=dist_min, distancemax=dist_max, resolution=res_str)
        dr.add_to_model()
        dr.set_label('pom152_%d_%d' % (z[0], z[1]))
        dr.set_weight(dr_weight)
        outputobjects.append(dr)
        print(dr.get_output())


#####################################################
# Restraints setup
# Cross-link restraints using the whole NPC DSS XL data
#####################################################
if (False):
    sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
    print "\nEVAL 0 : ", sf.evaluate(False), " (before applying the XL restraint) - ", rank

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
                                                        '../data/XL.csv',
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
#### optimize a bit before adding the EM restraint
#####################################################
if (inputs.rmf_input is None) :
    simo.optimize_floppy_bodies(30000)


#####################################################
# Restraints setup
# Sampling Boundary Restraint for the inner ring
#####################################################
if (True):
    resdensities = bm1.get_density_hierarchies([t[1] for t in domains])
    print ("resdensities=", resdensities)       ####  make sure resdensities are correct

    mass = sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))
    print ("Total mass for the Sampling Boundary EM restraint = ", mass)
    sbr = IMP.pmi.restraints.em.GaussianEMRestraint(resdensities,
                                                    '../data/pom152_relion_s40.gmm.50.txt',
                                                    target_mass_scale=mass,
                                                    slope=0.01,
                                                    #slope=0.0000001,
                                                    target_radii_scale=3.0)
    sbr.add_to_model()
    sbr.set_weight(0.5)        # play with the weight
    sbr.set_label("Sampling_Boundary")
    #sbr.center_model_on_target_density(simo)
    outputobjects.append(sbr)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
    print "\nEVAL 2 : ", sf.evaluate(False), " (after applying the Sampling Boundary EM restraint) - ", rank


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
                p = IMP.pmi.tools.select_by_tuple(representation, selection_tuple, resolution=1)
                #print(IMP.atom.Residue.get_is_setup(p[0]))
                ps1.append(p[0])
            print (ps1)

            ps2 = []
            for selection_tuple in t2:
                #print (selection_tuple)
                p = IMP.pmi.tools.select_by_tuple(representation, selection_tuple, resolution=1)
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
        IMP.pmi.tools.add_restraint_to_model(self.m, self.rs)

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
            ps1.append((z,z,"pom152"))
        ps2=[]
        for z in t2:
            ps2.append((z,z,"pom152"))
        print (ps1, ps2)

        pdr = PlaneDihedralRestraint(simo, [ps1, ps2], angle=90.0, k=1.)
        pdr.add_to_model()
        pdr.set_label('pom152_pdr_%d_%d' % (t1[0], t2[0]))
        pdr.set_weight(pdr_weight)
        outputobjects.append(pdr)
        print(pdr.get_output())


#####################################################
# Metropolis Monte Carlo sampling with Replica Exchange (PRE-SAMPLING)
#####################################################
sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
print "\nEVAL 4 : ", sf.evaluate(False), " (initial) - ", rank

if (False):
    initial_nframes = 1000
    mc1=IMP.pmi.macros.ReplicaExchange0(m,
                                        simo,
                                        monte_carlo_sample_objects = sampleobjects,
                                        output_objects = outputobjects,
                                        #crosslink_restraints=[xl1,xl2],
                                        crosslink_restraints = XL_restraints,
                                        monte_carlo_temperature = 1.0,
                                        replica_exchange_minimum_temperature = 1.0,
                                        replica_exchange_maximum_temperature = 2.5,
                                        number_of_best_scoring_models = 5,
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
                                        global_output_directory = "pre-EM_output",
                                        rmf_dir = "rmfs/",
                                        best_pdb_dir = "pdbs/",
                                        replica_stat_file_suffix = "stat_replica")
    mc1.execute_macro()
    rex1 = mc1.get_replica_exchange_object()
    print "\nEVAL 5 : ", sf.evaluate(False), " (after performing the pre-sampling) - ", rank
else:
    rex1 = None
    print "\n>> NO pre-sampling"


#####################################################
# Restraints setup
# EM 3D restraint using GMM  (Electron Microscopy Restraint)
#  The GaussianEMRestraint uses a density overlap function to compare model to data
#   First the EM map is approximated with a Gaussian Mixture Model (done separately)
#   Second, the componets of the model are represented with Gaussians (forming the model GMM)
#   Other options: scale_to_target_mass ensures the total mass of model and map are identical
#                  slope: nudge model closer to map when far away
#                  weight: experimental, needed becaues the EM restraint is quasi-bayesian
#####################################################
if (True):
    gem = IMP.pmi.restraints.em.GaussianEMRestraint(resdensities,
                                                    '../data/pom152_relion_s40.gmm.50.txt',
                                                    target_mass_scale=mass,
                                                    slope=0.0000001,
                                                    target_radii_scale=3.0)
    gem.add_to_model()
    gem.set_weight(5000.0)
    gem.set_label("GaussianEMRestraint_Relion")
    #gem.center_model_on_target_density(simo)
    outputobjects.append(gem)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi.tools.get_restraint_set(m))
    print "\nEVAL 6 : ", sf.evaluate(False), " (after applying the 3D EM restraint) - ", rank


#####################################################
# Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
# TODO: Ask how to save pdb files in the correct sequence order
mc2=IMP.pmi.macros.ReplicaExchange0(m,
                                    simo,
                                    monte_carlo_sample_objects = sampleobjects,
                                    output_objects = outputobjects,
                                    #crosslink_restraints=[xl1,xl2],
                                    crosslink_restraints = XL_restraints,
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 2.5,
                                    number_of_best_scoring_models = 5,
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
print "\nEVAL 7 : ", sf.evaluate(False), " (final evaluation) - ", rank


