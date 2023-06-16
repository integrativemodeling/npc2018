#!/usr/bin/env python
#####################################################
# Last Update: June 30th, 2016 by Seung Joong Kim
# Andrej Sali group, University of California San Francisco (UCSF)
#####################################################
from __future__ import print_function
import collections
import RMF
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container

import ihm.location
import ihm.dataset
try:
    import ihm.reference
except ImportError:
    pass
try:
    import ihm.citations
except ImportError:
    pass
import IMP.pmi1.mmcif
import IMP.pmi1.restraints.crosslinking
import IMP.pmi1.restraints.stereochemistry
import IMP.pmi1.restraints.em
import IMP.pmi1.restraints.em2d
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
import sys
import re
import math
import glob

sys.path.append('../util/')
import make_archive

#####################################################
# Parsing parameter inputs
#####################################################
import argparse

parser = argparse.ArgumentParser(description='Performing the INITIAL/REFINEMENT Monte Carlo job, with crosslinks and selected/ALL domain mapping data. Example of usage: setup_environment.sh python ./sj_SEA_XLDM.py -f models_1877.rmf -n 0')
parser.add_argument('--test', action='store_true',
                    help="Run in test mode (fewer steps of sampling)")
parser.add_argument('--dry-run', action='store_true',
                    help="Dry run (do not do any sampling)")
parser.add_argument('--one-spoke', action='store_true',
                    help="Speed up sampling by modeling only one spoke (not 3)")
parser.add_argument('--mmcif', action='store', type=str, default=None,
                    help="Record modeling protocol in a named mmCIF file")
parser.add_argument('--no-symmetry', action='store_false', dest="symmetry",
                    help="Don't add symmetry to make all 8 spokes")
parser.add_argument('-copy', action="store", dest="ncopy", help="copy numbers (stoichiometry) for SEA4 and Seh1" )
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
print(inputs)


#####################################################
# setting up topology and parameters
#####################################################
m = IMP.Model()
#s = IMP.pmi1.topology.System(m)
#st = s.create_state()
simo = IMP.pmi1.representation.Representation(m,upperharmonic=True,disorderedlength=False)

# Record additional software used
# Detection of templates for comparative modeling (HHPred)
s = ihm.Software(
          name='HHpred', classification='protein homology detection',
          description='Protein homology detection by HMM-HMM comparison',
          version='2.0.16',
          location='https://toolkit.tuebingen.mpg.de/hhpred')
if hasattr(ihm, 'citations'):
    s.citation = ihm.citations.hhpred
simo.add_metadata(s)

# Secondary structure prediction (PSIPRED)
s = ihm.Software(
          name='PSIPRED', classification='secondary structure prediction',
          description='Protein secondary structure prediction based on '
                      'position-specific scoring matrices',
          version='4.0',
          location='http://bioinf.cs.ucl.ac.uk/psipred/')
if hasattr(ihm, 'citations'):
    s.citation = ihm.citations.psipred
simo.add_metadata(s)

# Disordered region prediction (DISOPRED)
s = ihm.Software(
          name='DISOPRED', classification='disorder prediction',
          description='prediction of protein disorder', version=3,
          location='http://bioinf.cs.ucl.ac.uk/psipred/?disopred=1')
if hasattr(ihm, 'citations'):
    s.citation = ihm.citations.disopred
simo.add_metadata(s)

# Domain boundary prediction (DomPred)
s = ihm.Software(
          name='DomPred', classification='domain boundary prediction',
          description='identification of putative domain boundaries',
          location='http://bioinf.cs.ucl.ac.uk/dompred')
if hasattr(ihm, 'citations'):
    s.citation = ihm.Citation(
        pmid='17430199',
        title='Computer-assisted protein domain boundary prediction '
              'using the DomPred server.',
        journal='Curr Protein Pept Sci', volume=8, page_range=('181', '188'),
        year=2007, authors=['Bryson K', 'Cozzetto D', 'Jones DT'],
        doi='10.2174/138920307780363415')
simo.add_metadata(s)

# Coiled-coil region prediction (COILS/PCOILS)
s = ihm.Software(
          name='COILS/PCOILS', classification='coiled-coil prediction',
          description='prediction of coiled-coil structure',
          location='https://toolkit.tuebingen.mpg.de/#/tools/pcoils')
if hasattr(ihm, 'citations'):
    s.citation = ihm.Citation(
        pmid='2031185',
        title='Predicting coiled coils from protein sequences.',
        journal='Science', volume=252, page_range=('1162', '1164'),
        year=1991, authors=['Lupas A', 'Van Dyke M', 'Stock J'],
        doi='10.1126/science.252.5009.1162')
simo.add_metadata(s)

# EM particle picking (EMAN2) & generation of Nic96 class averages (ISAC)
s = ihm.Software(
          name='EMAN2', classification='image processing',
          description='processing of data from transmission electron '
                      'microscopes',
          version='2.2',
          location='http://blake.bcm.edu/emanwiki/EMAN2')
if hasattr(ihm, 'citations'):
    s.citation = ihm.Citation(
        pmid='16859925',
        title='EMAN2: an extensible image processing suite for electron '
              'microscopy.', journal='J Struct Biol', volume=157,
        page_range=('38', '46'), year=2007,
        authors=['Tang G', 'Peng L', 'Baldwin PR', 'Mann DS', 'Jiang W',
                 'Rees I', 'Ludtke SJ'], doi='10.1016/j.jsb.2006.05.009')
simo.add_metadata(s)

# Parallel tomographic analysis (Relion)
s = ihm.Software(
          name='RELION', classification='image processing',
          description='refinement of (multiple) 3D reconstructions or '
                      '2D class averages in electron cryo-microscopy',
          version='1.4',
          location='https://www2.mrc-lmb.cam.ac.uk/relion/')
if hasattr(ihm, 'citations'):
    s.citation = ihm.citations.relion
simo.add_metadata(s)

# Prediction of transmembrane domains (SGD)
s = ihm.Software(
          name='SGD', classification='database',
          description='biological information for the budding yeast '
                      'Saccharomyces cerevisiae along with search and '
                      'analysis tools to explore these data',
          location='https://www.yeastgenome.org/')
if hasattr(ihm, 'citations'):
    s.citation = ihm.Citation(
        pmid='22110037',
        title='Saccharomyces Genome Database: the genomics resource of '
              'budding yeast.', journal='Nucleic Acids Res', volume=40,
        page_range=('D700', 'D705'), year=2012,
        authors=['Cherry JM', 'Hong EL', 'Amundsen C', 'Balakrishnan R',
                 'Binkley G', 'Chan ET', 'Christie KR', 'Costanzo MC',
                 'Dwight SS', 'Engel SR', 'Fisk DG', 'Hirschman JE',
                 'Hitz BC', 'Karra K', 'Krieger CJ', 'Miyasato SR',
                 'Nash RS', 'Park J', 'Skrzypek MS', 'Simison M', 'Weng S',
                 'Wong ED'], doi='10.1093/nar/gkr1029')
simo.add_metadata(s)

# Prediction of membrane binding motifs (HeliQuest)
s = ihm.Software(
          name='HeliQuest', classification='helix prediction',
          description='prediction of helix content from primary sequence',
          location='http://heliquest.ipmc.cnrs.fr/')
if hasattr(ihm, 'citations'):
    s.citation = ihm.Citation(
        pmid='18662927',
        title='HELIQUEST: a web server to screen sequences with specific '
              'alpha-helical properties.', journal='Bioinformatics', volume=24,
        page_range=('2101', '2102'), year=2008,
        authors=['Gautier R', 'Douguet D', 'Antonny B', 'Drin G'],
        doi='10.1093/bioinformatics/btn392')
simo.add_metadata(s)

simo.dry_run = inputs.dry_run

try:
    from mpi4py import MPI
except ImportError:
    MPI = None

if MPI:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
else:
    rank = 0
print("rank = ", rank)

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

use_neighboring_spokes = not inputs.one_spoke
#Stopwatch_None_delta_seconds  ~22   (1 spoke for OR / IR + 3 spokes for others, 3.0G memory) with XL
#Stopwatch_None_delta_seconds  ~25   (1 spoke for OR / IR + 3 spokes for others, 3.0G memory) with XL + EM
#Stopwatch_None_delta_seconds  ~30   (1 spoke for OR / IR + 3 spokes for others, 3.0G memory) with XL + EM + EV
use_shuffle = True
use_ExcludedVolume = True
use_Immuno_EM = False
use_Composite = False
use_Distance_to_Point = False
use_end_to_end_157_170 = False
use_FG_anchor = False
use_sampling_boundary = False
use_XL = True
use_EM3D = True

class _ChainOfFGs(object):
    """Override IMP.pmi1.mmcif._Chain to substitute coordinates for FGs"""
    def __init__(self, orig_chain, fg_repeats, fg_spheres):
        self.orig_chain = orig_chain
        self.fg_repeats = fg_repeats
        self.fg_spheres = fg_spheres
    def _get_fg_spheres(self, fg_anchor, npc_anchor, residue_range, fg_spheres):
        # Make sure FG anchor matches NPC scaffold anchor, if the NPC is
        # unstructured (cannot compare to resolution=1 bead)
        # npc_anchor[3] is the list of residue indexes, which is None for
        # resolution=1 beads
        xyz0 = IMP.algebra.Vector3D(fg_anchor)
        xyz1 = IMP.algebra.Vector3D(npc_anchor[0])
        if npc_anchor[3] is not None and (xyz0 - xyz1).get_magnitude() > 0.1:
            raise ValueError("FG anchor coordinates mismatch; %s in FG "
                             "output RMF vs %s in NPC scaffold output RMF "
                             "for %s" % (xyz0, xyz1, self.comp))
        resstart = residue_range[0]
        for xyz, radius in fg_spheres:
            resend = min(resstart + self.fg_repeats.beadsize,
                          residue_range[1] + 1)
            yield (xyz, "BEA", resstart + resend // 2,
                   list(range(resstart, resend)), radius)
            resstart = resend
    def __get_spheres(self):
        fg_spheres = self.fg_spheres.get(self.comp, None)
        if fg_spheres:
            nup = self.comp.split('.')[0].split('@')[0]
            rng = self.fg_repeats.ranges[nup]
            orig_spheres = list(self.orig_chain.spheres)
            if rng[1] < rng[0]:
                for f in self._get_fg_spheres(fg_spheres[0][0],
                                              orig_spheres[0],
                                              (rng[1], rng[0]),
                                              reversed(fg_spheres[1:])):
                    yield f
            else:
                for f in self._get_fg_spheres(fg_spheres[0][0],
                                              orig_spheres[-1],
                                              rng, fg_spheres[1:]):
                    yield f
    spheres = property(__get_spheres)

    chain_id = property(lambda self: self.orig_chain.chain_id)
    # No atomic coordinates for FGs
    atoms = []
    asym_unit = property(lambda self: self.orig_chain.asym_unit)
    comp = property(lambda self: self.orig_chain.comp)
    orig_comp = property(lambda self: self.orig_chain.orig_comp)


class FGRepeats(object):
    """Handle FG repeat regions"""
    beadsize = 20
    Bead = collections.namedtuple('Bead', ['start', 'end', 'len'])
    anchor = {'Nsp1': (601,636),
              'Nup1': (301,350),
              'Nup49': (201,269),
              'Nup57': (201,286),
              'Nup60': (351,398),
              'Nup100': (551,575),
              'Nup116': (751,775),
              'Nup145': (201,225),
              'Nup159': (1082,1116)}
    ranges = {'Nsp1': (600,1),
              'Nup1': (352,1049),
              'Nup49': (200,1),
              'Nup57': (200,1),
              'Nup60': (399,498),
              'Nup100': (550,11),
              'Nup116': (750,11),
              'Nup145': (200,1),
              'Nup159': (1081,482)}
    fasta_fn = {'Nsp1': f_n82+"Nsp1.txt",
                'Nup1': f_npc+"Nup1.txt",
                'Nup49': f_n96+"Nup49.txt",
                'Nup57': f_n96+"Nup57.txt",
                'Nup60': f_npc+"Nup60.txt",
                'Nup100': f_npc+"Nup100.txt",
                'Nup116': f_n82+"Nup116.txt",
                'Nup145': f_npc+"Nup145.txt",
                'Nup159': f_n82+"Nup159.txt"}
    fasta_id = {'Nup100': 'YKL068W',
                'Nup145': 'YGL092W',
                'Nup49': 'YGL172W',
                'Nup60': 'YAR002W',
                'Nup57': 'YGR119C',
                'Nup1': 'YOR098C'}
    copies = {'Nsp1': ('.1', '.2', '.3', '.4', '.3@11', '.4@11'),
              'Nup1': ('',),
              'Nup49': ('.1', '.2', '.1@11', '.2@11'),
              'Nup57': ('.1', '.2', '.1@11', '.2@11'),
              'Nup60': ('.1', '.2'),
              'Nup100': ('.1', '.2'),
              'Nup116': ('.1', '.2'),
              'Nup145': ('.1', '.2'),
              'Nup159': ('.1', '.2')}

    def __init__(self, po):
        self.representation = po.create_representation("FG repeats")

    def get_assembly(self):
        """Get all components that constitute FG repeats"""
        compdict = {}
        for nup, rng in self.ranges.items():
            r = tuple(sorted(rng))
            for copy in self.get_copy_names(nup):
                compdict[copy] = r
        return compdict

    def get_number_of_beads(self, name):
        rng = self.ranges[name]
        return (max(rng) - min(rng) + self.beadsize) // self.beadsize

    def get_beads(self, name, start, end):
        """If an anchor point is passed, return the FG repeat beads"""
        # Remove copy number
        comp = name.split('.')[0].split('@')[0]
        if inputs.symmetry and self.anchor.get(comp, None) == (start, end):
            rng = self.ranges[comp]
            return self.Bead(start=min(rng), end=max(rng),
                             len=self.get_number_of_beads(comp))

    def create_assembly(self, po):
        self.assembly = po._get_subassembly(
                              self.get_assembly(), name="FG repeats",
                              description="All FG repeats (unstructured regions"
                                          " in the center of the pore)")

    def add_bead_coordinates(self, fname, model):
        """Read in FG bead coordinates.
           The model already contains coordinates for the scaffold. Replace
           these with coordinates for the FG repeats, checking that the
           anchor regions match."""
        if not inputs.symmetry:
            return
        spheres = self.read_rmf(fname)

        # Monkey patch the Model.all_chains() method to use FG coordinates
        orig_all_chains = model.all_chains
        def patched_all_chains(slf, simo):
            for c in orig_all_chains(simo):
                yield _ChainOfFGs(c, self, spheres)
        model.all_chains = patched_all_chains.__get__(model, type(model))

    def read_rmf(self, fname):
        r = RMF.open_rmf_file_read_only(fname)
        r.set_current_frame(RMF.FrameID(0))
        node_for_nup = {}
        for node in r.get_root_node().get_children():
            if node.get_name() in self.anchor:
                node_for_nup[node.get_name()] = node
        assert(len(node_for_nup) == len(self.anchor))

        spheres = {}
        for nup in self.anchor:
            self._read_nup(nup, r, node_for_nup, spheres)
        return spheres

    def get_copy_names(self, name):
        """Yield all names of Nup copies, including symmetry units.
           These names should match the ordering of FG repeats in the RMF
           file."""
        for c in self.copies[name]:
            if c.endswith('@11'):
                # NPC symmetry copies alternate positive/negative rotations,
                # while Nups are always positive rotations
                for symm in ('@11', '@12', '@14', '@16', '@18', '@17',
                             '@15', '@13'):
                    yield(name + c[:-3] + symm)
            else:
                for symm in ('', '@2', '@4', '@6', '@8', '@7', '@5', '@3'):
                    yield(name + c + symm)

    def _read_nup(self, name, rh, nodemap, spheres):
        rff = RMF.ReferenceFrameConstFactory(rh)
        pf = RMF.ParticleConstFactory(rh)
        node = nodemap[name]
        node_copies = node.get_children()
        assert(len(node_copies) == len(self.copies[name]) * 8)
        for nc, copyname in zip(node_copies, self.get_copy_names(name)):
            node_beads = nc.get_children()
            # First bead is anchor
            assert(len(node_beads) == self.get_number_of_beads(name) + 1)
            spheres[copyname] = s = []
            for b in node_beads:
                assert(pf.get_is(b))
                assert(rff.get_is(b))
                radius = pf.get(b).get_radius()
                coord = rff.get(b).get_translation()
                s.append((coord, radius))


class ProtocolOutput(IMP.pmi1.mmcif.ProtocolOutput):
    def __init__(self, *args, **keys):
        super(ProtocolOutput, self).__init__(*args, **keys)
        self.fgs = FGRepeats(self)

    def add_bead_element(self, state, name, start, end, num, hier):
        super(ProtocolOutput, self).add_bead_element(state, name, start,
                                                     end, num, hier)
        # Add representation for FG repeat beads
        beads = self.fgs.get_beads(name, start, end)
        if beads:
            super(ProtocolOutput, self).add_bead_element(state, name,
                                      beads.start, beads.end, beads.len, hier,
                                      self.fgs.representation)


if inputs.mmcif:
    # Record the modeling protocol to an mmCIF file
    po = ProtocolOutput(open(inputs.mmcif, 'w'))
    simo.add_protocol_output(po)
    # Exclude domains that weren't well resolved in the modeling from mmCIF
    po.system.comments.append(
                   "Missing residues: the following residues were poorly "
                   "constrained by the available experimental data, and as "
                   "such coordinates are not available in this model: "
                   "Mlp1 717-1875; Mlp2 691-1670; Nup42 1-430; Gle1 121-538.")
    if inputs.symmetry:
        subtitle = 'eight spokes'
    elif inputs.one_spoke:
        subtitle = 'a single spoke'
    else:
        subtitle = 'three spokes'
    po.system.title = ('Integrative structure and functional anatomy of '
                       '%s of a nuclear pore complex' % subtitle)

    # Add publication
    po.system.citations.append(ihm.Citation.from_pubmed_id(29539637))

    if use_neighboring_spokes:
        suffixes = ['', '@2', '@3']
    else:
        suffixes = ['']
    for suffix in suffixes:
        po.exclude_coordinates('Mlp1'+suffix, (717,1875))
        po.exclude_coordinates('Mlp2'+suffix, (691,1679))
        po.exclude_coordinates('Nup42'+suffix, (364,430))
        po.exclude_coordinates('Gle1'+suffix, (121,538))

    # Point to repositories where files are deposited
    simo.add_metadata(ihm.location.Repository(
        doi="10.5281/zenodo.1194533", root="npc_fg_2018",
        url="https://zenodo.org/record/1194533/files/npc_fg_2018-master.zip",
        top_directory="npc_fg_2018-master"))
    for subdir, zipname in make_archive.ARCHIVES.items():
        simo.add_metadata(ihm.location.Repository(
              doi="10.5281/zenodo.1194547", root="../%s" % subdir,
              url="https://zenodo.org/record/1194547/files/%s.zip" % zipname,
              top_directory=os.path.basename(subdir)))
        simo.add_metadata(ihm.location.Repository(
              doi="10.5281/zenodo.1194547", root="..",
              url="https://zenodo.org/record/1194547/files/npc2018-master.zip",
              top_directory="npc2018-master"))


#####################################################
# REPRESENTATION
#####################################################
# comp_name, hier_name, color, fasta_file, fasta_id, pdb_name, chain_id, res_range, read_em_files, bead_size, rb, super_rb, em_num_components, em_txt_file_name, em_mrc_file_name, chain_of_super_rb, keep_gaussian_on_flexible_beads
domains = []
if (use_EM3D):  gmm = True
else:           gmm = None
if (use_neighboring_spokes):
    clones_range_A = list(range(2,4)) + list(range(11,14))
    clones_range_B = list(range(2,4))
else:
    clones_range_A = list(range(11,12))
    clones_range_B = []
##########################
# Nup84 complex
##########################
if (is_n84):
    #n84_rb = 84;    n133_rb = 133;  n120_rb = 120
    n84_rb = None;    n133_rb = None;  n120_rb = None
    domains.append(('Nup84',  "Nup84",    0.0,  n84_fastafile,   "Nup84",   n84_pdb, "K",  (  1, 726,0),  gmm,  beadsize,  n84_rb,  None,    3,  " ",   " ",  None, False))
    domains.append(('Nup85',  "Nup85_1",  0.2,  n85_fastafile,   "Nup85",   n84_pdb, "L",  (  1, 492,0),  gmm,  beadsize,  n84_rb,  None,    3,  " ",   " ",  None, False))
    domains.append(('Nup85',  "Nup85_2",  0.25, n85_fastafile,   "Nup85",   n84_pdb, "L",  (493, 744,0),  gmm,  beadsize,  n84_rb,  None,    2,  " ",   " ",  None, False))
    domains.append(('Nup120', "Nup120_1", 0.35, n120_fastafile,  "Nup120",  n84_pdb, "M",  (  1, 714,0),  gmm,  beadsize,  n120_rb, None,    3,  " ",   " ",  None, False))
    domains.append(('Nup120', "Nup120_2", 0.4,  n120_fastafile,  "Nup120",  n84_pdb, "M",  (715,1037,0),  gmm,  beadsize,  n84_rb,  None,    2,  " ",   " ",  None, False))
    domains.append(('Nup133', "Nup133_1", 0.5,  n133_fastafile,  "Nup133",  n84_pdb, "N",  (  1, 489,0),  gmm,  beadsize,  n133_rb, None,    2,  " ",   " ",  None, False))
    domains.append(('Nup133', "Nup133_2", 0.55, n133_fastafile,  "Nup133",  n84_pdb, "N",  (490,1157,0),  gmm,  beadsize,  n84_rb,  None,    3,  " ",   " ",  None, False))
    domains.append(('Nup145c',"Nup145c_1",0.65, n145c_fastafile, "Nup145c", n84_pdb, "O",  (  1, 125,0),  gmm,  beadsize,  n84_rb,  None,    1,  " ",   " ",  None, False))
    domains.append(('Nup145c',"Nup145c_2",0.7,  n145c_fastafile, "Nup145c", n84_pdb, "O",  (126, 712,0),  gmm,  beadsize,  n84_rb,  None,    3,  " ",   " ",  None, False))
    domains.append(('Seh1',   "Seh1",     0.8,  seh1_fastafile,  "Seh1",    n84_pdb, "P",  (  1, 349,0),  gmm,  beadsize,  n84_rb,  None,    2,  " ",   " ",  None, False))
    domains.append(('Sec13',  "Sec13",    0.95, sec13_fastafile, "Sec13",   n84_pdb, "Q",  (  1, 297,0),  gmm,  beadsize,  n84_rb,  None,    2,  " ",   " ",  None, False))
    for i in clones_range_A:
    #for i in range(11,12):
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
    #n82_rb = 82
    n82_rb = None
    domains.append(("Dyn2.1",  "Dyn2.1",      0.48,  f_n82+"Dyn2.txt",   "Dyn2",   n82_pdb,  "A", (   1,  92,0),  gmm,   beadsize,   n82_rb,  None,    1,  " ",   " ",  None, False))
    domains.append(("Dyn2.2",  "Dyn2.2",      0.65,  f_n82+"Dyn2.txt",   "Dyn2",   n82_pdb,  "B", (   1,  92,0),  gmm,   beadsize,   n82_rb,  None,    1,  " ",   " ",  None, False))
    domains.append(("Nup82.1", "Nup82.1_1",   0.0,   f_n82+"Nup82.txt",  "Nup82",  n82_pdb,  "C", (   1, 452,0),  gmm,   beadsize,   n82_rb,  None,    2,  " ",   " ",  None, False))
    domains.append(("Nup82.1", "Nup82.1_2",   0.0,   f_n82+"Nup82.txt",  "Nup82",  n82_pdb,  "C", ( 453, 713,0),  gmm,   beadsize,   n82_rb,  None,    4,  " ",   " ",  None, False))
    domains.append(("Nup82.2", "Nup82.2_1",   0.15,  f_n82+"Nup82.txt",  "Nup82",  n82_pdb,  "D", (   1, 452,0),  gmm,   beadsize,   n82_rb,  None,    2,  " ",   " ",  None, False))
    domains.append(("Nup82.2", "Nup82.2_2",   0.15,  f_n82+"Nup82.txt",  "Nup82",  n82_pdb,  "D", ( 453, 713,0),  gmm,   beadsize,   n82_rb,  None,    4,  " ",   " ",  None, False))
    if (is_FG):
        domains.append(("Nup159.1","Nup159.1_1",  1.0,   f_n82+"Nup159.txt", "Nup159", n82_pdb,  "E", (   1, 381,0),  None,  beadsize,     1159,  None,    0,  None,  None, None))
        domains.append(("Nup159.1","Nup159.1_10", 1.0,   f_n82+"Nup159.txt", "Nup159", "BEADS",  " ", ( 382,1116,0),  None,  beadsize100,  1159,  None,    0,  None,  None, None))
        domains.append(("Nup159.2","Nup159.2_1",  0.9,   f_n82+"Nup159.txt", "Nup159", n82_pdb,  "F", (   1, 381,0),  None,  beadsize,     2159,  None,    0,  None,  None, None))
        domains.append(("Nup159.2","Nup159.2_10", 0.9,   f_n82+"Nup159.txt", "Nup159", "BEADS",  " ", ( 382,1116,0),  None,  beadsize100,  2159,  None,    0,  None,  None, None))
        domains.append(("Nsp1.1",  "Nsp1.1_10",   0.3,   f_n82+"Nsp1.txt",   "Nsp1",   "BEADS",  " ", (   1, 636,0),  None,  beadsize100,  1062,  None,    0,  None,  None, None))
        domains.append(("Nsp1.2",  "Nsp1.2_10",   0.38,  f_n82+"Nsp1.txt",   "Nsp1",   "BEADS",  " ", (   1, 636,0),  None,  beadsize100,  2062,  None,    0,  None,  None, None))
        #domains.append(("Nup116.1","Nup116.1_10", 0.75,  f_n82+"Nup116.txt", "Nup116", "BEADS",  " ", (   1, 750,0),  None,  beadsize100,  1116,  None,    0,  None,  None, None))
        #domains.append(("Nup116.2","Nup116.2_10", 0.8,   f_n82+"Nup116.txt", "Nup116", "BEADS",  " ", (   1, 750,0),  None,  beadsize100,  2116,  None,    0,  None,  None, None))
    else:
        domains.append(("Nup159.1","Nup159.1_10", 1.0,   f_n82+"Nup159.txt", "Nup159", "BEADS",  " ", (1082,1116,0),  None,  beadsize100, n82_rb, None,    0,  None,  None, None))
        domains.append(("Nup159.2","Nup159.2_10", 0.9,   f_n82+"Nup159.txt", "Nup159", "BEADS",  " ", (1082,1116,0),  None,  beadsize100, n82_rb, None,    0,  None,  None, None))
        domains.append(("Nsp1.1",  "Nsp1.1_10",   0.3,   f_n82+"Nsp1.txt",   "Nsp1",   "BEADS",  " ", ( 601, 636,0),  None,  beadsize100, n82_rb, None,    0,  None,  None, None))
        domains.append(("Nsp1.2",  "Nsp1.2_10",   0.38,  f_n82+"Nsp1.txt",   "Nsp1",   "BEADS",  " ", ( 601, 636,0),  None,  beadsize100, n82_rb, None,    0,  None,  None, None))
    domains.append(("Nup159.1","Nup159.1",    1.0,   f_n82+"Nup159.txt", "Nup159", n82_pdb,  "E", (1117,1460,0),  gmm,   beadsize,   n82_rb,  None,    6,  " ",   " ",  None, False))
    domains.append(("Nup159.2","Nup159.2",    0.9,   f_n82+"Nup159.txt", "Nup159", n82_pdb,  "F", (1117,1460,0),  gmm,   beadsize,   n82_rb,  None,    6,  " ",   " ",  None, False))
    domains.append(("Nsp1.1",  "Nsp1.1",      0.3,   f_n82+"Nsp1.txt",   "Nsp1",   n82_pdb,  "G", ( 637, 823,0),  gmm,   beadsize,   n82_rb,  None,    4,  " ",   " ",  None, False))
    domains.append(("Nsp1.2",  "Nsp1.2",      0.38,  f_n82+"Nsp1.txt",   "Nsp1",   n82_pdb,  "H", ( 637, 823,0),  gmm,   beadsize,   n82_rb,  None,    4,  " ",   " ",  None, False))
    #domains.append(("Nup116.1","Nup116.1",    0.75,  f_n82+"Nup116.txt", "Nup116", n82_pdb,  "I", ( 751,1113,0),  gmm,   beadsize25, n82_rb,  None,    1,  " ",   " ",  None, False))
    #domains.append(("Nup116.2","Nup116.2",    0.8,   f_n82+"Nup116.txt", "Nup116", n82_pdb,  "J", ( 751,1113,0),  gmm,   beadsize25, n82_rb,  None,    1,  " ",   " ",  None, False))
    for i in clones_range_B:
    #for i in []:
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
            #domains.append(("Nup116.1@%d"%i,"Nup116.1_10@%d"%i, 0.75,  f_n82+"Nup116.txt", "Nup116", "BEADS",  " ", (   1, 750,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
            #domains.append(("Nup116.2@%d"%i,"Nup116.2_10@%d"%i, 0.8,   f_n82+"Nup116.txt", "Nup116", "BEADS",  " ", (   1, 750,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
        else:
            domains.append(("Nup159.1@%d"%i,"Nup159.1_10@%d"%i, 1.0,   f_n82+"Nup159.txt", "Nup159", "BEADS",  " ", (1082,1116,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
            domains.append(("Nup159.2@%d"%i,"Nup159.2_10@%d"%i, 0.9,   f_n82+"Nup159.txt", "Nup159", "BEADS",  " ", (1082,1116,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
            domains.append(("Nsp1.1@%d"%i,  "Nsp1.1_10@%d"%i,   0.3,   f_n82+"Nsp1.txt",   "Nsp1",   "BEADS",  " ", ( 601, 636,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
            domains.append(("Nsp1.2@%d"%i,  "Nsp1.2_10@%d"%i,   0.38,  f_n82+"Nsp1.txt",   "Nsp1",   "BEADS",  " ", ( 601, 636,0),  None, beadsize100, None,  None,  0,  None,  None,  None))
        domains.append(("Nup159.1@%d"%i,"Nup159.1@%d"%i,    1.0,   f_n82+"Nup159.txt", "Nup159", n82_pdb,  "E", (1117,1460,0),  None, beadsize,    None,  None,  6,  None,  None,  None))
        domains.append(("Nup159.2@%d"%i,"Nup159.2@%d"%i,    0.9,   f_n82+"Nup159.txt", "Nup159", n82_pdb,  "F", (1117,1460,0),  None, beadsize,    None,  None,  6,  None,  None,  None))
        domains.append(("Nsp1.1@%d"%i,  "Nsp1.1@%d"%i,      0.3,   f_n82+"Nsp1.txt",   "Nsp1",   n82_pdb,  "G", ( 637, 823,0),  None, beadsize,    None,  None,  4,  None,  None,  None))
        domains.append(("Nsp1.2@%d"%i,  "Nsp1.2@%d"%i,      0.38,  f_n82+"Nsp1.txt",   "Nsp1",   n82_pdb,  "H", ( 637, 823,0),  None, beadsize,    None,  None,  4,  None,  None,  None))
        #domains.append(("Nup116.1@%d"%i,"Nup116.1@%d"%i,    0.75,  f_n82+"Nup116.txt", "Nup116", n82_pdb,  "I", ( 751,1113,0),  None, beadsize25,  None,  None,  1,  None,  None,  None))
        #domains.append(("Nup116.2@%d"%i,"Nup116.2@%d"%i,    0.8,   f_n82+"Nup116.txt", "Nup116", n82_pdb,  "J", ( 751,1113,0),  None, beadsize25,  None,  None,  1,  None,  None,  None))

##########################
# Nic96 complex
##########################
if (is_nic96):
    n96_rb = 196
    domains.append(("Nic96.1",  "Nic96.1_1",  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (  1,106,0),  gmm,  beadsize,    n96_rb, [n96_rb], 2,  " ",   " ",  None, False))
    domains.append(("Nic96.1",  "Nic96.1_2",  0.25, f_n96+"Nic96.txt", "YFR002W", "BEADS",  " ",  (107,204,0),  None, beadsize25,  n96_rb, [n96_rb], 0,  None,  None, None))
    #domains.append(("Nic96.1",  "Nic96.1_2",  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (107,204,0),  gmm,  beadsize,    1097,   [n96_rb], 1,  " ",   " ",  None, False))
    domains.append(("Nic96.1",  "Nic96.1",    0.25, f_n96+"Nic96.txt", "YFR002W", n961_pdb, "A",  (205,839,0),  gmm,  beadsize,    None,   None,     3,  " ",   " ",  None, False))
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
    domains.append(("Nic96.2",  "Nic96.2_2",  0.25, f_n96+"Nic96.txt", "YFR002W", "BEADS",  " ",  (107,204,0),  None, beadsize25,  n96_rb, [n96_rb], 0,  None,  None, None))
    #domains.append(("Nic96.2",  "Nic96.2_2",  0.25, f_n96+"Nic96.txt", "YFR002W", n96_pdb,  "A",  (107,204,0),  gmm,  beadsize,    2097,   [n96_rb], 1,  gmm_f+"Nic96.1_2.txt", gmm_f+"Nic96.1_2.mrc", None, False))
    domains.append(("Nic96.2",  "Nic96.2",    0.25, f_n96+"Nic96.txt", "YFR002W", n962_pdb, "A",  (205,839,0),  gmm,  beadsize,    None,   None,     3,  " ",   " ",  None, False))
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
        domains.append(("Nic96.2@%d"%i,  "Nic96.2@%d"%i,    0.25, f_n96+"Nic96.txt", "YFR002W", n962_pdb, "A",  (205,839,0),  gmm_c,  beadsize,    None, None, 3,  gmm_f+"Nic96.2.txt",   gmm_f+"Nic96.2.mrc",   None, False))
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
    domains.append(("Nup157",    "Nup157n",       0.0,  f_npc+"Nup157.txt", "YER105C", n157N_pdbfile,  "A",  (  1, 892,0),  gmm,  beadsize25,  None, None, 3, " ",   " ",  None, False))
    domains.append(("Nup157",    "Nup157c",       0.0,  f_npc+"Nup157.txt", "YER105C", n157C_pdbfile,  "A",  (893,1391,0),  gmm,  beadsize25,  None, None, 3, " ",   " ",  None, False))
    domains.append(("Nup170",    "Nup170n",       0.1,  f_npc+"Nup170.txt", "Nup170",  n170N_pdbfile,  "A",  (  1, 992,0),  gmm,  beadsize25,  None, None, 3, " ",   " ",  None, False))
    domains.append(("Nup170",    "Nup170c",       0.1,  f_npc+"Nup170.txt", "Nup170",  n170C_pdbfile,  "A",  (993,1502,0),  gmm,  beadsize25,  None, None, 3, " ",   " ",  None, False))
    domains.append(("Nup188",    "Nup188",        0.85, f_npc+"Nup188.txt", "YML103C", n188_pdbfile,   "A",  (  1,1655,0),  gmm,  beadsize25,  None, None, 6, " ",   " ",  None, False))
    domains.append(("Nup192",    "Nup192",        0.75, f_npc+"Nup192.txt", "YJL039C", n192_pdbfile,   "A",  (  1,1683,0),  gmm,  beadsize25,  None, None, 6, " ",   " ",  None, False))
    for i in clones_range_A:
    #for i in range(11,12):
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
    domains.append(("Nup53",     "Nup53",         0.0,  f_npc+"Nup53.txt",  "YMR153W", n53_pdbfile,  "A", (   1, 475,0),  gmm,  beadsize50,  None, None,  2,  " ",   " ",  None, False))
    domains.append(("Nup59",     "Nup59",         0.66, f_npc+"Nup59.txt",  "YDL088C", n59_pdbfile,  "A", (   1, 528,0),  gmm,  beadsize50,  None, None,  2,  " ",   " ",  None, False))
    domains.append(("Ndc1",      "Ndc1",          0.8,  f_npc+"Ndc1.txt",   "YML031W", "BEADS",      " ", (   1, 655,0),  gmm,  beadsize100, None, None,  0,  None,  None, None, False))
    domains.append(("Pom34",     "Pom34",         0.9,  f_npc+"Pom34.txt",  "YLR018C", "BEADS",      " ", (   1, 299,0),  gmm,  beadsize50,  None, None,  0,  None,  None, None, False))
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
        domains.append(("Nup116.1", "Nup116.1_10", 0.75,f_n82+"Nup116.txt", "Nup116",  "BEADS",      " ", (  1, 750,0), None, beadsize100, 1116, [1116], 0,  None,  None, None, False))
        domains.append(("Nup116.2", "Nup116.2_10", 0.8, f_n82+"Nup116.txt", "Nup116",  "BEADS",      " ", (  1, 750,0), None, beadsize100, 2116, [2116], 0,  None,  None, None, False))
        domains.append(("Nup42",    "Nup42_10",    0.6, f_npc+"Nup42.txt",  "YDR192C", "BEADS",      " ", (  1, 363,0), None, beadsize100, 42,   [42],   0,  None,  None, None, False))
    domains.append(("Nup100.1", "Nup100.1",    0.2, f_npc+"Nup100.txt", "YKL068W", n100_pdbfile, "A", (551, 959,0), None, beadsize25,  1100, [1100],2,  " ",   " ",  None, False))
    domains.append(("Nup100.2", "Nup100.2",    0.4, f_npc+"Nup100.txt", "YKL068W", n100_pdbfile, "B", (551, 959,0), None, beadsize25,  2100, [2100],2,  " ",   " ",  None, False))
    domains.append(("Nup116.1", "Nup116.1",    0.75,f_n82+"Nup116.txt", "Nup116",  n116_pdbfile, "I", (751,1113,0), None, beadsize25,  None,  None, 1,  " ",   " ",  None, False))
    domains.append(("Nup116.2", "Nup116.2",    0.8, f_n82+"Nup116.txt", "Nup116",  n116_pdbfile, "J", (751,1113,0), None, beadsize25,  None,  None, 1,  " ",   " ",  None, False))
    domains.append(("Nup42",    "Nup42",       0.6, f_npc+"Nup42.txt",  "YDR192C", "BEADS",      " ", (364, 430,0), None, beadsize50,  42,   [42],  0,  None,  None, None, False))
    domains.append(("Gle1",     "Gle1_10",     0.8, f_npc+"Gle1.txt",   "YDL207W", Gle1N_pdbfile,"A", (  1, 239,0), None, beadsize50,  610,  [611], 2,  " ",   " ",  None, False))
    domains.append(("Gle1",     "Gle1",        0.8, f_npc+"Gle1.txt",   "YDL207W", Gle1C_pdbfile,"B", (240, 538,0), None, beadsize25,  611,  [611], 2,  " ",   " ",  None, False))
    #domains.append(("Gle2.1",   "Gle2.1",      0.9, f_npc+"Gle2.txt",   "YER107C", Gle2_pdbfile, "A", (  1, 365,0), None, beadsize25,  1612, [1612],2,  " ",                 " ",                None, False))
    #domains.append(("Gle2.2",   "Gle2.2",      1.0, f_npc+"Gle2.txt",   "YER107C", Gle2_pdbfile, "A", (  1, 365,0), None, beadsize25,  2612, [2612],2,  gmm_f+"Gle2.1.txt",  gmm_f+"Gle2.1.mrc", None, False))
    for i in clones_range_B:
        if (is_FG):
            domains.append(("Nup100.1@%d"%i, "Nup100.1_10@%d"%i, 0.2, f_npc+"Nup100.txt", "YKL068W", "BEADS",      " ", (  1, 550,0), None, beadsize100, None, None, 0,  None,  None, None))
            domains.append(("Nup100.2@%d"%i, "Nup100.2_10@%d"%i, 0.4, f_npc+"Nup100.txt", "YKL068W", "BEADS",      " ", (  1, 550,0), None, beadsize100, None, None, 0,  None,  None, None))
            domains.append(("Nup116.1@%d"%i, "Nup116.1_10@%d"%i, 0.75,f_n82+"Nup116.txt", "Nup116",  "BEADS",      " ", (  1, 750,0), None, beadsize100, None, None, 0,  None,  None, None))
            domains.append(("Nup116.2@%d"%i, "Nup116.2_10@%d"%i, 0.8, f_n82+"Nup116.txt", "Nup116",  "BEADS",      " ", (  1, 750,0), None, beadsize100, None, None, 0,  None,  None, None))
            domains.append(("Nup42@%d"%i,    "Nup42_10@%d"%i,    0.6, f_npc+"Nup42.txt",  "YDR192C", "BEADS",      " ", (  1, 363,0), None, beadsize100, None, None, 0,  None,  None, None))
        domains.append(("Nup100.1@%d"%i, "Nup100.1@%d"%i,    0.2, f_npc+"Nup100.txt", "YKL068W", n100_pdbfile, "A", (551, 959,0), None, beadsize25,  None, None, 2,  None,  None, None))
        domains.append(("Nup100.2@%d"%i, "Nup100.2@%d"%i,    0.4, f_npc+"Nup100.txt", "YKL068W", n100_pdbfile, "B", (551, 959,0), None, beadsize25,  None, None, 2,  None,  None, None))
        domains.append(("Nup116.1@%d"%i, "Nup116.1@%d"%i,    0.75,f_n82+"Nup116.txt", "Nup116",  n116_pdbfile, "I", (751,1113,0), None, beadsize25,  None, None, 1,  None,  None, None))
        domains.append(("Nup116.2@%d"%i, "Nup116.2@%d"%i,    0.8, f_n82+"Nup116.txt", "Nup116",  n116_pdbfile, "J", (751,1113,0), None, beadsize25,  None, None, 1,  None,  None, None))
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
        domains.append(("Nup1",     "Nup1_10",     0.8, f_npc+"Nup1.txt",   "YOR098C", "BEADS",      " ", (352,1076,0), None, beadsize100, 1,    None,   0,  None,  None, None, False))
        domains.append(("Nup60.1",  "Nup60.1_10",  0.86,f_npc+"Nup60.txt",  "YAR002W", "BEADS",      " ", (399, 539,0), None, beadsize100, 60,   None,   0,  None,  None, None, False))
        domains.append(("Nup60.2",  "Nup60.2_10",  0.93,f_npc+"Nup60.txt",  "YAR002W", "BEADS",      " ", (399, 539,0), None, beadsize100, 60,   None,   0,  None,  None, None, False))
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
    domains.append((    "Mlp1",      "Mlp1",      0.0, f_npc+"Mlp1.txt",   "YKR095W", "BEADS",   " ", (238, 716,0), gmm,  beadsize50,  9191,  None,   0,  None,  None, None))
    domains.append((    "Mlp2",      "Mlp2",      0.2, f_npc+"Mlp2.txt",   "YIL149C", "BEADS",   " ", (215, 690,0), gmm,  beadsize50,  9191,  None,   0,  None,  None, None))
    for i in clones_range_B:
        domains.append(("Mlp1@%d"%i, "Mlp1@%d"%i, 0.0, f_npc+"Mlp1.txt",   "YKR095W", "BEADS",   " ", (238, 716,0), None, beadsize50,  None,  None,   0,  None,  None, None))
        domains.append(("Mlp2@%d"%i, "Mlp2@%d"%i, 0.2, f_npc+"Mlp2.txt",   "YIL149C", "BEADS",   " ", (215, 690,0), None, beadsize50,  None,  None,   0,  None,  None, None))

#####################################################
# Model Building
#####################################################
bm1 = IMP.pmi1.macros.BuildModel1(simo)
bm1.set_gmm_models_directory(gmm_f)

if (True):
    if (is_n84):
        n84=['Nup84', 'Nup85', 'Nup120', 'Nup133', 'Nup145c', 'Seh1', 'Sec13']
        for d in list(n84):
            bm1.set_rmf_file(d, "../data_npc/Outer_ring_rmfs/OR_833_0_best.rmf3", 0)
            #bm1.set_rmf_file(d, "../data_nup82/rmfs/B_8_1-95.rmf3", 0)

    if (is_n82):
        n82=['Dyn2.1', 'Dyn2.2', 'Nup82.1', 'Nup82.2', 'Nup159.1', 'Nup159.2', 'Nsp1.1', 'Nsp1.2']
        for d in list(n82):
            bm1.set_rmf_file(d, "../data_npc/Outer_ring_rmfs/OR_833_0_best.rmf3", 0)
            #if (is_FG): bm1.set_rmf_file(d, "../data_nup82/rmfs/B_8_1-95.rmf3", 0)
            #else:       bm1.set_rmf_file(d, "../data_nup82/rmfs/B_8_1-95_FGtruncated.rmf3", 0)

    if (is_nic96):
        nic96=['Nic96.1', 'Nsp1.3', 'Nup49.1', 'Nup57.1', 'Nic96.2', 'Nsp1.4', 'Nup49.2', 'Nup57.2']
        for d in list(nic96):
            bm1.set_rmf_file(d, "../data_npc/Inner_ring_rmfs/Nic96complex_initial.rmf3", 0)

    if (is_inner_ring):
        inner_ring=['Nup157', 'Nup170', 'Nup188', 'Nup192']
        for d in list(inner_ring):
            bm1.set_rmf_file(d, "../data_npc/Inner_ring_rmfs/IR_865_0_final.rmf3", 0)

    if (is_membrane):
        inner_ring=['Nup53', 'Nup59', 'Ndc1', 'Pom34']
        for d in list(inner_ring):
            bm1.set_rmf_file(d, "../data_npc/Inner_ring_rmfs/IR_865_0_final.rmf3", 0)
        bm1.set_rmf_file('Pom152', "../data_npc/Pom152_rmfs/Pom152_0_final.rmf3", 0)
    """
    if (is_cytoplasm):
        n116=['Nup116.1', 'Nup116.2']
        for d in list(n116):
            if (is_FG): bm1.set_rmf_file(d, "../data_nup82/rmfs/B_8_1-95.rmf3", 0)
            else:       bm1.set_rmf_file(d, "../data_nup82/rmfs/B_8_1-95_FGtruncated.rmf3", 0)
    """
    if is_basket:
        Mlps=['Mlp1', 'Mlp2']
        for d in list(Mlps):
            bm1.set_rmf_file(d, "../data_npc/Mlps_1.rmf3", 0)
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
bm1.scale_bead_radii(100, 0.6)

if sys.version_info[0] >= 3:
    def open_csv(fname):
        return open(fname, newline='')
else:
    def open_csv(fname):
        return open(fname, 'rb')

class SAXSFits(object):
    """Parse the SAXS csv file and add suitable fit data to the mmCIF file"""
    saxs_dir = '../input_data_files/SAXS'
    seqrange_re = re.compile('(\d+)\s*\-\s*(\d+)')
    sasbdb_re = re.compile('/data/(SASDB\w+)')

    def __init__(self, po):
        self.po = po

    def add_from_csv(self, model):
        import csv
        with open_csv(os.path.join(self.saxs_dir, 'Table6_SAXS.csv')) as fh:
            for row in csv.DictReader(fh):
                if row['FoXS fit score']:
                    self._add_one(model, row)

    def _add_one(self, model, row):
        protid = row['Protein ID']
        m = self.sasbdb_re.search(row['Notes'])
        if m:
            l = ihm.location.SASBDBLocation(m.group(1))
            detail = None
        else:
            profile = (glob.glob('%s/*/%s_*.sub' % (self.saxs_dir, protid))
                     + glob.glob('%s/*/%s_*.dat' % (self.saxs_dir, protid)))[0]
            l = ihm.location.InputFileLocation(profile,
                                 details = row['Notes'] if row['Notes']
                                                        else None)
        dataset = ihm.dataset.SASDataset(location=l)
        m = self.seqrange_re.match(row['Sequence coverage'])
        seqrange = (int(m.group(1)), int(m.group(2)))
        protein = row['Protein'].capitalize()
        # Exclude fit in the file with invalid seqrange (Nup82 has 713 residues)
        if protein == 'Nup82' and seqrange == (791,899):
            return
        # For proteins with multiple copies, map to the first
        if protein in ('Nup145', 'Nup100', 'Nup116', 'Nup82'):
            protein += '.1'

        self.po._add_foxs_restraint(model, protein, seqrange, dataset,
                                    row['Rg'], row['FoXS fit score'], None)


def create_8_spokes(simo, proteins, both_half_spokes, one_spoke=False):
    """Expand the 1 or 3 spoke model to encompass all 8 spokes.
       This is done by creating 7 or 5 more symmetry-related copies of each
       protein. If both_half_spokes is True, the protein exists in both
       half spokes, so create copies of both half spoke primaries."""
    def get_rotation(i):
        axis = IMP.algebra.Vector3D(0, 0, 1.0)
        if i % 2 == 0:
            angle = 0.25 * math.pi * ((i+2)//2)
        else:
            angle = -0.25 * math.pi * ((i+1)//2)
        return IMP.algebra.get_rotation_about_axis(axis, angle)
    # 5 copies of primary protein
    for i in range(2 if one_spoke else 4, 9):
        for p in proteins:
            simo.create_transformed_component("%s@%d" % (p, i), p,
                                              get_rotation(i - 2))
    if both_half_spokes:
        for i in range(12 if one_spoke else 14, 19):
            for p in proteins:
                simo.create_transformed_component("%s@%d" % (p, i),
                                                  p+ "@11",
                                                  get_rotation(i - 12))

#####################################################
# apply the rotational symmetry
#####################################################
if (use_neighboring_spokes):
    if (is_n84):
        proteins = ['Nup84', 'Nup85', 'Nup120', 'Nup133', 'Nup145c',
                    'Seh1', 'Sec13']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=True)
        for protein in proteins:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
            #simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
            #simo.create_rotational_symmetry(protein+'@11', [protein+'@%d'%i for i in range(12,14)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
    if (is_n82):
        proteins = ['Dyn2.1', 'Dyn2.2', 'Nup82.1', 'Nup82.2', 'Nup159.1',
                    'Nup159.2', 'Nsp1.1', 'Nsp1.2']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=False)
#       for protein in proteins:
#           simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
    if (is_nic96):
        proteins = ['Nic96.1', 'Nic96.2', 'Nsp1.3', 'Nsp1.4', 'Nup49.1',
                    'Nup49.2', 'Nup57.1', 'Nup57.2']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=True)
        for protein in proteins:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
            simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
            simo.create_rotational_symmetry(protein+'@11', [protein+'@%d'%i for i in range(12,14)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
    if (is_inner_ring):
        proteins = ['Nup157', 'Nup170', 'Nup188', 'Nup192']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=True)
        for protein in proteins:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
            #simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
            #simo.create_rotational_symmetry(protein+'@11', [protein+'@%d'%i for i in range(12,14)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
    if (is_membrane):
        proteins = ['Nup53', 'Nup59', 'Ndc1', 'Pom34', 'Pom152']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=True)
        for protein in proteins:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
            simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
            simo.create_rotational_symmetry(protein+'@11', [protein+'@%d'%i for i in range(12,14)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
    if (is_cytoplasm):
        #for protein in ['Nup100.1', 'Nup100.2', 'Nup116.1', 'Nup116.2', 'Nup42', 'Gle1', 'Gle2.1', 'Gle2.2']:
        proteins = ['Nup100.1', 'Nup100.2', 'Nup116.1', 'Nup116.2', 'Nup42',
                    'Gle1']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=False)
        for protein in proteins:
            simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
    if (is_nucleoplasm):
        proteins = ['Nup145.1', 'Nup145.2', 'Nup60.1', 'Nup60.2', 'Nup1']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=False)
        for protein in proteins:
            simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
    if (is_basket):
        proteins = ['Mlp1', 'Mlp2']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=False)
        for protein in proteins:
            simo.create_rotational_symmetry(protein, [protein+'@%d'%i for i in range(2,4)], rotational_axis=IMP.algebra.Vector3D(0, 0, 1.0), nSymmetry=8, skip_gaussian_in_clones=True)
else:
    if (is_n84):
        proteins = ['Nup84', 'Nup85', 'Nup120', 'Nup133', 'Nup145c', 'Seh1',
                    'Sec13']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=True,
                            one_spoke=True)
        for protein in proteins:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
    if is_n82:
        proteins = ['Dyn2.1', 'Dyn2.2', 'Nup82.1', 'Nup82.2', 'Nup159.1',
                    'Nup159.2', 'Nsp1.1', 'Nsp1.2']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=False,
                            one_spoke=True)
    if (is_nic96):
        proteins = ['Nic96.1', 'Nic96.2', 'Nsp1.3', 'Nsp1.4', 'Nup49.1',
                    'Nup49.2', 'Nup57.1', 'Nup57.2']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=True,
                            one_spoke=True)
        for protein in proteins:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
    if (is_inner_ring):
        proteins = ['Nup157', 'Nup170', 'Nup188', 'Nup192']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=True,
                            one_spoke=True)
        for protein in proteins:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
    if (is_membrane):
        proteins = ['Nup53', 'Nup59', 'Ndc1', 'Pom34', 'Pom152']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=True,
                            one_spoke=True)
        for protein in proteins:
            simo.create_rotational_symmetry(protein, [protein+'@11'], rotational_axis=IMP.algebra.Vector3D(1.0, 0, 0))
    if (is_cytoplasm):
        proteins = ['Nup100.1', 'Nup100.2', 'Nup116.1', 'Nup116.2', 'Nup42',
                    'Gle1']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=False,
                            one_spoke=True)
    if (is_nucleoplasm):
        proteins = ['Nup145.1', 'Nup145.2', 'Nup60.1', 'Nup60.2', 'Nup1']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=False,
                            one_spoke=True)
    if (is_basket):
        proteins = ['Mlp1', 'Mlp2']
        if inputs.mmcif and inputs.symmetry:
            create_8_spokes(simo, proteins, both_half_spokes=False,
                            one_spoke=True)


#####################################################
# rigidify floppy bodies
#####################################################
rigid_tuples = []
#for protein in ['Nup84', 'Nup85', (1,711,'Nup120'),(715,1037,'Nup120'), (1,480,'Nup133'),(490,1157,'Nup133'), 'Nup145c', 'Seh1', 'Sec13']:
for protein in ['Nup84', 'Nup85', 'Nup120', 'Nup133', 'Nup145c', 'Seh1', 'Sec13']:
    rigid_tuples.append(protein)
for protein in ['Dyn2.1', 'Dyn2.2', 'Nup82.1', 'Nup82.2']:
    rigid_tuples.append(protein)
#for protein in [(1117,1460,'Nup159.1'),(1117,1460,'Nup159.2'), (637,823,'Nsp1.1'),(637,823,'Nsp1.2'), (966,1113,'Nup116.1'),(966,1113,'Nup116.2')]:
for protein in ['Nup159.1', 'Nup159.2', 'Nsp1.1', 'Nsp1.2']:
    rigid_tuples.append(protein)
#for protein in [(1,56,'Nic96.1'),(205,839,'Nic96.1'), (637,823,'Nsp1.3'), (270,472,'Nup49.1'), (287,541,'Nup57.1')]:
for protein in [(637,823,'Nsp1.3'), (270,472,'Nup49.1'), (287,541,'Nup57.1')]:
    rigid_tuples.append(protein)
#for protein in [(1,56,'Nic96.2'),(205,839,'Nic96.2'), (637,823,'Nsp1.4'), (270,472,'Nup49.2'), (287,541,'Nup57.2')]:
for protein in [(637,823,'Nsp1.4'), (270,472,'Nup49.2'), (287,541,'Nup57.2')]:
    rigid_tuples.append(protein)
#for protein in [(88,892,'Nup157'),(900,1391,'Nup157'), (98,992,'Nup170'),(1000,1502,'Nup170'), 'Nup188', 'Nup192']:
for protein in ['Nup157', 'Nup170', 'Nup188', 'Nup192']:
    rigid_tuples.append(protein)
#for protein in [(201,410,'Nup53'), (251,402,'Nup59'), 'Mlp1', 'Mlp2']:
for protein in [(248,475,'Nup53'), (266,528,'Nup59'), 'Ndc1', 'Pom34', (379,1337,'Pom152'), 'Mlp1', 'Mlp2']:
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
    #simo.shuffle_configuration(bounding_box=((350, -150, 50), (1050, 200, 550)), ignore_initial_coordinates=True, cutoff=1.0, niterations=1000)
    simo.shuffle_configuration(bounding_box=((150, -150, 0), (550, 150, 350)), ignore_initial_coordinates=True, cutoff=1.0, niterations=1000)


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
            low_res = IMP.pmi1.tools.select(self.simo, resolution=res_ev,
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

if (use_Composite):
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
                  for r in (rot, rot.get_inverse())]
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
        IMP.pmi1.tools.add_restraint_to_model(simo.m, rsr)

# Add absolute position restraints
if inputs.mmcif:
    # Use final restraint values, not those for the initial modeling
    use_FG_anchor = True
    use_MembraneExclusion = True
    # Don't add Mlp1/2 restraints, since these weren't ultimately used
    is_basket = False
    exec(open("positional_restraints_final.py").read())
else:
    exec(open("positional_restraints_initial.py").read())

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
                                                        #inner_slope = 0.02,
                                                        inner_slope = 0.01,
                                                        filelabel = "wtDSS",
                                                        label = "wtDSS",
                                                        attributes_for_label = ["XLUniqueID"],
                                                        csvfile = True)
    # Point to the raw mass spec data and peaklists used to derive the
    # crosslinks.
    r = ihm.location.Repository(doi="10.5281/zenodo.1149746",
        url='https://zenodo.org/record/1149746/files/TODO')
    l = ihm.location.InputFileLocation(repo=r, path='TODO',
                         details='All raw mass spectrometry files and '
                                 'peaklists used in the study')
    d = ihm.dataset.MassSpecDataset(location=l)
    xl1.dataset.add_primary(d)

    xl1.add_to_model()
    xl1.set_weight(10.0)        # play with the weight
    sampleobjects.append(xl1)
    outputobjects.append(xl1)
    xl1.set_psi_is_sampled(False)
    psi2 = xl1.get_psi(1.0)[0]
    psi2.set_scale(0.05)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
    print("\nEVAL 1 : ", sf.evaluate(False), " (after applying the XL restraint) - ", rank)
    XL_restraints = [xl1]

    if inputs.mmcif:
        # Also add Mlp-specific crosslinks to the final deposition
        xl2 = IMP.pmi1.restraints.crosslinking.ISDCrossLinkMS(simo,
                 '../data_npc/XL_Merged_wholeNPC_MLPs.csv', length = 26.0,
                 slope = 0.00, columnmapping = columnmap, ids_map = ids_map,
                 resolution = 1.0, inner_slope = 0.01, filelabel = "wtDSS",
                 label = "wtDSS", attributes_for_label = ["XLUniqueID"],
                 csvfile = True)
        xl2.dataset.add_primary(d)
        xl2.add_to_model()
        xl2.set_weight(10.0)
        sampleobjects.append(xl2)
        outputobjects.append(xl2)
        xl2.set_psi_is_sampled(False)
        psi2 = xl2.get_psi(1.0)[0]
        psi2.set_scale(0.05)
        XL_restraints.append(xl2)
else:
    XL_restraints = None


#####################################################
# Restraints setup
# Distance restraints for XL cliques involving the membrane nups
#####################################################
if (is_inner_ring and is_membrane):
    db = IMP.pmi1.tools.get_db_from_csv('../data_npc/XL_cliques.csv')
    dist_max = 30.0

    for nxl, entry in enumerate(db):
        #print(nxl, entry)
        mol1 = entry["Protein 1"]
        res1 = int(entry["Residue 1"])
        mol2 = entry["Protein 2"]
        res2 = int(entry["Residue 2"])

        dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo, (res1,res1,mol1), (res2,res2,mol2), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        temp_label = mol1 + "_" + str(res1) + "-" + mol2 + "_" + str(res2)
        dr.set_label(temp_label)
        dr.add_to_model()
        dr.set_weight(dr_weight)
        outputobjects.append(dr)
        print(dr.get_output())

    print("\nDistance Restraints applied for XL cliques !!")
    print("weight = ", dr_weight, "dist_min = ", dist_min, "dist_max = ", dist_max, "\n")


#####################################################
# Restraints setup
# Sampling Boundary Restraint for the inner ring
#####################################################
if (use_sampling_boundary):
    main_spoke = [];  other_spokes = [];    main_spoke_hier_name = []
    for entry in domains:
        if 'Nup170n@11' in entry[1]:
            main_spoke.append(entry[0])
            main_spoke_hier_name.append(entry[1])
        elif 'Nup170c@11' in entry[1]:
            other_spokes.append(entry[0])
        elif 'Nup170n' in entry[1]:
            other_spokes.append(entry[0])
        elif 'Nup170c' in entry[1]:
            main_spoke.append(entry[0])
            main_spoke_hier_name.append(entry[1])
        elif '@' in entry[0]:
            other_spokes.append(entry[0])
        elif 'Ndc1' in entry[0]:
            other_spokes.append(entry[0])
        elif 'Pom34' in entry[0]:
            other_spokes.append(entry[0])
        elif 'Pom152' in entry[0]:
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
    #mass *= 1.2 * 2.0           # 1.2 for adjustment of the GMM (after removing flexible GMMs) and 2.0 for approximation of the NPC spoke mass
    mass *= 1.2           # 1.2 for adjustment of the GMM (after removing flexible GMMs)
    print ("Total mass for the Sampling Boundary EM restraint = ", mass)
    sbr = IMP.pmi1.restraints.em.GaussianEMRestraint(resdensities,
                                                    '../data_npc/em_gmm_model/SJ_SamplingBoundary.gmm.15.txt',
                                                    target_mass_scale=mass,
                                                    slope=0.01,
                                                    #slope=0.0000001,
                                                    target_radii_scale=3.0)
    sbr.set_label("Sampling_Boundary")
    sbr.add_to_model()
    sbr.set_weight(2.0)        # play with the weight
    #sbr.center_model_on_target_density(simo)
    outputobjects.append(sbr)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
    print("\nEVAL 2 : ", sf.evaluate(False), " (after applying the Sampling Boundary EM restraint) - ", rank)


#####################################################
# 1st Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
print("\nEVAL 3 : ", sf.evaluate(False), " (initial) - ", rank)

if not inputs.dry_run:
    simo.optimize_floppy_bodies(300)
print("\nEVAL 4 : ", sf.evaluate(False), " (after calling optimize_floppy_bodies(300)) - ", rank)

mc1 = IMP.pmi1.macros.ReplicaExchange0(m,
                                    simo,
                                    monte_carlo_sample_objects = sampleobjects,
                                    output_objects = outputobjects,
                                    crosslink_restraints = XL_restraints,
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 5.0,
                                    number_of_best_scoring_models = 0,
                                    monte_carlo_steps = 10,
                                    number_of_frames = 10 if inputs.test
                                                          else 500,
                                    write_initial_rmf = True,
                                    initial_rmf_name_suffix = "initial",
                                    stat_file_name_suffix = "stat",
                                    best_pdb_name_suffix = "model",
                                    do_clean_first = True,
                                    do_create_directories = True,
                                    global_output_directory = "1_pre_output",
                                    rmf_dir = "rmfs/",
                                    best_pdb_dir = "pdbs/",
                                    test_mode=simo.dry_run,
                                    replica_stat_file_suffix = "stat_replica")
mc1.execute_macro()
rex1 = mc1.get_replica_exchange_object()
print("\nEVAL 5 : ", sf.evaluate(False), " (after performing the pre_sampling) - ", rank)
#exit(0)


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
    mass *= 1.2 * 2.0           # 1.2 for adjustment of the GMM (after removing flexible GMMs) and 2.0 for approximation of the NPC spoke mass
    print ("Total mass for the EM restraint = ", mass)
    gem = IMP.pmi1.restraints.em.GaussianEMRestraint(resdensities,
                                                    '../data_npc/em_gmm_model/SJ_cropped_sym8_avg_monomer_final_rotated_adjusted90.gmm.400.txt',
                                                    target_mass_scale=mass,
                                                    #slope=0.0000005,
                                                    slope=0.0000001,
                                                    target_radii_scale=3.0,
                                                    representation=simo)
    gem.add_to_model()
    gem.set_weight(1000.0)        # play with the weight

    # Point to the original map in EMDB
    l = ihm.location.EMDBLocation('EMD-7321')
    emdb = ihm.dataset.EMDensityDataset(location=l)
    gem.dataset.add_primary(emdb)

    # Point back to Cryo-ET raw data (tilt series)
    l = ihm.location.EMPIARLocation('EMPIAR-10155')
    d = ihm.dataset.EMMicrographsDataset(location=l)
    emdb.add_primary(d)

    #gem.center_model_on_target_density(simo)
    outputobjects.append(gem)

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
    print("\nEVAL 6 : ", sf.evaluate(False), " (after applying the EM 3D restraint) - ", rank)


#####################################################
# 2nd Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
mc2 = IMP.pmi1.macros.ReplicaExchange0(m,
                                    simo,
                                    monte_carlo_sample_objects = sampleobjects,
                                    output_objects = outputobjects,
                                    crosslink_restraints = XL_restraints,
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 5.0,
                                    number_of_best_scoring_models = 0,
                                    monte_carlo_steps = 10,
                                    number_of_frames = 30 if inputs.test
                                                          else 3000,
                                    write_initial_rmf = True,
                                    initial_rmf_name_suffix = "initial",
                                    stat_file_name_suffix = "stat",
                                    best_pdb_name_suffix = "model",
                                    do_clean_first = True,
                                    do_create_directories = True,
                                    global_output_directory = "2_XL_EM_output",
                                    rmf_dir = "rmfs/",
                                    best_pdb_dir = "pdbs/",
                                    test_mode=simo.dry_run,
                                    #replica_stat_file_suffix = "stat_replica")
                                    replica_stat_file_suffix = "stat_replica",
                                    replica_exchange_object = rex1)
mc2.execute_macro()
rex2 = mc2.get_replica_exchange_object()
print("\nEVAL 7 : ", sf.evaluate(False), " (after performing the XL_EM_sampling) - ", rank)
#exit(0)


#####################################################
# Restraints setup
# Excluded Volume restraint for components in the main spoke
#####################################################
if (use_ExcludedVolume):
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

    ev1 = IMP.pmi1.restraints.stereochemistry.ExcludedVolumeSphere(simo,
                                                                 included_objects = included_objects,
                                                                 #other_objects = other_objects,
                                                                 resolution = res_ev)
    ev1.set_label('main_spoke')
    ev1.add_to_model()
    ev1.set_weight(0.3)
    outputobjects.append(ev1)
    print(ev1.get_output())
    print("ExcludedVolumeSphere1 for the main spoke !!\n")

    if (use_neighboring_spokes):
        ev2 = IMP.pmi1.restraints.stereochemistry.ExcludedVolumeSphere(simo,
                                                                     included_objects = included_objects,
                                                                     other_objects = other_objects,
                                                                     resolution = res_ev)
        ev2.set_label('bipartite')
        ev2.add_to_model()
        ev2.set_weight(0.3)
        outputobjects.append(ev2)
        print(ev2.get_output())
        print("ExcludedVolumeSphere2 between the main spoke and the neighboring spokes !!\n")

    sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
    print("\nEVAL 8 : ", sf.evaluate(False), " (after applying the Excluded Volume restraint) - ", rank)

"""
#####################################################
# 3rd Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
mc3 = IMP.pmi1.macros.ReplicaExchange0(m,
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
                                    test_mode=simo.dry_run,
                                    replica_stat_file_suffix = "stat_replica",
                                    replica_exchange_object = rex2)
mc3.execute_macro()
rex3 = mc3.get_replica_exchange_object()
print("\nEVAL 9 : ", sf.evaluate(False), " (after performing the EV_sampling) - ", rank)
"""


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
                                    test_mode=simo.dry_run,
                                    replica_stat_file_suffix = "stat_replica",
                                    replica_exchange_object = rex2)
                                    #replica_exchange_object = rex3)
mc4.execute_macro()
print("\nEVAL 10 : ", sf.evaluate(False), " (final evaluation) - ", rank)
#exit(0)

if inputs.mmcif:
    # Link entities to UniProt
    if hasattr(ihm, 'reference'):
        for subunit, accession, db_begin in (
                ('Nup84', 'P52891', 1), ('Nup85', 'P46673', 1),
                ('Nup120', 'P35729', 1), ('Nup133', 'P36161', 1),
                ('Nup145c', 'P49687', 606), ('Seh1', 'P53011', 1),
                ('Sec13', 'Q04491', 1), ('Dyn2.1', 'Q02647', 1),
                ('Nup82.1', 'P40368', 1), ('Nup159.1', 'P40477', 1),
                ('Nsp1.1', 'P14907', 1), ('Nic96.1', 'P34077', 1),
                ('Nup49.1', 'Q02199', 1), ('Nup57.1', 'P48837', 1),
                ('Nup157', 'P40064', 1), ('Nup170', 'P38181', 1),
                ('Nup188', 'P52593', 1), ('Nup192', 'P47054', 1),
                ('Nup53', 'Q03790', 1), ('Nup59', 'Q05166', 1),
                ('Ndc1', 'P32500', 1), ('Pom34', 'Q12445', 1),
                ('Pom152', 'P39685', 1), ('Nup100.1', 'Q02629', 1),
                ('Nup116.1', 'Q02630', 1), ('Nup42', 'P49686', 1),
                ('Gle1', 'Q12315', 1), ('Nup145.1', 'P49687', 1),
                ('Nup1', 'P20676', 1), ('Nup60.1', 'P39705', 1),
                ('Mlp1', 'Q02455', 1), ('Mlp2', 'P40457', 1)):
            ref = ihm.reference.UniProtSequence.from_accession(accession)
            ref.alignments.append(ihm.reference.Alignment(db_begin=db_begin))
            e = po.asym_units[subunit].entity.references.append(ref)

    if inputs.one_spoke:
        framework_rmf = '../results/RMF_files/cluster0_47-35_1spoke.rmf3'
    else:
        framework_rmf = '../results/RMF_files/cluster0_47-35_3spokes.rmf3'
    fgs_rmf = 'npc_fg_2018/RepresentativeEnsemble/modelN11_101.rmf'

    # Point to ChimeraX command scripts to display aspects of the mmCIF model
    if inputs.symmetry:
        l = ihm.location.VisualizationFileLocation(
                     path='../results/pdb-dev/chimerax/show-fg-ribbons.cxc',
                     details='Show FG repeats as ribbons')
        simo.add_metadata(l)
    l = ihm.location.VisualizationFileLocation(
                 path='../results/pdb-dev/chimerax/show-nic96-em2d.cxc',
                 details='Show fit of Nic96 complex against EM class averages')
    simo.add_metadata(l)

    # todo: fill in correct numbers
    pp = po._add_simple_postprocessing(num_models_begin=15000,
                                       num_models_end=2267)
    for c in simo.get_component_names():
        # Strip out 10-per-bead representation and Gaussians, since this
        # was previously done for the results RMF files (for speed)
        for reptop in simo.hier_dict[c].get_children():
            if 'Res:10' in reptop.get_name():
                IMP.atom.destroy(reptop)
        simo.set_coordinates_from_rmf(c, framework_rmf, 0, force_rigid_update=True, skip_gaussian_in_representation=True)
    den = {}
    den_nup = {}
    if inputs.symmetry:
        prefix = '8spokes-C8/C8_'
    elif inputs.one_spoke:
        prefix = '1spoke-C1/'
    else:
        prefix = '3spokes-C3/C3_'
    for copy in po.all_modeled_components:
        # One MRC file covers all symmetry copies of a Nup
        nup = copy.split('@')[0]
        # No coordinates for Nup42
        if nup == 'Nup42': continue
        if nup not in den_nup:
            den_nup[nup] = ihm.location.OutputFileLocation(
                          path='../results/localization_density_files_MRC/'
                               '%s%s.mrc' % (prefix, nup),
                          details="Localization density for %s" % nup)
        den[copy] = den_nup[nup]
    if inputs.one_spoke and not inputs.symmetry:
        # Only include full scaffold ensemble for 1-spoke model due to size
        r = ihm.location.Repository(doi="10.5281/zenodo.1194547",
             url='https://zenodo.org/record/1194547/files/scaffold-1spoke.dcd')
        f = ihm.location.OutputFileLocation(repo=r, path='.',
               details="All ensemble structures for scaffold")
    else:
        f = None
    c = po._add_simple_ensemble(pp, name="Scaffold cluster 1", num_models=5,
                                drmsd=1.0, num_models_deposited=1,
                                localization_densities=den, ensemble_file=f)
    scaffold_model = po.add_model(c.model_group)

    f = SAXSFits(po)
    f.add_from_csv(scaffold_model)

    # Add ensemble for FG repeats
    if inputs.symmetry:
        den = {}
        den_nup = {}
        for copy in po.all_modeled_components:
            # One MRC file covers all instances of a Nup (copies and symmetries)
            nup = copy.split('.')[0].split('@')[0]
            if nup in po.fgs.ranges:
                if nup not in den_nup:
                    den_nup[nup] = ihm.location.OutputFileLocation(
                                    path='npc_fg_2018/Densities/%s.mrc' % nup,
                                    details="Localization density for %s" % nup)
                den[copy] = den_nup[nup]
        po.fgs.create_assembly(po)
        # todo: add brownian dynamics script file, check numbers
        po._add_protocol()
        po._add_simple_dynamics(num_models_end=1000,
                                method="Brownian dynamics")
        pp = po._add_no_postprocessing(num_models=1000)
        f = ihm.location.OutputFileLocation(
               path='npc_fg_2018/RepresentativeEnsemble/fg_repeat_ensemble.dcd',
               details="All ensemble structures for FG repeats")
        c = po._add_simple_ensemble(pp, name="FG ensemble", num_models=1000,
                                    drmsd=1.0, num_models_deposited=1,
                                    localization_densities=den,
                                    ensemble_file=f)
        m = po.add_model(c.model_group, assembly=po.fgs.assembly,
                         representation=po.fgs.representation)
        # No restraints act on this model
        m._is_restrained = False
        po.fgs.add_bead_coordinates(fgs_rmf, m)

    # Add 2DEM restraint (note that we didn't use this for modeling, only
    # validation, and we actually used a different algorithm, but this
    # distinction isn't important for mmCIF output)
    pixel_size = 2.03
    image_size = 128
    n96_dir = '../validation/nic96_em2d'
    image_numbers = [6,25]
    images = ['%s/Images/Image-%d.pgm' % (n96_dir, num)
              for num in image_numbers]
    # todo: fill in correct # of micrographs
    em2d = IMP.pmi1.restraints.em2d.ElectronMicroscopy2D(simo, images,
                                                    resolution=1.0,
                                                    pixel_size = pixel_size,
                                                    image_resolution = 35.0,
                                                    projection_number = 10000,
                                                    micrographs_number = 800)
    # Point to the raw micrographs in EMPIAR
    l = ihm.location.EMPIARLocation('EMPIAR-10162')
    micrographs = ihm.dataset.EMMicrographsDataset(l)
    for d in em2d.datasets:
        d.add_primary(micrographs)

    em2d.add_to_model()
    # Add CCC and transformation to model (as if stored in a PMI stat file)
    sys.path.append('%s/Model_2B' % n96_dir)
    from get_transformations import get_transformations, get_centroid
    cccs = {}
    pat = re.compile('Score (\S+) (\d+) ccc= ([\d.]+)')
    with open('%s/Model_2B/C1_logs_35.txt' % n96_dir) as fh:
        for match in pat.findall(fh.read()):
            cccs[int(match[1])] = float(match[2])
            pdb_file = match[0]
    centroid = get_centroid('%s/Model_2B/%s' % (n96_dir, pdb_file))
    transforms = {}
    for num, trans in get_transformations(
                            "%s/Model_2B/Registration-Parameters" % n96_dir,
                            centroid, pixel_size, image_size):
        transforms[num] = trans
    stats = {}
    for i, image_num in enumerate(image_numbers):
        prefix = 'ElectronMicroscopy2D_None_Image%d' % (i+1)
        stats[prefix + '_CCC'] = cccs[image_num]
        rot = transforms[image_num].get_rotation().get_quaternion()
        for n in range(4):
            stats[prefix + '_Rotation%d' % n] = rot[n]
        trn = transforms[image_num].get_translation()
        for n in range(3):
            stats[prefix + '_Translation%d' % n] = trn[n]
    scaffold_model.em2d_stats = stats

    po.flush()
