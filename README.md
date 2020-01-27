[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1194547.svg)](https://doi.org/10.5281/zenodo.1194547)

These scripts demonstrate the use of [IMP](http://salilab.org/imp) in the modeling of the yeast NPC complex using diverse types of data as described in Seung Joong Kim, et al.'s 2018 NPC article published in Nature.

The scripts work with the [IMP](http://salilab.org/imp) (version 2.6).
A default build of IMP compiled with the IMP::npc module should work, but for most effective sampling, it should
be built with [MPI](http://integrativemodeling.org/nightly/doc/html/namespaceIMP_1_1mpi.html) so that replica exchange can be used.

## List of files and directories:

- `input_data_files`		            contains all relevant data

  EMBD_final-yNPC_map_28A.mrc.gz  : Cryo-ET density map
  
  Table1_crosslinks.xlsx : a table of all 3,077 chemical cross-links
  
  SAXS.zip  : SAXS source data (147 SAXS profiles for 18 different Nups) file used for assessment of the NPC structure 
  
- `data_npc`		            contains representation (PDB or rmf3 formats), sequence, MODELLER files

   protein_fasta.*.txt : Sequence for each Nups
   
   *.pdb (or *.rmf3): Representation pdb (or rmf3) files for each Nups
   
   emd_5556.map : EM 3D density map for Nup192
   
   XL_optimized_ambiguity.csv : A subset of the chemical cross-links used for the refinement
   
   XL_Merged_wholeNPC_MLPs : A subset of the chemical cross-links used for the basket components of Mlp1 and Mlp2
   
- `data_npc\Inner_ring_rmfs`		  Representation for the inner ring components
- `data_npc\Outer_ring_rmfs`		  Representation for the outer ring components
- `data_npc\Pom152_rmfs`		  Representation for the membrane ring (Pom152) component
- `data_npc\WholeNPC_rmfs`		  Intermediate files
- `data_npc\XL_backup`		  Backup of Chemical cross-links data
- `data_npc\em_gmm_model`		  Gaussian Mixture Model (GMM) for each Nups
- `data_npc\generating_em_gmm`		  Python scripts to generate Gaussian Mixture Model (GMM) for each Nups
- `data_npc\input_contact_frequencies`		  Contact frequency data from the 2007 NPC topological model
- `data_npc\input_density_maps`		  Localization probability densities from the 2007 NPC topological model

- `data_nic96`		            contains Nic96 complex-specific representation (PDB or rmf3 formats), sequence (*.txt), cross-links (XLs folder), and EM2D (Nic96complex_classes.hdf) files
- `data_nup82`		            contains Nup82 complex-specific representation (PDB or rmf3 formats), sequence (*.txt), cross-links (*.csv) files 
- `data_nup84_2016`		            contains Nup84 complex-specific representation (PDB or rmf3 formats), sequence (*.txt), cross-links (*.csv), hhpred (*.pdf or *.webarchive), and EM2D (EM_image.png) files 
- `pom152`		            Integrative structure determination of the Pom152 membrane ring - needed to be separate.

- `validation` contains data not used in modeling
  - `nic96_em2d` fit of final Nic96 structure against EM class averages

- `results`		                      contains resulting structures and output files
- `template`			                  contains modeling scripts

  1_modeling_wholeNPC.py  : Initial modeling script

  2_modeling_wholeNPC_FG_anchor_EV.py : Intermediate modeling script

  3_modeling_wholeNPC_refinement.py : Refinement script

  4_modeling_wholeNPC_refinement.py : Final refinement script
  
  Clustering_RMSD_GPU.py : Clustering script using GPU

  modeling_Nic96complex_initial.py : Initial modeling script for the Nic96 complex components
   
  modeling_MLPs.py : Refinement script for the basket components
   
  modeling_pom152.py : Initial modeling script for the membrane-ring (pom152) component
  
  modeling_pom152_PlaneDihedral.py : Refinement script for the membrane-ring (pom152) component
  


## Compiling IMP with NPC-specific module:
- Clone IMP version 2.6
- Clone the parent repository (npc) into imp/modules/npc/.
- Compile IMP

## Running the IMP scripts for the NPC complex:

1. Inner-ring components (Nup157, Nup170, Nup188, Nup192, Nic96, Nup53, Nup59, Ndc1, Pom34, and Pom152 NTD)
- Initial: template/inner_ring/job_IR501-510.sh (running script for template/inner_ring/modeling_inner_ring_initial.py), generates "3IR_502_0.rmf3" (an initial model)
- Refinement : template/inner_ring/job_IR860-869_refinement.sh (running script for template/inner_ring/modeling_inner_ring_refinement.py, which reads "3IR_502_0.rmf3 for the starting coordinate), generates "IR_865_0_final.rmf3" (a refined model).

2. Outer-ring components (Nup82 and Nup84 complexes)
- Initial : 
- Refinement :

- `cd template`
- `python XX.py & > XX.out` (on a single processor; prepend `mpirun -np 6` or similar if you built IMP with MPI support)

or job_test4.sh

## Information

_Author(s)_: Seung Joong Kim

_License_: [CC-BY-SA-4.0](https://creativecommons.org/licenses/by-sa/4.0/legalcode).
This work is freely available under the terms of the Creative Commons
Attribution-ShareAlike 4.0 International License.

_Last known good IMP version_: [![build info](https://integrativemodeling.org/systems/28/badge.svg?branch=master)](https://integrativemodeling.org/systems/) [![build info](https://integrativemodeling.org/systems/28/badge.svg?branch=develop)](https://integrativemodeling.org/systems/)

_Publications_:
- Seung Joong Kim\*, Javier Fernandez-Martinez\*, Ilona Nudelman\*, Yi Shi\*, Wenzhu Zhang\*, et al., [Integrative structure and Functional Anatomy of a Nuclear Pore Complex](https://www.nature.com/articles/nature26003), Nature 555, 475-482, 2018.

