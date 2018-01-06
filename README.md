# The entire 552-protein yeast NPC complex

These scripts demonstrate the use of [IMP](http://salilab.org/imp) in the modeling of the yeast NPC complex using diverse types of data as described in Seung Joong Kim, et al.'s 2018 NPC article published in Nature.

The scripts work with the [IMP](http://salilab.org/imp) (version 2.6).
A default build of IMP compiled with the IMP::npc module should work, but for most effective sampling, it should
be built with [MPI](http://integrativemodeling.org/nightly/doc/html/namespaceIMP_1_1mpi.html) so that replica exchange can be used.

## List of files and directories:

- `input_data_files`		            contains all relevant data

  EMBD_final-yNPC_map_28A.mrc.gz  : Cryo-ET density map
  
  Table1_crosslinks.xlsx : 3,077 chemical cross-links
  
  SAXS.zip  : SAXS source data (147 SAXS profiles for 18 different Nups) file used for assessment of the NPC structure 
  
- `data_npc`		            contains representation files (PDB format)
- `results`		                      contains resulting structures and output files
- `template`			                  contains modeling scripts

  1_modeling_wholeNPC.py  : Initial modeling script

  2_modeling_wholeNPC_FG_anchor_EV.py : Intermediate modeling script

  3_modeling_wholeNPC_refinement.py : Refinement modeling script

  4_modeling_wholeNPC_refinement.py : Final refinement modeling script

## Compiling IMP with NPC-specific module:
- Clone IMP version 2.6
- Clone the parent repository (npc) into imp/modules/npc/.
- Compile IMP

## Running the IMP scripts for the NPC complex:
- `cd scripts`
- `python test_NPC_scoring_functions.py & > test_NPC_scoring_functions.out` (on a single processor; prepend `mpirun -np 6` or similar if you built IMP with MPI support)

## Information

_Author(s)_: Seung Joong Kim

_License_: [CC-BY-SA-4.0](https://creativecommons.org/licenses/by-sa/4.0/legalcode).
This work is freely available under the terms of the Creative Commons
Attribution-ShareAlike 4.0 International License.

_Publications_:
- Seung Joong Kim\*, Javier Fernandez-Martines\*, Ilona Nudelman\*, Yi Shi\*, Wenzhu Zhang\*, et al., [Integrative structure and Functional Anatomy of a Nuclear Pore Complex](http://www.nature.com/nature/journal/), Nature , 2018, in press.

