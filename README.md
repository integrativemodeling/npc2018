# The entire 552-protein yeast NPC complex

These scripts demonstrate the use of [IMP](http://salilab.org/imp) in the modeling of the yeast NPC complex using data as described in Seung Joong Kim, et al.'s 2018 NPC article published in Nature. 

The scripts work with the [IMP](http://salilab.org/imp) (version 2.6).
A default build of IMP compiled with the IMP::npc module should work, but for most effective sampling, it should
be built with [MPI](http://integrativemodeling.org/nightly/doc/html/namespaceIMP_1_1mpi.html) so that replica exchange can be used.

## List of files and directories:

- `input_data_files`		            contains all relevant data
- `results`		                      contains resulting structures and output files
- `template`			                  contains modeling scripts

## Compiling IMP with NPC-specific module:
- Clone IMP version 385a178
- Clone this repository into imp/modules/npc/.
- Compile IMP

## Running the IMP scripts for the NPC complex:
- `cd scripts`
- `python test_NPC_scoring_functions.py & > test_NPC_scoring_functions.out` (on a single processor; prepend `mpirun -np 6` or similar if you built IMP with MPI support)

## Information

_Author(s)_: Seung Joong Kim

_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:

- Seung Joong Kim\*, Javier Fernandez-Martines\*, Ilona Nudelman\*, et al, [Integrative structure and Functional Anatomy of a Nuclear Pore Complex](http://www.nature.com/nature/journal/), Nature , 2018, in press.

