To generate mmCIF files suitable for deposition in PDB-Dev, use the
`../../template/1_modeling_wholeNPC.py` script with the following options:

 - For 1 spoke: `--mmcif=test.cif --dry-run --one-spoke --no-symmetry`
 - For 3 spokes: `--mmcif=test.cif --dry-run --no-symmetry`
 - For 8 spokes: `--mmcif=test.cif --dry-run`

The 8 spoke model also includes FG repeats, which are read from the
`npc_fg_2018` repository. To produce or visualize these modules, you will need
to get a copy of that repository - do that by running

    `git submodule init && git submodule update`

Note that since paths to external files (such as localization densities) are
relative to the script directory, for now in order to visualize the files you
will need to copy them from the current directory into `../../template` and
open them there.

Note that the scripts in **this** directory generate only very basic mmCIF
files, with only coordinates (no experimental information).
