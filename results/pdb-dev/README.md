To generate mmCIF files suitable for deposition in PDB-Dev, use the
`../../template/1_modeling_wholeNPC.py` script with the following options:

 - For 1 spoke: `--mmcif=test.cif --dry-run --one-spoke --no-symmetry`
 - For 3 spokes: `--mmcif=test.cif --dry-run --no-symmetry`
 - For 8 spokes: `--mmcif=test.cif --dry-run`

The 8 spoke model also includes FG repeats, which are read from the
`npc_fg_2018` subdirectory - clone the `https://github.com/salilab/npc_fg_2018/`
repository into this subdirectory first.
