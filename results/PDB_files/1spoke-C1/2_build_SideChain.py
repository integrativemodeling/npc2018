###################################################################################
## README:  Please run this script on Linux or any platform that is case-sensitive
## README:  Windows and OSX are not case-sensitive - so don't run it on them
###################################################################################
from math import sqrt
from Bio.PDB import *
#from Bio.PDB.PDBParser import PDBParser
import argparse
import string
import os

parser = argparse.ArgumentParser(description='Extract PDB chains')
parser.add_argument('-pdb', action="store", dest="pdbfile", help="the pdb file" )
#parser.add_argument('-n_Chain', action="store", dest="n_Chain", help="number of Chains" )
parser=parser.parse_args()

pdbparser = PDBParser()
structure = pdbparser.get_structure(parser.pdbfile, parser.pdbfile)
model = structure[0]
print(parser.pdbfile)


####################################################
## Extract each chain in the PDB file
####################################################
io = PDBIO()
for chain in structure.get_chains():
    print (chain.get_id())
    io.set_structure(chain)
    pdb_name = 'chain' + chain.get_id() + '.pdb'
    io.save(pdb_name)

    # Reconstruct side chain atoms
    os.system("pulchra " + pdb_name)
    os.remove(pdb_name)


####################################################
## recover chain ID 
####################################################
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb
import fnmatch
log.verbose()
env = environ()

for chain in structure.get_chains():
    pdb_name = 'chain' + chain.get_id() + '.rebuilt.pdb'    
    try:
        mdl = model(env, file=pdb_name)
        mdl.chains[0].name = chain.get_id()
        mdl.write(pdb_name)
    except:
        continue


####################################################
## combine individual chains into a single PDB file
####################################################
outfile_name = (structure.get_id()).replace(".pdb", "_rebuilt_temp.pdb")
print (outfile_name)

outfile=open(outfile_name, "w")
io = PDBIO(1)
for chain in structure.get_chains():
    pdb_name = 'chain' + chain.get_id() + '.rebuilt.pdb'
    structure_chain = pdbparser.get_structure(pdb_name, pdb_name)
    io.set_structure(structure_chain)
    io.save(outfile, write_end=True)
    os.remove(pdb_name)
outfile.close()
os.system('grep -vwE "(TER|END|ENDMDL|MODEL      0)" ' + outfile_name + " > " + (structure.get_id()).replace(".pdb", "_rebuilt.pdb"))
os.remove(outfile_name)


"""
io.set_structure(structure)
for i in range(0, int(parser.n_Chain)):
    if (i < 26):
        chainID = string.uppercase[i]
    elif (i < 52):
        chainID = string.lowercase[i-26]
    elif (i < 62):
        chainID = str(i-52)
    else:
        continue
    print (chainID, type(chainID))

    class chain_Select(Select):
        def accept_chain(self, chain):
            if (chain.get_id() == chainID):
                return 1
            else:
                return 0

    io.save('chain' + chainID + '.pdb', chain_Select())
    del chain_Select

exit(1)
"""
