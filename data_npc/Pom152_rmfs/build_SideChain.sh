MODELLER="python"
#MODELLER="mod9.15"

for FILES in *.pdb
do
    ## Extract chains separately, build side chain atoms, then combine each chain into a single PDB file
    python ./build_SideChain.py -pdb $FILES
done
