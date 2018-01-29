#MODELLER=mod9.15
MODELLER=python
modelname=Dbp5_3rrm_A_91_482
chain=A


respos=327
resname=LEU
$MODELLER mutate_model.py $modelname $respos $resname $chain
mv ${modelname}${resname}${respos}.pdb ${modelname}.pdb
rm -rf mutate_model.log

