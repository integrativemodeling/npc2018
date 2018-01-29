#MODELLER=mod9.15
MODELLER=python
modelname=Gle1_3rrm_B_244_538
chain=B


respos=337
resname=HIS
$MODELLER mutate_model.py $modelname $respos $resname $chain
mv ${modelname}${resname}${respos}.pdb ${modelname}.pdb
rm -rf mutate_model.log

