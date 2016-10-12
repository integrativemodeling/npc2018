test -r /Applications/EMAN/init.EMAN.1.9.sh && . /Applications/EMAN/init.EMAN.1.9.sh
test -r /Applications/EMAN2/eman2.bashrc && source /Applications/EMAN2/eman2.bashrc

HDF_FILE=pom152_T4_isac_052715
HDF_CLASSES=24

e2proc2d.py ../$HDF_FILE.hdf ./$HDF_FILE.hdf --split $HDF_CLASSES

: '
for FILE in *.hdf; do
    NEW_FILE=${FILE/.hdf/}
    echo "e2proc2d.py $FILE ${NEW_FILE}.png"
    e2proc2d.py $FILE ${NEW_FILE}.png
done
exit -1
'
#voxel size 2.03

rm -rf *.png; rm -rf *.pgm

for FILE in *.hdf; do
    NEW_FILE=${FILE/.hdf/}
    echo "e2proc2d.py $FILE ${NEW_FILE}a.pgm"
    e2proc2d.py $FILE ${NEW_FILE}a.pgm

    convert ${NEW_FILE}a.pgm -compress none ${NEW_FILE}.pgm 
done

rm -rf *a.pgm

