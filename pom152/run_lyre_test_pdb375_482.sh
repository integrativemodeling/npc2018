#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -r n
#$ -j y
##$ -l mem_free=4G
#$ -l arch=linux-x64
##$ -l netapp=5G,scratch=5G
##$ -l netappsali=5G
##$ -l scrapp=500G
##$ -l scrapp2=500G
#$ -l h_rt=335:59:59
#$ -R y
#$ -V
##$ -q lab.q
##$ -l hostname="i*"                 #-- anything based on Intel cores
##$ -l hostname="!opt*"			    #-- anything but opt*
##$ -m e                            #-- uncomment to get email when the job finishes
#$ -pe ompi 2
#$ -t 1
##$ -t 1-10                        #-- specify the number of tasks
#$ -N MVV_1-10
#########################################

#lyre usage : nohup ./job_test.sh 20000 output > job_test.log &
NSLOTS=8    ## Should be an "EVEN number" or 1
#NSLOTS=1    ## Should be an "EVEN number" or 1
SGE_TASK_ID=1

# load MPI modules
#module load openmpi-1.6-nodlopen
#module load sali-libraries
#mpirun -V

export IMP=setup_environment.sh
MODELING_SCRIPT=modeling_pdb375_482.py
SAXS_FILE=SAXS_26996_merged_q27.dat
XL_FILE=XL.csv
RMF_FILE=../data/rmfs/0.rmf3
RMF_FRAME=0
EM2D_FILE=../data/em2d_WT/no_bg.0.pgm
EM2D_WEIGHT=1000.0
#XL_FILE=XL_MVV_061515_confirmed.csv
#export PYTHONPATH=$PYTHONPATH:/netapp/sali/etjioe/IMP/local_modules/lib64/python/:/netapp/sali/etjioe/IMP/local_modules/lib/python/

for (( SGE_TASK_ID = 1; SGE_TASK_ID <= 10; SGE_TASK_ID++ )); do
    # Parameters
    if [ -z $1 ]; then
        REPEAT="10000"
    else
        REPEAT="$1"
    fi
    echo "number of REPEATs = $REPEAT"

    if [ -z $2 ]; then
        OUTPUT="output"
    else
        OUTPUT="$2"
    fi
    echo "OUTPUT foler = $OUTPUT"

    echo "SGE_TASK_ID = $SGE_TASK_ID"
    echo "JOB_ID = $JOB_ID"
    echo "NSLOTS = $NSLOTS"

    # write hostname and starting time
    hostname
    date

    let "SLEEP_TIME=$SGE_TASK_ID*2"
    #sleep $SLEEP_TIME

    PWD_PARENT=$(pwd)

    i=$(expr $SGE_TASK_ID)
    DIR=modeling$i

    rm -rf $DIR
    if [ ! -d $DIR ]; then
        mkdir $DIR
        cp -pr template/$MODELING_SCRIPT $DIR
        cp -pr template/representation_pom152.py $DIR
        #cp -pr template/topology.txt $DIR
    fi
    cd $DIR

    PWD=$(pwd)
    echo $PWD_PARENT : $PWD

    if [ $PWD_PARENT != $PWD ]; then
        # run the job
        #mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -sym False -r $REPEAT -out $OUTPUT -refine True -w 50.0 -x ../data/$XL_FILE
        
        #echo "mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT"
        #mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT

        echo "mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -rmf $RMF_FILE -rmf_n $RMF_FRAME"
        mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -rmf $RMF_FILE -rmf_n $RMF_FRAME

        #echo "mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -em2d $EM2D_FILE -weight $EM2D_WEIGHT"
        #mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -em2d $EM2D_FILE -weight $EM2D_WEIGHT

        cd ..
    fi

    hostname
    date
    continue

    #cp -pr data/$SAXS_FILE $DIR/$OUTPUT/pdbs

    cd $DIR/$OUTPUT/pdbs
    ../../../template/analysis_extract_chains2.sh
    #$IMP python ~/bin/MSE.py

    #rm -rf foxs.log
    rm -rf foxs_p.log
    for FILE in *r.pdb
    do
        #foxs $FILE $SAXS_FILE >> foxs.log
        foxs $FILE -p >> foxs_p.log
    done
    #gnuplot *.plt


    # done
    hostname
    date
    : '
    echo "copying SAXS dat files in $DIR for MES ..."
    for FILES in model.*r_SAXS_b2gp1a_490_559.dat; do
        NEW_FILES=${FILES/model./}
        NEW_FILES=${NEW_FILES/_SAXS_b2gp1a_490_559/}
        #echo "cp -pr $FILES ../../../mes/${DIR}_${NEW_FILES}"
        cp -pr $FILES ../../../mes/dat/${DIR}_${NEW_FILES}
    done
    '
    echo "copying SAXS dat files in $DIR for Multi_Foxs ..."
    for FILES in model.*r.pdb.dat; do
        NEW_FILES=${FILES/model./}
        NEW_FILES=${NEW_FILES/.pdb/}
        #echo "cp -pr $FILES ../../../mes/${DIR}_${NEW_FILES}"
        cp -pr $FILES ../../../multi_foxs/dat/${DIR}_${NEW_FILES}
    done

    echo "copying PDB files in $DIR for MES ..."
    for FILES in model.*r.pdb; do
        NEW_FILES=${FILES/model./}
        #echo "cp -pr $FILES ../../../pdb/${DIR}_${NEW_FILES}"
        cp -pr $FILES ../../../pdb/${DIR}_${NEW_FILES}
    done

    cd ../../../pdb
    echo "tar czf $DIR.tar.gz ${DIR}_*.pdb"
    tar czf $DIR.tar.gz ${DIR}_*.pdb
    rm -rf ${DIR}_*.pdb


    cd ..
    echo "tar czf $DIR.tar.gz $DIR"
    tar czf $DIR.tar.gz $DIR


    echo "cleaning up $DIR..."
    rm -rf $DIR

    echo ""
done

# done
hostname
date

#cd $OUTPUT
#process_output.py -f stat.0.out -n 1

#process_output.py -f stat.0.out -p

### searching for correlationo of the pdf file with rmf file / output log file
#process_output.py -f stat.5.out -s SimplifiedModel_Total_Score_None ISDCrossLinkMS_Data_Score_DSS rmf_file rmf_frame_index

#process_output.py -f stat.5.out -s SimplifiedModel_Total_Score_None rmf_file rmf_frame_index | grep 103.147
