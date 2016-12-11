#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -r n
#$ -j y
#$ -l mem_free=4G
#$ -l arch=linux-x64
##$ -l netapp=5G,scratch=5G
##$ -l netappsali=5G
##$ -l scrapp=500G
##$ -l scrapp2=500G
#$ -l h_rt=335:59:59
#$ -R y
#$ -V
##$ -q lab.q
#$ -l hostname="i*"                 #-- anything based on Intel cores
##$ -l hostname="!opt*"			    #-- anything but opt*
##$ -m e                            #-- uncomment to get email when the job finishes
#$ -pe ompi 6
##$ -t 1
#$ -t 1-20                        #-- specify the number of tasks
#$ -N npc_IR
#########################################

#: '#lyre usage : nohup ./job_test.sh 20000 output > job_test.log &
#NSLOTS=20   ## Should be an "EVEN number" or 1
NSLOTS=4    ## Should be an "EVEN number" or 1
SGE_TASK_ID=9191
#'
# load MPI modules
#module load openmpi-1.6-nodlopen
#module load sali-libraries
#mpirun -V

export IMP=setup_environment.sh
MODELING_SCRIPT=modeling_MLPs.py
SAXS_FILE=SAXS.dat
XL_FILE=XL.csv
RMF_FRAME=0
EM2D_FILE=../data/em2d/2.pgm
EM2D_WEIGHT=10000.0


# Parameters
if [ -z $1 ]; then
    #REPEAT="10000"
    REPEAT="1000"
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

for (( INDEX=2100; INDEX <= 2120; INDEX++ )); do
    # write hostname and starting time
    echo "INDEX = $INDEX"
    hostname
    date

    DIR=modeling${INDEX}
    RMF_FILE=../prefilter/${INDEX}kmeans_10_1/all_models.9/0_1spoke.rmf3

    rm -rf $DIR
    if [ ! -d $DIR ]; then
        mkdir $DIR
        cp -pr template/$MODELING_SCRIPT $DIR
        #cp -pr template/em2d_nup82.py $DIR
        #cp -pr template/representation_nup82.py $DIR
        #cp -pr template/crosslinking_nup82.py $DIR
    fi
    sleep 1
    cd $DIR

    PWD=$(pwd)
    echo $PWD_PARENT : $PWD

    if [ $PWD_PARENT != $PWD ]; then
        # run the job
        #echo "mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT"
        #mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT

        #mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -sym False -r $REPEAT -out $OUTPUT -refine True -w 50.0 -x ../data/$XL_FILE
        #echo "mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -em2d $EM2D_FILE -weight $EM2D_WEIGHT"
        #mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -em2d $EM2D_FILE -weight $EM2D_WEIGHT
        echo "mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -rmf $RMF_FILE -rmf_n $RMF_FRAME"
        mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -rmf $RMF_FILE -rmf_n $RMF_FRAME > $DIR.log
        cd ..
    fi

    # done
    hostname
    date
done


for (( INDEX=2200; INDEX <= 2220; INDEX++ )); do
    # write hostname and starting time
    echo "INDEX = $INDEX"
    hostname
    date

    DIR=modeling${INDEX}
    RMF_FILE=../prefilter/${INDEX}kmeans_10_1/all_models.9/0_1spoke.rmf3

    rm -rf $DIR
    if [ ! -d $DIR ]; then
        mkdir $DIR
        cp -pr template/$MODELING_SCRIPT $DIR
        #cp -pr template/em2d_nup82.py $DIR
        #cp -pr template/representation_nup82.py $DIR
        #cp -pr template/crosslinking_nup82.py $DIR
    fi
    sleep 1
    cd $DIR

    PWD=$(pwd)
    echo $PWD_PARENT : $PWD

    if [ $PWD_PARENT != $PWD ]; then
        # run the job
        #echo "mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT"
        #mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT

        #mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -sym False -r $REPEAT -out $OUTPUT -refine True -w 50.0 -x ../data/$XL_FILE
        #echo "mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -em2d $EM2D_FILE -weight $EM2D_WEIGHT"
        #mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -em2d $EM2D_FILE -weight $EM2D_WEIGHT
        echo "mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -rmf $RMF_FILE -rmf_n $RMF_FRAME"
        mpirun -np $NSLOTS $IMP python ./$MODELING_SCRIPT -r $REPEAT -out $OUTPUT -rmf $RMF_FILE -rmf_n $RMF_FRAME > $DIR.log
        cd ..
    fi

    # done
    hostname
    date
done
