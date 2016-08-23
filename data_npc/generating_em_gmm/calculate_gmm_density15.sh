#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -r n
#$ -j y
#$ -l mem_free=24G
#$ -l arch=linux-x64
##$ -l netapp=5G,scratch=5G
##$ -l netappsali=5G
##$ -l scrapp=500G
##$ -l scrapp2=500G
#$ -l h_rt=335:59:59
#$ -R y
#$ -V
#$ -q lab.q
#$ -l hostname="i*"                 #-- anything based on Intel cores
##$ -l hostname="!opt*"             #-- anything but opt*
##$ -m e                            #-- uncomment to get email when the job finishes
##$ -pe ompi 4
#$ -t 1
#$ -N GMM200
#########################################

#mapname=SJ_SamplingBoundary
mapname=SJ_Pom152
#mapname=SJ_outer_ring
#ngaussians=1
voxelsize=6
THRESHOLD=0.0012
NUM_ITER=5000
NUM_SAMPLES=10000000

# write hostname and starting time
hostname
date

for (( ngaussians=15; ngaussians<=15; ngaussians++ ))
do
    setup_environment.sh python ~/imp_git/imp/modules/isd_emxl/utility/create_gmm.py $mapname.mrc $ngaussians $mapname.gmm.$ngaussians.txt -m $mapname.gmm.$ngaussians.mrc -a $voxelsize -i $NUM_ITER -n $NUM_SAMPLES -s $THRESHOLD

    hostname
    date
done

# done
hostname
date
