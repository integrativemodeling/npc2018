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

mapname=SJ_SamplingBoundary
ngaussians=10
voxelsize=6
THRESHOLD=0.0012
NUM_ITER=10000
NUM_SAMPLES=20000000

# write hostname and starting time
hostname
date

#setup_environment.sh python ~/imp_git/imp/modules/isd_emxl/utility/create_gmm.py $mapname.mrc $ngaussians $mapname.gmm.$ngaussians.txt -m $mapname.gmm.$ngaussians.mrc -a $voxelsize -i 1000 -n 50000000 -s $THRESHOLD
setup_environment.sh python ~/imp_git/imp/modules/isd_emxl/utility/create_gmm.py $mapname.mrc $ngaussians $mapname.gmm.$ngaussians.txt -m $mapname.gmm.$ngaussians.mrc -a $voxelsize -i $NUM_ITER -n $NUM_SAMPLES -s $THRESHOLD
#setup_environment.sh python ~/imp_git/imp/modules/isd_emxl/utility/create_gmm.py $mapname.mrc $ngaussians $mapname.gmm.$ngaussians.txt -s 0.00079 -m $mapname.gmm.$ngaussians.mrc  -a 3.23 -i 200

# done
hostname
date
