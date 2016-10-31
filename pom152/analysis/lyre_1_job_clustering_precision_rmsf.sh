###############################################################################
# setting up a PREFILTER value
###############################################################################
if [ -z $1 ]; then
    PREFILTER="1300.0"
else
    PREFILTER="$1"
fi
echo "PREFILTER = $PREFILTER"

if [ -z $2 ]; then
    NMODS="1000"
else
    NMODS="$2"
fi
echo "NMODS = $NMODS"


###############################################################################
# submit a job_clustering_init job in parallel and write JobID into a temp file
###############################################################################
echo; echo "./lyre_job_clustering_init.sh $PREFILTER $NMODS"
./lyre_job_clustering_init.sh $PREFILTER $NMODS

: '
qsub ./job_clustering_init.sh $PREFILTER $NMODS | while read line
do
  bar=(`echo $line | tr '.' ' '`)
  jid=${bar[2]}
  echo $jid > "job_clustering_init_jobid.txt"
done


#read Job ID from the temp file
read jid < "job_clustering_init_jobid.txt"
rm -rf job_clustering_init_jobid.txt
'


###############################################################################
# submmit the job_clustering jobs (HOLDED) and write JobID into a temp file
###############################################################################
echo; echo "./lyre_job_clustering.sh $PREFILTER $NMODS"
./lyre_job_clustering.sh $PREFILTER $NMODS

: '
qsub -hold_jid $jid ./job_clustering.sh $PREFILTER $NMODS | while read line
do
  bar=(`echo $line | tr '.' ' '`)
  jid=${bar[2]}
  echo $jid > "job_clustering_jobid.txt"
done


#read Job ID from the temp file
read jid < "job_clustering_jobid.txt"
rm -rf job_clustering_jobid.txt
'


###############################################################################
# submmit the job_precision_rmsf jobs (HOLDED)
###############################################################################
echo; echo "./lyre_job_precision_rmsf.sh $NMODS"
#./lyre_job_precision_rmsf.sh $NMODS

nohup ./lyre_job_precision_rmsf1.sh $NMODS > ./lyre_job_precision_rmsf1.log &
nohup ./lyre_job_precision_rmsf2.sh $NMODS > ./lyre_job_precision_rmsf2.log &
nohup ./lyre_job_precision_rmsf3.sh $NMODS > ./lyre_job_precision_rmsf3.log &
nohup ./lyre_job_precision_rmsf4.sh $NMODS > ./lyre_job_precision_rmsf4.log &
nohup ./lyre_job_precision_rmsf5.sh $NMODS > ./lyre_job_precision_rmsf5.log &
