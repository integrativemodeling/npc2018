m_start=0
m_end=19
FILE="Nup59"

for ((i=$m_start; i<=$m_end; i++)); do
    cd ~/npc/data_npc/em_gmm_model
    rm -rf ${i}_$FILE.mrc
    rm -rf ${i}_$FILE.txt

    cd ~/npc
    ./job_test4.sh

    cd ~/npc/data_npc/em_gmm_model
    mv $FILE.mrc ${i}_$FILE.mrc
    mv $FILE.txt ${i}_$FILE.txt
done
