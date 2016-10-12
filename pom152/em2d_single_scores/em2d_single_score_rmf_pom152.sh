#rmf_slice ../back_up_pdb375_482/modeling1/output/rmfs/3.rmf3 test.rmf3 -f 1800 -s 1000000

setup_environment.sh python ./em2d_single_score_rmf_pom152.py -rmf ./test.rmf3 -rmf_n 0 -em2d ./pom152_wt_resized_classes.0.pgm -weight 10000
#setup_environment.sh python ./em2d_single_score_rmf_pom152.py -rmf ./test.rmf3 -rmf_n 0 -em2d ../data/em2d/pom152_wt_resized_classes.0.spi -weight 10000

#python ./em2d_single_scores/process_for_em.py -pdb ./test_pdb_writing.pdb -out temp_test_pdb_writing.pdb

#grep -vwE "(END|ENDMDL|MODEL      0)" temp_test_pdb_writing.pdb > tr_test_pdb_writing.pdb
#rm -rf chain*.pdb
#rm -rf temp_test_pdb_writing.pdb
