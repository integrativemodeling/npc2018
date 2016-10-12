#FILE=model.0.pdb
FILE=0.pdb

class_start=12
class_end=0
class_name=Pom152_WT_ali_classes_120315_2
#class_name=pom152_T4_isac_052715

#n_projection=500
n_projection=100

#n_resolution=5
n_resolution=40
n_pixel_size=2.03
#n_pixel_size=2.93

for ((i=$class_end; i<$class_start; i++)); do
    #PGM+="../data/em2d_T2_truncation_1-936/$class_name.${i}.pgm "
    PGM+="../data/em2d_WT/$class_name.${i}.pgm "
done

cd em2d_single_scores
#DIR=../analysis/kmeans_1000_1/cluster.0
#echo "copying from $DIR/$FILE"
#cp -pr $DIR/$FILE .

#em2d_single_score -r $n_resolution -s $n_pixel_size -n $n_projection $FILE $PGM
#exit -1

for ((i=$class_start; i>=$class_end; i--)); do
    echo ""
    echo "######################"
    echo "## EM 2D class $i   ##"
    echo "######################"

    if [ "$i" != "$class_start" ]; then
        PGM=../data/em2d_WT/$class_name.${i}.pgm
    fi

    #DIR=../modeling1/output/pdbs
    DIR=../analysis/kmeans_1000_1/cluster.0
    echo "copying from $DIR/$FILE"
    cp -pr $DIR/$FILE .
    
    em2d_single_score -r $n_resolution -s $n_pixel_size -n $n_projection -l 1 $FILE $PGM
    
    mv best_projections.pgm $class_name.${i}_best_projections_1.pgm
    mv images.pgm no_bg.${i}_1.pgm
done

