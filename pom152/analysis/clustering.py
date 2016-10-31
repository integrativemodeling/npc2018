import IMP
import IMP.pmi
import IMP.pmi.macros
import sys

#####################################################
# Parsing parameter inputs
#####################################################
import argparse

parser = argparse.ArgumentParser(description='generate clusters of the RMF solutions')
parser.add_argument('-mpi', action="store", dest="is_mpi", help="mpi enabled")
parser.add_argument('-preload', action="store", dest="load_distance_matrix_file", help="skip the matrix calcuklation and read the precalculated matrix")
parser.add_argument('-nmods', action="store", dest="nbestscoringmodels", help="number of models to be clustered")
parser.add_argument('-nclusters', action="store", dest="nclusters", help="number of clusters to be used by kmeans algorithm")
parser.add_argument('-prefilter', action="store", dest="prefiltervalue", help="prefilter the models by score")
inputs = parser.parse_args()

# mpi enabled
if (inputs.is_mpi=="True") or (inputs.is_mpi=="true") or (inputs.is_mpi=="Yes") or (inputs.is_mpi=="yes") :
    inputs.is_mpi = True
else:
    inputs.is_mpi = False

# skip the matrix calcuklation and read the precalculated matrix
if (inputs.load_distance_matrix_file=="True") or (inputs.load_distance_matrix_file=="true") or (inputs.load_distance_matrix_file=="Yes") or (inputs.load_distance_matrix_file=="yes") :
    inputs.load_distance_matrix_file = True
else:
    inputs.load_distance_matrix_file = False

# number of models to be clustered
if inputs.nbestscoringmodels==None:
    inputs.nbestscoringmodels = 500

# number of clusters to be used by kmeans algorithm
if inputs.nclusters==None:
    inputs.nclusters = 5

# prefilter the models by score
if inputs.prefiltervalue==None:
    inputs.prefiltervalue = 565.0

print inputs

is_mpi = inputs.is_mpi                                          # mpi enabled
load_distance_matrix_file = inputs.load_distance_matrix_file    # skip the matrix calcuklation and read the precalculated matrix
nbestscoringmodels = int(inputs.nbestscoringmodels)             # number of models to be clustered
nclusters = int(inputs.nclusters)                               # number of clusters to be used by kmeans algorithm
prefiltervalue = float(inputs.prefiltervalue)                   # prefilter the models by score


#####################################################
# initialize the macro
#####################################################
import macros_pom152

model=IMP.Model()

mc=macros_pom152.AnalysisReplicaExchange0(model,
                                        stat_file_name_suffix="stat",     # don't change
                                        #merge_directories=["../modeling1"], # change this list splitting the runs or adding new runs
                                        merge_directories=["../modeling1",
                                                           "../modeling2",
                                                           "../modeling3",
                                                           "../modeling4",
                                                           "../modeling5",
                                                           "../modeling6",
                                                           "../modeling7",
                                                           "../modeling8",
                                                           "../modeling9",
                                                           "../modeling10"],
                                        global_output_directory="output/")

# fields that have to be extracted for the stat file

feature_list=[
              #"ISDCrossLinkMS_Distance_intrarb",
              #"ISDCrossLinkMS_Distance_interrb",
              #"ISDCrossLinkMS_Data_Score",
              #"ISDCrossLinkMS_Psi",
              #"ISDCrossLinkMS_Sigma"
              "GaussianEMRestraint_GaussianEMRestraint_Relion",
              "GaussianEMRestraint_Sampling_Boundary",
              "SimplifiedModel_Linker_Score_None",
              "ExcludedVolumeSphere_None",
              "PlaneDihedralRestraint_pom152_pdr",
              "DistanceRestraint",
             ]

# Dictionary of densities to be calculated
# the key is the name of the file and the value if the selection
# example:
#              {"med17-CTD":[(200,300,"med17")],"med17-CTD.med14":[(200,300,"med17"),"med14"]   }


reduced_density_dict={"pom152":["pom152"],
                      "pom152_NTD":[(1,378,"pom152")],
                      "pom152_1":[(379,472,"pom152")],
                      "pom152_2":[(520,611,"pom152")],
                      "pom152_3":[(616,714,"pom152")],
                      "pom152_4":[(722,818,"pom152")],
                      "pom152_5":[(824,918,"pom152")],
                      "pom152_6":[(931,1026,"pom152")],
                      "pom152_7":[(1036,1141,"pom152")],
                      "pom152_8":[(1150,1229,"pom152")],
                      "pom152_9":[(1244,1337,"pom152")]}


# list of component names needed to calculate the RMSD for the clustering

components_names={"pom152":"pom152"
                  #"EloB":("EloB"),
                 }


mc.clustering("SimplifiedModel_Total_Score_None",  # don't change, field where to find the score
              "rmf_file",                          # don't change, field where to find the path for the rmf_file
              "rmf_frame_index",                   # don't change, field for the frame index
              prefiltervalue=prefiltervalue,               # prefilter the models by score
              number_of_best_scoring_models=nbestscoringmodels,   # number of models to be clustered
              alignment_components=None,           # don't change, (list of proteins you want to use for structural alignment
              #alignment_components=components_names,         # don't change, (list of proteins you want to use for structural alignment
              rmsd_calculation_components=components_names,  # list of proteins used to calculated the rmsd
              distance_matrix_file="distance.rawmatrix.pkl", # save the distance matrix
              outputdir="kmeans_"+str(nbestscoringmodels)+"_"+str(nclusters)+"/",  # directory name for the clustering
              feature_keys=feature_list,                     # extract these fields from the stat file
              load_distance_matrix_file=load_distance_matrix_file,                # skip the matrix calcuklation and read the precalculated matrix
              skip_clustering=False,                         # skip clustering
              display_plot=True,                            # display the heat map plot of the distance matrix
              exit_after_display=False,                      # exit after having displayed the distance matrix plot
              get_every=1,                                   # skip structures for faster computation
              #is_mpi=is_mpi,                                 # mpi enabled
              number_of_clusters=nclusters,                  # number of clusters to be used by kmeans algorithm
              voxel_size=5.0,                                # voxel size of the mrc files
              density_custom_ranges=reduced_density_dict)    # setup the list of densities to be calculated
