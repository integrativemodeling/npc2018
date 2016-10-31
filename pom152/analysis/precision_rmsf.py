import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.atom
import glob
import itertools


#####################################################
# Parsing parameter inputs
#####################################################
import argparse

parser = argparse.ArgumentParser(description='Performing calculations of the precision and the rmsf for each of the clusters')
parser.add_argument('-test', action="store", dest="test_mode", help="test_mode")
parser.add_argument('-dir', action="store", dest="root_cluster_directory", help="root_cluster_directory")
parser.add_argument('-mpi', action="store", dest="is_mpi", help="is_mpi")
inputs = parser.parse_args()

# runs on the first 10 structures to test if it runs smoothly
if (inputs.test_mode=="True") or (inputs.test_mode=="true") or (inputs.test_mode=="Yes") or (inputs.test_mode=="yes") :
    inputs.test_mode = True
else:
    inputs.test_mode = False

# specify the cluster directory to be analysed
if inputs.root_cluster_directory==None:
    inputs.root_cluster_directory = "kmeans_500_1"

# is_mpi ?
if (inputs.is_mpi=="True") or (inputs.is_mpi=="true") or (inputs.is_mpi=="Yes") or (inputs.is_mpi=="yes") :
    inputs.is_mpi = True
else:
    inputs.is_mpi = False
print inputs

test_mode = inputs.test_mode
root_cluster_directory = inputs.root_cluster_directory
is_mpi = inputs.is_mpi


#####################################################################
# choose whatever selection for the precision calculation
#####################################################################
selection_dictionary={"pom152":["pom152"],
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


#####################################################################
# don't change anything below
#####################################################################
rmfs=[]
frames=[]
clusterdirectories=glob.glob(root_cluster_directory+'/cluster.*/')

if test_mode:
  # runs on the first 10 structures to test if it runs smoothly
  for clusterdirectory in clusterdirectories:
      rmfs.append(glob.glob(clusterdirectory+'/*.rmf3')[0::3])
      #rmfs.append(glob.glob(clusterdirectory+'/*.rmf3')[0::10])
      frames.append([0]*len(rmfs[-1]))
else:
  for clusterdirectory in clusterdirectories:
      rmfs.append(glob.glob(clusterdirectory+'/*.rmf3')[0::1])
      frames.append([0]*len(rmfs[-1]))

model=IMP.Model()
pr=IMP.pmi.analysis.Precision(model, resolution=1, selection_dictionary=selection_dictionary)
pr.set_precision_style('pairwise_rmsd')

for n in range(len(rmfs)):
    pr.add_structures(zip(rmfs[n],frames[n]),clusterdirectories[n]) #,is_mpi=is_mpi)


for pair in itertools.product(range(len(rmfs)), repeat=2):
    clus1=pair[0]
    clus2=pair[1]
    outfile=root_cluster_directory+"/precision."+str(clus1)+"."+str(clus2)+".out"
    pr.get_precision(clusterdirectories[clus1],
                     clusterdirectories[clus2],
                     outfile=outfile,
                     #is_mpi=is_mpi,
                     skip=1)

for n in range(len(rmfs)):
    outdir=clusterdirectories[n]+"/"
    pr.get_rmsf(clusterdirectories[n],
                outdir=outdir,
                #is_mpi=is_mpi,
                skip=1)
    #pr.get_rmsf(clusterdirectories[n],clusterdirectories[n]+"/",is_mpi=is_mpi,skip=1,set_plot_yaxis_range=(0,100.0))
