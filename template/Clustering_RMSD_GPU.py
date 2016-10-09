from __future__ import print_function
import pyRMSD.RMSDCalculator
from pyRMSD.matrixHandler import MatrixHandler
import time
import resource
import numpy as np
import sys, os, glob

from sklearn.cluster import KMeans

import scipy as sp
from scipy import spatial

import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import os
import sys

from math import sqrt
import IMP.rmf
import RMF

conform = []
name = []
num = 0

for file in glob.glob("all_models.499/*.rmf3"):
   print("loading model: ", file, num)
   name.append(file.split("/")[-1])
   m = IMP.Model()
   inf = RMF.open_rmf_file_read_only(file)
   h = IMP.rmf.create_hierarchies(inf, m)[0]
   particle2s = IMP.core.get_leaves(h)
   IMP.rmf.load_frame(inf, 0)
   partcoord = []
   for p in IMP.core.XYZs(particle2s):
      if "Nic96" in p.get_name():
         partcoord.append(p.get_coordinates())
      elif "Nup192" in p.get_name():
         partcoord.append(p.get_coordinates())
      elif "Nup188" in p.get_name():
         partcoord.append(p.get_coordinates())
      elif "Nup170" in p.get_name():
         partcoord.append(p.get_coordinates())
      elif "Nup157" in p.get_name():
         partcoord.append(p.get_coordinates())
      elif "Nup59" in p.get_name():
         if "pdb" in p.get_name():
            partcoord.append(p.get_coordinates())
      elif "Nup53" in p.get_name():
         if "pdb" in p.get_name():
            partcoord.append(p.get_coordinates())
   
   conform.append(partcoord)
   partcoord = []
   num = num + 1

a = np.array(conform)
print(a.shape)


mHandler = MatrixHandler()
matrix = mHandler.createMatrix(a,"QCP_CUDA_MEM_CALCULATOR")
rmsd_matrix = mHandler.getMatrix()
inner_data = rmsd_matrix.get_data()

mHandler.saveMatrix("ex_computed_rmsd.data")
X = sp.spatial.distance.squareform(inner_data)

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as pl
from scipy.cluster import hierarchy as hrc

dendrogram = hrc.dendrogram(
   hrc.linkage(X),
   color_threshold=7,
   no_labels=True)
leaves_order = dendrogram['leaves']

fig = pl.figure(figsize=(4,4)) 
ax2 = fig.add_subplot(111)
cax = ax2.imshow(X[leaves_order,:][:,leaves_order],interpolation='nearest')

cb = fig.colorbar(cax)
cb.set_label('RMSD [Angstroms]')
ax2.set_xlabel('Model')
ax2.set_ylabel('Model')
figurename="RMSD_Matrix.pdf"
pl.savefig(figurename, dpi=300)
pl.close(fig)

from sklearn.cluster import KMeans

for j in range(1, 5):
   kmeans = KMeans(n_clusters=j)
   kmeans.fit_predict(X)
   f1 = open('identities_%s.txt' % j, 'w')

   for i in range(len(kmeans.labels_)):
      f1.write("%s %s %s\n" % (i, name[i], kmeans.labels_[i]))
      #f1.write("\n %s \n" % leaves_order)   
   f1.close()
conform = []

print("")
print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000)
print("")
