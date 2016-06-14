#!/usr/bin/env python

from __future__ import print_function
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.em2d
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.proteomics
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import IMP.npc
import IMP.npc.npc_restraints
import os
import sys

from math import sqrt
import IMP.em2d as em2d
import IMP.rmf
import RMF
import logging
import linecache
import time


####################################################
#  Remove FG nups in the input rmf file
####################################################
m = IMP.Model()
rh = RMF.open_rmf_file_read_only("./B_8_1-95.rmf3")
h = IMP.rmf.create_hierarchies(rh, m)
ps = IMP.atom.get_leaves(h[0])
print (len(ps))
for component in h[0].get_children():
    print ("component = ", component)
    
    if component.get_name()[0:6] == 'Nup159':
        for rep in component.get_children():
            print ("rep = ", rep)
            """
            # Remove Gaussian beads
            if rep.get_name() == 'Densities':
                IMP.atom.destroy(rep)
            """
            for p in rep.get_children():
                self_names=(p.get_name()).replace("-","_").split("_")
                if (len(self_names) > 2) :
                    residue_no = int(self_names[2])
                else :
                    residue_no = int(self_names[0])
                if (residue_no < 1082) :        #if (residue_no < 1117) :
                    print ("p = ", p, "self_names = ", self_names, "residue_no = ", residue_no)
                    IMP.atom.destroy(p)
    
    if component.get_name()[0:4] == 'Nsp1':
        for rep in component.get_children():
            print ("rep = ", rep)
            for p in rep.get_children():
                self_names=(p.get_name()).replace("-","_").split("_")
                if (len(self_names) > 2) :
                    residue_no = int(self_names[2])
                else :
                    residue_no = int(self_names[0])
                if (residue_no < 601) :        #if (residue_no < 637) :
                    print ("p = ", p, "self_names = ", self_names, "residue_no = ", residue_no)
                    IMP.atom.destroy(p)

    if component.get_name()[0:6] == 'Nup116':
        for rep in component.get_children():
            print ("rep = ", rep)
            for p in rep.get_children():
                self_names=(p.get_name()).replace("-","_").split("_")
                if (len(self_names) > 2) :
                    residue_no = int(self_names[2])
                else :
                    residue_no = int(self_names[0])
                if (residue_no < 751) :
                    print ("p = ", p, "self_names = ", self_names, "residue_no = ", residue_no)
                    IMP.atom.destroy(p)

ps = IMP.atom.get_leaves(h[0])
print (len(ps))
rh = RMF.create_rmf_file("./B_8_1-95_FGtruncated.rmf3")
IMP.rmf.add_hierarchies(rh, h)
IMP.rmf.save_frame(rh)
