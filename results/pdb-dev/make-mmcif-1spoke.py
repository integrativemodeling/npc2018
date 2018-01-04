# Generate an mmCIF file, for deposition in PDB-Dev

import IMP.mmcif
from rmf_mangle import RMFFrame
  
s = IMP.mmcif.System()
state = IMP.mmcif.State(s)
e = IMP.mmcif.Ensemble(state, name='Cluster 0')
frame = RMFFrame("../RMF_files/cluster0_47-35_1spoke.rmf3", frame=0,
                 name="Cluster representative")
e.add_frame(frame)
s.add_software(name="Integrative Modeling Platform (IMP)",
               classification="integrative model building",
               version=frame._producer,
               description=None, url="https://integrativemodeling.org")

s.write(fname='npc-1spoke.cif')
