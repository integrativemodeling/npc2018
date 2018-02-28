# Validation against electron microscopy class averages

The final model was validated against electron microscopy class averages
as follows:

1. Class averages were obtained (with pixel size 2.03Å) and deposited
   in the `Images` subdirectory.

2. The second copy of the Nic96 complex from the primary half spoke
   (consisting of components Nic96.2, Nsp1.4, Nup49.2, and Nup57.2) was
   extracted from the final model
   (`../../results/RMF_files/cluster0_47-35_1spoke.rmf3`) and stored in
   the `Model_2B` directory as `Nic96complex_C2.pdb`.

3. Flexible beads were removed from the structure to yield
   `Nic96complex_C2_NoBeads.pdb` (as a large amount of flexibility is seen
   in these regions, they do not fit well against the class averages).

4. The script `Model_2B/EM2D_Filter.py` was used to fit the PDB file against
   each of the class averages (using a resolution of 35Å and 10000 projections).
   This reported the cross correlation coefficient (ccc) for the fit in the
   log file `C1_logs_35.txt` and the transformation between the model and the
   image in the file `Registration-Parameters`.

5. Fits for class averages 6 and 25 (with CCCs of 0.85 and 0.80 respectively)
   were reported in the publication in Extended Data Figure 6, panel K.
