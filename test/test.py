#!/usr/bin/env python

import unittest
import os
import sys
import subprocess
import ihm.reader

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))

def get_mock_env():
    # Potentially override methods that need network access
    env = os.environ.copy()
    env['PYTHONPATH'] = os.path.join(TOPDIR, 'test', 'mock') \
                        + ':' + env.get('PYTHONPATH', '')
    return env

class Tests(unittest.TestCase):
    def run_modeller_script(self, script_dir, script_name, model_name, resrng,
                            chains=['A']):
        """Run a Modeller script and test the output model"""
        os.chdir(os.path.join(TOPDIR, script_dir))
        # Run script
        p = subprocess.check_call(["python", script_name, "--test"])
        # Make sure PDB was produced with the requested residue range
        with open('%s.B99990001.pdb' % model_name) as fh:
            pdb_lines = [x for x in fh.readlines() if x.startswith('ATOM')]
        rng = (int(pdb_lines[0][22:26]), int(pdb_lines[-1][22:26]))
        actual_chains = frozenset(line[21] for line in pdb_lines)
        self.assertEqual(rng, resrng)
        self.assertEqual(actual_chains, frozenset(chains))

    def test_nup100(self):
        """Test generation of comparative model for Nup100"""
        self.run_modeller_script(os.path.join('data_npc', 'MODELLER', 'Nup100'),
                                 'Nup100.py', 'Nup100', (803, 958), ['A', 'B'])

    def test_nic96_combined(self):
        """Test generation of combined comparative model for all of Nic96"""
        self.run_modeller_script(os.path.join('data_nic96', 'MODELLER',
                                              'combined'),
                                 'Build-combined.py', 'Nic96_complex',
                                 (20, 540), ['A', 'B', 'C', 'D'])

    def test_nic96_ntd(self):
        """Test generation of comparative model for Nic96 NTD"""
        self.run_modeller_script(os.path.join('data_nic96', 'MODELLER',
                                              'Nic96_NTD'),
                                 'Build-Nic96N.py', 'Nic96_N',
                                 (20, 56), ['F'])

    def test_nic96_nsp1(self):
        """Test generation of comparative model for Nic96 Nsp1"""
        self.run_modeller_script(os.path.join('data_nic96', 'MODELLER',
                                              'Nsp1'),
                                 'Build-Nsp1.py', 'NS3',
                                 (788, 823), ['C'])

    def test_nic96_nup49(self):
        """Test generation of comparative model for Nic96 Nup49"""
        self.run_modeller_script(os.path.join('data_nic96', 'MODELLER',
                                              'Nup49'),
                                 'Build-N49.py', 'N49',
                                 (270, 472), ['D'])

    def test_nic96_nup57(self):
        """Test generation of comparative model for Nic96 Nup57"""
        self.run_modeller_script(os.path.join('data_nic96', 'MODELLER',
                                              'Nup57'),
                                 'Build-N57.py', 'N57',
                                 (287, 541), ['E'])

    def test_mmcif_1spoke(self):
        """Test generation of 1 spoke mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'template'))
        if os.path.exists("npc-1spoke.cif"):
            os.unlink("npc-1spoke.cif")
        p = subprocess.check_call(
                ["python", "1_modeling_wholeNPC.py", "--mmcif=npc-1spoke.cif",
                 "--dry-run", "--one-spoke", "--no-symmetry"],
                env=get_mock_env())
        # Check output file
        self._check_mmcif_file('npc-1spoke.cif', 1)

    def test_mmcif_3spoke(self):
        """Test generation of 3 spoke mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'template'))
        if os.path.exists("npc-3spoke.cif"):
            os.unlink("npc-3spoke.cif")
        p = subprocess.check_call(
                ["python", "1_modeling_wholeNPC.py", "--mmcif=npc-3spoke.cif",
                 "--dry-run", "--no-symmetry"],
                env=get_mock_env())
        # Check output file
        self._check_mmcif_file('npc-3spoke.cif', 3)

    def test_mmcif_8spoke(self):
        """Test generation of 8 spoke mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'template'))
        if os.path.exists("npc-8spoke.cif"):
            os.unlink("npc-8spoke.cif")
        p = subprocess.check_call(
                ["python", "1_modeling_wholeNPC.py", "--mmcif=npc-8spoke.cif",
                 "--dry-run", "--one-spoke"],
                env=get_mock_env())
        # Check output file
        self._check_mmcif_file('npc-8spoke.cif', 8)

    def _check_mmcif_file(self, fname, num_spokes):
        with open(fname) as fh:
            s, = ihm.reader.read(fh)
        self.assertEqual(len(s.citations), 1)
        self.assertEqual(s.citations[0].doi, '10.1038/nature26003')
        self.assertEqual(len(s.software), 13)
        if num_spokes == 3:
            self.assertEqual(len(s.orphan_starting_models), 62 * 3)
        else:
            self.assertEqual(len(s.orphan_starting_models), 62)
        # Should be a single model in a single state for 1 or 3 spokes,
        # 2 models for 8 spokes
        self.assertEqual(len(s.state_groups), 1)
        self.assertEqual(len(s.state_groups[0]), 1)
        if num_spokes == 8:
            self.assertEqual(len(s.state_groups[0][0]), 2)
        else:
            self.assertEqual(len(s.state_groups[0][0]), 1)
        model1 = s.state_groups[0][0][0][0]
        self.assertEqual(len(model1._atoms), 0)
        self.assertEqual(len(model1._spheres), 29273 * num_spokes)

        # Check FG model
        if num_spokes == 8:
            model2 = s.state_groups[0][0][1][0]
            self.assertEqual(len(model2._atoms), 0)
            self.assertEqual(len(model2._spheres), 4104)
        # Should be 1 ensemble (cluster) for 1,3 spoke; 2 for the 8-spoke model
        # Full ensembles for 1-spoke scaffold and 8-spoke FG should be present
        if num_spokes == 8:
            self.assertEqual([e.num_models for e in s.ensembles], [5, 1000])
            self.assertNotEqual(s.ensembles[1].file, None)
        else:
            self.assertEqual([e.num_models for e in s.ensembles], [5])
            if num_spokes == 1:
                self.assertNotEqual(s.ensembles[0].file, None)

        # Check localization densities
        n_densities = {1:[68], 3:[204], 8:[544,200]}
        self.assertEqual([len(e.densities) for e in s.ensembles],
                         n_densities[num_spokes])
        self.assertEqual([len(e.sequence) for e in s.entities],
                         [726, 744, 1037, 1157, 712, 349, 297, 92, 713, 1460,
                          823, 839, 472, 541, 1391, 1502, 1655, 1683, 475,
                          528, 655, 299, 1337, 959, 1113, 430, 538, 1317,
                          1076, 539, 1875, 1679])
        self.assertEqual(len(s.asym_units), 69 * num_spokes)
        # 103 restraints
        self.assertEqual(len(s.restraints), 103)
        # 2 crosslink restraints
        xl_rsr = [r for r in s.restraints
                  if isinstance(r, ihm.restraint.CrossLinkRestraint)]
        self.assertEqual(len(xl_rsr), 2)
        xl1, xl2 = xl_rsr
        self.assertEqual(xl1.linker_type, 'DSS')
        self.assertEqual(xl2.linker_type, 'DSS')
        if num_spokes == 1:
            self.assertEqual(len(xl1.experimental_cross_links), 505)
        else:
            self.assertEqual(len(xl1.experimental_cross_links), 508)
        if num_spokes == 3:
            self.assertEqual(len(xl1.cross_links), 970)
        else:
            # todo: check 8-spoke model for these missing 9 XLs
            self.assertEqual(len(xl1.cross_links), 961)
        self.assertEqual(len(xl2.experimental_cross_links), 509)
        self.assertEqual(len(xl2.cross_links), 112)

        # 61 geometric restraints
        geom_rsr = [r for r in s.restraints
                    if isinstance(r, ihm.restraint.GeometricRestraint)]
        self.assertEqual(len(geom_rsr), 61)

        # 37 SAXS restraints
        sas_rsr = [r for r in s.restraints
                   if isinstance(r, ihm.restraint.SASRestraint)]
        self.assertEqual(len(sas_rsr), 37)
        self.assertEqual(sas_rsr[0].segment, False)
        self.assertEqual(sas_rsr[0].fitting_method, 'FoXS')
        self.assertAlmostEqual(sas_rsr[0].radius_of_gyration, 27.9, places=1)

        # 2 EM2D restraints
        em2d_rsr = [r for r in s.restraints
                    if isinstance(r, ihm.restraint.EM2DRestraint)]
        self.assertEqual(len(em2d_rsr), 2)
        self.assertAlmostEqual(em2d_rsr[0].image_resolution, 35.0, places=1)
        self.assertEqual(em2d_rsr[0].number_raw_micrographs, 800)
        self.assertEqual(len(em2d_rsr[0].fits), 1)

        # 1 EM3D restraint
        em3d_rsr = [r for r in s.restraints
                    if isinstance(r, ihm.restraint.EM3DRestraint)]
        self.assertEqual(len(em3d_rsr), 1)
        em3d, = em3d_rsr
        self.assertEqual(em3d.dataset.location.path,
                'em_gmm_model/SJ_cropped_sym8_avg_monomer_final_rotated_'
                'adjusted90.gmm.400.txt')
        self.assertEqual(em3d.dataset.parents[0].location.access_code,
                         'EMD-7321')


if __name__ == '__main__':
    unittest.main()
