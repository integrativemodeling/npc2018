#!/usr/bin/env python

import unittest
import os
import sys
import subprocess

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
        # Check size of output file
        with open("npc-1spoke.cif") as fh:
            wcl = len(fh.readlines())
        self.assertEqual(wcl, 152780)

    def test_mmcif_3spoke(self):
        """Test generation of 3 spoke mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'template'))
        if os.path.exists("npc-3spoke.cif"):
            os.unlink("npc-3spoke.cif")
        p = subprocess.check_call(
                ["python", "1_modeling_wholeNPC.py", "--mmcif=npc-3spoke.cif",
                 "--dry-run", "--no-symmetry"],
                env=get_mock_env())
        # Check size of output file
        with open("npc-3spoke.cif") as fh:
            wcl = len(fh.readlines())
        self.assertEqual(wcl, 389389)

    def test_mmcif_8spoke(self):
        """Test generation of 8 spoke mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'template'))
        if os.path.exists("npc-8spoke.cif"):
            os.unlink("npc-8spoke.cif")
        p = subprocess.check_call(
                ["python", "1_modeling_wholeNPC.py", "--mmcif=npc-8spoke.cif",
                 "--dry-run", "--one-spoke"],
                env=get_mock_env())
        # Check size of output file
        with open("npc-8spoke.cif") as fh:
            wcl = len(fh.readlines())
        self.assertEqual(wcl, 785763)


if __name__ == '__main__':
    unittest.main()
