#!/usr/bin/env python

import unittest
import os
import sys
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))

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


if __name__ == '__main__':
    unittest.main()
