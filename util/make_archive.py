#!/usr/bin/env python

"""
Simple script to make .zip archives of this repository, suitable for
upload to Zenodo or a similar DOI-providing service. Unlike Zenodo's
own GitHub integration (or 'git archive') the repository is split into
several smaller archives, so that users don't have to download enormous
archives containing trajectories just to get a small file like a Python
script or alignment.
"""

from __future__ import print_function, division
import sys
import os
import hashlib
import shutil
import subprocess

# Put larger directories (keys) in uniquely-named zipfiles (values)
ARCHIVES = {
  'validation': 'validation',
  'template': 'template',
  'input_data_files/SAXS': 'SAXS',
  'data_nic96': 'data_nic96',
  'pom152': 'pom152',
  'results/localization_density_files_MRC/1spoke-C1/': 'MRC-1spoke',
  'results/localization_density_files_MRC/3spokes-C3/': 'MRC-3spoke',
  'results/localization_density_files_MRC/8spokes-C8/': 'MRC-8spoke',
  'data_nup84_2016': 'data_nup84_2016',
  'data_npc/MODELLER': 'data_MODELLER',
  'data_npc/em_gmm_model': 'em_gmm_model',
  'data_npc/generating_em_gmm': 'generating_em_gmm',
  'data_nup82': 'data_nup82'
}

REPO="npc2018"

def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def format_size(b):
    suffixes = ['B', 'KiB', 'MiB', 'GiB']
    for s in suffixes:
        if b < 1024:
            return "%.2f %s" % (b, s)
        b /= 1024
    return "%.2f TiB" % b

class Archiver(object):
    ARCHIVE_DIR = "for_archival"

    def __init__(self, tag):
        if os.path.exists(self.ARCHIVE_DIR):
            raise ValueError("The %s directory already exists - please "
                             "delete it first" % self.ARCHIVE_DIR)
        self.tag = tag
        os.mkdir(self.ARCHIVE_DIR)
        self.topdir = os.path.join(self.ARCHIVE_DIR, "%s-%s" % (REPO, self.tag))

    def get_all_files(self):
        print("Extracting all files from %s at %s" % (REPO, self.tag))
        subprocess.check_call('git archive --format=tar --prefix=util/%s/ %s '
                              '| tar -xf -' % (self.topdir, self.tag),
                              shell=True, cwd='..')

    def zip_subdir(self, subdir, zipname):
        base = os.path.basename(subdir)
        subdir_full = os.path.join(self.topdir, subdir)
        cwd = os.getcwd()
        print("Archiving %s" % subdir)
        outzip_full = os.path.join(cwd, self.ARCHIVE_DIR, zipname)
        subprocess.check_call(['zip', '-r', outzip_full, base],
                              cwd=os.path.join(subdir_full, '..'))
        shutil.rmtree(subdir_full)
        os.mkdir(subdir_full)
        with open(os.path.join(subdir_full, 'README.txt'), 'w') as fh:
            fh.write("""
The files in this directory can be found in the %s file
at the same DOI where this archive is available.
""" % zipname)

    def zip_toplevel(self):
        print("Archiving top level")
        dirname = "%s-%s" % (REPO, self.tag)
        subprocess.check_call(['zip', '-r', dirname + '.zip', dirname],
                              cwd=self.ARCHIVE_DIR)
        shutil.rmtree(self.topdir)

    def summarize(self):
        for fname in sorted(os.listdir(self.ARCHIVE_DIR)):
            fullname = os.path.join(self.ARCHIVE_DIR, fname)
            sz = os.stat(fullname).st_size
            print("%s %-10s %s" % (md5(fullname), format_size(sz), fname))
        print("zip files created in %s. Upload them and then"
              % self.ARCHIVE_DIR)
        print("delete that directory.")


def main():
    if len(sys.argv) != 2:
        print("Usage: %s tag" % sys.argv[0], file=sys.stderr)
        sys.exit(1)
    tag = sys.argv[1]
    a = Archiver(tag)
    a.get_all_files()
    for subdir, zipname in ARCHIVES.items():
        a.zip_subdir(subdir, zipname + '.zip')
    a.zip_toplevel()
    a.summarize()

if __name__ == '__main__':
    main()
