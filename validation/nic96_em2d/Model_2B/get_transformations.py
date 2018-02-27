#!/usr/bin/env python

"""Convert the reported EM2D registrations (file "Registration-Parameters")
   for each image into a set of transformations, each of which transforms
   the model coordinates so that they superpose on the corresponding class
   average (assuming that the class average image lies in the xy plane
   with its origin at (0,0,0)). Print each transformation in a form suitable
   for including in an IHM mmCIF file."""

from __future__ import print_function
import IMP.algebra
import IMP.atom
import math

def get_centroid(fname):
    """Return the coordinates of the centroid of the given PDB file"""
    m = IMP.Model()
    mh = IMP.atom.read_pdb(fname, m, IMP.atom.CAlphaPDBSelector())
    return IMP.core.get_centroid(IMP.atom.get_leaves(mh))

def get_registrations(fname):
    """Yield information from the given registration parameters file"""
    with open(fname) as fh:
        for line in fh:
            if line.startswith('#') or '|' not in line:
                continue
            spl = line.split('|')
            # Image #, quaternion plus x and y shift
            yield(int(spl[0]), [float(spl[i]) for i in range(5,9)],
                  [float(spl[i]) for i in range(9,11)])

def get_transformations(registrations, centroid, pixel_size, image_size):
    """Yield a set of transformations from the given registration file"""
    # See do_project_particles() in IMP's modules/em2d/src/project.cpp for
    # the underlying algorithm used here.

    # Reported rotation quaternion assumes both the model centroid and the
    # center of the class average are at the origin, so we need to translate
    # accordingly
    to_centroid = IMP.algebra.Transformation3D(IMP.algebra.Vector3D(-centroid))
    image_origin = pixel_size * image_size * 0.5

    # em2d code uses the unusual convention of x=row and y=column, so correct
    # for this with a 90 degree rotation about the origin
    xy_rot = IMP.algebra.get_rotation_about_axis(
                           IMP.algebra.Vector3D(0., 0., 1.), -0.5 * math.pi)
    cleanup = IMP.algebra.Transformation3D(xy_rot,
                          IMP.algebra.Vector3D(image_origin, image_origin, 0.))

    for num, quat, shift in get_registrations(registrations):
        rot = IMP.algebra.Rotation3D(quat)
        # Reported shift is in pixels, so convert to angstroms
        tnsl = IMP.algebra.Vector3D(pixel_size * shift[0],
                                    pixel_size * shift[1], 0.)
        trans = cleanup * IMP.algebra.Transformation3D(rot, tnsl) * to_centroid
        yield num, trans

def main():
    # Experiment parameters
    pixel_size = 2.03
    image_size = 128

    centroid = get_centroid('Nic96complex_C2_NoBeads.pdb')
    for num, trans in get_transformations('Registration-Parameters', centroid,
                                          pixel_size, image_size):
        rot = trans.get_rotation()
        shift = trans.get_translation()
        rm = [["%.6f" % e for e in rot.get_rotation_matrix_row(i)]
              for i in range(3)]
        print("%d %s %s %s %s %s %s %s %s %s %.3f %.3f %.3f"
              % (num, rm[0][0], rm[1][0], rm[2][0], rm[0][1], rm[1][1],
                 rm[2][1], rm[0][2], rm[1][2], rm[2][2], shift[0],
                 shift[1], shift[2]))

if __name__ == '__main__':
    main()
