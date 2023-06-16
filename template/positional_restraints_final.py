#####################################################
# Restraints setup - Immuno-EM
# Supplementary Table 7. Upper and lower bounds on R-radial restraints of C-terminal bead of nups
# NupType : (min R value, max R value) (in Angstrom)
# Supplementary Table 7. Upper and lower bounds on Z-axial restraints of C-terminal bead of nups
# NupType : (min Z value, max Z value) (in Angstrom)
#####################################################
radial_weight = zaxial_weight = yaxial_weight = 10.0
if (use_Immuno_EM):
    nup_list = [entry[0] for entry in domains if not '@' in entry[0]]
    nup_list_unique = sorted(list(set(nup_list)))   # Make a unique list

    RADIAL = {
        "Ndc1"     : [106, 448],
        "Nic96.1"  : [  0, 633],
        "Nic96.2"  : [  0, 633],
        "Nup120"   : [190, 550],
        "Nup133"   : [190, 550],
        "Nup145c"  : [190, 550],
        "Nup157"   : [106, 448],
        "Nup159.1" : [ 43, 637],
        "Nup159.2" : [ 43, 637],
        "Nup170"   : [106, 448],
        "Nup188"   : [152, 368],
        "Nup192"   : [152, 368],
        "Nup49.1"  : [  0, 633],
        "Nup49.2"  : [  0, 633],
        "Nup53"    : [106, 448],
        "Nup57.1"  : [  0, 633],
        "Nup57.2"  : [  0, 633],
        "Nup59"    : [106, 448],
        "Nup82.1"  : [ 43, 637],
        "Nup82.2"  : [ 43, 637],
        "Nup84"    : [190, 550],
        "Nup85"    : [190, 550],
        "Pom34"    : [106, 448],
        "Pom152"   : [266, 734]
    }
    print("\nXYRadialPositionRestraint !!")
    for protein, r in RADIAL.items():
        if (protein not in nup_list_unique):
            continue
        xyr = IMP.pmi1.restraints.npc.XYRadialPositionRestraint(simo, protein, lower_bound=r[0], upper_bound=r[1], consider_radius=False, sigma=1.0)
        xyr.set_label('Lower_%d_Upper_%d_%s' % (r[0], r[1], protein))
        xyr.set_weight(radial_weight)
        xyr.add_to_model()
        outputobjects.append(xyr)
        print(xyr.get_output())

    ZAXIAL = {
        "Ndc1"     : [-70, 152],
        "Nic96.1"  : [-60, 235],
        "Nic96.2"  : [-60, 235],
        "Nup120"   : [ 38, 240],
        "Nup133"   : [ 38, 240],
        "Nup145c"  : [ 38, 240],
        "Nup157"   : [-70, 152],
        "Nup159.1" : [ 72, 355],
        "Nup159.2" : [ 72, 355],
        "Nup170"   : [-70, 152],
        "Nup188"   : [ 16, 124],
        "Nup192"   : [-12, 132],
        "Nup49.1"  : [-60, 235],
        "Nup49.2"  : [-60, 235],
        "Nup53"    : [-70, 152],
        "Nup57.1"  : [-60, 235],
        "Nup57.2"  : [-60, 235],
        "Nup59"    : [-70, 152],
        "Nup82.1"  : [ 72, 355],
        "Nup82.2"  : [ 72, 355],
        "Nup84"    : [ 38, 240],
        "Nup85"    : [ 38, 240],
        "Pom34"    : [-70, 152],
        "Pom152"   : [-68, 148]
    }
    print("\nZAxialPositionRestraint !!")
    for protein, z in ZAXIAL.items():
        if (protein not in nup_list_unique):
            continue
        zax = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(simo, protein, lower_bound=z[0], upper_bound=z[1], consider_radius=False, sigma=1.0)
        zax.set_label('Lower_%d_Upper_%d_%s' % (z[0], z[1], protein))
        zax.set_weight(zaxial_weight)
        zax.add_to_model()
        outputobjects.append(zax)
        print(zax.get_output())


#####################################################
# Restraints setup - Cytoplasm / Nucleoplasm / Basket
#####################################################
if (is_cytoplasm or is_nucleoplasm or is_basket):
    nup_list = [entry[0] for entry in domains if not '@' in entry[0]]
    nup_list_unique = sorted(list(set(nup_list)))   # Make a unique list

    ZAXIAL = {
        "Nup100.1" : [   0, 400],
        "Nup100.2" : [   0, 400],
        "Nup116.1" : [   0, 400],
        "Nup116.2" : [   0, 400],
        "Nup42"    : [   0, 400],
        "Gle1"     : [   0, 400],
        "Gle2.1"   : [   0, 400],
        "Gle2.2"   : [   0, 400],
        "Nup145.1" : [-400,   0],
        "Nup145.2" : [-400,   0],
        "Nup1"     : [-400,   0],
        "Nup60.1"  : [-400,   0],
        "Nup60.2"  : [-400,   0],
        "Mlp1"     : [-700,-200],
        "Mlp2"     : [-700,-200]
    }
    print("\nZAxialPositionRestraints for Cytoplasm / Nucleoplasm / Basket !!")
    for protein, z in ZAXIAL.items():
        if (protein not in nup_list_unique):
            continue
        # applying Z-position restraints on the C-termini of each protein
        zax = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(simo, protein, lower_bound=z[0], upper_bound=z[1], consider_radius=False, sigma=1.0)
        zax.set_label('Lower_%d_Upper_%d_%s' % (z[0], z[1], protein))
        zax.set_weight(zaxial_weight)
        zax.add_to_model()
        outputobjects.append(zax)
        print(zax.get_output())

    # localizing the Mlp1-Mlp2 hetero-dimer near the outer-ring on the nucleoplasm
    dist_min = -400.0
    dist_max = -200.0
    zax = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(simo, (716,716,"Mlp1"), lower_bound=dist_min, upper_bound=dist_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (dist_min, dist_max, "Mlp1_716"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print(zax.get_output())

    zax = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(simo, (690,690,"Mlp2"), lower_bound=dist_min, upper_bound=dist_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (dist_min, dist_max, "Mlp2_690"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print(zax.get_output())


#####################################################
# Restraints setup - FG anchor restraints
#####################################################
if (use_FG_anchor and is_nic96 and not is_FG):
    nup_list = [entry[0] for entry in domains if not '@' in entry[0]]
    nup_list_unique = sorted(list(set(nup_list)))   # Make a unique list

    RADIAL = {
        "Nsp1.3"   : [180, 225],
        "Nup49.1"  : [180, 225],
        "Nup57.1"  : [180, 225],
        "Nsp1.4"   : [180, 225],
        "Nup49.2"  : [180, 225],
        "Nup57.2"  : [180, 225],
        #"Nup1"     : [180, 350],
        #"Nup60.1"  : [180, 350],
        #"Nup60.2"  : [180, 350]
    }
    print("\nFG_XYRadialPositionRestraint !!")
    radial_weight = 10.0
    for protein, r in RADIAL.items():
        if (protein not in nup_list_unique):
            continue
        if (protein in ['Nup1', 'Nup60.1', 'Nup60.2']):
            xyr = IMP.pmi1.restraints.npc.XYRadialPositionRestraint(simo, protein, lower_bound=r[0], upper_bound=r[1], consider_radius=False, sigma=1.0, term='C')
        else:
            xyr = IMP.pmi1.restraints.npc.XYRadialPositionRestraint(simo, protein, lower_bound=r[0], upper_bound=r[1], consider_radius=False, sigma=1.0, term='N')
        xyr.set_label('Lower_%d_Upper_%d_%s' % (r[0], r[1], protein))
        xyr.set_weight(radial_weight)
        xyr.add_to_model()
        outputobjects.append(xyr)
        print(xyr.get_output())


#####################################################
# Restraints setup
# Distance restraints
#####################################################
dist_min = 3.0
dr_weight = 10.0

# Nup145n - Nup145c
if (is_n84 and is_nucleoplasm):
    dist_max = 10.0
    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(605,605,"Nup145.1"), (1,1,"Nup145c@11"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.set_label("Nup145n-Nup145c@11")
    dr.set_weight(dr_weight)
    dr.add_to_model()
    outputobjects.append(dr)
    print(dr.get_output())

# Nup60 homo-dimer cross-link
if (is_nucleoplasm):
    dist_max = 25.0
    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(151,151,"Nup60.1"), (151,151,"Nup60.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.set_label("Nup60.1_151-Nup60.2_151")
    dr.set_weight(dr_weight)
    dr.add_to_model()
    outputobjects.append(dr)
    print(dr.get_output())

"""
# Nup120 - Nup133 to form the outer ring  (Seo et al, PNAS 2009) ; Not sure if it is real
if (is_n84 and use_neighboring_spokes):
    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(11,11,"Nup133"), (641,641,"Nup120@2"), distancemin=3.0, distancemax=35.0, resolution=1.0)
    dr.set_label("Nup133-Nup120@2")
    dr.set_weight(10.0)
    dr.add_to_model()
    outputobjects.append(dr)
    print(dr.get_output())
"""

if (is_inner_ring):
    dist_max = 15.0
    # connection between NTD and CTD
    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(892,892,"Nup157"), (900,900,"Nup157"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.set_label("Nup157N-C")
    dr.set_weight(dr_weight)
    dr.add_to_model()
    outputobjects.append(dr)
    print(dr.get_output())

    # connection between NTD and CTD
    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(992,992,"Nup170"), (1000,1000,"Nup170"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.set_label("Nup170N-C")
    dr.set_weight(dr_weight)
    dr.add_to_model()
    outputobjects.append(dr)
    print(dr.get_output())

if (is_basket):
    # end-to-end distance of Mlps (230-350 angstrom)
    dist_min = 330.0    #dist_min = 230.0
    dist_max = 450.0    #dist_max = 350.0
    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(500,500,"Mlp1"), (1875,1875,"Mlp1"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.set_label("Mlp1_end-to-end")
    dr.set_weight(dr_weight)
    dr.add_to_model()
    outputobjects.append(dr)
    print(dr.get_output())

    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(500,500,"Mlp2"), (1679,1679,"Mlp2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.set_label("Mlp2_end-to-end")
    dr.set_weight(dr_weight)
    dr.add_to_model()
    outputobjects.append(dr)
    print(dr.get_output())

    # Distal Ring
    dist_min = 1.35*130.0
    dist_max = 1.35*170.0
    xyr = IMP.pmi1.restraints.npc.XYRadialPositionRestraint(simo, (1875,1875,"Mlp1"), lower_bound=dist_min, upper_bound=dist_max, consider_radius=False, sigma=1.0, term='M')
    xyr.set_label('Lower_%d_Upper_%d_%s' % (dist_min, dist_max, "Mlp1_1875"))
    xyr.set_weight(radial_weight)
    xyr.add_to_model()
    outputobjects.append(xyr)
    print(xyr.get_output())

    xyr = IMP.pmi1.restraints.npc.XYRadialPositionRestraint(simo, (1679,1679,"Mlp2"), lower_bound=dist_min, upper_bound=dist_max, consider_radius=False, sigma=1.0, term='M')
    xyr.set_label('Lower_%d_Upper_%d_%s' % (dist_min, dist_max, "Mlp2_1679"))
    xyr.set_weight(radial_weight)
    xyr.add_to_model()
    outputobjects.append(xyr)
    print(xyr.get_output())

    # localizing the Mlp1-Mlp2 hetero-dimer near the x-axis
    dist_min = 0.0
    dist_max = 40.0
    yax = IMP.pmi1.restraints.npc.YAxialPositionRestraint(simo, (1875,1875,"Mlp1"), lower_bound=dist_min, upper_bound=dist_max, consider_radius=False, sigma=1.0, term='M')
    yax.set_label('Lower_%d_Upper_%d_%s' % (dist_min, dist_max, "Mlp1_1875"))
    yax.set_weight(yaxial_weight)
    yax.add_to_model()
    outputobjects.append(yax)
    print(yax.get_output())

    yax = IMP.pmi1.restraints.npc.YAxialPositionRestraint(simo, (1679,1679,"Mlp2"), lower_bound=dist_min, upper_bound=dist_max, consider_radius=False, sigma=1.0, term='M')
    yax.set_label('Lower_%d_Upper_%d_%s' % (dist_min, dist_max, "Mlp2_1679"))
    yax.set_weight(yaxial_weight)
    yax.add_to_model()
    outputobjects.append(yax)
    print(yax.get_output())

    """
    # Mlp1-Mlp2 hetero-dimer
    dist_min = 3.0
    dist_max = 15.0
    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(1875,1875,"Mlp1"), (1679,1679,"Mlp2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.set_label("Mlp1-Mlp2")
    dr.set_weight(dr_weight)
    dr.add_to_model()
    outputobjects.append(dr)
    print(dr.get_output())
    """


#####################################################
# Restraints setup - Membrane Localization + ALPS Motif
# Campelo et al, PLOS CompBio, 2014 (PMC3983069)
#####################################################
tor_th      = 45.0
tor_th_ALPS = 12.0
tor_R       = 390.0 + 150.0
tor_r       = 150.0 - tor_th/2.0
tor_r_ALPS  = 150.0 - tor_th_ALPS/2.0
msl_sigma   = 1.0
msl_weight  = 1000.0

# Transmembrane domains
if (is_membrane):
    print("\nMembraneSurfaceLocationRestraint !!")
    msl = IMP.pmi1.restraints.npc.MembraneSurfaceLocationRestraint(simo, (111,194,'Pom152'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Pom152_111_194')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print(msl.get_output())

    # Pom152 201-1337 should be outside of the lipid bilayer
    msl = IMP.pmi1.restraints.npc.MembraneSurfaceLocationRestraint(simo, (201,1337,'Pom152'), tor_R=tor_R, tor_r=52.5, tor_th=105.0, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Pom152_201_1337')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print(msl.get_output())

    msl = IMP.pmi1.restraints.npc.MembraneSurfaceLocationRestraint(simo, (29,247,'Ndc1'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Ndc1')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print(msl.get_output())

    msl = IMP.pmi1.restraints.npc.MembraneSurfaceLocationRestraint(simo, (64,150,'Pom34'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Pom34')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print(msl.get_output())

    msl = IMP.pmi1.restraints.npc.MembraneSurfaceLocationRestraint(simo, (475,475,'Nup53'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup53')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print(msl.get_output())

    msl = IMP.pmi1.restraints.npc.MembraneSurfaceLocationRestraint(simo, (528,528,'Nup59'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup59')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print(msl.get_output())

# The Pom152 ring
if (is_membrane):
    dist_min = 3.0
    dist_max = 25.0

    # same residue cross-link of Pom152 62-62
    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(62,62,"Pom152"), (62,62,"Pom152@11"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.set_label("Pom152-Pom152@11")
    dr.set_weight(dr_weight)
    dr.add_to_model()
    outputobjects.append(dr)
    print(dr.get_output())

    """
    if (use_neighboring_spokes):
        # TODO: Pom152 orientation?  (clockwise or counter-clockwise?)
        #dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(351,351,"Pom152"), (1337,1337,"Pom152@12"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(351,351,"Pom152"), (1337,1337,"Pom152@13"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        #dr.set_label("Pom152-Pom152@12")
        dr.set_label("Pom152-Pom152@13")
        dr.set_weight(dr_weight)
        dr.add_to_model()
        outputobjects.append(dr)
        print(dr.get_output())

    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(351,351,"Pom152"), (351,351,"Pom152@11"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.set_label("Pom152-Pom152@11")
    dr.set_weight(dr_weight)
    dr.add_to_model()
    outputobjects.append(dr)
    print(dr.get_output())

    xyr = IMP.pmi1.restraints.npc.XYRadialPositionRestraint(simo, (859,859,"Pom152"), lower_bound=515, upper_bound=550, consider_radius=False, sigma=1.0, term='M')
    xyr.set_label('Lower_%d_Upper_%d_%s' % (515, 550, "Pom152_859"))
    xyr.set_weight(radial_weight)
    xyr.add_to_model()
    outputobjects.append(xyr)
    print(xyr.get_output())

    # TODO: Pom152 orientation?  (clockwise or counter-clockwise?)
    #yax = IMP.pmi1.restraints.npc.YAxialPositionRestraint(simo, (859,859,"Pom152"), lower_bound=180, upper_bound=205, consider_radius=False, sigma=1.0, term='M')
    yax = IMP.pmi1.restraints.npc.YAxialPositionRestraint(simo, (859,859,"Pom152"), lower_bound=-205, upper_bound=-180, consider_radius=False, sigma=1.0, term='M')
    #yax.set_label('Lower_%d_Upper_%d_%s' % (180, 205, "Pom152_859"))
    yax.set_label('Lower_%d_Upper_%d_%s' % (-205, -180, "Pom152_859"))
    yax.set_weight(yaxial_weight)
    yax.add_to_model()
    outputobjects.append(yax)
    print(yax.get_output())

    pom152_min = 5;     pom152_max = 55;
    zax = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(simo, (379,379,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_379"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print(zax.get_output())

    zax = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(simo, (520,520,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_520"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print(zax.get_output())

    zax = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(simo, (616,616,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_616"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print(zax.get_output())

    zax = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(simo, (722,722,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_722"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print(zax.get_output())

    zax = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(simo, (824,824,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_824"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print(zax.get_output())

    zax = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(simo, (931,931,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_931"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print(zax.get_output())

    zax = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(simo, (1036,1036,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_1036"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print(zax.get_output())

    zax = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(simo, (1150,1150,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_1150"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print(zax.get_output())

    zax = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(simo, (1244,1244,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_1244"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print(zax.get_output())
    """

# ALPS Motifs
if (is_nucleoplasm):
    msl = IMP.pmi1.restraints.npc.MembraneSurfaceLocationRestraint(simo, (1,32,'Nup1'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup1')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print(msl.get_output())

    msl = IMP.pmi1.restraints.npc.MembraneSurfaceLocationRestraint(simo, (27,47,'Nup60.1'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup60.1')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print(msl.get_output())

    msl = IMP.pmi1.restraints.npc.MembraneSurfaceLocationRestraint(simo, (27,47,'Nup60.2'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup60.2')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print(msl.get_output())

# ALPS Motifs
if (is_inner_ring):
    msl = IMP.pmi1.restraints.npc.MembraneSurfaceLocationRestraint(simo, (310,334,'Nup157'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup157')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print(msl.get_output())

    yax = IMP.pmi1.restraints.npc.YAxialPositionRestraint(simo, (310,334,'Nup157'), lower_bound=-60, upper_bound=-10, consider_radius=False, sigma=1.0, term='M')
    yax.set_label('Lower_%d_Upper_%d_%s' % (-60, -10, "Nup157_ALPS"))
    yax.set_weight(yaxial_weight)
    yax.add_to_model()
    outputobjects.append(yax)
    print(yax.get_output())

    zax = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(simo, (310,334,'Nup157'), lower_bound=-15, upper_bound=35, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (-15, 35, "Nup157_ALPS"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print(zax.get_output())


    msl = IMP.pmi1.restraints.npc.MembraneSurfaceLocationRestraint(simo, (320,344,'Nup170'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup170')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print(msl.get_output())

    yax = IMP.pmi1.restraints.npc.YAxialPositionRestraint(simo, (320,344,'Nup170'), lower_bound=-170, upper_bound=-120, consider_radius=False, sigma=1.0, term='M')
    yax.set_label('Lower_%d_Upper_%d_%s' % (-170, -120, "Nup170_ALPS"))
    yax.set_weight(yaxial_weight)
    yax.add_to_model()
    outputobjects.append(yax)
    print(yax.get_output())

    zax = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(simo, (320,344,'Nup170'), lower_bound=-60, upper_bound=-10, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (-60, -10, "Nup170_ALPS"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print(zax.get_output())

# ALPS Motifs
if (is_n84):
    msl = IMP.pmi1.restraints.npc.MembraneSurfaceLocationRestraint(simo, (252,270,'Nup133'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup133')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print(msl.get_output())

    msl = IMP.pmi1.restraints.npc.MembraneSurfaceLocationConditionalRestraint(simo, protein1=(135,152,'Nup120'), protein2=(197,216,'Nup120'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma, resolution = res_ev)
    msl.set_label('Nup120_135-152_197-216')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print(msl.get_output())


#####################################################
# Restraints setup - Membrane Exclusion
#####################################################
tor_th      = 150.0 - tor_th_ALPS
tor_R       = 390.0 + 150.0
tor_r       = tor_th/2.0
mex_sigma   = 0.2
mex_weight  = 1000.0

if (use_MembraneExclusion):
    nup_list = [entry[0] for entry in domains if not '@' in entry[0]]
    nup_list_unique = sorted(list(set(nup_list)))   # Make a unique list

    MEX_LIST = [
        [0, 0, "Nup84"],
        #[0, 0, "Nup85"],
        [0, 0, "Nup120"],
        [0, 0, "Nup133"],
        [0, 0, "Nup145c"],
        #[0, 0, "Seh1"],
        #[0, 0, "Sec13"],

        #[0, 0, "Nup82.1"],
        #[0, 0, "Nup82.2"],
        #[0, 0, "Nsp1.1"],
        #[0, 0, "Nsp1.2"],
        #[0, 0, "Nup159.1"],
        #[0, 0, "Nup159.2"],
        [1, 965, "Nup116.1"],
        [1, 965, "Nup116.2"],
        #[0, 0, "Gle1"],
        #[0, 0, "Nup42"],

        #[0, 0, "Nic96.1"],
        #[0, 0, "Nic96.2"],
        #[0, 0, "Nsp1.3"],
        #[0, 0, "Nsp1.4"],
        #[0, 0, "Nup57.1"],
        #[0, 0, "Nup57.2"],
        #[0, 0, "Nup49.1"],
        #[0, 0, "Nup49.2"],

        [1, 110, "Pom152"],
        [248, 655, "Ndc1"],
        [151, 299, "Pom34"],
        [0, 0, "Nup157"],
        [0, 0, "Nup170"],
        [0, 0, "Nup53"],
        [0, 0, "Nup59"],
        #[0, 0, "Nup188"],
        #[0, 0, "Nup192"],

        [0, 0, "Nup1"],
        [0, 0, "Nup60.1"],
        [0, 0, "Nup60.2"],
        [0, 0, "Nup145.1"],
        [0, 0, "Nup145.2"],
        [0, 0, "Nup100.1"],
        [0, 0, "Nup100.2"],
        [0, 0, "Mlp1"],
        [0, 0, "Mlp2"],
    ]
    print("\nMembraneExclusionRestraint !!")
    for z in MEX_LIST:
        if (z[2] not in nup_list_unique):
            continue
        if (z[0] > 0):
            mex = IMP.pmi1.restraints.npc.MembraneExclusionRestraint(simo, (z[0], z[1], z[2]), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=mex_sigma, resolution = res_ev)
        else:
            mex = IMP.pmi1.restraints.npc.MembraneExclusionRestraint(simo, z[2], tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=mex_sigma, resolution = res_ev)
        mex.set_label('%s_mex_%d_%d' % (z[2], z[0], z[1]))
        mex.set_weight(mex_weight)
        mex.add_to_model()
        outputobjects.append(mex)
        print(mex.get_output())
