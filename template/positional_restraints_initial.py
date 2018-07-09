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
        "Gle1"     : [200, 320],
        #"Gle2.1"   : [150, 330],
        #"Gle2.2"   : [150, 330],
        "Mlp1"     : [150, 550],    # no immuno-EM identified for MLPs
        #"Mlp2"     : [150, 550],
        "Ndc1"     : [280, 410],    #"Ndc1"     : [280, 400],   (Table S2 and Table S7 are different)
        "Nic96.1"  : [255, 350],    #"Nic96.1"  : [255, 525],
        "Nic96.2"  : [255, 350],    #"Nic96.2"  : [255, 525],
        "Nsp1.1"   : [155, 425],
        "Nsp1.2"   : [155, 425],
        "Nsp1.3"   : [155, 425],
        "Nsp1.4"   : [155, 425],
        "Nup1"     : [200, 360],    # redundant with the FG nup restraint below
        "Nup100.1" : [230, 330],
        "Nup100.2" : [230, 330],
        "Nup116.1" : [200, 350],    #"Nup116.1" : [250, 350],
        "Nup116.2" : [200, 350],    #"Nup116.2" : [250, 350],
        "Nup120"   : [250, 450],    #"Nup120"   : [250, 370],   (Table S2 and Table S7 are different)
        "Nup133"   : [300, 540],    #"Nup133"   : [300, 420],
        "Nup145c"  : [270, 520],    #"Nup145c"  : [270, 470],
        "Nup145.1" : [125, 395],
        "Nup145.2" : [125, 395],
        "Nup157"   : [190, 350],
        "Nup159.1" : [200, 430],    #"Nup159.1" : [250, 430],
        "Nup159.2" : [200, 430],    #"Nup159.2" : [250, 430],
        "Nup170"   : [170, 330],
        "Nup188"   : [200, 320],    #(Table S2 and Table S7 are different)
        "Nup192"   : [200, 320],
        "Nup42"    : [220, 400],
        "Nup49.1"  : [200, 300],
        "Nup49.2"  : [200, 300],
        "Nup53"    : [280, 410],    #"Nup53"    : [280, 380],
        "Nup57.1"  : [ 80, 300],
        "Nup57.2"  : [ 80, 300],
        "Nup59"    : [250, 410],    #"Nup59"    : [250, 370],
        "Nup60.1"  : [240, 400],
        "Nup60.2"  : [240, 400],
        "Nup82.1"  : [175, 505],
        "Nup82.2"  : [175, 505],
        "Nup84"    : [290, 520],    #"Nup84"    : [290, 450],
        "Nup85"    : [300, 520],    #"Nup85"    : [300, 420],
        "Pom34"    : [280, 410],    #"Pom34"    : [280, 380],
        "Seh1"     : [250, 370],
        "Pom152"   : [470, 630]     #"Pom152"   : [370, 630]     #(Table S2 and Table S7 are different)
    }
    print("\nXYRadialPositionRestraint !!")
    for protein, r in RADIAL.items():
        if (protein not in nup_list_unique):
            continue
        xyr = IMP.npc.npc_restraints.XYRadialPositionRestraint(simo, protein, lower_bound=r[0], upper_bound=r[1], consider_radius=False, sigma=1.0)
        xyr.set_label('Lower_%d_Upper_%d_%s' % (r[0], r[1], protein))
        xyr.set_weight(radial_weight)
        xyr.add_to_model()
        outputobjects.append(xyr)
        print (xyr.get_output())

    ZAXIAL = {
        #"Dyn2.1"   : [ 130, 400],
        #"Dyn2.2"   : [ 130, 400],
        "Gle1"     : [ 120, 180],       #(Table S2 and Table S7 are different)
        #"Gle2.1"   : [  20, 120],
        #"Gle2.2"   : [  20, 120],
        "Mlp1"     : [-400,-150],       # no immuno-EM identified for MLPs
        #"Mlp2"     : [-400,-150],
        "Ndc1"     : [ -10,  70],       #"Ndc1"     : [   0,  90],  (Table S2 and Table S7 are different)
        "Nic96.1"  : [  25, 175],
        "Nic96.2"  : [  25, 175],
        "Nsp1.1"   : [ 130, 400],       #"Nsp1.1"   : [   0, 120],
        "Nsp1.2"   : [ 130, 400],       #"Nsp1.2"   : [   0, 120],
        "Nsp1.3"   : [   0, 120],
        "Nsp1.4"   : [   0, 120],
        "Nup1"     : [-220,-140],
        "Nup100.1" : [  40, 220],       #"Nup100.1" : [  40, 120],
        "Nup100.2" : [  40, 220],       #"Nup100.2" : [  40, 120],
        "Nup116.1" : [  70, 330],       #"Nup116.1" : [  70, 150],
        "Nup116.2" : [  70, 330],       #"Nup116.2" : [  70, 150],
        "Nup120"   : [  65, 400],       #"Nup120"   : [  70, 150],
        "Nup133"   : [ 185, 400],       #"Nup133"   : [ 100, 200],
        "Nup145c"  : [  65, 400],       #"Nup145c"  : [  70, 150],
        "Nup145.1" : [-320, -50],       #"Nup145.1" : [-170, -50],
        "Nup145.2" : [-320, -50],       #"Nup145.2" : [-170, -50],
        "Nup157"   : [  40, 165],       #"Nup157"   : [   0,  95],  (Table S2 and Table S7 are different)
        "Nup159.1" : [ 130, 400],       #"Nup159.1" :  [120, 240],
        "Nup159.2" : [ 130, 400],       #"Nup159.2" :  [120, 240],
        "Nup170"   : [   0,  75],       #"Nup170"   : [   0, 100],  (Table S2 and Table S7 are different)
        "Nup188"   : [  40, 100],
        "Nup192"   : [  20, 100],
        "Nup42"    : [  70, 150],
        "Nup49.1"  : [   0,  50],       #"Nup49.1"  : [  40, 100],
        "Nup49.2"  : [ -50,   0],       #"Nup49.2"  : [  40, 100],
        "Nup53"    : [  20, 100],
        "Nup57.1"  : [ -20,  80],       #"Nup57.1"  : [   0,  80],  (Table S2 and Table S7 are different)
        "Nup57.2"  : [ -80,  20],       #"Nup57.2"  : [   0,  80],
        "Nup59"    : [  40, 120],
        "Nup60.1"  : [-200,-100],
        "Nup60.2"  : [-200,-100],
        "Nup82.1"  : [ 130, 400],       #"Nup82.1"  : [ 145, 295],
        "Nup82.2"  : [ 130, 400],       #"Nup82.2"  : [ 145, 295],
        "Nup84"    : [ 185, 400],       #"Nup84"    : [ 150, 170],
        "Nup85"    : [  65, 400],       #"Nup85"    : [ 140, 200],
        "Pom34"    : [ -10,  50],       #"Pom34"    : [   0,  65],  (Table S2 and Table S7 are different)
        "Seh1"     : [  65, 400],       #"Seh1"     : [  50, 170],
        "Pom152"   : [   5,  60]        #"Pom152"   : [   0,  95]   (Table S2 and Table S7 are different)
    }
    print("\nZAxialPositionRestraint !!")
    for protein, z in ZAXIAL.items():
        if (protein not in nup_list_unique):
            continue
        zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, protein, lower_bound=z[0], upper_bound=z[1], consider_radius=False, sigma=1.0)
        zax.set_label('Lower_%d_Upper_%d_%s' % (z[0], z[1], protein))
        zax.set_weight(zaxial_weight)
        zax.add_to_model()
        outputobjects.append(zax)
        print (zax.get_output())


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
        "Mlp1"     : [-400,-200],
        "Mlp2"     : [-400,-200]
    }
    print("\nZAxialPositionRestraints for Cytoplasm / Nucleoplasm / Basket !!")
    for protein, z in ZAXIAL.items():
        if (protein not in nup_list_unique):
            continue

        # distal basket restraints (10)
        zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, protein, lower_bound=z[0], upper_bound=z[1], consider_radius=False, sigma=1.0)
        zax.set_label('Lower_%d_Upper_%d_%s' % (z[0], z[1], protein))
        zax.set_weight(zaxial_weight)
        zax.add_to_model()
        outputobjects.append(zax)
        print (zax.get_output())


#####################################################
# Restraints setup - FG anchor restraints
#####################################################
if (use_FG_anchor and is_nic96 and not is_FG):
    nup_list = [entry[0] for entry in domains if not '@' in entry[0]]
    nup_list_unique = sorted(list(set(nup_list)))   # Make a unique list

    RADIAL = {
        "Nsp1.3"   : [180, 230],
        "Nup49.1"  : [180, 230],
        "Nup57.1"  : [180, 230],
        "Nsp1.4"   : [180, 230],
        "Nup49.2"  : [180, 230],
        "Nup57.2"  : [180, 230],
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
            xyr = IMP.npc.npc_restraints.XYRadialPositionRestraint(simo, protein, lower_bound=r[0], upper_bound=r[1], consider_radius=False, sigma=1.0, term='C')
        else:
            xyr = IMP.npc.npc_restraints.XYRadialPositionRestraint(simo, protein, lower_bound=r[0], upper_bound=r[1], consider_radius=False, sigma=1.0, term='N')
        xyr.set_label('Lower_%d_Upper_%d_%s' % (r[0], r[1], protein))
        xyr.set_weight(radial_weight)
        xyr.add_to_model()
        outputobjects.append(xyr)
        print (xyr.get_output())


#####################################################
# Restraints setup
# Distance restraints
#####################################################
dist_min = 3.0
dr_weight = 10.0

# Nup145n - Nup145c
if (is_n84 and is_nucleoplasm):
    dist_max = 10.0
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(605,605,"Nup145.1"), (1,1,"Nup145c@11"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.set_label("Nup145n-Nup145c@11")
    dr.add_to_model()
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

# Nup60 homo-dimer cross-link
if (is_nucleoplasm):
    dist_max = 25.0
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(151,151,"Nup60.1"), (151,151,"Nup60.2"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.set_label("Nup60.1_151-Nup60.2_151")
    dr.add_to_model()
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

"""
# Nup120 - Nup133 to form the outer ring  (Seo et al, PNAS 2009) ; Not sure if it is real
if (is_n84 and use_neighboring_spokes):
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(11,11,"Nup133"), (641,641,"Nup120@2"), distancemin=3.0, distancemax=35.0, resolution=1.0)
    dr.set_label("Nup133-Nup120@2")
    dr.add_to_model()
    dr.set_weight(10.0)
    outputobjects.append(dr)
    print(dr.get_output())

if (is_inner_ring):
    dist_max = 15.0
    # connection between NTD and CTD
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(892,892,"Nup157"), (900,900,"Nup157"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.set_label("Nup157N-C")
    dr.add_to_model()
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    # connection between NTD and CTD
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(992,992,"Nup170"), (1000,1000,"Nup170"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.set_label("Nup170N-C")
    dr.add_to_model()
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    # Enforcing a binary interface Between Nup157 #976 and Nup170 #1475 (based on a cross-link ID #221)
    if (False):
        dist_max = 30.0
        dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(976,976,"Nup157"), (1475,1475,"Nup170"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr.set_label("Nup157C-Nup170C")
        dr.add_to_model()
        dr.set_weight(dr_weight)
        outputobjects.append(dr)
        print(dr.get_output())
"""


#####################################################
# Restraints setup - Membrane Localization + ALPS Motif
#####################################################
tor_th      = 45.0
tor_th_ALPS = 5.0
tor_R       = 390.0 + 150.0
tor_r       = 150.0 - tor_th/2.0
tor_r_ALPS  = 150.0 - tor_th_ALPS/2.0
msl_sigma   = 1.0
msl_weight  = 1.0

# Transmembrane domains
if (is_membrane):
    print("\nMembraneSurfaceLocationRestraint !!")
    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (111,194,'Pom152'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma)
    msl.set_label('Pom152_101_200')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (29,247,'Ndc1'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma)
    msl.set_label('Ndc1')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (64,150,'Pom34'), tor_R=tor_R, tor_r=tor_r, tor_th=tor_th, sigma=msl_sigma)
    msl.set_label('Pom34')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (475,475,'Nup53'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma)
    msl.set_label('Nup53')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (528,528,'Nup59'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma)
    msl.set_label('Nup59')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

# The Pom152 ring
if (is_membrane):
    dist_max = 25.0

    # same residue cross-link of Pom152 62-62
    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(62,62,"Pom152"), (62,62,"Pom152@11"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.set_label("Pom152-Pom152@11")
    dr.add_to_model()
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    if (use_neighboring_spokes):
        # TODO: Pom152 orientation?  (clockwise or counter-clockwise?)
        #dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(351,351,"Pom152"), (1337,1337,"Pom152@12"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(351,351,"Pom152"), (1337,1337,"Pom152@13"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
        #dr.set_label("Pom152-Pom152@12")
        dr.set_label("Pom152-Pom152@13")
        dr.add_to_model()
        dr.set_weight(dr_weight)
        outputobjects.append(dr)
        print(dr.get_output())

    dr = IMP.pmi.restraints.basic.DistanceRestraint(simo,(351,351,"Pom152"), (351,351,"Pom152@11"), distancemin=dist_min, distancemax=dist_max, resolution=1.0)
    dr.set_label("Pom152-Pom152@11")
    dr.add_to_model()
    dr.set_weight(dr_weight)
    outputobjects.append(dr)
    print(dr.get_output())

    """
    xyr = IMP.npc.npc_restraints.XYRadialPositionRestraint(simo, (859,859,"Pom152"), lower_bound=515, upper_bound=550, consider_radius=False, sigma=1.0, term='M')
    xyr.set_label('Lower_%d_Upper_%d_%s' % (515, 550, "Pom152_859"))
    xyr.set_weight(radial_weight)
    xyr.add_to_model()
    outputobjects.append(xyr)
    print (xyr.get_output())

    # TODO: Pom152 orientation?  (clockwise or counter-clockwise?)
    #yax = IMP.npc.npc_restraints.YAxialPositionRestraint(simo, (859,859,"Pom152"), lower_bound=180, upper_bound=205, consider_radius=False, sigma=1.0, term='M')
    yax = IMP.npc.npc_restraints.YAxialPositionRestraint(simo, (859,859,"Pom152"), lower_bound=-205, upper_bound=-180, consider_radius=False, sigma=1.0, term='M')
    #yax.set_label('Lower_%d_Upper_%d_%s' % (180, 205, "Pom152_859"))
    yax.set_label('Lower_%d_Upper_%d_%s' % (-205, -180, "Pom152_859"))
    yax.set_weight(yaxial_weight)
    yax.add_to_model()
    outputobjects.append(yax)
    print (yax.get_output())

    pom152_min = 5;     pom152_max = 55;
    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (379,379,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_379"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (520,520,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_520"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (616,616,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_616"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (722,722,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_722"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (824,824,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_824"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (931,931,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_931"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (1036,1036,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_1036"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (1150,1150,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_1150"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (1244,1244,"Pom152"), lower_bound=pom152_min, upper_bound=pom152_max, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (pom152_min, pom152_max, "Pom152_1244"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())
    """

# ALPS Motifs
if (is_nucleoplasm):
    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (1,32,'Nup1'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma)
    msl.set_label('Nup1')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (27,47,'Nup60.1'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma)
    msl.set_label('Nup60.1')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (27,47,'Nup60.2'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma)
    msl.set_label('Nup60.2')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

"""
# ALPS Motifs
if (is_n84):
    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (252,270,'Nup133'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma)
    msl.set_label('Nup133')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (135,152,'Nup120'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma)
    msl.set_label('Nup120_1')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (197,216,'Nup120'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma)
    msl.set_label('Nup120_2')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

# ALPS Motifs
if (is_inner_ring):
    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (310,334,'Nup157'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma)
    msl.set_label('Nup157')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    yax = IMP.npc.npc_restraints.YAxialPositionRestraint(simo, (310,334,'Nup157'), lower_bound=-60, upper_bound=-10, consider_radius=False, sigma=1.0, term='M')
    yax.set_label('Lower_%d_Upper_%d_%s' % (-60, -10, "Nup157_ALPS"))
    yax.set_weight(yaxial_weight)
    yax.add_to_model()
    outputobjects.append(yax)
    print (yax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (310,334,'Nup157'), lower_bound=-15, upper_bound=35, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (-15, 35, "Nup157_ALPS"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())


    msl = IMP.npc.npc_restraints.MembraneSurfaceLocationRestraint(simo, (320,344,'Nup170'), tor_R=tor_R, tor_r=tor_r_ALPS, tor_th=tor_th_ALPS, sigma=msl_sigma)
    msl.set_label('Nup170')
    msl.set_weight(msl_weight)
    msl.add_to_model()
    outputobjects.append(msl)
    print (msl.get_output())

    yax = IMP.npc.npc_restraints.YAxialPositionRestraint(simo, (320,344,'Nup170'), lower_bound=-170, upper_bound=-120, consider_radius=False, sigma=1.0, term='M')
    yax.set_label('Lower_%d_Upper_%d_%s' % (-170, -120, "Nup170_ALPS"))
    yax.set_weight(yaxial_weight)
    yax.add_to_model()
    outputobjects.append(yax)
    print (yax.get_output())

    zax = IMP.npc.npc_restraints.ZAxialPositionRestraint(simo, (320,344,'Nup170'), lower_bound=-60, upper_bound=-10, consider_radius=False, sigma=1.0, term='M')
    zax.set_label('Lower_%d_Upper_%d_%s' % (-60, -10, "Nup170_ALPS"))
    zax.set_weight(zaxial_weight)
    zax.add_to_model()
    outputobjects.append(zax)
    print (zax.get_output())
"""

