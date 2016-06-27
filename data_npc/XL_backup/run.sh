#!/bin/bash
: '####################################
# 1. Wild Type ScNup82 complex DSS XLs
#######################################
setup_environment.sh python ./crosslinks_dataset_standardizer.py -f XL_wtNup82_DSS.csv
mv cross-links_standardized.csv XL_wtNup82_DSS_standardized.csv

setup_environment.sh python ./crosslinks_dataset_2copies.py -f XL_wtNup82_DSS_standardized.csv
mv cross-links_standardized_2copies.csv XL_wtNup82_DSS_standardized_2copies.csv
cp -pr XL_wtNup82_DSS_standardized_2copies.csv ..
'

: '# Homo-dimer XLs
Nup159.1 1384 Nup159.2 1384 1.0 171
Nup159.1 1414 Nup159.2 1414 1.0 557
Nup159.1 1417 Nup159.2 1417 1.0 22

Nup159.1 1387 Nup159.2 1387 1.0 589
Nup159.1 1432 Nup159.2 1432 1.0 153
Nup82.1 517 Nup82.2 517 1.0 399
'


: '####################################
# 2. Wild Type SkNup82 complex DSS XLs
#######################################
setup_environment.sh python ./crosslinks_dataset_standardizer.py -f XL_skNup82_DSS.csv
mv cross-links_standardized.csv XL_skNup82_DSS_standardized.csv

# Replace Diploid82 with Nup82, subtracting 1 residue of Diploid82 for skNup82 residues 359-472
setup_environment.sh python ./crosslinks_dataset_equivalent.py -f XL_skNup82_DSS_standardized.csv
mv cross-links_standardized_equiv.csv XL_skNup82_DSS_standardized_equiv.csv

setup_environment.sh python ./crosslinks_dataset_2copies.py -f XL_skNup82_DSS_standardized_equiv.csv
mv cross-links_standardized_2copies.csv XL_skNup82_DSS_standardized_equiv_2copies.csv
cp -pr XL_skNup82_DSS_standardized_equiv_2copies.csv ..
'

: '# remove redundant XLs and make unique IDs (+1000) - uniqueness test
Nsp1    	795	Diploid82	632	1	36
Nsp1	    795	Nup82	    632	1	36

Diploid82	671	Diploid82	632	1	77
Diploid82	671	Nup82	632	1	77

Diploid82	649	Diploid82	692	1	87
Diploid82	649	Nup82	692	1	271

Diploid82	632	Diploid82	660	1	133
Nup82	632	Diploid82	660	1	133

Nup159	1397	Diploid82	632	1	153
Nup159	1397	Nup82	    632	1	153

Diploid82	497	Diploid82	541	1	168
Nup82	497	Nup82	    541	1	238

Diploid82	269	Diploid82	129	1	197
Diploid82	269	Nup82	129	1	197

Nup82	246	Nup159	1447	    1	231
Diploid82	246	Nup159	1447	1	275

Diploid82	365	Diploid82	129	1	265
Diploid82	365	Nup82	129	1	265

Diploid82	632	Diploid82	621	1	313
Nup82	632	Diploid82	621	1	313

Diploid82	616	Diploid82	632	1	321
Diploid82	616	Nup82	632	1	321

Diploid82	246	Diploid82	632	1	323
Diploid82	246	Nup82	632	1	323

Diploid82	362	Diploid82	129	1	335
Diploid82	362	Nup82	129	1	335
'


: '####################################
# 3. Wild Type ScNup82 complex EDC XLs
#######################################
setup_environment.sh python ./crosslinks_dataset_standardizer_EDC.py -f XL_wtNup82_EDC.csv
mv cross-links_standardized.csv XL_wtNup82_EDC_standardized.csv

setup_environment.sh python ./crosslinks_dataset_2copies.py -f XL_wtNup82_EDC_standardized.csv
mv cross-links_standardized_2copies.csv XL_wtNup82_EDC_standardized_2copies.csv
cp -pr XL_wtNup82_EDC_standardized_2copies.csv ..
'

setup_environment.sh python ./crosslinks_dataset_Ambiguity.py -f XL_optimized.csv
