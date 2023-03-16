# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is to run several times the main file.
"""

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

import main

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

overlap_L = [10]
chi_L = [0.2]
kappa_c_L = [50]

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

for overlap in overlap_L :
    for chi in chi_L :
        for kappa_c in kappa_c_L :
            main.open_main(overlap, chi, kappa_c)
