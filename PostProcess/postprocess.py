# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

The goal of this file is to plot results from the pressure solution simulations
"""
#------------------------------------------------------------------------------
# Librairies
#------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import math
import pickle
import os
from pathlib import Path

#own
import Grain
import Report

#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# User
#------------------------------------------------------------------------------

overlap_L = [5, 10]
chi_L = [0.05, 0.1, 0.2]
kappa_c_L = [5, 50]

fontsize_title = 44
fontsize_axis = 36
linewidth_plot = 6
fontsize_legend = 26
fontsize_ticks = 30

#------------------------------------------------------------------------------
# Load data
#------------------------------------------------------------------------------

L_L_vertical_strain = []
L_L_t = []
L_L_log_vertical_strain = []
L_L_log_time = []
L_L_fit = []
L_L_ed_che = []
L_L_ed_mec = []
L_L_ed = []
L_L_area_sphericity_g0 = []
L_L_diameter_sphericity_g0 = []
L_L_circle_ratio_sphericity_g0 = []
L_L_perimeter_sphericity_g0 = []
L_L_width_to_length_ratio_sphericity_g0 = []

i_figure = 0

for overlap in overlap_L :
    for chi in chi_L :
        for kappa_c in kappa_c_L :
            if Path('../../Data_2G_ACS_Parametric_Study/'+str(overlap)+'_'+str(int(100*chi))+'_'+str(kappa_c)+'_run_1').exists():

                toload = open('../../Data_2G_ACS_Parametric_Study/'+str(overlap)+'_'+str(int(100*chi))+'_'+str(kappa_c)+'_run_1_save','rb')
                dict_save = pickle.load(toload,encoding = 'bytes')
                toload.close()
                dict_tracker = dict_save['tracker']

                #------------------------------------------------------------------------------
                # Save data
                #------------------------------------------------------------------------------

                L_L_t.append(dict_tracker['L_t'])
                L_L_area_sphericity_g0.append(dict_tracker['L_area_sphericity_g0'])
                L_L_diameter_sphericity_g0.append(dict_tracker['L_diameter_sphericity_g0'])
                L_L_circle_ratio_sphericity_g0.append(dict_tracker['L_circle_ratio_sphericity_g0'])
                L_L_perimeter_sphericity_g0.append(dict_tracker['L_perimeter_sphericity_g0'])
                L_L_width_to_length_ratio_sphericity_g0.append(dict_tracker['L_width_to_length_ratio_sphericity_g0'])

                #------------------------------------------------------------------------------
                # Work and save data
                #------------------------------------------------------------------------------

                #vertical strain
                L_vertical_strain = []
                for i in range(len(dict_tracker['L_int_displacement'])):
                    L_vertical_strain.append(dict_tracker['L_int_displacement'][i]/400)
                L_L_vertical_strain.append(L_vertical_strain)

                #Andrade fit
                L_log_vertical_strain = []
                L_log_time = []
                log_k_mean = 0
                n_mean = 0
                for i in range(1,len(dict_tracker['L_t'])):
                    L_log_vertical_strain.append(math.log(dict_tracker['L_int_displacement'][i]/400,10))
                    L_log_time.append(math.log(dict_tracker['L_t'][i],10))
                    if i > 0.2*len(dict_tracker['L_t']):
                        log_k_mean = log_k_mean + L_log_vertical_strain[-1] - 1/3*L_log_time[-1]
                        n_mean = n_mean + 1
                L_L_log_vertical_strain.append(L_log_vertical_strain)
                L_L_log_time.append(L_log_time)

                #compute Andrade fit
                log_k_mean = log_k_mean / n_mean
                L_fit = []
                for i in range(len(L_log_time)):
                    L_fit.append(L_log_time[i]*1/3+log_k_mean)
                L_L_fit.append(L_fit)

                #energy distribution (chemical, mechanical, total)
                L_ed_che = []
                for i in range(len(dict_tracker['c_at_the_center'])):
                    L_ed_che.append(dict_tracker['c_at_the_center'][i]*dict_save['sollicitation']['chi'])
                L_L_ed_che.append(L_ed_che)

                L_ed_mec = []
                for i in range(len(dict_tracker['sum_min_etai_L'])):
                    L_ed_mec.append(dict_save['sollicitation']['alpha']/dict_tracker['sum_min_etai_L'][i])
                L_L_ed_mec.append(L_ed_mec)

                L_ed = []
                for i in range(len(dict_tracker['sum_min_etai_L'])):
                    L_ed.append(L_ed_mec[i]-L_ed_che[i])
                L_L_ed.append(L_ed)

#------------------------------------------------------------------------------
# Data not available
#------------------------------------------------------------------------------

            else :

                L_L_vertical_strain.append([])
                L_L_t.append([])
                L_L_log_vertical_strain.append([])
                L_L_log_time.append([])
                L_L_fit.append([])
                L_L_ed_che.append([])
                L_L_ed_mec.append([])
                L_L_ed.append([])
                L_L_area_sphericity_g0.append([])
                L_L_diameter_sphericity_g0.append([])
                L_L_circle_ratio_sphericity_g0.append([])
                L_L_perimeter_sphericity_g0.append([])
                L_L_width_to_length_ratio_sphericity_g0.append([])

#------------------------------------------------------------------------------
# Print data (influence of overlap)
#------------------------------------------------------------------------------

#Influence of overlap
if not Path('overlap_influence').exists():
    os.mkdir('overlap_influence')

#times - strain
plt.figure(1,figsize=(16,9))
plt.suptitle('Times (s) - Strain (-)')
for i_chi in range(len(chi_L)):
    for i_kappa_c in range(len(kappa_c_L)):
        i_subplot = i_chi*len(kappa_c_L) + i_kappa_c
        plt.subplot(len(chi_L),len(kappa_c_L),i_subplot+1)
        plt.title(r'$\chi$ = '+str(chi_L[i_chi])+r' $\kappa_c$ = '+str(kappa_c_L[i_kappa_c]))
        for i_overlap in range(len(overlap_L)):
            i_list = i_overlap*len(chi_L)*len(kappa_c_L) + i_chi*len(kappa_c_L) + i_kappa_c
            plt.plot(L_L_t[i_list], L_L_vertical_strain[i_list], label='overlap = '+str(overlap_L[i_overlap]))
        plt.legend()
plt.savefig('overlap_influence/times_strain.png')
plt.close()

#log(times) - log(strain) (Andrade fit)
plt.figure(1,figsize=(16,9))
plt.suptitle('log(Times) - log(Strain)')
for i_chi in range(len(chi_L)):
    for i_kappa_c in range(len(kappa_c_L)):
        i_subplot = i_chi*len(kappa_c_L) + i_kappa_c
        plt.subplot(len(chi_L),len(kappa_c_L),i_subplot+1)
        plt.title(r'$\chi$ = '+str(chi_L[i_chi])+r' $\kappa_c$ = '+str(kappa_c_L[i_kappa_c]))
        for i_overlap in range(len(overlap_L)):
            i_list = i_overlap*len(chi_L)*len(kappa_c_L) + i_chi*len(kappa_c_L) + i_kappa_c
            plt.plot(L_L_log_time[i_list], L_L_log_vertical_strain[i_list], label='overlap = '+str(overlap_L[i_overlap]))
        plt.legend()
plt.savefig('overlap_influence/log_times_log_strain.png')
plt.close()

#strain - diameter sphericity
plt.figure(1,figsize=(16,9))
plt.suptitle('Strain (-) - Diameter sphericity (-)')
for i_chi in range(len(chi_L)):
    for i_kappa_c in range(len(kappa_c_L)):
        i_subplot = i_chi*len(kappa_c_L) + i_kappa_c
        plt.subplot(len(chi_L),len(kappa_c_L),i_subplot+1)
        plt.title(r'$\chi$ = '+str(chi_L[i_chi])+r' $\kappa_c$ = '+str(kappa_c_L[i_kappa_c]))
        for i_overlap in range(len(overlap_L)):
            i_list = i_overlap*len(chi_L)*len(kappa_c_L) + i_chi*len(kappa_c_L) + i_kappa_c
            plt.plot(L_L_vertical_strain[i_list], L_L_diameter_sphericity_g0[i_list], label='overlap = '+str(overlap_L[i_overlap]))
        plt.legend()
plt.savefig('overlap_influence/strain_diameter_sphericity.png')
plt.close()

#strain - circle sphericity
plt.figure(1,figsize=(16,9))
plt.suptitle('Strain (-) - Circle sphericity (-)')
for i_chi in range(len(chi_L)):
    for i_kappa_c in range(len(kappa_c_L)):
        i_subplot = i_chi*len(kappa_c_L) + i_kappa_c
        plt.subplot(len(chi_L),len(kappa_c_L),i_subplot+1)
        plt.title(r'$\chi$ = '+str(chi_L[i_chi])+r' $\kappa_c$ = '+str(kappa_c_L[i_kappa_c]))
        for i_overlap in range(len(overlap_L)):
            i_list = i_overlap*len(chi_L)*len(kappa_c_L) + i_chi*len(kappa_c_L) + i_kappa_c
            plt.plot(L_L_vertical_strain[i_list], L_L_circle_ratio_sphericity_g0[i_list], label='overlap = '+str(overlap_L[i_overlap]))
        plt.legend()
plt.savefig('overlap_influence/strain_circle_sphericity.png')
plt.close()

#------------------------------------------------------------------------------
# Print data (influence of chi)
#------------------------------------------------------------------------------

#Influence of overlap
if not Path('chi_influence').exists():
    os.mkdir('chi_influence')

#times - strain
plt.figure(1,figsize=(16,9))
plt.suptitle('Times (s) - Strain (-)')
for i_overlap in range(len(overlap_L)):
    for i_kappa_c in range(len(kappa_c_L)):
        i_subplot = i_overlap*len(kappa_c_L) + i_kappa_c
        plt.subplot(len(overlap_L),len(kappa_c_L),i_subplot+1)
        plt.title('overlap = '+str(overlap_L[i_overlap])+r' $\kappa_c$ = '+str(kappa_c_L[i_kappa_c]))
        for i_chi in range(len(chi_L)):
            i_list = i_overlap*len(chi_L)*len(kappa_c_L) + i_chi*len(kappa_c_L) + i_kappa_c
            plt.plot(L_L_t[i_list], L_L_vertical_strain[i_list], label=r'$\chi$ = '+str(chi_L[i_chi]))
        plt.legend()
plt.savefig('chi_influence/times_strain.png')
plt.close()

#log(times) - log(strain) (Andrade fit)
plt.figure(1,figsize=(16,9))
plt.suptitle('log(Times) - log(Strain)')
for i_overlap in range(len(overlap_L)):
    for i_kappa_c in range(len(kappa_c_L)):
        i_subplot = i_overlap*len(kappa_c_L) + i_kappa_c
        plt.subplot(len(overlap_L),len(kappa_c_L),i_subplot+1)
        plt.title('overlap = '+str(overlap_L[i_overlap])+r' $\kappa_c$ = '+str(kappa_c_L[i_kappa_c]))
        for i_chi in range(len(chi_L)):
            i_list = i_overlap*len(chi_L)*len(kappa_c_L) + i_chi*len(kappa_c_L) + i_kappa_c
            plt.plot(L_L_log_time[i_list], L_L_log_vertical_strain[i_list], label=r'$\chi$ = '+str(chi_L[i_chi]))
        plt.legend()
plt.savefig('chi_influence/log_times_log_strain.png')
plt.close()

#strain - diameter sphericity
plt.figure(1,figsize=(16,9))
plt.suptitle('Strain (-) - Diameter sphericity (-)')
for i_overlap in range(len(overlap_L)):
    for i_kappa_c in range(len(kappa_c_L)):
        i_subplot = i_overlap*len(kappa_c_L) + i_kappa_c
        plt.subplot(len(overlap_L),len(kappa_c_L),i_subplot+1)
        plt.title('overlap = '+str(overlap_L[i_overlap])+r' $\kappa_c$ = '+str(kappa_c_L[i_kappa_c]))
        for i_chi in range(len(chi_L)):
            i_list = i_overlap*len(chi_L)*len(kappa_c_L) + i_chi*len(kappa_c_L) + i_kappa_c
            plt.plot(L_L_vertical_strain[i_list], L_L_diameter_sphericity_g0[i_list], label=r'$\chi$ = '+str(chi_L[i_chi]))
        plt.legend()
plt.savefig('chi_influence/strain_diameter_sphericity.png')
plt.close()

#strain - circle sphericity
plt.figure(1,figsize=(16,9))
plt.suptitle('Strain (-) - Circle sphericity (-)')
for i_overlap in range(len(overlap_L)):
    for i_kappa_c in range(len(kappa_c_L)):
        i_subplot = i_overlap*len(kappa_c_L) + i_kappa_c
        plt.subplot(len(overlap_L),len(kappa_c_L),i_subplot+1)
        plt.title('overlap = '+str(overlap_L[i_overlap])+r' $\kappa_c$ = '+str(kappa_c_L[i_kappa_c]))
        for i_chi in range(len(chi_L)):
            i_list = i_overlap*len(chi_L)*len(kappa_c_L) + i_chi*len(kappa_c_L) + i_kappa_c
            plt.plot(L_L_vertical_strain[i_list], L_L_circle_ratio_sphericity_g0[i_list], label=r'$\chi$ = '+str(chi_L[i_chi]))
        plt.legend()
plt.savefig('chi_influence/strain_circle_sphericity.png')
plt.close()

#------------------------------------------------------------------------------
# Print data (influence of kappa_c)
#------------------------------------------------------------------------------

#Influence of overlap
if not Path('kappa_c_influence').exists():
    os.mkdir('kappa_c_influence')

#times - strain
plt.figure(1,figsize=(16,9))
plt.suptitle('Times (s) - Strain (-)')
for i_overlap in range(len(overlap_L)):
    for i_chi in range(len(chi_L)):
        i_subplot = i_overlap*len(chi_L) + i_chi
        plt.subplot(len(overlap_L),len(chi_L),i_subplot+1)
        plt.title('overlap = '+str(overlap_L[i_overlap])+r' $\chi$ = '+str(chi_L[i_chi]))
        for i_kappa_c in range(len(kappa_c_L)):
            i_list = i_overlap*len(chi_L)*len(kappa_c_L) + i_chi*len(kappa_c_L) + i_kappa_c
            plt.plot(L_L_t[i_list], L_L_vertical_strain[i_list], label=r'$\kappa_c$ = '+str(kappa_c_L[i_kappa_c]))
        plt.legend()
plt.savefig('kappa_c_influence/times_strain.png')
plt.close()

#log(times) - log(strain) (Andrade fit)
plt.figure(1,figsize=(16,9))
plt.suptitle('log(Times) - log(Strain)')
for i_overlap in range(len(overlap_L)):
    for i_chi in range(len(chi_L)):
        i_subplot = i_overlap*len(chi_L) + i_chi
        plt.subplot(len(overlap_L),len(chi_L),i_subplot+1)
        plt.title('overlap = '+str(overlap_L[i_overlap])+r' $\chi$ = '+str(chi_L[i_chi]))
        for i_kappa_c in range(len(kappa_c_L)):
            i_list = i_overlap*len(chi_L)*len(kappa_c_L) + i_chi*len(kappa_c_L) + i_kappa_c
            plt.plot(L_L_log_time[i_list], L_L_log_vertical_strain[i_list], label=r'$\kappa_c$ = '+str(kappa_c_L[i_kappa_c]))
        plt.legend()
plt.savefig('kappa_c_influence/log_times_log_strain.png')
plt.close()

#strain - diameter sphericity
plt.figure(1,figsize=(16,9))
plt.suptitle('Strain (-) - Diameter sphericity (-)')
for i_overlap in range(len(overlap_L)):
    for i_chi in range(len(chi_L)):
        i_subplot = i_overlap*len(chi_L) + i_chi
        plt.subplot(len(overlap_L),len(chi_L),i_subplot+1)
        plt.title('overlap = '+str(overlap_L[i_overlap])+r' $\chi$ = '+str(chi_L[i_chi]))
        for i_kappa_c in range(len(kappa_c_L)):
            i_list = i_overlap*len(chi_L)*len(kappa_c_L) + i_chi*len(kappa_c_L) + i_kappa_c
            plt.plot(L_L_vertical_strain[i_list], L_L_diameter_sphericity_g0[i_list], label=r'$\kappa_c$ = '+str(kappa_c_L[i_kappa_c]))
        plt.legend()
plt.savefig('kappa_c_influence/strain_diameter_sphericity.png')
plt.close()

#strain - circle sphericity
plt.figure(1,figsize=(16,9))
plt.suptitle('Strain (-) - Circle sphericity (-)')
for i_overlap in range(len(overlap_L)):
    for i_chi in range(len(chi_L)):
        i_subplot = i_overlap*len(chi_L) + i_chi
        plt.subplot(len(overlap_L),len(chi_L),i_subplot+1)
        plt.title('overlap = '+str(overlap_L[i_overlap])+r' $\chi$ = '+str(chi_L[i_chi]))
        for i_kappa_c in range(len(kappa_c_L)):
            i_list = i_overlap*len(chi_L)*len(kappa_c_L) + i_chi*len(kappa_c_L) + i_kappa_c
            plt.plot(L_L_vertical_strain[i_list], L_L_circle_ratio_sphericity_g0[i_list], label=r'$\kappa_c$ = '+str(kappa_c_L[i_kappa_c]))
        plt.legend()
plt.savefig('chi_influence/strain_circle_sphericity.png')
plt.close()
