# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the main file.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

from pathlib import Path
from datetime import datetime
import numpy as np
import os
import shutil
import math

#Own functions and classes
import Grain
import Owntools
import User
import Report

#-------------------------------------------------------------------------------
#Plan simulation
#-------------------------------------------------------------------------------

if Path('Input').exists():
    shutil.rmtree('Input')
os.mkdir('Input')
if Path('Output').exists():
    shutil.rmtree('Output')
os.mkdir('Output')
if Path('Data').exists():
    shutil.rmtree('Data')
os.mkdir('Data')
if Path('Debug').exists():
    shutil.rmtree('Debug')
os.mkdir('Debug')
os.mkdir('Debug/Configuration')
os.mkdir('Debug/Comparison_Init_Current')

#-------------------------------------------------------------------------------
#Create a simulation
#-------------------------------------------------------------------------------

#create a simulation report
simulation_report = Report.Report('Debug/Report',datetime.now())
simulation_report.tic_tempo(datetime.now())

#general parameters
dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
if dict_algorithm['SaveData']:
    if not Path('../'+dict_algorithm['foldername']).exists():
        os.mkdir('../'+dict_algorithm['foldername'])
    #tempo save of the user file
    shutil.copy('User.py','../'+dict_algorithm['foldername']+'/User_'+dict_algorithm['namefile']+'_tempo.txt')

#create the two grains
User.Add_2grains(dict_material,dict_sample)
#Compute initial sum_eta
Owntools.Compute_sum_eta(dict_sample)
#Compute the surface of the contact initially
User.Add_S0(dict_sample, dict_sollicitation)
#Compute the sphericity initially for the first grain
dict_sample['L_g'][0].geometric_study(dict_sample)
dict_sample['L_g'][0].Compute_sphericity(dict_algorithm)
#create the solute
User.Add_solute(dict_sample)

#plot
Owntools.Plot_config(dict_algorithm, dict_sample)

simulation_report.tac_tempo(datetime.now(),'Initialisation')

#-------------------------------------------------------------------------------
#trackers
#-------------------------------------------------------------------------------

dict_tracker = {
'L_displacement' : [0],
'L_sum_solute' : [0],
'L_sum_eta' : [dict_sample['sum_eta']],
'L_sum_total' : [dict_sample['sum_eta']],
'L_area_sphericity_g0' : [dict_sample['L_g'][0].area_sphericity],
'L_diameter_sphericity_g0' : [dict_sample['L_g'][0].diameter_sphericity],
'L_circle_ratio_sphericity_g0' : [dict_sample['L_g'][0].circle_ratio_sphericity],
'L_perimeter_sphericity_g0' : [dict_sample['L_g'][0].perimeter_sphericity],
'L_width_to_length_ratio_sphericity_g0' : [dict_sample['L_g'][0].width_to_length_ratio_sphericity]
}

#-------------------------------------------------------------------------------
#main
#-------------------------------------------------------------------------------

dict_algorithm['i_PFDEM'] = 0
while not User.Criteria_StopSimulation(dict_algorithm):

    #prepare iteration
    simulation_report.tic_tempo(datetime.now())
    dict_algorithm['i_PFDEM'] = dict_algorithm['i_PFDEM'] + 1
    simulation_report.write_and_print(f"\nITERATION {dict_algorithm['i_PFDEM']} / {dict_algorithm['n_t_PFDEM']}\n\n",f"\nITERATION {dict_algorithm['i_PFDEM']} / {dict_algorithm['n_t_PFDEM']}\n")
    os.mkdir('Output/Ite_'+str(dict_algorithm['i_PFDEM']))

    #move grain
    Grain.Compute_overlap_2_grains(dict_sample)
    Grain.Apply_overlap_target(dict_material,dict_sample,dict_sollicitation,dict_tracker)

    #compute the intersection surface
    Owntools.Compute_S_int(dict_sample)

    #plot
    Owntools.Plot_config(dict_algorithm, dict_sample)

    #write data
    Owntools.Write_eta_txt(dict_algorithm, dict_sample)
    Owntools.Write_solute_txt(dict_algorithm, dict_sample)
    Owntools.Write_kc_txt(dict_algorithm, dict_material, dict_sample)
    Owntools.Write_ep_txt(dict_algorithm, dict_sample)

    #create i
    Owntools.Create_i(dict_algorithm,dict_sample,dict_material)

    simulation_report.tac_tempo(datetime.now(),f"Iteration {dict_algorithm['i_PFDEM']}: preparation of the pf simulation")
    simulation_report.tic_tempo(datetime.now())

    #---------------------------------------------------------------------------
    #PF simulation
    #---------------------------------------------------------------------------

    #run
    os.system('mpiexec -n '+str(dict_algorithm['np_proc'])+' ~/projects/moose/modules/combined/combined-opt -i '+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'.i')

    simulation_report.tac_tempo(datetime.now(),f"Iteration {dict_algorithm['i_PFDEM']}: pf simulation")
    simulation_report.tic_tempo(datetime.now())

    #sorting files
    j_str = Owntools.Sort_Files(dict_algorithm)

    #---------------------------------------------------------------------------
    #PF to DEM
    #---------------------------------------------------------------------------

    #look for the new grains shape
    for grain in dict_sample['L_g']:
        grain.PFtoDEM_Multi('Output/Ite_'+str(dict_algorithm['i_PFDEM'])+'/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str,dict_algorithm,dict_sample)
        grain.geometric_study(dict_sample)
    #look for the new solute shape
    Owntools.solute_PFtoDEM_Multi('Output/Ite_'+str(dict_algorithm['i_PFDEM'])+'/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str,dict_algorithm,dict_sample)

    #plot
    Owntools.Plot_config(dict_algorithm, dict_sample)
    Owntools.Plot_init_current_shape(dict_algorithm, dict_sample)

    #---------------------------------------------------------------------------
    #postprocess
    #---------------------------------------------------------------------------

    #Compute the sphericity for the first grain
    dict_sample['L_g'][0].Compute_sphericity(dict_algorithm)

    #compute the mass of grain
    Owntools.Compute_sum_eta(dict_sample)

    #compute the mass of the solute
    Owntools.Compute_sum_c(dict_sample)

    #---------------------------------------------------------------------------
    #tracker
    #---------------------------------------------------------------------------

    dict_tracker['L_sum_solute'].append(dict_sample['sum_c'])
    dict_tracker['L_sum_eta'].append(dict_sample['sum_eta'])
    dict_tracker['L_sum_total'].append(dict_sample['sum_c']+dict_sample['sum_eta'])
    dict_tracker['L_area_sphericity_g0'].append(dict_sample['L_g'][0].area_sphericity)
    dict_tracker['L_diameter_sphericity_g0'].append(dict_sample['L_g'][0].diameter_sphericity)
    dict_tracker['L_circle_ratio_sphericity_g0'].append(dict_sample['L_g'][0].circle_ratio_sphericity)
    dict_tracker['L_perimeter_sphericity_g0'].append(dict_sample['L_g'][0].perimeter_sphericity)
    dict_tracker['L_width_to_length_ratio_sphericity_g0'].append(dict_sample['L_g'][0].width_to_length_ratio_sphericity)

    #Plot trackers
    Owntools.Plot_trackers(dict_tracker)

    #---------------------------------------------------------------------------
    #tempo save
    #---------------------------------------------------------------------------

    if dict_algorithm['SaveData']:
        Owntools.save_dicts_tempo(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker)
        shutil.copy('Debug/Report.txt','../'+dict_algorithm['foldername']+'/Report_'+dict_algorithm['namefile']+'_tempo.txt')

    simulation_report.tac_tempo(datetime.now(),f"Iteration {dict_algorithm['i_PFDEM']}: from pf to dem")

#-------------------------------------------------------------------------------
#close simulation
#-------------------------------------------------------------------------------

simulation_report.end(datetime.now())

#final save
if dict_algorithm['SaveData']:
    Owntools.save_dicts_final(dict_algorithm, dict_material, dict_sample, dict_sollicitation, dict_tracker)
    name_actual_folder = os.path.dirname(os.path.realpath(__file__))
    shutil.copytree(name_actual_folder, '../'+dict_algorithm['foldername']+'/'+dict_algorithm['namefile'])
    os.remove('../'+dict_algorithm['foldername']+'/User_'+dict_algorithm['namefile']+'_tempo.txt')
    os.remove('../'+dict_algorithm['foldername']+'/Report_'+dict_algorithm['namefile']+'_tempo.txt')
