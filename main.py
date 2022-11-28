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
#Compute the surface of the contact initially
User.Add_S0(dict_sample, dict_sollicitation)
#create the solute
User.Add_solute(dict_sample)

#plot
Owntools.Plot_config(dict_sample)

simulation_report.tac_tempo(datetime.now(),'Initialisation')

#-------------------------------------------------------------------------------
#trackers
#-------------------------------------------------------------------------------

dict_tracker = {
'L_displacement' : [],
'L_sum_solute' : [0],
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
    Owntools.Plot_config(dict_sample)

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
    #postprocess
    #---------------------------------------------------------------------------

    #compute the mass of each grain
    for grain in dict_sample['L_g']:
        sum_eta = 0
        for l in range(len(dict_sample['y_L'])):
            for c in range(len(dict_sample['x_L'])):
                sum_eta = sum_eta + grain.etai_M[l][c]
        simulation_report.write(f"Grain {grain.id} -> sum of etai = {sum_eta}\n")

    #compute the mass of the solute
    sum_c = 0
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            sum_c = sum_c + dict_sample['solute_M'][l][c]
    simulation_report.write(f"Solute -> sum of etai = {sum_c}\n")

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
    Owntools.Plot_config(dict_sample)

    #---------------------------------------------------------------------------
    #postprocess
    #---------------------------------------------------------------------------

    #compute the mass of each grain
    for grain in dict_sample['L_g']:
        sum_eta = 0
        for l in range(len(dict_sample['y_L'])):
            for c in range(len(dict_sample['x_L'])):
                sum_eta = sum_eta + grain.etai_M[l][c]
        simulation_report.write(f"Grain {grain.id} -> sum of etai = {sum_eta}\n")

    #compute the mass of the solute
    sum_c = 0
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            sum_c = sum_c + dict_sample['solute_M'][l][c]
    simulation_report.write(f"Solute -> sum of etai = {sum_c}\n")

    #---------------------------------------------------------------------------
    #tracker
    #---------------------------------------------------------------------------

    dict_tracker['L_sum_solute'].append(sum_c)

    #Plot trackers
    Owntools.Plot_trackers(dict_tracker)

    #---------------------------------------------------------------------------
    #tempo save
    #---------------------------------------------------------------------------

    if dict_algorithm['SaveData']:
        Owntools.save_dicts_tempo(dict_algorithm, dict_material, dict_sample, dict_sollicitation)
        shutil.copy('Debug/Report.txt','../'+dict_algorithm['foldername']+'/Report_'+dict_algorithm['namefile']+'_tempo.txt')

    simulation_report.tac_tempo(datetime.now(),f"Iteration {dict_algorithm['i_PFDEM']}: from pf to dem")

#-------------------------------------------------------------------------------
#close simulation
#-------------------------------------------------------------------------------

simulation_report.end(datetime.now())

#final save
if dict_algorithm['SaveData']:
    Owntools.save_dicts_final(dict_algorithm, dict_material, dict_sample, dict_sollicitation)
    name_actual_folder = os.path.dirname(os.path.realpath(__file__))
    shutil.copytree(name_actual_folder, '../'+dict_algorithm['foldername']+'/'+dict_algorithm['namefile'])
    os.remove('../'+dict_algorithm['foldername']+'/User_'+dict_algorithm['namefile']+'_tempo.txt')
    os.remove('../'+dict_algorithm['foldername']+'/Report_'+dict_algorithm['namefile']+'_tempo.txt')
