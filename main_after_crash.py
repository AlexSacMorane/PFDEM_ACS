# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the file to restart a simulation after a crash.
There is a save at the end of each PFDEM iteration.
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
import pickle

#Own functions and classes
import Grain
import Owntools
import User
import Report

#-------------------------------------------------------------------------------
#Plan the simulation
#-------------------------------------------------------------------------------

#create a simulation report
simulation_report = Report.Report('Debug/Report_after_crash',datetime.now())

#-------------------------------------------------------------------------------
#load data
#-------------------------------------------------------------------------------

toload = open('../Data_2G_ACS/PF_AC_PS_1_save_tempo','rb')
dict_save = pickle.load(toload,encoding = 'bytes')
toload.close()
dict_algorithm = dict_save['algorithm']
dict_material = dict_save['material']
dict_sample = dict_save['sample']
dict_sollicitation = dict_save['sollicitation']
dict_tracker = dict_save['tracker']

#-------------------------------------------------------------------------------
#main
#-------------------------------------------------------------------------------

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

    #---------------------------------------------------------------------------
    #postprocess
    #---------------------------------------------------------------------------

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
