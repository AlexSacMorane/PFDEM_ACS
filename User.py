# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the file where the user can change the different parameters for the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

import math
import numpy as np
from pathlib import Path

#Own functions and classes
import Grain

#-------------------------------------------------------------------------------
#User
#-------------------------------------------------------------------------------

def All_parameters():
    '''
    This function is called in main.py to have all the parameters needed in the simulation.

        Input :
            Nothing
        Output :
            an algorithm dictionnary (a dictionnary)
            a material dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
            a sollicitation dictionnary (a dictionnary)
    '''
    #---------------------------------------------------------------------------
    #Sample parameters

    #spatial discretisation
    x_min = -210
    x_max = 210
    nx = 180
    x_L = np.linspace(x_min,x_max,nx)
    y_min = -130
    y_max = 130
    ny = 100
    y_L = np.linspace(y_min,y_max,ny)

    #approximatively the number of vertices for one grain during DEM simulation
    grain_discretisation = 80

    dict_sample = {
    'x_L' : x_L,
    'y_L' : y_L,
    'grain_discretisation' : grain_discretisation,
    'Emec_M' : np.array(np.zeros((len(y_L),len(x_L)))),
    'Eche_M' : np.array(np.zeros((len(y_L),len(x_L)))),
    'Ed_M' : np.array(np.zeros((len(y_L),len(x_L))))
    }

    #---------------------------------------------------------------------------
    #Algorithm parameters

    np_proc = 4 #number of processor used
    n_t_PFDEM = 200 #number of cycle PF-DEM

    #Time step for phase field
    n_t_PF = 6
    dt_PF_init = 0.3
    dt_PF_level1 = dt_PF_init/2
    dt_PF_level2 = dt_PF_level1/2
    dt_PF_level3 = dt_PF_level2/2

    #criteria to switch level
    Ed_level1 = 3
    Ed_level2 = 5
    Ed_level3 = 10

    #Visual parameters
    c_min = 0
    c_max = 0.1

    #Margin for sphericity study
    sphericity_margin = 0.05
    #Discretisation to find the inscribing (number of nodes in one direction)
    n_spatial_inscribing = 100

    #structural element for dilation
    struct_element = np.array(np.ones((10,6)), dtype = bool)

    #List of plot to do
    # Config, C_at_P, Diff_Solute, dt, Ed, Eta_c, Init_Current_Shape, Kc, Movie (need Config to work), Sint_MinEtai, Sphericity, sum_Ed
    L_flag_plot = ['Config', 'C_at_P', 'dt', 'Sphericity', 'Movie', 'sum_Ed']

    #Find a simulation ame
    template = 'PS_Long_Run' #template of the name of the simulation
    SaveData = False #Save data or not
    foldername = 'Data_2G_ACS' #name of the folder where data are saved
    if SaveData :
        i_run = 1
        folderpath = Path('../'+foldername+'/'+template+'_'+str(i_run))
        while folderpath.exists():
            i_run = i_run + 1
            folderpath = Path('../'+foldername+'/'+template+'_'+str(i_run))
        namefile = template+'_'+str(i_run)
    else :
        namefile = template

    dict_algorithm = {
    'struct_element' : struct_element,
    'c_min' : c_min,
    'c_max' : c_max,
    'sphericity_margin' : sphericity_margin,
    'n_spatial_inscribing' : n_spatial_inscribing,
    'np_proc' : np_proc,
    'SaveData' : SaveData,
    'namefile' : namefile,
    'dt_PF' : dt_PF_init,
    'dt_PF_init' : dt_PF_init,
    'dt_PF_level1' : dt_PF_level1,
    'dt_PF_level2' : dt_PF_level2,
    'dt_PF_level3' : dt_PF_level3,
    'Ed_level1' : Ed_level1,
    'Ed_level2' : Ed_level2,
    'Ed_level3' : Ed_level3,
    'n_t_PF' : n_t_PF,
    'foldername' : foldername,
    'n_t_PFDEM' : n_t_PFDEM,
    'L_flag_plot' : L_flag_plot
    }

    #---------------------------------------------------------------------------
    #External sollicitation parameters

    overlap_target = 10 #overlap verified before each phase field iteration
    chi = 0.2 #coefficient applied to the chemical energy

    dict_sollicitation = {
    'overlap_target' : overlap_target,
    'chi' : chi
    }

    #---------------------------------------------------------------------------
    #Material parameters

    #phase field
    width_int = math.sqrt((x_L[4]-x_L[0])**2+(y_L[4]-y_L[0])**2)
    Mobility = 1 #L1, L2 in .i
    #Diffusion of etai
    kappa_eta = 0.01
    #Energy barrier of etai
    Energy_barrier = 20*kappa_eta/(width_int)**2
    #Diffusion of c, the solute generated by the dissolution
    kappa_c = 50
    method_to_compute_kc = 'dilation' #dilation or interpolation
    #Definition of the evolution of the penalty term inside a grain (for interpolation)
    tau_kappa_c = 5 #kc = kc0 e(-(R_i-d)/R_i/tau_kappa_c)
    #Define the spring value
    Y = 70*10**9*10**(-6*2)
    nu = 0.3


    dict_material = {
    'w' : width_int,
    'M' : Mobility,
    'kappa_eta' : kappa_eta,
    'kappa_c' : kappa_c,
    'tau_kappa_c' : tau_kappa_c,
    'method_to_compute_kc' : method_to_compute_kc,
    'Energy_barrier' : Energy_barrier,
    'Y' : Y,
    'nu' : nu
    }

    #---------------------------------------------------------------------------

    return dict_algorithm, dict_material, dict_sample, dict_sollicitation

#-------------------------------------------------------------------------------

def Add_2grains(dict_material,dict_sample):
    '''
    Generate two grains in the sample.

        Input :
            a material dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
        Output :
            Nothing but the sample dictionnary gets a new grains configuration (a list)
    '''
    #grain 1
    radius = 100
    center = np.array([np.mean(dict_sample['x_L'])-radius,np.mean(dict_sample['y_L'])])
    grain_1 = Grain.Grain(1,radius,center,dict_material,dict_sample)

    #grain 2
    radius = 100
    center = np.array([np.mean(dict_sample['x_L'])+radius,np.mean(dict_sample['y_L'])])
    grain_2 = Grain.Grain(2,radius,center,dict_material,dict_sample)

    #add element in dict
    dict_sample['L_g'] = [grain_1, grain_2]

#-------------------------------------------------------------------------------

def Add_S0(dict_sample, dict_sollicitation):
    '''
    Compute the initial intersection surface between two disk particles.
    See https://calculis.net/q/aire-intersection-disques-35

        Input :
            a sample dictionnary (a dictionnary)
        Output :
            Nothing but the sample dictionnary gets a new entry (a float)
    '''
    #notation from the forum
    R = dict_sample['L_g'][0].r_mean
    R_prime = dict_sample['L_g'][1].r_mean
    D = R + R_prime - dict_sollicitation['overlap_target']
    d  = (R**2 + D**2 - R_prime**2)/2/D
    d_prime = D - d
    #Part from particle 1
    S1 = R**2*math.acos(d/R)-d*math.sqrt(R**2-d**2)
    #Part from particle 2
    S2 = R_prime**2*math.acos(d_prime/R_prime)-d_prime*math.sqrt(R_prime**2-d_prime**2)

    #add element in dictionnary
    dict_sample['S_int_0'] = S1 + S2

#-------------------------------------------------------------------------------

def Add_solute(dict_sample):
    '''
    Generate solute in the sample.

        Input :
            a sample dictionnary (a dictionnary)
        Output :
            Nothing but the sample dictionnary gets a new solute configuration (a n_y x n_x numpy array)
    '''
    #solute
    solute_M = np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L'])))

    #add element in dict
    dict_sample['solute_M'] = solute_M

#-------------------------------------------------------------------------------

def Add_alpha_emec(dict_sample, dict_sollicitation):
    '''
    Compute the coefficient alpha applied to the mechanical energy.

        Input :
            a sample dictionnary (a dictionnary)
            a sollicitation dictionnary (a dictionnary)
        Output :
            Nothing but the sollicitation dictionnary gets a new attribut (a float)
    '''
    alpha =  0.2*dict_sample['S_int_0']

    #update dict
    dict_sollicitation['alpha'] = alpha

#-------------------------------------------------------------------------------

def Criteria_StopSimulation(dict_algorithm):
    '''
    Define a stop criteria for the PFDEM simulation.

    The simulation stops if the number of iterations reaches a maximum value.

        Input :
            an algorithm dictionnary (a dictionnary)
        Output :
            The result depends on the fact if the stop criteria is reached or not (a bool)
    '''
    Criteria_Verified = False
    if dict_algorithm['i_PFDEM'] >= dict_algorithm['n_t_PFDEM']:
        Criteria_Verified = True
    return Criteria_Verified
