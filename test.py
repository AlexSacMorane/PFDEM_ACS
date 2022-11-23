# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This is the test file to verify if the code is running well.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

from pathlib import Path
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import unittest
import os
import shutil
import math

#own functions and classes
import User
import Owntools
import Grain
import Report

#-------------------------------------------------------------------------------
#Test
#-------------------------------------------------------------------------------

class TestExample(unittest.TestCase):

    def test_example(self):
        self.assertTrue(True)

#-------------------------------------------------------------------------------

class TestReport(unittest.TestCase):
#test functions from Report.py

    def test_Report(self):
        #try to create a report.txt file
        simulation_report = Report.Report('Report',datetime.now())
        #check if the .txt has been created
        self.assertTrue(Path('Report.txt').is_file(),"The file Report.txt has not been created by the function Report.Report()!")
        os.remove('Report.txt')

#-------------------------------------------------------------------------------

class TestUser(unittest.TestCase):
#test functions from User.py

    def test_All_parameters(self):
        #try to acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #check if all data are dictionnaries
        self.assertTrue((type(dict_algorithm) is dict) and (type(dict_material) is dict) and (type(dict_sample) is dict) and (type(dict_sollicitation) is dict),'Outputs from User.All_parameters() are not dictionaries!')

    #---------------------------------------------------------------------------

    def test_Add_2grains(self):
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #try to create 2 grains
        User.Add_2grains(dict_sample,dict_material)
        #check if there are 2 grains
        self.assertTrue(len(dict_sample['L_g'])==2,'The function User.Add_2grains does not create 2 grains!')

    #---------------------------------------------------------------------------

    def test_Add_solute(self):
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #try to create 2 grains
        User.Add_solute(dict_sample)
        #check if there are 2 grains
        self.assertTrue('solute_M' in dict_sample.keys(),'The function User.Add_solute does not create a solute!')

#-------------------------------------------------------------------------------

class TestOwntools(unittest.TestCase):
#test functions from Owntools.py

    def test_is_PF_ACS_base_here(self):
        #Check if the file PF_ACS_base.i is in the directory
        self.assertTrue(Path('PF_ACS_base.i').is_file(),"The file PF_CH_AC_base.i should exists!")

    #---------------------------------------------------------------------------

    def test_Create_i(self):
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        dict_algorithm['i_PFDEM'] = 0
        #try to create a .i file
        Owntools.Create_i(dict_algorithm,dict_sample,dict_material)
        #Check if the file PF_CH_AC_base.i is in the directory
        self.assertTrue(Path(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'.i').is_file(),"The file namefile.i has not been created!")
        os.remove(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'.i')

    #---------------------------------------------------------------------------

    def test_index_to_str(self):
        #check if the function works well in different configurations
        self.assertTrue(Owntools.index_to_str(7)=='007','The conversion index_to_str() seems to do not for 00x...')
        self.assertTrue(Owntools.index_to_str(26)=='026','The conversion index_to_str() seems to do not for 0xx...')
        self.assertTrue(Owntools.index_to_str(666)=='666','The conversion index_to_str() seems to do not for xxx...')

    #---------------------------------------------------------------------------

    def test_Plot_config(self):
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #create the two grains
        User.Add_2grains(dict_sample,dict_material)
        #create a folder
        if Path('Debug').exists():
            shutil.rmtree('Debug')
        os.mkdir('Debug')
        #try to plot the configuration
        Owntools.Plot_config(dict_sample)
        #check if the .png has been created
        self.assertTrue(Path('Debug/Configuration_0.png').is_file(),"The image Debug/Configuration_0.png has not been created by the function Owntools.Plot_config()!")
        shutil.rmtree('Debug')

    #---------------------------------------------------------------------------

    def test_Cosine_Profile(self):
        #check if the function works well in different configurations
        self.assertTrue(Owntools.Cosine_Profile(1,0,0.5)==1,'The Owntools.Cosine_Profile() seems to do not for a point inside the grain...')
        self.assertTrue(Owntools.Cosine_Profile(1,1,0.5)==0.5*(1 + math.cos(math.pi/2)),'The Owntools.Cosine_Profile() seems to do not for a point inside the interface...')
        self.assertTrue(Owntools.Cosine_Profile(1,2,0.5)==0,'The Owntools.Cosine_Profile() seems to do not for a point outside the grain...')

    #---------------------------------------------------------------------------

    def test_Write_eta_txt(self):
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        dict_algorithm['i_PFDEM'] = 0
        #create the two grains
        User.Add_2grains(dict_sample,dict_material)
        #create a folder
        if Path('Data').exists():
            shutil.rmtree('Data')
        os.mkdir('Data')
        #try to create .txt files
        Owntools.Write_eta_txt(dict_algorithm, dict_sample)
        #Check if the files are in the directory
        self.assertTrue(Path('Data/eta1_0.txt').is_file(),"The file Data/eta1_0.txt has not been created!")
        self.assertTrue(Path('Data/eta2_0.txt').is_file(),"The file Data/eta2_0.txt has not been created!")
        shutil.rmtree('Data')

    #---------------------------------------------------------------------------

    def test_Write_solute_txt(self):
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        dict_algorithm['i_PFDEM'] = 0
        #create the two grains
        User.Add_solute(dict_sample)
        #create a folder
        if Path('Data').exists():
            shutil.rmtree('Data')
        os.mkdir('Data')
        #try to create .txt files
        Owntools.Write_solute_txt(dict_algorithm, dict_sample)
        #Check if the files are in the directory
        self.assertTrue(Path('Data/c_0.txt').is_file(),"The file Data/c_0.txt has not been created!")
        shutil.rmtree('Data')

    #---------------------------------------------------------------------------

    def test_Write_ep_txt(self):
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        dict_algorithm['i_PFDEM'] = 0
        #create the two grains
        User.Add_2grains(dict_sample,dict_material)
        #create a folder
        if Path('Data').exists():
            shutil.rmtree('Data')
        os.mkdir('Data')
        #try to create .txt files
        Owntools.Write_ep_txt(dict_algorithm, dict_sample)
        #Check if the files are in the directory
        self.assertTrue(Path('Data/ep_0.txt').is_file(),"The file Data/ep_0.txt has not been created!")
        shutil.rmtree('Data')

    #---------------------------------------------------------------------------

    def test_Write_kc_txt(self):
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        dict_algorithm['i_PFDEM'] = 0
        #create the two grains
        User.Add_2grains(dict_sample,dict_material)
        #create a folder
        if Path('Data').exists():
            shutil.rmtree('Data')
        os.mkdir('Data')
        #try to create .txt files
        Owntools.Write_kc_txt(dict_algorithm, dict_material, dict_sample)
        #Check if the files are in the directory
        self.assertTrue(Path('Data/kc_0.txt').is_file(),"The file Data/kc_0.txt has not been created!")
        shutil.rmtree('Data')

#-------------------------------------------------------------------------------

class TestGrain(unittest.TestCase):
#test functions from Grain.py

    def test_geometric_study(self):
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #Create one grain
        grain = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L']),np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        #try to do the geometric study of the grain
        grain.geometric_study(dict_sample)
        #define margine of the surface and the center (the study is done with a Monte Carlo method, some noise can be introduced)
        margin_center = 0.03*10 #10 is the radius see line upper
        margin_surface = 0.03*math.pi*10**2 #10 is the radius see line upper
        #check if the center computed is near the analytical one
        self.assertTrue(np.linalg.norm(np.array([np.mean(dict_sample['x_L']),np.mean(dict_sample['y_L'])])-grain.center)<margin_center,'The estimation of the center position seems false (because of the Monte Carlo Method try to rerun or increase the margin)...')
        #check if the surface computed is near the analytical one
        self.assertTrue(abs(math.pi*10**2-grain.surface)<margin_surface,'The estimation of the surface seems false (because of the Monte Carlo Method try to rerun or increase the margin)...')

    #---------------------------------------------------------------------------

    def test_P_is_inside(self):
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #Create one grain
        grain = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L']),np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        #check if the function works well in different configurations
        self.assertFalse(grain.P_is_inside(np.array([np.mean(dict_sample['x_L'])+11,np.mean(dict_sample['y_L'])])),'An outside point is detected as inside by Grain.P_is_inside()...')
        self.assertTrue(grain.P_is_inside(np.array([np.mean(dict_sample['x_L']),np.mean(dict_sample['y_L'])+1])),'An inside point is detected as outside by Grain.P_is_inside()...')

    #---------------------------------------------------------------------------

    def test_move_grain_rebuild(self):
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #Create one grain
        grain = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L']),np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        #try to move the grain
        grain.move_grain_rebuild(np.array([5,0]),dict_material,dict_sample)
        #check if the center computed is near the analytical one
        self.assertTrue(np.linalg.norm(np.array([np.mean(dict_sample['x_L'])+5,np.mean(dict_sample['y_L'])])-grain.center)==0,'The grain has not been well moved!')
        #Study the geometric of the grain
        grain.geometric_study(dict_sample)
        #define margine of the new center (the study is done with a Monte Carlo method, some noise can be introduced)
        margin_center = 0.03*10 #10 is the radius see line upper
        #check if the center computed is near the analytical one
        self.assertTrue(np.linalg.norm(np.array([np.mean(dict_sample['x_L'])+5,np.mean(dict_sample['y_L'])])-grain.center)<margin_center,'The displacement of the grain seems false after the etai_M rebuild (because of the Monte Carlo Method try to rerun or increase the margin)...')

    #---------------------------------------------------------------------------

    def test_move_grain_interpolation(self):
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #Create one grain
        grain = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L']),np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        #try to move the grain
        grain.move_grain_interpolation(np.array([5,0]),dict_material,dict_sample)
        #check if the center computed is near the analytical one
        self.assertTrue(np.linalg.norm(np.array([np.mean(dict_sample['x_L'])+5,np.mean(dict_sample['y_L'])])-grain.center)==0,'The grain has not been well moved!')
        #Study the geometric of the grain
        grain.geometric_study(dict_sample)
        #define margine of the new center (the study is done with a Monte Carlo method, some noise can be introduced)
        margin_center = 0.03*10 #10 is the radius see line upper
        #check if the center computed is near the analytical one
        self.assertTrue(np.linalg.norm(np.array([np.mean(dict_sample['x_L'])+5,np.mean(dict_sample['y_L'])])-grain.center)<margin_center,'The displacement of the grain seems false after the etai_M rebuild (because of the Monte Carlo Method try to rerun or increase the margin)...')

    #---------------------------------------------------------------------------

    def test_Compute_overlap_2_grains(self):
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        dict_sample['grain_discretisation'] = 20
        #Create two grain
        g1 = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L'])-5,np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        g2 = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L'])+5,np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        dict_sample['L_g'] = [g1,g2]
        #Check if there is an overlap between those grains
        Grain.Compute_overlap_2_grains(dict_sample)
        self.assertTrue(dict_sample['overlap']==10,'The overlap between two grains in contact is not well computed!')
        #Create two grain
        g1 = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L'])-11,np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        g2 = Grain.Grain(0,10,np.array([np.mean(dict_sample['x_L'])+11,np.mean(dict_sample['y_L'])]),dict_material,dict_sample)
        dict_sample['L_g'] = [g1,g2]
        #Check if there is not an overlap between those grains
        Grain.Compute_overlap_2_grains(dict_sample)
        self.assertTrue(dict_sample['overlap']==-2,'The overlap between two grains not in contact is not well computed!')

    #---------------------------------------------------------------------------

    def test_Apply_overlap_target(self):
        #Acquire data
        dict_algorithm, dict_material, dict_sample, dict_sollicitation = User.All_parameters()
        #Create two grain
        User.Add_2grains(dict_sample,dict_material)
        #Compute the initial overlap
        Grain.Compute_overlap_2_grains(dict_sample)
        #try to apply a target overlap
        Grain.Apply_overlap_target(dict_material,dict_sample,dict_sollicitation)
        #Compute the current overlap
        Grain.Compute_overlap_2_grains(dict_sample)
        #Check if the overlap target is well applied
        self.assertTrue(dict_sample['overlap']==dict_sollicitation['overlap_target'],'The target overlap between two grains is not well applied!')

#-------------------------------------------------------------------------------
#main
#-------------------------------------------------------------------------------

unittest.main(verbosity = 2)
