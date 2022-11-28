# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions used in the simulation.
    ''      Compute the current overlap between two grains.
        It is assumed the sample is composed  by only 2 grains.

        Input :
            a sample dictionnary (a dictionnary)
    Output :
        Nothing but the sample dictionnary gets an updated value (a float)
  ''
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import os
import math
import pickle

#-------------------------------------------------------------------------------

def Create_i(dict_algorithm,dict_sample,dict_material):
  '''
  Create the .i file to run MOOSE simulation.

  The file is generated from a template nammed PF_ACS_base.i

    Input :
        a algorithm dictionnary (a dictionnary)
        a sample dictionnary (a dictionnary)
        a material dictionnary (a dictionnary)
    Output :
        Nothing but a txt .i file is generated (a file)
  '''
  file_to_write = open(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'.i','w')
  file_to_read = open('PF_ACS_base.i','r')
  lines = file_to_read.readlines()
  file_to_read.close()

  j = 0
  for line in lines :
    j = j + 1
    if j == 4:
      line = line[:-1] + ' ' + str(len(dict_sample['x_L'])-1)+'\n'
    elif j == 5:
      line = line[:-1] + ' ' + str(len(dict_sample['y_L'])-1)+'\n'
    elif j == 7:
      line = line[:-1] + ' ' + str(min(dict_sample['x_L']))+'\n'
    elif j == 8:
      line = line[:-1] + ' ' + str(max(dict_sample['x_L']))+'\n'
    elif j == 9:
      line = line[:-1] + ' ' + str(min(dict_sample['y_L']))+'\n'
    elif j == 10:
      line = line[:-1] + ' ' + str(max(dict_sample['y_L']))+'\n'
    elif j == 116:
      line = line[:-1] + "'"+str(dict_material['M'])+' '+str(dict_material['kappa_eta'])+"'\n"
    elif j == 136:
      line = line[:-1] + ' ' + str(dict_material['Energy_barrier'])+"'\n"
    elif j == 175 or j == 179 or j == 183 or j == 187 or j == 191:
      line = line[:-1] + str(dict_algorithm['i_PFDEM']) + '.txt\n'
    elif j == 221:
      line = line[:-1] + ' ' + str(dict_algorithm['dt_PF']*dict_algorithm['n_t_PF']) +'\n'
    elif j == 225:
      line = line[:-1] + ' ' + str(dict_algorithm['dt_PF']) +'\n'
    file_to_write.write(line)

  file_to_write.close()

#-------------------------------------------------------------------------------

def index_to_str(j):
  '''
  An integer is converted to a float with 3 components

    Input :
        a number (a float)
    Output :
        a string with 3 components (a string)
  '''
  if j < 10:
      j_str = '00'+str(j)
  elif 10 <= j and j < 100:
      j_str = '0'+str(j)
  else :
      j_str = str(j)
  return j_str

#-------------------------------------------------------------------------------

def Sort_Files(dict_algorithm):
     '''
     Sort files generated by MOOSE to different directories

        Input :
            an algorithm dictionnary (a dict)
        Output :
            Nothing but files are sorted
     '''

     os.rename(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_out.e','Output/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_out.e')
     os.rename(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'.i','Input/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'.i')
     j = 0
     j_str = index_to_str(j)
     filepath = Path(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'.pvtu')
     while filepath.exists():
         for i_proc in range(dict_algorithm['np_proc']):
            os.rename(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'_'+str(i_proc)+'.vtu','Output/Ite_'+str(dict_algorithm['i_PFDEM'])+'/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'_'+str(i_proc)+'.vtu')
         os.rename(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'.pvtu','Output/Ite_'+str(dict_algorithm['i_PFDEM'])+'/'+dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'.pvtu')
         j = j + 1
         j_str = index_to_str(j)
         filepath = Path(dict_algorithm['namefile']+'_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'.pvtu')

     return index_to_str(j-1)

#-------------------------------------------------------------------------------

def Compute_S_int(dict_sample):
    '''
    Searching Surface,
    #Monte Carlo Method
    #A box is defined, we take a random point and we look if it is inside or outside the grain
    #Properties are the statistic times the box properties


        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing but the dictionnary gets an updated value for the intersection surface (a float)
    '''
    #Defining the limit
    box_min_x = min(dict_sample['L_g'][1].l_border_x)
    box_max_x = max(dict_sample['L_g'][0].l_border_x)
    box_min_y = min(dict_sample['L_g'][0].l_border_y)
    box_max_y = max(dict_sample['L_g'][0].l_border_y)

    #Compute the intersection surface
    N_MonteCarlo = 5000 #The larger it is, the more accurate it is
    sigma = 1
    M_Mass = 0

    for i in range(N_MonteCarlo):
        P = np.array([random.uniform(box_min_x,box_max_x),random.uniform(box_min_y,box_max_y)])
        if dict_sample['L_g'][0].P_is_inside(P) and dict_sample['L_g'][1].P_is_inside(P):
            M_Mass = M_Mass + sigma

    Mass = (box_max_x-box_min_x)*(box_max_y-box_min_y)/N_MonteCarlo*M_Mass
    Surface = Mass/sigma

    #Update element in dictionnary
    dict_sample['S_int'] = Surface

#-------------------------------------------------------------------------------

def Plot_config(dict_sample):
    '''
    Plot the sample configuration.

        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing but a .png file is generated (a file)
    '''
    #look for the name of the new plot
    template_name = 'Debug/Configuration_'
    j = 0
    plotpath = Path(template_name+str(j)+'.png')
    while plotpath.exists():
        j = j + 1
        plotpath = Path(template_name+str(j)+'.png')
    name = template_name+str(j)+'.png'

    #plot
    plt.figure(1,figsize=(16,9))
    for i in range(len(dict_sample['L_g'])):
        plt.plot(dict_sample['L_g'][i].l_border_x,dict_sample['L_g'][i].l_border_y)
    plt.axis('equal')
    plt.xlim(min(dict_sample['x_L']),max(dict_sample['x_L']))
    plt.savefig(name)
    plt.close(1)

#-------------------------------------------------------------------------------

def Cosine_Profile(R,r,w):
    '''
    Compute the phase field variable at some point.

    A cosine profile is assumed (see https://mooseframework.inl.gov/source/ics/SmoothCircleIC.html).

        Input :
            the radius R of the grain in the direction (a float)
            the distance r between the current point and the center (a float)
            the width w of the interface (a float)
        Output :
            the value of the phase field variable (a float)
    '''
    #inside the grain
    if r<R-w/2:
        return 1
    #outside the grain
    elif r>R+w/2:
        return 0
    #inside the interface
    else :
        return 0.5*(1 + np.cos(math.pi*(r-R+w/2)/w))

#-------------------------------------------------------------------------------

def Write_eta_txt(dict_algorithm, dict_sample):
    '''
    Write a .txt file needed for MOOSE simulation.

    The variables etai are transmitted to the MOOSE simulation.
    It is assumed the sample is composed by only two grains.

        Input :
            an algorithm dictionnary (a dict)
            an sample dictionnary (a dict)
        Output :
            Nothing but a .txt file is generated (a file)
    '''
    grain1 = dict_sample['L_g'][0]
    grain2 = dict_sample['L_g'][1]

    file_to_write_1 = open('Data/eta1_'+str(dict_algorithm['i_PFDEM'])+'.txt','w')
    file_to_write_2 = open('Data/eta2_'+str(dict_algorithm['i_PFDEM'])+'.txt','w')
    file_to_write_1.write('AXIS X\n')
    file_to_write_2.write('AXIS X\n')
    line = ''
    for x in dict_sample['x_L']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write_1.write(line)
    file_to_write_2.write(line)

    file_to_write_1.write('AXIS Y\n')
    file_to_write_2.write('AXIS Y\n')
    line = ''
    for y in dict_sample['y_L']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write_1.write(line)
    file_to_write_2.write(line)

    file_to_write_1.write('DATA\n')
    file_to_write_2.write('DATA\n')
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):

            #grain 1
            file_to_write_1.write(str(grain1.etai_M[-1-l][c])+'\n')
            #grain 2
            file_to_write_2.write(str(grain2.etai_M[-1-l][c])+'\n')

    file_to_write_1.close()
    file_to_write_2.close()

#-------------------------------------------------------------------------------

def Write_solute_txt(dict_algorithm, dict_sample):
    '''
    Write a .txt file needed for MOOSE simulation.

    The variable c is transmitted to the MOOSE simulation.

        Input :
            an algorithm dictionnary (a dict)
            an sample dictionnary (a dict)
        Output :
            Nothing but a .txt file is generated (a file)
    '''
    file_to_write = open('Data/c_'+str(dict_algorithm['i_PFDEM'])+'.txt','w')
    file_to_write.write('AXIS X\n')
    line = ''
    for x in dict_sample['x_L']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write.write(line)

    file_to_write.write('AXIS Y\n')
    line = ''
    for y in dict_sample['y_L']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write.write(line)

    file_to_write.write('DATA\n')
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            file_to_write.write(str(dict_sample['solute_M'][-1-l][c])+'\n')

    file_to_write.close()

#-------------------------------------------------------------------------------

def Write_ep_txt(dict_algorithm, dict_sample):
    '''
    Write a .txt file needed for MOOSE simulation.

    The variable ep is transmitted to the MOOSE simulation.
    This variable is the dissolution field due to mechanical energy at the contact level.

        Input :
            an algorithm dictionnary (a dict)
            an sample dictionnary (a dict)
        Output :
            Nothing but a .txt file is generated (a file)
    '''
    #compute the variable e_mec
    e_mec = 0.2*dict_sample['S_int_0']/dict_sample['S_int']

    #write data
    file_to_write = open('Data/ep_'+str(dict_algorithm['i_PFDEM'])+'.txt','w')
    file_to_write.write('AXIS X\n')
    line = ''
    for x in dict_sample['x_L']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write.write(line)

    file_to_write.write('AXIS Y\n')
    line = ''
    for y in dict_sample['y_L']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write.write(line)

    file_to_write.write('DATA\n')
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):
            file_to_write.write(str(e_mec*min(dict_sample['L_g'][0].etai_M[-1-l][c],dict_sample['L_g'][1].etai_M[-1-l][c]))+'\n')

    file_to_write.close()

#-------------------------------------------------------------------------------

def Write_kc_txt(dict_algorithm, dict_material, dict_sample):
    '''
    Write a .txt file needed for MOOSE simulation.

    The variable kc is transmitted to the MOOSE simulation.
    This variable is the diffusion coefficient of the solute.
    It takes the value 0 if the point is inside one grain and not in the other.
    Else it takes an user defined value.

        Input :
            an algorithm dictionnary (a dict)
            an material dictionnary (a dict)
            an sample dictionnary (a dict)
        Output :
            Nothing but a .txt file is generated (a file)
    '''

    file_to_write = open('Data/kc_'+str(dict_algorithm['i_PFDEM'])+'.txt','w')
    file_to_write.write('AXIS X\n')
    line = ''
    for x in dict_sample['x_L']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write.write(line)

    file_to_write.write('AXIS Y\n')
    line = ''
    for y in dict_sample['y_L']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write.write(line)

    file_to_write.write('DATA\n')
    for l in range(len(dict_sample['y_L'])):
        for c in range(len(dict_sample['x_L'])):

            #inside g1 and not g2
            if dict_sample['L_g'][0].etai_M[-1-l][c] > 0.9 and dict_sample['L_g'][1].etai_M[-1-l][c] < 0.1:
                file_to_write.write('0\n')
            #inside g2 and not g1
            elif dict_sample['L_g'][0].etai_M[-1-l][c] < 0.1 and dict_sample['L_g'][1].etai_M[-1-l][c] > 0.9:
                file_to_write.write('0\n')
            #at the contact or outside of grains
            else:
                file_to_write.write(str(dict_material['kappa_c'])+'\n')

    file_to_write.close()

#---------------------------------------------------------------------------

def solute_PFtoDEM_Multi(FileToRead,dict_algorithm,dict_sample):
    '''
    Read file from MOOSE simulation to reconstruct the phase field of the solute.

        Input :
            the name of the file to read (a string)
            an algorithm dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
        Output :
            Nothing but the sample dictionnary gets an updated attribute (a n_y x n_x numpy array)
    '''

    #---------------------------------------------------------------------------
    #Global parameters
    #---------------------------------------------------------------------------

    dict_sample['solute_M'] = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))

    id_L = None
    c_selector_len = len('        <DataArray type="Float64" Name="c')
    end_len = len('        </DataArray>')
    XYZ_selector_len = len('        <DataArray type="Float64" Name="Points"')
    data_jump_len = len('          ')

    for i_proc in range(dict_algorithm['np_proc']):

        L_Work = [[], #X
                  [], #Y
                  []] #c

    #---------------------------------------------------------------------------
    #Reading file
    #---------------------------------------------------------------------------

        f = open(f'{FileToRead}_{i_proc}.vtu','r')
        data = f.read()
        f.close
        lines = data.splitlines()

        #iterations on line
        for line in lines:

            if line[0:c_selector_len] == '        <DataArray type="Float64" Name="c':
                id_L = 2

            elif line[0:XYZ_selector_len] == '        <DataArray type="Float64" Name="Points"':
                id_L = 0

            elif (line[0:end_len] == '        </DataArray>' or  line[0:len('          <InformationKey')] == '          <InformationKey') and id_L != None:
                id_L = None

            elif line[0:data_jump_len] == '          ' and id_L == 2: #Read c
                line = line[data_jump_len:]
                c_start = 0
                for c_i in range(0,len(line)):
                    if line[c_i]==' ':
                        c_end = c_i
                        L_Work[id_L].append(float(line[c_start:c_end]))
                        c_start = c_i+1
                L_Work[id_L].append(float(line[c_start:]))

            elif line[0:data_jump_len] == '          ' and id_L == 0: #Read [X, Y, Z]
                line = line[data_jump_len:]
                XYZ_temp = []
                c_start = 0
                for c_i in range(0,len(line)):
                    if line[c_i]==' ':
                        c_end = c_i
                        XYZ_temp.append(float(line[c_start:c_end]))
                        if len(XYZ_temp)==3:
                            L_Work[0].append(XYZ_temp[0])
                            L_Work[1].append(XYZ_temp[1])
                            XYZ_temp = []
                        c_start = c_i+1
                XYZ_temp.append(float(line[c_start:]))
                L_Work[0].append(XYZ_temp[0])
                L_Work[1].append(XYZ_temp[1])

        #Adaptating data and update of solute_M
        for i in range(len(L_Work[0])):
            #Interpolation method
            L_dy = []
            for y_i in dict_sample['y_L'] :
                L_dy.append(abs(y_i - L_Work[1][i]))
            L_dx = []
            for x_i in dict_sample['x_L'] :
                L_dx.append(abs(x_i - L_Work[0][i]))
            dict_sample['solute_M'][-1-list(L_dy).index(min(L_dy))][list(L_dx).index(min(L_dx))] = L_Work[2][i]


#-------------------------------------------------------------------------------

def save_dicts_tempo(dict_algorithm, dict_material, dict_sample, dict_sollicitation):
    '''
    Save dictionnaries at the end of each PFDEM interations.

        Input :
            an algorithm dictionnary (a dictionnary)
            a material dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
            a sollicitation dictionnary (a dictionnary)
        Output :
            Nothing but a save file is generated (a file)
    '''
    outfile = open('../'+dict_algorithm['foldername']+'/'+dict_algorithm['namefile']+'_save_tempo','wb')
    dict_save = {}
    dict_save['algorithm'] = dict_algorithm
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitation'] = dict_sollicitation
    pickle.dump(dict_save,outfile)
    outfile.close()

#-------------------------------------------------------------------------------

def save_dicts_final(dict_algorithm, dict_material, dict_sample, dict_sollicitation):
    '''
    Save dictionnaries at the end of the simulation.

        Input :
            an algorithm dictionnary (a dictionnary)
            a material dictionnary (a dictionnary)
            a sample dictionnary (a dictionnary)
            a sollicitation dictionnary (a dictionnary)
        Output :
            Nothing but a save file is generated (a file)
    '''
    os.remove('../'+dict_algorithm['foldername']+'/'+dict_algorithm['namefile']+'_save_tempo')
    outfile = open('../'+dict_algorithm['foldername']+'/'+dict_algorithm['namefile']+'_save','wb')
    dict_save = {}
    dict_save['algorithm'] = dict_algorithm
    dict_save['material'] = dict_material
    dict_save['sample'] = dict_sample
    dict_save['sollicitation'] = dict_sollicitation
    pickle.dump(dict_save,outfile)
    outfile.close()
