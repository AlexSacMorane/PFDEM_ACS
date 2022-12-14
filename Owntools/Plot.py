# -*- coding: utf-8 -*-
"""
@author: Alexandre Sac--Morane
alexandre.sac-morane@uclouvain.be

This file contains the different functions to plot used in the simulation.
"""

#-------------------------------------------------------------------------------
#Librairy
#-------------------------------------------------------------------------------

from pathlib import Path
from scipy.ndimage import binary_dilation
import numpy as np
import matplotlib.pyplot as plt
import os
import imageio
from Owntools import index_to_str

#-------------------------------------------------------------------------------

def Plot_Diffusion_Solute(dict_algorithm, dict_material, dict_sample):
    '''
    Plot the delta of solute concentration in the case of a pure diffusion problem with the same initial configuration as phase field simulation.

        Input :
            an algorithm dictionnary (a dict)
            a sample dictionnary (a dict)
        Output :
            Nothing but a .png file is generated (a file)
    '''
    #---------------------------------------------------------------------------
    #create .i
    #---------------------------------------------------------------------------
    file_to_write = open('Debug_Diff_Solute_'+str(dict_algorithm['i_PFDEM'])+'.i','w')
    file_to_read = open('Owntools/Debug_Diff_Solute_base.i','r')
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
        elif j == 40:
          line = line[:-1] + ' ' + str(dict_material['M'])+'\n'
        elif j == 53 or j == 57:
          line = line[:-1] + str(dict_algorithm['i_PFDEM']) + '.txt\n'
        elif j == 87:
          line = line[:-1] + ' ' + str(dict_algorithm['dt_PF']*dict_algorithm['n_t_PF']) +'\n'
        elif j == 91:
          line = line[:-1] + ' ' + str(dict_algorithm['dt_PF']) +'\n'
        file_to_write.write(line)

    file_to_write.close()

    #---------------------------------------------------------------------------
    #run simulation
    #---------------------------------------------------------------------------

    os.system('mpiexec -n '+str(dict_algorithm['np_proc'])+' ~/projects/moose/modules/combined/combined-opt -i '+'Debug_Diff_Solute_'+str(dict_algorithm['i_PFDEM'])+'.i')

    #---------------------------------------------------------------------------
    #sort files
    #---------------------------------------------------------------------------

    os.rename('Debug_Diff_Solute_'+str(dict_algorithm['i_PFDEM'])+'_out.e','Debug/Diff_Solute/Ite_'+str(dict_algorithm['i_PFDEM'])+'/Debug_Diff_Solute_'+str(dict_algorithm['i_PFDEM'])+'_out.e')
    os.rename('Debug_Diff_Solute_'+str(dict_algorithm['i_PFDEM'])+'.i','Debug/Diff_Solute/Ite_'+str(dict_algorithm['i_PFDEM'])+'/Debug_Diff_Solute_'+str(dict_algorithm['i_PFDEM'])+'.i')
    j = 0
    j_str = index_to_str(j)
    filepath = Path('Debug_Diff_Solute_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'.pvtu')
    while filepath.exists():
        for i_proc in range(dict_algorithm['np_proc']):
            os.rename('Debug_Diff_Solute_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'_'+str(i_proc)+'.vtu','Debug/Diff_Solute/Ite_'+str(dict_algorithm['i_PFDEM'])+'/Debug_Diff_Solute_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'_'+str(i_proc)+'.vtu')
        os.rename('Debug_Diff_Solute_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'.pvtu','Debug/Diff_Solute/Ite_'+str(dict_algorithm['i_PFDEM'])+'/Debug_Diff_Solute_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'.pvtu')
        j = j + 1
        j_str = index_to_str(j)
        filepath = Path('Debug_Diff_Solute_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'.pvtu')
    j_str = index_to_str(j-1)

    #---------------------------------------------------------------------------
    #read files
    #---------------------------------------------------------------------------

    solute_diff_M = np.array(np.zeros((len(dict_sample['y_L']),len(dict_sample['x_L']))))

    id_L = None
    c_selector_len = len('        <DataArray type="Float64" Name="c')
    end_len = len('        </DataArray>')
    XYZ_selector_len = len('        <DataArray type="Float64" Name="Points"')
    data_jump_len = len('          ')

    for i_proc in range(dict_algorithm['np_proc']):

        L_Work = [[], #X
                  [], #Y
                  []] #c

        f = open('Debug/Diff_Solute/Ite_'+str(dict_algorithm['i_PFDEM'])+'/Debug_Diff_Solute_'+str(dict_algorithm['i_PFDEM'])+'_other_'+j_str+'_'+str(i_proc)+'.vtu','r')
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
            solute_diff_M[-1-list(L_dy).index(min(L_dy))][list(L_dx).index(min(L_dx))] = L_Work[2][i]

    #---------------------------------------------------------------------------
    #Compare with initial value
    #---------------------------------------------------------------------------

    delta_solute = solute_diff_M - dict_sample['solute_M']

    #---------------------------------------------------------------------------
    #Plot result
    #---------------------------------------------------------------------------

    #look for the name of the new plot
    template_name = 'Debug/Diff_Solute/Diffusion_solute_'
    j = 0
    plotpath = Path(template_name+str(j)+'.png')
    while plotpath.exists():
        j = j + 1
        plotpath = Path(template_name+str(j)+'.png')
    name = template_name+str(j)+'.png'

    #plot
    plt.figure(1,figsize=(16,9))

    plt.subplot(311)
    #diffusion coefficient
    im = plt.imshow(dict_sample['kc_M'],interpolation='nearest', extent=[min(dict_sample['x_L']),max(dict_sample['x_L']),min(dict_sample['y_L']),max(dict_sample['y_L'])])
    plt.colorbar(im)
    #etai
    for i in range(len(dict_sample['L_g'])):
        plt.plot(dict_sample['L_g'][i].l_border_x,dict_sample['L_g'][i].l_border_y)
    plt.title('Diffusion coefficient solute and grains')
    plt.axis('equal')
    plt.xlim(min(dict_sample['x_L']),max(dict_sample['x_L']))

    plt.subplot(312)
    #diffusion coefficient
    im = plt.imshow(delta_solute,interpolation='nearest', extent=[min(dict_sample['x_L']),max(dict_sample['x_L']),min(dict_sample['y_L']),max(dict_sample['y_L'])])
    plt.colorbar(im)
    #etai
    for i in range(len(dict_sample['L_g'])):
        plt.plot(dict_sample['L_g'][i].l_border_x,dict_sample['L_g'][i].l_border_y)
    plt.title('Delta solute (= diffused - initial)')
    plt.axis('equal')
    plt.xlim(min(dict_sample['x_L']),max(dict_sample['x_L']))

    grad_y, grad_x = np.gradient(-delta_solute)
    plt.subplot(325)
    im = plt.imshow(grad_x,interpolation='nearest', extent=[min(dict_sample['x_L']),max(dict_sample['x_L']),min(dict_sample['y_L']),max(dict_sample['y_L'])])
    plt.colorbar(im)
    #etai
    for i in range(len(dict_sample['L_g'])):
        plt.plot(dict_sample['L_g'][i].l_border_x,dict_sample['L_g'][i].l_border_y)
    plt.title('Grad_x(Delta solute)')
    plt.axis('equal')
    plt.xlim(min(dict_sample['x_L']),max(dict_sample['x_L']))
    plt.subplot(326)
    im = plt.imshow(grad_y,interpolation='nearest', extent=[min(dict_sample['x_L']),max(dict_sample['x_L']),min(dict_sample['y_L']),max(dict_sample['y_L'])])
    plt.colorbar(im)
    #etai
    for i in range(len(dict_sample['L_g'])):
        plt.plot(dict_sample['L_g'][i].l_border_x,dict_sample['L_g'][i].l_border_y)
    plt.title('Grad_y(Delta solute)')
    plt.axis('equal')
    plt.xlim(min(dict_sample['x_L']),max(dict_sample['x_L']))

    plt.savefig(name)
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_config(dict_algorithm, dict_sample):
    '''
    Plot the sample configuration.

        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing but a .png file is generated (a file)
    '''
    #look for the name of the new plot
    template_name = 'Debug/Configuration/Configuration_'
    j = 0
    plotpath = Path(template_name+str(j)+'.png')
    while plotpath.exists():
        j = j + 1
        plotpath = Path(template_name+str(j)+'.png')
    name = template_name+str(j)+'.png'

    #plot
    plt.figure(1,figsize=(16,9))

    plt.subplot(211)
    #solute
    im = plt.imshow(dict_sample['solute_M'],interpolation='nearest', extent=[min(dict_sample['x_L']),max(dict_sample['x_L']),min(dict_sample['y_L']),max(dict_sample['y_L'])], vmin = dict_algorithm['c_min'], vmax = dict_algorithm['c_max'])
    plt.colorbar(im)
    #etai
    for i in range(len(dict_sample['L_g'])):
        plt.plot(dict_sample['L_g'][i].l_border_x,dict_sample['L_g'][i].l_border_y)
    plt.title('Solute c and grains')
    plt.axis('equal')
    plt.xlim(min(dict_sample['x_L']),max(dict_sample['x_L']))

    #eta1_M
    plt.subplot(223)
    im = plt.imshow(dict_sample['L_g'][0].etai_M,interpolation='nearest', extent=[min(dict_sample['x_L']),max(dict_sample['x_L']),min(dict_sample['y_L']),max(dict_sample['y_L'])], vmin = 0, vmax = 1)
    plt.colorbar(im)
    plt.plot(dict_sample['L_g'][0].l_border_x,dict_sample['L_g'][0].l_border_y,'r')
    plt.title(r'$\eta$1')

    #eta2_M
    plt.subplot(224)
    im = plt.imshow(dict_sample['L_g'][1].etai_M,interpolation='nearest', extent=[min(dict_sample['x_L']),max(dict_sample['x_L']),min(dict_sample['y_L']),max(dict_sample['y_L'])], vmin = 0, vmax = 1)
    plt.colorbar(im)
    plt.plot(dict_sample['L_g'][1].l_border_x,dict_sample['L_g'][1].l_border_y,'r')
    plt.title(r'$\eta$2')

    plt.savefig(name)
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_mp4(template_name,name_video):
    '''The goal of this function is to create a movie with pictures.

    from https://www.blog.pythonlibrary.org/2021/06/23/creating-an-animated-gif-with-python/

        Input :
            the template of the pictures used (a string)
            the name of the video (a string)
        Output :
            a movie file (a .mp4)
    '''
    #look for the largest iteration
    i_f = 0
    plotpath = Path(template_name+str(i_f)+'.png')
    while plotpath.exists():
        i_f = i_f + 1
        plotpath = Path(template_name+str(i_f)+'.png')

    fileList = []
    for i in range(0,i_f):
        fileList.append(template_name+str(i)+'.png')

    duration_movie  = 10 #sec
    writer = imageio.get_writer(name_video, fps=int(i_f/duration_movie))
    for im in fileList:
        writer.append_data(imageio.imread(im))
    writer.close()

#-------------------------------------------------------------------------------

def Plot_Ed(dict_sample):
    '''
    Plot the energy source configuration at the start of the simultion.

        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing but a .png file is generated (a file)
    '''
    #look for the name of the new plot
    template_name = 'Debug/Ed/Ed_'
    j = 0
    plotpath = Path(template_name+str(j)+'.png')
    while plotpath.exists():
        j = j + 1
        plotpath = Path(template_name+str(j)+'.png')
    name = template_name+str(j)+'.png'

    #plot
    plt.figure(1,figsize=(16,9))

    #Ed_M
    plt.subplot(211)
    im = plt.imshow(dict_sample['Ed_M'],interpolation='nearest', extent=[min(dict_sample['x_L']),max(dict_sample['x_L']),min(dict_sample['y_L']),max(dict_sample['y_L'])])
    plt.colorbar(im)
    #etai
    for i in range(len(dict_sample['L_g'])):
        plt.plot(dict_sample['L_g'][i].l_border_x,dict_sample['L_g'][i].l_border_y)
    plt.axis('equal')
    plt.xlim(min(dict_sample['x_L']),max(dict_sample['x_L']))
    plt.title('Ed = Emec - Eche')

    #Emec_M
    plt.subplot(223)
    im = plt.imshow(dict_sample['Emec_M'],interpolation='nearest', extent=[min(dict_sample['x_L']),max(dict_sample['x_L']),min(dict_sample['y_L']),max(dict_sample['y_L'])])
    plt.colorbar(im)
    #etai
    for i in range(len(dict_sample['L_g'])):
        plt.plot(dict_sample['L_g'][i].l_border_x,dict_sample['L_g'][i].l_border_y)
    plt.axis('equal')
    plt.xlim(min(dict_sample['x_L']),max(dict_sample['x_L']))
    plt.title('Emec')

    #Eche_M
    plt.subplot(224)
    im = plt.imshow(dict_sample['Eche_M'],interpolation='nearest', extent=[min(dict_sample['x_L']),max(dict_sample['x_L']),min(dict_sample['y_L']),max(dict_sample['y_L'])])
    plt.colorbar(im)
    #etai
    for i in range(len(dict_sample['L_g'])):
        plt.plot(dict_sample['L_g'][i].l_border_x,dict_sample['L_g'][i].l_border_y)
    plt.axis('equal')
    plt.xlim(min(dict_sample['x_L']),max(dict_sample['x_L']))
    plt.title('Eche')

    plt.savefig(name)
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_kc(dict_sample):
    '''
    Plot the distribution of the diffusion coefficient for the solute.

        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing but a .png file is generated (a file)
    '''
    #look for the name of the new plot
    template_name = 'Debug/Kc/Kc_'
    j = 0
    plotpath = Path(template_name+str(j)+'.png')
    while plotpath.exists():
        j = j + 1
        plotpath = Path(template_name+str(j)+'.png')
    name = template_name+str(j)+'.png'

    #plot
    plt.figure(1,figsize=(16,9))
    #kc_M
    im = plt.imshow(dict_sample['kc_M'],interpolation='nearest', extent=[min(dict_sample['x_L']),max(dict_sample['x_L']),min(dict_sample['y_L']),max(dict_sample['y_L'])])
    plt.colorbar(im)
    #plot g1 and g2 boundaries
    plt.plot(dict_sample['L_g'][0].l_border_x,dict_sample['L_g'][0].l_border_y)
    plt.plot(dict_sample['L_g'][1].l_border_x,dict_sample['L_g'][1].l_border_y)

    plt.axis('equal')
    plt.xlim(min(dict_sample['x_L']),max(dict_sample['x_L']))
    plt.title('Diffusion coefficient of the solute')
    plt.savefig(name)
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_init_current_shape(dict_sample):
    '''
    Plot the comparison between initial and current shape for the grain 1.

        Input :
            a sample dictionnary (a dict)
        Output :
            Nothing but a .png file is generated (a file)
    '''
    #look for the name of the new plot
    template_name = 'Debug/Comparison_Init_Current/Init_Current_Shape_'
    j = 0
    plotpath = Path(template_name+str(j)+'.png')
    while plotpath.exists():
        j = j + 1
        plotpath = Path(template_name+str(j)+'.png')
    name = template_name+str(j)+'.png'

    #prepare plot
    L_border_x_init = []
    L_border_y_init = []
    for i in range(len(dict_sample['L_g'][0].l_border_x_init)):
        L_border_x_init.append(dict_sample['L_g'][0].l_border_x_init[i] - dict_sample['L_g'][0].center_init[0])
        L_border_y_init.append(dict_sample['L_g'][0].l_border_y_init[i] - dict_sample['L_g'][0].center_init[1])
    L_border_x = []
    L_border_y = []
    for i in range(len(dict_sample['L_g'][0].l_border_x)):
        L_border_x.append(dict_sample['L_g'][0].l_border_x[i] - dict_sample['L_g'][0].center[0])
        L_border_y.append(dict_sample['L_g'][0].l_border_y[i] - dict_sample['L_g'][0].center[1])
    #plot
    plt.figure(1,figsize=(16,9))
    plt.plot(L_border_x_init,L_border_y_init,'k',label='Initial')
    plt.plot(L_border_x,L_border_y,label='Current')
    plt.savefig(name)
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_sum_eta_c(dict_tracker):
    '''
    Plot the trackers.

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing but a .png files are generated (files)
    '''
    #plot Displacement and sum of c and etai
    plt.figure(1,figsize=(16,9))

    plt.subplot(211)
    plt.plot(dict_tracker['L_t'], dict_tracker['L_int_displacement'])
    plt.title('Total displacement done')

    plt.subplot(234)
    plt.plot(dict_tracker['L_t'], dict_tracker['L_sum_eta'])
    plt.title('Sum of etai')

    plt.subplot(235)
    plt.plot(dict_tracker['L_t'], dict_tracker['L_sum_solute'])
    plt.title('Sum of c')

    plt.subplot(236)
    plt.plot(dict_tracker['L_t'], dict_tracker['L_sum_total'])
    plt.title('Sum of etai and c')

    plt.savefig('Debug/Trackers.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_sphericity(dict_tracker):
    '''
    Plot the sphericity of the grain 1.

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing but a .png files are generated (files)
    '''
    plt.figure(1,figsize=(16,9))

    plt.plot(dict_tracker['L_t'], dict_tracker['L_area_sphericity_g0'],label='Area sphericity')
    plt.plot(dict_tracker['L_t'], dict_tracker['L_diameter_sphericity_g0'],label='Diameter sphericity')
    plt.plot(dict_tracker['L_t'], dict_tracker['L_circle_ratio_sphericity_g0'],label='Circle sphericity')
    plt.plot(dict_tracker['L_t'], dict_tracker['L_perimeter_sphericity_g0'],label='Perimeter sphericity')
    plt.plot(dict_tracker['L_t'], dict_tracker['L_width_to_length_ratio_sphericity_g0'],label='Width to length ratio sphericity')
    plt.legend()
    plt.title('2D sphericity of grain 1')

    plt.savefig('Debug/Sphericity_g_1.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_c_at_p(dict_tracker):
    '''
    Plot the value of the solute concentration at the point defined.

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing but a .png files are generated (files)
    '''
    plt.figure(1,figsize=(16,9))

    plt.plot(dict_tracker['L_t'], dict_tracker['c_at_the_center'])
    plt.title('Value of the solute concentration at the center')
    plt.savefig('Debug/Solute_Contration_Center.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_sum_Ed(dict_tracker):
    '''
    Plot the value of the total energy in the sample.

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing but a .png files are generated (files)
    '''
    plt.figure(1,figsize=(16,9))

    plt.subplot(331)
    plt.plot(dict_tracker['L_t'][:-1], dict_tracker['sum_Ed_che_L'])
    plt.title('Total chemical energy Ed_che')

    plt.subplot(334)
    plt.plot(dict_tracker['L_t'][:-1], dict_tracker['sum_Ed_mec_L'])
    plt.title('Total mechanical energy Ed_mec')

    plt.subplot(337)
    plt.plot(dict_tracker['L_t'][:-1], dict_tracker['sum_ed_L'])
    plt.title('Total Energy Ed = Ed_mec - Ed_che')

    plt.subplot(132)
    plt.plot(dict_tracker['L_t'][:-1], dict_tracker['sum_ed_plus_L'], label = 'Ed+')
    plt.plot(dict_tracker['L_t'][:-1], dict_tracker['sum_ed_minus_L'], label = 'Ed-')
    plt.title('Repartition of the energy in a + and a -  terms')

    plt.subplot(133)
    plt.plot(dict_tracker['L_t'], dict_tracker['L_int_displacement'])
    plt.title('Total displacement done')

    plt.savefig('Debug/Evolution_sum_Ed.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_Sint_SumMinEtai(dict_tracker):
    '''
    Plot the evolution of the intersection surface and the sum of min of etai.

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing but a .png files are generated (files)
    '''
    plt.figure(1,figsize=(16,9))

    plt.subplot(211)
    plt.plot(dict_tracker['L_t'][:-1], dict_tracker['S_int_L'])
    plt.title('Intersection surface')

    plt.subplot(212)
    plt.plot(dict_tracker['L_t'][:-1], dict_tracker['sum_min_etai_L'])
    plt.title('Sum of minimum etai')

    plt.savefig('Debug/Evolution_Sint_SumMinEtai.png')
    plt.close(1)

#-------------------------------------------------------------------------------

def Plot_dt_used(dict_tracker):
    '''
    Plot the evolution of the time step used in the phase field simulation.

        Input :
            a tracker dictionnary (a dict)
        Output :
            Nothing but a .png files are generated (files)
    '''
    plt.figure(1,figsize=(16,9))
    plt.plot(dict_tracker['L_dt'])
    plt.title('Evolution of the time step used')
    plt.savefig('Debug/Evolution_dt_used.png')
    plt.close(1)
