#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#=======================pyDAMPF==============================
#The code written below is part of the pyDAMPF project and has
#been developed as a joint project between the IJS and UMSA. This 
#software will enable students and researchers from the 
#AFM community to easily plan AFM experiments.
#======================CREDITS==============================
#
# PI Horacio V. Guzman Ph.D.
# B.Sc. Willy Menacho N.
#
#We thank you in advance for sending your feedback and/or 
#suggestions to:
#             horacio.v.g@gmail.com
#
#======================pyDAMPF LICENSE==============================

#Copyright (C) 2022  Horacio V. Guzman and Willy Menacho N.

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# What we observe is not nature itself, but nature exposed to
# our method of questioning. 
#                                  W. Heisenberg
############################################################
#====================================================================
import numpy as np
import os
import shutil

def permutador(lst, data):
    y, x = np.shape(data)
    arr_aux = np.array([0])
    lst_arr = []
    for i in range(y):
        arr_aux = data[i]
        arr_tmp = np.insert(arr_aux, x, lst[0])
        arr_family = arr_tmp
        for j in lst[1:]:
            arr_tmp = np.insert(arr_aux, x, j)
            arr_family = np.concatenate((arr_family, arr_tmp))
        arr_family = np.reshape(arr_family, (len(lst), x+1))
        lst_arr.append(arr_family)
    final_arr = np.concatenate(lst_arr)
    y, x = np.shape(final_arr)
    final_arr[:,0] = np.arange(1,y+1)
    return final_arr

def gran_permutador(lst, data2):
    aux = permutador(lst[0], data2)
    for i in lst[1:]:
        aux = permutador(i, aux)
    y, x = np.shape(aux)
    aux[:,0] = np.arange(1,y+1)
    return aux


#### values for the numerical integrator
nper = [float(input('Number of periods  ' ))]
npp = [float(input('Number of points per period  ' ))]
nperfin = [float(input('Interval of periods  '))]
zcmin = [float(input('Minimum distance [nm]  '))]
dzc = [float(input('Average tip–surface distance step [nm]  '))]
fixzc = [float(input('Distance fixed [nm]  '))]

#### essential values 
a0 = [float(input('Amplitude [nm]  '))]
Es = [float(input('Young modulus sample [MPa]  '))]
q = [float(input('Quality factor [Adim]  '))]
RH = [float(input('Relative Humidity [%] ' ))]

#### values for interactive forces
ham = [float(input('Hamaker Const [*10^-20 J]  '  ))]
eta = [float(input('Viscosity [Pa s]  '  ))]
epsilon2 = [float(input('Surface energy of solid-liquid  [mN/m]  ' ))]
epsilon1 = [float(input('LR Disp. [Adim]  ' ))]
eps = [float(input('relative permeability [C^2/N m^2]  ' ))]
eps0 = [float(input('vacuum permeability [C^2/N m^2]  ' ))] 
sigmat = [float(input('tip-surface loading density [C/m^2]  ' ))] 
sigmas = [float(input('sample surface charge density [C/m^2]  ' ))]
debye = [float(input('Debye lengths [nm]  ' ))]
ljmin = [float(input('Epsilon force LJ [nN]  ' ))]
lenght = [float(input('Sigma force LJ [nm]  ' ))]
#mtip = [float(input('Tip magn [A m2] ' ))]
#msample = [float(input('Surface magn [A m2] ' ))]  

data = np.genfromtxt("cantilevers_data.txt")

######IF THIS VALUES ARE 0 = INDETERMINATION ##### FOR FIX THIS:
if debye == [0]:
    debye = [1]
if eps == [0]:
    eps = [1]
if eps0 == [0]:
    eps0 = [1]
        
def inputs_processor (nper,npp,nperfin,zcmin,dzc,fixzc,a0,Es,q,ham,eta,epsilon2,epsilon1,eps,eps0,sigmat,sigmas,debye,ljmin,lenght,RH,data):
    a, b = np.shape(data)
    #list
    big_lst = [nper,npp,nperfin,zcmin,dzc,fixzc,a0,Es,q,ham,eta,epsilon2,epsilon1,eps,eps0,sigmat,sigmas,debye,ljmin,lenght,RH]
    final = gran_permutador(big_lst, data)
    #### SAVE DATA ####
    f_name = "tempall.txt"
    np.savetxt(f_name, final)
    print('THE FILE HAS BEEN CREATED WITH ALL THE DATA AND VARIABLES!')
    directory = os.getcwd()
    shutil.copy(directory+'/tempall.txt' , directory+'/EXECUTE_pyDAMPF/')
    shutil.copy(directory+'/tempall.txt' , directory+'/EXECUTE_pyDAMPF/pyDAMPF_BASE/nrun/')
    shutil.copy(directory+'/tempall.txt' , directory+'/EXECUTE_pyDAMPF/pyDAMPF_BASE/nrun/runa/')
    
if __name__ == '__main__':
    inputs_processor()
    