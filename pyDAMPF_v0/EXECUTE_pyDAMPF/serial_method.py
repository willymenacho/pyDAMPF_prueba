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

import os
import shutil 
import numpy as np
import re
################argumento####################################
tempall = np.genfromtxt('tempall.txt')
tcases = len(tempall[:,0])
factor=1

###############funciones#######################################
def gen_limites(tcases, factor):
    paralleliz = int(tcases/factor)
    aux = tcases - (paralleliz*factor) #casos sobrantes
    lim_inf = 0
    lim_sup = lim_inf + (paralleliz)
    if aux>0:
        lim_sup = lim_sup + 1
        aux = aux -1
    lst = [(lim_inf, lim_sup)]
    for i in range(1,factor):
        lim_inf = lim_sup
        lim_sup = lim_inf + (paralleliz)
        if aux>0:
            lim_sup = lim_sup + 1
            aux = aux - 1
        lst.append((lim_inf, lim_sup))
    return lst

def change_dir():
    my_path = os.getcwd()
    ls = os.listdir(my_path)
    if "SERIALBASIC_0" in ls:
        r = re.compile("SERIALBASIC_[0-9]")
        lst = list(filter(r.match, ls))
        nums = []
        lst.sort(reverse = True)
        for i in lst:
            aux = i.split("_")[1]
            nums.append(int(aux))
            new_dir = "SERIALBASIC_" + str(int(aux) + 1)
            os.rename(i, new_dir)
        os.makedirs("SERIALBASIC_0")
    else:
        os.makedirs("SERIALBASIC_0")
    
    
    
#################creacion de carpetas, copia de contenido y cambio paralelize#######
def serial_method (tcases , factor, tempall):
	change_dir()
	for i in range (1,factor+1):
	    direc = os.getcwd()
	    direc2 = direc+'/pyDAMPF_BASE/'
	    direc3 = direc+'/SERIALBASIC_0/'+str(i)+'/'
	    shutil.copytree(direc2, direc3)
	os.chdir(direc+'/SERIALBASIC_0/1/nrun/')  
	exec(open("generate_cases.py").read()) 
	 
if __name__ == '__main__':
    serial_method (tcases , factor, tempall)
