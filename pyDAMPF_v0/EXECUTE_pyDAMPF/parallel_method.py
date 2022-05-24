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

from sys import argv
import os
import shutil 
from tempfile import mkstemp
import numpy as np
import re
################argumento####################################
script, factor = argv
tempall = np.genfromtxt('tempall.txt')
tcases = len(tempall[:,0])
factor = int(factor)

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

def reemplazo(source_file_path, pattern, substring):
    fh, target_file_path = mkstemp()
    with open(target_file_path, 'w') as target_file:
        with open(source_file_path, 'r') as source_file:
            for line in source_file:
                target_file.write(line.replace(pattern, substring))
    #print(source_file_path)
    os.remove(source_file_path)
    shutil.move(target_file_path, source_file_path)
    
def change_dir():
    my_path = os.getcwd()
    ls = os.listdir(my_path)
    if "PARALLELBASIC_0" in ls:
        r = re.compile("PARALLELBASIC_[0-9]")
        lst = list(filter(r.match, ls))
        nums = []
        lst.sort(reverse = True)
        for i in lst:
            aux = i.split("_")[1]
            nums.append(int(aux))
            new_dir = "PARALLELBASIC_" + str(int(aux) + 1)
            os.rename(i, new_dir)
        os.makedirs("PARALLELBASIC_0")
    else:
        os.makedirs("PARALLELBASIC_0")
        
 #################creacion de carpetas, copia de contenido y cambio paralelize#######
def parallel_method(tcases, factor,tempall):
    lst = gen_limites(tcases, factor)
    change_dir()
    for i in range (1,factor+1):
        lim_inferior, lim_superior = lst[i-1][0], lst[i-1][1]
        direc = os.getcwd()
        direc2 = direc+'/pyDAMPF_BASE/'
        direc3 = direc+'/PARALLELBASIC_0/'+str(i)+'/'
        shutil.copytree(direc2, direc3)
        factorantiguo='factor=1'
        factornuevo='factor='+str(factor)
        rangoantiguo='(0,paraleliz)'
        rangonuevo='('+str(lim_inferior)+","+str(lim_superior)+')'
        casoantiguo= 'tcases =11'
        casonuevo= 'tcases ='+(str(lim_superior-lim_inferior))
        map0antiguo='[0:tcases]'
        map0nuevo='['+str(lim_inferior)+":"+str(lim_superior)+']'
        rangoposantiguo='(0,tcases)'
        os.chdir(direc+'/PARALLELBASIC_0/'+str(i))
        pyname='nrun/generate_cases.py'
        posname='nrun/runa/pmPos.py'
        newpath=direc+'/PARALLELBASIC_0/'+str(i)+'/'+pyname
        reemplazo(newpath, factorantiguo, factornuevo)
        reemplazo(newpath, rangoantiguo, rangonuevo)
        os.chdir(direc)
        newpathpos=direc+'/PARALLELBASIC_0/'+str(i)+'/'+posname
        reemplazo(newpathpos, casoantiguo, casonuevo)
        reemplazo(newpathpos, map0antiguo, map0nuevo)
        reemplazo(newpathpos, rangoposantiguo, rangonuevo)
        os.chdir(direc)
    print(str(factor)+'PARALLEL CASE FOLDERS HAVE BEEN CREATED')
    
if __name__ == '__main__':
    parallel_method(tcases, factor,tempall)
    
