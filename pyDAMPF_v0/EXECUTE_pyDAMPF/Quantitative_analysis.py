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

import pandas as pd
import numpy as np
import os
from IPython.display import display


tempall = np.genfromtxt(open('tempall.txt',"r"))
model = pd.read_csv("cantileversdata_model.cvs")
casos = np.array(len(tempall[:,0]))


def findfile(name, path):
    for dirpath, dirname, filename in os.walk(path):
        if name in filename:
            return os.path.join(dirpath, name)
#dominio temporal

def cauntitative (tempall, model):
    direc = os.getcwd()
    direc = direc+'/PARALLELBASIC_0/'
    forcet=[];maxforcet=[];vdwt=[];hertzt=[];viscot=[];capt=[];ljt=[];dlvot=[]
         
    for i in range(casos):
        filepatht = findfile("tdom"+str(i+1)+".dfo", direc)
        tdom = np.genfromtxt(open(str(filepatht),"r"))        
        forcet.append(tdom[:,4])
        maxforcet.append(tdom[:,5])
        vdwt.append(tdom[:,6])
        hertzt.append(tdom[:,7])
        viscot.append(tdom[:,8])
        capt.append(tdom [:,9]) 
        ljt.append(tdom[:,10])
        dlvot.append(tdom [:,11])    
        
    forcetM=[];vdwtM=[];hertztM=[];viscotM=[];captM=[];ljtM=[];dlvotM=[]
    for j in range(casos):
        forcetM.append([max(forcet[j])]) 
        vdwtM.append([max(vdwt[j])]) 
        hertztM.append([max(hertzt[j])]) 
        viscotM.append([max(viscot[j])]) 
        captM.append([max(capt[j])]) 
        ljtM.append([max(ljt[j])]) 
        dlvotM.append([max(dlvot[j])]) 
    
    df = model.assign(ID=tempall[:,0] , LENGTH=tempall[:,1], WIDTH=tempall[:,2], THICKNESS=tempall[:,3] , F0=tempall[:,4] , R_tip=tempall[:,5] , Kc=tempall[:,6] , VOLUME=tempall[:,8] , A0=tempall[:,15] , RH=tempall[:,29] , Es=tempall[:,16], Q=tempall[:,17], F_T_tT=forcetM , F_Hertz=hertztM , F_vdw=vdwtM , F_vis=viscotM , F_cap=captM , F_LJ=ljtM , F_DLVO=dlvotM )
    
    pd.set_option('display.max_columns', None) 
    print(df)
if __name__ == '__main__':
    cauntitative (tempall, model)