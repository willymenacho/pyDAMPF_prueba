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
# IMPORT MODULES
import os
import shutil
# IMPORT NEWLY VISUALIZATION MODULES
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

#======================Prev Definitions=======
tcases =10
aspcases=tcases*4
dminS=np.zeros(aspcases)
fmaxS=np.zeros(aspcases)
nfound=0
ki=0
aca=0
aux=99999.99

#------------------Copy to Inputs template--------------------------
# Open
map0 = np.genfromtxt('tempall.txt')
map0 = map0[0:tcases]

#creamos un array con los Asp 
a, b = np.shape(map0)
arr_aux = np.array([0])
Asp = [0.6500, 0.7500, 0.8500, 0.9500]

#funcion que permite permutar una sola lista######
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

#nuevo txt temporal con Asp
final_arr= permutador(Asp, map0)
namefile = 'pmPosAsp.txt'
np.savetxt(namefile,final_arr) 
print('se ha guardado pmPosAsp.txt')
id0 = final_arr[:,0]
r = final_arr[:,1]
q = final_arr[:,2]
es = final_arr[:,3]
a0 = final_arr[:,4]
k = final_arr[:,5]
fo = final_arr[:,6]
cap = final_arr[:,7]
vis = final_arr[:,8]
ham = final_arr[:,9]
asp = final_arr[:,10]

#for i in range(0,tcases):
for i in range(0,3):     ####borrar depues#### prueba!!!!!!!!
    numplot=(i+1)
    idout='zcdom'+str(i+1)+'.dfo'
    #idout='zcdom1.dfo'
    zcd = np.genfromtxt(idout)
    amp1d = zcd[:,1]
    fase1d = zcd[:,2]
    dmind = zcd[:,3]
    fmaxd = zcd[:,5]
    
    #-------Once gathered keep on seeking---------
    amp1r=amp1d.round(3)
    
    #a0w1=max(amp1d)
    #w#a0w1=round(a0w1)
    #aspv1=[ 0.600*a0w1, 0.700*a0w1, 0.800*a0w1, 0.900*a0w1 ]
    a0w=int(a0[0])
    aspv=[]
    for i in Asp:
        aspv.append(i*a0w)
    fmaxR=[]
    dminR=[]
    for j in range(len(aspv)):
        asp=round(aspv[j],3)
        #print (asp)
        ind=np.nonzero(amp1r<asp)
        #print (ind)
        print ('el min indice que tenemos es:')
        reind=ind[0]
        p=np.size(ind)
        ind=np.resize(reind,(1,1))
        print (p)
        if p!=0:
            #print (dmind[ind])
            #ja=(i+1)*4
            dminS[ki]=dmind[ind]
            fmaxS[ki]=fmaxd[ind]
            #dminS[0:ja]=dmind[ind]
            #fmaxS[0:ja]=fmaxd[ind]
           #jb=ja
        else:
            nfound=nfound+1
            dminS[ki]=aux
            fmaxS[ki]=aux
            ki=ki+1
        dminR.append(dminS[ki])
        fmaxR.append(fmaxS[ki])
        
    ########guardamos los dmin y fmax de cada Asp en un txt#######
    
    #creamos una lista de los calores dmin fmax a cierto Asp
    big_lst = [Asp, dminR, fmaxR]
    
    #funcion que permite permutar listas multiples#####
    def permu_multi(lst, data):
        y, x = np.shape(data)
        y_orig = y
        arr_aux = np.array([0])
        lst_arr = []
        len_lst = len(lst)
        big_add = []
        for i in range(len(lst[0])):
            lst_add = []
            for j in lst:
                lst_add.append(j[i])
            big_add.append(lst_add)
        for i in range(y):
            arr_aux = data[i]
            arr_tmp = np.insert(arr_aux, x, big_add[0])
            arr_family = arr_tmp
            for j in big_add[1:]:
                arr_tmp = np.insert(arr_aux, x, j)
                arr_family = np.concatenate((arr_family, arr_tmp))
            lst_arr.append(arr_family)
        final_arr = np.concatenate(lst_arr)
        y_new = y_orig*len(lst[0])
        x_new = x + len_lst
        final_arr = np.reshape(final_arr, (y_new, x_new))
        y_final,x_final=np.shape(final_arr)
        final_arr[:,0]=np.arange(1,y_final + 1)
        return final_arr
    
    final = permu_multi(big_lst, map0)
    namefile = 'pmPosAspFmaxDmin.txt'
    np.savetxt(namefile,final) 
    print('se ha guardado pmPosAspFmaxDmin.txt')
    
    ##### graficas Asp #####
    sid = map0[:,0]
    fmaxA = final[:,12]
    dminA = final[:,11]
    aspA = final[:,10]
    
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
    
    casosgraf = len(final[:,0])
    num =len(Asp)
    parametro = int(casosgraf/num)
    print(parametro)
    lstR=gen_limites(casosgraf,tcases)
    
    #for i in range (0,parametro):
    #    print(i)
    lim_inferior, lim_superior = lstR[0][0], lstR[0][1]
    #print(aspA[lim_inferior:lim_superior])
    #print(fmaxA[lim_inferior:lim_superior]) 
    #def fmaxasp():
    #========Amp 1 vs. zc=================
    plt.figure(1)
    plt.grid(True)
    plt.title('$Fmax$($Asp$) zcdom'+str(numplot))
    plt.ylabel('$Fmax \; [nN]$',weight='bold',size='x-large')
    plt.xlabel('$Asp\; [nm]$',weight='bold',size='x-large')
    plt.plot(aspA[lim_inferior:lim_superior],fmaxA[lim_inferior:lim_superior],'r.--')
    plt.savefig('fmaxVasp_zcdom'+str(numplot)+'.png')
    plt.show()
    plt.close()
    
    #def dminasp():
    #========Amp 1 vs. zc=================
    plt.figure(2)
    plt.grid(True)
    plt.title('$dmin$($Asp$) zcdom'+str(numplot))
    plt.ylabel('$dmin \; [nm]$',weight='bold',size='x-large')
    plt.xlabel('$Asp\; [nm]$',weight='bold',size='x-large')
    plt.plot(aspA[lim_inferior:lim_superior],dminA[lim_inferior:lim_superior],'r.--')
    plt.savefig('dminVasp_zcdom'+str(numplot)+'.png')
    plt.show()
    plt.close()
print('thats all')
