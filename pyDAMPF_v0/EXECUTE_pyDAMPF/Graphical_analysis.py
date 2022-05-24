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
import sys
import numpy as np
from sys import argv
import matplotlib.pyplot as plt


####### with argv ########
'''
#for casos_a in sys.argv:
#    casos_a = []

#script, escenario ,*casos_a = argv
#escenario,*casos_a =input('prueba').split()
escenario,*casos_a =input('prueba')
casos_a1 = str(casos_a)
casos_a2 = np.array(casos_a)
esc = escenario
#esc = str(escenario)
#inputs = int(casos_a)
#xcasos = map(float, inputs.strip('[]').split())
xcasos = map(float, casos_a1.strip('[]').split())
#xcasos=input('Which cases would you like to compare?  ').split()
#xcasos = np.array([int(x) for x in xcasos])
xcasos = np.array([int(x) for x in casos_a])
'''

#### with inputs####
casos_a =input('Which cases would you like to compare?  ').split()
escenario = input ('Serial method (s) or Parallel method (p)  ')


def findfile(name, path):
    for dirpath, dirname, filename in os.walk(path):
        if name in filename:
            return os.path.join(dirpath, name)
        
def plotter (casos_a,escenario):
    esc = escenario
    xcasos = np.array([int(x) for x in casos_a])       
    direc = os.getcwd()
    if esc=='p' or esc == 'P'or esc == '-P' or esc == '-p' :
        direc = direc+'/PARALLELBASIC_0/'
    elif esc == 's' or esc == 'S'or esc == '-s'or esc == '-S':
        direc = direc+'/SERIALBASIC_0/'
    else:
        print('Select parallel mode (-p) or serial mode (-s) in the first argument' )
        sys.exit()
    #print(direc)  
      
    
    
    #dominio temporal
    tt=[]; zt=[];vt=[];forcet=[];maxforcet=[];vdwt=[];hertzt=[];viscot=[];capt=[];ljt=[];dlvot=[]
    #dominio espacial
    zc=[]; amp1d=[];fase1d=[];dmind=[];dmaxd=[];fmaxd=[];defld=[];ets1d=[];vts1d=[];pts1d=[];fmedd=[];tcd=[]        
    for i in xcasos:
        filepatht = findfile("tdom"+str(i)+".dfo", direc)
        #print(filepatht)
        filepathzc = findfile("zcdom"+str(i)+".dfo", direc)
        #print(filepathzc)       
        tdom = np.genfromtxt(open(str(filepatht),"r"))        
        tt.append(tdom[:,0])
        zt.append(tdom[:,1])
        vt.append(tdom[:,2])
        forcet.append(tdom[:,4])
        maxforcet.append(tdom[:,5])
        vdwt.append(tdom[:,6])
        hertzt.append(tdom[:,7])
        viscot.append(tdom[:,8])
        capt.append(tdom [:,9]) 
        ljt.append(tdom [:,10])
        dlvot.append(tdom [:,11])       
        zcdom = np.genfromtxt(open(filepathzc,"r"))
        zc.append(zcdom[:,0])
        amp1d.append(zcdom[:,1])
        fase1d.append(zcdom[:,2])
        dmind.append(zcdom[:,3])
        dmaxd.append(zcdom[:,4])
        fmaxd.append(zcdom[:,5])
        defld.append(zcdom[:,9])
        ets1d.append(zcdom[:,6])
        vts1d.append(zcdom[:,7])
        pts1d.append(zcdom[:,11])
        fmedd.append(zcdom[:,8])
        tcd.append(zcdom[:,10])        
            
            
    Asp = [6.500, 7.500, 8.500, 9.500]
    colores = ['blue', 'orange', 'green', 'red', 'khaki', 'darkviolet', 'brown', 'pink', 'grey', 'gold', 'skyblue', 'cyan', 'darksalmon', 'palegreen', 'darkblue', 'blueviolet','honeydew', 'lightslategray', 'magenta', 'seagreen', 'royalblue', 'deeppink', 'yellowgreen', 'mediumaquamarine', 'rebeccapurple', 'linen', 'mistyrose', 'papayawhip', 'moccasin', 'navy']
    
    #======================PLOT========================
        #||||||||||||||||||||||||||||||||||||||||||
    
    cplots = direc+"/COMPARATIVE_PLOTS/"
    try:
        os.stat(cplots)
    except:
        os.mkdir(cplots)
    os.chdir(cplots)
    
    #def amp1():
    #========Amp 1 vs. zc=================
    plt.figure(1)
    plt.grid(True)
    plt.title('$A_1$($z_c$)')
    plt.ylabel('$A_1 \; [nm]$',weight='bold',size='x-large')
    plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(zc[i],amp1d[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(0.5,max(amp1d[0]))
    plt.legend()
    plt.savefig('amp1Vzc_caso'+str(xcasos)+'.png')
    plt.show()
    
    #def pha1():
    #========Phase 1 vs. zc=================
    plt.figure(2)
    plt.grid(True)
    plt.title('$\phi_1$($z_c$)')
    plt.ylabel('$\phi_1 \; [deg]$',weight='bold',size='x-large')
    plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(zc[i],fase1d[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.legend()
    plt.xlim(0.5,max(amp1d[0]))
    plt.savefig('phi1Vzc_caso'+str(xcasos)+'.png')
    plt.show()
    '''
    #def den1():
    #========Diss Energy 1st vs. zc=================
    plt.figure(3)
    plt.grid(True)
    plt.title('First disp. energy vs. Average distance')
    plt.ylabel('$E_{ts1} * 10^{-20} \; [J]$',weight='bold',size='x-large')
    plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(zc[i],ets1d[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(0.5,max(amp1d[0]))
    plt.legend()
    plt.savefig('ets1Vzc_caso'+str(xcasos)+'.png')
    plt.show()
    
    #def dpo1():
    #========Diss Power 1st vs. zc=================
    plt.figure(4)
    plt.grid(True)
    plt.title('First disp. power vs. Average distance')
    plt.ylabel('$P_{ts1} * 10^{-20} \; [W]$',weight='bold',size='x-large')
    plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(zc[i],pts1d[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(0.5,max(amp1d[0]))
    plt.legend()
    plt.savefig('pts1Vzc_caso'+str(xcasos)+'.png')
    plt.show()
    '''
    #def vir1():
    #========Virial 1st vs. zc=================
    plt.figure(5)
    plt.grid(True)
    plt.title('Virial first vs. Average distance')
    plt.ylabel('$V_{ts1} * 10^{-20} \; [J]$',weight='bold',size='x-large')
    plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(zc[i],vts1d[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(0.5,max(amp1d[0]))
    plt.legend()
    plt.savefig('vts1Vzc_caso'+str(xcasos)+'.png')
    plt.show()
    
    #def def1():
    #========Deflection vs. zc=================
    plt.figure(6)
    plt.grid(True)
    plt.title('Deflection vs. Average distance')
    plt.ylabel('$Defl \; [nm]$',weight='bold',size='x-large')
    plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(zc[i],defld[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(0.5,max(amp1d[0]))
    plt.legend()
    plt.savefig('deflVzc_caso'+str(xcasos)+'.png')
    plt.show()
    
    #def dmin1():
    #========Dmin vs. zc=================
    plt.figure(7)
    plt.grid(True)
    plt.title('Minimum Distance vs. Average distance')
    plt.ylabel('$dmin \; [nm]$',weight='bold',size='x-large')
    plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(zc[i],dmind[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(0.5,max(amp1d[0]))
    plt.legend()
    plt.savefig('dminVzc_caso'+str(xcasos)+'.png')
    plt.show()
    
    #def maf():
    #========Fmax vs. zc=================
    plt.figure(8)
    plt.grid(True)
    plt.title('Maximum Force vs. Average distance')
    plt.ylabel('$F_{max} \; [nN]$',weight='bold',size='x-large')
    plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(zc[i],fmaxd[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(0.5,max(amp1d[0]))
    plt.legend()
    plt.savefig('fmaxVzc_caso'+str(xcasos)+'.png')
    plt.show()
    
    #def tco():
    #========Contact time vs. zc=================
    plt.figure(9)
    plt.grid(True)
    plt.title('Contact Time (tc/T) vs. Average distance (nm)')
    plt.ylabel('$Contact Time \; [t]$',weight='bold',size='x-large')
    plt.xlabel('$z_c\; [nm]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(zc[i],tcd[i], label='Case '+str(xcasos[i]) ,color=colores[(xcasos[i])-1])
    plt.xlim(0.5,max(amp1d[0]))
    plt.legend()
    plt.savefig('ctimeVzc_caso'+str(xcasos)+'.png')
    plt.show()
    
    #----------------------------------------------------------------
    #===================Time Domain==================================
    #----------------------------------------------------------------
    #def ztt():
    #========Instantaneous position vs. time=================
    plt.figure(10)
    plt.grid(True)
    plt.title("$z_T$($t$)")
    plt.ylabel('$z_T \; [nm]$',weight='bold',size='x-large')
    plt.xlabel('$[t/T]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(tt[i],zt[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(max(tt[0])-3,max(tt[0]))
    plt.legend()
    plt.savefig('totalposTime_caso'+str(xcasos)+'.png')
    plt.show()
    
    #def vtt():
    #========Instantaneous velocity vs. time=================
    plt.figure(11)
    plt.grid(True)
    plt.title("$v_T$($t$)")
    plt.ylabel('$v_T \; [nm/\mu s]$',weight='bold',size='x-large')
    plt.xlabel('$[t/T]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(tt[i-1],vt[i-1], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1] )
    plt.xlim(max(tt[0])-3,max(tt[0]))
    plt.legend()
    plt.savefig('totalveloTime_caso'+str(xcasos)+'.png')
    plt.show()
    
    #def ftt():
    #========Instantaneous force vs. time=================
    plt.figure(12)
    plt.grid(True)
    plt.title('$f_T$($t$)')
    plt.ylabel('$f_T \; [nN]$',weight='bold',size='x-large')
    plt.xlabel('$[t/T]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(tt[i],forcet[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(max(tt[0])-3,max(tt[0]))
    plt.legend()
    plt.savefig('totalforceTime_caso'+str(xcasos)+'.png')
    plt.show()
    
    #def ftdt():
    #=========max force vs. distance (ida-vuelta)===============
    plt.figure(13)
    plt.grid(True)
    plt.title('$F_{max} ($d$)\ TOTAL$')
    plt.ylabel('$F_{max} \; [nN]$',weight='bold',size='x-large')
    plt.xlabel('$[d]$',weight='bold',size='x-large')
    limx=[]
    for i in range(len(xcasos)):
        limx.append(min(maxforcet[i][0:1024]))
        plt.plot(maxforcet[i][0:1024],forcet[i][0:1024], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(min(limx),3)
    plt.legend()
    plt.savefig('totalmaxforceTimeidavuelta_caso'+str(xcasos)+'.png')
    plt.show()
    
    #def ftdi():
    #=========max force vs. distance (ida)===============
    plt.figure(14)
    plt.grid(True)
    plt.title('$F_{max}(d)\ one\ way\ $')
    plt.ylabel('$F_{max} \; [nN]$',weight='bold',size='x-large')
    plt.xlabel('$d [nm]$',weight='bold',size='x-large')
    limx=[]
    for i in range(len(xcasos)):
        limx.append(min(maxforcet[i][0:1024]))
        aux=np.where(forcet[i][0:1024]==forcet[i][0:1024].max())
        aux=int(aux[0][0])
        plt.plot(maxforcet[i][0:aux],forcet[i][0:aux], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(min(limx),3)
    plt.legend()
    plt.savefig('TmaxforceTimeida_caso'+str(xcasos)+'.png')
    plt.show()
    
    #def ftdv():
    #=========max force vs. distance (vuelta)===============
    plt.figure(15)
    
    plt.grid(True)
    plt.title('$F_{max}(d)\ way\ back$')
    plt.ylabel('$F_{max} \; [nN]$',weight='bold',size='x-large')
    plt.xlabel('$[d]$',weight='bold',size='x-large')
    limx=[]
    for i in range(len(xcasos)):
        limx.append(min(maxforcet[i][0:1024]))
        aux=np.where(forcet[i][0:1024]==forcet[i][0:1024].max())
        aux=int(aux[0][0])
        plt.plot(maxforcet[i][aux:1024],forcet[i][aux:1024], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(min(limx),3)
    plt.legend()
    plt.savefig('totalmaxforceTimevuelta_caso'+str(xcasos)+'.png')
    plt.show()
    
    #def ftvdw():
    #=========vdw force vs. time============================
    plt.figure(16)
    plt.grid(True)
    plt.title('$F_{VdW}$($d$)')
    plt.ylabel('$F_{VdW} \; [nN]$',weight='bold',size='x-large')
    plt.xlabel('$[t/T]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(tt[i],vdwt[i], label='Case '+str(xcasos[i]) ,color=colores[(xcasos[i])-1])
    plt.xlim(max(tt[0])-3,max(tt[0]))
    plt.legend()
    plt.savefig('vdwforce_caso'+str(xcasos)+'.png')
    plt.show()
    
    #def fthertz():
    #=========hertz force vs. time=========================
    plt.figure(17)
    plt.grid(True)
    plt.title('$F_{Hertz}$($d$)')
    plt.ylabel('$F_{Hertz} \; [nN]$',weight='bold',size='x-large')
    plt.xlabel('$[t/T]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(tt[i],hertzt[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(max(tt[0])-3,max(tt[0]))
    plt.legend()
    plt.savefig('hertzforce_caso'+str(xcasos)+'.png')
    plt.show()
    
    #def ftvisco():
    #=========viscoelastic force vs. time==================
    plt.figure(18)
    plt.grid(True)
    plt.title('$F_{viscosa}$($d$)')
    plt.ylabel('$F_{viscosa} \; [nN]$',weight='bold',size='x-large')
    plt.xlabel('$[t/T]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(tt[i],viscot[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1] )
    plt.xlim(max(tt[0])-3,max(tt[0]))
    plt.legend()
    plt.savefig('viscoforce_caso'+str(xcasos)+'.png')
    plt.show()
    
    
    #def ftcap():
    #=========capilar force vs. time============================
    plt.figure(19)
    plt.grid(True)
    plt.title('$F_{capilar}$($d$)')
    plt.ylabel('$F_{cap} \; [nN]$',weight='bold',size='x-large')
    plt.xlabel('$[t/T]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(tt[i],capt[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(max(tt[0])-3,max(tt[0]))
    plt.legend()
    plt.savefig('capforce_caso'+str(xcasos)+'.png')
    plt.show()
    
    
    #def ftlj():
    #=========lennard jones force vs. time============================
    plt.figure(20)
    plt.grid(True)
    plt.title('$F_{Lennard-Jones}$($t/T$)')
    plt.ylabel('$F_{LJ} \; [nN]$',weight='bold',size='x-large')
    plt.xlabel('$[t/T]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(tt[i],ljt[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(max(tt[0])-3,max(tt[0]))
    plt.legend()
    plt.savefig('ljforce_caso'+str(xcasos)+'.png')
    plt.show()
    
    
    #def ftdlvo():
    #=========DLVO force vs. time============================
    plt.figure(21)
    plt.grid(True)
    plt.title('$F_{DLVO}$($t/T$)')
    plt.ylabel('$F_{DLVO} \; [nN]$',weight='bold',size='x-large')
    plt.xlabel('$[t/T]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(tt[i],dlvot[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1])
    plt.xlim(max(tt[0])-3,max(tt[0]))
    plt.legend()
    plt.savefig('dlvoforce_caso'+str(xcasos)+'.png')
    plt.show()
    
    
    #def p1a1():
    #========Phase vs. Amp =================
    plt.figure(22)
    plt.grid(True)
    plt.title("$\phi$ ($A$)")
    plt.ylabel('$\phi \; [deg]$',weight='bold',size='x-large')
    plt.xlabel('$A\; [nm]$',weight='bold',size='x-large')
    for i in range(len(xcasos)):
        plt.plot(amp1d[i],fase1d[i], label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1] )
    plt.legend()
    plt.savefig('phi1Vamp1_caso'+str(xcasos)+'.png')
    plt.show()
    
    #def asp():
    #========Fmax vs. Asp=================
    Aspn = [0.65,0.75,0.85,0.95]
    for i in range(len(xcasos)):
        plt.figure(23)
        plt.grid(True)
        plt.title("$F_{max}$($Asp$)")
        plt.ylabel('$F_{max} [nN]$',weight='bold',size='x-large')
        plt.xlabel('$Asp [nm]$',weight='bold',size='x-large')
        fmaxdasp=[]
        for j in range (len(Asp)):
            aux1 = np.where(amp1d[i] == amp1d[i].flat[np.abs(amp1d[i] - Asp[j]).argmin()])
            fmaxdasp.append(fmaxd[i][aux1])
        plt.plot(Aspn,fmaxdasp,'o--',label='Case '+str(xcasos[i]),color=colores[(xcasos[i])-1] )
    plt.legend()
    plt.savefig('FmaxVAsp_caso'+str(xcasos)+'.png')
    plt.show() 
    
    
    os.chdir(direc)

if __name__ == '__main__':
    plotter (casos_a,escenario)