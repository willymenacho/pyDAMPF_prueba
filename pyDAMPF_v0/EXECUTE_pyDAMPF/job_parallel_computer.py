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

##################################################################################
###################   STRONG SCALING SUBMISSION SCRIPT     
#######################
import os
from sys import argv 
script, factor = argv
factor = int(factor)

#factor = 10
def compute (factor):
	direc = os.getcwd()
	for i in range(1,factor+1):  #change
	    os.chdir(direc+'/PARALLELBASIC_0/'+str(i)+'/nrun')
	    os.system('python3 generate_cases.py &')
	    os.chdir(direc)
	    
	print('pyDAMPF execution is finished ')
if __name__ == '__main__':
    compute (factor)   
