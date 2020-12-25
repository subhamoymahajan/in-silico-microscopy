#    Copyright 2020 SUBHAMOY MAHAJAN 
#    
#    This file is part of InSilicoMicroscopy software.
# 
#    InSilicoMicroscopy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.)

import sys
import os
import numpy as np
lam=np.zeros(10,dtype=int)
lam_col_hue=np.zeros(10)
lam_I0=np.zeros(10)

for i in range(2,len(sys.argv),2):
    if sys.argv[i-1]=='-f':
        filename=sys.argv[i] #Starting filename without lambda or fs values
    elif sys.argv[i-1]=='-p':
        paramfile=sys.argv[i]
    elif sys.argv[i-1]=='-tmax':
        tmax=int(sys.argv[i]) 
    elif sys.argv[i-1]=='-t0':
        t0=int(sys.argv[i]) 
    elif sys.argv[i-1]=='-tdiff':
        tdiff=int(sys.argv[i]) 

f=open(paramfile)
for lines in f:
    foo=lines.split('=')
    foo1=foo[1].split()
    if foo[0]=='T':
        T=int(foo1[0])
    elif foo[0]=='fs':
        fs=int(foo1[0])
    elif foo[0][0:3]=='lam':
        if foo[0][-4:]=='_hue':
            lam_col_hue[int(foo[0][3])-1]=float(foo1[0])/360.0
        elif foo[0][-3:]=='_I0':
            print(foo1)
            lam_I0[int(foo[0][3])-1]=float(foo1[0])    
        else:
            lam[int(foo[0][3])-1]=int(foo1[0])

Istring=''
for i in range(10):
    if lam_I0[i]>1E-10:
        Istring=Istring+'_'+str(round(lam_I0[i],3))
print('T = '+str(T)+'\nfs = '+str(fs)+'\ntmax = '+str(tmax)+'\nI0 = '+str(lam_I0)+'\ntdiff = '+str(tdiff))    
    
import cv2
vidname=filename+'_fs'+str(fs)+'_T'+str(T)+'_I'+str(Istring)+'.avi'
filename=filename+str(t0)'_fs'+str(fs)+'_T'+str(T)+'_I'+str(Istring)+'.png'
if os.path.exists(filename):
    img0=cv2.imread(filename)
else:
    print(filename+' does not exist')
    sys.exit()
h,w,l=img0.shape
video=cv2.VideoWriter(vidname,0,5,(w,h))

for i in range(t0,tmax,tdiff):
    filename=filename+str(i)+'_fs'+str(fs)+'_T'+str(T)+'_I'+str(Istring)+'.png'
    if os.path.exists(filename):
        video.write(cv2.imread(filename))
    else:
        print (filename+' does not exists')
        sys.exit()
cv2.destroyAllWindows()
video.release()

