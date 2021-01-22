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

# Creates a colored in-silico microscopy image from in-silico monochrome image 

import numpy as np
import sys
import os
import colorsys

lam=np.zeros(10,dtype=int)
lam_col_hue=np.zeros(10)
lam_I0=np.zeros(10)

for i in range(2,len(sys.argv),2):
    if sys.argv[i-1]=='-f':
        filename=sys.argv[i] #Starting filename without lambda or fs values
    elif sys.argv[i-1]=='-p':
        paramfile=sys.argv[i]
    elif sys.argv[i-1]=='-t':
        ti=int(sys.argv[i]) 

f=open(paramfile)
for lines in f:
    foo=lines.split('=')
    if len(foo)<2:
        continue
    foo1=foo[1].split()
    if foo[0]=='T':
        T=int(foo1[0])
    elif foo[0]=='fs':
        fs=int(foo1[0])
    elif foo[0]=='I0':
        I0=float(foo1[0])
    elif foo[0]=='scale':
        scale=float(foo1[0])
    elif foo[0]=='size':
        size=float(foo1[0])
    elif foo[0][0:3]=='lam':
        if foo[0][-4:]=='_hue':
            lam_col_hue[int(foo[0][3])-1]=float(foo1[0])/360.0
        elif foo[0][-3:]=='_I0':
            lam_I0[int(foo[0][3])-1]=float(foo1[0])    
        else:
            lam[int(foo[0][3])-1]=int(foo1[0])

print('T = '+str(T)+'\nfs = '+str(fs)+'\nt = '+str(ti)+'\n')

for i in range(10):
    if lam[i]<1: # no lambda
        break
nlam=i 
print('nlam = '+str(nlam))
print('lam_I0 = ', lam_I0[:nlam])    

Img_dat=[]
T_dat=[]
for n in range(T):
    for i in range(nlam):
        if ti<0:#Just a single file with no time dimension
            fname=filename+'_lam'+str(lam[i])+'_fs'+str(fs)+'.dat'
        else:#Multiple files for different time whose timestep is mentioned right after 'filename'
            fname=filename+str(ti+n)+'_lam'+str(lam[i])+'_fs'+str(fs)+'.dat'

        print('Opening file: '+fname)       
        f=open(fname,'r')
        if n==0:
            I=[]
            Tini=[]
            for lines in f:
                foo=lines.split()
                if foo[0][0]=='#':
                    continue
                Irow=[]
                Trow=[]
                for j in range(len(foo)):
                    Irow.append(float(foo[j]))
                    if float(foo[j])>0:
                        Trow.append(1)
                    else:
                        Trow.append(0)
                I.append(Irow)
                Tini.append(Trow)
            Img_dat.append(I)
            T_dat.append(Tini)
        else:
            j=0
#Img_dat < 0 implies that the box size is smaller than expected. Should be rendered white (background color)             
            for lines in f:
                foo=lines.split()
                if foo[0][0]=='#':
                    continue
                for k in range(len(foo)):
                    if float(foo[k])>0:
                        if Img_dat[i][j][k]>0:
                            Img_dat[i][j][k]+=float(foo[k])
                        else:
                            Img_dat[i][j][k]=float(foo[k])
                        T_dat[i][j][k]+=1
                j+=1    

#Average Intensity for each wavelength
for i in range(len(Img_dat)):
    for j in range(len(Img_dat[0])):
        for k in range(len(Img_dat[0][0])):
            if Img_dat[i][j][k]>0 and T_dat[i][j][k]>1:
                Img_dat[i][j][k]=Img_dat[i][j][k]/float(T_dat[i][j][k])
          

Img_dat1=np.array(Img_dat)
IMG=np.zeros((len(Img_dat[0]),len(Img_dat[0][0]),3))
for i in range(len(Img_dat[0])):
    for j in range(len(Img_dat[0][0])):
        foo=0
        for lam_id in range(nlam):
            if Img_dat[lam_id][i][j]>-1E-10:
                foo=+1
        if foo==0:
            #if all Img_dat are -1 then res is -1
            #Making the color blue now, later make it white
            IMG[i][j][0]=1      
            IMG[i][j][1]=1
            IMG[i][j][2]=1
        else:
            xres=0
            yres=0
            Is=[]
            ncol=0
            for lam_id in range(nlam):
                if Img_dat[lam_id][i][j]>1E-10:
                    xres+=Img_dat[lam_id][i][j]*np.cos(lam_col_hue[lam_id]*2*np.pi)*lam_I0[lam_id]    
                    yres+=Img_dat[lam_id][i][j]*np.sin(lam_col_hue[lam_id]*2*np.pi)*lam_I0[lam_id]    
                    Is.append(Img_dat[lam_id][i][j]*lam_I0[lam_id])
                    ncol+=1
            if ncol==0:
                continue
            Is=sorted(Is) 
            hres=np.arctan2(yres,xres)/(2*np.pi)
            if hres<0:
                hres+=1
            vres=Is[-1]
            if vres>1:
                vres=1
            if ncol>2:
                sres=1-Is[-3]/Is[-1]
            else:
                sres=1

            rgb=list(colorsys.hsv_to_rgb(hres,sres,vres))
            IMG[i][j][0]=rgb[0]
            IMG[i][j][1]=rgb[1]
            IMG[i][j][2]=rgb[2]

L=int(scale/size*len(IMG[0]))
Llast=int(0.9*len(IMG[0]))
Hlast=int(0.9*len(IMG))
wid=int(len(IMG)*0.005)
for i in range(Llast-L,Llast+1):
    for j in range(Hlast-wid,Hlast+wid+1):
        IMG[j][i][0]=1.0
        IMG[j][i][1]=1.0
        IMG[j][i][2]=1.0

img_width=len(IMG[0])*3/float(len(IMG))
fontS=12
if img_width<3:
    fontS=img_width*72/18.0

import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,1,figsize=(img_width,3))
ax.imshow(IMG)
ax.set_xticks([])
#ax.text(Llast-L,0.97*len(IMG),str(scale)+' nm',fontsize=fontS,transform=ax.transData,color='w') #Uncomment to display scalebar text on image
ax.set_yticks([])
plt.axis('off')
plt.tight_layout(pad=0)
if ti<0:
    cnt=0
    while os.path.exists(filename+str(cnt)+'.png'):
        cnt+=1
#    plt.show()
    print('Output: '+filename+str(cnt)+'.png')
    
    plt.savefig(filename+str(cnt)+'.png',dpi=600)
else:
#    plt.show()
    Istring=''
    for i in range(10):
        if lam_I0[i]>1E-10:
            Istring=Istring+'_'+str(round(lam_I0[i],3))

    print('Output: '+filename+str(ti)+'_fs'+str(fs)+'_T'+str(T)+'_I'+str(Istring)+'.png')
    plt.savefig(filename+str(ti)+'_fs'+str(fs)+'_T'+str(T)+'_I'+str(Istring)+'.png',dpi=600)
plt.close()



