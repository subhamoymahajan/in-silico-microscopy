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

lam=np.zeros(10,dtype=int)
lam_col_hue=np.zeros(10)
lam_I0=np.zeros(10)

import sys
def get_img(filename,scale,size,I0,lam,T,ti,fs):
    if ti<0 and T>1:
        print('Error: T>1 but ti<0')
        sys.exit(-1)
    IMG=[]
    Timg=[]
    for i in range(T):
        if ti>=0:
            fname=filename+str(i+ti)+'_lam'+str(lam)+'_fs'+str(fs)+'.dat'
        else:
            fname=filename+'_lam'+str(lam)+'_fs'+str(fs)+'.dat'
        
        f=open(fname,'r')
        if i==0:
            for lines in f:
                foo=lines.split()
                if foo[0][0]=='#':
                    continue
                img_row=[]
                t_row=[]
                for j in range(len(foo)):
                    I=float(foo[j])*I0
                    if I>1:#max intensity value is 1
                        img_row.append(1.0)
                        t_row.append(1)
                    elif I<0: #white background image
                        img_row.append(0.0)
                        t_row.append(0)
                    else:
                        img_row.append(I)
                        t_row.append(1)
                IMG.append(img_row)
                Timg.append(t_row)
        if i>0:
            j=0
            for lines in f:
                foo=lines.split()
                if foo[0][0]=='#':
                    continue
                for k in range(len(foo)):
                    I=float(foo[k])*I0
                    if I>1:#max intensity value is 1
                        IMG[j][k]=IMG[j][k]+1
                        Timg[j][k]=Timg[j][k]+1
                    elif I>0: #microscopy image
                        IMG[j][j]=IMG[j][k]+I
                        Timg[j][k]=Timg[j][k]+1
                j+=1
    for j in range(len(IMG)):
        for k in range(len(IMG[0])):
            if Timg[j][k]>0:
                IMG[j][k]=IMG[j][k]/Timg[j][k]
            else:
                IMG[j][k]=1 #white background image

    L=int(scale/size*len(IMG[0]))
    Llast=int(0.9*len(IMG[0]))
    Hlast=int(0.9*len(IMG))
    wid=int(len(IMG)*0.005)
    for i in range(Llast-L,Llast+1):
        for j in range(Hlast-wid,Hlast+wid+1):
            IMG[j][i]=1.0
    return IMG

# read parameters
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

import matplotlib.pyplot as plt
for i in range(len(lam)):
    if lam[i]>0.5:
        IMG=get_img(filename,scale,size,lam_I0[i],lam[i],T,ti,fs)
        img_width=len(IMG[0])*3/float(len(IMG))
        fig,ax=plt.subplots(1,1,figsize=(img_width,3))
        ax.imshow(IMG,vmin=0, vmax=1,cmap='gray')
        fontS=12
        if img_width<3:
            fontS=img_width*72/18.0
        L=int(scale/size*len(IMG[0]))
        Llast=int(0.9*len(IMG[0]))
        #ax.text(Llast-L,0.97*len(IMG),str(scale)+' nm',fontsize=fontS,transform=ax.transData,color='w') #Uncomment to display scalebar text on image
        ax.set_xticks([])
        ax.set_yticks([])
        plt.axis('off')
        plt.tight_layout(pad=0)
        #plt.show()
        if ti>=0:
            plt.savefig('mono_'+filename+str(ti)+'_lam'+str(lam[i])+'_fs'+str(fs)+'_T'+str(T)+'_I'+str(lam_I0[i])+'.png',dpi=600)
        else:
            plt.savefig('mono_'+filename+'_lam'+str(lam[i])+'_fs'+str(fs)+'_T'+str(T)+'_I'+str(lam_I0[i])+'.png',dpi=600)


