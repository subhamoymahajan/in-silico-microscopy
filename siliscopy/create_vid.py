#    Copyright 2020,2021 SUBHAMOY MAHAJAN 
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
import cv2

def gen_vid(begin_name,end_name,vid_ext,fps,t0,tmax,tdiff,fourcc):
    """ Generates a video from given begining and ending names of JPEG images.
    
    Parameters
    ----------
    """
    vidname=begin_name+end_name+vid_ext
    fname=begin_name+str(t0)+end_name+'.jpeg'
    if not os.path.exists(fname):
        raise Exception("File "+fname+" does not exist")
    img0=cv2.imread(fname)
    h,w,l=img0.shape
    fourcc=cv2.VideoWriter_fourcc(*fourcc)
    video=cv2.VideoWriter(vidname,fourcc,fps,(w,h))
    for i in range(t0,tmax,tdiff):
        fname=begin_name+str(i)+end_name+'.jpeg'
        if not os.path.exists(fname):
            raise Exception("File "+fname+" does not exist")
        video.write(cv2.imread(fname))
    video.release()
    print("Writing: "+vidname)

def gen_vid_mono(filename,t0,tmax,tdiff,T,fs,lam_I0s,lams,vid_ext,fps,fourcc):
    for i in range(len(lams)):
        end_name='_lam'+str(lams[i])+'_fs'+str(fs)+'_T'+str(T)+'_I'+str(lam_I0s[i])
        gen_vid(filename,end_name,vid_ext,fps,t0,tmax,tdiff,fourcc)
    
def gen_vid_col(filename,t0,tmax,tdiff,T,fs,lam_I0s,vid_ext,fps,fourcc):
    Istring=''
    for i in range(len(lam_I0s)):
        Istring+='_'+str(lam_I0s[i])
    
    end_name='_fs'+str(fs)+'_T'+str(T)+'_I'+Istring
    gen_vid(filename,end_name,vid_ext,fps,t0,tmax,tdiff,fourcc)


def gen_vid_data(datafile,outname,fps,fourcc):
    f=open(datafile,'r')
    i=0
    fourcc=cv2.VideoWriter_fourcc(*fourcc)

    for lines in f:
        fname=lines.strip()
        if not os.path.exists(fname):
            raise Exception("File "+fname+" does not exist")
        img=cv2.imread(fname)         
        if i==0:
            h,w,l=img.shape
            video=cv2.VideoWriter(outname,fourcc,fps,(w,h))
        video.write(img)
        i+=1
    video.release
    print("Writing: "+outname)
