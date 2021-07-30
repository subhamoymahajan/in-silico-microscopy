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

def gen_vid(begin_name,end_name,vid_ext,fps,tbegin,tmax,tdiff,fourcc):
    """ Generates a video from given begining and ending names of JPEG or PNG
        images.
    
    Parameters
    ----------
    begin_name: str
        Begining name of the video and image files
    end_name: str
        Ending name of images and video file.
    vid_ext: str
        Format of output video.
    fps: int
        Frames per second for the output video.
    tbegin: int
        index of first timestep, which is included.
    tmax: int
        index of maximum timestep, tmax is not included.
    tdiff: int 
        difference between timesteps to generate the video.
    fourcc:
        Codec for the output video
   
    Writes
    ------
    A video file from given images.
    """
    vidname=begin_name+end_name+vid_ext
    fname=begin_name+str(tbegin)+end_name
    img_ext=''
    if os.path.exists(fname+'.jpeg'):
        img_ext='.jpeg'
    elif os.path.exists(fname+'.png'):
        img_ext='.png'
    else:
        raise Exception("File "+fname+" (.png or jpeg) does not exist")
    img0=cv2.imread(fname+img_ext)
    h,w,l=img0.shape
    fourcc=cv2.VideoWriter_fourcc(*fourcc)
    video=cv2.VideoWriter(vidname,fourcc,fps,(w,h))
    for i in range(tbegin,tmax,tdiff):
        fname=begin_name+str(i)+end_name+img_ext
        if not os.path.exists(fname):
            raise Exception("File "+fname+" does not exist")
        video.write(cv2.imread(fname))
    video.release()
    print("Writing: "+vidname)

def gen_vid_mono(filename,tbegin,tmax,tdiff,T,fs,lam_I0s,lams,vid_ext,fps,
    fourcc):
    """ Generates monochrome videos. see ''gen_vid'' for more details. filename
        is same as begin_name in ''gen_vid''
    
    Parameters
    ----------
    T: int
        Number of timesteps to perform an average.
    fs: int
        Scaling factor for wave vector or MS position coordinates.
    lam_I0s: array of floats
        The maximum image intensities of all fluorophore types
    lams: array of int
        The wavelength of all fluorophore types

    Writes
    ------
    A video file from given images.
    """
    for i in range(len(lams)):
        end_name='_lam'+str(lams[i])+'_fs'+str(fs)+'_T'+str(T)+'_I'+str(lam_I0s[i])
        gen_vid(filename,end_name,vid_ext,fps,tbegin,tmax,tdiff,fourcc)
    
def gen_vid_col(filename,tbegin,tmax,tdiff,T,fs,lam_I0s,vid_ext,fps,fourcc):
    """ Generates color videos. see ''gen_vid_mono'' for more details.    
 
    """
    Istring=''
    for i in range(len(lam_I0s)):
        Istring+='_'+str(lam_I0s[i])
    
    end_name='_fs'+str(fs)+'_T'+str(T)+'_I'+Istring
    gen_vid(filename,end_name,vid_ext,fps,tbegin,tmax,tdiff,fourcc)


def gen_vid_data(datafile,outname,fps,fourcc):
    """ Generate videos from images specified in a file

    Parameters
    ----------
    datafile: str
        filename of the file containing all image filenames.
    outname: str
        output name of the video file
    fps:
        see ''gen_vid''
    fourcc:
        see ''gen_vid''
 
    Writes
    ------
    A video file from a list of images.
    """
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
