import unittest
import warnings
import siliscopy
import numpy as np
import cffi
import os
import importlib
import ctypes
import copy
import cv2
import tifffile as tif

class TestPSFGeneration(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter("ignore")
        #beta=1 rad
        self.NA=1.
        self.meu=1./np.sin(1.)
        self.dlmn=[1,1,1]
        self.Plmn=[1,1,1]
        self.fs=1
        #meu*k'*(r'=1)=5
        self.lam=2*np.pi*self.meu/5.
        #Gibson Lanni
        self.lam_gl=2*np.pi/5.
        self.meu0=self.meu
        self.t0=1
        self.meug=1.522
        self.meug0=1.522
        self.tg=1
        self.tg0=1
        self.meus=1.33
        self.tsO=0
    def tearDown(self):
        if os.path.exists('test.dat'):
            os.remove('test.dat')

    def test_Gandy_r(self):
        PSF0=[]
        PSF0.append(1.)
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append((0.0381952*1.5/(1-np.cos(1.)**1.5))**2)
        PSF0.append((0.0396359*1.5/(1-np.cos(1.)**1.5))**2)
        siliscopy.gen_psf.psf_gandy_sep(self.NA, self.meu, self.lam, self.dlmn,
            self.Plmn, self.fs, 'test.dat','w',0)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+ \
                str(PSF0[i]))
            i+=1 
        f.close()

    def test_Gandy_n(self): #n' = 1
        PSF0=[]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append((1.5/(1-np.cos(1.)**1.5))**2*(0.229135**2+0.223351**2))
        PSF0.append((1.5/(1-np.cos(1.)**1.5))**2*(0.0829956**2+0.0268196**2))
        PSF0.append((1.5/(1-np.cos(1.)**1.5))**2*(0.0402964**2+0.0287536**2))
        siliscopy.gen_psf.psf_gandy_sep(self.NA, self.meu, self.lam, self.dlmn,
            self.Plmn, self.fs, 'test.dat','w',1)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+ \
                str(PSF0[i]))
            i+=1 
        f.close()
        
    def test_GL1991_r(self): #OPD = 0
        PSF0=[1.]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0171693)  
        PSF0.append(0.0000219836)
        siliscopy.gen_psf.psf_GL1991_sep(self.NA, self.meu, self.meu0, self.t0,
            self.tsO, self.meus, self.tg, self.tg0, self.meug, self.meug0, 
            self.lam_gl, self.dlmn, self.Plmn, self.fs, 'test.dat','w',0)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+ \
                str(PSF0[i]))
            i+=1
        f.close()

    def test_GL1991_n(self): #OPD only due to n'
        PSF0=[0.643715]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0382171)  
        PSF0.append(0.00481841)
        siliscopy.gen_psf.psf_GL1991_sep(self.NA, self.meu, self.meu0, self.t0,
            self.tsO, self.meus, self.tg, self.tg0, self.meug, self.meug0, 
            self.lam_gl, self.dlmn, self.Plmn, self.fs, 'test.dat','w',1)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+ \
                str(PSF0[i]))
            i+=1
        f.close()

    def test_GL1991_tsO(self): #OPD only due to tsO
        PSF0=[0.997761]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0171865)  
        PSF0.append(0.0000705946)
        siliscopy.gen_psf.psf_GL1991_sep(self.NA, self.meu, self.meu0, self.t0,
            1, self.meus, self.tg, self.tg0, self.meug, self.meug0, 
            self.lam_gl, self.dlmn, self.Plmn, self.fs, 'test.dat','w',0)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+ \
                str(PSF0[i]))
            i+=1
        f.close()

    def test_GL1991_meug(self): #OPD only due to tsO
        PSF0=[0.999982]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0171697)  
        PSF0.append(0.0000223409)
        siliscopy.gen_psf.psf_GL1991_sep(self.NA, self.meu, self.meu0, self.t0,
            self.tsO, self.meus, self.tg, self.tg0, 1.6, self.meug0, 
            self.lam_gl, self.dlmn, self.Plmn, self.fs, 'test.dat','w',0)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+ \
                str(PSF0[i]))
            i+=1
        f.close()

    def test_GL1991_tg(self): #OPD only due to tsO
        PSF0=[0.997761]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0171865)  
        PSF0.append(0.0000705946)
        siliscopy.gen_psf.psf_GL1991_sep(self.NA, self.meu, self.meu0, self.t0,
            self.tsO, self.meus, self.tg0+1, self.tg0, self.meus, self.meus, 
            self.lam_gl, self.dlmn, self.Plmn, self.fs, 'test.dat','w',0)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+ \
                str(PSF0[i]))
            i+=1
        f.close()

    def test_GL1991_toil(self): #OPD only due to tsO
        PSF0=[0.997761]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0171865)  
        PSF0.append(0.0000705946)
        siliscopy.gen_psf.psf_GL1991_sep(self.NA, self.meu, self.meus, -1,
            self.tsO, self.meus, self.tg0, self.tg0, self.meug, self.meug0, 
            self.lam_gl, self.dlmn, self.Plmn, self.fs, 'test.dat','w',0)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+ \
                str(PSF0[i]))
            i+=1
        f.close()

##############################################
    def test_Mod_Gandy_r(self): #OPD = 0
        PSF0=[1.]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0222045)  
        PSF0.append(0.000233725)
        siliscopy.gen_psf.psf_Mod_Gandy_sep(self.NA, self.meu, self.meu0, 
            self.t0, self.tsO, self.meus, self.tg, self.tg0, self.meug, 
            self.meug0, self.lam_gl, self.dlmn, self.Plmn, self.fs, 
            'test.dat','w',0)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+ \
                str(PSF0[i]))
            i+=1
        f.close()

    def test_Mod_Gandy_n(self): #OPD only due to n'
        PSF0=[0.638634]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0410877)  
        PSF0.append(0.00547789)
        siliscopy.gen_psf.psf_Mod_Gandy_sep(self.NA, self.meu, self.meu0, 
            self.t0, self.tsO, self.meus, self.tg, self.tg0, self.meug, 
            self.meug0, self.lam_gl, self.dlmn, self.Plmn, self.fs,
            'test.dat','w',1)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+ \
                str(PSF0[i]))
            i+=1
        f.close()

    def test_Mod_Gandy_tsO(self): #OPD only due to tsO
        PSF0=[0.997601]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0222105)  
        PSF0.append(0.000287965)
        siliscopy.gen_psf.psf_Mod_Gandy_sep(self.NA, self.meu, self.meu0, 
            self.t0, 1, self.meus, self.tg, self.tg0, self.meug, self.meug0, 
            self.lam_gl, self.dlmn, self.Plmn, self.fs, 'test.dat','w',0)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+ \
                str(PSF0[i]))
            i+=1
        f.close()

    def test_Mod_Gandy_meug(self): #OPD only due to tsO
        PSF0=[0.999982]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0222048)  
        PSF0.append(0.000234126)
        siliscopy.gen_psf.psf_Mod_Gandy_sep(self.NA, self.meu, self.meu0, 
            self.t0, self.tsO, self.meus, self.tg, self.tg0, 1.6, self.meug0, 
            self.lam_gl, self.dlmn, self.Plmn, self.fs, 'test.dat','w',0)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+ \
                str(PSF0[i]))
            i+=1
        f.close()

    def test_Mod_Gandy_tg(self): #OPD only due to tg
        PSF0=[0.997601]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0222105)  
        PSF0.append(0.000287965)
        siliscopy.gen_psf.psf_Mod_Gandy_sep(self.NA, self.meu, self.meu0, 
            self.t0, self.tsO, self.meus, self.tg0+1, self.tg0, self.meus, 
            self.meus, self.lam_gl, self.dlmn, self.Plmn, self.fs, 
            'test.dat','w',0)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+ \
                str(PSF0[i]))
            i+=1
        f.close()

    def test_Mod_Gandy_toil(self): #OPD only due to toil
        PSF0=[0.997601]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0222105)  
        PSF0.append(0.000287965)
        siliscopy.gen_psf.psf_Mod_Gandy_sep(self.NA, self.meu, self.meus, -1,
            self.tsO, self.meus, self.tg, self.tg0, self.meug, self.meug0, 
            self.lam_gl, self.dlmn, self.Plmn, self.fs, 'test.dat','w',0)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+ \
                str(PSF0[i]))
            i+=1
        f.close()

def readIMG(filename):
    f=open(filename,'r')
    data=[]
    for lines in f:
        if lines[0]=='#':
            continue
        foo=lines.split()
        foo=[float(x) for x in foo ] 
        data.append(foo)
    data=np.array(data)
    f.close()
    return data

def writePSF(filename,PSF):
    f=open(filename,'w')
    f.write('#Comment message\n')
    for i in range(len(PSF)):
        f.write(str(PSF[i,0]) +' '+str(PSF[i,1]) +' '+str(PSF[i,2]) +' '+ \
                str(PSF[i,3]) +'\n')
    f.close()

def write_param(filename,params):
    f=open(filename,'w')
    for key in params:
        if type(params[key]) in [int, float, str]:
            f.write(key+'='+str(params[key]))        
        elif type(params[key])==list:
            f.write(key+'=')
            for i in range(len(params[key])):
                f.write(str(params[key][i])+' ')
        f.write('\n')
    f.close()

def write_pos(filename,pos,atoms,box):
    f=open(filename,'w')
    f.write('Comment : Molecule t= 0.12\n '+str(len(pos))+'\n')
    for i in range(len(pos)):
        f.write('    1'+'DC   '+'   '+str(atoms[i])+' '+str(i).ljust(5)+ \
                str(int(pos[i,0]*1000)/1000).ljust(8)+ \
                str(int(pos[i,1]*1000)/1000).ljust(8)+ \
                str(int(pos[i,2]*1000)/1000).ljust(8)+'\n')
    f.write(' '+str(box[0])+' '+str(box[1])+' '+str(box[2]))
    f.close()

class TestMonoGeneration(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        pass

    @classmethod
    def tearDownClass(self):
        pass

    def setUp(self):
        self.params={'frame_col': 1.0, 'mix_type': 'mt', 'psf_type': 0, 'T': 1,
                     'add_n': 1, 'fourcc': "'mp4v'", 'vid_ext': '.mov',  
                     'fps': 1, 'tdiff': 1000, 'tmax': 21000, 'tbegin': 0, 
                     'dpi': 600, 'scale': 5, 'lam_I0_1': 0.13, 'lam_I0_2': 0.27,
                     'lam_hue1': 255, 'lam_hue2': 60, 'tsO': 0, 'meus': 1.33,
                     'meu0': 1.51, 'meu': 1.51, 'NA': 1.3, 'pbc': 'xyz', 
                     'Plmn': [1.0, 1.0, 1.0], 'dlmn': [0.2, 0.2, 0.2], 'fs': 1, 
                     'lam1': 1, 'lam2': 2, 'lam_names1': 'A', 'lam_names2': 'B',
                     'maxlen': [1.2, 1.2, 1.2], 'focus_cor': 0, 'opt_axis': 2}

    def tearDown(self):
        os.system('rm -f out*.dat psf*.dat param.dat pos.gro')

    def test_img1(self): #In-focus image for psf_type=0
        #Also Checks 'white image frame', pbc = in xy with z as optical axis

        #Write PSF
        PSF1=np.array([[0,0,0,1.],[0.2,0,0,0.9],[0.2,0.2,0,0.8]])
        PSF2=np.array([[0,0,0,0.3],[0.2,0,0,0.2],[0.2,0.2,0,0.1]])
       
        writePSF('psf_lam1_fs1.dat',PSF1)
        writePSF('psf_lam2_fs1.dat',PSF2)

        #Write Parameter file
        write_param('param.dat',self.params)

        #Write positions file
        atoms,pos,box=[np.array(['A']),np.array([[0.,0.,0.]]),[1.,1.,1.]]
        write_pos('pos.gro',pos,atoms,box)
     
        #Generate image data
        siliscopy.gen_mono.gen_mono_c(['pos.gro','param.dat','psf','out'], 
            silent=True)

        IMG1=readIMG('out_lam1_fs1.dat')
        IMG2=readIMG('out_lam2_fs1.dat')
        self.assertEqual(IMG1.shape,(6,6))
        self.assertEqual(IMG2.shape,(6,6))
        
        IMG1_0=np.zeros((6,6))
        IMG1_0[5,:]=-1.
        IMG1_0[:,5]=-1.
        IMG1_0[:5,0]=[1,0.9,0.,0.,0.9]
        IMG1_0[:5,1]=[0.9,0.8,0.,0.,0.8]
        IMG1_0[:5,4]=[0.9,0.8,0.,0.,0.8]

        IMG2_0=np.zeros((6,6))
        IMG2_0[5,:]=-1.
        IMG2_0[:,5]=-1.
        self.assertTrue((abs(IMG1-IMG1_0)<1E-6).all())
        self.assertTrue((abs(IMG2-IMG2_0)<1E-6).all())
 
    def test_img2(self): 
        # Checks if correct PSF is read when dlmn in parameter file is coarser 
        # than the dlmn used to generate PSF. 
        pass
        #Write PSF
        PSF1=np.array([[0,0,0,1.],[0.2,0,0,0.9],[0.2,0.2,0,0.8],
             [0,0,0.1,0.6],[0.2,0,0.1,0.5],[0.2,0.2,0.1,0.4],
             [0,0,0.2,0.3],[0.2,0,0.2,0.2],[0.2,0.2,0.2,0.1]])

        PSF2=np.array([[0,0,0,0.9],[0.2,0,0,0.8],[0.2,0.2,0,0.7],
             [0,0,0.1,0.5],[0.2,0,0.1,0.4],[0.2,0.2,0.1,0.3],
             [0,0,0.2,0.2],[0.2,0,0.2,0.1],[0.2,0.2,0.2,0.05]])
       
        writePSF('psf_lam1_fs1.dat',PSF1)
        writePSF('psf_lam2_fs1.dat',PSF2)

        write_param('param.dat',self.params)

        #Write positions file
        atoms,pos,box=[np.array(['A','B']),np.array([[0.15,0.15,0.06], 
            [0.09,0.71,0.11]]),[1.0,1.0,1.0]]
        write_pos('pos.gro',pos,atoms,box)
     
        #Generate image data
        siliscopy.gen_mono.gen_mono_c(['pos.gro','param.dat','psf','out'],
            silent=True)

        IMG1=readIMG('out_lam1_fs1.dat')
        IMG2=readIMG('out_lam2_fs1.dat')
        self.assertEqual(IMG1.shape,(6,6))
        self.assertEqual(IMG2.shape,(6,6))
        
        IMG1_0=np.zeros((6,6))
        IMG1_0[5,:]=-1.
        IMG1_0[:,5]=-1.
        IMG1_0[0,:3]=[0.8,0.9,0.8]
        IMG1_0[1,:3]=[0.9,1.0,0.9]
        IMG1_0[2,:3]=[0.8,0.9,0.8]

        IMG2_0=np.zeros((6,6))
        IMG2_0[5,:]=-1.
        IMG2_0[:,5]=-1.
        IMG2_0[0,:5]=[0.1,0.,0.,0.1,0.2]
        IMG2_0[1,:5]=[0.05,0.,0.,0.05,0.1]
        IMG2_0[4,:5]=[0.05,0.,0.,0.05,0.1]

        self.assertTrue((abs(IMG1-IMG1_0)<1E-6).all())
        self.assertTrue((abs(IMG2-IMG2_0)<1E-6).all())



    def test_img3(self): #Multiple out-of-focus image for psf_type=0
        #Also checks location of l,m, addition of intensities
        #Write PSF
        PSF1=np.array([[0,0,0,1.],[0.2,0,0,0.9],[0.2,0.2,0,0.8],
                       [0,0,0.2,0.6],[0.2,0,0.2,0.5],[0.2,0.2,0.2,0.4],
                       [0,0,0.4,0.3],[0.2,0,0.4,0.2],[0.2,0.2,0.4,0.1]])
        PSF2=np.array([[0,0,0,0.3],[0.2,0,0,0.2],[0.2,0.2,0,0.1]])
       
        writePSF('psf_lam1_fs1.dat',PSF1)
        writePSF('psf_lam2_fs1.dat',PSF2)

        #Write Parameter file
        write_param('param.dat',self.params)

        #Write positions file
        atoms,pos,box=[np.array(['A','A']),np.array([[0.15,0.15,-0.15],
            [0.49,0.69,0.31]]), [1.0,1.0,1.0]]
        write_pos('pos.gro',pos,atoms,box)
     
        #Generate image data
        siliscopy.gen_mono.gen_mono_c(['pos.gro','param.dat','psf','out'],
            silent=True)

        IMG1=readIMG('out_lam1_fs1.dat')
        IMG2=readIMG('out_lam2_fs1.dat')
        self.assertEqual(IMG1.shape,(6,6))
        self.assertEqual(IMG2.shape,(6,6))
        
        IMG1_0=np.zeros((6,6))
        IMG1_0[5,:]=-1.
        IMG1_0[:,5]=-1.
        IMG1_0[0,:3]=[0.4,0.5,0.4]
        IMG1_0[1,:5]=[0.5,0.6,0.6,0.2,0.1]
        IMG1_0[2,:5]=[0.4,0.5,0.6,0.3,0.2]
        IMG1_0[3,2:5]=[0.1,0.2,0.1]

        IMG2_0=np.zeros((6,6))
        IMG2_0[5,:]=-1.
        IMG2_0[:,5]=-1.
        
        
        self.assertTrue((abs(IMG1-IMG1_0)<1E-6).all())
        self.assertTrue((abs(IMG2-IMG2_0)<1E-6).all())
    
    def test_img4(self): #Multiple in and out-of-focus image for psf_type=1
        #Checks location of l,m center.
        #Write PSF
        PSF1=np.array([[0,0,0,1.],[0.2,0,0,0.9],[0.2,0.2,0,0.8],
                       [0,0,0.2,0.6],[0.2,0,0.2,0.5],[0.2,0.2,0.2,0.4],
                       [0,0,-0.2,0.3],[0.2,0,-0.2,0.2],[0.2,0.2,-0.2,0.1]])

        PSF2=np.array([[0,0,0,0.9],[0.2,0,0,0.8],[0.2,0.2,0,0.7],
                       [0,0,0.2,0.5],[0.2,0,0.2,0.4],[0.2,0.2,0.2,0.3],
                       [0,0,-0.2,0.2],[0.2,0,-0.2,0.1],[0.2,0.2,-0.2,0.05]])
       
        writePSF('psf_tsO'+"%g"%self.params['tsO']+'_lam1_fs1.dat',PSF1)
        writePSF('psf_tsO'+"%g"%self.params['tsO']+'_lam2_fs1.dat',PSF2)

        #Write Parameter file
        self.params['psf_type']=1
        write_param('param.dat',self.params)

        #Write positions file
        atoms,pos,box=[np.array(['A','B']),np.array([[0.15,0.15,0.21],
            [0.49,0.69,-0.15]]), [1.0,1.0,1.0]]
        write_pos('pos.gro',pos,atoms,box)
     
        #Generate image data
        siliscopy.gen_mono.gen_mono_c(['pos.gro','param.dat','psf','out'],
            silent=True)
        IMG1=readIMG('out_tsO'+"%g"%self.params['tsO']+'_lam1_fs1.dat')
        IMG2=readIMG('out_tsO'+"%g"%self.params['tsO']+'_lam2_fs1.dat')
        self.assertEqual(IMG1.shape,(6,6))
        self.assertEqual(IMG2.shape,(6,6))
        
        IMG1_0=np.zeros((6,6))
        IMG1_0[5,:]=-1.
        IMG1_0[:,5]=-1.
        IMG1_0[:3,0]=[0.1,0.2,0.1]
        IMG1_0[:3,1]=[0.2,0.3,0.2]
        IMG1_0[:3,2]=[0.1,0.2,0.1]

        IMG2_0=np.zeros((6,6))
        IMG2_0[5,:]=-1.
        IMG2_0[:,5]=-1.
        IMG2_0[1:4,2]=[0.3,0.4,0.3]
        IMG2_0[1:4,3]=[0.4,0.5,0.4]
        IMG2_0[1:4,4]=[0.3,0.4,0.3]
        self.assertTrue((abs(IMG1-IMG1_0)<1E-6).all())
        self.assertTrue((abs(IMG2-IMG2_0)<1E-6).all())

 
    def test_img5(self): #Check if changing optical axis works.
        #Write PSF
        PSF1=np.array([[0,0,0,1.],[0.2,0,0,0.9],[0.2,0.2,0,0.8],
             [0,0,0.2,0.6],[0.2,0,0.2,0.5],[0.2,0.2,0.2,0.4],
             [0,0,-0.2,0.3],[0.2,0,-0.2,0.2],[0.2,0.2,-0.2,0.1]])

        PSF2=np.array([[0,0,0,0.3],[0.2,0,0,0.2],[0.2,0.2,0,0.1]])
       
        writePSF('psf_tsO'+"%g"%self.params['tsO']+'_lam1_fs1.dat',PSF1)
        writePSF('psf_tsO'+"%g"%self.params['tsO']+'_lam2_fs1.dat',PSF2)

        #Write Parameter file
        self.params['psf_type']=1
        self.params['opt_axis']=1
        write_param('param.dat',self.params)

        #Write positions file
        atoms,pos,box=[np.array(['A','B']),np.array([[0.15,0.15,-0.21],
            [0.09,0.7,0.25]]), [1.0,1.0,1.0]]
        write_pos('pos.gro',pos,atoms,box)
     
        #Generate image data
        siliscopy.gen_mono.gen_mono_c(['pos.gro','param.dat','psf','out'],
            silent=True)

        IMG1=readIMG('out_tsO'+"%g"%self.params['tsO']+'_lam1_fs1.dat')
        IMG2=readIMG('out_tsO'+"%g"%self.params['tsO']+'_lam2_fs1.dat')
        self.assertEqual(IMG1.shape,(6,6))
        self.assertEqual(IMG2.shape,(6,6))
        
        IMG1_0=np.zeros((6,6))
        IMG1_0[5,:]=-1.
        IMG1_0[:,5]=-1.
        IMG1_0[0,:3]=[0.1,0.2,0.1]
        IMG1_0[3,:3]=[0.1,0.2,0.1]
        IMG1_0[4,:3]=[0.2,0.3,0.2]

        IMG2_0=np.zeros((6,6))
        IMG2_0[5,:]=-1.
        IMG2_0[:,5]=-1.
        self.assertTrue((abs(IMG1-IMG1_0)<1E-6).all())
        self.assertTrue((abs(IMG2-IMG2_0)<1E-6).all())

    def test_img6(self): #Check if changing optical axis works with pbc
        #Write PSF
        PSF1=np.array([[0,0,0,1.],[0.2,0,0,0.9],[0.2,0.2,0,0.8],
             [0,0,0.2,0.6],[0.2,0,0.2,0.5],[0.2,0.2,0.2,0.4],
             [0,0,-0.2,0.3],[0.2,0,-0.2,0.2],[0.2,0.2,-0.2,0.1]])

        PSF2=np.array([[0,0,0,0.9],[0.2,0,0,0.8],[0.2,0.2,0,0.7],
             [0,0,0.2,0.5],[0.2,0,0.2,0.4],[0.2,0.2,0.2,0.3],
             [0,0,-0.2,0.2],[0.2,0,-0.2,0.1],[0.2,0.2,-0.2,0.05]])
       
        writePSF('psf_tsO'+"%g"%self.params['tsO']+'_lam1_fs1.dat',PSF1)
        writePSF('psf_tsO'+"%g"%self.params['tsO']+'_lam2_fs1.dat',PSF2)

        #Write Parameter file
        self.params['psf_type']=1
        self.params['opt_axis']=1
        self.params['pbc']='xy'
        write_param('param.dat',self.params)

        #Write positions file
        atoms,pos,box=[np.array(['A','B']),np.array([[0.15,0.15,-0.21],
            [0.09,0.71,0.25]]), [1.0,1.0,1.0]]
        write_pos('pos.gro',pos,atoms,box)
     
        #Generate image data
        siliscopy.gen_mono.gen_mono_c(['pos.gro','param.dat','psf','out'],
            silent=True)

        IMG1=readIMG('out_tsO'+"%g"%self.params['tsO']+'_lam1_fs1.dat')
        IMG2=readIMG('out_tsO'+"%g"%self.params['tsO']+'_lam2_fs1.dat')
        self.assertEqual(IMG1.shape,(6,6))
        self.assertEqual(IMG2.shape,(6,6))
        
        IMG1_0=np.zeros((6,6))
        IMG1_0[5,:]=-1.
        IMG1_0[:,5]=-1.
        IMG1_0[0,:3]=[0.1,0.2,0.1]

        IMG2_0=np.zeros((6,6))
        IMG2_0[5,:]=-1.
        IMG2_0[:,5]=-1.
        IMG2_0[0,:5]=[0.4,0.3,0.,0.,0.3]
        IMG2_0[1,:5]=[0.5,0.4,0.,0.,0.4]
        IMG2_0[2,:5]=[0.4,0.3,0.,0.,0.3]

        self.assertTrue((abs(IMG1-IMG1_0)<1E-6).all())
        self.assertTrue((abs(IMG2-IMG2_0)<1E-6).all())

def writeIMG(filename,I):
    f=open(filename,'w')
    f.write('# Test here\n')
    for i in range(len(I)):
        for j in range(len(I[0])):
            f.write(str(I[i,j])+' ')
        f.write('\n')
    f.close()

class TestPlotImage(unittest.TestCase):
    def setUp(self):
        self.I1=np.ones((6,6))*-1
        self.I1[2,2]=0.3
        self.I2=np.ones((6,6))*-1
        self.I2[2:4,2:4]=[[0.1, 0.2],[0.4,0.5]]
        self.I3=np.ones((6,6))*-1
        self.I3[1:4,1:4]=[[0.15, 0.25, 0.35],[0.45,0.55,0.65],[0.75,0.85,0.95]]

    def tearDown(self):
        os.system('rm -f img*.dat')
        os.system('rm -f *.tiff')
        os.system('rm -f *.jpeg')

    def test_noise(self):
        I1=copy.deepcopy(self.I1)
        I=siliscopy.plot_image.add_noise(I1,0,0)
        self.assertTrue((abs(I-self.I1)<1E-6).all())
        I1=copy.deepcopy(self.I1)
        I=siliscopy.plot_image.add_noise(I1,1,1)
        for i in range(6):
            for j in range(6):
                if i==2 and j==2:
                    self.assertNotEqual(I[2,2],self.I1[2,2])
                else:
                    self.assertEqual(I[i,j],-1)

        I2=copy.deepcopy(self.I2)
        I=siliscopy.plot_image.add_noise(I2,0,0)
        self.assertTrue((abs(I-self.I2)<1E-6).all())
        I2=copy.deepcopy(self.I2)
        I=siliscopy.plot_image.add_noise(I2,1,1)
        for i in range(6):
            for j in range(6):
                if i>=2 and i<4 and j>=2 and j<4:
                    self.assertNotEqual(I[i,j],self.I2[i,j])
                else:
                    self.assertEqual(I[i,j],-1)

        I3=copy.deepcopy(self.I3)
        I=siliscopy.plot_image.add_noise(I3,0,0)
        self.assertTrue((abs(I-self.I3)<1E-6).all())
        I3=copy.deepcopy(self.I3)
        I=siliscopy.plot_image.add_noise(I3,1,1)
        for i in range(6):
            for j in range(6):
                if i>=1 and i<4 and j>=1 and j<4:
                    self.assertNotEqual(I[i,j],self.I3[i,j])
                else:
                    self.assertEqual(I[i,j],-1)

    def test_get_grey1(self): 
        #no time
        writeIMG('img_lam100_fs1.dat',self.I2)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,-1,1,[6,6],
            frame=True)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img_lam100_fs1.dat')   
        #one time frame
        writeIMG('img10_lam100_fs1.dat',self.I2)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,10,1,[6,6],
            frame=True)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img10_lam100_fs1.dat')   
        #three time frames 
        writeIMG('img10_lam100_fs1.dat',self.I1)
        writeIMG('img11_lam100_fs1.dat',self.I2)
        writeIMG('img12_lam100_fs1.dat',self.I3)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,3,10,1,[6,6],
            frame=True)
        Itrue=np.ones((6,6))*-1
        Itrue[1:4,1:4]=[[0.15,0.25,0.35],[0.45,0.3166666666,0.425],
            [0.75,0.625,0.725]]
        self.assertTrue((abs(IMG-Itrue)<1E-6).all())
        for i in range(10,13):
            os.remove('img'+str(i)+'_lam100_fs1.dat')   

        #test 3D slices x          
        writeIMG('img_x20_lam100_fs1.dat',self.I2)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,-1,1,[6,6],
            frame=True,nidx=20,opt_axis=0)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img_x20_lam100_fs1.dat')   
        #one time frame
        writeIMG('img10_x20_lam100_fs1.dat',self.I2)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,10,1,[6,6],
            frame=True,nidx=20,opt_axis=0)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img10_x20_lam100_fs1.dat')   
        #three time frames 
        writeIMG('img10_x20_lam100_fs1.dat',self.I1)
        writeIMG('img11_x20_lam100_fs1.dat',self.I2)
        writeIMG('img12_x20_lam100_fs1.dat',self.I3)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,3,10,1,[6,6],
            frame=True,nidx=20,opt_axis=0)
        Itrue=np.ones((6,6))*-1
        Itrue[1:4,1:4]=[[0.15,0.25,0.35],[0.45,0.3166666666,0.425],
            [0.75,0.625,0.725]]
        self.assertTrue((abs(IMG-Itrue)<1E-6).all())
        for i in range(10,13):
            os.remove('img'+str(i)+'_x20_lam100_fs1.dat')  
 
        #test 3D slices y          
        writeIMG('img_y20_lam100_fs1.dat',self.I2)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,-1,1,[6,6],
            frame=True,nidx=20,opt_axis=1)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img_y20_lam100_fs1.dat')   
        #one time frame
        writeIMG('img10_y20_lam100_fs1.dat',self.I2)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,10,1,[6,6],
            frame=True,nidx=20,opt_axis=1)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img10_y20_lam100_fs1.dat')   
        #three time frames 
        writeIMG('img10_y20_lam100_fs1.dat',self.I1)
        writeIMG('img11_y20_lam100_fs1.dat',self.I2)
        writeIMG('img12_y20_lam100_fs1.dat',self.I3)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,3,10,1,[6,6],
            frame=True,nidx=20,opt_axis=1)
        Itrue=np.ones((6,6))*-1
        Itrue[1:4,1:4]=[[0.15,0.25,0.35],[0.45,0.3166666666,0.425],
            [0.75,0.625,0.725]]
        self.assertTrue((abs(IMG-Itrue)<1E-6).all())
        for i in range(10,13):
            os.remove('img'+str(i)+'_y20_lam100_fs1.dat')  
 
        #test 3D slices z          
        writeIMG('img_z20_lam100_fs1.dat',self.I2)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,-1,1,[6,6],
            frame=True,nidx=20,opt_axis=2)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img_z20_lam100_fs1.dat')   
        #one time frame
        writeIMG('img10_z20_lam100_fs1.dat',self.I2)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,10,1,[6,6],
            frame=True,nidx=20,opt_axis=2)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img10_z20_lam100_fs1.dat')   
        #three time frames 
        writeIMG('img10_z20_lam100_fs1.dat',self.I1)
        writeIMG('img11_z20_lam100_fs1.dat',self.I2)
        writeIMG('img12_z20_lam100_fs1.dat',self.I3)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,3,10,1,[6,6],
            frame=True,nidx=20,opt_axis=2)
        Itrue=np.ones((6,6))*-1
        Itrue[1:4,1:4]=[[0.15,0.25,0.35],[0.45,0.3166666666,0.425],
            [0.75,0.625,0.725]]
        self.assertTrue((abs(IMG-Itrue)<1E-6).all())

    def test_add_color(self):
        IMGs=np.zeros((3,2,5))
        #Check hue-value addition
        IMGs[0,0,:]=[1,1,0,0,0] #red+yellow
        IMGs[0,1,:]=[1,0,0,0,1] #red+violet
        #Check saturation
        IMGs[1,0,:]=[0.5,0,0,0.1,0] #red+cyan
        #3-color addition
        IMGs[1,1,:]=[0.1,0,0.3,0,0.5] #red+green+violet
        #4-color addition
        IMGs[2,0,:]=[0.9,0.6,0.6,0,0.3] 
        #5-color addition
        IMGs[2,1,:]=[0.5,0.5,0.5,0.9,0.8]

        hues=[0,60,120,180,270] #red, yellow, green, cyan, violet
        newIMGs=siliscopy.plot_image.add_color(IMGs,hues)
        #red+yellow=orange
        self.assertTrue((abs(newIMGs[0,0,:]-np.array([1,0.5,0]))<1E-6).all())
        #red+violet
        self.assertTrue((abs(newIMGs[0,1,:]-np.array([1,0,0.75]))<1E-6).all())
        #red+cyan
        self.assertTrue((abs(newIMGs[1,0,:]-np.array([0.5,0,0]))<1E-6).all())
        #red+green+violet
        self.assertTrue((abs(newIMGs[1,1,:]-np.array([0.221605828,0.1,0.5]))< \
            1E-6).all())
        #4-color
        self.assertTrue((abs(newIMGs[2,0,:]-np.array([0.9,0.7969928058,0.6]))< \
            1E-6).all())
        #5-color
        self.assertTrue((abs(newIMGs[2,1,:]-np.array([0.5,0.9,0.8375138173]))< \
            1E-6).all())

    def test_plot_ism(self): #Checks the shape of images produced.
        #Grey 2d, dynamic
        img0=np.random.random((4,8)) #XY
        sha0=img0.shape
        siliscopy.plot_image.plot_ism(img0,[0.1],[2],1,3,4,filename='out',dpi=2)
        img1=cv2.imread('out3_lam2_fs4_T1_I0.1.jpeg')
        sha1=img1.shape
        self.assertTrue(abs(sha0[1]/sha0[0]-sha1[1]/sha1[0])<1E-6) 

    def test_mono2d(self):
        Img=np.random.rand(4,6)
        Img[0,:]=-1.0
        Img[-1,:]=-1.0
        Img[:,0]=-1.0
        Img[:,-1]=-1.0
        writeIMG('img0_lam1_fs1.dat',Img)
        siliscopy.plot_image.plot_grey_img('img',[1],[1],1,0,1,[4,6],1,0.1,
            otype='tiff8',dlmn=[0.1,0.1,0.1],pbc=[1,1,1])
        data=tif.TiffFile('img0_lam1_fs1_T1_I1.tiff')
        foo=data.pages[0].tags['XResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        foo=data.pages[0].tags['YResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        self.assertEqual(data.imagej_metadata['unit'],'nm')
        self.assertEqual(data.series[0].axes,'YX')
        self.assertEqual(data.imagej_metadata['pbc'],str([1,1,1]))
        self.assertEqual(data.imagej_metadata['bounds'],str([[1,5,1,3]]))
        self.assertEqual(str(data.series[0].shape),str(tuple(Img.shape)))        

    def test_mono2dt(self):
        Img=np.random.rand(2,4,6)
        Img[0,0,:]=-1.0
        Img[0,-1,:]=-1.0
        Img[0,:,0]=-1.0
        Img[0,:,-1]=-1.0

        Img[1,:2,:]=-1.0
        Img[1,-1,:]=-1.0
        Img[1,:,:2]=-1.0
        Img[1,:,-1]=-1.0

        writeIMG('img0_lam1_fs1.dat',Img[0,:,:])
        writeIMG('img10_lam1_fs1.dat',Img[1,:,:])
        siliscopy.plot_image.plot_grey_2dtimg('img',[1],[1],1,0,11,10,1,[4,6],
            [0.1,0.1,0.1],0.5,otype='tiff8')

        data=tif.TiffFile('img0-11_lam1_fs1_T1_I1.tiff')
        foo=data.pages[0].tags['XResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        foo=data.pages[0].tags['YResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        self.assertEqual(data.imagej_metadata['unit'],'nm')
        self.assertEqual(data.series[0].axes,'TYX')
        self.assertEqual(data.imagej_metadata['pbc'],str([1,1,1]))
        self.assertEqual(data.imagej_metadata['bounds'],str([[1,5,1,3],
            [2,5,2,3]]))
        self.assertEqual(data.imagej_metadata['finterval'],0.5)
        self.assertEqual(data.imagej_metadata['funit'],'ns')
        self.assertEqual(str(data.series[0].shape),str(tuple(Img.shape)))        

    def test_mono3d(self):
        Img=np.random.rand(3,4,6)
        Img[:,0,:]=-1.0
        Img[:,-1,:]=-1.0
        Img[:,:,0]=-1.0
        Img[:,:,-1]=-1.0
        Img[0,:,:]=-1.0
        Img[-1,:,:]=-1.0
        for i in range(3):
            writeIMG('img0_z'+str(i)+'_lam1_fs1.dat',Img[i,:,:])

        siliscopy.plot_image.plot_grey_3dimg('img',[1],[1],1,0,1,[4,6],[0.1,0.1,0.2],3,2,otype='tiff8')

        data=tif.TiffFile('img0_z_lam1_fs1_T1_I1.tiff')
        foo=data.pages[0].tags['XResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        foo=data.pages[0].tags['YResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        self.assertEqual(data.imagej_metadata['spacing'],0.2)
        self.assertEqual(data.imagej_metadata['unit'],'nm')
        self.assertEqual(data.series[0].axes,'ZYX')
        self.assertEqual(data.imagej_metadata['pbc'],str([1,1,1]))
        self.assertEqual(data.imagej_metadata['bounds'],str([[1,2,1,5,1,3]]))
        self.assertEqual(str(data.series[0].shape),str(tuple(Img.shape)))        


    def test_mono3dt(self):
        Img=np.random.rand(2,3,4,6)
        Img[:,:,0,:]=-1.0
        Img[:,:,-1,:]=-1.0
        Img[:,:,:,0]=-1.0
        Img[:,:,:,-1]=-1.0
        Img[:,0,:,:]=-1.0
        Img[:,-1,:,:]=-1.0
        for j in range(2):
            for i in range(3):
                writeIMG('img'+str(j*10)+'_z'+str(i)+'_lam1_fs1.dat',Img[j,i,:,:])

        siliscopy.plot_image.plot_grey_3dtimg('img',[1],[1],1,0,11,10,1,[4,6],[0.1,0.1,0.2],3,2,0.5, otype='tiff8')

        data=tif.TiffFile('img0-11_z_lam1_fs1_T1_I1.tiff')
        foo=data.pages[0].tags['XResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        foo=data.pages[0].tags['YResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        self.assertEqual(data.imagej_metadata['spacing'],0.2)
        self.assertEqual(data.imagej_metadata['finterval'],0.5)
        self.assertEqual(data.imagej_metadata['funit'],'ns')
        self.assertEqual(data.imagej_metadata['unit'],'nm')
        self.assertEqual(data.series[0].axes,'TZYX')
        self.assertEqual(data.imagej_metadata['pbc'],str([1,1,1]))
        self.assertEqual(data.imagej_metadata['bounds'],str([[1,2,1,5,1,3],[1,2,1,5,1,3]]))
        self.assertEqual(str(data.series[0].shape),str(tuple(Img.shape)))        

    def test_color2d(self):
        Img=np.random.rand(2,4,6)
        Img[:,0,:]=-1.0
        Img[:,-1,:]=-1.0
        Img[:,:,0]=-1.0
        Img[:,:,-1]=-1.0
        writeIMG('img0_lam1_fs1.dat',Img[0,:,:])
        writeIMG('img0_lam2_fs1.dat',Img[1,:,:])

        siliscopy.plot_image.plot_col_img('img',[1,1],[1,2],[0,120],1,0,1,[4,6],[0.1,0.1,0.1],1,0.1,otype='tiff8',pbc=[1,1,1])

        data=tif.TiffFile('img0_fs1_T1_I_1_1.tiff')
        foo=data.pages[0].tags['XResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        foo=data.pages[0].tags['YResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        self.assertEqual(data.imagej_metadata['unit'],'nm')
        self.assertEqual(data.series[0].axes,'CYX')
        self.assertEqual(data.imagej_metadata['pbc'],str([1,1,1]))
        self.assertEqual(data.imagej_metadata['bounds'],str([[1,5,1,3]]))
        sha=list(Img.shape)
        sha[0]=3
        self.assertEqual(str(data.series[0].shape),str(tuple(sha)))        

    def test_nomix2d(self):
        Img=np.random.rand(2,4,6)
        Img[:,0,:]=-1.0
        Img[:,-1,:]=-1.0
        Img[:,:,0]=-1.0
        Img[:,:,-1]=-1.0
        writeIMG('img0_lam1_fs1.dat',Img[0,:,:])
        writeIMG('img0_lam2_fs1.dat',Img[1,:,:])

        siliscopy.plot_image.plot_col_img('img',[1,1],[1,2],[0,120],1,0,1,[4,6],
            [0.1,0.1,0.1],1,0.1,otype='tiff8',pbc=[1,1,1],mix_type='nomix')

        data=tif.TiffFile('img0_fs1_T1_I_1_1.tiff')
        foo=data.pages[0].tags['XResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        foo=data.pages[0].tags['YResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        self.assertEqual(data.imagej_metadata['unit'],'nm')
        self.assertEqual(data.series[0].axes,'CYX')
        self.assertEqual(data.imagej_metadata['pbc'],str([1,1,1]))
        self.assertEqual(data.imagej_metadata['bounds'],str([[1,5,1,3]]))
        self.assertEqual(str(data.series[0].shape),str(tuple(Img.shape)))        

    def test_color2dt(self):
        Img=np.random.rand(2,2,4,6)
        Img[:,0,0,:]=-1.0
        Img[:,0,-1,:]=-1.0
        Img[:,0,:,0]=-1.0
        Img[:,0,:,-1]=-1.0

        Img[:,1,:2,:]=-1.0
        Img[:,1,-1,:]=-1.0
        Img[:,1,:,:2]=-1.0
        Img[:,1,:,-1]=-1.0

        writeIMG('img0_lam1_fs1.dat',Img[0,0,:,:])
        writeIMG('img10_lam1_fs1.dat',Img[0,1,:,:])
        writeIMG('img0_lam2_fs1.dat',Img[1,0,:,:])
        writeIMG('img10_lam2_fs1.dat',Img[1,1,:,:])
        siliscopy.plot_image.plot_col_2dtimg('img',[1,1],[1,2],[0,120],1,0,11,
            10,1,[4,6],0.5,[0.1,0.1,0.1],otype='tiff8')

        data=tif.TiffFile('img0-11_fs1_T1_I_1_1.tiff')
        foo=data.pages[0].tags['XResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        foo=data.pages[0].tags['YResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        self.assertEqual(data.imagej_metadata['unit'],'nm')
        self.assertEqual(data.series[0].axes,'TCYX')
        self.assertEqual(data.imagej_metadata['pbc'],str([1,1,1]))
        self.assertEqual(data.imagej_metadata['bounds'],str([[1,5,1,3],
            [2,5,2,3]]))
        self.assertEqual(data.imagej_metadata['finterval'],0.5)
        self.assertEqual(data.imagej_metadata['funit'],'ns')
        sha=list(Img.shape)
        sha[1]=3
        self.assertEqual(str(data.series[0].shape),str(tuple(sha)))        

    def test_color3d(self):
        Img=np.random.rand(2,5,4,6) #CZXY
        Img[:,:,0,:]=-1.0
        Img[:,:,-1,:]=-1.0
        Img[:,:,:,0]=-1.0
        Img[:,:,:,-1]=-1.0
        Img[:,0,:,:]=-1.0
        Img[:,-1,:,:]=-1.0
        for i in range(5):
            writeIMG('img0_z'+str(i)+'_lam1_fs1.dat',Img[0,i,:,:])
            writeIMG('img0_z'+str(i)+'_lam2_fs1.dat',Img[1,i,:,:])

        siliscopy.plot_image.plot_col_3dimg('img',[1,1],[1,2],[0,120],1,0,1,
            [4,6],[0.1,0.1,0.2],5,2,otype='tiff8')

        data=tif.TiffFile('img0_z_fs1_T1_I_1_1.tiff')
        foo=data.pages[0].tags['XResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        foo=data.pages[0].tags['YResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        self.assertEqual(data.imagej_metadata['spacing'],0.2)
        self.assertEqual(data.imagej_metadata['unit'],'nm')
        self.assertEqual(data.series[0].axes,'ZCYX')
        self.assertEqual(data.imagej_metadata['pbc'],str([1,1,1]))
        self.assertEqual(data.imagej_metadata['bounds'],str([[1,4,1,5,1,3]]))
        sha=[5,3,4,6]
        self.assertEqual(str(data.series[0].shape),str(tuple(sha)))        


    def test_color3dt(self):
        Img=np.random.rand(2,2,3,4,6)
        Img[:,:,:,0,:]=-1.0
        Img[:,:,:,-1,:]=-1.0
        Img[:,:,:,:,0]=-1.0
        Img[:,:,:,:,-1]=-1.0
        Img[:,:,0,:,:]=-1.0
        Img[:,:,-1,:,:]=-1.0
        for j in range(2):
            for i in range(3):
                writeIMG('img'+str(j*10)+'_z'+str(i)+'_lam1_fs1.dat',
                    Img[0,j,i,:,:])
                writeIMG('img'+str(j*10)+'_z'+str(i)+'_lam2_fs1.dat',
                    Img[1,j,i,:,:])

        siliscopy.plot_image.plot_col_3dtimg('img',[1,1],[1,2],[0,120],1,0,11,
            10,1,[4,6],[0.1,0.1,0.2],3,2,0.5, otype='tiff8')

        data=tif.TiffFile('img0-11_z_fs1_T1_I_1_1.tiff')
        foo=data.pages[0].tags['XResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        foo=data.pages[0].tags['YResolution'].value 
        self.assertAlmostEqual(foo[1]/foo[0],0.1,4)
        self.assertEqual(data.imagej_metadata['spacing'],0.2)
        self.assertEqual(data.imagej_metadata['finterval'],0.5)
        self.assertEqual(data.imagej_metadata['funit'],'ns')
        self.assertEqual(data.imagej_metadata['unit'],'nm')
        self.assertEqual(data.series[0].axes,'TZCYX')
        self.assertEqual(data.imagej_metadata['pbc'],str([1,1,1]))
        self.assertEqual(data.imagej_metadata['bounds'],str([[1,2,1,5,1,3],
            [1,2,1,5,1,3]]))
        sha=[2,3,3,4,6]
        self.assertEqual(str(data.series[0].shape),str(tuple(sha)))        

class TestConvert(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        pass

    @classmethod
    def tearDownClass(self):
        pass

    def setUp(self):
        self.PSF=np.array([[0.,0.,0.,1.0],[0.1,0.,0.,0.829401],
             [0.1,0.1,0.,0.683437],[0.,0.,0.2,0.92737],[0.1,0.,0.2,0.76942],
             [0.1,0.1,0.2,0.634507]])
        self.PSF2=np.array([[0.,0.,-0.2,0.90737],[0.1,0.,-0.2,0.70942],
             [0.1,0.1,-0.2,0.604507],[0.,0.,0.,1.0],[0.1,0.,0.,0.829401],
             [0.1,0.1,0.,0.683437],[0.,0.,0.2,0.92737],[0.1,0.,0.2,0.76942],
             [0.1,0.1,0.2,0.634507]])

    def tearDown(self):
        pass

    def test_dat2tiff(self):
        writePSF('test.dat',self.PSF)
        siliscopy.convert.psf_dat2tiff('test.dat','test.tiff',[0.1,0.1,0.2],
                 [0.1,0.1,0.2],dtype='uint8',psf_type=0)
        psf_dat=tif.imread('test.tiff')
        foo_PSF=(self.PSF*np.iinfo('uint8').max).astype('uint8')
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    x,y=sorted([abs(i-1),abs(j-1)],reverse=True)
                    z=abs(k-1)
                    self.assertEqual(psf_dat[k,j,i],foo_PSF[3*z+x+y,3])
        os.system('rm test.dat')
        os.system('rm test.tiff')

        writePSF('test.dat',self.PSF2)
        siliscopy.convert.psf_dat2tiff('test.dat','test.tiff',[0.1,0.1,0.2],
                 [0.1,0.1,0.2],dtype='uint8',psf_type=1)
        psf_dat=tif.imread('test.tiff')
        foo_PSF=(self.PSF2*np.iinfo('uint8').max).astype('uint8')
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    x,y=sorted([abs(i-1),abs(j-1)],reverse=True)
                    self.assertEqual(psf_dat[k,j,i],foo_PSF[3*k+x+y,3])
        os.system('rm test.dat')
        os.system('rm test.tiff')

    def test_dat2tiff2(self):
        pass #Not tested.

    def test_n2tiff2(self):
        #Check YX -> ZYX
        imgs=np.random.randint(0,255,(4,4,4),dtype='uint8')
        for i in range(4):
            tif.imsave('test'+str(i)+'.tif', imgs[i,:,:], resolution=(0.1,0.1), 
                imagej=True, metadata={'axes':'YX', 'unit':'nm',
                'pbc':[1,1,1], 'bounds':[[0,4,0,4]]})
        f=open('fname.dat','w')
        for i in range(4): 
            f.write('test'+str(i)+'.tif\n')
        f.close()
        siliscopy.convert.nstack2tiff('fname.dat','test_all.tif',2)
        
        img_dat=tif.TiffFile('test_all.tif')
        self.assertEqual(img_dat.series[0].axes,'ZYX')
        self.assertEqual(img_dat.series[0].shape,(4,4,4))
        self.assertEqual(img_dat.imagej_metadata['unit'],'nm')
        self.assertEqual(img_dat.imagej_metadata['spacing'],2)
        foo=(img_dat.imagej_metadata['pbc']=='[1, 1, 1]') or (img_dat.imagej_metadata['pbc']=='[1 1 1]')
        self.assertTrue(foo)
        self.assertEqual(img_dat.imagej_metadata['bounds'],str([[0,4,0,4,0,4]]))
        self.assertEqual(img_dat.pages[0].tags['XResolution'].value,(1,10))
        self.assertEqual(img_dat.pages[0].tags['YResolution'].value,(1,10))
        imgs_r=tif.imread('test_all.tif')
        self.assertTrue((imgs_r==imgs).all())

        #Check CYX -> ZCYX
        imgs=np.random.randint(0,255,(4,3,4,4),dtype='uint8')
        for i in range(4):
            tif.imsave('test'+str(i)+'.tif', imgs[i,:,:,:], resolution=(0.1,0.1), 
                imagej=True, metadata={'axes':'CYX', 'unit':'nm',
                'pbc':[1,1,1], 'bounds':[[0,4,0,4]]})
        f=open('fname.dat','w')
        for i in range(4): 
            f.write('test'+str(i)+'.tif\n')
        f.close()
        siliscopy.convert.nstack2tiff('fname.dat','test_all.tif',1)
        
        img_dat=tif.TiffFile('test_all.tif')
        self.assertEqual(img_dat.series[0].axes,'ZCYX')
        self.assertEqual(img_dat.series[0].shape,(4,3,4,4))
        self.assertEqual(img_dat.imagej_metadata['unit'],'nm')
        self.assertEqual(img_dat.imagej_metadata['spacing'],1)
        foo=(img_dat.imagej_metadata['pbc']=='[1, 1, 1]') or (img_dat.imagej_metadata['pbc']=='[1 1 1]')
        self.assertTrue(foo)
        self.assertEqual(img_dat.imagej_metadata['bounds'],str([[0,4,0,4,0,4]]))
        self.assertEqual(img_dat.pages[0].tags['XResolution'].value,(1,10))
        self.assertEqual(img_dat.pages[0].tags['YResolution'].value,(1,10))
        imgs_r=tif.imread('test_all.tif')
        self.assertTrue((imgs_r==imgs).all())

        #Check TCYX -> TZCYX
        imgs=np.random.randint(0,255,(5,4,3,4,4),dtype='uint8')
        for i in range(4):
            tif.imsave('test'+str(i)+'.tif', imgs[:,i,:,:,:], resolution=(0.1,0.1), 
                imagej=True, metadata={'axes':'TCYX', 'unit':'nm',
                'pbc':[1,1,1], 'bounds':[[0,4,0,4]]*5, 'finterval':0.5, 
                'funit':'ns'})
        f=open('fname.dat','w')
        for i in range(4): 
            f.write('test'+str(i)+'.tif\n')
        f.close()
        siliscopy.convert.nstack2tiff('fname.dat','test_all.tif',3)
        
        img_dat=tif.TiffFile('test_all.tif')
        self.assertEqual(img_dat.series[0].axes,'TZCYX')
        self.assertEqual(img_dat.series[0].shape,(5,4,3,4,4))
        self.assertEqual(img_dat.imagej_metadata['unit'],'nm')
        self.assertEqual(img_dat.imagej_metadata['spacing'],3)
        foo=(img_dat.imagej_metadata['pbc']=='[1, 1, 1]') or (img_dat.imagej_metadata['pbc']=='[1 1 1]')
        self.assertTrue(foo)
        self.assertEqual(img_dat.imagej_metadata['finterval'],0.5)
        self.assertEqual(img_dat.imagej_metadata['funit'],'ns')
        self.assertEqual(img_dat.imagej_metadata['bounds'],str([[0,4,0,4,0,4]]*5))
        self.assertEqual(img_dat.pages[0].tags['XResolution'].value,(1,10))
        self.assertEqual(img_dat.pages[0].tags['YResolution'].value,(1,10))
        imgs_r=tif.imread('test_all.tif')
        self.assertTrue((imgs_r==imgs).all())
        os.system('rm *.tif')
        os.system('rm *.dat')

    def test_t2tiff2(self):
        #Check YX -> TYX
        imgs=np.random.randint(0,255,(5,4,4),dtype='uint8')
        for i in range(5):
            tif.imsave('test'+str(i)+'.tif', imgs[i,:,:], resolution=(0.1,0.1), 
                imagej=True, metadata={'spacing':1, 'axes':'YX', 'unit':'nm',
                'pbc':[1,1,1], 'bounds':[[0,4,0,4]]})
        f=open('fname.dat','w')
        for i in range(5): 
            f.write('test'+str(i)+'.tif\n')
        f.close()
        siliscopy.convert.tstack2tiff('fname.dat',0.5,'test_all.tif')
        
        img_dat=tif.TiffFile('test_all.tif')
        self.assertEqual(img_dat.series[0].axes,'TYX')
        self.assertEqual(img_dat.series[0].shape,(5,4,4))
        self.assertEqual(img_dat.imagej_metadata['finterval'],0.5)
        self.assertEqual(img_dat.imagej_metadata['funit'],'ns')
        self.assertEqual(img_dat.imagej_metadata['unit'],'nm')
        foo=(img_dat.imagej_metadata['pbc']=='[1, 1, 1]') or (img_dat.imagej_metadata['pbc']=='[1 1 1]')
        self.assertTrue(foo)
        self.assertEqual(img_dat.imagej_metadata['bounds'],str([[0,4,0,4]]*5))
        self.assertEqual(img_dat.pages[0].tags['XResolution'].value,(1,10))
        self.assertEqual(img_dat.pages[0].tags['YResolution'].value,(1,10))

        imgs_r=tif.imread('test_all.tif')
        self.assertTrue((imgs_r==imgs).all())
        
        #Check CYX -> TCYX
        imgs=np.random.randint(0,255,(5,3,4,4),dtype='uint8')
        for i in range(5):
            tif.imsave('test'+str(i)+'.tif', imgs[i,:,:,:], resolution=(0.1,0.1), 
                imagej=True, metadata={'spacing':1, 'axes':'CYX', 'unit':'nm',
                'pbc':[1,1,1], 'bounds':[[0,4,0,4]]})
        f=open('fname.dat','w')
        for i in range(5): 
            f.write('test'+str(i)+'.tif\n')
        f.close()
        siliscopy.convert.tstack2tiff('fname.dat',0.5,'test_all.tif')
        
        img_dat=tif.TiffFile('test_all.tif')
        self.assertEqual(img_dat.series[0].axes,'TCYX')
        self.assertEqual(img_dat.series[0].shape,(5,3,4,4))
        self.assertEqual(img_dat.imagej_metadata['finterval'],0.5)
        self.assertEqual(img_dat.imagej_metadata['funit'],'ns')
        self.assertEqual(img_dat.imagej_metadata['unit'],'nm')
        foo=(img_dat.imagej_metadata['pbc']=='[1, 1, 1]') or (img_dat.imagej_metadata['pbc']=='[1 1 1]')
        self.assertTrue(foo)
        self.assertEqual(img_dat.imagej_metadata['bounds'],str([[0,4,0,4]]*5))
        self.assertEqual(img_dat.pages[0].tags['XResolution'].value,(1,10))
        self.assertEqual(img_dat.pages[0].tags['YResolution'].value,(1,10))
        imgs_r=tif.imread('test_all.tif')
        self.assertTrue((imgs_r==imgs).all())

        #Check ZCYX -> TZCYX
        imgs=np.random.randint(0,255,(5,4,3,4,4),dtype='uint8')
        for i in range(5):
            tif.imsave('test'+str(i)+'.tif', imgs[i,:,:,:,:], resolution=(0.1,0.1), 
                imagej=True, metadata={'spacing':1, 'axes':'ZCYX', 'unit':'nm',
                'pbc': [1,1,1], 'bounds':[[0,4,0,4,0,4]]})
        f=open('fname.dat','w')
        for i in range(5): 
            f.write('test'+str(i)+'.tif\n')
        f.close()
        siliscopy.convert.tstack2tiff('fname.dat',0.5,'test_all.tif')
        
        img_dat=tif.TiffFile('test_all.tif')
        self.assertEqual(img_dat.series[0].axes,'TZCYX')
        self.assertEqual(img_dat.series[0].shape,(5,4,3,4,4))
        self.assertEqual(img_dat.imagej_metadata['finterval'],0.5)
        self.assertEqual(img_dat.imagej_metadata['spacing'],1)
        self.assertEqual(img_dat.imagej_metadata['funit'],'ns')
        self.assertEqual(img_dat.imagej_metadata['unit'],'nm')
        foo=(img_dat.imagej_metadata['pbc']=='[1, 1, 1]') or (img_dat.imagej_metadata['pbc']=='[1 1 1]')
        self.assertTrue(foo)
        self.assertEqual(img_dat.imagej_metadata['bounds'],str([[0,4,0,4,0,4]]*5))
        self.assertEqual(img_dat.pages[0].tags['XResolution'].value,(1,10))
        self.assertEqual(img_dat.pages[0].tags['YResolution'].value,(1,10))

        imgs_r=tif.imread('test_all.tif')
        os.system('rm *.tif')
        os.system('rm *.dat')
        self.assertTrue((imgs_r==imgs).all())
    def test_imgs2color(self): 
        # Not tested.
        pass

class TestProp(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        pass

    @classmethod
    def tearDownClass(self):
        pass

    def setUp(self):
        self.I=np.random.random((5,5))
        self.I[0,:]=-np.ones(5)
        self.I[4,:]=-np.ones(5) 
        self.I[:,0]=-np.ones(5)
        self.I[:,4]=-np.ones(5)

    def tearDown(self):
        pass

    def test_maxI(self):
        writeIMG('test_lam1_fs2.dat',self.I)
        Imax=siliscopy.prop.get_maxI('test',[1],2)
        self.assertAlmostEqual(Imax,np.amax(self.I))
        os.system('rm test*.dat')

    def test_I0s(self):
        pass 
    
    def test_hist(self):
        pass
  
    def test_numarea_1(self):
        #Check diagonals are not considered the same
        #Simple particle count is correct
        IMG=np.zeros((10,10),dtype='uint8')
        for i in range(10):
            IMG[i,i]=200
        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'YX', 'bounds':[[0,10,0,10]], 'pbc':[1,1,1], 
                      'unit':'nm'})
        Bin,area=siliscopy.prop.get_num_area('test.tif',0.5,outname='test0')
        for i in range(10):
            self.assertAlmostEqual(area[0][i],0.02)
        os.system('rm *.tif')

    def test_numarea_2(self):
        #Check simple connections
        #Check threshold
        IMG=np.zeros((8,8),dtype='uint8')
        IMG[1,1:3]=200
        IMG[1:3,4]=210
        IMG[1,6:8]=215
        IMG[2,7]=213
        IMG[3,1:3]=230
        IMG[4,1]=225
        IMG[4,4:7]=211
        IMG[6,1:4]=100
        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'YX', 'bounds':[[0,8,0,8]], 'pbc':[1,1,1], 
                      'unit':'nm'})
        Bin,areas=siliscopy.prop.get_num_area('test.tif',0.5,outname='test0')
        Bin_0=np.zeros((8,8),dtype='int')
        Bin_0[1,1:3]=92
        Bin_0[1:3,4]=133
        Bin_0[1,6:8]=174
        Bin_0[2,7]=174
        Bin_0[3,1:3]=215
        Bin_0[4,1]=215
        Bin_0[4,4:7]=255
        self.assertTrue((Bin[0]==Bin_0).all()) 
        self.assertTrue((abs(np.array(areas[0])-np.array([0.04,0.04,0.06,0.06,
            0.06]))<1E-6).all())

    def test_numarea_3(self):
        #Check complex connections
        #Check minimum pixels
        #Check pbc
        IMG=np.zeros((10,9),dtype='uint8')
        IMG[0:2,4]=200
        IMG[3,0:2]=210
        IMG[3,3:6]=215
        IMG[3,7:9]=213
        IMG[5:8,1]=212
        IMG[4:7,5]=211
        IMG[5,3]=210
        IMG[5,8]=210
        IMG[6,3:5]=210
        IMG[7,0]=209
        IMG[7,8]=209
        IMG[8:10,4]=209
        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'YX', 'bounds':[[0,9,0,10]], 'pbc':[1,1,1], 
                      'unit':'nm'})
        Bin,areas=siliscopy.prop.get_num_area('test.tif',0.5,outname='test0',
            min_pix=2)
        Bin_0=np.array([[0,0,0,0,102,0,0,0,0],
                        [0,0,0,0,102,0,0,0,0],
                        [0,0,0,0,0,0,0,0,0],
                        [153,153,0,204,204,204,0,153,153],
                        [0,0,0,0,0,204,0,0,0],
                        [0,255,0,204,0,204,0,0,0],
                        [0,255,0,204,204,204,0,0,0],
                        [255,255,0,0,0,0,0,0,255],
                        [0,0,0,0,102,0,0,0,0],
                        [0,0,0,0,102,0,0,0,0]])
        self.assertTrue((Bin[0]==Bin_0).all()) 
        self.assertTrue((abs(np.array(areas[0])-np.array([0.08,0.08,0.18,
            0.1]))<1E-6).all())

    def test_numarea_4(self):
        #Check complex connections
        #Check minimum pixels
        #Check pbc
        IMG=np.zeros((10,9),dtype='uint8')
        IMG[0:2,4]=200
        IMG[3,0:2]=210
        IMG[3,3:6]=215
        IMG[3,7:9]=213
        IMG[5:8,1]=212
        IMG[4:7,5]=211
        IMG[5,3]=210
        IMG[5,8]=210
        IMG[6,3:5]=210
        IMG[7,0]=209
        IMG[7,8]=209
        IMG[8:10,4]=209
        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'YX', 'bounds':[[0,9,0,10]], 'pbc':[0,0,0],
                      'unit':'nm'})
        Bin,areas=siliscopy.prop.get_num_area('test.tif',0.5,outname='test0',
            min_pix=2)
        Bin_0=np.array([[0,0,0,0,85,0,0,0,0],
                        [0,0,0,0,85,0,0,0,0],
                        [0,0,0,0,0,0,0,0,0],
                        [119,119,0,153,153,153,0,187,187],
                        [0,0,0,0,0,153,0,0,0],
                        [0,221,0,153,0,153,0,0,0],
                        [0,221,0,153,153,153,0,0,0],
                        [221,221,0,0,0,0,0,0,0],
                        [0,0,0,0,255,0,0,0,0],
                        [0,0,0,0,255,0,0,0,0]])
        self.assertTrue((Bin[0]==Bin_0).all()) 
        self.assertTrue((abs(np.array(areas[0])-np.array([0.04,0.04,0.18,0.04,
            0.08,0.04]))<1E-6).all())
        os.system('rm *.tif')

    def test_numarea_5(self):
        # Check different tiff file axes
        # CYX
        IMG=np.zeros((3,10,10),dtype='uint8')
        for i in range(10):
            IMG[0,i,i]=200
        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'CYX', 'bounds':[[0,10,0,10]], 'pbc':[1,1,1], 
                      'unit':'nm'})
        Bin,area=siliscopy.prop.get_num_area('test.tif',0.5,outname='test0',
            col_channel=0)
        for i in range(10):
            self.assertAlmostEqual(area[0][i],0.02)
        Bin_0=np.zeros((10,10),dtype='int')
        for i in range(10):
            Bin_0[i,i]=min((i+1)*int(255*0.8/10.+0.5)+int(0.2*255+0.5),255)
        self.assertTrue((Bin[0]==Bin_0).all())
        os.system('rm *.tif')

        # ZYX
        IMG=np.zeros((10,10,10),dtype='uint8')
        for i in range(10):
            IMG[5,i,i]=200
        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'ZYX', 'bounds':[[0,10,0,10,0,10]], 'pbc':[1,1,1],
                      'unit':'nm', 'spacing':0.1})

        Bin,area=siliscopy.prop.get_num_area('test.tif',0.5,outname='test0',
            ncoor=0.5)
        for i in range(10):
            self.assertAlmostEqual(area[0][i],0.02)
        Bin_0=np.zeros((10,10),dtype='int')
        for i in range(10):
            Bin_0[i,i]=min((i+1)*int(255*0.8/10.+0.5)+int(0.2*255+0.5),255)
        self.assertTrue((Bin[0]==Bin_0).all())
        os.system('rm *.tif')

        # CZYX
        IMG=np.zeros((10,3,10,10),dtype='uint8')
        for i in range(10):
            IMG[5,0,i,i]=200
        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'ZCYX', 'bounds':[[0,10,0,10,0,10]], 'pbc':[1,1,1],
                      'unit':'nm', 'spacing':0.1})

        Bin,area=siliscopy.prop.get_num_area('test.tif',0.5,outname='test0',
            ncoor=0.5,col_channel=0)
        for i in range(10):
            self.assertAlmostEqual(area[0][i],0.02)
        Bin_0=np.zeros((10,10),dtype='int')
        for i in range(10):
            Bin_0[i,i]=min((i+1)*int(255*0.8/10.+0.5)+int(0.2*255+0.5),255)
        self.assertTrue((Bin[0]==Bin_0).all())
        os.system('rm *.tif')

        # TCYX
        IMG=np.zeros((5,3,10,10),dtype='uint8')
        for i in range(10):
            for t in range(5):
                IMG[t,0,i,i]=200
        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'TCYX', 'bounds':[[0,10,0,10]]*5, 'pbc':[1,1,1], 
                      'unit':'nm', 'finterval':0.5, 'funit':'ns'})
        Bin,area=siliscopy.prop.get_num_area('test.tif',0.5,outname='test0',
            col_channel=0)
        Bin_0=np.zeros((5,10,10),dtype='int')
        for t in range(5):
            for i in range(10):
                self.assertAlmostEqual(area[t][i],0.02)
                Bin_0[t,i,i]=min((i+1)*20+51,255)
        self.assertTrue((Bin==Bin_0).all())
        os.system('rm *.tif')

        # TZYX
        IMG=np.zeros((5,10,10,10),dtype='uint8')
        for i in range(10):
            for t in range(5):
                IMG[t,5,i,i]=200
        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'TZYX', 'bounds':[[0,10,0,10,0,10]]*5, 
                      'finterval':0.5, 'funit': 'ns', 'pbc':[1,1,1], 
                      'unit':'nm', 'spacing':0.1})

        Bin,area=siliscopy.prop.get_num_area('test.tif',0.5,outname='test0',
            ncoor=0.5)
        Bin_0=np.zeros((5,10,10),dtype='int')
        for t in range(5):
            for i in range(10):
                self.assertAlmostEqual(area[t][i],0.02)
                Bin_0[t,i,i]=min((i+1)*20+51,255)
        self.assertTrue((Bin==Bin_0).all())
        os.system('rm *.tif')

        # TZCYX
        IMG=np.zeros((5,10,3,10,10),dtype='uint8')
        for i in range(10):
            for t in range(5):
                IMG[t,5,0,i,i]=200
        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'TZCYX', 'bounds':[[0,10,0,10,0,10]]*5, 
                      'finterval':0.5, 'funit':'ns', 'pbc':[1,1,1], 'unit':'nm', 
                      'spacing':0.1})

        Bin,area=siliscopy.prop.get_num_area('test.tif',0.5,outname='test0',
            ncoor=0.5,col_channel=0)
        Bin_0=np.zeros((5,10,10),dtype='int')
        for t in range(5):
            for i in range(10):
                self.assertAlmostEqual(area[t][i],0.02)
                Bin_0[t,i,i]=min((i+1)*20+51,255)
        self.assertTrue((Bin==Bin_0).all())
        os.system('rm *.tif')



    def test_numvol_1(self):
        #Check diagonals are not considered the same
        #Simple particle count is correct
        IMG=np.zeros((10,10,10),dtype='uint8')
        for i in range(10):
            IMG[i,i,i]=200
        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'ZYX', 'bounds':[[0,10,0,10,0,10]], 'pbc':[1,1,1],
                      'unit':'nm', 'spacing':0.3})
        Bin,vol=siliscopy.prop.get_num_vol('test.tif',0.5,outname='test0')
        for i in range(10):
            self.assertAlmostEqual(vol[0][i],0.006)
        os.system('rm *.tif')

    def test_numvol_2(self):
        #Check simple connections
        #Check threshold
        IMG=np.zeros((5,8,8),dtype='uint8')
        IMG[1,1,1:3]=200 #3 voxels 
        IMG[2:3,1,2]=200
        IMG[1,1:3,4]=210 #3 voxels
        IMG[2:3,2,4]=209
        IMG[1,1,6:8]=215 #4 voxels
        IMG[1:3,2,7]=213
        IMG[1,3,1:3]=230 #4 voxels
        IMG[1:3,4,1]=225
        IMG[1,4,4:7]=211 #5 voxels
        IMG[2:4,4,6]=210
        IMG[1,6,1:4]=100 #5 voxels
        IMG[2:4,6,3]=90

        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'ZYX', 'bounds':[[0,5,0,8,0,8]], 'pbc':[1,1,1], 
                      'unit':'nm', 'spacing':0.3})
        Bin,vols=siliscopy.prop.get_num_vol('test.tif',0.5,outname='test0')
        Bin_0=np.zeros((5,8,8),dtype='int')
        Bin_0[1,1,1:3]=92
        Bin_0[2:3,1,2]=92
        Bin_0[1,1:3,4]=133
        Bin_0[2:3,2,4]=133
        Bin_0[1,1,6:8]=174
        Bin_0[1:3,2,7]=174
        Bin_0[1,3,1:3]=215
        Bin_0[1:3,4,1]=215
        Bin_0[1,4,4:7]=255
        Bin_0[2:4,4,6]=255
        self.assertTrue((Bin[0]==Bin_0).all()) 
        self.assertTrue((abs(np.array(vols[0])-np.array([0.018,0.018,0.024,
            0.024,0.03]))<1E-6).all())

    def test_numvol_3(self):
        #Check complex connections
        #Check minimum pixels
        #Check pbc
        IMG=np.zeros((5,10,9),dtype='uint8') #ZXY
        IMG[0,0:2,4]=200 # Particle 1
        IMG[0,8:10,4]=209
        IMG[4,0:2,4]=208
        IMG[0,3,0:2]=210 # Particle 2
        IMG[0,3,7:9]=213
        IMG[3:5,3,8]=212
        IMG[0,3,3:6]=215 # Particle 3
        IMG[0,4:7,5]=211
        IMG[0,5,3]=210
        IMG[0,6,3:5]=210
        IMG[4,5,3:5]=210
        IMG[0,5:8,1]=212 # Particle 4
        IMG[0,7,0]=209 
        IMG[0,7,8]=209
        IMG[1,7,1:4]=208
        IMG[0,5,8]=210 # Particle x
        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'ZYX', 'bounds':[[0,5,0,9,0,10]], 'pbc':[1,1,1], 
                      'unit':'nm', 'spacing':0.3})
        Bin,vols=siliscopy.prop.get_num_vol('test.tif',0.5,outname='test0', 
            min_pix=2)
        Bin_0=np.zeros((5,10,9),dtype='int')
        Bin_0[0,0:2,4]=102 # Particle 1
        Bin_0[0,8:10,4]=102
        Bin_0[4,0:2,4]=102
        Bin_0[0,3,0:2]=153 # Particle 2
        Bin_0[0,3,7:9]=153
        Bin_0[3:5,3,8]=153
        Bin_0[0,3,3:6]=204 # Particle 3
        Bin_0[0,4:7,5]=204
        Bin_0[0,5,3]=204
        Bin_0[0,6,3:5]=204
        Bin_0[4,5,3:5]=204
        Bin_0[0,5:8,1]=255 # Particle 4
        Bin_0[0,7,0]=255
        Bin_0[0,7,8]=255
        Bin_0[1,7,1:4]=255
        self.assertTrue((Bin==Bin_0).all()) 
        self.assertTrue((abs(np.array(vols[0])-np.array([0.036,0.036,0.066,
            0.048]))<1E-6).all())

    def test_numvol_4(self):
        #Check pbc
        IMG=np.zeros((5,5,5),dtype='uint8')
        IMG[0,:,4]=200
        IMG[0,2,4]=0
        IMG[1,2,:]=200
        IMG[1,2,2]=0
        IMG[:,3,2]=200
        IMG[2,3,2]=0

        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'ZYX', 'bounds':[[0,5,0,5,0,5]], 'pbc':[0,0,9], 
                      'unit':'nm', 'spacing':0.3})
        Bin,vols=siliscopy.prop.get_num_vol('test.tif',0.5,outname='test0',
            min_pix=2)
        Bin_0=np.zeros((5,5,5),dtype='int')
        Bin_0[0,:2,4]=85
        Bin_0[:2,3,2]=119
        Bin_0[0,3:,4]=153
        Bin_0[1,2,:2]=187
        Bin_0[1,2,3:]=221
        Bin_0[3:,3,2]=255
        self.assertTrue((Bin==Bin_0).all()) 
        self.assertTrue((abs(np.array(vols[0])-np.array([0.012,0.012,0.012, 
            0.012,0.012,0.012]))<1E-6).all())
        os.system('rm *.tif')

    def test_numvol_5(self):
        # Check different tiff file axes

        # CZYX
        IMG=np.zeros((10,3,10,10),dtype='uint8')
        for i in range(10):
            IMG[i,0,i,i]=200
        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'ZCYX', 'bounds':[[0,10,0,10,0,10]], 
                      'pbc':[1,1,1], 'unit':'nm', 'spacing':0.1})

        Bin,vol=siliscopy.prop.get_num_vol('test.tif',0.5,outname='test0',
            col_channel=0)
        for i in range(10):
            self.assertAlmostEqual(vol[0][i],0.002)
        Bin_0=np.zeros((10,10,10),dtype='int')
        for i in range(10):
            Bin_0[i,i,i]=min((i+1)*20+51,255)
        self.assertTrue((Bin==Bin_0).all())
        os.system('rm *.tif')

        # TZYX
        IMG=np.zeros((5,10,10,10),dtype='uint8')
        for i in range(10):
            for t in range(5):
                IMG[t,i,i,i]=200
        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'TZYX', 'bounds':[[0,10,0,10,0,10]]*5, 
                      'finterval':0.5, 'funit':'ns', 'pbc':[1,1,1], 
                      'unit':'nm', 'spacing':0.1})

        Bin,vol=siliscopy.prop.get_num_vol('test.tif',0.5,outname='test0')
        Bin_0=np.zeros((5,10,10,10),dtype='int')
        for t in range(5):
            for i in range(10):
                self.assertAlmostEqual(vol[t][i],0.002)
                Bin_0[t,i,i,i]=min((i+1)*20+51,255)
        self.assertTrue((Bin==Bin_0).all())
        os.system('rm *.tif')

        # TZCYX
        IMG=np.zeros((5,10,3,10,10),dtype='uint8')
        for i in range(10):
            for t in range(5):
                IMG[t,i,0,i,i]=200
        tif.imwrite('test.tif',IMG,imagej=True,resolution=(1./0.1,1./0.2),
            metadata={'axes':'TZCYX', 'bounds':[[0,10,0,10,0,10]]*5, 
                      'finterval':0.5, 'funit':'ns', 'pbc':[1,1,1], 
                      'unit':'nm', 'spacing':0.1})

        Bin,vol=siliscopy.prop.get_num_vol('test.tif',0.5,outname='test0',
                                           col_channel=0)
        Bin_0=np.zeros((5,10,10,10),dtype='int')
        for t in range(5):
            for i in range(10):
                self.assertAlmostEqual(vol[t][i],0.002)
                Bin_0[t,i,i,i]=min((i+1)*20+51,255)
        self.assertTrue((Bin==Bin_0).all())
        os.system('rm *.tif')



        
    
if __name__ == '__main__':
    unittest.main()
