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

class TestPSFGeneration(unittest.TestCase):
    def setUp(self):
        warnings.simplefilter('ignore', category=ResourceWarning)
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
        print(self.NA)
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
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+str(PSF0[i]))
            i+=1 

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
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+str(PSF0[i]))
            i+=1 
        
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
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+str(PSF0[i]))
            i+=1

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
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+str(PSF0[i]))
            i+=1

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
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+str(PSF0[i]))
            i+=1

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
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+str(PSF0[i]))
            i+=1

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
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+str(PSF0[i]))
            i+=1

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
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+str(PSF0[i]))
            i+=1

##############################################
    def test_Mod_Gandy_r(self): #OPD = 0
        PSF0=[1.]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0222045)  
        PSF0.append(0.000233725)
        siliscopy.gen_psf.psf_Mod_Gandy_sep(self.NA, self.meu, self.meu0, self.t0,
            self.tsO, self.meus, self.tg, self.tg0, self.meug, self.meug0, 
            self.lam_gl, self.dlmn, self.Plmn, self.fs, 'test.dat','w',0)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+str(PSF0[i]))
            i+=1

    def test_Mod_Gandy_n(self): #OPD only due to n'
        PSF0=[0.638634]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0410877)  
        PSF0.append(0.00547789)
        siliscopy.gen_psf.psf_Mod_Gandy_sep(self.NA, self.meu, self.meu0, self.t0,
            self.tsO, self.meus, self.tg, self.tg0, self.meug, self.meug0, 
            self.lam_gl, self.dlmn, self.Plmn, self.fs, 'test.dat','w',1)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+str(PSF0[i]))
            i+=1

    def test_Mod_Gandy_tsO(self): #OPD only due to tsO
        PSF0=[0.997601]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0222105)  
        PSF0.append(0.000287965)
        siliscopy.gen_psf.psf_Mod_Gandy_sep(self.NA, self.meu, self.meu0, self.t0,
            1, self.meus, self.tg, self.tg0, self.meug, self.meug0, 
            self.lam_gl, self.dlmn, self.Plmn, self.fs, 'test.dat','w',0)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+str(PSF0[i]))
            i+=1

    def test_Mod_Gandy_meug(self): #OPD only due to tsO
        PSF0=[0.999982]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0222048)  
        PSF0.append(0.000234126)
        siliscopy.gen_psf.psf_Mod_Gandy_sep(self.NA, self.meu, self.meu0, self.t0,
            self.tsO, self.meus, self.tg, self.tg0, 1.6, self.meug0, 
            self.lam_gl, self.dlmn, self.Plmn, self.fs, 'test.dat','w',0)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+str(PSF0[i]))
            i+=1

    def test_Mod_Gandy_tg(self): #OPD only due to tg
        PSF0=[0.997601]
        #Calculated the definite integral for PSF using WolframAlpha
        PSF0.append(0.0222105)  
        PSF0.append(0.000287965)
        siliscopy.gen_psf.psf_Mod_Gandy_sep(self.NA, self.meu, self.meu0, self.t0,
            self.tsO, self.meus, self.tg0+1, self.tg0, self.meus, self.meus, 
            self.lam_gl, self.dlmn, self.Plmn, self.fs, 'test.dat','w',0)
        f=open('test.dat','r')
        i=0 
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split()
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+str(PSF0[i]))
            i+=1

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
            self.assertAlmostEqual(float(foo[-1]),PSF0[i],6,foo[-1]+' '+str(PSF0[i]))
            i+=1

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
        f.write('    1'+'DC   '+'   '+str(atoms[i])+' '+str(i).ljust(5)+\
                str(int(pos[i,0]*1000)/1000).ljust(8)+\
                str(int(pos[i,1]*1000)/1000).ljust(8)+\
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
        siliscopy.gen_mono.gen_mono_c(['pos.gro','param.dat','psf','out'],silent=True)

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
        atoms,pos,box=[np.array(['A','B']),np.array([[0.15,0.15,0.06],[0.09,0.71,0.11]]),
                       [1.0,1.0,1.0]]
        write_pos('pos.gro',pos,atoms,box)
     
        #Generate image data
        siliscopy.gen_mono.gen_mono_c(['pos.gro','param.dat','psf','out'],silent=True)

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

        print('IMG1',IMG1)
        print('IMG1_0',IMG1_0)
        print('IMG2',IMG2)
        print('IMG2_0',IMG2_0)

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
        atoms,pos,box=[np.array(['A','A']),np.array([[0.15,0.15,-0.15],[0.49,0.69,0.31]]),
                       [1.0,1.0,1.0]]
        write_pos('pos.gro',pos,atoms,box)
     
        #Generate image data
        siliscopy.gen_mono.gen_mono_c(['pos.gro','param.dat','psf','out'],silent=True)

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
        atoms,pos,box=[np.array(['A','B']),np.array([[0.15,0.15,0.21],[0.49,0.69,-0.15]]),
                       [1.0,1.0,1.0]]
        write_pos('pos.gro',pos,atoms,box)
     
        #Generate image data
        siliscopy.gen_mono.gen_mono_c(['pos.gro','param.dat','psf','out'],silent=True)
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
        atoms,pos,box=[np.array(['A','B']),np.array([[0.15,0.15,-0.21],[0.09,0.7,0.25]]),
                       [1.0,1.0,1.0]]
        write_pos('pos.gro',pos,atoms,box)
     
        #Generate image data
        siliscopy.gen_mono.gen_mono_c(['pos.gro','param.dat','psf','out'],silent=True)

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
        atoms,pos,box=[np.array(['A','B']),np.array([[0.15,0.15,-0.21],[0.09,0.71,0.25]]),
                       [1.0,1.0,1.0]]
        write_pos('pos.gro',pos,atoms,box)
     
        #Generate image data
        siliscopy.gen_mono.gen_mono_c(['pos.gro','param.dat','psf','out'],silent=True)

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
    @classmethod
    def setUpClass(self):
        pass

    @classmethod
    def tearDownClass(self):
        pass

    def setUp(self):
        self.I1=np.ones((6,6))*-1
        self.I1[2,2]=0.3
        self.I2=np.ones((6,6))*-1
        self.I2[2:4,2:4]=[[0.1, 0.2],[0.4,0.5]]
        self.I3=np.ones((6,6))*-1
        self.I3[1:4,1:4]=[[0.15, 0.25, 0.35],[0.45,0.55,0.65],[0.75,0.85,0.95]]

    def tearDown(self):
        pass

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
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,-1,1,[6,6],frame=True)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img_lam100_fs1.dat')   
        #one time frame
        writeIMG('img10_lam100_fs1.dat',self.I2)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,10,1,[6,6],frame=True)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img10_lam100_fs1.dat')   
        #three time frames 
        writeIMG('img10_lam100_fs1.dat',self.I1)
        writeIMG('img11_lam100_fs1.dat',self.I2)
        writeIMG('img12_lam100_fs1.dat',self.I3)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,3,10,1,[6,6],frame=True)
        Itrue=np.ones((6,6))*-1
        Itrue[1:4,1:4]=[[0.15,0.25,0.35],[0.45,0.3166666666,0.425],[0.75,0.625,0.725]]
        self.assertTrue((abs(IMG-Itrue)<1E-6).all())
        for i in range(10,13):
            os.remove('img'+str(i)+'_lam100_fs1.dat')   

        #test 3D slices x          
        writeIMG('img_x20_lam100_fs1.dat',self.I2)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,-1,1,[6,6],frame=True,nidx=20,opt_axis=0)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img_x20_lam100_fs1.dat')   
        #one time frame
        writeIMG('img10_x20_lam100_fs1.dat',self.I2)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,10,1,[6,6],frame=True,nidx=20,opt_axis=0)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img10_x20_lam100_fs1.dat')   
        #three time frames 
        writeIMG('img10_x20_lam100_fs1.dat',self.I1)
        writeIMG('img11_x20_lam100_fs1.dat',self.I2)
        writeIMG('img12_x20_lam100_fs1.dat',self.I3)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,3,10,1,[6,6],frame=True,nidx=20,opt_axis=0)
        Itrue=np.ones((6,6))*-1
        Itrue[1:4,1:4]=[[0.15,0.25,0.35],[0.45,0.3166666666,0.425],[0.75,0.625,0.725]]
        self.assertTrue((abs(IMG-Itrue)<1E-6).all())
        for i in range(10,13):
            os.remove('img'+str(i)+'_x20_lam100_fs1.dat')  
 
        #test 3D slices y          
        writeIMG('img_y20_lam100_fs1.dat',self.I2)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,-1,1,[6,6],frame=True,nidx=20,opt_axis=1)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img_y20_lam100_fs1.dat')   
        #one time frame
        writeIMG('img10_y20_lam100_fs1.dat',self.I2)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,10,1,[6,6],frame=True,nidx=20,opt_axis=1)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img10_y20_lam100_fs1.dat')   
        #three time frames 
        writeIMG('img10_y20_lam100_fs1.dat',self.I1)
        writeIMG('img11_y20_lam100_fs1.dat',self.I2)
        writeIMG('img12_y20_lam100_fs1.dat',self.I3)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,3,10,1,[6,6],frame=True,nidx=20,opt_axis=1)
        Itrue=np.ones((6,6))*-1
        Itrue[1:4,1:4]=[[0.15,0.25,0.35],[0.45,0.3166666666,0.425],[0.75,0.625,0.725]]
        self.assertTrue((abs(IMG-Itrue)<1E-6).all())
        for i in range(10,13):
            os.remove('img'+str(i)+'_y20_lam100_fs1.dat')  
 
        #test 3D slices z          
        writeIMG('img_z20_lam100_fs1.dat',self.I2)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,-1,1,[6,6],frame=True,nidx=20,opt_axis=2)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img_z20_lam100_fs1.dat')   
        #one time frame
        writeIMG('img10_z20_lam100_fs1.dat',self.I2)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,1,10,1,[6,6],frame=True,nidx=20,opt_axis=2)
        self.assertTrue((abs(IMG-self.I2)<1E-6).all())
        os.remove('img10_z20_lam100_fs1.dat')   
        #three time frames 
        writeIMG('img10_z20_lam100_fs1.dat',self.I1)
        writeIMG('img11_z20_lam100_fs1.dat',self.I2)
        writeIMG('img12_z20_lam100_fs1.dat',self.I3)
        IMG=siliscopy.plot_image.get_grey_img('img',1,100,3,10,1,[6,6],frame=True,nidx=20,opt_axis=2)
        Itrue=np.ones((6,6))*-1
        Itrue[1:4,1:4]=[[0.15,0.25,0.35],[0.45,0.3166666666,0.425],[0.75,0.625,0.725]]
        self.assertTrue((abs(IMG-Itrue)<1E-6).all())
        for i in range(10,13):
            os.remove('img'+str(i)+'_z20_lam100_fs1.dat')  

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
        print(newIMGs)
        #red+yellow=orange
        self.assertTrue((abs(newIMGs[0,0,:]-np.array([1,0.5,0]))<1E-6).all())
        #red+violet
        self.assertTrue((abs(newIMGs[0,1,:]-np.array([1,0,0.75]))<1E-6).all())
        #red+cyan
        self.assertTrue((abs(newIMGs[1,0,:]-np.array([0.5,0,0]))<1E-6).all())
        #red+green+violet
        self.assertTrue((abs(newIMGs[1,1,:]-np.array([0.221605828,0.1,0.5]))<1E-6).all())
        #4-color
        self.assertTrue((abs(newIMGs[2,0,:]-np.array([0.9,0.7969928058,0.6]))<1E-6).all())
        #5-color
        self.assertTrue((abs(newIMGs[2,1,:]-np.array([0.5,0.9,0.8375138173]))<1E-6).all())

    def test_plot_ism(self): #Checks the shape of images produced.
        #Grey 2d, dynamic
        img0=np.random.random((2,4))
        sha0=img0.shape
        siliscopy.plot_image.plot_ism(img0,[0.1],[2],1,3,4,filename='out',dpi=1)
        #
        #def plot_ism(IMG, lam_I0, lam, T, ti, fs, img_hei=3.0, filename=None, 
        #   show=False, gcolmap='gray', dpi=600, otype='jpeg', psf_type=0, tsO=None,
        #   frame_col=1.0):
        #
        img1=cv2.imread('out3_lam2_fs4_T1_I0.1.jpeg')
        sha1=img1.shape
        print(sha0,sha1)
        self.assertTrue(abs(sha0[1]/sha0[0]-sha1[1]/sha1[0])<1E-6) 
        os.system('rm -f *.jpeg')
    
    #    #Grey 2d, static: Expected to be the same. Removed to keep test.py fast.
    #    siliscopy.plot_image.plot_ism(img0,[0.1],[2],1,-1,4,filename='out',dpi=1)
    #    img1=cv2.imread('out_lam2_fs4_I0.1.jpeg')
    #    sha1=img1.shape
    #    print(sha0,sha1)
    #    self.assertTrue(abs(sha0[1]/sha0[0]-sha1[1]/sha1[0])<1E-6) 
    #    os.system('rm -f *.jpeg')
    #    #Color 2d, dynamic
    #    siliscopy.plot_image.plot_ism(img0,[0.1,0.2,0.3],[2,3,4],1,5,6,filename='out',dpi=1)
    #    img1=cv2.imread('out5_fs6_T1_I_0.1_0.2_0.3.jpeg')
    #    sha1=img1.shape
    #    print(sha0,sha1)
    #    self.assertTrue(abs(sha0[1]/sha0[0]-sha1[1]/sha1[0])<1E-6) 
    #    os.system('rm -f *.jpeg')
    #    #Color 2d, static
    #    siliscopy.plot_image.plot_ism(img0,[0.1,0.2,0.3],[2,3,4],1,-1,6,filename='out',dpi=1)
    #    img1=cv2.imread('out_fs6_I_0.1_0.2_0.3.jpeg')
    #    sha1=img1.shape
    #    print(sha0,sha1)
    #    self.assertTrue(abs(sha0[1]/sha0[0]-sha1[1]/sha1[0])<1E-6) 
    #    os.system('rm -f *.jpeg')
        #Grey 2d, dynamic, tiff
        siliscopy.plot_image.plot_ism(img0,[0.1],[2],1,3,4,filename='out',dpi=1, otype='tif8')
        img1=cv2.imread('out3_lam2_fs4_T1_I0.1.tif')
        sha1=img1.shape
        print(sha0,sha1)
        self.assertTrue(abs(sha0[1]/sha0[0]-sha1[1]/sha1[0])<1E-6) 
        os.system('rm -f *.tif')

        #color 2d, dynamic, tiff
        siliscopy.plot_image.plot_ism(img0,[0.1,0.2],[2,3],1,4,5,filename='out',dpi=1, otype='tif8')
        img1=cv2.imread('out4_fs5_T1_I_0.1_0.2.tif')
        sha1=img1.shape
        print(sha0,sha1)
        self.assertTrue(abs(sha0[1]/sha0[0]-sha1[1]/sha1[0])<1E-6) 
        os.system('rm -f *.tif')
    
if __name__ == '__main__':
    unittest.main()
