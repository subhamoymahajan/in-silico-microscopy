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

from optparse import OptionParser
from siliscopy.gen_psf import *
from siliscopy.gen_mono import *
from siliscopy.plot_image import *
from siliscopy.create_vid import *
from siliscopy.prop import *
from siliscopy.convert import *
small=1E-10
import numpy as np

def main():
    parser = OptionParser()
    parser.add_option('-f', '--file', dest="filename", metavar="GRO", 
                      type="str",
                      help="Provide Gromacs containing 1 timestep.")
    parser.add_option('-p', '--paramfile', dest="pfile", metavar="FILE", 
                      type="str", 
                      help="Provide files containing all the parameters")
    parser.add_option('-q', '--psf', dest="psfheader", metavar="FILE", 
                      type="str", help="Starting name of the PSF file")
    parser.add_option('-o', '--output', dest="outname", metavar="FILE", 
                      type="str", help="Provide output filename")
    parser.add_option('-s', '--multiprocess', dest="mprocess", 
                      action="store_true", default=False, 
                      help="Use for multiprocessing")
    parser.add_option('-m', '--method', dest="method",type="str",
                      help="Method of PSF calculation")
    parser.add_option('-c', '--calc', dest="calc", type="str",
                      help="Calculate for a 'specific' n' or 'all'")
    parser.add_option('-d', '--data', dest="data", type="str", default=None,
                      help="Additional data file")
    parser.add_option('-n', '--nid', dest="nid", type="int",
                      help="The index of n' coordinate. int(n'/dn)")
    parser.add_option('-t', '--timestep', dest="timestep", type="int",
                      help="The timestep of the simulation.")
    parser.add_option('-a','--threshold', dest="threshold", type="float",
                      help="Thresold for intensity")    
    parser.add_option('-b','--lambdaID', dest="lambda_ID", type="int",
                       help="Index of the wavelength")    
    parser.add_option('-e','--type', dest="type", type="str",
                       help="Data type", default='jpeg')    
    parser.add_option('-i','--iterations', dest="iters", type="int",
                       help="Number of iterations", default=1)    


    options, remainder = parser.parse_args()
    params={}
    #################DEFAULT VALUES ####################
    params['lam']=np.zeros(10,dtype=int)
    params['hue']=np.zeros(10)
    params['I0']=np.zeros(10)
    params['fourcc']='mp4v'
    params['fpns']=1.0
    params['vid_ext']='.mov'
    params['meu']=1.515
    params['meu0']=1.515
    params['t0']=300.0 #nm
    params['meug']=1.522
    params['meug0']=1.522
    params['tg']=320.0 #nm
    params['tg0']=320.0 #nm
    params['tsO']=0 #nm
    params['meus']=1.33 
    params['NA']=1
    params['add_n']=1
    params['fs']=530
    params['psf_type']=0
    params['poi']=None
    params['gauss']=None
    params['frame_col']=1.0
    params['mix_type']='mt'
    params['scale']=None
    params['tbegin']=None
    params['tmax']=None
    params['tdiff']=None
    #####################################################
    #Read parameters
    if options.pfile !=None:
        f=open(options.pfile,'r')
        for lines in f:
            if '=' not in lines:
                continue
            foo=lines.split('=')
            varname=foo[0].strip()
            val_string=foo[1].split('/')[0].strip()
            val_string=val_string.split()
            
            if varname in ["fs", "T", "dpi", "tbegin", "tmax", "tdiff", \
                "opt_axis", "add_n", "psf_type", "min_pix", "pos_prec"]:
                params[varname]=int(val_string[0])
            elif varname in ["NA", "meu", "scale", "meu0", "t0", "meug", \
                "meug0", "tg", "tg0", "tsO", "meus", "poi", "gauss",\
                "frame_col", "focus_cor", "sig_r", "sig_n", "fpns", "fcs_tmax"]:
                params[varname]=float(val_string[0])
            elif varname == 'fccs_cols':
                val_string=[int(x) for x in val_string]
                params[varname]=np.array(val_string)
            elif varname in ["dlmn", "Plmn","maxlen", "fcsfit_bounds"]:
                val_string=[float(x) for x in val_string]
                params[varname]=np.array(val_string)
            elif varname[0:3]=='lam':
                if varname[3].isdigit():
                    num=int(varname[3:])
                    params['lam'][num-1]=int(val_string[0])
                if varname[3:7]=='_hue':
                    num=int(varname[7:])
                    params['hue'][num-1]=float(val_string[0])
                if varname[3:7]=='_I0_':
                    num=int(varname[7:])
                    params['I0'][num-1]=float(val_string[0])
            elif varname in ["vid_ext", "mix_type", "pp_file"]:
                params[varname]=val_string[0].strip()
            elif varname=='fourcc':
                foo=val_string[0].strip()
                params[varname]=foo[1:-1]
            elif varname=='pbc':
                params[varname]=[0,0,0]
                if 'x' in val_string[0]:
                    params[varname][0]=1
                if 'y' in val_string[0]:
                    params[varname][1]=1
                if 'z' in val_string[0]:
                    params[varname][2]=1
        f.close()
        #Reduce the lam variable
        for i in range(10):
            if params['lam'][i]<1:
                break
        params['lam']=params['lam'][:i]
        params['hue']=params['hue'][:i]
        params['I0']=params['I0'][:i]
        if 'maxlen' in params:
            MaxBox=[0,0]
            MaxBox[0]=int(params['maxlen'][(params['opt_axis']+1)%3]/ \
                          params['dlmn'][0]+0.5)
            MaxBox[1]=int(params['maxlen'][(params['opt_axis']+2)%3]/ \
                          params['dlmn'][1]+0.5)
            params['nmax']=int(params['maxlen'][params['opt_axis']]/ \
                               params['dlmn'][2]+0.5)
    
    if remainder[0]=='gen_psf':

        if options.method in ["gandy", "Gandy", "GL1991", "GL1992", 
            "GibsonLanni", "Mod_Gandy"]: 
            print('NA = '+str(params['NA']))
            print('meu immersion oil = '+str(params['meu']))
            print('fs = '+str(params['fs']))

        if options.method in ["GL1991", "GL1992", "GibsonLanni", "Mod_Gandy"]: 
            print("meu immersion oil (design) = "+str(params["meu0"])) 
            print("thickness immersion oil (design) = "+str(params["t0"])+" nm")
            print("meu coverslip = "+str(params["meug"])) 
            print("meu coverslip (design) = "+str(params["meug0"])) 
            print("thickness coverslip = "+str(params["tg"])+" nm")
            print("thickness coverslip (design) = "+str(params["tg0"])+" nm")
            print("distance between object focal plane and coverslip = "+ \
                str(params["tsO"])+" nm")

        for key in ["lam", "dlmn", "Plmn"]:
            print(key+' = '+str(params[key]))
            if len(params[key])==0:
                raise Exception('Empty array')

        if options.method in ["gauss", "Gauss"]:
            print('sigma_r = '+str(params['sig_r']))
            print('sigma_n = '+str(params['sig_n']))

            for lambd in params['lam']:
                psf_gauss(params['sig_r'], params['sig_n'], params['dlmn'], 
                    params['Plmn'], options.outname, lambd, params['fs'])
        elif options.method in ["gandy", "Gandy"]:
            if options.calc=='all': #Calculates for all n' coordiantes
                if options.mprocess==True: #Multiprocessing. 
                    for lambd in params['lam']:
                        psf_gandy_mp(params['NA'], params['meu'], lambd, 
                            params['dlmn'], params['Plmn'], params['fs'],
                            options.outname)
                else: # Calculates Serially 
                    for lambd in params['lam']:
                        psf_gandy(params['NA'], params['meu'], lambd, 
                                  params['dlmn'], params['Plmn'], 
                                  params['fs'], options.outname)
            elif options.calc in ["specific", "spec"]: #Calculates for one n'
                for lambd in params['lam']:
                    outname = options.outname + '_lam' + str(lambd) + \
                              '_fs' + str(params['fs']) + '_' + \
                              str(options.nid) + '.dat'
                    psf_gandy_sep(params['NA'], params['meu'], lambd, 
                        params['dlmn'], params['Plmn'], params['fs'],
                        outname+'_lam'+str(lambd),'w',options.nid)

        elif options.method in ["GL1991", "GL1992", "GibsonLanni"]:
            if options.calc=='all': #Calculates for all n' coordiantes
                if options.mprocess==True: #Multiprocessing. 
                    for lambd in params['lam']:
                        psf_GL1991_mp(params['NA'], params['meu'], 
                            params['meu0'], params['t0'], params['tsO'], 
                            params['meus'], params['tg'], params['tg0'], 
                            params['meug'], params['meug0'], lambd, 
                            params['dlmn'], params['Plmn'], params['fs'],
                            options.outname)
                else: # Calculates Serially 
                    for lambd in params['lam']:
                        psf_GL1991(params['NA'], params['meu'], 
                            params['meu0'], params['t0'], params['tsO'], 
                            params['meus'], params['tg'], params['tg0'], 
                            params['meug'], params['meug0'], lambd, 
                            params['dlmn'], params['Plmn'], params['fs'], 
                            options.outname)
            elif options.calc in ["specific", "spec"]: #Calculates for one n'
                for lambd in params['lam']:
                    outname = options.outname + '_lam' + str(lambd) + \
                              '_fs' + str(params['fs']) + '_' + \
                              str(options.nid) + '.dat'
                    psf_GL1991_sep(params['NA'], params['meu'], 
                        params['meu0'], params['t0'], params['tsO'], 
                        params['meus'], params['tg'], params['tg0'], 
                        params['meug'], params['meug0'], lambd, params['dlmn'], 
                        params['Plmn'], params['fs'],
                        options.outname+'_lam'+str(lambd),'w',options.nid)

        elif options.method in ["Mod_Gandy"]:
            if options.calc=='all': #Calculates for all n' coordiantes
                if options.mprocess==True: #Multiprocessing. 
                    for lambd in params['lam']:
                        psf_Mod_Gandy_mp(params['NA'], params['meu'], 
                            params['meu0'], params['t0'], params['tsO'], 
                            params['meus'], params['tg'], params['tg0'], 
                            params['meug'], params['meug0'], lambd, 
                            params['dlmn'], params['Plmn'], params['fs'],
                            options.outname)
                else: # Calculates Serially 
                    for lambd in params['lam']:
                        psf_Mod_Gandy(params['NA'], params['meu'], 
                            params['meu0'], params['t0'], params['tsO'], 
                            params['meus'], params['tg'], params['tg0'], 
                            params['meug'], params['meug0'], lambd, 
                            params['dlmn'], params['Plmn'], params['fs'],
                            options.outname)
            elif options.calc in ["specific", "spec"]: #Calculates for one n'
                for lambd in params['lam']:
                    outname = options.outname + '_lam' + str(lambd) + \
                              '_fs' + str(params['fs']) + '_' + \
                              str(options.nid) + '.dat'
                    psf_Mod_Gandy_sep(params['NA'], params['meu'], 
                        params['meu0'], params['t0'], params['tsO'], 
                        params['meus'], params['tg'], params['tg0'], 
                        params['meug'], params['meug0'], lambd, params['dlmn'],
                        params['Plmn'], params['fs'],
                        options.outname+'_lam'+str(lambd),'w', options.nid)

    elif remainder[0]=='gen_spm':
        print("Precision of positions = "+str(params['pos_prec']))
        print("Photophysics file = "+str(params['pp_file']))
        os.system(os.path.dirname(__file__) +'/photo_phys -f '+options.data+' -p '+options.pfile+' -pp '+params['pp_file']+' -o '+options.outname)
         
    elif remainder[0]=='gen_mono':
        if options.method == None:
            options.method='slice'
        if options.data!=None: #Data available in a file
            if options.mprocess==True: #Use multiprocessing
                # siliscopy gen_mono -d [data] -s
                gen_mono_c_mp(options.data,True,False)
            else:#Run serially
                # siliscopy gen_mono -d [data]
                gen_mono_c_serial(options.data,True,False)
        elif options.method=='volume':
            gen_mono_c_vol([options.filename, options.pfile, 
                options.psfheader,options.outname],True, False,
                 params['maxlen'],params['opt_axis'], params['dlmn'], 
                add_n=params['add_n'], mprocess=options.mprocess)
        elif options.method=='slice': #Generates one image slice.
            # siliscopy gen_mono -f [gro] -p [param] -o [out]
            gen_mono_c([options.filename, options.pfile, options.psfheader, 
                        options.outname],False,False)
    elif remainder[0]=='gen_mono_pp':
        if options.method == None:
            options.method='slice'
        if options.data!=None: #Data available in a file
            if options.mprocess==True: #Use multiprocessing
                # siliscopy gen_mono_pp -d [data] -s
                gen_mono_c_mp(options.data,True, True)
            else:#Run serially
                # siliscopy gen_mono_pp -d [data]
                gen_mono_c_serial(options.data,False, True)
        elif options.method=='volume':
            gen_mono_c_vol([options.filename, options.pfile, 
                options.psfheader,options.outname], True, True, 
                params['maxlen'], params['opt_axis'], params['dlmn'], 
                add_n=params['add_n'], mprocess=options.mprocess)
        elif options.method=='slice': #Generates one image slice.
            # siliscopy gen_mono -f [gro] -p [param] -o [out]
            gen_mono_c([options.filename, options.pfile, options.psfheader, 
                        options.outname], False, True)
    
    elif remainder[0]=='plot':
        filename=options.filename
        foo=list(params['pbc'])*2
        pbc_lmn=np.array(foo[params['opt_axis']+1:params['opt_axis']+4])
        for key in ["lam","I0","T","fs","dlmn","scale","dpi"]:
            if key not in params:
                continue
            print(key+" = "+str(params[key]))

        if options.outname==None:
            outname=options.filename
        else:
            outname=options.outname
 
        print("maxlen = "+str(params['maxlen']))
        Bm=params['maxlen'][(params['opt_axis']+1)%3] 
        noise=False
        if 'dpi' not in params:
            params['dpi']=600
            
        if options.method in ["mono", "grey", "gray", "noise_mono", \
            "noise_grey", "noise_gray"]:
            if options.method[:6]=="noise_":
                noise=True
            if options.calc=='show': #show specific
                plot_grey_img(filename, params['I0'], params['lam'],
                    params['T'], options.timestep, params['fs'], MaxBox, Bm, 
                    params['scale'], params['dpi'], noise=noise, 
                    poi=params['poi'], gauss=params['gauss'], 
                    psf_type=params['psf_type'], tsO=params['tsO'],
                    frame_col=params['frame_col'], dlmn=params['dlmn'],
                    pbc=pbc_lmn)
    
            elif options.calc in ["specific", "spec"]: #Specific output
                plot_grey_img(filename, params['I0'], params['lam'],
                    params['T'], options.timestep, params['fs'], MaxBox, Bm, 
                    params['scale'], params['dpi'], outname, noise=noise, 
                    poi=params['poi'], gauss=params['gauss'], 
                    otype=options.type, psf_type=params['psf_type'], 
                    tsO=params['tsO'], frame_col=params['frame_col'],
                    dlmn=params['dlmn'],pbc=pbc_lmn)
            
            elif options.calc == 'all':
                print('tbegin = '+str(params['tbegin']))
                print('tmax = '+str(params['tmax']))
                print('tdiff = '+str(params['tdiff']))
                if options.mprocess==True: #parallel
                    plot_grey_mp(filename, params['I0'], params['lam'],
                        params['T'], params['tbegin'], params['tmax'],
                        params['tdiff'], params['fs'], MaxBox, Bm, 
                        params['scale'], params['dpi'], outname, 
                        params['frame_col'], noise, params['poi'],
                        params['gauss'], options.type, params['psf_type'],
                        params['tsO'], params['dlmn'],pbc_lmn)
                else: #Serial
                    plot_grey_serial(filename, params['I0'], params['lam'],
                        params['T'], params['tbegin'], params['tmax'], 
                        params['tdiff'], params['fs'], MaxBox, Bm, 
                        params['scale'], params['dpi'], outname, 
                        params['frame_col'], noise, params['poi'], 
                        params['gauss'], options.type, params['psf_type'],
                        params['tsO'], params['dlmn'], pbc_lmn)

        elif options.method in ["col", "color", "noise_col", "noise_color"]:
            print("hue = "+str(params["hue"]))
            if options.method[:6]=='noise_':
                noise=True
            if options.calc=='show': #show specific
                plot_col_img(filename, params['I0'], params['lam'],
                    params['hue'], params['T'], options.timestep, params['fs'],
                    MaxBox, params['dlmn'], Bm, params['scale'], params['dpi'], 
                    frame_col=params['frame_col'], mix_type=params['mix_type'],
                    noise=noise, poi=params['poi'], tsO=params['tsO'], 
                    gauss=params['gauss'], psf_type=params['psf_type'],
                    pbc=pbc_lmn) 
                                                 
            elif options.calc in ["specific", "s pec"]: #Save Specific
                plot_col_img(filename, params['I0'], params['lam'],
                    params['hue'], params['T'], options.timestep, params['fs'],
                    MaxBox, params['dlmn'], Bm, params['scale'], params['dpi'], 
                    outfile=outname, otype=options.type, noise=noise, 
                    poi=params['poi'], gauss=params['gauss'], 
                    psf_type=params['psf_type'], tsO=params['tsO'],
                    mix_type=params['mix_type'], frame_col=params['frame_col'],
                    pbc=pbc_lmn)
            
            elif options.calc == 'all':
                print('tbegin = '+str(params['tbegin']))
                print('tmax = '+str(params['tmax']))
                print('tdiff = '+str(params['tdiff']))
                if options.mprocess==True: #parallel
                    plot_col_mp(filename, params['I0'], params['lam'],
                        params['hue'], params['T'], params['tbegin'], 
                        params['tmax'], params['tdiff'], params['fs'], MaxBox, 
                        Bm, params['scale'], params['dpi'], outname, 
                        params['frame_col'], params['mix_type'],noise, 
                        params['poi'], params['gauss'], options.type,
                        params['psf_type'], params['tsO'], params['dlmn'],
                        pbc_lmn)

                else: #Serial
                    plot_col_serial(filename, params['I0'], params['lam'],
                        params['hue'], params['T'], params['tbegin'], 
                        params['tmax'], params['tdiff'], params['fs'], MaxBox,
                        Bm, params['scale'], params['dpi'], outname, 
                        params['frame_col'], params['mix_type'], noise, 
                        params['poi'], params['gauss'], options.type, 
                        params['psf_type'], params['tsO'], params['dlmn'],
                        pbc_lmn)

        elif options.method in ["mono2dt", "grey2dt", "gray2dt", \
            "noise_mono2dt", "noise_grey2dt", "noise_gray2dt"]:
            print('tbegin = '+str(params['tbegin']))
            print('tmax = '+str(params['tmax']))
            print('tdiff = '+str(params['tdiff']))
             
            if options.method[:6]=="noise_":
                noise=True
            plot_grey_2dtimg(filename, params['I0'], params['lam'], params['T'],
                params['tbegin'], params['tmax'], params['tdiff'], params['fs'],
                MaxBox, params['dlmn'], params['fpns'], outfile=outname, 
                mprocess=options.mprocess, noise=noise, poi=params['poi'], 
                otype=options.type, gauss=params['gauss'], 
                psf_type=params['psf_type'], tsO=params['tsO'], pbc=pbc_lmn)

        elif options.method in ["col2dt", "color2dt", "noise_col2dt", 
            "noise_color2dt"]:
            print('tbegin = '+str(params['tbegin']))
            print('tmax = '+str(params['tmax']))
            print('tdiff = '+str(params['tdiff']))
             
            if options.method[:6]=="noise_":
                noise=True
            plot_col_2dtimg(filename, params['I0'], params['lam'],params['hue'],
                params['T'], params['tbegin'], params['tmax'], params['tdiff'], 
                params['fs'], MaxBox, params['fpns'], params['dlmn'], 
                noise=noise, outfile=outname, otype=options.type, pbc=pbc_lmn,
                mprocess=options.mprocess, poi=params['poi'], tsO=params['tsO'],
                gauss=params['gauss'], psf_type=params['psf_type'], 
                mix_type=params['mix_type'], )


        elif options.method in ["mono3d", "grey3d", "gray3d", "noise_mono3d", \
            "noise_grey3d", "noise_gray3d"]:
             
            if options.method[:6]=="noise_":
                noise=True
            plot_grey_3dimg(filename, params['I0'], params['lam'], params['T'],
                options.timestep, params['fs'], MaxBox, params['dlmn'], 
                params['nmax'], params['opt_axis'], add_n=params['add_n'], 
                outfile=outname, poi=params['poi'], tsO=params['tsO'],
                otype=options.type, mprocess=options.mprocess, noise=noise,
                gauss=params['gauss'], psf_type=params['psf_type'], pbc=pbc_lmn) 

        elif options.method in ["col3d", "color3d", "noise_col3d", 
            "noise_color3d"]:
             
            if options.method[:6]=="noise_":
                noise=True
            plot_col_3dimg(filename, params['I0'], params['lam'], params['hue'],
                params['T'], options.timestep, params['fs'], MaxBox, 
                params['dlmn'], params['nmax'], params['opt_axis'], 
                add_n=params['add_n'], outfile=outname, noise=noise, 
                otype=options.type, mprocess=options.mprocess, 
                poi=params['poi'], gauss=params['gauss'],
                psf_type=params['psf_type'], tsO=params['tsO'],pbc=pbc_lmn)

        elif options.method in ["mono3dt", "grey3dt", "gray3dt", 
            "noise_mono3dt", "noise_grey3dt", "noise_gray3dt"]:
            print('tbegin = '+str(params['tbegin']))
            print('tmax = '+str(params['tmax']))
            print('tdiff = '+str(params['tdiff']))
             
            if options.method[:6]=="noise_":
                noise=True
            plot_grey_3dtimg(filename, params['I0'], params['lam'], params['T'],
                params['tbegin'], params['tmax'], params['tdiff'], params['fs'],
                MaxBox, params['dlmn'], params['nmax'], params['opt_axis'],
                params['fpns'], add_n=params['add_n'], noise=noise, 
                outfile=outname, otype=options.type, mprocess=options.mprocess, 
                poi=params['poi'], gauss=params['gauss'], 
                psf_type=params['psf_type'], tsO=params['tsO'], pbc=pbc_lmn)

        elif options.method in ["col3dt", "color3dt", "noise_col3dt", 
            "noise_color3dt"]:
            print('tbegin = '+str(params['tbegin']))
            print('tmax = '+str(params['tmax']))
            print('tdiff = '+str(params['tdiff']))
             
            if options.method[:6]=="noise_":
                noise=True
            plot_col_3dtimg(filename, params['I0'], params['lam'], 
                params['hue'], params['T'], params['tbegin'], params['tmax'], 
                params['tdiff'], params['fs'], MaxBox, params['dlmn'],
                params['nmax'], params['opt_axis'], params['fpns'], 
                add_n=params['add_n'], outfile=outname, noise=noise, 
                otype=options.type, mprocess=options.mprocess, pbc=pbc_lmn,
                poi=params['poi'], gauss=params['gauss'], 
                psf_type=params['psf_type'], tsO=params['tsO'])
    
        elif options.method == 'region':
            print("hue = "+str(params["hue"]))
            if options.calc=='show': #show specific
                plot_region(filename, params['I0'], params['lam'],
                             params['hue'], params['T'], options.timestep, 
                             params['fs'], MaxBox, Bm, params['scale'], 
                             params['dpi'])
    
            elif options.calc in ["specific", "spec"]: #Save Specific
                plot_region(filename, params['I0'], params['lam'],
                             params['hue'], params['T'], options.timestep, 
                             params['fs'], MaxBox, Bm, params['scale'], 
                             params['dpi'], outfile=outname)
            
            elif options.calc == 'all':
                print('Not implemented!')
 
        elif options.method == 'lumin':
            print("hue = "+str(params["hue"]))
            if options.calc=='show': #show specific
                plot_lumin(filename, params['I0'], params['lam'],
                             params['hue'], params['T'], options.timestep, 
                             params['fs'], MaxBox)
    
            elif options.calc in ["specific", "spec"]: #Save Specific
                print('Not implemented!')
            
            elif options.calc == 'all':
                print('Not implemented!')
        else:
            print('Method Not implemented') 

    elif remainder[0]=='video':
        print("fpns = "+str(params['fpns']))
        print("fourcc = "+str(params['fourcc']))
        if options.method!='data':
            for key in ["tbegin", "tmax", "tdiff", "T", "fs", "I0", "lam"]:
                print(key+" = "+str(params[key]))
            print('vid_ext = '+params['vid_ext'])

        if options.method == 'data':
            gen_vid_data(options.data,options.outname,
                         params['fpns'],params['fourcc'])       
        elif options.method in ["mono", "grey", "gray"]:
            gen_vid_mono(options.filename, params['tbegin'], params['tmax'],
                         params['tdiff'], params['T'], params['fs'], 
                         params['I0'], params['lam'], params['vid_ext'],
                         params['fpns'], params['fourcc'])
        elif options.method in ["color", "col"]:
            gen_vid_col(options.filename, params['tbegin'], params['tmax'], 
                        params['tdiff'], params['T'], params['fs'],
                        params['I0'], params['vid_ext'], params['fpns'], 
                        params['fourcc'])
        else:
            print('Method implemented')
 
    elif remainder[0]=='prop':
        if options.method == 'maxI':
            Imax=get_maxI(options.filename, params['lam'], params['fs'],
                tstep=options.timestep)
            print(Imax)
        elif options.method == 'predI0':
            if options.timestep!=None: #Using timestep as number of iterations
                get_I0s(options.filename, params['lam'],params['fs'],
                        iterations=options.iters, tstep=options.timestep)
            else:
                get_I0s(options.filename, params['lam'], params['fs'])

        elif options.method == 'hist':
            nor=False
            if options.calc=='norm':
                nor=True
            Imax=get_maxI(options.filename, params['lam'], params['fs'],
                tstep=options.timestep)
            get_hist(options.filename, params['lam'], params['fs'],
                options.outname, maxI=max(Imax), dI=0.1, norm=nor,
                tstep=options.timestep)
        elif options.method == 'num_area':
            if options.outname==None:
                write=False
            else:
                write=True
            get_num_area(options.filename, options.threshold, options.outname, 
                min_pix=params['min_pix'], col_channel=options.lambda_ID, 
                ncoor=params['focus_cor'], write=write) 
        elif options.method == 'num_vol':
            if options.outname==None:
                write=False
            else:
                write=True
            get_num_vol(options.filename, options.threshold, options.outname, 
                min_pix=params['min_pix'], col_channel=options.lambda_ID, 
                write=write)
        elif options.method == 'fcs':
            nidx=int(params['focus_cor']/params['dlmn'][2]+0.5)
            get_fcs(options.filename, options.outname, params['fcs_tmax'],
                    nidx=nidx) 
        elif options.method == 'fccs':
            nidx=int(params['focus_cor']/params['dlmn'][2]+0.5)
            get_fccs(options.filename, options.outname, params['fcs_tmax'], 
                     params['fccs_cols'][0], params['fccs_cols'][1], nidx=nidx)

        elif options.method == 'max_Nf':
            max_Nf(options.filename, params['tbegin'], params['tmax'],
                   params['tdiff'],options.mprocess)

        else:
            print('Method not implemented')
    
    elif remainder[0]=="convert":
        if options.method == "psf2tiff":
            if options.calc== None:
                dtype='uint8'
            else:
                dtype=options.calc
            psf_dat2tiff(options.filename, options.outname, params['Plmn'], 
                params['dlmn'], dtype, psf_type=params['psf_type'])
        if options.method == "psf2tiff2":
            if options.calc== None:
                dtype='uint8'
            else:
                dtype=options.calc
            psf_dat2tiff2(options.filename, options.outname, params['Plmn'], 
                params['dlmn'], dtype, psf_type=params['psf_type'])
        
        elif options.method == "nstack2tiff":
            nstack2tiff(options.data, options.outname, 
                params['dlmn'][2]*params['add_n']) 
            #Threshold -a is used as the spacing

        elif options.method == "tstack2tiff":
            tstack2tiff(options.data, params['fpns'], options.outname)

        elif options.method == "img2color":
            imgs2color(options.data,options.outname,params['mix_type'],
                params['hue'])
    else:
        print("Function not specified or implemented") 

if __name__=='__main__':
    main() 
