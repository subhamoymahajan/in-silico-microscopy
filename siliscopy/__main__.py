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

def main():
    import numpy as np
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


    options, remainder = parser.parse_args()
    #Read parameters
    if options.pfile !=None:
        params={}
        #################DEFAULT VALUES ####################
        params['lam']=np.zeros(10,dtype=int)
        params['hue']=np.zeros(10)
        params['I0']=np.zeros(10)
        params['fourcc']='mp4v'
        params['fps']=1
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
        params['fs']=530
        params['psf_type']=0
        params['poi_a']=None
        params['gauss_b']=None
        params['frame_col']=1.0
        params['mix_type']='mt'
        #####################################################
        f=open(options.pfile,'r')
        for lines in f:
            foo=lines.split('=')
            varname=foo[0].strip()
            val_string=foo[1].split('/')[0].strip()
            val_string=val_string.split()
            
            if varname in ["fs", "T", "dpi", "tbegin", "tmax", "tdiff", "opt_axis",\
                "add_n", "psf_type", "nmax"]:
                params[varname]=int(val_string[0])
            elif varname in ["NA", "meu", "scale", "meu0", "t0", "meug", \
                "meug0", "tg", "tg0", "tsO", "meus", "poi", "gauss",\
                "frame_col"]:
                params[varname]=float(val_string[0])
            elif varname in ["dlmn", "Plmn","maxlen"]:
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
            elif varname in ["vid_ext", "mix_type"]:
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
        #Reduce the lam variable
        for i in range(10):
            if params['lam'][i]<1:
                break
        params['lam']=params['lam'][:i]
        params['hue']=params['hue'][:i]
        params['I0']=params['I0'][:i]
        MaxBox=[0,0]
        MaxBox[0]=int(params['maxlen'][(params['opt_axis']+1)%3]/ \
                      params['dlmn'][0]+small)
        MaxBox[1]=int(params['maxlen'][(params['opt_axis']+2)%3]/ \
                      params['dlmn'][1]+small)
    
    if remainder[0]=='gen_psf':
 
        print('NA = '+str(params['NA']))
        print('meu immersion oil = '+str(params['meu']))

        for key in ["lam", "dlmn", "Plmn"]:
            print(key+' = '+str(params[key]))
            if len(params[key])==0:
                raise Exception('Empty array')

        print('fs = '+str(params['fs']))

        if options.method in ["gandy", "Gandy"]:
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
            print("meu immersion oil (design) = "+str(params["meu0"])) 
            print("thickness immersion oil (design) = "+str(params["t0"])+" nm")
            print("meu coverslip = "+str(params["meug"])) 
            print("meu coverslip (design) = "+str(params["meug0"])) 
            print("thickness coverslip = "+str(params["tg"])+" nm")
            print("thickness coverslip (design) = "+str(params["tg0"])+" nm")
            print("distance between object focal plane and coverslip = "+str(params["tsO"])+" nm")
   
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
            print("meu immersion oil (design) = "+str(params["meu0"])) 
            print("thickness immersion oil (design) = "+str(params["t0"])+" nm")
            print("meu coverslip = "+str(params["meug"])) 
            print("meu coverslip (design) = "+str(params["meug0"])) 
            print("thickness coverslip = "+str(params["tg"])+" nm")
            print("thickness coverslip (design) = "+str(params["tg0"])+" nm")
            print("distance between object focal plane and coverslip = "+str(params["tsO"])+" nm")

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
    
    elif remainder[0]=='gen_mono':
        psf_header=options.psfheader
        tsO=None
        if options.method == None:
            options.method='slice'
        if params['psf_type']==1:
            psf_header+='_tsO'+str(params['tsO'])
            tsO=params['tsO']
        if options.data!=None: #Data available in a file
            if options.mprocess==True: #Use multiprocessing
                # siliscopy gen_mono -d [data] -s
                gen_mono_c_mp(options.data, params['maxlen'], 
                    params['opt_axis'], params['dlmn'], add_n=params['add_n'],
                    tsO=tsO)
            else:#Run serially
                # siliscopy gen_mono -d [data]
                gen_mono_c_serial(options.data, prams['maxlen'],
                    params['opt_axis'], parmas['dlmn'], add_n=params['add_n'],
                    tsO=tsO)
        elif options.method=='volume':
            gen_mono_c_vol([options.filename, options.pfile, 
                psf_header,options.outname], params['maxlen'],
                params['opt_axis'], params['dlmn'], 
                add_n=params['add_n'], mprocess=options.mprocess)
        elif options.method=='slice': #Generates one image slice.
            # siliscopy gen_mono -f [gro] -p [param] -o [out]
            gen_mono_c([options.filename, options.pfile, psf_header, 
                        options.outname])
    
    elif remainder[0]=='plot':
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
                plot_grey_img(options.filename, params['I0'], params['lam'],
                              params['T'], options.timestep, params['fs'],
                              MaxBox, Bm, params['scale'], params['dpi'], 
                              noise=noise, poi_a=params['poi_a'], 
                              gauss_b=params['gauss_b'])
    
            elif options.calc in ["specific", "spec"]: #Specific output
                plot_grey_img(options.filename, params['I0'], params['lam'],
                              params['T'], options.timestep, params['fs'],
                              MaxBox, Bm, params['scale'], params['dpi'], 
                              outname, noise=noise, poi_a=params['poi_a'],
                              gauss_b=params['gauss_b'], otype=options.type)
            
            elif options.calc == 'all':
                print('tbegin = '+str(params['tbegin']))
                print('tmax = '+str(params['tmax']))
                print('tdiff = '+str(params['tdiff']))
                if options.mprocess==True: #parallel
                    plot_grey_mp(options.filename, params['I0'], params['lam'],
                                 params['T'], params['tbegin'], params['tmax'],
                                 params['tdiff'], params['fs'], MaxBox, Bm, 
                                 params['scale'], params['dpi'], outname, 1.0,
                                 noise, params['poi_a'], params['gauss_b'], 
                                 options.type)
                else: #Serial
                    plot_grey_serial(options.filename, params['I0'], 
                                     params['lam'], params['T'], params['tbegin'],
                                     params['tmax'], params['tdiff'], 
                                     params['fs'], MaxBox, Bm, params['scale'], 
                                     params['dpi'], outname, 1.0, noise, 
                                     params['poi_a'], params['gauss_b'], 
                                     options.type)

        elif options.method in ["col", "color", "noise_col", "noise_color"]:
            print("hue = "+str(params["hue"]))
            if options.method[:6]=='noise_':
                noise=True
            if options.calc=='show': #show specific
                plot_col_img(options.filename, params['I0'], params['lam'],
                             params['hue'], params['T'], options.timestep, 
                             params['fs'], MaxBox, Bm, params['scale'], 
                             params['dpi'], frame_col=params['frame_col'],
                             mix_type=params['mix_type'], otype=options.type,
                             noise=noise, poi_a=params['poi_a'], 
                             gauss_b=params['gauss_b'])
    
            elif options.calc in ["specific", "spec"]: #Save Specific
                plot_col_img(options.filename, params['I0'], params['lam'],
                             params['hue'], params['T'], options.timestep, 
                             params['fs'], MaxBox, Bm, params['scale'], 
                             params['dpi'], outfile=outname, otype=options.type,
                             noise=noise, poi_a=params['poi_a'], 
                             gauss_b=params['gauss_b'])
            
            elif options.calc == 'all':
                print('tbegin = '+str(params['tbegin']))
                print('tmax = '+str(params['tmax']))
                print('tdiff = '+str(params['tdiff']))
                if options.mprocess==True: #parallel
                    plot_col_mp(options.filename, params['I0'], params['lam'],
                                params['hue'], params['T'], params['tbegin'], 
                                params['tmax'], params['tdiff'], params['fs'], 
                                MaxBox, Bm, params['scale'], params['dpi'], 
                                outname, params['frame_col'], 
                                params['mix_type'],noise, params['poi_a'],
                                params['gauss_b'], options.type)

                else: #Serial
                    plot_col_serial(options.filename, params['I0'], 
                                    params['lam'], params['hue'], params['T'],
                                    params['tbegin'], params['tmax'],
                                    params['tdiff'], params['fs'], MaxBox, Bm,
                                    params['scale'], params['dpi'],
                                    outname, params['frame_col'],
                                    params['mix_type'], noise, params['poi_a'],
                                    params['gauss_b'], options.type)

        elif options.method in ["mono3d", "grey3d", "gray3d", "noise_mono3d", \
            "noise_grey3d", "noise_gray3d"]:
             
            if options.method[:6]=="noise_":
                noise=True
            plot_grey_3dimg(options.filename, params['I0'], params['lam'],
                            params['T'], options.timestep, params['fs'], MaxBox,
                            params['dlmn'], params['nmax'], params['opt_axis'],
                            add_n=params['add_n'], outfile=outname, noise=noise,
                            frame_col=1.0, otype=options.type, 
                            mprocess=options.mprocess, poi_a=params['poi_a'],
                            gauss_b=params['gauss_b'])

        elif options.method in ["col3d", "colo3d", "noise_col3d", 
            "noise_color3d"]:
             
            if options.method[:6]=="noise_":
                noise=True
            plot_col_3dimg(options.filename, params['I0'], params['lam'],
                           params['hue'], params['T'], options.timestep, 
                           params['fs'], MaxBox, params['dlmn'], params['nmax'], 
                           params['opt_axis'], add_n=params['add_n'], 
                           outfile=outname, noise=noise, frame_col=1.0, 
                           otype=options.type, mprocess=options.mprocess, 
                           poi_a=params['poi_a'], gauss_b=params['gauss_b'])

        elif options.method in ["mono3dt", "grey3dt", "gray3dt", "noise_mono3dt",\
            "noise_grey3dt", "noise_gray3dt"]:
             
            if options.method[:6]=="noise_":
                noise=True
            plot_grey_3dtimg(options.filename, params['I0'], params['lam'],
                             params['T'], params['tbegin'], params['tmax'], 
                             params['tdiff'], params['fs'], MaxBox, 
                             params['dlmn'], params['nmax'], params['opt_axis'],
                             params['fps'], add_n=params['add_n'], noise=noise,
                             outfile=outname, frame_col=1.0, otype=options.type,  
                             mprocess=options.mprocess, poi_a=params['poi_a'],
                             gauss_b=params['gauss_b'])

        elif options.method in ["col3dt", "colo3dt", "noise_col3dt", 
            "noise_color3dt"]:
             
            if options.method[:6]=="noise_":
                noise=True
            plot_col_3dimg(options.filename, params['I0'], params['lam'],
                           params['hue'], params['T'], options.timestep, 
                           params['fs'], MaxBox, params['dlmn'], params['nmax'], 
                           params['opt_axis'], params['fps'], frame_col=1.0,
                           add_n=params['add_n'], outfile=outname, noise=noise, 
                           otype=options.type, mprocess=options.mprocess, 
                           poi_a=params['poi_a'], gauss_b=params['gauss_b'])
    
        elif options.method == 'region':
            print("hue = "+str(params["hue"]))
            if options.calc=='show': #show specific
                plot_region(options.filename, params['I0'], params['lam'],
                             params['hue'], params['T'], options.timestep, 
                             params['fs'], MaxBox, Bm, params['scale'], 
                             params['dpi'])
    
            elif options.calc in ["specific", "spec"]: #Save Specific
                plot_region(options.filename, params['I0'], params['lam'],
                             params['hue'], params['T'], options.timestep, 
                             params['fs'], MaxBox, Bm, params['scale'], 
                             params['dpi'], outfile=outname)
            
            elif options.calc == 'all':
                print('Not implemented!')
 
        elif options.method == 'lumin':
            print("hue = "+str(params["hue"]))
            if options.calc=='show': #show specific
                plot_lumin(options.filename, params['I0'], params['lam'],
                             params['hue'], params['T'], options.timestep, 
                             params['fs'], MaxBox)
    
            elif options.calc in ["specific", "spec"]: #Save Specific
                print('Not implemented!')
            
            elif options.calc == 'all':
                print('Not implemented!')
        else:
            print('Method Not implemented') 

    elif remainder[0]=='video':
        print("fps = "+str(params['fps']))
        print("fourcc = "+str(params['fourcc']))
        if options.method!='data':
            for key in ["tbegin", "tmax", "tdiff", "T", "fs", "I0", "lam"]:
                print(key+" = "+str(params[key]))
            print('vid_ext = '+params['vid_ext'])

        if options.method == 'data':
            gen_vid_data(options.data,options.outname,
                         params['fps'],params['fourcc'])       
        elif options.method in ["mono", "grey", "gray"]:
            gen_vid_mono(options.filename, params['tbegin'], params['tmax'],
                         params['tdiff'], params['T'], params['fs'], 
                         params['I0'], params['lam'], params['vid_ext'],
                         params['fps'], params['fourcc'])
        elif options.method in ["color", "col"]:
            gen_vid_col(options.filename, params['tbegin'], params['tmax'], 
                        params['tdiff'], params['T'], params['fs'],
                        params['I0'], params['vid_ext'], params['fps'], 
                        params['fourcc'])
        else:
            print('Method implemented') 
    elif remainder[0]=='prop':
        if options.method == 'maxI':
            Imax=get_maxI(options.filename, params['lam'], params['fs'])
            print(Imax)
        elif options.method == 'predI0':
            if options.timestep!=None: #Using timestep as number of iterations
                get_I0s(options.filename, params['lam'],params['fs'],
                        iterations=options.timestep)
            else:
                get_I0s(options.filename, params['lam'], params['fs'])

        elif options.method == 'hist':
            nor=False
            if options.calc=='norm':
                nor=True
            Imax=get_maxI(options.filename, params['lam'], params['fs'])
            get_hist(options.filename, params['lam'], params['fs'],
                     options.outname, maxI=max(Imax), dI=0.1, norm=nor)
        elif options.method == 'num_area':
            foo=params['pbc']
            pbc=[foo[(params['opt_axis']+1)%3], foo[(params['opt_axis']+2)%3]]
            sh=False
            white =True
            if options.calc=='show':
                sh = True
            elif options.calc=='show-test':
                sh = True
                white = False
            elif options.calc=='test':
                white = False
            if options.threshold<1:
                thres=int(options.threshold*255)
            else:
                thres=int(options.threshold)
            get_num_area(options.filename, params['I0'][options.lambda_ID], \
                          params['lam'][options.lambda_ID], params['T'], \
                          options.timestep, params['fs'], \
                          thres, MaxBox, params['dlmn'], \
                          pbc, options.outname, sh, white ) 
        else:
            print('Method not implemented')
    elif remainder[0]=="convert":
        if options.method == "psf2tiff":
            if options.calc== None:
                dtype='uint16'
            else:
                dtype=options.calc
            psf_dat2tiff(options.filename, options.outname, params['Plmn'], \
                params['dlmn'], dtype)
        
        elif options.method == "zstack2tiff":
            if options.calc== None:
                dtype='uint16'
            else:
                dtype=options.calc
            zmax=int(params['maxlen'][params['opt_axis']]/params['dlmn'][2]+0.5)
            zstack_dat2tiff(options.filename, options.outname, params['lam'], \
                params['I0'], params['fs'], options.timestep, params['T'], \
                MaxBox, zmax, params['opt_axis'], dtype, params['add_n']) 

        elif options.method == "tstack2tiff":
            if options.calc== None:
                dtype='uint16'
            else:
                dtype=options.calc
            tstack_dat2tiff(options.filename, options.outname, params['lam'], \
                params['I0'], params['fs'], params['tbegin'], params['tmax'], \
                params['tdiff'], params['T'], MaxBox, dtype) 
        elif options.method == "img2color":
            print('not implemented')  
    else:
        print("Function not specified or implemented") 
if __name__=='__main__':
    main() 
