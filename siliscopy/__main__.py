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

from optparse import OptionParser
from siliscopy.gen_psf import *
from siliscopy.gen_mono import *
from siliscopy.plot_image import *
from siliscopy.create_vid import *
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
    
    
    options, remainder = parser.parse_args()
    #Read parameters
    if options.pfile !=None:
        params={}
        params['lam']=np.zeros(10,dtype=int)
        params['hue']=np.zeros(10)
        params['I0']=np.zeros(10)
        params['fourcc']='mp4v'
        params['fps']=1
        params['vid_ext']='.mov'
        f=open(options.pfile,'r')
        for lines in f:
            foo=lines.split('=')
            varname=foo[0].strip()
            val_string=foo[1].split('/')[0].strip()
            val_string=val_string.split()
            if varname in ["fs", "T", "dpi", "t0", "tmax", "tdiff", "opt_axis"]:
                params[varname]=int(val_string[0])
            elif varname in ["NA", "meu", "beta", "scale"]:
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
            elif varname=='vid_ext':
                params[varname]=val_string[0].strip()
            elif varname=='fourcc':
                foo=val_string[0].strip()
                params[varname]=foo[1:-1]
             
        #Reduce the lam variable
        for i in range(10):
            if params['lam'][i]<1:
                break
        params['lam']=params['lam'][:i]
        params['hue']=params['hue'][:i]
        params['I0']=params['I0'][:i]
    
    
    if remainder[0]=='gen_psf':
 
        beta=None
        if 'NA' in params and 'meu' in params:
            beta=np.arcsin(params['NA']/params['meu'])
        if 'beta' in params:
            beta=params['beta']
        if beta==None:
            raise Exception("Beta not provided!")
        
        print('beta = '+str(beta))
        for key in ["lam", "dlmn", "Plmn", "fs"]:
            print(key+' = '+str(params[key]))

        if options.method in ["gandy", "Gandy"]:#Currently the only method
            if options.calc=='all': #Calculates for all n' coordiantes
                if options.mprocess==True: #Multiprocessing. 
                    for lambd in params['lam']:
                        psf_gandy_mp(beta, lambd, params['dlmn'], 
                            params['Plmn'], params['fs'],options.outname)
                else: # Calculates Serially 
                    for lambd in params['lam']:
                        psf_gandy(beta, lambd, params['dlmn'], params['Plmn'], 
                                  params['fs'], options.outname)
            elif options.calc in ["specific", "spec"]: #Calculates for one n'
                for lambd in params['lam']:
                    psf_gandy_sep(beta, lambd, params['dlmn'], params['Plmn'], 
                        params['fs'],options.outname,'w',options.nid)
    
    elif remainder[0]=='gen_mono':
        if options.data!=None: #Data available in a file
            if options.mprocess==True: #Use multiprocessing
                # siliscopy gen_mono -d [data] -s
                gen_mono_c_mp(options.data)
            else:#Run serially
                # siliscopy gen_mono -d [data]
                gen_mono_c_serial(options.data)
        else: #No data file
            # siliscopy gen_mono -f [gro] -p [param] -o [out]
            gen_mono_c([options.filename, options.pfile, options.psfheader, 
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
        MaxBox=[0,0]
        MaxBox[0]=int(params['maxlen'][(params['opt_axis']+1)%3]/ \
                      params['dlmn'][0]+small)
        MaxBox[1]=int(params['maxlen'][(params['opt_axis']+2)%3]/ \
                      params['dlmn'][1]+small)
        Bm=params['maxlen'][(params['opt_axis']+1)%3] 

        if 'dpi' not in params:
            params['dpi']=600
            
        if options.method in ["mono", "grey", "gray"]:
            if options.calc=='show': #show specific
                plot_grey_img(options.filename, params['I0'], params['lam'],
                              params['T'], options.timestep, params['fs'],
                              MaxBox, Bm, params['scale'], params['dpi'])
    
            elif options.calc in ["specific", "spec"]: #Specific output
                plot_grey_img(options.filename, params['I0'], params['lam'],
                              params['T'], options.timestep, params['fs'],
                              MaxBox, Bm, params['scale'], params['dpi'], 
                              outname)
            
            elif options.calc == 'all':
                print('t0 = '+str(params['t0']))
                print('tmax = '+str(params['tmax']))
                print('tdiff = '+str(params['tdiff']))
                if options.mprocess==True: #parallel
                    plot_grey_mp(options.filename, params['I0'], params['lam'],
                                 params['T'], params['t0'], params['tmax'],
                                 params['tdiff'], params['fs'], MaxBox, Bm, 
                                 params['scale'], params['dpi'], outname)
                else: #Serial
                    plot_grey_serial(options.filename, params['I0'], 
                                     params['lam'], params['T'], params['t0'],
                                     params['tmax'], params['tdiff'], 
                                     params['fs'], MaxBox, Bm, params['scale'], 
                                     params['dpi'], outname)
        elif options.method in ["col", "color"]:
            print("hue = "+str(params["hue"]))
            if options.calc=='show': #show specific
                plot_col_img(options.filename, params['I0'], params['lam'],
                             params['hue'], params['T'], options.timestep, 
                             params['fs'], MaxBox, Bm, params['scale'], 
                             params['dpi'])
    
            elif options.calc in ["specific", "spec"]: #Save Specific
                plot_col_img(options.filename, params['I0'], params['lam'],
                             params['hue'], params['T'], options.timestep, 
                             params['fs'], MaxBox, Bm, params['scale'], 
                             params['dpi'], outfile=outname)
            
            elif options.calc == 'all':
                print('t0 = '+str(params['t0']))
                print('tmax = '+str(params['tmax']))
                print('tdiff = '+str(params['tdiff']))
                if options.mprocess==True: #parallel
                    plot_col_mp(options.filename, params['I0'], params['lam'],
                                params['hue'], params['T'], params['t0'], 
                                params['tmax'], params['tdiff'], params['fs'], 
                                MaxBox, Bm, params['scale'], params['dpi'], 
                                outname)
                else: #Serial
                    plot_col_serial(options.filename, params['I0'], 
                                    params['lam'], params['hue'], params['T'],
                                    params['t0'], params['tmax'],
                                    params['tdiff'], params['fs'], MaxBox, Bm,
                                    params['scale'], params['dpi'],outname) 

    elif remainder[0]=='video':
        print("fps = "+str(params['fps']))
        print("fourcc = "+str(params['fourcc']))
        if options.method!='data':
            for key in ["t0", "tmax", "tdiff", "T", "fs", "I0", "lam"]:
                print(key+" = "+str(params[key]))
            print('vid_ext = '+params['vid_ext'])

        if options.method == 'data':
            gen_vid_data(options.data,options.outname,
                         params['fps'],params['fourcc'])       
        elif options.method in ["mono", "grey", "gray"]:
            gen_vid_mono(options.filename, params['t0'], params['tmax'],
                         params['tdiff'], params['T'], params['fs'], 
                         params['I0'], params['lam'], params['vid_ext'],
                         params['fps'], params['fourcc'])
        elif options.method in ["color", "col"]:
            gen_vid_col(options.filename, params['t0'], params['tmax'], 
                        params['tdiff'], params['T'], params['fs'],
                        params['I0'], params['vid_ext'], params['fps'], 
                        params['fourcc'])
         
if __name__=='__main__':
    main() 
