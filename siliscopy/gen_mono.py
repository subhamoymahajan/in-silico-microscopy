import os
import multiprocessing as mp
def gen_mono_c(data):
    os.system(os.path.dirname(__file__)+'/gen_mono -f '+data[0]+' -p '+data[1]+' -psf '+data[2]+' -o '+data[3])

def read_data(datafile):
    Arguments=[]
    f=open(datafile,'r')
    for lines in f: #Read data 
        foo=lines.split(',')
        Arguments.append(foo)
    return Arguments 

def gen_mono_c_mp(datafile):
    Arguments=read_data(datafile)
    pool=mp.Pool(mp.cpu_count())
    results=pool.map(gen_mono_c,Arguments)

def gen_mono_c_serial(datafile):
    Arguments=read_data(datafile)
    for i in range(len(Arguments)):
        gen_mono_c(Arguments[i])
