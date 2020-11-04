#!/usr/bin/env python3.5
import glob, sys, string, time
import os, fnmatch, re
import shutil, getpass 
import stat, platform, math
import multiprocessing
import matplotlib
import optparse
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import numpy.random as rnd
from matplotlib import cm as CM
import subprocess as subp

## Create a general dictionary of script parameters. Easy to expand later
author = 'prabu@biochem.mpg.de'
quick2d_parameters = {}
quick2d_parameters['datadir'] = os.getcwd()
quick2d_parameters['user'] = getpass.getuser()
quick2d_parameters['user_uid'] = os.getuid()
quick2d_parameters['micrograph_name_pattern'] = '*'
quick2d_parameters['micrograph_name_suffix'] = 'mrcs'
quick2d_parameters['micrograph_name_exclude'] = '\#'
quick2d_parameters['pixel_size'] = 0
quick2d_parameters['z'] = 0
quick2d_parameters['kv'] = 300
quick2d_parameters['cs'] = 2.7
quick2d_parameters['diameter'] = 0
quick2d_parameters['box'] = 0
quick2d_parameters['cccutoff'] = 0.2
quick2d_parameters['mask'] = 0
quick2d_parameters['dose'] = 0
quick2d_parameters['Tval'] = 2
quick2d_parameters['nclasses'] = 15
quick2d_parameters['ignorectf'] = ''
quick2d_parameters['negative'] = 0
quick2d_parameters['contrast'] = ' --invert_contrast '
quick2d_parameters['partbin'] = ' --scale 96 '
quick2d_parameters['pickcontrast'] = '  '
quick2d_parameters['preexposure'] = 0
quick2d_parameters['time'] =  time.strftime("%d_%b_%Y_%I_%M_%S%p")
quick2d_parameters['timestamp'] =  'QUICK2D_' + quick2d_parameters['time']
quick2d_parameters['workdir'] = os.getcwd() + '/' + quick2d_parameters['timestamp'] + '/'
quick2d_parameters['power_users'] = 0
quick2d_parameters['write_movies'] = 'NO'
quick2d_parameters['ngpu'] = 0
quick2d_parameters['bin'] = 1
quick2d_parameters['ncpu'] = 0
quick2d_parameters['gain'] = ' '
quick2d_parameters['input_commands'] = ' '.join(sys.argv)
quick2d_parameters['rescut'] = 0
quick2d_parameters['dose_filter'] = 'NO' 
quick2d_parameters['template'] = '' 
quick2d_parameters['templatepix'] = '' 
quick2d_parameters['3dmodel'] = '' 
quick2d_parameters['3dmodelpix'] = '' 
quick2d_parameters['lavgmin'] = ' --lave_min -1 '
quick2d_parameters['lavgmax'] = ' --lave_max 1.2 '



### Command line and help text
### Argparse
from optparse import OptionParser

usage = "usage %prog "+ """ [options] 

        Examples:
        
        Minimal input: pixel size and diameter
        quick2d.com --pix=1.99 --dia=200
    
        Different kv and cs
        quick2d.com --pix=1.99 --dia=200 --kv=200 --cs=2.0

        For negative stain data ( default is cryo data )
        quick2d.com --pix=1.99 --dia=200 --negative

        By defaults particles will be binned to 96 px.. If you want to to keep original particle size  
        quick2d.com --pix=1.99 --dia=200 --nopartbin

        Give custom box size ( in pixels ) and mask size ( A )  
        quick2d.com --pix=1.99 --dia=200 --box=200 --mask=250

        Give template and its pixels size for autopicking 
        quick2d.com --pix=1.99 --dia=200 --template=class.mrcs --Tpix=2.37

        Give 3D Model and its pixels size, projections are made and used for autopicking 
        quick2d.com --pix=1.99 --dia=200 --3dmodel=run_it025.mrc --3dmodelpix=2.37

        For the dose wighting 
        quick2d.com --pix=1.99 --dia=200 --dose=22

Read below for much finer control of individual modules
        """


parser = OptionParser(usage=usage)
must = optparse.OptionGroup(parser, 'Mandatory Options')
must.add_option("--pix",               dest="pixel_size",        help="Pixel size of the Micrograph eg --pix=1.42 No default")
must.add_option("--dia",               dest="diameter",          help="Particle dimeter in A eg --dia=200 No default")
parser.add_option_group(must)


gen = optparse.OptionGroup(parser, 'General Options')
gen.add_option("--negative",        dest="negative",          action='store_true', default=False,  help="For negative stain data Default  false")
gen.add_option("--nopartbin",       dest="partbin",           action='store_false', default=True,  help="by default all particles will be binned to 96 px (for speed), use this option if you dont want to bin it")
parser.add_option_group(gen)

mic = optparse.OptionGroup(parser, 'Microscopy')
mic.add_option("--kv",                 dest="kv",                help="Voltage of the Microscope  eg  --kv=300  default 300 ")
mic.add_option("--cs",                 dest="cs",                help="Spherical aberation of the Microscope  eg  --cs=2.7  default  2.7  ")
parser.add_option_group(mic)

select = optparse.OptionGroup(parser, 'Micrograph Selection')
select.add_option("--mic",             dest="pattern",           help="Common unique pattern in the micrograph name ( to identify/select)  eg  --m=_frames default all mrc files in the directory")
select.add_option("--exc",             dest="exclude_pattern",   help="Micrograph name that you want to exclude    eg  --exc=_DW No default ")
select.add_option("--suf",             dest="suffix",            help="Micrograph suffix  eg --suf=mrc  default mrcs or mrc ")
parser.add_option_group(select)

movie = optparse.OptionGroup(parser, 'Movies')
movie.add_option("--dose",             dest="dose",              help="Exposure per frame e/A^2 --dose=22 No default ")
movie.add_option("--gain",             dest="gain",              help="Camera gain file for correction eg --gain=K2gain.em  No default")
movie.add_option("--bin",              dest="bin",               help="Binning factor for micrographs default - no binning")
movie.add_option("--stack",            dest="stack",             action="store_true", default=False,  help="Write alinged movie stacks eg --stack Default false")
movie.add_option("--preexp",           dest="preexposure",       help=" Pre Exposure per frame e/A^2 No default")
parser.add_option_group(movie)

picking = optparse.OptionGroup(parser, 'Particle picking')
picking.add_option('--box',            dest="box",               help="Box Size (in pixels) for windowing the particles eg --box=128 default - 1.8 x particle diameter")
picking.add_option('--cc',             dest="cccutoff",          help="cc cutoff to use for particle picking eg --cc=0.3 Default  0.2") 
picking.add_option('--lave_min ',      dest="lavgmin",           help="autopicking option for ice detection  eg --lave_min  Default  -1.0") 
picking.add_option('--lave_max ',      dest="lavgmax",           help="autopicking option for ice detection  eg --lave_max  Default   1.2") 
picking.add_option("--template",       dest="template",          help="template for particle picking Default  false")
picking.add_option("--tpix",           dest="templatepix",       help="template pixel size")
picking.add_option("--3dmodel",        dest="model",             help="Projections will be made using this 3D model and used for particle picking Default  false")
picking.add_option("--3dmodelpix",     dest="modelpix",          help="Pixel size of the 3D model  Default nothing")
parser.add_option_group(picking)

class2d = optparse.OptionGroup(parser, '2D Classification')
class2d.add_option("--mask",           dest="mask",              help="Mask diameter eg --mask=200 Default  - 1.1 x particle diameter")
class2d.add_option("--class",          dest="nclasses",          help="Number of 2D classes to calculate eg --class=30 Default  15")
class2d.add_option("--Tval",           dest="Tval",              help="T regularization parameter eg --Tval=1 Default  1")
class2d.add_option("--ignorectf",      dest="ignorectf",         action='store_true', default=False,  help="Ignore ctf till first peak eg --ignorectf Default  false")
parser.add_option_group(class2d)

computing = optparse.OptionGroup(parser, 'Computing')
computing.add_option("--cpu",          dest="ncpu",              help="Number of cpus to use eg --cpu=20 Default use all available")
computing.add_option("--gpu",          dest="ngpu",              help="Number of gpus to use eg --gpu=2 defaults use all available")
computing.add_option("--devel",        dest="power",             action="store_true", default=False, help="For developers....")
computing.add_option("--rcut",         dest="rescut",            help="only works with --devel")
parser.add_option_group(computing)

(options, args) = parser.parse_args()

if options.kv :
    quick2d_parameters['kv'] = options.kv

if options.pixel_size :
    quick2d_parameters['pixel_size'] = options.pixel_size

if options.dose :
    quick2d_parameters['dose'] = options.dose

if options.cccutoff :
    quick2d_parameters['cccutoff'] = options.cccutoff

if options.box :
    quick2d_parameters['box'] = options.box

if options.Tval :
    quick2d_parameters['Tval'] = options.Tval

if options.template :
    quick2d_parameters['template'] = options.template

if options.templatepix :
    quick2d_parameters['templatepix'] = options.templatepix

if options.model :
    quick2d_parameters['3dmodel'] = options.model

if options.nclasses :
    quick2d_parameters['nclasses'] = options.nclasses

if options.lavgmin :
    quick2d_parameters['lavgmin'] = ' --lave_min ' + options.lavgmin

if options.lavgmax :
    quick2d_parameters['lavgmax'] = ' --lave_max ' + options.lavgmax

if options.modelpix :
    quick2d_parameters['3dmodelpix'] = options.modelpix

if options.mask :
    quick2d_parameters['mask'] = options.mask

if options.diameter :
    quick2d_parameters['diameter'] = options.diameter

if options.cs :
    quick2d_parameters['cs'] = options.cs

if options.preexposure :
    quick2d_parameters['preexposure'] = options.preexposure

if options.pattern :
    quick2d_parameters['micrograph_name_pattern'] = options.pattern

if options.exclude_pattern :
    quick2d_parameters['micrograph_name_exclude'] = options.exclude_pattern

if options.gain :
    quick2d_parameters['gain'] = options.gain

if options.bin :
    quick2d_parameters['bin'] = options.bin

if options.suffix :
    quick2d_parameters['micrograph_name_suffix'] = options.suffix

if options.ncpu :
    quick2d_parameters['ncpu'] = options.ncpu

if options.ngpu :
    quick2d_parameters['ngpu'] = options.ngpu

if options.rescut :
    quick2d_parameters['rescut'] = options.rescut

if str(options.stack) == "True" :
    quick2d_parameters['write_movies'] = 'YES'

if str(options.ignorectf) == "True" :
    quick2d_parameters['ignorectf'] = ' --ctf_intact_first_peak '

if str(options.partbin) == "False" :
    quick2d_parameters['partbin'] = '  '

if str(options.negative) == "True" :
    quick2d_parameters['contrast'] = ' '
    quick2d_parameters['negative'] = 1
    quick2d_parameters['pickcontrast'] = ' --dont_invertT '

if str(options.power) == "True" :
    quick2d_parameters['power_users'] = 1

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)


if quick2d_parameters['pixel_size'] == 0 :
    print ( "\nI need pixel size as a input, I kind of assume defaults for other values" )
    print ( "please provide pixel size For eg. --pix=1.43 or -p 1.43\n" )
    print ( "If you need help please use quick2d.com -h" )
    quit ()

if quick2d_parameters['template'] != '' :
    file = quick2d_parameters['datadir'] + '/' + quick2d_parameters['template']
    if os.path.exists (file) != True :
        print ("\nI cant fine the template file for particle picking")
        quit ()
    elif quick2d_parameters['templatepix'] == '' :
        print ("\n You have given template but not the pixel size" )
        print ("\n For eg --tpix=2.54 ")
        quit ()

if quick2d_parameters['3dmodel'] != '' :
    file = quick2d_parameters['datadir'] + '/' + quick2d_parameters['3dmodel']
    if os.path.exists (file) != True :
        print ("\nI cant fine the template file for particle picking")
        quit ()
    elif quick2d_parameters['3dmodelpix'] == '' :
        print ("\n You have given 3dmodel but not the pixel size" )
        print ("\n For eg --3dmodelpix=2.54 ")
        quit ()

if quick2d_parameters['box'] != 0 :
    quick2d_parameters['partbin'] = ''
    print ("\n You have given box size, so no particles binning" )
elif quick2d_parameters['mask'] != 0 :
    quick2d_parameters['partbin'] = ''
    print ("\n You have given mask size, so no particles binning" )


if quick2d_parameters['diameter'] == 0 :
    print ( "\nI need particle diameter as a input " )
    print ( "please provide particle dimaeter For eg. --dia=128 or -d 128\n" )
    print ( "If you need help please use quick2d.com -h" )
    quit ()

if quick2d_parameters['dose'] != 0 :
    print ( "\nI will do dose weighting" )
    quick2d_parameters['dose_filter'] = 'TRUE' 
    if quick2d_parameters['kv'] == 0 :
        print ( "please provide KV of the microscope For eg. --kv=200 or -k 200\n" )
        print ( "If you need help please use quick2d.com -h" )

#if quick2d_parameters['gain'] == '' :
#    print ( "\nI need camera gain, I kind of assume defaults for other values" )
#    print ( "please provide the gain For eg. --g=Gain.em or -g Gain.em\n" )
#    print ( "If you need help please use quick2d.com -h" )
    #quit ()

## CHECK FOR THE PRESENSE OF EXECUTABLES
def check_execs() :

    exec_list = ['header', 'Gctf', 'relion_refine_mpi', 'relion_preprocess_mpi', 'MotionCor2', 'Gautomatch','e2proc2d.py','relion_display']

    for i in range (len(exec_list)):
        if not check_for_executables(exec_list[i]) :
            quick2d_parameters.setdefault('missing_execs', []).append(exec_list[i])

    if 'missing_execs' in quick2d_parameters :
        print ('\n\nThe following exeucutables are missing', quick2d_parameters['missing_execs'] )
        print ('I cant work without them.. Please install them and retry it')
        print ('I will quit\n')
        quit ()

def check_for_executables(program):
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file


    return None

check_execs()

class bcolors:
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

sg = bcolors.OKGREEN
eb = bcolors.ENDC


def count_micrographs(pattern, micdir, exclude, printflag) :
    micrographs_list = glob.glob(os.path.join(micdir,pattern ))
 
    for name in glob.glob(os.path.join(micdir,pattern )) :
        if exclude in name :
            micrographs_list.remove(name)
    micrographs_total = len( micrographs_list )

    if micrographs_total != 0 :
        print ('\nYou have',micrographs_total, 'micrographs')
        if micrographs_total > 25 :
            micrographs_list = micrographs_list[0:25] 
            micrographs_total = len( micrographs_list )
            print ('\nI have taken only',micrographs_total, 'micrographs ( quick 2d )')
        
    return micrographs_total,micrographs_list

def number_of_gpu( micrographs_list ) :
    gpus = len(glob.glob('/proc/driver/nvidia/gpus/*'))
    if gpus == 0 :
        print ('I cant find any usable GPU')
        print ('I will quit')
        quit ()
     
    if len (micrographs_list) < gpus :
        quick2d_parameters['ngpu']  = len (micrographs_list)

    quick2d_parameters['ngpu']  = int( quick2d_parameters['ngpu'] )
    if quick2d_parameters['ngpu'] == 0 :
        print (' Available GPUS: ',gpus )
    else :
        if ( gpus != 0 ) :
            print (' Available GPUS: ',gpus )
            print (' Using ', quick2d_parameters['ngpu'], 'GPUS' )
            gpus = quick2d_parameters['ngpu']
    return gpus

def number_of_cpu() :
    cpus = multiprocessing.cpu_count()
    quick2d_parameters['ncpu']  = int (quick2d_parameters['ncpu'] )
    if quick2d_parameters['ncpu'] == 0 :
        print (' Available CPUS: ',cpus )
    else :
        if ( cpus != 0 ) :
            cpus = quick2d_parameters['ncpu']
            print (' Using ', quick2d_parameters['ncpu'], 'CPUS' )
    return cpus

def get_movie_frame_number ( micrographs_list ) :
    args = 'header -size -i ' + micrographs_list[0] 
    proc = subp.Popen( args,  shell = True,  stdout=subp.PIPE ) 
    proc.wait
    z  =  int ( proc.communicate()[0].decode("utf-8").split()[2] )
    return z

def grep(pattern,fileObj):
    r=[]
    linenumber=0
    for line in fileObj:
        linenumber +=1
        if re.search(pattern,line):
            r.append((line))
    return r

def round_up_to_even(f):
    return math.ceil(f / 2.) * 2

def make_dir_soft_links ( workdir, micrographs_list) :
    dest_dir  = workdir + '/'
    mic_dir  = workdir + '/micrographs'
    os.makedirs(workdir)
    os.makedirs(mic_dir)
    #for file in glob.glob (pattern) :
     #   if exclude not in file :
      #      file_with_dir = micdir + '/' + file
       #     file_with_mic_dir = mic_dir + '/' + file
        #    os.symlink ( file_with_dir, file_with_mic_dir )
    for file in micrographs_list :
        micname =  file.split('/')[-1]
        file_with_mic_dir = mic_dir + '/' + micname
        os.symlink ( file, file_with_mic_dir )

    os.chdir ( mic_dir )

def make_dir_soft_links1 ( workdir, micrographs_list ) :
    dest_dir  = workdir + '/'
    movie_dir  = workdir + '/movies'
    os.makedirs(workdir)
    os.makedirs(movie_dir)
    number_of_gpus = quick2d_parameters['ngpu']
    split_by_gpu = np.array_split(micrographs_list,number_of_gpus)
    i = 0
    while ( i < number_of_gpus ) :
        gpu_dir = movie_dir + '/gpu' + str(i)
        os.makedirs(gpu_dir)
        for file in split_by_gpu[i] :
            micname =  file.split('/')[-1]
            file_with_gpu_dir = gpu_dir + '/' + micname
            os.symlink ( file, file_with_gpu_dir )
        i += 1

def write_unblur_input ( micname, frames, pixel_size , dose_filter, exposure, kv, preexposure, movies ) :
    name = ''.join ( micname.split('.')[:-1] )
    suf = micname.split('.')[-1]
    real = '../../' + micname
    dest = name + '_in.mrc' 
    os.symlink ( real, dest )
    unblur_input = []
    unblur_input.append ( '#!/usr/bin/csh -f' )
    unblur_input.append ('unblur  > ' +  name + '.unblur.log << EOF' )
    unblur_input.append ( name + '_in.mrc' )
    unblur_input.append ( frames )
    unblur_input.append ( '../unblur_sum/' + name + '.mrc' )
    unblur_input.append ( name + '_shifts.txt' )
    unblur_input.append ( pixel_size )
    if dose_filter == 'YES' :
        unblur_input.append ( 'YES' )
        unblur_input.append ( exposure )
        unblur_input.append ( kv )
        unblur_input.append ( preexposure )
    else :
        unblur_input.append ( 'NO' )
    if movies == 'YES' :
        unblur_input.append ( movies )
        unblur_input.append ( '../unblur_movie/' + name + '_movie.mrc' )
    else :
        unblur_input.append ( movies )
    unblur_input.append ( 'NO' )
    unblur_input.append ( 'EOF' )
    i = 0
    out = name + '.com'
    com = open (out, "w")
    while (i < len (unblur_input) ):
        com.write( str (unblur_input[i] ) ) 
        com.write("\n"  ) 
        i += 1
    st = os.stat(out)
    os.chmod(out, st.st_mode | stat.S_IEXEC)
    com.close()

def unblur_align ( micrographs_list) :
    number_of_cpus = number_of_cpu()
    os.makedirs('unblur_log')
    os.makedirs('unblur_sum')
    os.makedirs('unblur_movie')
    os.makedirs('unblur_ctf')
    os.chdir   ('unblur_log')

    the_dot_with_suff = '.' + quick2d_parameters['micrograph_name_suffix']
    for mics in micrographs_list :
       micname =     mics.split('/')[-1]
       write_unblur_input( micname, quick2d_parameters['z'], quick2d_parameters['pixel_size'], quick2d_parameters['dose_filter'], \
                  quick2d_parameters['dose'], quick2d_parameters['kv'], quick2d_parameters['preexposure'], quick2d_parameters['write_movies']  )
    com_list = glob.glob("*.com" )
    if  len (micrographs_list) < number_of_cpus :
        number_of_cpus = len (micrographs_list)
    split_by_cpu = np.array_split(com_list,number_of_cpus)
    for i in range (number_of_cpus) :
        out = 'cpu'+ str (i)
        outfile = open ( out, 'w')
        for j in split_by_cpu[i] :
            outfile.write( j ) 
            outfile.write("\n"  ) 
        outfile.close()
        st = os.stat(out)
        os.chmod(out, st.st_mode | stat.S_IEXEC)

    ps = []
    output = open("unblur.log", 'w')
    for i in range(number_of_cpus) :
        args = 'cpu'+ str (i)
        proc = subp.Popen( args,  shell = True,  stdout = output)
        time.sleep(0.1)
        ps.append(proc)

    #Keep the user informed
    #shifts.txt = len(glob.glob("*_shifts.txt" ))
    b = len (micrographs_list )
    c = 0
    a = 0 
    print ('\n Aligning using Unblur \n')
    while ( a < b ) :
        if ( a > c ) :
            print ( a,'/',b, 'dat file are converted to mrc(s)' )
            c = a
        a = len(glob.glob("*_shifts.txt" ))

    for p in ps:
        p.wait()

def MotionCor2_align ( number_of_gpus ) :
    output_dir = quick2d_parameters['workdir'] + '/aligned/'
    os.makedirs(output_dir)
    os.chdir (output_dir)
    if quick2d_parameters['gain'] != "":
        gain = ' -Gain ' + quick2d_parameters['gain'] + ' '
    else :
        gain = ""

    if quick2d_parameters['micrograph_name_suffix'] == "tif":
        input_stem = ' -InTiff '
    else :
        input_stem = ' -InMrc '

    output = open("MotionCor2.log", 'w')
    ps = []
    for i in range(number_of_gpus) :
        input_dir  = quick2d_parameters['workdir'] + '/movies/gpu' + str (i) + '/'
        outlog = 'gpu' + str (i) + '.log'
        output = open(outlog, 'w')
        args = 'MotionCor2 ' + input_stem  + input_dir  + ' -OutMrc ' + output_dir + ' -patch 5 5  -Tol 0.5 -Serial 1  -Gpu ' + str (i )  \
                 + ' -FtBin ' + str(quick2d_parameters['bin']) + gain + ' -InitDose ' + str (quick2d_parameters['preexposure'])  + ' -FmDose ' + str(quick2d_parameters['dose']) \
                 + ' -PixSize '+ str (quick2d_parameters['pixel_size']) + ' -kV ' + str (quick2d_parameters['kv']) + ' -Throw 1 '
        output.writelines ( str(args ) + '\n' )
        proc = subp.Popen( args,  shell = True,  stdout = output)
        #print (args)
        ps.append(proc)

    #Keep the user informed       
    if quick2d_parameters['dose'] == 0 :
        #list = output_dir + '/*' + quick2d_parameters['micrograph_name_suffix']
        list = output_dir + '/*.mrc'
    else :
        #list = output_dir + '/*_DW.' + quick2d_parameters['micrograph_name_suffix']
        list = output_dir + '/*_DW.mrc'
    a = len(glob.glob(list))
    b = len (micrographs_list )
    c = 0
    print ('\nAligning the frames using MotionCor2\n')
    while ( a < b ) :
        if ( a > c ) :
            #pbar = ProgressBar(maxval=b).start()
            #time.sleep(0.01)
            #pbar.update(a+1)
            print ( a,'/',b, 'Micrographs are done' )
            c = a
        a = len(glob.glob(list ))
    #pbar.finish()
    for p in ps:
        p.wait()

    if quick2d_parameters ['micrograph_name_suffix'] == 'mrcs' :
        mrcs_dir_name = output_dir + '/*mrcs'
        #names_list = [w.replace(suf, '') for w in micrographs_list]
        for mrcs in glob.glob(mrcs_dir_name) :
            name = mrcs.replace('mrcs','')
            mrc = name + 'mrc'
            os.rename ( mrcs, mrc)
    print (sg + '\nMovie alignment done')
    print (eb)

def GCTF ( pixel_size, cs, kv, micrographs_list, negative ) :
    #output_dir = quick2d_parameters['workdir'] + '/ctf_Gctf/'
    #os.makedirs(output_dir)
    #os.chdir (output_dir) 
    number_of_gpus = quick2d_parameters['ngpu']
    split_by_gpu = np.array_split(micrographs_list,number_of_gpus)
    import subprocess as subp
    import time
    if float (pixel_size)  > 3 :
        resH = 2 * float (pixel_size)
        B_resH = 9.00
    else :
        resH = 4  
        B_resH = 7.00
    if ( negative == 1 ) :
        box = ' '
        ac = ' --ac 0.35 '
    else :
        box = ' --boxsize 512 ' 
        ac = ' --ac 0.07 '
    output = open("gctf.log", 'w')
    ps = []
    for i in range(number_of_gpus) :
        star = 'gpu' + str (i) + '.star'
        args = 'Gctf --apix '+ str(pixel_size)+ ' --kV ' + str(kv) +' --cs ' + str(cs) + ' ' + ' '.join(split_by_gpu[i]) \
               + ' --gid ' + str(i) + ' --do_validation' + ' --resH ' + str(resH) + ' --B_resH ' + str(B_resH) \
               + box + ac + ' --ctfstar ' + star
               #+ ' --ac 0.07 --do_EPA  --boxsize 512 --do_Hres_ref --Href_resL 20'  
        #print (args)
        proc = subp.Popen( args,  shell = True,  stdout = output)
        ps.append(proc)

    #Keep the user informed       
    a = len(glob.glob("*gctf.log" ))
    b = len (micrographs_list )
    c = 0
    print ('\nCalculate CTF for your micrographs using Gctf\n')
    while ( a < b ) :
        if ( a > c ) :
            print ( a,'/',b, 'Micrographs are done' )
            c = a
        a = len(glob.glob("*gctf.log" ))

    for p in ps:
        p.wait()
         
    star = 'micrographs.star'
    g = open (star, "w")
    g.writelines ( '\n' + 'data_' + '\n\n' + 'loop_' + '\n' + '_rlnMicrographName #1' + '\n' )
    for mic in micrographs_list :
        mic = mic.split('/')[-1] 
        g.writelines (  mic + '\n' )
    g.close()
    

    output = open("star.log", 'w')
    args = 'relion_run_ctffind --i micrographs.star --only_make_star --use_gctf --o . '
    proc = subp.Popen( args,  shell = True, stdout = output  )
    proc.wait
    output.close()

    #input_star  = open ("micrographs_ctf.star", 'r' )
    #output_star = open ("micrographs_DW.star",  'w' )
    #for line in input_star:
        #output_star.write(line.replace('.mrc', '_DW.mrc')) 
    #    print (line)

    defocus1_list = []
    defocus2_list = []
    defocus_angle_list = []
    CCC_list = []
    resolution_list = []
    for file in glob.glob("*_gctf.log") :
        myfile = open (file, "r" )
        b  =  ' '.join(grep("Final",myfile ))
        defocus1_list.append(b.split()[0])
        defocus2_list.append(b.split()[1])
        defocus_angle_list.append(b.split()[2])
        CCC_list.append(b.split()[3])
    for file in glob.glob("*_gctf.log") :
        myfile = open (file, "r" )
        b  =  ' '.join(grep("RES_LIMIT",myfile ))
        resolution_list.append(b.split()[-1])
    defocus1_list = list(map(float, defocus1_list))
    defocus2_list = list(map(float, defocus2_list))
    defocus_angle_list = list(map(float, defocus_angle_list))
    CCC_list = list(map(float, CCC_list))
    resolution_list = list(map(float, resolution_list))
    myfile.close()
    mean_defocus1 = round ( sum (defocus1_list) / len (defocus1_list), 1 )
    mean_defocus2 = round ( sum (defocus2_list) / len (defocus2_list), 1 )
    mean_defocus_angle = round ( sum (defocus_angle_list) / len (defocus_angle_list), 1 )
    mean_CCC = round ( sum (CCC_list) / len (CCC_list), 1 )
    mean_resolution = round ( sum (resolution_list) / len (resolution_list),1 )
    print (sg)
    print ( '\nMean Defocus    :', mean_defocus1 )   
    print ( 'Mean FOM        :', mean_CCC )   
    print ( 'Mean Resolution :', mean_resolution )   
    print (eb)

def make_template_from3D ( model ) :
    args = 'e2project3d.py --outfile model.mrcs ' + model + ' --orientgen=eman:delta=45  --sym=c1 --postprocess=filter.lowpass.gauss:cutoff_freq=0.025'
    output = open("maketemplate.log", 'w')
    proc = subp.Popen( args,  shell = True,  stdout = output)
    while not os.path.exists('model.mrcs'):
        time.sleep(1)
    proc.wait
    output.close()


def Gautomatch ( miccrographs_list, pixel_size, diameter, cccutoff, contrast, template, lavgmin, lavgmax ) :
    #output_dir = quick2d_parameters['workdir'] + '/particles_Gautomatch/'
    #os.makedirs(output_dir)
    #os.chdir (output_dir) 
    number_of_gpus = quick2d_parameters['ngpu']
    split_by_gpu = np.array_split(micrographs_list,number_of_gpus)
    import subprocess as subp
    import time
    output = open("gautomatch.log", 'w')
    ps = []
    for i in range(number_of_gpus) :
        args = 'Gautomatch --apixM ' +  str(pixel_size)+ ' --diameter ' + str(diameter)  + ' --speed 1 ' + ' --lsigma_cutoff 1.2 ' + lavgmin + ' ' +  lavgmax + ' ' \
               + ' --cc_cutoff ' + str(cccutoff) + ' '  + ' '.join(split_by_gpu[i]) + ' --gid ' + str (i) + ' ' + contrast + ' ' + template
        #print (args )
        proc = subp.Popen( args,  shell = True,  stdout = output)
        ps.append(proc)

    #Keep the user informed       
    a = len(glob.glob("*_automatch.star" ))
    b = len (micrographs_list )
    c = 0
    print ('\nPicking particles using Gautomatch\n')
    while ( a < b ) :
        if ( a > c ) :
            print ( a,'/',b, 'Micrographs are done' )
            c = a
        a = len(glob.glob("*_automatch.star" ))

    for p in ps:
        p.wait()
   
    particles = 0 
    for file in (glob.glob("*_automatch.box" )) :
        particles = particles +  sum(1 for line in open(file))
    print (sg + '\nTotal number of picked particles:', particles  )
    print (eb )

def relion_particles_extract (  box_size, contrast, starfile, diameter, pixel_size, scale ) :
    output_dir = quick2d_parameters['workdir'] + 'Particles/'
    if scale != "" :
        bgradius =  36
    else :
        bgradius = int ( round ( (float (box_size) * 0.75) / 2  ) ) 
    os.makedirs(output_dir)
    args = 'mpirun -np 16 relion_preprocess_mpi --i  ' + starfile +' --coord_dir ./ --coord_suffix _automatch.star --part_star ' + output_dir + 'particles.star --part_dir ../Particles/ --extract --extract_size ' \
            + str(box_size) +  ' --norm --bg_radius ' +  str(bgradius) + ' --white_dust -1 --black_dust -1 ' + contrast + ' ' +  scale
    #print (args)
    print ('\nExtracting particles using Relion\n')
    output = open("extract.log", 'w')
    #print (args)
    proc = subp.Popen( args,  shell = True,  stdout = output)

    #Keep the user informed       
    a = len(glob.glob("../Particles/*_extract.star" ))
    b = len (micrographs_list )
    c = 0
    while ( a < b ) :
        if ( a > c ) :
            print ( a,'/',b, 'Micrographs are done' )
            c = a
        a = len(glob.glob("../Particles/*_extract.star" ))
    print (sg + '\nParticle Extraction done')
    print (eb)
    proc.wait()

def relion_2dclass ( diameter, nclasses, ignorectf, Tval, run_number, particles  ) :
    output_dir = quick2d_parameters['workdir'] + 'Class2D/'
    if os.path.isdir (output_dir) != True :
        os.makedirs(output_dir)
    os.chdir(output_dir)
    args = 'mpirun -np 9 relion_refine_mpi --o ./run' + str(run_number) + ' --i ' + particles + ' --dont_combine_weights_via_disc --pool 100 --ctf  --pad 1  --iter 25 --tau2_fudge ' + str(Tval) + ' --K ' + str(nclasses)  \
           + ' --particle_diameter ' + str (diameter) + ' --flatten_solvent  --zero_mask  --oversampling 1 --psi_step 10 --offset_range 5 --offset_step 2  --dont_check_norm --scale  --j 2 --gpu ' + ignorectf
    #print (args)
    output = open("relion2d.log", 'w')
    print ('\n2D Classification using Relion\n')
    proc = subp.Popen( args,  shell = True,  stdout = output)

    #Keep the user informed       
    name = 'run' + str(run_number) + '*model.star'
    a = len(glob.glob( name ))
    b = 26
    c = 0
    while ( a < b ) :
        if ( a > c ) :
            print ( a,'/',b-1, 'Iterations are done' )
            c = a
        a = len(glob.glob(name ))

    proc.wait()

def resize_3dmodel ( inpix, outpix, outbox, inmrc ) :
    ratio = float (inpix) / float (outpix)
    args =  'sxprocess.py ' + inmrc + ' temp1.mrc --changesize --ratio=' + ratio 
    proc = subp.Popen( args,  shell = True,  stdout = output)
    while not os.path.exists('temp1.mrc'):
        time.sleep(1)
    proc.wait()

    args = 'e2proc3d.py --clip=' + outbox + ' temp1.mrc  resized.mrc '
    while not os.path.exists('resized.mrc'):
        time.sleep(1)
    proc.wait()

def make_clean ( datadir, workdir ) :
    remove_dir = datadir + '/' + workdir
    tempdir = datadir +'/'+workdir+'/*'
    print ('I am cleaning up the temp files' )
    for file in glob.glob(tempdir) :
        os.remove(file)
    #time.sleep(3)
    os.chdir ( datadir )
    shutil.rmtree( remove_dir, ignore_errors=True )
    

        
pattern = '*'+quick2d_parameters['micrograph_name_pattern'] +'*'+ quick2d_parameters['micrograph_name_suffix']
micrographs_total, micrographs_list = count_micrographs(pattern, quick2d_parameters['datadir'], quick2d_parameters['micrograph_name_exclude'], 1  )
if micrographs_total == 0 :
    quick2d_parameters['micrograph_name_suffix'] = 'mrc'
    pattern = '*'+quick2d_parameters['micrograph_name_pattern'] +'*'+ quick2d_parameters['micrograph_name_suffix']
    micrographs_total, micrographs_list = count_micrographs(pattern, quick2d_parameters['datadir'], quick2d_parameters['micrograph_name_exclude'], 1  )
    if micrographs_total == 0 :
        print ('I cant find the micrographs in the directory')
        print ('I will quit')
        quit ()

    
quick2d_parameters['z'] = get_movie_frame_number ( micrographs_list )
quick2d_parameters['ngpu'] = number_of_gpu ( micrographs_list  )

if quick2d_parameters['z'] == 1 :
    print ('\nYou have given micrographs rather than movies.. I wil start with CTF estimation')
    #make_dir_soft_links(quick2d_parameters['workdir'],quick2d_parameters['datadir'] , pattern, quick2d_parameters['micrograph_name_exclude']    )
    make_dir_soft_links(quick2d_parameters['workdir'], micrographs_list )
    micrographs_total = len(glob.glob(os.path.join(quick2d_parameters['workdir'],'micrographs','*' )))
    micrographs_list = glob.glob(os.path.join(quick2d_parameters['workdir'],'micrographs','*' ))
else :
    make_dir_soft_links1 ( quick2d_parameters['workdir'], micrographs_list )
    MotionCor2_align (quick2d_parameters['ngpu'] )
    aligned_micrographs_dir = quick2d_parameters['workdir'] + '/aligned/' 
    aligned_micrographs_suf = 'mrc'
    pattern = '*'+quick2d_parameters['micrograph_name_pattern'] +'*'+ aligned_micrographs_suf
    micrographs_total, micrographs_list = count_micrographs( pattern , aligned_micrographs_dir, 'DW', 0 )
    if quick2d_parameters['bin'] != 1 :
        quick2d_parameters['pixel_size'] = float ( quick2d_parameters['pixel_size'] ) * float ( quick2d_parameters['bin'] )

GCTF( quick2d_parameters['pixel_size'], quick2d_parameters['cs'], quick2d_parameters['kv'], micrographs_list, quick2d_parameters['negative'])

if quick2d_parameters['dose'] != 0 :
    pattern = '*'+quick2d_parameters['micrograph_name_pattern'] +'*DW*'+ aligned_micrographs_suf
    print (pattern)
    micrographs_total, micrographs_list = count_micrographs( pattern , aligned_micrographs_dir, quick2d_parameters['micrograph_name_exclude'], 0 )

#if quick2d_parameters['template'] != '' :
#    template = ' --T ' + quick2d_parameters['datadir'] + '/' + quick2d_parameters['template']
if  quick2d_parameters['3dmodel'] != '' :
    make_template_from3D ( quick2d_parameters['datadir'] + '/' + quick2d_parameters['3dmodel'] )
    template = ' --T  model.mrcs --apixT ' + quick2d_parameters['3dmodelpix'] 
elif quick2d_parameters['template'] != '' :
    template = ' --T ' + quick2d_parameters['datadir'] + '/' + quick2d_parameters['template']  + ' --apixT ' + quick2d_parameters['templatepix']
elif quick2d_parameters['template'] == '' :
    template = ' '

Gautomatch ( micrographs_list, quick2d_parameters['pixel_size'] , quick2d_parameters['diameter'], quick2d_parameters['cccutoff'], quick2d_parameters['pickcontrast'], template, quick2d_parameters['lavgmin'], quick2d_parameters['lavgmax'] )

if quick2d_parameters['dose'] != 0 :
    input_star  = open ("micrographs_ctf.star", 'r' )
    output_star = open ("micrographs_ctf_DW.star",  'w' )
    for line in input_star:
        output_star.write(line.replace('.mrc', '_DW.mrc')) 
    starfile = 'micrographs_ctf_DW.star'
    output_star.close()
    input_star.close()
else :
    starfile = 'micrographs_ctf.star'

if quick2d_parameters['box'] == 0 :
    quick2d_parameters['box'] =  round_up_to_even ( ( float ( quick2d_parameters['diameter'])   * 1.8 ) / float ( quick2d_parameters['pixel_size'] ) )

if quick2d_parameters['mask'] == 0 :
    quick2d_parameters['mask'] = float ( quick2d_parameters['diameter'] ) * 1.1
relion_particles_extract ( quick2d_parameters['box'], quick2d_parameters['contrast'], starfile, quick2d_parameters['diameter'], quick2d_parameters['pixel_size'], quick2d_parameters['partbin'] )

run_number = 1
particles = '../Particles/particles.star'
relion_2dclass ( quick2d_parameters['mask'],  quick2d_parameters['nclasses'], quick2d_parameters['ignorectf'], quick2d_parameters['Tval'], run_number, particles ) 



mail_author = 'echo " " | mail -s QUICK2D ' + author
os.system( mail_author )

args = 'relion_display --class  --i  run' + str (run_number) + '_it025_model.star --scale 1 --sort rlnClassDistribution --reverse --allow_save  \
           --fn_imgs ./classes' + str (run_number) + '.star  --fn_parts particles' + str(run_number) + '.star'
proc = os.system ( args )

print (sg + '\nPrcoessing is done.. ')
print ('\nResults are in:    ' +  quick2d_parameters['timestamp'] + '\n' )
print (eb)
