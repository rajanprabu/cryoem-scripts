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
author = 'test@test.com'
quick3d_parameters = {}
quick3d_parameters['datadir'] = os.getcwd()
quick3d_parameters['user'] = getpass.getuser()
quick3d_parameters['user_uid'] = os.getuid()
quick3d_parameters['micrograph_name_pattern'] = '*'
quick3d_parameters['micrograph_name_suffix'] = 'mrcs'
quick3d_parameters['micrograph_name_exclude'] = '\#'
quick3d_parameters['pixel_size'] = 0
quick3d_parameters['z'] = 0
quick3d_parameters['kv'] = 300
quick3d_parameters['cs'] = 2.62
quick3d_parameters['diameter'] = 0
quick3d_parameters['box'] = 0
quick3d_parameters['cccutoff'] = 0.2
quick3d_parameters['mask'] = 0
quick3d_parameters['dose'] = 0
quick3d_parameters['Tval'] = 2
quick3d_parameters['nclasses'] = 50
quick3d_parameters['ignorectf'] = ''
quick3d_parameters['negative'] = 0
quick3d_parameters['contrast'] = ' --invert_contrast '
quick3d_parameters['partbin'] = ' --scale 96 '
quick3d_parameters['pickcontrast'] = '  '
quick3d_parameters['preexposure'] = 0
quick3d_parameters['time'] =  time.strftime("%d_%b_%Y_%I_%M_%S%p")
quick3d_parameters['timestamp'] =  'processing_' + quick3d_parameters['time']
quick3d_parameters['workdir'] = os.getcwd() + '/' + quick3d_parameters['timestamp'] + '/'
quick3d_parameters['power_users'] = 0
quick3d_parameters['write_movies'] = 'NO'
quick3d_parameters['ngpu'] = 0
quick3d_parameters['bin'] = 0
quick3d_parameters['ncpu'] = 0
quick3d_parameters['sparxcpu'] = 0
quick3d_parameters['gain'] = ''
quick3d_parameters['input_commands'] = ' '.join(sys.argv)
quick3d_parameters['rescut'] = 0
quick3d_parameters['dose_filter'] = 'NO' 
quick3d_parameters['template'] = '' 
quick3d_parameters['templatepix'] = '' 
quick3d_parameters['3dmodel'] = '' 
quick3d_parameters['3dmodelpix'] = '' 
quick3d_parameters['lavgmin'] = ' --lave_min -1 '
quick3d_parameters['lavgmax'] = ' --lave_max 1.2 '
quick3d_parameters['rotgain'] = 0
quick3d_parameters['qcpu'] = 0
quick3d_parameters['sparx'] = 0
quick3d_parameters['relion_to_sparx'] = 0
quick3d_parameters['reextract'] = 0
quick3d_parameters['bin_pixel_size'] = 0
quick3d_parameters['gfact'] = 1
quick3d_parameters['ask_pixel_size'] = 0
quick3d_parameters['krios'] = 0
quick3d_parameters['talos'] = 0
quick3d_parameters['listfile'] = ''
quick3d_parameters['3dmask'] = ''
quick3d_parameters['dwcombi'] = 0
quick3d_parameters['3d'] = ''
quick3d_parameters['auto'] = 0
quick3d_parameters['auto_2d_iter'] = 2
quick3d_parameters['relion_2d_iter'] = 25
quick3d_parameters['relion_3d_iter'] = 25
quick3d_parameters['number_of_3dclasses'] = 3

mail_author = 'echo " " | mail -s QUICK3D ' + author
os.system( mail_author )

if quick3d_parameters['user_uid'] not in ( 4421, 3623, 3336, 5021 ) :
    print ("\nThis script is in beta testing. Please contact rajan if you want to test it\n") 
    quit ()


### Command line and help text
### Argparse
from optparse import OptionParser

usage = "usage %prog "+ """ [options] 

        Examples:
        
        Minimal input: pixel size and diameter
        quick3d.com --pix=1.99 --dia=200
    
        Different kv and cs
        quick3d.com --pix=1.99 --dia=200 --kv=200 --cs=2.0

        For negative stain data ( default is cryo data )
        quick3d.com --pix=1.99 --dia=200 --negative

        By defaults particles will be binned to 96 px.. If you want to to keep original particle size  
        quick3d.com --pix=1.99 --dia=200 --nopartbin

        Give custom box size ( in pixels ) and mask size ( A )  
        quick3d.com --pix=1.99 --dia=200 --box=200 --mask=250

        Give template and its pixels size for autopicking 
        quick3d.com --pix=1.99 --dia=200 --template=class.mrcs --Tpix=2.37

        Give 3D Model and its pixels size, projections are made and used for autopicking 
        quick3d.com --pix=1.99 --dia=200 --3dmodel=run_it025.mrc --3dmodelpix=2.37

        For the dose wighting 
        quick3d.com --pix=1.99 --dia=200 --dose=22

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
gen.add_option("--dwcombi",         dest="dwcombi",           action='store_true', default=False,  help="estimate ctf from micrographs but extrat from dose weighted ones")
parser.add_option_group(gen)

mic = optparse.OptionGroup(parser, 'Microscopy')
mic.add_option("--kv",                 dest="kv",                help="Voltage of the Microscope  eg  --kv=300  default 300 ")
mic.add_option("--cs",                 dest="cs",                help="Spherical aberation of the Microscope  eg  --cs=2.7  default  2.7  ")
mic.add_option("--krios",              dest="krios",             action="store_true", default=False, help="will set kv cs and RotGain")
mic.add_option("--talos",              dest="talos",             action="store_true", default=False, help="will set kv cs and RotGain")
parser.add_option_group(mic)

select = optparse.OptionGroup(parser, 'Micrograph Selection')
select.add_option("--mic",             dest="pattern",           help="Common unique pattern in the micrograph name ( to identify/select)  eg  --m=_frames default all mrc files in the directory")
select.add_option("--exc",             dest="exclude_pattern",   help="Micrograph name that you want to exclude    eg  --exc=_DW No default ")
select.add_option("--suf",             dest="suffix",            help="Micrograph suffix  eg --suf=mrc  default mrcs or mrc ")
select.add_option("--list",            dest="list",              help="Micrograph list file  eg --list=file1.list ")
parser.add_option_group(select)

movie = optparse.OptionGroup(parser, 'Movies')
movie.add_option("--dose",             dest="dose",              help="Exposure per frame e/A^2 --dose=2 No default ")
movie.add_option("--gain",             dest="gain",              help="Camera gain file for correction eg --gain=K2gain.em  No default")
movie.add_option("--bin",              dest="bin",               help="Binning factor for micrographs default - no binning")
movie.add_option("--stack",            dest="stack",             action="store_true", default=False,  help="Write alinged movie stacks eg --stack Default false")
movie.add_option("--preexp",           dest="preexposure",       help=" Pre Exposure per frame e/A^2 No default")
movie.add_option("--rotgain",          dest="rotgain",           help=" Gain rotation No default")
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

class3d = optparse.OptionGroup(parser, '3D Classification')
class3d.add_option("--3d",              dest="model3d",             help="3D Model  eg --3d=m.mrc - No Default ")
class3d.add_option("--3dmask",          dest="mask3d",              help="3D Mask  eg --3d=m.mrc - No Default ")
class2d.add_option("--3dclass",         dest="nclasses3d",          help="Number of 3D classes to calculate eg --3dclass=3 Default  5")
parser.add_option_group(class3d)

computing = optparse.OptionGroup(parser, 'Computing')
computing.add_option("--cpu",          dest="ncpu",              help="Number of cpus to use eg --cpu=20 Default use all available")
computing.add_option("--qcpu",          dest="qcpu",             help="Number of cpus to use for quesu submission eg --qcpu=40 Default 40")
computing.add_option("--gpu",          dest="ngpu",              help="Number of gpus to use eg --gpu=2 defaults use all available")
computing.add_option("--devel",        dest="power",             action="store_true", default=False, help="For developers....")
computing.add_option("--sparx",        dest="sparx",             action="store_true", default=False, help="For developers....")
computing.add_option("--rel2sparx",    dest="rel2sparx",         action="store_true", default=False, help="For developers....")
computing.add_option("--auto",         dest="auto",              action='store_true', default=False,  help="auto pilot mode")
computing.add_option("--rcut",         dest="rescut",            help="only works with --devel")
parser.add_option_group(computing)

(options, args) = parser.parse_args()

if options.kv :
    quick3d_parameters['kv'] = options.kv

if options.pixel_size :
    quick3d_parameters['pixel_size'] = options.pixel_size
    quick3d_parameters['full_pixel_size'] = options.pixel_size

if options.dose :
    quick3d_parameters['dose'] = options.dose

if options.cccutoff :
    quick3d_parameters['cccutoff'] = options.cccutoff

if options.box :
    quick3d_parameters['box'] = options.box

if options.Tval :
    quick3d_parameters['Tval'] = options.Tval

if options.template :
    quick3d_parameters['template'] = options.template

if options.templatepix :
    quick3d_parameters['templatepix'] = options.templatepix

if options.list :
    quick3d_parameters['listfile'] = options.list

if options.model :
    quick3d_parameters['3dmodel'] = options.model

if options.nclasses :
    quick3d_parameters['nclasses'] = options.nclasses

if options.lavgmin :
    quick3d_parameters['lavgmin'] = ' --lave_min ' + options.lavgmin

if options.lavgmax :
    quick3d_parameters['lavgmax'] = ' --lave_max ' + options.lavgmax

if options.modelpix :
    quick3d_parameters['3dmodelpix'] = options.modelpix

if options.model3d :
    quick3d_parameters['3d'] = options.model3d

if options.mask3d :
    quick3d_parameters['3dmask'] = options.mask3d

if options.mask :
    quick3d_parameters['mask'] = options.mask

if options.rotgain :
    quick3d_parameters['rotgain'] = options.rotgain

if options.diameter :
    quick3d_parameters['diameter'] = options.diameter

if options.cs :
    quick3d_parameters['cs'] = options.cs

if options.preexposure :
    quick3d_parameters['preexposure'] = options.preexposure

if options.pattern :
    quick3d_parameters['micrograph_name_pattern'] = options.pattern

if options.exclude_pattern :
    quick3d_parameters['micrograph_name_exclude'] = options.exclude_pattern

if options.gain :
    quick3d_parameters['gain'] = options.gain

if options.bin :
    quick3d_parameters['bin'] = options.bin

if options.suffix :
    quick3d_parameters['micrograph_name_suffix'] = options.suffix

if options.ncpu :
    quick3d_parameters['ncpu'] = options.ncpu

if options.qcpu :
    quick3d_parameters['qcpu'] = options.qcpu

if options.ngpu :
    quick3d_parameters['ngpu'] = options.ngpu

if options.rescut :
    quick3d_parameters['rescut'] = options.rescut

if str(options.stack) == "True" :
    quick3d_parameters['write_movies'] = 'YES'

if str(options.ignorectf) == "True" :
    quick3d_parameters['ignorectf'] = ' --ctf_intact_first_peak '

if str(options.dwcombi) == "True" :
    quick3d_parameters['dwcombi'] = 1

if str(options.auto) == "True" :
    quick3d_parameters['auto'] = 1

if str(options.partbin) == "False" :
    quick3d_parameters['partbin'] = ''

if str(options.negative) == "True" :
    quick3d_parameters['contrast'] = ' '
    quick3d_parameters['negative'] = 1
    quick3d_parameters['pickcontrast'] = ' --dont_invertT '

if str(options.power) == "True" :
    quick3d_parameters['power_users'] = 1

if str(options.sparx) == "True" :
    quick3d_parameters['sparx'] = 1

if str(options.rel2sparx) == "True" :
    quick3d_parameters['relion_to_sparx'] = 1

if str(options.krios) == "True" :
    quick3d_parameters['krios'] = 1

if str(options.talos) == "True" :
    quick3d_parameters['talos'] = 1

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

if quick3d_parameters['ask_pixel_size'] == 1 :
    if quick3d_parameters['pixel_size'] == 0 :
        print ( "\nI need pixel size as a input, I kind of assume defaults for other values" )
        print ( "please provide pixel size For eg. --pix=1.43 or -p 1.43\n" )
        print ( "If you need help please use quick3d.com -h" )
        quit ()

if quick3d_parameters['template'] != '' :
    file = quick3d_parameters['datadir'] + '/' + quick3d_parameters['template']
    if os.path.exists (file) != True :
        print ("\nI cant fine the template file for particle picking")
        quit ()
    elif quick3d_parameters['templatepix'] == '' :
        print ("\n You have given template but not the pixel size" )
        print ("\n For eg --tpix=2.54 ")
        quit ()

if quick3d_parameters['gain'] != '' :
    file = quick3d_parameters['datadir'] + '/' + quick3d_parameters['gain']
    if os.path.exists (file) != True :
        print ("\nI cant fine the Gain reference quick3d_parameters['gain']")
        quit ()
    else :
        quick3d_parameters['gain'] = file


if quick3d_parameters['gain'] == '' and  quick3d_parameters['krios'] == 1 :
        gain_list = glob.glob('*.dm4')
        if len (gain_list) == 1 :
            gain = '\n'.join(str(p) for p in gain_list)
            quick3d_parameters['gain'] = input ('\nPossible Gain available in this directory is  ( ' + gain +   ' )  :')
            quick3d_parameters['gain'] =  quick3d_parameters['datadir'] + '/' + gain
        else :
            print ("\nI cant find the Gain reference quick3d_parameters['gain']")
            quit ()

if quick3d_parameters['3dmodel'] != '' :
    file = quick3d_parameters['datadir'] + '/' + quick3d_parameters['3dmodel']
    if os.path.exists (file) != True :
        print ("\nI cant fine the template file for particle picking")
        quit ()
    elif quick3d_parameters['3dmodelpix'] == '' :
        print ("\n You have given 3dmodel but not the pixel size" )
        print ("\n For eg --3dmodelpix=2.54 ")
        quit ()

#if quick3d_parameters['box'] != 0 :
#    quick3d_parameters['partbin'] = ''
#    print ("\n You have given box size, so no particles binning" )
#elif quick3d_parameters['mask'] != 0 :
#    quick3d_parameters['partbin'] = ''
#    print ("\n You have given mask size, so no particles binning" )


if quick3d_parameters['krios'] == 1 :
    quick3d_parameters['cs'] = 2.62
    quick3d_parameters['kv'] = 300
    quick3d_parameters['rotgain'] = 0

if quick3d_parameters['talos'] == 1 :
    quick3d_parameters['cs'] = 2.62
    quick3d_parameters['kv'] = 200

if quick3d_parameters['diameter'] == 0 :
    print ( "\nI need particle diameter as a input " )
    print ( "please provide particle dimaeter For eg. --dia=128 or -d 128\n" )
    print ( "If you need help please use quick3d.com -h" )
    quit ()

if quick3d_parameters['dose'] != 0 :
    print ( "\nI will do dose weighting" )
    quick3d_parameters['dose_filter'] = 'TRUE' 
    if quick3d_parameters['kv'] == 0 :
        print ( "please provide KV of the microscope For eg. --kv=200 or -k 200\n" )
        print ( "If you need help please use quick3d.com -h" )

#if quick3d_parameters['gain'] == '' :
#    print ( "\nI need camera gain, I kind of assume defaults for other values" )
#    print ( "please provide the gain For eg. --g=Gain.em or -g Gain.em\n" )
#    print ( "If you need help please use quick3d.com -h" )
    #quit ()

## CHECK FOR THE PRESENSE OF EXECUTABLES
def check_execs() :

    exec_list = ['header', 'Gctf', 'relion_refine_mpi', 'relion_preprocess_mpi', 'MotionCor2', 'Gautomatch','e2proc2d.py','relion_display']

    for i in range (len(exec_list)):
        if not check_for_executables(exec_list[i]) :
            quick3d_parameters.setdefault('missing_execs', []).append(exec_list[i])

    if 'missing_execs' in quick3d_parameters :
        print ('\n\nThe following exeucutables are missing', quick3d_parameters['missing_execs'] )
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
        #print ('I cant find the micrographs in the directory')
        #print ('I will quit')
        #quit ()
    #elif printflag == 1 :
        print ('\nYou have',micrographs_total, 'micrographs')
    return micrographs_total,micrographs_list

def count_micrographs_list(listfile) :
    micrographs_list = open(listfile).read().splitlines()
    micrographs_total = len( micrographs_list )
    print ('\nYou have',micrographs_total, 'micrographs')
    return micrographs_total,micrographs_list


def number_of_gpu( micrographs_list ) :
    gpus = len(glob.glob('/proc/driver/nvidia/gpus/*'))
    if gpus == 0 :
        print ('I cant find any usable GPU')
        print ('I will quit')
        quit ()
     
    if len (micrographs_list) < gpus :
        quick3d_parameters['ngpu']  = len (micrographs_list)

    quick3d_parameters['ngpu']  = int( quick3d_parameters['ngpu'] )
    if quick3d_parameters['ngpu'] == 0 :
        print (' Available GPUS: ',gpus )
    else :
        if ( gpus != 0 ) :
            print (' Available GPUS: ',gpus )
            print (' Using ', quick3d_parameters['ngpu'], 'GPUS' )
            gpus = quick3d_parameters['ngpu']
    return gpus

def number_of_cpu() :
    cpus = multiprocessing.cpu_count()
    quick3d_parameters['ncpu']  = int (quick3d_parameters['ncpu'] )
    if quick3d_parameters['ncpu'] == 0 :
        print (' Available CPUS: ',cpus )
    else :
        if ( cpus != 0 ) :
            cpus = quick3d_parameters['ncpu']
            print (' Using ', quick3d_parameters['ncpu'], 'CPUS' )
    return cpus

def get_pixel_size ( micrographs_list ) :
    args = 'header -pixel -i ' + micrographs_list[0]
    proc = subp.Popen( args,  shell = True,  stdout=subp.PIPE )
    proc.wait
    p  =  float ( proc.communicate()[0].decode("utf-8").split()[2] )
    return p

def get_movie_frame_number ( micrographs_list ) :
    args = 'header -size -i ' + micrographs_list[0] 
    proc = subp.Popen( args,  shell = True,  stdout=subp.PIPE ) 
    proc.wait
    z  =  int ( proc.communicate()[0].decode("utf-8").split()[2] )
    return z

def make_qsub ( name, wdir, command ) :
    output_dir = quick3d_parameters['workdir'] + 'logs'
    if os.path.isdir (output_dir) != True :
        os.makedirs(output_dir)
    #os.chdir(quick3d_parameters['workdir'])

    qsub_input = []
    qsub_input.append ( '#!/bin/bash' )
    qsub_input.append ( '#$ -pe openmpi ' + quick3d_parameters['qcpu'] )
    qsub_input.append ( '#$ -N  ' + name )
    qsub_input.append ( '#$ -l h_rt=604800') 
    qsub_input.append ( '#$ -l h_vmem=45G ' )
    qsub_input.append ( '#$ -e ' + name + '.err' ) 
    qsub_input.append ( '#$ -o ' + name + '.out' )
    qsub_input.append ( '#$ -cwd ')
    qsub_input.append ( '#$ -S /bin/bash ' )
    qsub_input.append ( '#$ -j y ' )
    qsub_input.append ( '#$ -m n ')
    qsub_input.append ( 'cd ' + wdir )
    qsub_input.append ( 'source /fs/pool/pool-schulman-soft/setup/source.sh')
    qsub_input.append ( command )
    i = 0
    out = output_dir + '/' +  name + '.submit'
    com = open (out, "w")
    while (i < len (qsub_input) ):
        com.write( str (qsub_input[i] ) )
        com.write("\n"  )
        i += 1
    st = os.stat(out)
    os.chmod(out, st.st_mode | stat.S_IEXEC)
    com.close()

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
    if os.path.exists (workdir) != True :
        os.makedirs(workdir)
    if os.path.exists (mic_dir) != True :
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
    command = workdir + 'input'
    output = open(command, 'w')
    output.writelines ( quick3d_parameters['input_commands'] + '\n' )
    output.close


    os.chdir ( mic_dir )

def make_dir_soft_links1 ( workdir, micrographs_list ) :
    dest_dir  = workdir + '/'
    movie_dir  = workdir + '/movies'
    os.makedirs(workdir)
    os.makedirs(movie_dir)
    number_of_gpus = quick3d_parameters['ngpu']
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

def micrographs_list_nodir_name (micrographs_list ) :
    micrographs_list_simple = []
    for file in micrographs_list :
        micname =  file.split('/')[-1]
        micrographs_list_simple.append(micname)
    return micrographs_list_simple
 
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

    the_dot_with_suff = '.' + quick3d_parameters['micrograph_name_suffix']
    for mics in micrographs_list :
       micname =     mics.split('/')[-1]
       write_unblur_input( micname, quick3d_parameters['z'], quick3d_parameters['pixel_size'], quick3d_parameters['dose_filter'], \
                  quick3d_parameters['dose'], quick3d_parameters['kv'], quick3d_parameters['preexposure'], quick3d_parameters['write_movies']  )
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

def MotionCor2_align_script ( number_of_gpus, suffix, gain, gainrot, micrographs_list ) :
    output_dir = quick3d_parameters['workdir'] + '/aligned/'
    os.makedirs(output_dir)
    os.chdir (output_dir) 
    output = open("MotionCor2.log", 'w')
    jobs = number_of_gpus * int(quick3d_parameters['gfact'])
    split_by_jobs = np.array_split(micrographs_list, jobs)
    if quick3d_parameters['gain'] != "" :
        k2gain = ' -Gain ' + gain + ' '
    else :
        k2gain = ""

    if quick3d_parameters['bin'] != 0 :
        binning = ' -Ftbin ' + quick3d_parameters['bin'] + ' ' 
    else :
        binning = ''

    if gainrot != 0 :
       rotgain = ' -RotGain ' + str (gainrot)
    else :
       rotgain = ' -RotGain 0 '

    if suffix == 'tif' :
        input_stem = ' -InTiff '
    else :
        input_stem = ' -InMrc '
     
    if gain.split('.')[-1] == 'dm4' :
        g =  os.path.basename(gain)
        args = 'dm2mrc ' + gain + ' ../' +  g + '.mrc'
        proc = subp.Popen( args,  shell = True,  stdout = output)
        proc.wait
        k2gain = ' -Gain ../' + g + '.mrc '

    for i in range(jobs) :
        gpuid = int ( i / quick3d_parameters['gfact'] )
        out = 'motioncor2_gpu' + str(i) + '.sh'
        com = open (out, "w")
        com.write("#!/bin/bash \n"  )
        for j in (split_by_jobs[i]) :
            k =  os.path.basename(j).split('.')[0]
            args = 'MotionCor2 ' + input_stem  + j  + ' -OutMrc ' + output_dir + k + '.mrc' + ' -patch 5 5  -Tol 0.5 -Iter 10  -Gpu ' + str (gpuid)  \
                 + binning  + k2gain + rotgain + ' -InitDose ' + str (quick3d_parameters['preexposure'])  + ' -FmDose ' + str(quick3d_parameters['dose']) \
                 + ' -PixSize '+ str (quick3d_parameters['pixel_size']) + ' -kV ' + str (quick3d_parameters['kv']) + ' -Throw 1 '
            com.write(str(args))
            com.write("\n")
        st = os.stat(out)
        os.chmod(out, st.st_mode | stat.S_IEXEC)
        com.close()
    ps = []
    for i in range(jobs) :
        args = 'motioncor2_gpu'+ str (i) + '.sh'
        proc = subp.Popen( args,  shell = True,  stdout = output)
        time.sleep(0.1)
        ps.append(proc)


    #Keep the user informed       
    if quick3d_parameters['dose'] == 0 :
        #list = output_dir + '/*' + quick3d_parameters['micrograph_name_suffix']
        list = output_dir + '/*.mrc'
    else :
        #list = output_dir + '/*_DW.' + quick3d_parameters['micrograph_name_suffix']
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

    if quick3d_parameters ['micrograph_name_suffix'] == 'mrcs' :
        mrcs_dir_name = output_dir + '/*mrcs'
        #names_list = [w.replace(suf, '') for w in micrographs_list]
        for mrcs in glob.glob(mrcs_dir_name) :
            name = mrcs.replace('mrcs','')
            mrc = name + 'mrc'
            os.rename ( mrcs, mrc)
    print (sg + '\nMovie alignment done')
    print (eb)



def MotionCor2_align ( number_of_gpus ) :
    output_dir = quick3d_parameters['workdir'] + '/aligned/'
    os.makedirs(output_dir)
    os.chdir (output_dir) 
    if quick3d_parameters['gain'] != "":
        gain = ' -Gain ' + quick3d_parameters['gain'] + ' '
    else :
        gain = ""

    if quick3d_parameters['micrograph_name_suffix'] == "tif":
        input_stem = ' -InTiff '
    else :
        input_stem = ' -InMrc '

    output = open("MotionCor2.log", 'w')
    ps = []
    for i in range(number_of_gpus) :
        input_dir  = quick3d_parameters['workdir'] + '/movies/gpu' + str (i) + '/'
        outlog = 'gpu' + str (i) + '.log'
        output = open(outlog, 'w')
        args = 'MotionCor2 ' + input_stem  + input_dir  + ' -OutMrc ' + output_dir + ' -patch 5 5 20  -Tol 0.5 -Serial 1  -Gpu ' + str (i )  \
                 + ' -FtBin ' + str(quick3d_parameters['bin']) + gain + ' -InitDose ' + str (quick3d_parameters['preexposure'])  + ' -FmDose ' + str(quick3d_parameters['dose']) \
                 + ' -PixSize '+ str (quick3d_parameters['pixel_size']) + ' -kV ' + str (quick3d_parameters['kv']) + ' -Throw 1 '
        output.writelines ( str(args ) + '\n' )
        proc = subp.Popen( args,  shell = True,  stdout = output)
        #print (args)
        ps.append(proc)

    #Keep the user informed       
    if quick3d_parameters['dose'] == 0 :
        #list = output_dir + '/*' + quick3d_parameters['micrograph_name_suffix']
        list = output_dir + '/*.mrc'
    else :
        #list = output_dir + '/*_DW.' + quick3d_parameters['micrograph_name_suffix']
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

    if quick3d_parameters ['micrograph_name_suffix'] == 'mrcs' :
        mrcs_dir_name = output_dir + '/*mrcs'
        #names_list = [w.replace(suf, '') for w in micrographs_list]
        for mrcs in glob.glob(mrcs_dir_name) :
            name = mrcs.replace('mrcs','')
            mrc = name + 'mrc'
            os.rename ( mrcs, mrc)    
    print (sg + '\nMovie alignment done')
    print (eb)

def GCTF ( pixel_size, cs, kv, micrographs_list, negative ) :
    #output_dir = quick3d_parameters['workdir'] + '/ctf_Gctf/'
    #os.makedirs(output_dir)
    #os.chdir (output_dir) 
    number_of_gpus = quick3d_parameters['ngpu']
    jobs = number_of_gpus * int(quick3d_parameters['gfact'])
    split_by_jobs = np.array_split(micrographs_list, jobs)


    #split_by_gpu = np.array_split(micrographs_list,number_of_gpus)
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
        ac = ' --ac 0.1  '
    output = open("gctf.log", 'w')
    ps = []
    #for i in range(number_of_gpus) :
    for i in range(jobs) :
        gpuid = int ( i / quick3d_parameters['gfact'] )
        star = 'gpu' + str (i) + '.star'
        args = 'Gctf --apix '+ str(pixel_size)+ ' --kV ' + str(kv) +' --cs ' + str(cs) + ' ' + ' '.join(split_by_jobs[i]) \
               + ' --gid ' + str(gpuid) + ' --do_validation --do_EPA --astm 100  --resL 30 ' + ' --resH ' + str(resH) + ' --B_resH ' + str(B_resH) \
               + box + ac + ' --ctfstar ' + star
               #+ ' --ac 0.07 --do_EPA  --boxsize 512 --do_Hres_ref --Href_resL 20'  
        #print (args)
        output.writelines ( str(args ) + '\n' )
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
        g.writelines (   mic + '\n' )
    g.close()
	    

    output = open("star.log", 'w')
    args = 'relion_run_ctffind --i micrographs.star --only_make_star --use_gctf --o .'
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

def sparx_ctf ( pixel_size, cs, kv, micrographs_list, negative, cpu ) :
    if ( negative == 1 ) :
        ac = '35'
    else :
        ac = '10'
    output = open("sxcter.log", 'w')
    args = 'mpirun -np ' + str(cpu)  + ' sxcter.py "./*.' + quick3d_parameters['micrograph_name_suffix'] +  '" ctf --wn=512 --apix=' + str(pixel_size) + ' --Cs=' + str(cs) + '  --voltage=' + str(kv) + ' --ac=' + str(ac)
    output.writelines ( str(args ) + '\n' )
    if quick3d_parameters['qcpu'] != 0 :
        make_qsub( 'sxcter', os.getcwd() , args )
    else :
        proc = subp.Popen( args,  shell = True,  stdout = output)
        #Keep the user informed       
        a = len(glob.glob("ctf/pwrot/*txt" ))
        b = len (micrographs_list )
        c = 0
        print ('\nCalculating CTF using sxcter\n')
        while ( a < b ) :
            if ( a > c ) :
                print ( a,'/',b, 'Micrographs are done' )
                c = a
            a = len(glob.glob("ctf/pwrot/*txt" ))

        print (sg + '\nCTF Estimation done')
        print (eb)

        proc.wait
        output.close()

def make_template_from3D ( model ) :
    args = 'e2project3d.py --outfile model.mrcs ' + model + ' --orientgen=eman:delta=45  --sym=c1 --postprocess=filter.lowpass.gauss:cutoff_freq=0.025'
    output = open("maketemplate.log", 'w')
    output.writelines ( str(args ) + '\n' )
    proc = subp.Popen( args,  shell = True,  stdout = output)
    while not os.path.exists('model.mrcs'):
        time.sleep(1)
    proc.wait
    output.close()


def Gautomatch ( micrographs_list, pixel_size, diameter, cccutoff, contrast, template, lavgmin, lavgmax ) :
    #output_dir = quick3d_parameters['workdir'] + '/particles_Gautomatch/'
    #os.makedirs(output_dir)
    #os.chdir (output_dir) 
    number_of_gpus = quick3d_parameters['ngpu']
    split_by_gpu = np.array_split(micrographs_list,number_of_gpus)
    import subprocess as subp
    import time
    output = open("gautomatch.log", 'w')
    ps = []
    for i in range(number_of_gpus) :
        args = 'Gautomatch --apixM ' +  str(pixel_size)+ ' --diameter ' + str(diameter)  + ' --speed 2 ' + \
                ' --cc_cutoff ' + str(cccutoff) + ' '  + ' '.join(split_by_gpu[i]) + ' --gid ' + str (i) + ' ' + contrast + ' ' + template
           # ' --lsigma_cutoff 1.3 ' + lavgmin + ' ' +  lavgmax + ' ' \
        #print (args )
        output.writelines ( str(args ) + '\n' )
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
    output_dir = quick3d_parameters['workdir'] + 'Particles/'
    if scale != "" :
        bgradius =  36
    else :
        bgradius = int ( round ( (float (box_size) * 0.75) / 2  ) ) 
    os.makedirs(output_dir)
    args = 'mpirun -np 16 relion_preprocess_mpi --i  ' + starfile +' --coord_dir ./ --coord_suffix _automatch.star --part_star ' + output_dir + 'particles.star --part_dir ../Particles/ --extract --extract_size ' \
            + str(box_size) +  ' --norm --bg_radius ' +  str(bgradius) + ' --white_dust 3 --black_dust 3 ' + contrast + ' ' +  scale
    #print (args)
    print ('\nExtracting particles using Relion\n')
    output = open("extract.log", 'w')
    output.writelines ( str(args ) + '\n' )
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

def relion_particles_reextract (  box_size, contrast, diameter, run_number, dose ) :
    output_dir = quick3d_parameters['workdir'] + 'reextract/'
    bgradius = int ( round ( (float (box_size) * 0.75) / 2  ) )
    os.makedirs(output_dir)
    if dose == 0 :
        micrographs = 'micrographs_ctf.star'
    else :
        micrographs = 'micrographs_ctf_DW.star'
    args = 'mpirun -np 40 relion_preprocess_mpi --i ' +  micrographs + ' --reextract_data_star ../Class2D/particles' + str(run_number) + '.star --recenter --part_star ' + output_dir + 'particles.star --part_dir ../reextract/ --extract --extract_size ' \
            + str(box_size) +  ' --norm --bg_radius ' +  str(bgradius) + ' --white_dust 3 --black_dust 3 ' + contrast
    print ('\nRe-extracting particles using Relion\n')
    output = open("reextract.log", 'w')
    output.writelines ( str(args ) + '\n' )
    #print (args)
    proc = subp.Popen( args,  shell = True,  stdout = output)
    proc.wait()


def sparx_particles_extract (  box_size, contrast, cpu ) :
    output_dir = quick3d_parameters['workdir'] + 'sparx-particles/'
    os.makedirs(output_dir)
    os.chdir(output_dir)
    if  contrast != "" :
        contrast = ' --skip_invert'
    args = 'mpirun -np ' + str(cpu) + ' sxwindow.py "../micrographs/*.mrc"  "../micrographs/*_automatch.box" "../micrographs/ctf/partres.txt" "extract" --box_size=' + str(box_size) + contrast 
    #print (args)
    output = open("sxwindow.log", 'w')
    output.writelines ( str(args ) + '\n' )
    if quick3d_parameters['qcpu'] != 0 :
        make_qsub( 'sxwindow', os.getcwd() , args )
    else :
        print ('\nExtracting particles using sxwindow\n')
        proc = subp.Popen( args,  shell = True,  stdout = output)
        b = 0
        for file in (glob.glob("../micrographs/*_automatch.box" )) :
            if ( sum(1 for line in open(file)) ) != 0 :
                b = b + 1

        a = len(glob.glob("extract/mpi_proc_*" ))
        #b = len (micrographs_list )
        c = 0
        while ( a < b ) :
            if ( a > c ) :
                print ( a,'/',b, 'Micrographs are done' )
                c = a
            a = len(glob.glob("extract/mpi_proc_*" ))
        print (sg + '\nParticle Extraction done')
        print (eb)
        proc.wait()
        
    args = 'e2bdb.py  extract/mpi_proc_*  --makevstack=bdb:particles#data'
    if quick3d_parameters['qcpu'] != 0 :
        make_qsub( 'e2bdb', os.getcwd() , args )
    else :
        proc = subp.Popen( args,  shell = True,  stdout = output)
        proc.wait()

def sparx_2dclass ( diameter, particles, cpu  ) :
    #output_dir = quick3d_parameters['workdir'] + 'isac2D/'
    #if os.path.isdir (output_dir) != True :
    #    os.makedirs(output_dir)
    os.chdir(quick3d_parameters['workdir'])
    #args = 'sxrelion2sparx.py ' + particles +  ' --create_stack' 
    #print (args)
    #output = open("sxrelion2sparx.log", 'w')
    #output.writelines ( str(args ) + '\n' )
    #print ('\n2D Classification using Relion\n')
    #print ('\nConverting particles to sparx format\n')
    #proc = subp.Popen( args,  shell = True,  stdout = output)
    #proc.wait() 
    #stack = output_dir + 'work/sparx_stack.hdf'
    #while not os.path.exists(stack):
    #    time.sleep(1)
    #print (sg + '\nParticle sparx format conversion done')
    #print (eb)
 
    radius = round_up_to_even ( float (diameter) /  ( 2 * float ( quick3d_parameters['pixel_size'] )) )
    args = 'mpirun -np ' + str(cpu) + ' sxisac2.py bdb:sparx-particles/particles#data  isac2D --radius=' + str(radius) + ' --CTF'
    if quick3d_parameters['qcpu'] != 0 :
        make_qsub( 'sxisac', os.getcwd() , args )
    else :
        output = open("isac.log", 'w')
        output.writelines ( str(args ) + '\n' )
        proc = subp.Popen( args,  shell = True,  stdout = output)
        print ('\nRunning ISAC\n')
        proc.wait() 

def relion_to_sparx_2dclass ( diameter, particles, cpu  ) :
    output_dir = quick3d_parameters['workdir'] + 'sparx-particles'
    if os.path.isdir (output_dir) != True :
        os.makedirs(output_dir)
    os.chdir(output_dir)
    args = 'sxrelion2sparx.py ' + particles +  ' --create_stack --output_dir=extract --star_section=data_'
    #print (args)
    output = open("sxrelion2sparx.log", 'w')
    output.writelines ( str(args ) + '\n' )
    print ('\nConverting particles to sparx format\n')
    proc = subp.Popen( args,  shell = True,  stdout = output)
    proc.wait() 
    stack = output_dir + '/extract/sparx_stack.hdf'
    while not os.path.exists(stack):
        time.sleep(1)
    print (sg + '\nParticle sparx format conversion done')
    print (eb)
    
    os.chdir(quick3d_parameters['workdir'])
    print ('\n2D Classification using isac\n')
    radius = round_up_to_even ( float (diameter) /  ( 2 * float ( quick3d_parameters['pixel_size'] )) )
    args = 'mpirun -np ' + str(cpu) + ' sxisac2.py sparx-particles/extract/sparx_stack.hdf  isac2D --radius=' + str(radius) + ' --CTF'
    if quick3d_parameters['qcpu'] != 0 :
        make_qsub( 'sxisac', os.getcwd() , args )
    else :
        output = open("isac.log", 'w')
        output.writelines ( str(args ) + '\n' )
        proc = subp.Popen( args,  shell = True,  stdout = output)
        print ('\nRunning ISAC\n')
        proc.wait()

def relion_2dclass ( diameter, nclasses, ignorectf, Tval, run_number, particles  ) :
    output_dir = quick3d_parameters['workdir'] + 'Class2D/'
    if os.path.isdir (output_dir) != True :
        os.makedirs(output_dir)
    os.chdir(output_dir)
    args = 'mpirun -np 17 relion_refine_mpi --o ./run' + str(run_number) + ' --i ' + particles + ' --dont_combine_weights_via_disc --pool 100 --ctf  --pad 2 --fast_subsets --iter ' + str(quick3d_parameters['relion_2d_iter']) + ' --tau2_fudge ' + str(Tval) + ' --K ' + str(nclasses)  \
           + ' --particle_diameter ' + str (diameter) + ' --flatten_solvent  --zero_mask  --oversampling 1 --psi_step 10 --offset_range 5 --offset_step 2  --dont_check_norm --scale  --j 4 --gpu ' + ignorectf
    #print (args)
    output = open("relion2d.log", 'w')
    output.writelines ( str(args ) + '\n' )
    print ('\n2D Classification using Relion\n')
    proc = subp.Popen( args,  shell = True,  stdout = output)

    #Keep the user informed       
    name = 'run' + str(run_number) + '*model.star'
    a = len(glob.glob( name ))
    b = quick3d_parameters['relion_2d_iter'] + 1
    c = 0
    while ( a < b ) :
        if ( a > c ) :
            print ( a,'/',b-1, 'Iterations are done' )
            c = a
        a = len(glob.glob(name ))

    proc.wait()

def resize_3dmodel ( inpix, outpix, outbox, inmrc ) :
    output_dir = quick3d_parameters['workdir'] + 'Class3D/'
    if os.path.isdir (output_dir) != True :
        os.makedirs(output_dir)
    os.chdir(output_dir)

    print ('\nRe sizing the input 3D model\n')
    ratio = float (inpix) / float (outpix)
    output = open("resize.log",'w')
    args =  'sxprocess.py ' + inmrc + ' temp1.mrc --changesize --ratio=' + str (ratio)
    output.writelines ( str(args ) + '\n' )
    proc = subp.Popen( args,  shell = True,  stdout = output)
    while not os.path.exists('temp1.mrc'):
        time.sleep(1)
    proc.wait()

    args = 'e2proc3d.py --clip=' + str(outbox) + ' temp1.mrc  resized.mrc '
    output.writelines ( str(args ) + '\n' )
    proc = subp.Popen( args,  shell = True,  stdout = output)
    while not os.path.exists('resized.mrc'):
        time.sleep(1)
    proc.wait()

def relion_inital ( diameter, particles ) :
    output_dir = quick3d_parameters['workdir'] + 'InitialModel/'
    if os.path.isdir (output_dir) != True :
        os.makedirs(output_dir)
    os.chdir(output_dir)

    args = 'mpirun -np 9 relion_refine_mpi --o ./run1 --sgd  --subset_size 200 --strict_highres_sgd 20 --write_subsets 10 --denovo_3dref  \
            --i ' + particles + ' --ctf --sym C1 --zero_mask --dont_combine_weights_via_disc --pool 100 --iter 1 --particle_diameter ' + str(diameter) + \
             ' --oversampling 1 --healpix_order 1 --offset_range 10 --offset_step 4 --j 2 --gpu '
    output = open("relion-model.log", 'w')
    output.writelines ( str(args ) + '\n' )
    print ('\nInitial model using Relion\n')
    proc = subp.Popen( args,  shell = True,  stdout = output)

    proc.wait()

def relion_auto_select ( modelstar, datastar, run_number, maxp ) :
    ### should write more pythonic.. I am in a hurry.. sorry
    args = 'awk \'{if (NF>2 && $5 > 0 && $5 < 30 && $3 < 5 && $6 < 5 ) print $0}\' ' + modelstar + '| awk -F@ \'{print $1+0 }\' > list' + str(run_number)
    os.system ( args )
    #time.sleep(5) 
    if maxp != 0 :
        args = 'class.com -i ' + datastar + ' -o maxp_particles.star -l list' + str(run_number) + ' -m '  + str(maxp)
        outpart = 'maxp_particles.star'
    else :
        args = 'class.com -i ' + datastar + ' -o particles' + str(run_number) + '.star -l list' + str(run_number)   
        outpart = 'particles' + str(run_number) + '.star' 
    os.system ( args )
    #time.sleep(5) 
    return outpart

def relion_3dclass ( diameter, nclasses, ignorectf, Tval, run_number, particles, model3d, mask ) :
    output_dir = quick3d_parameters['workdir'] + 'Class3D/'
    if os.path.isdir (output_dir) != True :
        os.makedirs(output_dir)
    os.chdir(output_dir)
    if mask != '' :
        mask = ' --solvent_mask ' + mask
    else :
        mask = ''
    args = 'mpirun -np 17 relion_refine_mpi --o ./run' + str(run_number) + ' --i ' + particles + ' --ref ' + model3d + ' --firstiter_cc --ini_high 40 --dont_combine_weights_via_disc --pool 100 --ctf  --pad 2 --fast_subsets  --iter ' + str(quick3d_parameters['relion_2d_iter']) + ' --tau2_fudge ' + str(Tval) + ' --K ' + str(nclasses)  \
           + ' --particle_diameter ' + str (diameter) + ' --flatten_solvent  --zero_mask  --oversampling 1 --healpix_order 2  --offset_range 5 --offset_step 2 --sym C1 --norm --scale  --j 4 --gpu ' + ignorectf + mask
    #print (args)
    output = open("relion3d.log", 'w')
    output.writelines ( str(args ) + '\n' )
    print ('\n3D Classification using Relion\n')
    proc = subp.Popen( args,  shell = True,  stdout = output)

    #Keep the user informed       
    name = 'run' + str(run_number) + '*model.star'
    a = len(glob.glob( name ))
    b = quick3d_parameters['relion_3d_iter'] + 1
    c = 0
    while ( a < b ) :
        if ( a > c ) :
            print ( a,'/',b-1, 'Iterations are done' )
            c = a
        a = len(glob.glob(name ))

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
    

if quick3d_parameters['listfile'] == '' :
    pattern = '*'+quick3d_parameters['micrograph_name_pattern'] +'*'+ quick3d_parameters['micrograph_name_suffix']
    micrographs_total, micrographs_list = count_micrographs(pattern, quick3d_parameters['datadir'], quick3d_parameters['micrograph_name_exclude'], 1  )
    if micrographs_total == 0 :
        quick3d_parameters['micrograph_name_suffix'] = 'tif'
        pattern = '*'+quick3d_parameters['micrograph_name_pattern'] +'*'+ quick3d_parameters['micrograph_name_suffix']
        micrographs_total, micrographs_list = count_micrographs(pattern, quick3d_parameters['datadir'], quick3d_parameters['micrograph_name_exclude'], 1  )
        if micrographs_total == 0 :
            print ('I cant find the micrographs in the directory')
            print ('I will quit')
            quit ()
    if micrographs_total == 0 :
        quick3d_parameters['micrograph_name_suffix'] = 'mrc'
        pattern = '*'+quick3d_parameters['micrograph_name_pattern'] +'*'+ quick3d_parameters['micrograph_name_suffix']
        micrographs_total, micrographs_list = count_micrographs(pattern, quick3d_parameters['datadir'], quick3d_parameters['micrograph_name_exclude'], 1  )
        if micrographs_total == 0 :
            print ('I cant find the micrographs in the directory')
            print ('I will quit')
            quit ()
else :
    micrographs_total, micrographs_list = count_micrographs_list(quick3d_parameters['listfile'] )

    
quick3d_parameters['z'] = get_movie_frame_number ( micrographs_list )
if quick3d_parameters['ask_pixel_size'] == 0 and quick3d_parameters['pixel_size'] == 0:
    quick3d_parameters['pixel_size'] = float ( get_pixel_size ( micrographs_list ) )
    if quick3d_parameters['pixel_size'] == 1 :
        print ( 'The pixel size from the header is 1.0, I hope thats the correct one')
        print (' If not please rerun the script with correct pixel size --pixel=1.63 for example')

quick3d_parameters['ngpu'] = number_of_gpu ( micrographs_list  )

quick3d_parameters['ncpu'] = number_of_cpu ()
if quick3d_parameters['ncpu'] > len ( micrographs_list) :
    quick3d_parameters['sparxcpu'] = len ( micrographs_list)
elif quick3d_parameters['qcpu'] != 0 :
    if int ( quick3d_parameters['qcpu'])  > len ( micrographs_list) :
        quick3d_parameters['sparxcpu'] = len ( micrographs_list)
    else :
        quick3d_parameters['sparxcpu'] = quick3d_parameters['qcpu']
else :
    quick3d_parameters['sparxcpu'] = quick3d_parameters['ncpu']

if quick3d_parameters['z'] == 1 :
    print ('\nYou have given micrographs rather than movies.. I wil start with CTF estimation')
    #make_dir_soft_links(quick3d_parameters['workdir'],quick3d_parameters['datadir'] , pattern, quick3d_parameters['micrograph_name_exclude']    )
    make_dir_soft_links(quick3d_parameters['workdir'], micrographs_list )
    raw_list = micrographs_list
    micrographs_total = len(glob.glob(os.path.join(quick3d_parameters['workdir'],'micrographs','*' )))
    micrographs_list = glob.glob(os.path.join(quick3d_parameters['workdir'],'micrographs','*' ))
    if quick3d_parameters['dwcombi'] == 1 :
        quick3d_parameters['dose'] = 1
        dw_list = [x.replace('.mrc', '_DW.mrc') for x in raw_list]
        make_dir_soft_links(quick3d_parameters['workdir'], dw_list )
        aligned_micrographs_suf = 'mrc'
        aligned_micrographs_dir = quick3d_parameters['workdir'] + '/micrographs/' 
else :
    make_dir_soft_links1 ( quick3d_parameters['workdir'], micrographs_list )
    #MotionCor2_align_script (quick3d_parameters['ngpu'], micrographs_list )
    MotionCor2_align_script (quick3d_parameters['ngpu'], quick3d_parameters['micrograph_name_suffix'], quick3d_parameters['gain'],  quick3d_parameters['rotgain'], micrographs_list )
    #MotionCor2_align (quick3d_parameters['ngpu']  )
    aligned_micrographs_dir = quick3d_parameters['workdir'] + '/aligned/' 
    aligned_micrographs_suf = 'mrc'
    pattern = '*'+quick3d_parameters['micrograph_name_pattern'] +'*'+ aligned_micrographs_suf
    micrographs_total, micrographs_list = count_micrographs( pattern , aligned_micrographs_dir, 'DW', 0 )
    if quick3d_parameters['bin'] != 0 :
        quick3d_parameters['pixel_size'] = float ( quick3d_parameters['pixel_size'] ) * float ( quick3d_parameters['bin'] )

micrographs_list_simple = micrographs_list_nodir_name(micrographs_list)

if quick3d_parameters['sparx'] == 1 :
    sparx_ctf( quick3d_parameters['pixel_size'], quick3d_parameters['cs'], quick3d_parameters['kv'], micrographs_list, quick3d_parameters['negative'], quick3d_parameters['sparxcpu'] )
else :
    GCTF( quick3d_parameters['pixel_size'], quick3d_parameters['cs'], quick3d_parameters['kv'], micrographs_list_simple, quick3d_parameters['negative'])

if quick3d_parameters['dose'] != 0 :
    pattern = '*'+quick3d_parameters['micrograph_name_pattern'] +'*DW*'+ aligned_micrographs_suf
    #print (pattern)
    if quick3d_parameters['dwcombi'] == 1 :
        quick3d_parameters['micrograph_name_exclude'] = '\#'
    micrographs_total, micrographs_list = count_micrographs( pattern , aligned_micrographs_dir, quick3d_parameters['micrograph_name_exclude'], 0 )
    micrographs_list_simple = micrographs_list_nodir_name(micrographs_list)

#if quick3d_parameters['template'] != '' :
#    template = ' --T ' + quick3d_parameters['datadir'] + '/' + quick3d_parameters['template']
if  quick3d_parameters['3dmodel'] != '' :
    make_template_from3D ( quick3d_parameters['datadir'] + '/' + quick3d_parameters['3dmodel'] )
    template = ' --T  model.mrcs --apixT ' + quick3d_parameters['3dmodelpix'] 
elif quick3d_parameters['template'] != '' :
    template = ' --T ' + quick3d_parameters['datadir'] + '/' + quick3d_parameters['template']  + ' --apixT ' + quick3d_parameters['templatepix']
elif quick3d_parameters['template'] == '' :
    template = ' '

Gautomatch ( micrographs_list_simple, quick3d_parameters['pixel_size'] , quick3d_parameters['diameter'], quick3d_parameters['cccutoff'], quick3d_parameters['pickcontrast'], template, quick3d_parameters['lavgmin'], quick3d_parameters['lavgmax'] )

if quick3d_parameters['dose'] != 0 :
    input_star  = open ("micrographs_ctf.star", 'r' )
    output_star = open ("micrographs_ctf_DW.star",  'w' )
    for line in input_star:
        output_star.write(line.replace('.mrc', '_DW.mrc')) 
    starfile = 'micrographs_ctf_DW.star'
    output_star.close()
    input_star.close()
else :
    starfile = 'micrographs_ctf.star'

if quick3d_parameters['box'] == 0 :
    quick3d_parameters['box'] =  round_up_to_even ( ( float ( quick3d_parameters['diameter'])   * 1.5 ) / float ( quick3d_parameters['pixel_size'] ) )
    quick3d_parameters['full_box'] = quick3d_parameters['box'] 

if quick3d_parameters['mask'] == 0 :
    quick3d_parameters['mask'] = float ( quick3d_parameters['diameter'] ) * 1.1

if quick3d_parameters['sparx'] == 1 :
    sparx_particles_extract ( quick3d_parameters['box'], quick3d_parameters['contrast'], quick3d_parameters['sparxcpu'] )
else :
    relion_particles_extract ( quick3d_parameters['box'], quick3d_parameters['contrast'], starfile, quick3d_parameters['diameter'], quick3d_parameters['pixel_size'], quick3d_parameters['partbin'] )

run_number = 1
particles = '../Particles/particles.star'

if quick3d_parameters['sparx'] != 1 :
    if quick3d_parameters['partbin'] != '' :
        quick3d_parameters['bin_pixel_size'] = round_up_to_even ( float ( quick3d_parameters['pixel_size'] ) * float ( quick3d_parameters['box'] ) / 96  )
        quick3d_parameters['box'] = 96

if quick3d_parameters['sparx'] == 1 :
    if quick3d_parameters['qcpu'] == 0 :
        sparx_2dclass ( quick3d_parameters['diameter'], particles , quick3d_parameters['ncpu'] ) 
    else :
        sparx_2dclass ( quick3d_parameters['diameter'], particles , quick3d_parameters['qcpu'] ) 
elif quick3d_parameters['relion_to_sparx'] == 1 :
    if quick3d_parameters['qcpu'] == 0 :
        relion_to_sparx_2dclass ( quick3d_parameters['diameter'], particles , quick3d_parameters['ncpu'] ) 
    else :
        relion_to_sparx_2dclass ( quick3d_parameters['diameter'], particles , quick3d_parameters['qcpu'] ) 
elif quick3d_parameters['3d'] == ''  :
    relion_2dclass ( quick3d_parameters['mask'],  quick3d_parameters['nclasses'], quick3d_parameters['ignorectf'], quick3d_parameters['Tval'], run_number, particles ) 


if quick3d_parameters['sparx'] == 1 and quick3d_parameters['qcpu'] != 0 :
    output_dir = quick3d_parameters['workdir'] + 'logs'
    os.chdir(output_dir)
    qhold_input = []
    qhold_input.append ( '#!/bin/bash' )
    qhold_input.append ( 'qsub7 -N job1 sxcter.submit' )
    qhold_input.append ( 'qsub7 -hold_jid job1 -N job2 sxwindow.submit' )
    qhold_input.append ( 'qsub7 -hold_jid job2 -N job3 e2bdb.submit' )
    qhold_input.append ( 'qsub7 -hold_jid job3 -N job4 sxisac.submit' )
    i = 0
    out = output_dir + '/qhold.sh'
    com = open (out, "w")
    while (i < len (qhold_input) ):
        com.write( str (qhold_input[i] ) )
        com.write("\n"  )
        i += 1
    st = os.stat(out)
    os.chmod(out, st.st_mode | stat.S_IEXEC)
    com.close()
    os.system('qhold.sh')

elif quick3d_parameters['relion_to_sparx'] == 1 and quick3d_parameters['qcpu'] != 0 :
    output_dir = quick3d_parameters['workdir'] + 'logs'
    os.chdir(output_dir)
    qhold_input = []
    qhold_input.append ( '#!/bin/bash' )
    qhold_input.append ( 'qsub7 -N job1 sxisac.submit' )
    i = 0
    out = output_dir + '/qhold.sh'
    com = open (out, "w")
    while (i < len (qhold_input) ):
        com.write( str (qhold_input[i] ) )
        com.write("\n"  )
        i += 1
    st = os.stat(out)
    os.chmod(out, st.st_mode | stat.S_IEXEC)
    com.close()
    os.system('qhold.sh')




#if quick3d_parameters['3dmodel'] != '' :
#    resize_3dmodel ( quick3d_parameters['3dmodelpix'], quick3d_parameters['pixel_size'] , quick3d_parameters['box'] , quick3d_parameters['datadir'] + '/' + quick3d_parameters['3dmodel'] ) 
    #relion_2dclass ( quick3d_parameters['mask'],  quick3d_parameters['nclasses'], quick3d_parameters['ignorectf'], quick3d_parameters['Tval'], run_number, particles ) 
#    relion_3dclass ( quick3d_parameters['diameter'],  '3' , quick3d_parameters['ignorectf'], '4', run_number, particles, 'resized.mrc' ) 
#else :
#    sparx_2dclass ( quick3d_parameters['diameter'], particles , quick3d_parameters['ncpu'] ) 

    

print (sg + '\nPrcoessing is done.. ')
print ('\nResults are in:    ' +  quick3d_parameters['timestamp'] + '\n' )
print (eb)


## RAJAN
class2d = 'no'
if quick3d_parameters['3d'] != ''  :
    class2d = 'done'
    user_input = '3d'

while ( class2d != 'done' ) :
    if quick3d_parameters['auto'] != 1 :
        os.chdir (quick3d_parameters['workdir']+'/Class2D')
        if os.path.exists ('backup_selection.star') == True :
            os.remove ('backup_selection.star')
        args = 'relion_display --class  --i  run' + str (run_number) + '_it0' + str(quick3d_parameters['relion_2d_iter']) +'_model.star --scale 1 --sort rlnClassDistribution --reverse --allow_save  \
               --fn_imgs ./classes' + str (run_number) + '.star  --fn_parts particles' + str(run_number) + '.star'
        proc = os.system ( args )

    #if os.path.exists ( 'particles1.star' ) != True :
    #    quit ()
    #classes = 'classes' + str(run_number) + '.star'
    #auto = 1
    #if auto != 1 :
        if os.path.exists ( 'particles1.star' ) != True :
            quit ()
        classes = 'classes' + str(run_number) + '.star'
        user_input = input ('\nContinue (2d/3d/no/reextract) : ')
        while user_input != '2d' and user_input != '3d' and user_input != 'no' and user_input != 'reextract':
            user_input = input ('Continue (2d/3d/no/reextract) : ')
        if user_input == '3d' or user_input == 'reextract' :
            class2d = 'done'
        elif user_input == 'no' :
            quit ()

        if user_input == '2d' :
            class_input = input ('\nClasses (  ' + str(quick3d_parameters['nclasses']) +  ' )  :')
            if class_input != '' :
                quick3d_parameters['nclasses'] = class_input
            particles = 'particles' + str (run_number) + '.star'
            run_number += 1
            relion_2dclass ( quick3d_parameters['mask'],  quick3d_parameters['nclasses'], quick3d_parameters['ignorectf'], quick3d_parameters['Tval'], run_number, particles ) 
    else :
        datastar = 'run' + str(run_number) + '_it0' + str (quick3d_parameters['relion_2d_iter'])  + '_data.star' 
        modelstar = 'run' + str(run_number) + '_it0' + str (quick3d_parameters['relion_2d_iter']) + '_model.star'
        particles = relion_auto_select ( modelstar, datastar, run_number, 0 )
        run_number += 1
        quick3d_parameters['nclasses'] = 15 + quick3d_parameters['nclasses']
        if run_number == quick3d_parameters['auto_2d_iter'] + 1 :
            user_input = '3dmodel' 
            particles = '../Class2D/particles' + str(quick3d_parameters['auto_2d_iter']) + '.star'
            particles_for_3d = '../Class2D/' + relion_auto_select ( modelstar, datastar, run_number-1, 400 )
            run_number = run_number - 1
            class2d = 'done'
        else :
            relion_2dclass ( quick3d_parameters['mask'],  quick3d_parameters['nclasses'], quick3d_parameters['ignorectf'], quick3d_parameters['Tval'], run_number, particles ) 

if user_input == '3dmodel' :
    relion_inital ( quick3d_parameters['mask'] , particles_for_3d )
    resize_3dmodel ( quick3d_parameters['bin_pixel_size'], quick3d_parameters['full_pixel_size'] , quick3d_parameters['full_box'] , quick3d_parameters['workdir'] + '/' + 'InitialModel/run1_it001_class001.mrc' )
    if quick3d_parameters['auto'] == 1 :
        user_input = 'reextract'

if user_input == 'reextract' :
    if  quick3d_parameters['z'] != 1 :
        os.chdir ( quick3d_parameters['workdir'] + '/' + 'aligned' )
    else :
        os.chdir ( quick3d_parameters['workdir'] + '/' + 'micrographs' )
    if quick3d_parameters['auto'] == 0 :
        user_input = input ('Give the box size for extraction: ')
        while user_input.replace('.','').isdigit() != True :
            print ('Please give box size as a number/float/decimal')
            user_input = input ('Give the box size for extraction: ')
        quick3d_parameters['box'] = user_input
    else :
        quick3d_parameters['box'] = quick3d_parameters['full_box']
        run_number = quick3d_parameters['auto_2d_iter']  

    relion_particles_reextract ( quick3d_parameters['box'], quick3d_parameters['contrast'], quick3d_parameters['diameter'], run_number, quick3d_parameters['dose']  )
    quick3d_parameters['reextract'] = 1

if quick3d_parameters['reextract'] == 1 :
    particles = '../reextract/particles.star'
    if quick3d_parameters['auto'] == 0 :
        user_input = input ('\nContinue (3d/no) : ')
        while user_input != '3d'  and user_input != 'no' :
            user_input = input ('Continue (3d/no) : ')
    else :
        user_input = '3d'
else  :
    if quick3d_parameters['auto'] == 0 :
        particles = '../Class2D/particles' + str(run_number) + '.star' 
    if quick3d_parameters['partbin'] != '' :
        quick3d_parameters['pixel_size'] = quick3d_parameters['bin_pixel_size']

if quick3d_parameters['3d'] != '' :
    particles = '../Particles/particles.star'
    if quick3d_parameters['partbin'] != '' :
        quick3d_parameters['pixel_size'] = quick3d_parameters['bin_pixel_size']


if user_input == '3d' :
    if quick3d_parameters['auto'] == 0 :
        if quick3d_parameters['3d'] != '' :
            quick3d_parameters['3dmodel'] = quick3d_parameters['3d']
            #quick3d_parameters['3dmodelpix'] = get_pixel_size ([ quick3d_parameters['datadir'] + '/' + quick3d_parameters['3dmodel']] )
        if quick3d_parameters['3dmodel'] != '' :
            if quick3d_parameters['3d'] != '' :
                quick3d_parameters['3dmodel'] = quick3d_parameters['datadir'] + '/' + quick3d_parameters['3dmodel']
            else :
                resize_3dmodel ( quick3d_parameters['3dmodelpix'], quick3d_parameters['pixel_size'] , quick3d_parameters['box'] , quick3d_parameters['datadir'] + '/' + quick3d_parameters['3dmodel'] ) 
                quick3d_parameters['3dmodel'] = 'resized.mrc'
        else :
            user_input = input ('Give the 3D model full path: ') 
            while os.path.exists (user_input) != True :
                print ('model file does not exist')
                user_input = input ('Give the 3D model full path: ') 
            quick3d_parameters['3dmodel'] = user_input
            print ( quick3d_parameters['3dmodel'] )
            user_input = input ('Give the 3D model pixel size: ') 
            while user_input.replace('.','').isdigit() != True :
                print ('Please give pixel size as a number/float/decimal')
                user_input = input ('Give the 3D model pixel size: ') 
            quick3d_parameters['3dmodelpix'] = user_input
            resize_3dmodel ( quick3d_parameters['3dmodelpix'], quick3d_parameters['pixel_size'] , quick3d_parameters['box'] , quick3d_parameters['3dmodel'] ) 
            quick3d_parameters['3dmodel'] = 'resized.mrc'
    else :
        quick3d_parameters['3dmodel'] = 'resized.mrc'
        quick3d_parameters['number_of_3dclasses'] = 1
    #3dclasses = '5'
    Tval = '4'
    run_number = '1'
    if quick3d_parameters['3dmask'] != '' :
        mask = quick3d_parameters['datadir'] + '/' + quick3d_parameters['3dmask']
    else :
        mask = ''
     
    relion_3dclass ( quick3d_parameters['diameter'], quick3d_parameters['number_of_3dclasses'] , quick3d_parameters['ignorectf'], Tval, run_number, particles, quick3d_parameters['3dmodel'] , mask )
