#!/usr/bin/env python3.5
import glob, sys, string, time
import os, fnmatch, re
import shutil, getpass 
import stat, platform, math
import multiprocessing
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import numpy.random as rnd
from matplotlib import cm as CM
import subprocess as subp


## Create a general dictionary of script parameters. Easy to expand later
author = 'prabu@biochem.mpg.de'
drift_parameters = {}
drift_parameters['datadir'] = os.getcwd()
drift_parameters['user'] = getpass.getuser()
drift_parameters['user_uid'] = os.getuid()
drift_parameters['kv'] = 300
drift_parameters['pixel_size'] = 0
drift_parameters['cs'] = 2.7
drift_parameters['micrograph_name_pattern'] = '*'
drift_parameters['micrograph_name_suffix'] = 'mrc'
drift_parameters['micrograph_name_exclude'] = '\#'
drift_parameters['magnification'] = 'N/A'
drift_parameters['time'] =  time.strftime("%d_%b_%Y_%I_%M_%S%p")
drift_parameters['timestamp'] =  'EM_' + drift_parameters['time']
drift_parameters['workdir'] = os.getcwd() + '/' + drift_parameters['timestamp'] + '/'
drift_parameters['Microscope'] = 'N/A'
drift_parameters['gctf'] = 1
drift_parameters['ctffind4'] = 0
drift_parameters['negative'] = 0
drift_parameters['gctf_available'] = 1
drift_parameters['ctffind4_available'] = 1
drift_parameters['plot_only'] = 0
drift_parameters['power_users'] = 0
drift_parameters['ngpu'] = 0
drift_parameters['ncpu'] = 0
drift_parameters['output_text'] = []
drift_parameters['input_commands'] = ' '.join(sys.argv)
drift_parameters['rescut'] = 0


### Command line and help text
### Argparse
from optparse import OptionParser

usage = "usage: %prog [options] \n\n" \
        "makesum.com --pix=1.42\n" \
        "makesum.com --pix=1.42 --cs=2.7 --kv=300 --suf=mrc --mag=55000 \n"\
        "makesum.com --pix=1.42 --cs=2.7 --kv=300 --suf=mrc --mag=55000 --cem=Halos --mic=FoilHole\n"\
        "makesum.com --pix=1.42 --cs=2.7 --kv=300 --suf=mrc --mag=55000 --cem=Halos --mic=_aligned\n\n"\
        "The script calculates the CTF [using GCTF or CTFFIND4] of input micrographs and prints a nice summary\n\
         and some useful plots"

parser = OptionParser(usage=usage)
parser.add_option("--pix",   dest="pixel_size", help="Pixel size of the Micrograph eg -p 1.42 or --pix=1.42")
parser.add_option("--kv",    dest="kv", help="Voltage of the Microscope  eg -k 300 or --kv=300 default 300")
parser.add_option("--cs",    dest="spherical_aberration", help="Spherical aberration of the Microscope  eg -c 2.7 or --cs=2.7 default 2.7")
parser.add_option("--mic",   dest="pattern", help="Common unique pattern in the micrograph name ( to identify/select)  eg -m _frames or --m=_frames default all mrc files in the directory")
parser.add_option("--exc",   dest="exclude_pattern", help="Micrograph name that you want to exclude    eg -e _DW or --exc=_DW no default ")
parser.add_option("--mag",   dest="magnification", help="Magnification of the microscope  eg -x 55000 or --mag=55000 No default -only used for the summary ")
parser.add_option("--suf",   dest="suffix", help="Micrograph suffix  eg mrc  default mrc ")
parser.add_option("--cem",        dest="microscope_name", help="Microscope name  eg Halos  no default")
parser.add_option("--ctffind",    dest="ctffind", action="store_true", default=False, help="Use ctffind4 instead of Gctf")
parser.add_option("--negative",   dest="negative", action="store_true", default=False, help="Negative stain data")
parser.add_option("--plot",       dest="plot_only", action="store_true", default=False, help="plot results of precomputed ctf Gctf or/and ctffind4")
parser.add_option("--cpu",        dest="ncpu", help="Number of cpus to use")
parser.add_option("--gpu",        dest="ngpu", help="Number of gpus to use")
parser.add_option("--devel",      dest="power", action="store_true", default=False, help="For developers....")
parser.add_option("--rcut",       dest="rescut", help="only works with --devel")
(options, args) = parser.parse_args()

if options.kv :
    drift_parameters['kv'] = options.kv

if options.pixel_size :
    drift_parameters['pixel_size'] = options.pixel_size 

if options.spherical_aberration :
    drift_parameters['cs'] = options.spherical_aberration

if options.pattern :
    drift_parameters['micrograph_name_pattern'] = options.pattern

if options.exclude_pattern :
    drift_parameters['micrograph_name_exclude'] = options.exclude_pattern

if options.magnification :
    drift_parameters['magnification'] = options.magnification

if options.suffix :
    drift_parameters['micrograph_name_suffix'] = options.suffix

if options.microscope_name :
    drift_parameters['Microscope'] = options.microscope_name

if options.ncpu :
    drift_parameters['ncpu'] = options.ncpu

if options.ngpu :
    drift_parameters['ngpu'] = options.ngpu

if options.rescut :
    drift_parameters['rescut'] = options.rescut

if str(options.ctffind) == "True" :
    drift_parameters['ctffind4'] = 1
    drift_parameters['gctf'] = 0
    #print ('\nUsing ctffind4 instead of Gctf\n')

if str(options.plot_only) == "True" :
    drift_parameters['plot_only'] = 1
    print ('\nplotting the results, no ctf calculations\n')

if str(options.power) == "True" :
    drift_parameters['power_users'] = 1

if str(options.negative) == "True" :
    drift_parameters['negative'] = 1

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

if drift_parameters['pixel_size'] == 0 :
    print ( "\nI need pixel size as a input, I kind of assume defaults for other values" )
    print ( "please provide pixel size For eg. --pix=1.43 or -p 1.43\n" )
    print ( "If you need help please use makesum.com -h" )
    quit ()

text_print_list = ( drift_parameters['pixel_size'] , drift_parameters['magnification'] , drift_parameters['Microscope'] )

## CHECK FOR THE PRESENSE OF EXECUTABLES
def check_execs() :

    exec_list = ['Gctf', 'ctffind415']

    for i in range (len(exec_list)):
        if not check_for_executables(exec_list[i]) :
            drift_parameters.setdefault('missing_execs', []).append(exec_list[i])

    if 'missing_execs' in drift_parameters :
        if 'Gctf' in drift_parameters['missing_execs'] :
            drift_parameters['gctf_available'] = 0
        if 'ctffind415' in drift_parameters['missing_execs'] :
            drift_parameters['ctffind4_available'] = 0
    if ( drift_parameters['ctffind4_available'] == 0 and drift_parameters['gctf_available'] == 0 ) :
        print ('\nNo  Gctf and ctffind415 executables in the path')
        print ('I will quit')
        quit ()
    elif ( drift_parameters['gctf'] == 1 and drift_parameters['gctf_available'] == 0 ) :
        if ( drift_parameters['ctffind4_available'] == 1 ):
            print ('\nNo Gctf executables in the path, but ctffind415 is available')
            drift_parameters['ctffind4'] = 1
            drift_parameters['gctf'] = 0
    elif ( drift_parameters['ctffind4'] == 1 and drift_parameters['ctffind4_available'] == 0 ) :
        if ( drift_parameters['gctf_available'] == 1 ):
            print ('ctffind requested.. its not in the path But Gctf available' )
            print ('\n I will calculate with GCTF\n')
            time.sleep (3)
            drift_parameters['gctf'] = 1
            drift_parameters['ctffind4'] = 0

if (  drift_parameters['gctf_available'] == 1 and drift_parameters['ctffind4_available'] == 1   and drift_parameters['power_users'] == 1) :
    drift_parameters['ctffind4'] = 1
    drift_parameters['gctf'] = 1

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

#check_execs()

if drift_parameters['plot_only'] == 1 :
    drift_parameters['gctf'] = 0
    drift_parameters['ctffind4'] = 0
else :
    check_execs()



def count_micrographs(pattern, micdir, exclude) :
    micrographs_list = glob.glob(os.path.join(micdir,pattern ))
 
    for name in glob.glob(os.path.join(micdir,pattern )) :
        if exclude in name :
            micrographs_list.remove(name)
    micrographs_total = len( micrographs_list )

    if micrographs_total == 0 :
        print ('I cant find the micrographs in the directory')
        print ('I will quit')
        quit ()
    else:
        print ('\nYou have',micrographs_total, 'micrographs')
    return micrographs_total,micrographs_list

def make_dir_soft_links ( workdir, micdir, pattern, exclude) :
    dest_dir  = workdir + '/'
    os.makedirs(workdir)
    for file in glob.glob (pattern) :
        if exclude not in file :
            file_with_dir = micdir + '/' + file
            file_with_dest_dir = dest_dir + '/' + file
            os.symlink ( file_with_dir, file_with_dest_dir )
    os.chdir ( dest_dir )

def count_micrographs_workdir(pattern, micdir) :
    micrographs_total = len(glob.glob(os.path.join(micdir,pattern )))
    micrographs_list = glob.glob(os.path.join(micdir,pattern ))
    suf = '.' + drift_parameters['micrograph_name_suffix']
    names_list = [w.replace(suf, '') for w in micrographs_list]
    return micrographs_total,micrographs_list,names_list

def number_of_gpu() :
    gpus = len(glob.glob('/proc/driver/nvidia/gpus/*'))
    drift_parameters['ngpu']  = int( drift_parameters['ngpu'] )
    if drift_parameters['ngpu'] == 0 :
        print (' Available GPUS: ',gpus )
    else :
        if ( gpus != 0 ) :
            print (' Available GPUS: ',gpus )
            print (' Using ', drift_parameters['ngpu'], 'GPUS' )
            gpus = drift_parameters['ngpu']
    return gpus

def number_of_cpu() :
    cpus = multiprocessing.cpu_count()
    drift_parameters['ncpu']  = int (drift_parameters['ncpu'] )
    if drift_parameters['ncpu'] == 0 :
        print (' Available CPUS: ',cpus )
    else :
        if ( cpus != 0 ) :
            cpus = drift_parameters['ncpu']
            print (' Using ', drift_parameters['ncpu'], 'CPUS' )
    return cpus

def GCTF ( pixel_size, cs, kv, micrographs_list, negative) :
    number_of_gpus = number_of_gpu()
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
    a = len(glob.glob("*.ctf" ))
    b = len (micrographs_list ) 
    c = 0
    print ('\n I am going to calculate CTF for your micrographs using Gctf\n')
    while ( a < b ) :
        if ( a > c ) :
            print ( a,'/',b, 'Micrographs are done' )
            c = a
        a = len(glob.glob("*.ctf" ))
       
    for p in ps:
        p.wait()

def write_ctffind4_input ( micname,pixel,cs, kv ) :
    name = ''.join ( micname.split('.')[:-1] )
    suf = micname.split('.')[-1]
    ctf_input_list = [ pixel, kv, cs, '0.07','512','20','5','5000','50000','500','no','no','yes','100','no','no','EOF']
    ctf_input = []
    ctf_input.append ( '#!/usr/bin/csh' )
    ctf_input.append ('ctffind-4.1.5.exe > ' + name + '.log << EOF' )
    ctf_input.append ( name + '.' + suf )
    ctf_input.append ( name + '.ctf' )
    ctf_input =  ctf_input + ctf_input_list
    i = 0
    out = name + '.com'
    com = open (out, "w")
    while (i < 21 ):
        com.write( str (ctf_input[i] ) ) 
        com.write("\n"  ) 
        i += 1
    st = os.stat(out)
    os.chmod(out, st.st_mode | stat.S_IEXEC)
    com.close()

def CTFFIND4 ( pixel_size, cs, kv, micrographs_list) :
    number_of_cpus = number_of_cpu()

    the_dot_with_suff = '.' + drift_parameters['micrograph_name_suffix']
    for mics in micrographs_list :
       micname =     mics.split('/')[-1]
       write_ctffind4_input( micname,pixel_size ,cs, kv   )
    com_list = glob.glob("*.com" )
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
    output = open("test.log", 'w')
    for i in range(number_of_cpus) :
        args = 'cpu'+ str (i)
        proc = subp.Popen( args,  shell = True,  stdout = output)
        time.sleep(0.1)
        ps.append(proc)

    #Keep the user informed
    a = len(glob.glob("*.ctf" ))
    b = len (micrographs_list )
    c = 0
    print ('\n I am going to calculate CTF for your micrographs using ctffind4\n')
    while ( a < b ) :
        if ( a > c ) :
            print ( a,'/',b, 'Micrographs are done' )
            c = a
        a = len(glob.glob("*.ctf" ))

    for p in ps:
        p.wait()

        
def grep(pattern,fileObj):
    r=[]
    linenumber=0
    for line in fileObj:
        linenumber +=1
        if re.search(pattern,line):
            r.append((line))
    return r

def tail(file):
    with open(file, 'r') as o:
        return o.readlines()[-1]

def CTFFIND4_results_list () :
    defocus1_list = []
    defocus2_list = []
    defocus_angle_list = []
    CCC_list = []
    resolution_list = []

    for file in glob.glob("*.txt") :
        if 'avrot'not  in file :
            b  =  tail(file )
            defocus1_list.append(b.split()[1])
            defocus2_list.append(b.split()[2])
            defocus_angle_list.append(b.split()[3])
            CCC_list.append(b.split()[5])
            resolution_list.append(b.split()[6])
    for i in  range (len (resolution_list)) :
        if ( resolution_list[i] == 'inf' ) :
            resolution_list[i] = 50 
    defocus1_list = list(map(float, defocus1_list))
    defocus2_list = list(map(float, defocus2_list))
    defocus_angle_list = list(map(float, defocus_angle_list))
    CCC_list = list(map(float, CCC_list))
    resolution_list = list(map(float, resolution_list))
    return defocus1_list, defocus2_list, defocus_angle_list, CCC_list, resolution_list

def GCTF_results_list () :
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
    return defocus1_list, defocus2_list, defocus_angle_list, CCC_list, resolution_list

def plot_only_gctf_micrographs_list () :
    suf = '.' + drift_parameters['micrograph_name_suffix']
    gctf_list = glob.glob("*_gctf.log")
    micrographs_list = [w.replace('_gctf.log', suf) for w in gctf_list]
    micrographs_total = len (micrographs_list)
    return micrographs_total, micrographs_list

def plot_only_ctffind3_micrographs_list () :
    suf = '.' + drift_parameters['micrograph_name_suffix']
    gctf_list = glob.glob("*_ctffind3.log")
    micrographs_list = [w.replace('_ctffind3.log', suf) for w in gctf_list]
    micrographs_total = len (micrographs_list)
    return micrographs_total, micrographs_list

def plot_only_ctffind4_micrographs_list () :
    suf = '.' + drift_parameters['micrograph_name_suffix']
    gctf_list = glob.glob("*_ctffind4.log")
    micrographs_list = [w.replace('_ctffind4.log', suf) for w in gctf_list]
    micrographs_total = len (micrographs_list)
    return micrographs_total, micrographs_list

def power_calculation () :
    id_list = [ 3336, 491 ]
    extra_calculation = 0
    for id in id_list :
        if id == drift_parameters['user_uid'] :
            a = 1
    if a == 1 :
        extra_calculation = 1
    else :
        print ('Not implemented' )
        quit ()
    return extra_calculation

def plot_and_format_results (defocus1_list, defocus2_list, defocus_angle_list, CCC_list, resolution_list, text_print_list, time, gctf_flag ) :
    mail_author = 'echo " " | mail -s MAKESUM ' + author 
    os.system( mail_author )
    print ('\nCTF done\nI am plotting the results')
    mean_defocus1 = round ( sum (defocus1_list) / len (defocus1_list), 1 )
    mean_defocus2 = round ( sum (defocus2_list) / len (defocus2_list), 1 )
    mean_defocus_angle = round ( sum (defocus_angle_list) / len (defocus_angle_list), 1 )
    mean_CCC = round ( sum (CCC_list) / len (CCC_list), 1 )
    mean_resolution = round ( sum (resolution_list) / len (resolution_list),1 )
    csfont = {'fontname':'Sans-serif'}

    plt.style.use('ggplot')
    title_font_size = 12
    summary_font_size = 10
    ax = plt.axes([.65, .6, .4, .4])
    plt.hist(resolution_list, range(2,30), color='royalblue',alpha=0.75)
    #plt.hist(resolution_list,  color='royalblue',alpha=0.75)
    plt.title('Estimated Resolution', size=title_font_size, y = 1.05 )
    ax.set_xlabel( r'$Resolution ( {\AA})$')
    ax.set_ylabel(r'$No.\ of\ Micrographs$')
    #plt.tick_params( axis='x', which='both', bottom='on', top='off', labelbottom='on') 
    #plt.tick_params( axis='y', which='both', right='off', left='on', labelleft='on') 
    #ax.xaxis.set_tick_params(labeltop='on')
    #ax.xaxis.set_tick_params(labelbottom='off')

    ax1 = plt.axes([1.2, .6, .4, .4])
    plt.hist(CCC_list, color = 'crimson', alpha=0.75)
    plt.title('CTF Score', size=title_font_size, y = 1.04 )
    ax1.set_xlabel(r'$CCC$')
    ax1.set_ylabel(r'$No.\ of\ Micrographs$')

    with plt.style.context(('ggplot')):
        defocus1_list_in_micrometer = [ x / 10000 for x in defocus1_list ]
        defocus2_list_in_micrometer = [ x / 10000 for x in defocus2_list ]
        ax3 = plt.axes([0.1, 0.02, .4, .4])
        ax3.set_xlabel (  r'$Defocus-x ({\mu}m)$' )
        ax3.set_ylabel (  r'$Defocus-y ({\mu}m)$' )

        plt.hist2d( defocus1_list_in_micrometer, defocus2_list_in_micrometer, bins = 40, cmap=CM.Blues )
        plt.colorbar()
        plt.title('Defocus Spread', size=title_font_size, y = 1.05 )
        plt.grid(b=False)

    with plt.style.context(('ggplot')):
        ax = plt.axes([.1, 0.6, .4, .4 ])
        y = list(reversed(np.arange(0,1,0.15)))
        #text = [time.strftime("%c"),\
        text = ['Number of  Micrographs : '+ str(len(micrographs_list)),\
               'Mean Resolution estimate : '+str(mean_resolution)+  ' '+  r'$ \AA $',\
               'Mean Defocus estimate : '+str(mean_defocus1)+' '+  r'$ \AA $', \
               'Pixel size : '+ text_print_list[0] +  r'$\ \AA $', \
               'Magnification : '+ text_print_list[1],\
               'Microscope : '+ text_print_list[2]   ]
        for i in range(6):
            plt.text(0.05, y[i], text[i], alpha=1,clip_on=True,fontsize=summary_font_size)
            plt.xticks(())
            plt.yticks(())
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            #ax.spines['left'].set_visible(False)
            ax.set_axis_bgcolor([ 0.95, 0.95, 0.95, 1.0] )
            title = 'Summary: ' + time
            plt.title(title,size=title_font_size, y = 1.05 )

    from matplotlib.patches import Ellipse
    
    with plt.style.context(('ggplot')):
        ells= []
        xlimit = max(defocus1_list)/10000 + 0.5
        ylimit = max(defocus2_list)/10000 + 0.5
        for i, j, k in zip (defocus1_list,defocus2_list,defocus_angle_list ) :
            ells.append ( Ellipse( xy = [ xlimit/2 , ylimit/2 ] , width = i / 10000 , height = j / 10000 , angle = k, linewidth = 0.4 ) )
        ax = plt.axes([.7 ,  0.02, .3, .4])
        for e in ells:
            ax.add_artist(e)
            e.set_clip_box(ax.bbox)
            e.set_alpha(0.1)
            e.set_facecolor('none' )
            e.set_edgecolor( [0.4, 0.4, 0.4 ])
           
        ax.grid(b=False)
        ax.set_xlim(0, xlimit)
        ax.set_ylim(0, ylimit)
        ax.set_axis_bgcolor([ 0.9, 0.9, 0.9, 1.0] )    
        ax.set_xlabel (  r'$Defocus-x ({\mu}m)$' )
        ax.set_ylabel (  r'$Defocus-y ({\mu}m)$' )
        plt.title('Beam Astigmatism',size=title_font_size, y = 1.05 )
    
    if ( gctf_flag == 1 ) :
        ax = plt.axes([1.2 , 0.02, .4, .4])
        scores_twenty_to_eight,scores_fifteen_to_six,scores_twelve_to_five,scores_ten_to_four,scores_eight_to_three = ([] for i in range(5))
        bins_list = [ '20-08A', '15-06A', '12-05A','10-04A','08-03A' ]
        if ( float ( drift_parameters['pixel_size'] ) * 2  ) < 3 :
            bins_list = [ '20-08A', '15-06A', '12-05A','10-04A','08-03A' ]
        elif ( float ( drift_parameters['pixel_size'] ) * 2  ) < 4 :
            bins_list = [ '20-08A', '15-06A', '12-05A','10-04A' ]
            scores_eight_to_three.append(1)
        elif ( float ( drift_parameters['pixel_size'] ) * 2  ) < 5 :
            bins_list = [ '20-08A', '15-06A', '12-05A', ]
            scores_eight_to_three.append(1)
            scores_ten_to_four.append(1) 
        elif ( float ( drift_parameters['pixel_size'] ) * 2  ) < 6 :
            bins_list = [ '20-08A', '15-06A' ]
            scores_eight_to_three.append(1)
            scores_ten_to_four.append(1) 
            scores_twelve_to_five.append(1)

        for file in glob.glob("*_gctf.log") :
            for bin in bins_list :
                myfile = open (file, "r" )
                if   ( bin == '20-08A' ):
                    a  =  ' '.join(grep(bin,myfile ))
                    scores_twenty_to_eight.append ((a.split()[-1]))
                elif ( bin == '15-06A' ):
                    a  =  ' '.join(grep(bin,myfile ))
                    scores_fifteen_to_six.append ((a.split()[-1]))
                elif ( bin == '12-05A' ):
                    a  =  ' '.join(grep(bin,myfile ))
                    scores_twelve_to_five.append ((a.split()[-1]))
                elif ( bin == '10-04A' ):
                    a  =  ' '.join(grep(bin,myfile ))
                    scores_ten_to_four.append ((a.split()[-1]))
                elif ( bin == '08-03A' ):
                    a  =  ' '.join(grep(bin,myfile ))
                    scores_eight_to_three.append ((a.split()[-1]))
        myfile.close()
        scores_twenty_to_eight = (np.array(list(map(int, scores_twenty_to_eight ))))
        scores_fifteen_to_six  = (np.array(list(map(int, scores_fifteen_to_six ))))
        scores_twelve_to_five  = (np.array(list(map(int, scores_twelve_to_five ))))
        scores_ten_to_four     = (np.array(list(map(int, scores_ten_to_four ))))
        scores_eight_to_three  = (np.array(list(map(int, scores_eight_to_three ))))

        scores_list = ( scores_twenty_to_eight,scores_fifteen_to_six,scores_twelve_to_five,scores_ten_to_four,scores_eight_to_three )
        scores = []
        for name in scores_list :
            for i in range(1,6):
                scores.append(np.sum ( name == i ))
        j = 1
        x = []
        y = []
        for i in scores_list :
            x.extend ([j] * len (i)  )
            y.extend (i)
            j += 1
        plt.hist2d ( x, y, bins = 5, cmap=CM.Blues )
        plt.grid(b=False)

        x = [ 1.4, 2.2, 2.9 ,3.7, 4.5, 5.5 ]
        y = [ 1.3, 2.1, 2.9, 3.7, 4.5, 5.3 ]
        k=0 
        for i in range(5):
            for j in range(5):
                plt.text(x[i], y[j], scores[k], alpha=1,clip_on=True,fontsize=8)  
                k += 1

        x = [ 1.4, 2.3, 3.0 ,3.8, 4.6 ]
        y = [ 1, 2, 3, 4, 5 ]
        bin_names = [ '20-8', '15-6', '12-5', '10-4', '8-3' ]
        plt.xticks(x, bin_names  )
        plt.yticks(x, y  )
        ax.set_xlabel(r'$Resolution\ Bins\ ({\AA})$')
        ax.set_ylabel(r'$CTF\ \  Fitting\ \  Score$')
        #ax.set_xticks(np.arange(x[0])+0.5, minor=False)
        #ax.xaxis.set_tick_params(labeltop='on')
        #ax.xaxis.set_tick_params(labelbottom='off')
        plt.title('CTF Validation', size = title_font_size, y = 1.05)
        #plt.axes().set_aspect('equal')
        plt.margins(0.2)
        plt.subplots_adjust(bottom=0.15)

    ax3 = plt.axes([0.1, -0.2 , 1.2, .05])
    ax3.text(0.01, 0.6, drift_parameters['input_commands'], alpha=1,clip_on=True,fontsize=10) 
    plt.xticks(())
    plt.yticks(())
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax3.set_axis_bgcolor([ 1,1,1, 1.0] )

    plt.savefig('try.pdf',bbox_inches='tight')
    plt.gcf().clear()
    drift_parameters['output_text'] = [ len(micrographs_list) , mean_resolution, mean_defocus1, text_print_list[0] , text_print_list[1], text_print_list[2] ]

def make_output_pdf ( filename ) :
    output_file_name =  '../' + filename +'.pdf'
    if os.path.isfile("try.pdf"):
        shutil.move("try.pdf", output_file_name)
        #print ('\nCaculations are done.. Your output PDF is :',filename+'.pdf\n')
        pdf = filename+'.pdf'
    else :
        print ('something went wrong.. Inform Rajan')
    return pdf

def make_output_pdf_plot_only ( filename ) :
    output_file_name =  filename +'.pdf'
    if os.path.isfile("try.pdf"):
        shutil.move("try.pdf", output_file_name)
        #print ('\nCaculations are done.. Your output PDF is :',filename+'.pdf\n')
        pdf = filename+'.pdf'
    else :
        print ('something went wrong.. Inform Rajan')
    return pdf

def make_clean ( datadir, workdir, miclist ) :
    os.chdir ( datadir )
    remove_dir = datadir + '/' + workdir 
    #os.system( rmdir )
    for link in miclist :
        os.unlink(link) 
    tempdir = datadir +'/'+workdir+'/*'
    #print ('I am cleaning up the temp files' )
    for file in glob.glob(tempdir) :
        os.remove(file)
    #time.sleep(3)
    shutil.rmtree( remove_dir, ignore_errors=True )


## Call modules and get the job done
## check if all ok
extra_calculations = power_calculation()

if  drift_parameters['gctf'] == 1 :
    pattern = '*'+drift_parameters['micrograph_name_pattern'] +'*'+ drift_parameters['micrograph_name_suffix']
    micrographs_total, micrographs_list = count_micrographs(pattern, drift_parameters['datadir'], drift_parameters['micrograph_name_exclude']  )
    make_dir_soft_links(drift_parameters['workdir'],drift_parameters['datadir'] , pattern, drift_parameters['micrograph_name_exclude']    )
    micrographs_total, micrographs_list,mics_names_list = count_micrographs_workdir(pattern, drift_parameters['workdir'])
    GCTF( drift_parameters['pixel_size'], drift_parameters['cs'], drift_parameters['kv'], micrographs_list, drift_parameters['negative'])
    defocus1_list, defocus2_list, defocus_angle_list, CCC_list, resolution_list = GCTF_results_list()
    plot_and_format_results(defocus1_list, defocus2_list, defocus_angle_list, CCC_list, resolution_list , text_print_list, drift_parameters['time'] , 1 )
    pdf = make_output_pdf ( drift_parameters['timestamp']+'_GCTF'  )

if  drift_parameters['ctffind4'] == 1 :
    pattern = '*'+drift_parameters['micrograph_name_pattern'] +'*'+ drift_parameters['micrograph_name_suffix']
    micrographs_total, micrographs_list = count_micrographs(pattern, drift_parameters['datadir'], drift_parameters['micrograph_name_exclude']  )
    if drift_parameters['power_users'] == 0 :
        make_dir_soft_links(drift_parameters['workdir'],drift_parameters['datadir'] , pattern, drift_parameters['micrograph_name_exclude']    )
    micrographs_total, micrographs_list,mics_names_list = count_micrographs_workdir(pattern, drift_parameters['workdir'])
    CTFFIND4( drift_parameters['pixel_size'], drift_parameters['cs'], drift_parameters['kv'], micrographs_list)
    defocus1_list, defocus2_list, defocus_angle_list, CCC_list, resolution_list = CTFFIND4_results_list()
    plot_and_format_results(defocus1_list, defocus2_list, defocus_angle_list, CCC_list, resolution_list , text_print_list, drift_parameters['time'],  0 )
    pdf = make_output_pdf ( drift_parameters['timestamp']+'_CTFFIND4'  )


if  drift_parameters['plot_only'] == 1 :
    gctf_logs     = len(glob.glob("*_gctf.log"))
    ctffind3_logs = len(glob.glob("*_ctffind3.log"))
    ctffind4_logs = len(glob.glob("*_ctffind4.log"))
    if ( gctf_logs == 0 and ctffind3_logs == 0 and ctffind4_logs == 0 ) :
        print ('There are no files with *_gctf.log or *_ctffind4.log ')
        print ('Please run --plot where you have precomputed ctf log files are')
        quit () 

    if gctf_logs > 0 :
        micrographs_total, micrographs_list = plot_only_gctf_micrographs_list ()
        defocus1_list, defocus2_list, defocus_angle_list, CCC_list, resolution_list = GCTF_results_list()
        plot_and_format_results(defocus1_list, defocus2_list, defocus_angle_list, CCC_list, resolution_list , text_print_list, drift_parameters['time'] , 0 )
        pdf = make_output_pdf_plot_only ( drift_parameters['timestamp'] + '-GCTF' )
   
    if ctffind3_logs > 0 :
        micrographs_total, micrographs_list = plot_only_ctffind3_micrographs_list ()
        defocus1_list, defocus2_list, defocus_angle_list, CCC_list, resolution_list = CTFFIND4_results_list()
        plot_and_format_results(defocus1_list, defocus2_list, defocus_angle_list, CCC_list, resolution_list , text_print_list, drift_parameters['time'] , 0 )
        pdf = make_output_pdf_plot_only ( drift_parameters['timestamp'] + '-CTFFIND3' )
   
    if ctffind4_logs > 0 :
        micrographs_total, micrographs_list = plot_only_ctffind3_micrographs_list ()
        defocus1_list, defocus2_list, defocus_angle_list, CCC_list, resolution_list = CTFFIND4_results_list()
        plot_and_format_results(defocus1_list, defocus2_list, defocus_angle_list, CCC_list, resolution_list , text_print_list, drift_parameters['time'] , 0 )
        pdf = make_output_pdf_plot_only ( drift_parameters['timestamp'] + '_CTFFIND4' )

if drift_parameters['power_users'] == 1 and extra_calculations == 1:
    power_calculation()
    micrographs_total, micrographs_list,mics_names_list = count_micrographs_workdir(pattern, drift_parameters['workdir'])
    suf = '.' + drift_parameters['micrograph_name_suffix']
    names_list = [w.replace(suf, '') for w in micrographs_list]
    res_dic = {}
    for name in names_list :
        name = name.split('/')[-1]
        ctf_log  = '.txt'
        gctf_log = '_gctf.log'
        filec = name + ctf_log
        fileg = name + gctf_log
        myfile = open (fileg, "r" )
        g  =  float ( ( ' '.join(grep("RES_LIMIT",myfile )) ).split()[-1] )
        file = open (filec, "r" )
        b  =  tail(filec )
        c = float (b.split()[6])
        res_dic[name] = g,c
        myfile.close()

    good = 'micrographs_matching.star'
    bad  = 'micrographs_differ.star'
    g = open (good, "w")
    b = open (bad,  "w")
    g.writelines ( '\n' + 'data_' + '\n\n' + 'loop_' + '\n' + '_rlnMicrographName #1' + '\n')
    b.writelines ( '\n' + 'data_' + '\n\n' + 'loop_' + '\n' + '_rlnMicrographName #1' + '\n')
    for i in res_dic :
        if  abs ( res_dic [i][0] - res_dic [i][1] )  > 3   :
            b.writelines( i+suf+ '\n')
        else :
            if ( float (drift_parameters['rescut'] ) == 0 ) :
                g.write(i+suf+ '\n')    
            else :
                if res_dic [i][0] < float ( drift_parameters['rescut'] ):
                    g.write(i+suf+ '\n')    
    g.close()
    b.close()
    

## PRINT A NICE SUMMARY 
print ('\n') 
print ('######################################################################')
print ('Command Input              :', drift_parameters['input_commands'] )
print ('Output PDF file            :', pdf )
print ('Number of  Micrographs     :', drift_parameters['output_text'] [0] )
print ('Mean Resolution estimate   :', drift_parameters['output_text'] [1],'A' )
print ('Mean Defocus estimate      :', drift_parameters['output_text'] [2],'A' )
print ('Pixel size                 :', drift_parameters['output_text'] [3],'A' )
print ('Magnification              :', drift_parameters['output_text'] [4] )
print ('Microscope                 :', drift_parameters['output_text'] [5] )
if drift_parameters['power_users'] == 1 and extra_calculations == 1:
    print ('Matching/similar ctf       : micrographs_ctf_matching.star')
    shutil.move("micrographs_matching.star", "../micrographs_ctf_matching.star")
    print ('Differences in ctf         : micrographs_ctf_differ.star')
    shutil.move("micrographs_differ.star", "../micrographs_ctf_differ.star")
print ('######################################################################')
## Output

#make_output_pdf ( drift_parameters['timestamp'] )

## CLEAN UP

make_clean ( drift_parameters['datadir'] , drift_parameters['timestamp'], micrographs_list )
