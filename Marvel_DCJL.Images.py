import pynbody,pickle,sys,pymp,warnings,argparse,os
import numpy as np

def myprint(string,clear=False):
    if clear:
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K") 
    print(string)
warnings.filterwarnings("ignore")
parser = argparse.ArgumentParser()
parser.add_argument("-s","--simulation",choices=['cptmarvel','elektra','storm','rogue','h148','h229','h242','h329'],required=True)
parser.add_argument("-o","--overwrite",action='store_true')
parser.add_argument("-p","--path",required=True)
parser.add_argument("-n","--numproc",required=True,type=int)
args = parser.parse_args()


#########################################################################
#Change the name of the property you're collecting (used for file naming)
image_dir = args.path+f'Images/'
#########################################################################


#Load in Simulation Info (paths and halo numbers)
Sims = pickle.load(open('SimulationInfo.pickle','rb'))
#Load in the existing data file, or create new/overwrite if necessary
#ile_name = args.path+f'Marvel_DCJL.{property}.pickle'
#if args.overwrite:
#    os.system('rm '+file_name)
#    Data = {}
#    for sim in Sims:
#        Data[sim] = {}
#    print('Writing new Data File...')
#else:
#    try:
#        Data = pickle.load(open(file_name,'rb'))
#        print('Data File loaded.')
#    except:
#        print('No Data File found. Writing new one...')
#        Data = {}
#        for sim in Sims:
#            Data[sim] = {}


sim  = args.simulation
#Skip if this simulation if it has already been done (unless overwriting)
if False:
    print(f'{sim} already completed')
else:
    print(f'Loading {sim}...')
    #Create shared dictionary for multiprocessing data and shared 1-element
    # array to act as scalar for traking progress
    #SimData = pymp.shared.dict()
    prog=pymp.shared.array((1,),dtype=int)
    #Load in simulation
    s = pynbody.load(Sims[sim]['path'])
    s.physical_units()
    h = s.halos()
    myprint(f'{sim} Loaded.',clear=True)
    print('\tWriting: 0%')
    #Begin parallelized analysis
    with pymp.Parallel(args.numproc) as pl:
        for i in pl.xrange(len(Sims[sim]['halos'])):
            #Load in specific halo
            halonum = Sims[sim]['halos'][i]
            halo = h[halonum]

            pynbody.analysis.angmom.faceon(halo)
            R = pynbody.analysis.luminosity.half_light_r(halo)
            pynbody.plot.stars.render(halo,width=8*R,filename=image_dir+f'{sim}.{halonum}.Faceon.png')
            halo.rotate_x(90)
            pynbody.plot.stars.render(halo,width=8*R,filename=image_dir+f'{sim}.{halonum}.Sideon.png')   

            #Print progress to terminal
            myprint(f'\tWriting: {round(float(prog[0]+1)/float(len(Sims[sim]["halos"]))*100,2)}%',clear=True)
            
            #with pl.lock:  <----- the 'lock' block probably isn't necessary since the only thing that might
            #                    encounter a race condition is the progress counter, which is only for the UI
            #Insert the current halo's data into the shared dictionary
            
            prog[0]+=1


    print(f'{sim} Finished.')