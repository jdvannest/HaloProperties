import pynbody,pickle,sys,pymp,warnings,argparse,os
import numpy as np

def myprint(string,clear=False):
    if clear:
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K") 
    print(string)
warnings.filterwarnings("ignore")
parser = argparse.ArgumentParser()
parser.add_argument("-s","--simulation",choices=['cptmarvel','elektra','storm','rogue','h148','h229','h242','h329','storm_bubble'],required=True)
parser.add_argument("-o","--overwrite",action='store_true')
parser.add_argument("-p","--path",required=True)
parser.add_argument("-n","--numproc",required=True,type=int)
args = parser.parse_args()


#########################################################################
#Change the name of the property you're collecting (used for file naming)
property = 'DarkMatterProfiles'
image_dir = args.path+f'Images/DarkMatterProfiles/'
def CoreEinasto(r,rho_s,r_s,r_c,a):
    return( rho_s*np.exp(-(2/a)*(((r+r_c)/r_s)**a-1)) )
from pynbody.analysis.profile import Profile
from scipy.optimize import curve_fit
from scipy import polyfit
import matplotlib.pylab as plt

#########################################################################


#Load in Simulation Info (paths and halo numbers)
Sims = pickle.load(open('SimulationInfo.pickle','rb'))
#Load in the existing data file, or create new/overwrite if necessary
file_name = args.path+f'Marvel_DCJL.{property}.pickle'
if args.overwrite:
    os.system('rm '+file_name)
    Data = {}
    for sim in Sims:
        Data[sim] = {}
    print('Writing new Data File...')
else:
    try:
        Data = pickle.load(open(file_name,'rb'))
        print('Data File loaded.')
    except:
        print('No Data File found. Writing new one...')
        Data = {}
        for sim in Sims:
            Data[sim] = {}


sim  = args.simulation
#Skip if this simulation if it has already been done (unless overwriting)
if len(Data[sim])>0:
    print(f'{sim} already completed')
else:
    print(f'Loading {sim}...')
    #Create shared dictionary for multiprocessing data and shared 1-element
    # array to act as scalar for traking progress
    SimData = pymp.shared.dict()
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

            #The shared dictionary object does not allow writing to sub-dictionaries
            #So use 'current' as a 'middle-man' step that will be inserted as a whole
            #into the shared dictionary at the end
            current={}

            ########################################################
            #Replace this block with the desired property to collect

            f,ax = plt.subplots(1,1)
            current['rbins'] = np.nan
            current['density'] = np.nan
            current['fitpar'] = np.nan
            current['CoreSlope'] = np.nan
            try:
                pynbody.analysis.angmom.faceon(halo)
                prof = Profile(halo.d,ndim=3,min=0.25,max=10,type='log')
                current['rbins'] = prof['rbins']
                current['density'] = prof['density']
                x,y = np.log10(prof['rbins']),np.log10(prof['density'])
                ax.plot(x,y,c='k')

                par,err = curve_fit(CoreEinasto,x,y)
                current['fitpar'] = par
                ax.plot(x,CoreEinasto(x,par[0],par[1],par[2],par[3]),c='b')
                ax.axvline(np.log10(par[2]),color='k',linestyle='--')

                x_in = np.log10(prof['rbins'][(prof['rbins']<par[2])])
                y_in = np.log10(prof['density'][(prof['rbins']<par[2])])
                line =  polyfit(x_in,y_in,1)
                current['CoreSlope'] = line[0]
                ax.plot(x_in,line[0]*x_in+line[1],c='r')
            except:
                error = 1

            ax.set_xlabel(r'Log$_{10}$(r)')
            ax.set_ylabel(r'Log$_{10}$($\rho_{dm}$)')
            f.savefig(image_dir+f'{sim}.{halonum}.png',bbox_inches='tight',pad_inches=.1)
            plt.close()

            ########################################################
            

            #Print progress to terminal
            myprint(f'\tWriting: {round(float(prog[0]+1)/float(len(Sims[sim]["halos"]))*100,2)}%',clear=True)
            
            #with pl.lock:  <----- the 'lock' block probably isn't necessary since the only thing that might
            #                    encounter a race condition is the progress counter, which is only for the UI
            #Insert the current halo's data into the shared dictionary
            SimData[str(halonum)] = current
            prog[0]+=1
            del current

    #Transfer the data from the shared dictionary to the main Data File (saving the shared
    #dictionary just results in garbage data for some reason)
    for halo in Sims[sim]['halos']:
        Data[sim][str(halo)] = SimData[str(halo)]

    #Update the Data File
    out = open(file_name,'wb')
    pickle.dump(Data,out)
    out.close()
    print(f'File updated with {sim}.')
