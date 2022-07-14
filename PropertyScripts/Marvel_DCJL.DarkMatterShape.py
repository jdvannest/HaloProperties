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


Sims = pickle.load(open('SimulationInfo.pickle','rb'))


#Change the name of the property you're collecting (used for file naming)
property = 'DarkMatterShape'


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
if len(Data[sim])>0:
    print(f'{sim} already completed')
else:
    print(f'Loading {sim}...')
    SimData = pymp.shared.dict()
    s = pynbody.load(Sims[sim]['path'])
    s.physical_units()
    h = s.halos()
    myprint(f'{sim} Loaded.',clear=True)
    print('\tWriting: 0%')
    prog=pymp.shared.array((1,),dtype=int)
    with pymp.Parallel(args.numproc) as pl:
        for i in pl.xrange(len(Sims[sim]['halos'])):
            halonum = Sims[sim]['halos'][i]
            current={}
            halo = h[halonum]

            #Replace this block with the desired property to collect
            cont=True
            for i in [1,2,3]:
                current[f'b_{i}Reff'] = np.nan
                current[f'c_{i}Reff'] = np.nan
            current['b_Edge'] = np.nan
            current['c_Edge'] = np.nan
            current['Reff'] = np.nan
            current['b'] = np.nan
            current['c'] = np.nan
            current['rbins'] = np.nan
            try:
                pynbody.analysis.angmom.faceon(halo)
                Rhalf = pynbody.analysis.luminosity.half_light_r(halo,band='v')
                current['Reff'] = Rhalf
            except:
                cont=False
            if cont:
                try:
                    r,ba,ca,angle,Es = pynbody.analysis.halo.halo_shape(halo)
                    for i in [1,2,3]:
                        ind = np.where(np.abs(r-i*Rhalf)==min(np.abs(r-i*Rhalf)))[0][0]
                        current[f'b_{i}Reff'] = ba[ind]
                        current[f'c_{i}Reff'] = ca[ind]
                    current['b_Edge'] = ba[-1]
                    current['c_Edge'] = ca[-1]
                    current['b'] = ba
                    current['c'] = ca
                    current['rbins'] = r
                except:
                    err=1

            
            myprint(f'\tWriting: {round(float(prog[0]+1)/float(len(Sims[sim]["halos"]))*100,2)}%',clear=True)
            #with pl.lock:
            SimData[str(halonum)] = current
            prog[0]+=1
            del current


    for halo in Sims[sim]['halos']:
        Data[sim][str(halo)] = SimData[str(halo)]

    out = open(file_name,'wb')
    pickle.dump(Data,out)
    out.close()
    print(f'File updated with {sim}.')








            