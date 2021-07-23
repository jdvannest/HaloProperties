import pynbody,pickle,sys,pymp,warnings,argparse,os
import numpy as np
#from modules.Custom import halo_shape_stellar
def myprint(string,clear=False):
    if clear:
        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K") 
    print(string)
warnings.filterwarnings("ignore")
parser = argparse.ArgumentParser()
parser.add_argument("-s","--simulation",choices=['cptmarvel','elektra','storm','rogue','h148','h229','h242','h329'],required=True)
parser.add_argument("-o","--overwrite",action='store_true')
args = parser.parse_args()


marvel = '/myhome2/users/munshi/dwarf_volumes/'
dc = '/myhome2/users/munshi/e12gals/'
Sims = {
    'cptmarvel' : {
        'path' : marvel+'cptmarvel.cosmo25cmb.4096g5HbwK1BH/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096/cptmarvel.cosmo25cmb.4096g5HbwK1BH.004096',
        'halos' : [1,2,3,5,6,7,10,11,13,14,24]
    },
    'elektra' : {
        'path' : marvel+'elektra.cosmo25cmb.4096g5HbwK1BH/elektra.cosmo25cmb.4096g5HbwK1BH.004096/elektra.cosmo25cmb.4096g5HbwK1BH.004096',
        'halos' : [1,2,3,4,5,8,9,10,11,12,17,36,64]
    },
    'storm' : {
        'path' : marvel+'storm.cosmo25cmb.4096g5HbwK1BH/storm.cosmo25cmb.4096g5HbwK1BH.004096/storm.cosmo25cmb.4096g5HbwK1BH.004096',
        'halos' : [1,2,3,4,5,6,7,8,10,11,12,14,15,22,23,31,37,44,48,55,118]
    },
    'rogue' : {
        'path' : marvel+'rogue.cosmo25cmb.4096g5HbwK1BH/rogue.cosmo25cmb.4096g5HbwK1BH.004096/rogue.cosmo25cmb.4096g5HbwK1BH.004096',
        'halos' : [1,3,7,8,10,11,12,15,16,17,28,31,37,58,116]
    },
    'h148' : {
        'path' : dc+'h148.cosmo50PLK.3072g3HbwK1BH/h148.cosmo50PLK.3072g3HbwK1BH.004096/h148.cosmo50PLK.3072g3HbwK1BH.004096',
        'halos' : [ 2, 3, 4, 6, 7, 11, 12, 13, 15, 20, 23, 27, 28, 29, 33, 34, 37, 38, 41, 43, 51, 59, 65, 75, 86, 94, 109, 114, 122]
    },
    'h229' : {
        'path' : dc+'h229.cosmo50PLK.3072gst5HbwK1BH/h229.cosmo50PLK.3072gst5HbwK1BH.004096/h229.cosmo50PLK.3072gst5HbwK1BH.004096',
        'halos' : [ 2, 3, 6, 14, 15, 18, 20, 22, 25, 33, 47, 48, 49, 52, 57, 62, 89, 92, 127]
    },
    'h242' : {
        'path' : dc+'h242.cosmo50PLK.3072gst5HbwK1BH/h242.cosmo50PLK.3072gst5HbwK1BH.004096/h242.cosmo50PLK.3072gst5HbwK1BH.004096',
        'halos' : [8, 10, 21, 26, 30, 34, 38, 42, 44, 45, 63, 70, 81, 138]
    },
    'h329' : {
        'path' : dc+'h329.cosmo50PLK.3072gst5HbwK1BH/h329.cosmo50PLK.3072gst5HbwK1BH.004096/h329.cosmo50PLK.3072gst5HbwK1BH.004096',
        'halos' : [7, 29, 30, 37, 53, 92, 115, 117, 127]
    }
}


#Change the name of the property you're collecting (used for file naming)
property = 'DarkMatterShape'


file_name = f'/myhome2/users/vannest/Data/Marvel_DCJL.{property}.pickle'
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
    numproc = 9
    with pymp.Parallel(numproc) as pl:
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








            