import os,argparse,pickle

parser = argparse.ArgumentParser(description='This script iterates through'+
            'the Marvel+DCJL resolved dwarf halos and collects the desired'+
            'property. The data is written out as a dictionary in a .pickle'+
            'file, and is read as DataFile[simulation][halo][property]')
parser.add_argument("-p","--property",required=True,
                    help='the desired property to be collected')

args = parser.parse_args()

Data = pickle.load(open(f'DataFiles/Marvel_DCJL.{args.property}.pickle','rb'))
Marvel = {}
DCJL = {}
for sim in ['cptmarvel','elektra','storm','rogue']:
    Marvel[sim] = Data[sim]
out = open(f'DataFiles/Marvel.{args.property}.pickle','wb')
pickle.dump(Marvel,out)
out.close()
for sim in ['h148','h229','h242','h329']:
    DCJL[sim] = Data[sim]
out = open(f'DataFiles/DCJL.{args.property}.pickle','wb')
pickle.dump(DCJL,out)
out.close()

os.system(f'rm DataFiles/Marvel_DCJL.{args.property}.pickle')