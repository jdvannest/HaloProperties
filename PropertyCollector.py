import os,argparse,pickle

parser = argparse.ArgumentParser(description='This script iterates through'+
            'the Marvel+DCJL resolved dwarf halos and collects the desired'+
            'property. The data is written out as a dictionary in a .pickle'+
            'file, and is read as DataFile[simulation][halo][property]')
parser.add_argument("-p","--property",required=True,
                    help='the desired property to be collected')
parser.add_argument("-o","--overwrite",action='store_true',
                    help='Add flag to overwrite existing datafile')
parser.add_argument("-n","--numproc",type=int,default=1,
                    help='The number of processes to use (default=1)')
parser.add_argument("-b","--bubble",action='store_true',
                    help='Run Storm Superbubble instead of Marvel/DCJL')
args = parser.parse_args()


config = pickle.load(open('Config.pickle','rb'))

if args.bubble:
    o_flag = '-o' if args.overwrite else ''
    os.system(f'{config["python_path"]} PropertyScripts/Marvel_DCJL.{args.property}.py -s storm_bubble '+
            f'-p {config["output_path"]} -n {args.numproc} {o_flag}')

else:
    for sim in ['cptmarvel','elektra','storm','rogue','h148','h229','h242','h329']:
        o_flag = '-o' if sim=='cptmarvel' and args.overwrite else ''
        os.system(f'{config["python_path"]} PropertyScripts/Marvel_DCJL.{args.property}.py -s {sim} '+
                f'-p {config["output_path"]} -n {args.numproc} {o_flag}')