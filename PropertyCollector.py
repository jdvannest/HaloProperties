import os,argparse

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
args = parser.parse_args()

#Path to directory where datafiles are written
output_path = '/myhome2/users/vannest/Data/'
#Path to prefered python executable
python_path = '/myhome2/users/vannest/anaconda3/bin/python'

for sim in ['cptmarvel','elektra','storm','rogue','h148','h229','h242','h329']:
    o_flag = '-o' if sim=='cptmarvel' and args.overwrite else ''
    os.system(f'{python_path} Marvel_DCJL.{args.property}.py -s {sim} '+
              f'-p {output_path} -n {args.numproc} {o_flag}')