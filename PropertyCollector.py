import os,argparse

parser = argparse.ArgumentParser()
parser.add_argument("-p","--property",required=True)
parser.add_argument("-o","--overwrite",action='store_true')
args = parser.parse_args()

for sim in ['cptmarvel','elektra','storm','rogue','h148','h229','h242','h329']:
    if sim=='cptmarvel' and args.overwrite:
        os.system(f'/myhome2/users/vannest/anaconda3/bin/python /home/vannest/Code/Marvel_DCJL.{args.property}.py -s {sim} -o')
    else:
        os.system(f'/myhome2/users/vannest/anaconda3/bin/python /home/vannest/Code/Marvel_DCJL.{args.property}.py -s {sim}')