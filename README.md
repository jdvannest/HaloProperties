# HaloProperties

This is a set of parallelized scripts designed to collect properties from all of the resolved Marvel+DCJL dwarfs.

## Setup
### Config.py
This script initializes a dictionary file with all of the necessary simulation information, namely the paths to the files and the list of resolved halos. It also intializes a config dictionary file that is used by the other scripts.

The following variables should be edited depending on which machine the scripts are being run on

**output\_path**: the path to the desired directory for the data files (default is the DataFiles directory within the repository)

**python\_path**: the absolute path to your preferred python executable

**marvel\_path**: the absolute path to the Marvel suite

**dcjl\_path**: the absolute path to the DCJL suite

Config.py only needs to be run after a path is updated.

### Running in Parallel
In order to run in parallel, the python package [pymp](https://github.com/classner/pymp.git) is required. This can be installed with pip:
	
	pip install pymp-pypi
	
Further documentation can be found [here](https://github.com/classner/pymp.git).

## Usage
### Property Scripts already in the GitHub
To collect a property for which there is already a script (i.e. if there is a file called "Marvel\_DCJL.<*property*>.py"), then simply run the "PropertyCollector" script:

	python PropertyCollector.py -p Masses -n 5 -o

The PropertyCollector script takes the following flags:

**-p/--property**: The desired property to be collected (must match syntax of the associated Marvel_DCJL.<*property*>.py file)

**-n/--numproc** (optional): The number of processes to use (default = 1)

**-o/--overwrite** (optional): Add to overwrite an existing file for this property, otherwise already completed sims will be skipped (this way, a crash won't require a complete restart)

### New Property Scripts
To collect a new property, a "Marvel\_DCJL.<*property*>.py" must first be generated. To do this, simply copy and rename the "Marvel\_DCJL.TEMPLATE.py" script and edit the following two code blocks:

	#########################################################################
	#Change the name of the property you're collecting (used for file naming)
	property = 'TEST'
	#########################################################################

The above variable 'property' should be updated to the desired name for the property(ies) you want to collect, and the names of their associated files.

	########################################################
	#Replace this block with the desired property to collect
	current['test_property_1'] = 1
	current['test_property_2'] = [3,6,9]
	current['test_property_3'] = {'a':1,'b':[1,2,3]}
	########################################################

The above block should be replaced with the actual property calculations/gathering. The desired values should be written as keys to the 'current' dictionary object, and the halo object is already loaded into the variable 'halo'. For example, if you wanted to gather the stellar mass, you would do:

	current['Mstar'] = halo.s['mass'].sum()
	
The keys that you add here will be the keys in the resultant dictionary file.

Once the "Marvel\_DCJL.<*property*>.py" has been created, you can simply run the PropertyCollector script as above.

## Output
Data from the PropertyCollector script are written out as dictionary files named "Marvel\_DCJL.<*property*>.pickle" to the output directory specified in "Config.py". 

The dictionaries are structured as follows:

	DataFile[simulation][halo_number][property]
	
So in the stellar mass example above, it would access the data by:

	in [1]: import pickle
	
	in [2]: Data = pickle.load(open(<path_to_datafile>,'rb'))
	
	in [3]: Data['cptmarvel']['1']['Mstar']
	out[3]: 54344890.82587707
	