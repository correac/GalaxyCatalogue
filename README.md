GalaxyCatalogue
=========

A python package that makes visualisation plots of stars and dark matter from the TangoSIDM simulation suite.

Requirements
------------

The GalaxyCatalogue package requires:

+ `python3.6` or above
+ see requirements.txt


Installing
----------

To get started using the package you need to set up a python virtual environment. The steps are as follows:

Clone GalaxyCatalogue
```
git clone https://github.com/correac/GalaxyCatalogue.git

cd GalaxyCatalogue

python3 -m venv galaxycatalogue_env
```

Now activate the virtual environment.

```
source galaxycatalogue_env/bin/activate
```

Update pip just in case
```
pip install pip --upgrade

pip install -r requirements.txt
```

How to use it
-------------

To run the script type
```bash
 python3 visualization.py -d run_directory \
                          -s snapshot_name \
                          -c catalogue_name \
                          -n name_of_the_run \
			  -t simulation_type \
                          -o path_to_output_directory 
```

