#!/bin/bash

python3 visualization.py \
	-d /Users/camila/SimulationData/mahti/L025N376/Hydro/SigmaConstant00 \
	-s snapshot_0036.hdf5 \
	-c subhalo_0036.properties \
	-n RefModel2SigmaConstant00 \
	-t Hydro \
	-o ./output_data


