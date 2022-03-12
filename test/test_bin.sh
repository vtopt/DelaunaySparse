#!/bin/bash

# Run delsparses on 2d/4d VarSys test problems and analyze output
bin/delsparses data/varsys/sample_input2d.dat > sample_out2d.txt
if [[ `wc -l < sample_out2d.txt` == 710 ]]
then
	echo The command-line executables seem to be installed correctly.
	rm sample_out2d.txt
else
	echo There seems to be an issue with the CL install of delaunaysparses.
	echo See sample_out2d.txt for more information...
    exit 1
fi
bin/delsparses data/varsys/sample_input4d.dat > sample_out4d.txt
if [[ `wc -l < sample_out4d.txt` == 3027 ]]
then
	echo The command-line executables seem to be installed correctly.
	rm sample_out4d.txt
else
	echo There seems to be an issue with the CL install of delaunaysparses.
	echo See sample_out4d.txt for more information...
    exit 1
fi

# Run delsparsep on 2d/4d VarSys test problems and analyze output
export OMP_NUM_THREADS=2
bin/delsparsep data/varsys/sample_input2d.dat > sample_out2d.txt
if [[ `wc -l < sample_out2d.txt` == 710 ]]
then
	echo The command-line executables seem to be installed correctly.
	rm sample_out2d.txt
else
	echo There seems to be an issue with the CL install of delaunaysparsep.
	echo See sample_out2d.txt for more information...
    exit 1
fi
bin/delsparsep data/varsys/sample_input4d.dat > sample_out4d.txt
if [[ `wc -l < sample_out4d.txt` == 3027 ]]
then
	echo The command-line executables seem to be installed correctly.
	rm sample_out4d.txt
else
	echo There seems to be an issue with the CL install of delaunaysparsep.
	echo See sample_out4d.txt for more information...
    exit 1
fi
