#!/bin/bash

setup_PYTHIA() {
    export PYTHIA8LOCATION=/u/at/pnef/Work/Code/pythia8183/
    export PYTHIA8DATA=${PYTHIA8LOCATION}xmldoc/
    export LD_LIBRARY_PATH=${PYTHIA8LOCATION}lib/:$LD_LIBRARY_PATH
}

setup_ROOT() {
    source /u/at/pnef/Work/Code/root_v5.34.17/bin/thisroot.sh
}

setup_fastjet() {
	export FASTJETLOCATION=/u/at/pnef/Work/Code/TrackBasedGrooming/stable/fastjet-3.0.3/fastjet-install/
	#export FASTJETLOCATION=/u/at/pnef/Work/Code/TrackBasedGrooming/fastjet-3.0.3/fastjet-install/
    #export FASTJETLOCATION=/u/at/pnef/Work/Code/fastjet-install/
    export LD_LIBRARY_PATH=${FASTJETPATH}lib/:$LD_LIBRARY_PATH
}

setup_boost() {
    export BOOSTINCDIR=/usr/include
    export BOOSTLIBLOCATION=/usr/lib64
}

setup_ROOT
setup_PYTHIA
setup_fastjet
setup_boost

