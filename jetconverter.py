#!/usr/bin/env python
'''
file: jetconverter.py
author: Luke de Oliveira, July 2015 

This file takes files (*.root) produced by the event-gen
portion of the jet-simulations codebase and converts them into 
a more usable format. In particular, this produces a 
numpy record array with the following fields:
    
    * 'image'  : (25, 25) numpy arrays that have been
                   rotated using cubic spline interpolation
                   to have the leading and subleading 
                   subjets aligned. In additon, the images
                   are pooled/flipped so one side of the image 
                   holds the most energy.
    
    * 'signal' : {0, 1} for whether or not the sample 
                   matches signal or not. Invoke:
                   python jetconverter.py -h for more help

    * 'jet_{x}': x is {pt, eta, phi}, and is the value of x
                   for the leading subjet.

    * 'tau_{NM}': N-subjettiness.

'''


from argparse import ArgumentParser
import sys
import logging

import numpy as np

from jettools import plot_mean_jet, buffer_to_jet, is_signal


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def perfectsquare(n):
    '''
    I hope this is self explanatory...
    '''
    return n % n**0.5 == 0

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--verbose', 
                        action='store_true', 
                        help='Verbose output')

    parser.add_argument('--signal', 
                        default='Wprime',
                        help = 'String to search for in\
                         filenames to indicate a signal file')

    parser.add_argument('--save', 
                        default='output.npy', 
                        help = 'Filename to write out the data.')
    parser.add_argument('--plot',  
                        help = 'File prefix that\
                         will be part of plotting filenames.')

    parser.add_argument('files', nargs='*', help='Files to pass in')

    args = parser.parse_args()

    if len(args.files) < 1:
        logger.error('Must pass at least one file in -- terminating with error.')
        exit(1)

    signal_match = args.signal
    files = args.files
    savefile = args.save
    plt_prefix = ''
    if args.plot:
        plt_prefix = args.plot

    try:
        from rootpy.io import root_open
    except ImportError:
        raise ImportError('rootpy (www.rootpy.org) not installed\
         -- install, then try again!')

    pix_per_side = -999
    entries = []
    for fname in files:
        logger.info('working on file: {}'.format(fname))
        with root_open(fname) as f:
            df = f.EventTree.to_array()

            n_entries = df.shape[0]

            pix = df[0]['Intensity'].shape[0]

            if not perfectsquare(pix):
                raise ValueError('shape of image array must be square.')

            if (pix_per_side > 1) and (int(np.sqrt(pix)) != pix_per_side):
                raise ValueError('all files must have same sized images.')
            
            pix_per_side = int(np.sqrt(pix))

            tag = is_signal(fname, signal_match)
            for jet_nb, jet in enumerate(df):
                if jet_nb % 1000 == 0:
                    logger.info('processing jet {} of {} for file {}'.format(
                            jet_nb, n_entries, fname
                        ))
                entries.append(
                    buffer_to_jet(jet, tag, max_entry=2600, pix=pix_per_side)
                    )


    # -- datatypes for outputted file.
    _bufdtype = [('image', 'float64', (pix_per_side, pix_per_side)), 
                 ('signal', float),
                 ('jet_pt', float),
                 ('jet_eta', float),
                 ('jet_phi', float), 
                 ('jet_mass', float),
                 ('tau_32', float), 
                 ('tau_21', float)]

    df = np.array(entries, dtype=_bufdtype)
    logger.info('saving to file: {}'.format(savefile))
    np.save(savefile, df)

    if plt_prefix != '':
        logger.info('plotting...')
        plot_mean_jet(df[df['signal'] == 0], title="Average Jet Image, Background").savefig(plt_prefix + '_bkg.pdf')
        plot_mean_jet(df[df['signal'] == 1], title="Average Jet Image, Signal").savefig(plt_prefix + '_signal.pdf')
