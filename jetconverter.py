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
import array




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

    parser.add_argument('--dump', 
                        default='dump.root',
                        help = 'ROOT file to dump all this into (writes to TTree `images`)')

    parser.add_argument('--save', 
                        default='output.npy', 
                        help = 'Filename to write out the data.')
    parser.add_argument('--plot',  
                        help = 'File prefix that\
                         will be part of plotting filenames.')
    parser.add_argument('--ptmin', default=200.0, help = 'minimum pt to consider')
    parser.add_argument('--ptmax', default=400.0, help = 'maximum pt to consider')

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
        from rootpy.tree import Tree, TreeModel, FloatCol, FloatArrayCol
    except ImportError:
        raise ImportError('rootpy (www.rootpy.org) not installed\
         -- install, then try again!')

    class JetImage(TreeModel):
        '''
        Buffer for Jet Image
        '''
        image = FloatArrayCol(25 ** 2)
        signal = FloatCol()
        jet_pt = FloatCol()
        jet_eta = FloatCol()
        jet_phi = FloatCol()
        jet_delta_R = FloatCol()
        tau_32 = FloatCol()
        tau_21 = FloatCol()
        tau_1 = FloatCol()
        tau_2 = FloatCol()
        tau_3 = FloatCol()
                
            

    pix_per_side = -999
    entries = []
    with root_open(args.dump, "recreate") as ROOTfile:
        tree = Tree('images', model=JetImage)
        
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
                            )
                        )
                    if (np.abs(jet['LeadingEta']) < 2) & (jet['LeadingPt'] > float(args.ptmin)) & (jet['LeadingPt'] < float(args.ptmax)):
                        buf = buffer_to_jet(jet, tag, max_entry=100000, pix=pix_per_side)
                        
                        tree.image = buf[0].ravel().astype('float32')
                        tree.signal = buf[1]
                        tree.jet_pt = buf[2]
                        tree.jet_eta = buf[3]
                        tree.jet_phi = buf[4]
                        tree.jet_delta_R = buf[5]
                        tree.tau_32 = buf[6]
                        tree.tau_21 = buf[7]
                        tree.tau_1 = buf[8]
                        tree.tau_2 = buf[9]
                        tree.tau_3 = buf[10]
                        entries.append(buf)
                        tree.fill()
        tree.write()


    # -- datatypes for outputted file.
    _bufdtype = [('image', 'float32', (pix_per_side, pix_per_side)), 
                 ('signal', 'float32'),
                 ('jet_pt', 'float32'),
                 ('jet_eta', 'float32'),
                 ('jet_phi', 'float32'), 
                 ('jet_mass', 'float32'),
                 ('jet_delta_R', 'float32'),
                 ('tau_32', 'float32'), 
                 ('tau_21', 'float32'),
                 ('tau_1', 'float32'),
                 ('tau_2', 'float32'),
                 ('tau_3', 'float32')]

    df = np.array(entries, dtype=_bufdtype)
    logger.info('saving to file: {}'.format(savefile))
    np.save(savefile, df)

    if plt_prefix != '':
        logger.info('plotting...')
        plot_mean_jet(df[df['signal'] == 0], title="Average Jet Image, Background").savefig(plt_prefix + '_bkg.pdf')
        plot_mean_jet(df[df['signal'] == 1], title="Average Jet Image, Signal").savefig(plt_prefix + '_signal.pdf')
