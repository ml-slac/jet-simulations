'''
file: jetconverter.py
author: Luke de Oliveira, Mar. 2015 

This file takes files (*.root) produced by the event-gen
portion of the jet-images codebase and converts them into 
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



'''
import sys
import argparse
import numpy as np
from jettools import rotate_jet, flip_jet, plot_mean_jet


def buffer_to_jet(entry, tag = 0, side = 'r', max_entry = 2000, pix = 25):
    """
    Takes an entry from an ndarray, and a tag = {0, 1} indicating
    if its a signal entry or not. The parameter 'side' indicates 
    which side of the final jet image we want the highest energy.
    """

    image = flip_jet(rotate_jet(np.array(entry['Intensity']), -entry['RotationAngle'], normalizer=max_entry, dim=pix), side)
    return (image / np.linalg.norm(image), tag, 
        entry['LeadingPt'], entry['LeadingEta'], entry['LeadingPhi'], entry['LeadingM'], entry['Tau32'], entry['Tau21'])


def is_signal(f, matcher = 'Wprime'):
    """
    Takes as input a filename and a string to match. If the 
    'matcher' string is found in the filename, the file is 
    taken to be a signal file.
    """
    key = matcher.lower().replace(' ', '').replace('-', '')
    if key in f.lower().replace(' ', '').replace('-', ''):
        return 1.0
    return 0.0

def perfectsquare(n):
    return n % n**0.5 == 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
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

    # parser.add_argument('--size', 
    #                     default='output.npy', 
    #                     help = 'Filename to write out the data.')

    parser.add_argument('--plot',  
                        help = 'File prefix that will be part of plotting filenames.')

    parser.add_argument('files', nargs='*', help='Files to pass in')

    args = parser.parse_args()

    if len(args.files) < 1:
        sys.stderr.write('Must pass at least one file in.\n')
        exit(1)

    signal_match = args.signal
    files = args.files
    savefile = args.save
    plt_prefix = ''
    if args.plot:
        plt_prefix = args.plot


    from rootpy.io import root_open

    pix_per_side = -999
    entries = []
    for fname in files:
        if args.verbose:
            print 'Working on file: {}'.format(fname)
        with root_open(fname) as f:
            df = f.EventTree.to_array()

            pix = df[0]['Intensity'].shape[0]

            if not perfectsquare(pix):
                raise ValueError('shape of image array must be square.')

            if (pix_per_side > 1) and (int(np.sqrt(pix)) != pix_per_side):
                raise ValueError('all files must have same sized images.')
            
            pix_per_side = int(np.sqrt(pix))

            tag = is_signal(fname, signal_match)
            for jet in df:
                entries.append(buffer_to_jet(jet, tag, max_entry=2600, pix=pix_per_side))


    # datatypes for outputted file.
    _bufdtype = [('image', 'float64', (pix_per_side, pix_per_side)), 
             ('signal', float),
             ('jet_pt', float),
             ('jet_eta', float),
             ('jet_phi', float), 
             ('jet_mass', float),
             ('tau_32', float), 
             ('tau_21', float)]

    df = np.array(entries, dtype=_bufdtype)
    np.save(savefile, df)

    if plt_prefix != '':
        plot_mean_jet(df[df['signal'] == 0], title="Average Jet Image, Background").savefig(plt_prefix + '_bkg.pdf')
        plot_mean_jet(df[df['signal'] == 1], title="Average Jet Image, Signal").savefig(plt_prefix + '_signal.pdf')
