'''
processing.py
author: Luke de Oliveira, July 2015 

Simple utilities for processing the junk that comes out of the ntuple event generation.
'''

import numpy as np
from .jettools import rotate_jet, flip_jet, plot_mean_jet

def buffer_to_jet(entry, tag = 0, side = 'r', max_entry = None, pix = 25):
    """
    Takes an *single* entry from an structured ndarray, i.e., X[i], 
    and a tag = {0, 1} indicating if its a signal entry or not. 
    The parameter 'side' indicates which side of the final 
    jet image we want the highest energy.

    The `entry` must have the following fields (as produced by event-gen)
        * Intensity
        * RotationAngle
        * LeadingPt
        * LeadingEta
        * LeadingPhi
        * LeadingM
        * DeltaR
        * Tau32
        * Tau21
        * Tau{n} for n = 1, 2, 3 
    """

    image = flip_jet(rotate_jet(np.array(entry['Intensity']), -entry['RotationAngle'], normalizer=4000.0, dim=pix), side)
    e_norm = np.linalg.norm(image)
    return ((image / e_norm).astype('float32'), np.float32(tag), 
        np.float32(entry['LeadingPt']), np.float32(entry['LeadingEta']), 
        np.float32(entry['LeadingPhi']), np.float32(entry['LeadingM']), np.float32(entry['DeltaR']),
        np.float32(entry['Tau32']), np.float32(entry['Tau21']), np.float32(entry['Tau1']), np.float32(entry['Tau2']), np.float32(entry['Tau3']))


def is_signal(f, matcher = 'wprime'):
    """
    Takes as input a filename and a string to match. If the 
    'matcher' string is found in the filename, the file is 
    taken to be a signal file.
    """
    key = matcher.lower().replace(' ', '').replace('-', '')
    if key in f.lower().replace(' ', '').replace('-', ''):
        return 1.0
    return 0.0
