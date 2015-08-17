'''
processing.py
author: Luke de Oliveira, July 2015 

Simple utilities for processing the junk that comes out of the ntuple event generation.
'''

import numpy as np
from .jettools import rotate_jet, flip_jet, plot_mean_jet

def buffer_to_jet(entry, tag = 0, side = 'r', max_entry = 2000, pix = 25):
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
        * Tau32
        * Tau21
        * Tau32old
        * Tau21old
    """

    image = flip_jet(rotate_jet(np.array(entry['Intensity']), -entry['RotationAngle'], normalizer=max_entry, dim=pix), side)
    e_norm = np.linalg.norm(image)
    return (image / e_norm, tag, 
        entry['LeadingPt'], entry['LeadingEta'], 
        entry['LeadingPhi'], entry['LeadingM'], 
        entry['Tau32'], entry['Tau21'], 
        entry['Tau32old'], entry['Tau21old'], e_norm)


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
