from rootpy.tree import TreeChain
from rootpy.io import root_open
from jettools import rotate_jet, flip_jet
import numpy as np

import sys

_bufdtype = [('image', 'float64', (25, 25)), 
		     ('signal', float),
		     ('jet_pt', float),
		     ('jet_eta', float),
		     ('jet_phi', float)]

def buffer_to_jet(entry, tag = 0):
	"""
	Takes an entry from an ndarray, and a tag = {0, 1} indicating
	if its a signal entry or not.
	"""
	image = flip_jet(rotate_jet(np.array(entry['Intensity']), -entry['RotationAngle']))
	return (image, tag, entry['LeadingPt'], entry['LeadingEta'], entry['LeadingPhi'])


def is_signal(f, matcher = 'Wprime'):
	key = matcher.lower().replace(' ', '').replace('-', '')
	if key in f.lower().replace(' ', '').replace('-', ''):
		return 1.0
	return 0.0

if __name__ == '__main__':
	
	entries = []
	for fname in sys.argv[1:]:
		print 'Working on file: {}'.format(fname)
		with root_open(fname) as f:
			df = f.EventTree.to_array()
			tag = is_signal(fname, 'wprime')
			for jet in df:
				entries.append(buffer_to_jet(jet, tag))

	df = np.array(entries, dtype=_bufdtype)

	np.save('test.npy', x)
