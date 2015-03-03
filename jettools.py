import numpy as np
import skimage.transform.rotate as im_rotate

def rotate_jet(jet, angle, in_radians = True, normalizer = 1450, dim=25):
	"""
	Take an *flat* unrotated in the form of an array from MI.exe, and an angle, 
	and rotates the jet to that angle using a
	"""
	im = jet.reshape((dim, dim))
	im = np.flipud(im.T) / normalizer

	if in_radians:
		angle = np.rad2deg(angle)

	return im_rotate(im, angle, order = 3)


def flip_jet(jet, pool = 'r'):
	"""
	Takes a rotated jet (25, 25) from rotate_jet() and flips it across the
	vertical axis according to whether you want the right side
	(pool = {r, R, right, Right, ...}) or the 
	left side (pool = {l, L, Left, left, ...}) to contain the most energy.
	"""
	weight = jet.sum(axis = 0)

	halfway = jet.shape[0] / 2.
	l, r = np.int(np.floor(halfway)), np.int(np.ceil(halfway))
	l_weight, r_weight = np.sum(weight[:l]), np.sum(weight[r:])

	if ('r' in pool.lower()) and ('l' in pool.lower()):
		raise ValueError('Jet pooling side must have l -OR- r in the name.')
	if 'r' in pool.lower():
		if r_weight > l_weight:
			return jet
		return np.fliplr(jet)
	elif 'l' in pool.lower():
		if l_weight > r_weight:
			return jet
		return np.fliplr(jet)
	else:
		raise ValueError('Jet pooling side must have l -OR- r in the name.')




x = np.array((rotate_jet(j['Intensity'], -j['RotationAngle']), j['LeadingPt']) for j in wprime, dtype=[('image', 'float64', (25, 25)), ('jet_pt',float)])


new = [(rotate_jet(j['Intensity'], -j['RotationAngle']), j['LeadingPt']) for j in wprime]
