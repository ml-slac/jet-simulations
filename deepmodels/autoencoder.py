from .activation import *
from .layer import Layer
from .cost import BaseCost, SquaredErrorCost, CrossEntropyCost

import numpy as np
import sys


class Autoencoder(object):
	"""
	Class for Autoencoder
	"""
	def __init__(self, n_visible, n_hidden, activ = (sigmoid, identity), cost=SquaredErrorCost(), **kwargs):
		super(Autoencoder, self).__init__()
		self.n_visible = n_visible
		self.n_hidden = n_hidden
		self.activ = activ

		self.errors = []

		if cost.__class__ is type:
			raise TypeError('Cost must be passed as an\
			 instance of a class deriving from BaseCost')

		self.cost = cost

		# -- Loop through keyword arguments
		encoder_args, decoder_args = {}, {}

		for k in ['learning','momentum','l2_reg']:
			# -- look for particular parameter name in passed args
			try:
				parm = kwargs[k]
				# -- try to convert the arguments to a list...
				try:
					parm = list(parm)
					if len(parm) > 2:
						raise ValueError('Parameter {} can\
						 have a tuple of length at most 2'.format(k))

					encoder_args[k], decoder_args[k] = parm[0], parm[1]
				# -- ...if it's not a tuple, apply the parameter to enc and dec
				except TypeError:
					encoder_args[k], decoder_args[k] = parm, parm
			# -- couldn't find argument
			except KeyError:
				pass

		# -- unpack and initialize
		self.encoder = Layer(n_visible, n_hidden, activ=activ[0], **encoder_args)
		self.decoder = Layer(n_hidden, n_visible, activ=activ[1], **decoder_args)

	def encode(self, X):
		return self.encoder.predict(X)

	def decode(self, R):
		return self.decoder.predict(R)

	def reconstruct(self, X):
		return self.decoder.predict(self.encoder.predict(X))

	def fit(self, X, holdout=True, batch=5, epochs=10, batches_per_epoch=None):
		if holdout is True:
			split = int(0.75 * X.shape[0])
			ix = np.array(range(0, X.shape[0]))
			np.random.shuffle(ix)
			holdout = X[ix[split:]]
			X = X[ix[:split]]

			
		# -- how many batch-sized patches to pull from X?
		n_patches = batches_per_epoch
		if batches_per_epoch is None:
			n_patches = X.shape[0] / batch



		for ep in xrange(0, epochs):

			sys.stdout.write('\rEpoch {} of {}'.format(ep + 1, epochs))
			for patch in xrange(0, n_patches):
				ix = np.random.randint(0, X.shape[0], batch)

				X_ = self.reconstruct(X[ix])

				self.decoder.backpropagate(self.cost.cost_derivative(X_, X[ix]))
				self.encoder.backpropagate(self.decoder.pass_down().T)

			if holdout is not None:
				self.errors.append(self.cost.cost(self.reconstruct(holdout), holdout))

	

class NormalNoise(object):
    """docstring for NormalNoise"""
    def __init__(self, mean=0, sd=0.01):
        super(NormalNoise, self).__init__()
        self.mean = mean
        self.sd = sd
    def corrupt(self, X):
        return X + np.random.normal(self.mean, self.sd, X.shape)

class SaltAndPepper(object):
    """docstring for SaltAndPepper"""
    def __init__(self, p = 0.1):
        super(SaltAndPepper, self).__init__()
        self.p = p
    def corrupt(self, X):
        return np.multiply(X, 1.0 * (np.random.uniform(0, 1, X.shape) > self.p))

class DenoisingAutoencoder(Autoencoder):
	"""docstring for DenoisingAutoencoder"""
	def __init__(self, n_visible, n_hidden, activ = (sigmoid, identity), cost=SquaredErrorCost(), noise=NormalNoise(), **kwargs):
		super(DenoisingAutoencoder, self).__init__(n_visible, n_hidden, activ, cost, **kwargs)
		self.noise = noise

	def noisy_reconstruct(self, X):
		return self.decoder.predict(self.encoder.predict(self.noise.corrupt(X)))

	def fit(self, X, holdout=True, batch=5, epochs=10, batches_per_epoch=None):
		if holdout is True:
			split = int(0.75 * X.shape[0])
			ix = np.array(range(0, X.shape[0]))
			np.random.shuffle(ix)
			holdout = X[ix[split:]]
			X = X[ix[:split]]

			
		# -- how many batch-sized patches to pull from X?
		n_patches = batches_per_epoch
		if batches_per_epoch is None:
			n_patches = X.shape[0] / batch


		n_passes = epochs * n_patches

		prog = 1.0

		for ep in xrange(0, epochs):
			for patch in xrange(0, n_patches):
				sys.stdout.write('\r{}% complete.'.format(np.round(100 * (prog / n_passes), 2)))
				sys.stdout.flush()

				prog += 1.0

				ix = np.random.randint(0, X.shape[0], batch)
				X_ = self.noisy_reconstruct(X[ix])

				self.decoder.backpropagate(self.cost.cost_derivative(X_, X[ix]))
				self.encoder.backpropagate(self.decoder.pass_down().T)

			if holdout is not None:
				self.errors.append(self.cost.cost(self.reconstruct(holdout), holdout))

