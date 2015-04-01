# -- Global import
import numpy as np

# -- Local import 
from autoencoder import NormalNoise, SaltAndPepper, DenoisingAutoencoder




class DeepAE(object):
	"""docstring for DeepAE"""
	def __init__(self, structure = None, learning = None, momentum = None, ):
		super(DeepAE, self).__init__()
		self.structure = structure
