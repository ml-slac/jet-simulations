# -- Global import
import numpy as np

# -- Local import 
from autoencoder import NormalNoise, SaltAndPepper, DenoisingAutoencoder
from .activation import *
from collections import OrderedDict


model = [{'type' : 'Denoising autoencoder',
		  'inputs' : 25 ** 2, 
		  'outputs' : 26 ** 2, 
		  'activation' : (sigmoid, sigmoid),
		  'learning' : (GeometricLearningSchedule(0.1, 0.001, 0.95), 
                        GeometricLearningSchedule(0.01, 0.001, 0.95)),
		  'momentum' : (0.96, 0.97)}, 

		  
					 'layer 2' : {}, 
					 'layer 3' : {}
					 })



class DeepAE(object):
	"""docstring for DeepAE"""
	def __init__(self, structure):
		super(DeepAE, self).__init__()
		self.structure = structure
		
