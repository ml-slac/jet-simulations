# -- global imports
import numpy as np

# -- local imports
from .activation import *
from .cost import ConstantLearningSchedule, LinearLearningSchedule, GeometricLearningSchedule


USE_GPU = False


if USE_GPU:
	from gnumpy import garray as array
	from gnumpy import as_garray as asarray
	from gnumpy import dot as dotproduct
else:
	from numpy import array as array
	from numpy import asarray as asarray
	from numpy import dot as dotproduct


class Layer(object):
	"""
	simple building block for constructing more complicated layer types
	"""
	def __init__(self, inputs=0, outputs=0, activ=sigmoid, learning=0.02, momentum=0.9, l2_reg=0.001):
		super(Layer, self).__init__()

		# -- scalar parameters
		self.inputs = inputs
		self.outputs = outputs
		self.activ = activ
		if (type(learning) == float) or (type(learning) == int):
			self.learning = ConstantLearningSchedule(learning)
		else:
			self.learning = learning
		self.momentum = momentum
		self.l2_reg = l2_reg
 
		# -- weight stuff
		self.W = asarray(0.1 * np.random.normal(0, 1, (self.outputs, self.inputs)))  #/ np.sqrt(outputs + inputs)
		self.b = asarray(0.1 * np.random.normal(0, 1, (self.outputs, 1))) #/ np.sqrt(outputs + inputs)

		# -- store previous weights
		self._W = asarray(np.zeros((self.outputs, self.inputs)))
		self._b = asarray(np.zeros((self.outputs, 1)))

		# -- store weight gradient
		self._grad_W = asarray(np.zeros((self.outputs, self.inputs)))
		self._grad_b = asarray(np.zeros((self.outputs, 1)))


	def predict(self, X):
		if USE_GPU:
			self.X = asarray(X.T)
		else:
			self.X = X.T
		# -- Z = s(W * x + b)
		if USE_GPU:
			if self.activ is sigmoid:
				self.Z = (dotproduct(self.W, self.X) + self.b).logistic()
			elif self.activ is identity:
				self.Z = dotproduct(self.W, self.X) + self.b
		else:
			self.Z = self.activ(dotproduct(self.W, self.X) + self.b)
		return self.Z.T

	def calculate_derivatives(self, err):
		# -- calculate and store the derivatives and the values to pass down
		if USE_GPU:
			if self.activ is sigmoid:
				self.delta = err.T * self.Z * (1 - self.Z)
			elif self.activ is identity:
				self.delta = err.T
		else:
			self.delta = np.multiply(err.T, derivative[self.activ](self.Z))
		self._grad_W = dotproduct(self.delta, self.X.T) / err.shape[0] + self.l2_reg * self.W
		self._grad_b = (self.delta.sum(axis = 1) / err.shape[0]).reshape(self.b.shape)
		self._dump = dotproduct(self.W.T, self.delta)


	def backpropagate(self, err):
		self.calculate_derivatives(err)

		self._W *= self.momentum
		self._b *= self.momentum

		self._W -= self.learning.rate() * self._grad_W
		self._b -= self.learning.rate() * self._grad_b

		self._W -= self.learning.rate() * self.l2_reg * self.W;

		self.W += self._W
		self.b += self._b

	def pass_down(self):
		return self._dump

