# -- global imports
import numpy as np

# -- local imports
from .activation import *
from .cost import ConstantLearningSchedule, LinearLearningSchedule, GeometricLearningSchedule

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
		self.W = 0.1 * np.random.normal(0, 1, (self.outputs, self.inputs))  #/ np.sqrt(outputs + inputs)
		self.b = 0.1 * np.random.normal(0, 1, (self.outputs, 1))  #/ np.sqrt(outputs + inputs)

		# -- store previous weights
		self._W = np.zeros((self.outputs, self.inputs))
		self._b = np.zeros((self.outputs, 1))

		# -- store weight gradient
		self._grad_W = np.zeros((self.outputs, self.inputs))
		self._grad_b = np.zeros((self.outputs, 1))


	def predict(self, X):
		self.X = X.T

		# -- Z = s(W * x + b)
		self.Z = self.activ(np.dot(self.W, self.X) + self.b)
		return self.Z.T

	def calculate_derivatives(self, err):
		# -- calculate and store the derivatives and the values to pass down
		self.delta = np.multiply(err.T, derivative[self.activ](self.Z))
		self._grad_W = np.dot(self.delta, self.X.T) / err.shape[0] + self.l2_reg * self.W
		self._grad_b = (self.delta.sum(axis = 1) / err.shape[0]).reshape(self.b.shape)
		self._dump = np.dot(self.W.T, self.delta)


	def backpropagate(self, err):
		self.calculate_derivatives(err)

		self._W *= self.momentum
		self._b *= self.momentum

		self._W -= self.learning * self._grad_W
		self._b -= self.learning * self._grad_b

		self._W -= self.learning * self.l2_reg * self.W;

		self.W += self._W
		self.b += self._b

	def pass_down(self):
		return self._dump

