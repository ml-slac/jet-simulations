import numpy as np


# -- Cost functions!
class BaseCost(object):
	"""docstring for BaseCost"""
	def __init__(self):
		super(BaseCost, self).__init__()

	def cost(self, predicted, target):
		raise NotImplementedError('Not implemented in base class.')
		
	def cost_derivative(self, predicted, target):
		raise NotImplementedError('Not implemented in base class.')


class SquaredErrorCost(BaseCost):
	"""docstring for SquaredErrorCost"""
	def __init__(self, fraction=0.5):
		super(SquaredErrorCost, self).__init__()
		self.fraction = fraction
	
	def cost(self, predicted, target):
		return self.fraction * np.sum(predicted - target)

	def cost_derivative(self, predicted, target):
		return 2 * self.fraction * (predicted - target)

class CrossEntropyCost(BaseCost):
	"""docstring for CrossEntropyCost"""
	def __init__(self):
		super(CrossEntropyCost, self).__init__()
	
	def cost(self, predicted, target):
		return -np.sum(np.multiply(target, np.log(predicted)) + np.multiply((1 - target), np.log(1 - predicted)))

	def cost_derivative(self, predicted, target):
		E = predicted - target
		return np.divide(E, (predicted - predicted ** 2))


# -- Learning rate schedules
class BaseLearningSchedule(object):
	"""docstring for BaseLearningSchedule"""
	def __init__(self, initial = None, final = None):
		super(BaseLearningSchedule, self).__init__()
		self.initial = initial
		self.final = final
		self._rate = initial

		self.at_final = False

	def rate(self):
		raise NotImplementedError('Not implemented in base class.')

	def next(self):
		raise NotImplementedError('Not implemented in base class.')

	def set(self, rate):
		self.rate = rate

class ConstantLearningSchedule(BaseLearningSchedule):
	"""docstring for ConstantLearningSchedule"""
	def __init__(self, rate = 0.1):
		super(ConstantLearningSchedule, self).__init__()
		self._rate = rate
	def rate(self):
		return self._rate	
	def next(self):
		return


class LinearLearningSchedule(BaseLearningSchedule):
	"""docstring for LinearLearningSchedule"""
	def __init__(self, initial = 0.1, final = 0.0001, steps = 100):
		super(LinearLearningSchedule, self).__init__()
		self.initial = initial
		self.final = final
		self._steps = steps
		self._delta = (initial - final) / steps

		self._rate = initial

	def rate(self, decrement = False):
		if decrement:
			self.next()
		return self._rate

	def next(self):
		if self._rate > self.final:
			self._rate -= self._delta


class GeometricLearningSchedule(BaseLearningSchedule):
	"""docstring for GeometricLearningSchedule"""
	def __init__(self, initial = 0.1, final = 0.0001, decay = 0.9):
		super(GeometricLearningSchedule, self).__init__()
		self.initial = initial
		self.final = final
		self._decay = decay

		self._rate = initial


	def rate(self, decrement = False):
		if decrement:
			self.next()
		return self._rate

	def next(self):
		if self._rate > self.final:
			self._rate *= self._decay





		

		

		







