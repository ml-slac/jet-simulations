import numpy as np
from scipy.special import expit

def sigmoid(x):
	return expit(x)

def sigmoid_derivative(x):
	return np.multiply(x, 1 - x)

def identity(x):
	return x

def identity_derivative(x):
	return 1.0

derivative = {sigmoid : sigmoid_derivative, identity : identity_derivative}