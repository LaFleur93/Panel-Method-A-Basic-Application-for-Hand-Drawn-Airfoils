import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def parametric_airfoil(nodes):
	data = np.array(nodes)

	tck,u = interpolate.splprep(data.transpose(), s=0)
	unew = np.arange(0, 1.01, 0.01)#0.01
	out = interpolate.splev(unew, tck)

	return out, data