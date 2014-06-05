import plot as plt
import scipy.optimize as optimization
from numpy import *

k = 13

path = "h_125.dat"

data = array(plt.parse(path))

f = lambda x, a, b, c, d: a*exp(-b*x**c) + d

fig, ax = plt.new_plot(title="K-means, %d Clusters" % k,
                       xaxis="Time [ms]",
                       yaxis="K-means distortion")

p0 = [2.7, 1., 1.1, .25]
p, e = optimization.curve_fit(f, data[:,0], data[:,1], p0)

print p

X = linspace(0, max(data[:,0]), 100)
plt.dot(data[:,0]*.005, data[:,1], "H = 50 mT")
plt.add(X*.005, f(X, p[0], p[1], p[2], p[3]), "Exponential fit")
plt.save_plot(fig, ax, path[:-4]+".png")
