import plot as plt
import scipy.optimize as optimization
from numpy import *

k = 13

path = "h1.dat"

data = array(plt.parse(path))

f = lambda x, a, b, c:  a*exp(-x*b) + c

p0 = [.35, 1, .05]
p, e = optimization.curve_fit(f, data[:,0], data[:,1], p0)
print p

fig, ax = plt.new_plot(title="K-means, %d Clusters" % k,
                       xaxis="Time [ms]",
                       yaxis="K-means distortion")

X = linspace(0, max(data[:,0]), 100)
plt.dot(data[:,0], data[:,1], "H = 100 mT")
plt.add(X, f(X, p[0], p[1], p[2]), "Exponential fit")
plt.save_plot(fig, ax, "plt3.png")
