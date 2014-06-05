import plot as plt
from numpy import *

k = 13

path = "h1000.dat"

data = array(plt.parse(path))

blur = 15
smoothed = array(data[:,1])
for i in range(blur, len(smoothed)-blur):
    for j in range(-blur, blur):
        smoothed[i] += data[i+j,1]
    smoothed[i] /= (2*blur)

fig, ax = plt.new_plot(title="K-means, %d Clusters" % k,
                       xaxis="Time [ms]",
                       yaxis="K-means distortion")
plt.add(arange(len(smoothed))*.001, smoothed, 
        "K-Means Cluster Distortion Smoothed")
plt.dot(arange(len(data))*.001, data[:,1], 
        "K-Means Cluster Distortion")
plt.save_plot(fig, ax, "plt1.png")
