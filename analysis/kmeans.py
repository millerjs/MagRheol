import plot as plt
from scipy.cluster.vq import vq, kmeans, whiten
from numpy import *

n = 800
k = 13

def cluster(i):
    path = "../checkpoints/mag_%05d.dat" % i
    data = array(plt.parse(path))[:, 0:2]
    
    whitened = whiten(data)
    book = array((whitened[0],whitened[2]))
    
    centroids, distort = kmeans(whitened, k)
    print "%d\t%f" % (i, distort)
    return distort
    
y = []
for i in range(1,n):
    cluster(i)

