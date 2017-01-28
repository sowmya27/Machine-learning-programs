__author__ = 'annapurnaannadatha'
from scipy.cluster.vq import kmeans2,kmeans,vq
import numpy as np
from pylab import plot,show


data = np.array([[3.60,79],
              [1.800,54],
              [2.283,62],
              [3.333,74],
              [2.883,55],
              [4.533,85],
              [1.950,51],
              [1.833,54],
              [4.700,88],
              [3.600,85],
              [1.600,52],
              [4.350,85],
              [3.917,84],
              [4.200,78],
              [1.750,62],
              [1.800,51],
              [4.700,83],
              [2.167,52],
              [4.800,84],
              [1.750,47]
             ])

#print(kmeans2(X, 3, iter=10, thresh=1e-05, minit='random', missing='warn', check_finite=True))
'''
# computing K-Means with K = 2 (2 clusters)
centroids,_ = kmeans(data,2)
print(centroids)
# assign each sample to a cluster
idx,_ = vq(data,centroids)

# some plotting using numpy's logical indexing
plot(data[idx==0,0],data[idx==0,1],'ob',
     data[idx==1,0],data[idx==1,1],'or')
plot(centroids[:,0],centroids[:,1],'sg',markersize=8)
show()
'''
# now with K = 3 (3 clusters)
centroids,_ = kmeans(data,3)
print(centroids)
idx,_ = vq(data,centroids)
print(idx)

plot(data[idx==0,0],data[idx==0,1],'ob',
     data[idx==1,0],data[idx==1,1],'or',
     data[idx==2,0],data[idx==2,1],'og') # third cluster points
plot(centroids[:,0],centroids[:,1],'sm',markersize=8)
show()