"""
Created on Sun Feb 16 19:54:34 2020
@author: Andres Duque
"""

"""
DIFFUSION PSEUDOTIME WITH LANDMARKS     
"""


import warnings 
import numpy as np
import scipy as sci
import graphtools
from scipy.spatial import distance
from sklearn.neighbors import NearestNeighbors

class DPT():
    def __init__(self, k = 10,
                 a = 40, adaptativeK = False, landmarks = True, n_pca = 10,
                 n_landmarks = 2000):
        """
        k : number of neighbors  (default:10)
        a : alpha decay for the kernel (default:40)
        landmarks : Use landmarks (default: True)
        n_landmarks : Number of landmarks (default: 2000)
        n_pca : Number of principal components to use for large data sets (default:10)
        """
        self.k = k
        self.alpha = a
        self.adaptativeK = adaptativeK
        self.landmarks = landmarks
        self.n_pca = n_pca
        self.n_landmark = n_landmarks
        
    def compute_diffusion_opt(self):
        D = distance.squareform(distance.pdist(self.X, 'euclidean'))
        th = 1e-4
        #bth=-np.log(th)**(1/self.alpha)
        self.dimensions = self.X.shape[1]
        self.observations = self.X.shape[0]
        if self.adaptativeK is True:
            kdist= np.sort(D)
            epsilon = kdist[:, self.k+1]
        else:
            epsilon = np.median(D)
            epsilon = np.repeat(epsilon,1)#cambiar
        epsilon[epsilon == 0] = th   
        epsilon = epsilon.reshape((epsilon.shape[0],1))    
        epsilon = np.repeat(epsilon, D.shape[0], axis=1)
        K = np.exp(-(np.divide(D, epsilon))**self.alpha)
        K[K<=th] = 0
        K = K + np.transpose(K)
        self.diff_op = K/K.sum(axis=1)[:,None]
        return self
    
    def compute_diffusion_opt_landmarks(self):
        self.graph = graphtools.Graph(self.X, 
                                      n_landmark = self.n_landmark,
                                      knn=self.k,
                                      n_pca= self.n_pca)
        self.pc_coord = self.graph.transform(self.X)
        self.diff_op = self.graph.landmark_op
        # Compute the diffusion from the landmarks to the points 
        # Using spectral clustering 
        return self
    
    def compute_dpt(self, X):
        self.X = X
        if self.landmarks == True:
            self.compute_diffusion_opt_landmarks()
        else:    
            self.compute_diffusion_opt()
        #vec1, val1, ver1 = np.linalg.svd(self.diff_op)
        val, vec = sci.linalg.eig(self.diff_op)
        weights = np.abs(np.transpose((val/(1-val))[1:,None]))
        weights = (np.repeat(weights, self.diff_op.shape[0], axis=0))**(1/2)
        if self.landmarks == True:
            self.diffusion_time_coord = np.multiply(vec[:,1:], weights)
            self.diffusion_time_coord = self.graph.interpolate(self.diffusion_time_coord)
        else:
            self.diffusion_time_coord = np.multiply(vec[:,1:], weights)
#        nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(self.diffusion_time_coord)
#        distances, indices = nbrs.kneighbors(self.diffusion_time_coord)
        self.dpt = distance.squareform(distance.pdist(self.diffusion_time_coord, 'sqeuclidean'))
        return self        