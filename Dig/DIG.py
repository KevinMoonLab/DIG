
"""
Created on Sat Feb  8 13:36:58 2020
@author: Andres Duque
DIG implementation with Diffusion-PT
"""

"""
Computes the mahalanobis distance in the histograms space, using DPT
"""

import warnings 
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import euclidean_distances

class DIG():
    def __init__(self, n_bins = 20, L1 = 500, L2 = 10, L3=5,
                 histograms_distance = 'mahalanobis', 
                 emb_dimension = 2, dpt = None):
        
        """
		n_bins : Number of bins to create the histograms (default:20)
		L1 : Time window of data to compute the histograms (default:500)
		L2 : Time window of data to compute the covariance matrix for L2 histograms (default:10)
		L3 : Distance between "histogram centers" (default:5)
		histograms_distance : Distance between histograms (euclidean/mahalanobis) (default: mahalanobis)
		emb_dimension : (default:2)
		dpt : dpt matrix for non-ordered data
        """
        
        self.n_bins = n_bins
        self.L1 = L1
        self.L2 = L2
        self.L3 = L3
        self.histograms_distance = histograms_distance
        self.emb_dimension = emb_dimension
        self.dpt = dpt
        
    def compute_histogram(self):

        self.dimensions = self.X.shape[1]
        self.observations = self.X.shape[0]
        if self.dpt is None :
            warnings.warn("You are not using a DPT matrix")
            #cont = input("Is this time series Data? [y/n]:")
            self.ts = 'y'
            self.histogram_ts()
            if self.ts != 'y':
                raise ValueError("DIG can not continue with out time ordering")
        else:
            self.dpt = self.dpt
            self.histogram_dpt()
            self.ts = 'n'
        return self
                
    def histogram_ts(self):
         self.centers_histograms = np.arange(np.ceil(self.L1/2), 
                                             self.observations+np.ceil(self.L1/2), 
                                             self.L3)
         self.histograms = np.array((len(self.centers_histograms), 
                                     self.dimensions*self.n_bins),
                                     dtype = float)
         for i in range(0, self.dimensions):
             hist_dim = np.array([])
             bins_position = list(np.linspace(self.X.min(axis=0)[i], self.X.max(axis=0)[i], num=self.n_bins+1))
             for ji in self.centers_histograms:
                 data_int = self.X[int(ji-np.ceil(self.L1/2)):int(ji+np.ceil(self.L1/2)),i]
                 hist_dim = np.hstack((hist_dim, np.histogram(data_int, bins = bins_position)[0]))
             
             hist_dim =  np.divide(hist_dim, self.L1)  
             hist_dim = hist_dim.reshape(len(self.centers_histograms), 
                                  self.n_bins) 
             if i == 0:
                 self.histograms = hist_dim
             else:
                 self.histograms = np.hstack((self.histograms, hist_dim))                 
         return self
    
    def histogram_dpt(self):
         self.centers_histograms = np.arange(np.ceil(self.L1/2), 
                                             self.observations-np.ceil(self.L1/2), 
                                             self.L3).astype(int)-1
         self.histograms = np.array((len(self.centers_histograms), 
                                     self.dimensions*self.n_bins),
                                     dtype = float)
         for i in range(0, self.dimensions):
             hist_dim = np.array([])
             bins_position = list(np.linspace(self.X.min(axis=0)[i], self.X.max(axis=0)[i], num=self.n_bins+1))
             for ji in self.centers_histograms:
                 data_int = self.X[self.dpt[:self.L1, ji],i]
                 hist_dim = np.hstack((hist_dim, np.histogram(data_int, bins = bins_position)[0]))
             
             hist_dim =  np.divide(hist_dim, self.L1)  
             hist_dim = hist_dim.reshape(len(self.centers_histograms), 
                                  self.n_bins) 
             if i == 0:
                 self.histograms = hist_dim
             else:
                 self.histograms = np.hstack((self.histograms, hist_dim))                 
         return self
    
    def histogramMeans(self):
        histograms_mean = np.zeros((self.histograms.shape))
        if self.ts == 'y':
            for i in range(self.L2, len(self.centers_histograms)-self.L2):
                histograms_mean[i,:] = self.histograms[i-self.L2:i+self.L2,:].mean(axis=0)
        else: 
            Dpt_C = self.dpt[:, self.centers_histograms]
            for i in range(0, len(self.centers_histograms)):
                closests = Dpt_C[np.in1d(Dpt_C[:,i], self.centers_histograms), i]
                histograms_window = self.histograms[np.where(np.in1d(Dpt_C[0,:], closests[:self.L2]))[0],:] 
                histograms_mean[i,:] = histograms_window.mean(axis=0)
        self.histogram_means = histograms_mean
        return self
                 
    def histograms_covariance(self):
        inverse_covariance = np.zeros((self.histograms.shape[1],
                                       self.histograms.shape[1],
                                       self.histograms.shape[0]))
        if self.ts == 'y':
            for i in range(2*self.L2, len(self.centers_histograms)-2*self.L2):
                histograms_window = self.histogram_means[i-self.L2:i+self.L2, :]
                c = np.cov(histograms_window, rowvar = False)
                u, sd, v = np.linalg.svd(c)
                inverse_covariance[:,:,i] = np.matmul(np.matmul(u[:,:self.emb_dimension], 
                                                      np.linalg.inv(np.diag(sd))[:self.emb_dimension,:self.emb_dimension]),
                                                        v[:self.emb_dimension,:])
        else:
            Dpt_C = self.dpt[:, self.centers_histograms]
            for i in range(0, len(self.centers_histograms)):
                #closests = np.intersect1d(Dpt_C[:,i], self.centers_histograms)
                closests = Dpt_C[np.in1d(Dpt_C[:,i], self.centers_histograms), i]
                #histograms_window = self.histogram_means[np.where(np.in1d(Dpt_C[0,:], closests[:self.L2]))[0],:]
                histograms_window = self.histograms[np.where(np.in1d(Dpt_C[0,:], closests[:self.L2]))[0],:] 
                covm = np.cov(histograms_window, rowvar = False)
                u, sd, v = np.linalg.svd(covm)
                inverse_covariance[:,:,i] = np.matmul(np.matmul(u[:,:self.emb_dimension], 
                                                      np.linalg.inv(np.diag(sd))[:self.emb_dimension,:self.emb_dimension]),
                                                        v[:self.emb_dimension,:])
            
        self.inverse_cov = inverse_covariance
    
        return self
    
            
    def compute_hist_distances(self, X):
        self.X = X
        self.compute_histogram()
        self.histogramMeans()
        if self.histograms_distance == 'euclidean':
            Dis = euclidean_distances(self.histogram_means, self.histogram_means)
        else:   
            self.histograms_covariance()
            transpose_hm = np.transpose(self.histogram_means);
            Dis = np.zeros((len(self.centers_histograms), len(self.centers_histograms)))
            for j in range(0,len(self.centers_histograms)):
                tmp1 = np.matmul(self.inverse_cov[:,:,j], transpose_hm[:,j])
                a2 = np.matmul(self.histogram_means[j,:], tmp1)
                b2 = np.sum(np.multiply(self.histogram_means,
                            np.transpose(np.matmul(self.inverse_cov[:,:,j], transpose_hm))),
                            axis = 1)
                ab = np.matmul(self.histogram_means, tmp1)
                Dis[:,j] = np.repeat(a2, len(self.centers_histograms)) + b2 - 2*ab
        D = Dis + np.transpose(Dis)
        np.fill_diagonal(D, 0)
        self.histogram_distances = np.abs(D)
        return self                                           
            
                                                  
         
                 
    
                
            
    
                
                 
             
         
       
            
        












