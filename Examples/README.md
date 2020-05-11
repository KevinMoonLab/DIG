# DIG demos 


### How to use it 

#### Data with time ordering

In the case the data is naturally ordered by time, as in `EEG_demo.ipynb`. Where the rows are the measurements and the columns are channels/features/sensors.

```
dig = DIG(L1=3840, L3=3840)
hist = dig.compute_hist_distances(Data_eeg)
```

#### Data with no time ordering

When the data has no order, the neighborhood of each observation is compute by the [Diffusion Pseudo Time](https://www.nature.com/articles/nmeth.3971) distance. 
We implemented a landmark version in `DPT.py` in order to deal with big data sets, and it is recommended to preprocessing the data by applying  PCA to reduce the dimensionality first.

```
dpt_c = DPT()
dpt = dpt_c.compute_dpt(data_pca)

# Order the dpt for each observation 
dpt_or = np.argsort(dpt.dpt, axis=0)
```

Now you can use DIG using diffusion pseudo time

```
dig = DIG(dpt = dpt_or, n_bins = 10, L1 = 1000, L2 = 10, L3=10)
hist = dig.compute_hist_distances(data_pca)
```

Once you have the distances in the histograms space, you can apply PHATE with precomputed distance. 

```
phate_emb = phate.PHATE(knn_dist='precomputed_distance', knn = 20)
phate_fit = phate_emb.fit_transform(hist.histogram_distances)
```

### DIG parameters

```
n_bins : Number of bins to create the histograms (default:20)
L1 : Time window of data to compute the histograms (default:500)
L2 : Time window of data to compute the covariance matrix for L2 histograms (default:10)
L3 : Distance between "histogram centers" (default:5)
histograms_distance : Distance between histograms (euclidean/mahalanobis) (default: mahalanobis)
emb_dimension : (default:2)
dpt : dpt matrix for non-ordered data
```

### DATA sets

* EB data from [here](https://github.com/KrishnaswamyLab/PHATE/blob/master/data/EBdata.mat), preprocessed using PCA and retaining 10 dimensions. 
* EEG data 
