function  [D,  info] = EIGDis(data, L, ndim, nbins, L2)
[N, len] = size(data);
%% Concatenate 1D histograms (marginals) of each time series in short windows 
z_hist = [];
for dim=1:N
    z_hist_dim = [];
    hist_bins = linspace(min(data(dim,:)), max(data(dim,:)), nbins);
    for i=1:floor(len/L)
        % Take the data in this interval of time 
        interval = data(dim, 1+(i-1)*L:i*L);
        % create the histogram with the measurements taken above
        z_hist_dim(:, end+1) = (hist(interval, hist_bins))'./L;
    end
    z_hist = [z_hist z_hist_dim'];
end
%% Covariance estimation
% Store the mean histogram in each local neighborhood
z_mean = zeros(size(z_hist)); 
% Store the inverse covariance matrix of histograms in each local neighborhood
inv_c = zeros(N*length(hist_bins), N*length(hist_bins), length(z_hist)); 
% Store the covariance matrix of histograms in each local neighborhood
c = zeros(N*length(hist_bins), N*length(hist_bins), length(z_hist)); 
% HAY PROBLEMAS CUANDO LAS OBSERVACIONES SON MENOS Q EL TAMANO DEL
% HISTOGRAMA
for i=1+L2:length(z_hist)-L2
    % Estimate covariance in short time windows
    win = z_hist(i-L2:i+L2-1,:);
    c(:,:,i) = cov(win);
    % Denoise via projection on "known" # of dimensions
    [U, S, V] = svd(c(:,:,i));
    inv_c(:,:,i) = U(:,1:ndim) * inv(S(1:ndim,1:ndim)) * V(:,1:ndim)';
    %inv_c(:,:,i) = eye(size(win,2));
    % mean of each feature( histogram bin) in the interva;l
    z_mean(i, :) = mean(win,1);
end
info.z_mean = z_mean;
info.cov = c;
info.invcov = inv_c;
M = size(z_mean, 1);

%% Mahalanobis Distances 
Dis = zeros(M, M);
h = waitbar(0, 'Please wait');
for j = 1:M
    waitbar(j/M, h);
    tmp1 = inv_c(:,:,j) * z_mean(j,:)';
    a2 = z_mean(j,:) * tmp1;
    b2 = sum(z_mean .* (inv_c(:,:,j) * z_mean')',2);
    ab = z_mean * tmp1;
    Dis(:,j) = repmat(a2, M, 1) + b2 - 2*ab;
end
close(h)
D = abs(Dis+Dis');
D = D - diag(diag(D));
info.Hist_distances = D;
end 