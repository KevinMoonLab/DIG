function D = KLDis(data, L)


% Compute covariance matrix in each window 
[p, TimeSteps] = size(data);
NumWindows = floor(TimeSteps/L); 
for i = 1:NumWindows
    x = data(:, ((i-1)*L + 1):i*L);
    C{i} = cov(x');
end 

D = zeros(NumWindows, NumWindows);

% Compute KL divergence zero mean Multivariate Gaussians 
for i = 1:NumWindows
    for j = 1:NumWindows
    D(i,j)= 1/2 * (trace(C{i}\C{j}) - p + log(det(C{i})/det(C{j})));
    end 
end 

D = D + D';

end 