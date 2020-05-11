function D = InfGeoDis(data, L)


% Compute covariance matrix in each window 
[~, TimeSteps] = size(data);
NumWindows = floor(TimeSteps/L); 
for i = 1:NumWindows
    x = data(:, ((i-1)*L + 1):i*L);
    C{i} = cov(x');
end 

D = zeros(NumWindows, NumWindows);

% Compute Information Geometry Distance (C. Atkinson, A. Mitchell) 
for i = 1:NumWindows
    for j = i+1:NumWindows
    M = C{i}\C{j};
    Eig = eig(M);
    if all(Eig >= 0)
    D(i,j)= 1/2 * sum((log(Eig)).^2);
    else 
    Eig = max(1.0000e-04, Eig);
    D(i,j)= 1/2 * sum((log(Eig)).^2);
    end
    end 
end 


D = D + D';

end 