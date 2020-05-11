function Pmm = affinity_matrix(D,k,a, adaptativeK)

[kdist, idx] = sort(D,2);
th = 1e-4;
k_knn = k * 20;
bth=(-log(th))^(1/a);
if adaptativeK == true
    epsilon = kdist(:,k+1);
else 
    epsilon = median(D, 'all');
    epsilon = repmat(epsilon, 1, size(D, 1));
end 

below_thresh=kdist(end,:)>=bth*epsilon;
epsilon(epsilon == 0) = mean(epsilon);
idx_thresh=find(below_thresh);
idx_thresh=1:size(D,2);
if ~isempty(idx_thresh)
    K=exp(-(D./repmat(epsilon(idx_thresh)',size(D,1),1)).^a);
    K(K<=th)=0;
end

K = K + K';
Pmm = bsxfun(@rdivide, K, sum(K,2));

end 