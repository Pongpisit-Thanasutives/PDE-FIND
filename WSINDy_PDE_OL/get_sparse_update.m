function H = get_sparse_update(SRmeth)
    
    if isequal(SRmeth,'IVHT') % iterative (variable) hard thresholding
        H = @(x,lambda,LB,astep) x.*(abs(x)>=LB*astep*lambda);
    
    elseif isequal(SRmeth,'IHTs') % iterative hard thresholding on s largest entries
        H = @(x,s,LB,astep) get_largest(x,s);
        lambda_max = 0;
    
    elseif isequal(SRmeth,'LASSO') % prox. GD on L1-regularized least squares
        H = @(x,lambda,LB,astep) sign(x).*max(abs(x)-astep*lambda,0);

end