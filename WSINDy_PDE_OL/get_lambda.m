function lambda = get_lambda(lambda,lambda_step,lambda_max,obj_fcn_old,obj_fcn_new,sparsity_old,sparsity_new)
    if obj_fcn_new > obj_fcn_old
            if and(all( sparsity_new(:) <= sparsity_old(:) ), any( sparsity_new(:) < sparsity_old(:) )) % valuable term dropped, decrease lambda
                lambda = (1-lambda_step)*lambda;
            elseif and(all( sparsity_new(:) >= sparsity_old(:) ), any( sparsity_new(:) > sparsity_old(:) )) % unnecessary term added, increase lambda
                lambda = (1-lambda_step)*lambda + lambda_step*lambda_max;
            end
    elseif all(sparsity_new(:)==sparsity_old(:))        
            lambda = lambda*(1-lambda_step)+lambda_step*lambda_max;        
    end
end