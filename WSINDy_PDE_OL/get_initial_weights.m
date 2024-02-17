function W_0=get_initial_weights(G,b,gamma,ICmeth)
    if isequal(ICmeth,'LS') % least squares on initial snapshots
        W_0 = [G;gamma*eye(size(G,2))] \ [b;zeros(size(G,2),size(b,2))];

    elseif isequal(ICmeth,'zeros') % all zeros
        W_0 = zeros(size(G,2),size(b,2));
    
    elseif isequal(ICmeth,'randn') % random Gaussian entries
        W_0 = randn(size(G,2),size(b,2));
    
    end
end