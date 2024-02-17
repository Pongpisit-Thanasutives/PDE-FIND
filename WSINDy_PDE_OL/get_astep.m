function astep = get_astep(G,M,sparsity,aa_fac,astep0)
    if ~any(sparsity(:))
        astep=1; % if the zero vector, use 1.
    else
        if astep0 == 0
            if aa_fac>0
                astep = aa_fac/norm((M.*G')*(G(:,sum(sparsity,2)~=0).*M(sum(sparsity,2)~=0)'));
            else
                s=max(sum(sparsity));
                astep=2/sqrt(s)/min(sqrt(size(sparsity,1))-1,sqrt(size(sparsity,1)-s)+1);
            end
        else
            astep = astep0;
        end
    end
