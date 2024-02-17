function [U_obs,xs_obs] = coarsen_data(U_obs, coarse_data_pattern, xs_obs)
    if ~isempty(coarse_data_pattern)
        n = length(U_obs);
        dims = size(U_obs{1});
        dim = length(size(U_obs{1}));
        inds = cell(1,dim);
        if size(coarse_data_pattern,1)~=dim
            coarse_data_pattern=[coarse_data_pattern;[0 1 1]];
        end
        for j=1:dim
            inds{j} = 1+floor(coarse_data_pattern(j,1)*dims(j)):ceil(coarse_data_pattern(j,2)):ceil(coarse_data_pattern(j,3)*dims(j));
            xs_obs{j} = xs_obs{j}(inds{j});
        end
        for j=1:n
            U_obs{j} = U_obs{j}(inds{:});
        end
    end
end