function Theta_cell = get_lib_OL(...
    U,Theta_cell,lib_list,libtree,scales,sub_inds,Cfs_ffts,tt)

    dim = length(sub_inds);
    n = length(U);
    for i=1:length(libtree)
        fcn = ones(size(U{1}));
        tags_root = zeros(1,n);
        indout = libtree{i};
        while ~isempty(indout)
            tags = lib_list(indout(1),1:n);
            sametags = indout(ismember(lib_list(indout,1:n),tags,'rows'));
            for k=1:n
                if isreal(tags(k))
                    fcn = fcn.*(U{k} / scales(k)).^(tags(k)-tags_root(k));
                else
                    if imag(tags(k))<0
                        fcn = sin(abs(imag(tags(k)))*U{k});
                    else
                        fcn = cos(abs(imag(tags(k)))*U{k});
                    end
                end
            end
            for ind = sametags
                test_conv_cell = cell(dim,1);
                for k=1:dim
                    test_conv_cell{k} = Cfs_ffts{k}(lib_list(ind,n+k)+1,:);
                end
                fcn_conv = convNDfft(fcn,test_conv_cell,sub_inds,2);
                foo=repmat({':'},dim,1);
                Theta_cell{ind}(foo{:},tt) = fcn_conv;
            end
            indout = indout(~ismember(indout,sametags));
            tags_root = tags;
        end
    end

end