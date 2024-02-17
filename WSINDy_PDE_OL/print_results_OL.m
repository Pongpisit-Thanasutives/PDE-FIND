warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle')
[m,num_eq] = size(W);

str_wsindy = cell(num_eq,1);
for k=1:num_eq
    tags_pde_rdx = tags_pde(~ismember(1:m+num_eq,lhs_ind));
    str_wsindy{k} = print_pde(W(:,k),tags_pde_rdx,tags_pde{lhs_ind(k)});
end

if ~isempty(axi)
    if isequal(class(axi),'double')
        Tps = tpscore({W},axi);
        dW = wnorm({W},axi,Inf);
    elseif isequal(class(axi),'cell')
        Tps = tpscore({W},axi{(Kmem+1)/2+tOL-1});
        dW = wnorm({W},sum(Cfs_t(1,:)/norm(Cfs_t(1,:),1).*cell2mat(axi(tOL:Kmem+tOL-1)),2),Inf);
    end
else
    Tps=NaN;
    dW=[];
end

if ~isequal(print_loc,0)
    if ~isequal(print_loc,1)
        print_loc = fopen(print_loc,'a');
    end

    fprintf(print_loc,'\n-------------------online iteration number: %u \n',tOL);

    for k=1:num_eq
        fprintf(print_loc,['\nRecovered PDE: ',str_wsindy{k}]);
        fprintf(print_loc,'\nRelative Res: ||b-G*W||_2/||b||_2 = %.2e\n',norm(resid(:,k)/norm(b)));
        if ~isempty(dW)
            fprintf(print_loc,'\nMax Weight Error: max|W-W_{true}| = %.2e\n', dW(k));
        end
    end
    for nn=1:n
        cellfun(@(x,y)fprintf(print_loc,[x,' ||G*W||/||b|| = %3.8f \n'],y),tags_pde_G(W(:,nn)~=0),...
            num2cell(vecnorm(G(:,W(:,nn)~=0).*W(W(:,nn)~=0,nn)')/norm(b(:,nn))),'uni',0);
    end
    fprintf(print_loc,'TP Score = %1.2f\n', Tps);
    fprintf(print_loc,'\npolys = ');
    fprintf(print_loc,'%u ',polys);
    fprintf(print_loc,'\ntrigs = ');
    fprintf(print_loc,'%u ',trigs);
    fprintf(print_loc,'\nMax derivs [t x] = ');
    fprintf(print_loc,'%u ',[max_dt max_dx]);
    fprintf(print_loc,'\n[m_x m_t] = ');
    fprintf(print_loc,'%u ',[m_x m_t]);
    fprintf(print_loc,'\n[s_x s_t] = ');
    fprintf(print_loc,'%u ',[max(sub_inds{1}(end)-sub_inds{1}(max(end-1,1)),1) max(sub_inds{end}(end)-sub_inds{end}(max(end-1,1)),1)]);
    fprintf(print_loc,'\n[p_x p_t] = ');
    pps=zeros(1,2);
    ps=[p_x p_t];
    for i=1:2
        if isequal(phi_class{i},1)
            pps(i)=ps(i);
        elseif isequal(phi_class{i},2)
            pps(i)=1/ps(i);
        else
            pps(i)=NaN;
        end
    end
    fprintf(print_loc,'%u ',pps);
    fprintf(print_loc,'\nSize of 1 snapshot (Mb) = ');
    fprintf(print_loc,'%f',whos('U').bytes/10^6);
    fprintf(print_loc,'\nDimensions of Kmem snapshots = ');
    fprintf(print_loc,'%u ',[size(U_obs_gap{1}) n]);
    fprintf(print_loc,'\nDimensions of G = ');
    fprintf(print_loc,'%u ',size(G));
    fprintf(print_loc,'\nCondition number of G = ');
    fprintf(print_loc,'%u ',cond(G));
    fprintf(print_loc,'\nlambda = ');
    fprintf(print_loc,'%.3e ',lambda);
    fprintf(print_loc,'\nsigma_NR = ');
    fprintf(print_loc,'%.3e \n',sigma_NR);
    fprintf(print_loc,'\nET_load_data = %4.4f \n',ET_load_data);
    fprintf(print_loc,'ET_offline = %4.4f \n',ET_offline);
    fprintf(print_loc,'ET_online = %4.4f \n',ET_build_Gb+ET_online_iteration);
    fprintf(print_loc,'\nnumber of snapshots (Kmem) = %u \n',Kmem);
    fprintf(print_loc,'astep = %2.16f \n',astep);
    fprintf(print_loc,'||w||_0 = %u \n',length(find(W)));
    fprintf(print_loc,'lower. gap: %3.8f \n',lower_gap);
    fprintf(print_loc,'upper. gap: %3.8f \n',upper_gap);
    if ~all(print_loc==1)
        fclose(print_loc);
    end
end