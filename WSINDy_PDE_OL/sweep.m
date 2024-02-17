%% set parameters, load data 
clear all;
runs = 100;

tic;
wsindy_OL_inputs;
if pde_num>0
    pde_name = pde_names{pde_num};
    load([data_dr,pde_name],'U_exact','xs','lhs','true_nz_weights');
elseif ode_num>0
    ode_data_gen;
else
    try 
        pde_name='';
    catch
        disp(['no data found'])
    end
end
if ~exist('true_nz_weights','var')
    true_nz_weights=[];
end
[U_exact,xs_obs] = coarsen_data(U_exact, coarse_data_pattern, xs);
ET_load_data = toc;

Wscell = {};
for j=1:runs

wsindy_OL_inputs;
rng_seed = rng().Seed; rng(rng_seed);
[U_obs,noise,snr,sigma] = gen_noise(U_exact,sigma_NR,noise_dist,noise_alg,...
    rng_seed,0);

offline_computations;
for tOL=2:T+1

    %%%%%%%%%%%%%% build new linear system using new snapshot U
    tic;
    
    U = cellfun(@(x)x(space_inds{:},tOL+Kmem-1),U_obs,'uni',0); % get new snapshots
    U_obs_gap=cellfun(@(x)cat(dim,x(space_inds{:},2:Kmem),U{1}),U_obs_gap,'uni',0); % relevant snapshots, only for plotting
    xs_OL{end}={xs_obs{end}(tOL:tOL+Kmem-1)}; % update time grid
    
    Theta_cell=cellfun(@(x)x(space_inds{:},2:Kmem),Theta_cell,'uni',0); % remove oldest snapshot
    Theta_cell = get_lib_OL(...
        U,Theta_cell,lib_list,libtree,scales,sub_inds(1:end-1),Cfs_ffts,Kmem); % add weak-form terms from new snapshot
    [G,b,~]=build_linear_system(Theta_cell,Cfs_t_scaled,lib_list,lhs_ind,[]); % build new linear system
    
    ET_build_Gb = toc;

    %%%%%%%%%%%%%% step forward 
    tic;

    M=1./max(vecnorm(G)',eps); % 
    LB = max(norm(b).*M,1);

    astep = get_astep(G,M,sparsity,aa_fac,astep0);

    Z = W - astep*(M.^2.*(G'*(G*W-b)) -gamma^2*W);
    W = H(Z,lambda,LB,astep);
    W = get_proj_op(W,upper_bound);

    %%%%%%%%%%%%%% update lambda
    
    sparsity_old = sparsity;
    resid_old = resid;
    resid = b-G*W;
    sparsity = W~=0;
    obj_fcns(tOL) = norm(resid(:))^2+sum(sum((W~=0).*(LB./M*lambda).^2));
    lambda = get_lambda(lambda,lambda_step,lambda_max,obj_fcns(tOL-1),obj_fcns(tOL),sparsity_old,sparsity);
    
    ET_online_iteration = toc;

    %%%%%%%%%%%%% Record optimality gaps
    
    tic; 
    lower_gap = zeros(num_eq,1);
    upper_gap = zeros(num_eq,1);

    upper_lim = abs(W./LB/astep);
    lower_lim = abs((G'*resid).*(M.^2./LB));

    for ne=1:num_eq
        if all(~sparsity(:,ne))
            lower_gap(ne)=lambda-max(lower_lim(:,ne));
            upper_gap(ne)=-lambda;
        elseif all(sparsity(:,ne))
            lower_gap(ne)=lambda;
            upper_gap(ne)= min(upper_lim(:,ne)) - lambda;
        else 
            lower_gap(ne)= lambda-max(lower_lim(~sparsity(:,ne),ne));
            upper_gap(ne)= min(upper_lim(sparsity(:,ne),ne)) - lambda;
        end        
    end

    ET_save = toc;

    %%%%%%%%%%%%% save results

    if toggle_save
        Ws = [Ws, {W}];
        ETs = [ETs, ET_online_iteration+ET_build_Gb+ET_save];
    end

    %%%%%%%%%%%%% Display results

    if and(toggle_OL_print,mod(tOL,toggle_OL_print)==0)
        clc;
        print_results_OL;
        display_results;
    end


end

Wscell = [Wscell,{Ws}];
disp([j sum(ETs)])
end