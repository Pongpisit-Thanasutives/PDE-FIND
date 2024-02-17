%%%%%%%%%%%% load observation data from file

if pde_num>0
    pde_name = pde_names{pde_num};
    load([data_dr, 'kuramoto_sivishinky30_ksvdreg.mat'],'U_exact','xs','lhs');
    U_obs = cellfun(@real, U_exact, 'UniformOutput', false);
    U_obs = cellfun(@double, U_obs, 'UniformOutput', false);
    xs_obs = cellfun(@double, xs, 'UniformOutput', false);
    lhs = double(lhs);
elseif ode_num>0
    ode_data_gen;
else
    try 
        U_obs = U_exact;
        xs_obs = xs;
        pde_name='';
    catch
        disp(['no data found'])
    end
end

if ~exist('true_nz_weights','var')
    true_nz_weights=[];
end

%%%%%%%%%%%% coarsen data: rewrites coarsened versions to U_obs, xs_obs

[U_obs,xs_obs] = coarsen_data(U_obs, coarse_data_pattern, xs_obs);

%%%%%%%%%%%% add noise
% rng_seed = rng().Seed; rng(rng_seed);
% [U_obs,noise,snr,sigma] = gen_noise(U_obs,sigma_NR,noise_dist,noise_alg,...
%     rng_seed,0);