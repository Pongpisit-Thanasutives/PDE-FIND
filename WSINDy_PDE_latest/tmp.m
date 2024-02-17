%% Load data

clc;
% close all;
clear all;

pde_num = 3;
pde_names = {'burgers.mat','KdV.mat','KS.mat','NLS.mat','Sine_Gordon.mat','rxn_diff.mat','Nav_Stokes.mat','porous.mat','sod.mat'};
dr = ['datasets/',pde_names{pde_num}];
% dr = ['/home/danielmessenger/Dropbox/Boulder/research/data/WSINDy_PDE/datasets/',pde_names{pde_num}];

try
    %%% don't reload data if already loaded
    U_obs = U_exact;
    xs_obs = xs;
catch
    load(dr);
    U_obs = U_exact;
    xs_obs = xs;
end

%%% select subset of equations
eq = 1:length(U_obs);
lhs = lhs(unique(min(eq,end)),:);
true_nz_weights = true_nz_weights(unique(min(eq,end)));

dims = size(U_obs{1});
dim = length(dims);
n = length(U_obs);

%% Subsample data (if desired)

coarsen_data = repmat([0 1 1],dim,1);

%%% set row d of coarsen_data to [initial_frac inc final_frac] to subsample dth coordinate to 
%%% start at index initial_frac*L, where L is the number of points in dth coordinate
%%% end at index final_frac*L,
%%% skip every inc gridpoint.

coarsen_data(1:dim-1,:) = repmat([0 1 1],dim-1,1); 

[xs_obs,U_obs] = subsamp(xs_obs,U_obs,coarsen_data,dims);
dims = cellfun(@(x) length(x), xs_obs);

%% Add noise

sigma_NR = 0.25;
noise_dist = 0; 
noise_alg = 0;
rng(1);
rng_seed = rng().Seed;
 
rng(rng_seed);
[U_obs,noise,snr,sigma] = gen_noise(U_obs,sigma_NR,noise_dist,noise_alg,rng_seed,0);
