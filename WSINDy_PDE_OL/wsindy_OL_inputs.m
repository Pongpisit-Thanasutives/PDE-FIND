%%%%%%%%%%%%%%% data file specifications
%%% ------------------------------------
%%% data U_exact, xs, lhs, true_nz_weights will be added from 
%%% [data_dr,pde_names{pde_num}] 

data_dr = '/Users/pongpisit/Desktop/research/PDE-FIND/WSINDy_PDE_latest/inputs/KS/';
pde_names = {'burgers.mat',...          %1
    'KdV.mat',...                       %2
    'KS_largetime.mat',...              %3
    'NLS.mat',...                       %4
    'SG_largetime_uu3.mat',...          %5
    'rxn_diff.mat',...                  %6
    'Nav_Stokes.mat',...                %7
    'porous2.mat',...                   %8
    'KS.mat',...                        %9
    'Sine_Gordon.mat',...               %10
    'AC_temp.mat',...                   %11
    'qanr1Dcont.mat',...                %12
    'KS_1024.mat',...                   %13
    'wave3D.mat',...                    %14
    'SG_largetime',...                  %15
    'porous.mat',...                    %16
    'wave3D_N128.mat',...               %17
    'waveVCspeed.mat',....              %18
    'waveVCspeedNL.mat',...             %19
    'wave2Du3.mat',...                  %20
    'wave2DVCspeedu3.mat',...           %21    
    'wave2DVCspeedu3_c0.5.mat',...      %22
    'wave2DVCspeeduu3_c0.2.mat',...     %23
    'wave3D_13-Feb-2022_1.mat'...,      %24
    'wave2DVCspeedu3_c1.mat',...        %25
    'wave2DVCspeedu3_c0.9.mat',...      %26
    'wave2DVCspeedu3_c0.95.mat',...     %27
    'wave2DVCspeedu3_c0.95-19-May-2022.mat', ...    %28
    'wave2DVCspeedu3_c0.5-14-Aug-2022.mat', ... %29
    };

pde_num = 3;
ode_num = 0;

%%%%%%%%%%%%%%% coarsen data
%%% ------------------------
%%% each row k is a triple [L s R] corresponding to subsampling in
%%% coordinate k: xs{k} <- xs{k}(L:s:R)

coarse_data_pattern = [[0 1 1];[0 1 1];[0 1 1]]; 

%%%%%%%%%%%%%%% add noise to data
%%% -----------------------------
%%% sigma_NR=||ep||/||U_exact||, ep is the noise, U_exact is noise-free data
%%% comment out rng('shuffle') to reproduce results
%%% noise_dist chooses Gaussian (0) or uniform (1) noise
%%% noise_alg chooses additive (0) or multiplicative (1) noise

sigma_NR = 0.0;
% rng('shuffle');
noise_dist = 0;
noise_alg = 0;

%%%%%%%%%%%%%%% number of time snapshots allowed in memory
%%% ------------------------------------------------------

Kmem = floor(251/2);

%%%%%%%%%%%%%%% number of online iterations 
%%% ---------------------------------------
%%% (use inf to run on entire dataset)

maxOLits = inf;

%%%%%%%%%%%%%%% method of computing initial guess
%%% ---------------------------------------------
%%% options: 'LS' for least squares, 'zeros' for all zeros, 'randn' for standard Gaussian entries

ICmeth = 'LS';

%%%%%%%%%%%%%%% Tikhonov regularization 
%%% -----------------------------------
%%% adds 0.5*gamma*||w||^2 to cost functions
gamma = 0;

%%%%%%%%%%%%%%% gradient descent stepsize
%%% -------------------------------------
%%% either aa_fac/||G*G_S|| or astep0 (see get_astep.m for details) 

aa_fac = 1;
astep0 = 0;

%%%%%%%%%%%%%%% iterative sparse regression method
%%% ----------------------------------------------
%%% 'IVHT', 'IHTs', 'LASSO' (see get_sparse_update.m for details)

SRmeth = 'IVHT';
lambda_init = 0.0001;
lambda_step = 0.1;
lambda_max = 0.1;
upper_bound = 10^2;

%%%%%%%%%%%%%%% PDE library
%%% -----------------------

max_dx = 4; % use spatial derivatives up to this order
max_dt = 1; % use temporal derivatives up to this order
polys = [0:2]; % use polynomials of state variables 
trigs = [];  % use trig functions of state variables
use_all_dt =  0; % use all temporal derivatives 
use_cross_dx = 0; % use all cross derivatives
custom_add = []; % add custom terms to the library [poly_tag derivative_tag]
custom_remove = {}; % remove terms using conditions {@(mat) mat(1,:)==0} will remove terms that depend on first state variable  

%%%%%%%%%%%%%%% weak discretization
%%% -------------------------------

max_rows = 10000; % maximum allowed rows in matrix, sets Query points, s_x
phi_class = {1,1}; % 1=piecewise polynomial, 2=Gaussian, {space,time}

use_cornerpt = 1; % use cornerpoint algorithm applied to first Kmem snapshots
% tauhat_x = 2; tau_x = 10^-10;
tauhat_x = 5; tau_x = 10^-6; % Only for my kuramoto_sivishinky datasets
tauhat_t = 1; tau_t = 10^-10;

m_x = 10; % values chosen if use_cornerpt = 0
p_x = 9;
m_t = floor((Kmem-1)/2);
p_t = 9;

%%%%%%%%%%%%%%% Display results during online iterations
%%% ----------------------------------------------------

print_loc = 1;
toggle_OL_print=10;

%%%%%%%%%%%%%%% Plot results during online iterations
%%% -------------------------------------------------

toggle_plot_basis_fcn = 0;
toggle_plot_sol = 0;
toggle_plot_weights = 200;
toggle_plot_fft = 0;

%%%%%%%%%%%%%%% Save results during online iterations
%%% -------------------------------------------------
toggle_save = 1;
