%% get library terms

tic,
n = length(U_obs);
dims = cellfun(@(x) length(x), xs_obs);
dim = length(dims);
num_eq = size(lhs,1);

[tags_pde,lib_list,pdx_list,lhs_ind,max_dx,max_dt,polys] = get_lib_tags(...
    n,dim,lhs,max_dx,max_dt,polys,trigs,use_cross_dx,use_all_dt,...
    custom_remove,custom_add,true_nz_weights);
tags_pde_G = tags_pde(~ismember(1:length(tags_pde),lhs_ind));
lib_list_G = lib_list(~ismember(1:size(lib_list,1),lhs_ind),:);

%% set weak discretization using first Kmem snapshots

dx = mean(diff(xs_obs{1}));
dt = mean(diff(xs_obs{end}));

tic,
dimsgap = [cellfun(@(x)length(x),xs_obs(1:end-1)) Kmem];
xs_OL=[xs_obs(1:end-1),{xs_obs{end}(1:Kmem)}];
space_inds=repmat({':'},1, dim-1);
U_obs_gap=cellfun(@(x)x(space_inds{:},1:Kmem),U_obs,'uni',0);

mx_max = min(cellfun(@(x) floor((length(x)-1)/2), xs_obs(1:end-1)));
mt_max=floor((Kmem-1)/2/1);

if use_cornerpt>0
    tauhat_ind = unique(min(use_cornerpt,length(U_obs)));
    [corners,sig_est] = findcornerpts(U_obs_gap{tauhat_ind},xs_OL,0);
    k_x = floor(mean(cellfun(@(x)x(end),corners(1:dim-1))));
    k_t = corners{dim}(end);
    [~,m_x,p_x] = get_phi_handle(min(dims(1:end-1)),'phi_class',phi_class{1},...
        'tau',tau_x,'tauhat',tauhat_x,'k',k_x,'maxd',max_dx);
    [~,m_t,p_t] = get_phi_handle(Kmem,'phi_class',phi_class{2},...
        'tau',tau_t,'tauhat',tauhat_t,'k',k_t,'maxd',max_dt);
end
m_x=max(min(m_x,mx_max),max_dx*2);
m_t=max(min(m_t,mt_max),max_dt*2);
ps = -[p_x p_t]; ms = [m_x m_t];
%%% rescale data (not used here)
scales = ones(1,n+dim); 

s_t = 1;
s_x = ceil(((Kmem-2*m_t)*prod(cellfun(@(x)length(x),xs_obs(1:end-1))-2*m_x)/max_rows)^(1/(dim-1)));
[Cfs_x,Cfs_t,p_x,p_t,sub_inds] = get_testfcn_weights(dimsgap,max_dx,max_dt,m_x,m_t,s_x,s_t,ps,phi_class);
Cfs_ffts = cell(dim-1,1);
[mm,nn] = size(Cfs_x);
for k=1:dim-1
    Cfs_ffts{k} = [zeros(mm,dims(k)-nn) (m_x*dx*scales(n+k)).^(-(0:mm-1)').*Cfs_x];
    Cfs_ffts{k} = fft(Cfs_ffts{k},[],2);
end
Cfs_t_scaled=(m_t*dt*scales(n+dim)).^(-(0:size(Cfs_t,1)-1)').*Cfs_t;

%% get cell array of weak-formed library over first Kmem snapshots

Theta_cell = repmat({zeros([cellfun(@(x)length(x),sub_inds(1:dim-1)) Kmem])},size(lib_list,1),1);
n = length(U_obs_gap);
libtree = liblisttree(lib_list,n);
for tt=1:Kmem
    U=cellfun(@(x)x(space_inds{:},tt),U_obs_gap,'uni',0);
    Theta_cell = get_lib_OL(...
    U,Theta_cell,lib_list,libtree,scales,sub_inds(1:end-1),Cfs_ffts,tt);
end

%% build initial linear system

[G,b,~] = build_linear_system(Theta_cell,Cfs_t_scaled,lib_list,lhs_ind,[]);

%% set iterative method for solving sparse reg. prob

H = get_sparse_update(SRmeth);

%% set initial conditions

W_0=get_initial_weights(G,b,gamma,ICmeth);

T=min(length(xs_obs{end})-Kmem,maxOLits);
W = W_0;
lambda=lambda_init;
resid = b-G*W;
sparsity = W~=0;
M=1./max(vecnorm(G)',eps);
LB = max(norm(b)./max(vecnorm(G)',eps),1);
obj_fcns = zeros(1,T+1);
obj_fcns(1) = [norm(resid(:))^2+sum(sum((W_0~=0).*(max(1,LB)./M*lambda).^2))];

%% if available, track coefficient error

if contains(pde_name,{'wave','2D','VC','du3'})
    load([data_dr,pde_name],'c');
    axi=repmat({W_0*0},1,length(xs_obs{end}));
    speed_tags =ismember(tags_pde_G,{'u^{1}_{xx}','u^{1}_{yy}'});
    cvec=c(xs_obs{end});
    for j=1:length(xs_obs{end})
        axi{j}(ismember(tags_pde_G,{'u^{1}_{xx}','u^{1}_{yy}','u^{3}_{}'}))=[cvec(j) cvec(j) -1];
    end
    axi0={sum(Cfs_t(1,:)/norm(Cfs_t(1,:),1).*cell2mat(axi(1:Kmem)),2)};
else
    axi=tags2axi(true_nz_weights,lib_list(~ismember(lib_list,lhs,'rows'),:));
    axi0=axi;
end

Ws = {W_0};

ET_offline = toc;
ETs = ET_offline+ET_load_data;

fprintf(1,'offline precomputation: %4.4f sec \n',ET_offline);