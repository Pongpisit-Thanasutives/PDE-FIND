%% Korteweg-de Vries equation (two-soliton solution)
nn = 512;
steps = 500;
L = 20;
dom = [-L L]; x = chebfun('x',dom); tspan = linspace(0,40,steps+1);
S = spinop(dom, tspan);

% com=2 from WF
% v = -0.9103100649234885;
% S.lin = @(u) -0.9037167644398507*diff(u,3);
% S.nonlin = @(u) +v*0.5*diff(u.^2);
% com=2 from CWF
S.lin = @(u) -0.9950437*diff(u,3);
S.nonlin = @(u) -0.50073689*diff(u.^2)-1e-4;

% com=3 from WF
% The solution blew up. Require a too small time-step. (BUT Dedalus3 works)
% v = -0.9601505537264142;
% S.lin = @(u) -0.889691174996562*diff(u,3);
% S.nonlin = @(u) +v*0.5*diff(u.^2)-0.05551274586939359*(...);
% com=3 from CWF (tested)
% S.lin = @(u) -0.98987647*diff(u,3);
% S.nonlin = @(u) -0.50145825*diff(u.^2)-0.00490339*diff(u.^2,3)-1e-4;

% S.init = -sin(pi*x/20); % ground
% S.init = sin(-0.15645656*x); % pysr_params
S.init = sin(-0.1564577903993788*x); % jaxfit

u = spin(S,nn,1e-5,'plot','off'); % WF = 1e-4, CWF = 1e-5
usol = zeros(nn,steps+1);
for i = 1:steps+1
    usol(:,i) = u{i}.values;
end

x = linspace(-L,L,nn+1);
x = x(1:end-1);
t = tspan;
pcolor(t,x,usol); shading interp, axis tight, colormap(jet);
save('./KdV_sim/kdv_sim_CWF_BIAS_adjustedIC_com2.mat','t','x','usol');
% save('./KdV_sim/kdv_sim_CWF_adjustedIC_com3.mat','t','x','usol');