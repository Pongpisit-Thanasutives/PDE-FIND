%% Kuramoto-Sivashinsky equation and chaos
nn = 1024;
steps = 251;

maxx = 32*pi;
minx = maxx/nn;
maxx = maxx+minx;

dom = [minx maxx]; x = chebfun('x',dom); tspan = linspace(0,100,steps);
S = spinop(dom, tspan);

%% true
% S.lin = @(u) - diff(u,2) - diff(u,4);
% S.nonlin = @(u) - 0.5*diff(u.^2);

%% com=3
% S.lin = @(u) - 0.99514305*diff(u,2) - 0.99660169*diff(u,4);
% S.nonlin = @(u) - (0.99681542/2)*diff(u.^2);

%% com=4
S.lin = @(u) - (0.995202555)*diff(u,2) + (6.38896462e-04)*diff(u,3) - (0.996653643)*diff(u,4);
S.nonlin = @(u) - (0.996894189/2)*diff(u.^2);

%% com=5
% S.lin = @(u) + (1.11347422e-03)*u - (9.95333392e-01)*diff(u,2) - (9.96563189e-01)*diff(u,4);
% S.nonlin = @(u) - (0.997362570/2)*diff(u.^2) - (1.17039851e-04)*(...);

% ic = cos(x/16).*(1+sin(x/16));
% ic = cos(0.062976345*x).*(0.9813761+sin(0.062976345*x)); % pysr_params
ic = cos(0.06263812*x).*(0.9839311+sin(0.06286727*x)); % jaxfit
S.init = ic;

u = spin(S,nn,1e-4,'plot','off'); % 1e-3, use 1e-4 for ground truth generation
usol = zeros(nn,steps);
for i = 1:steps
    usol(:,i) = u{i}.values;
end

x = linspace(minx,maxx,nn+1);
x = x(1:end-1);
t = tspan;

% save('./KS_sim/ks_chaotic_simV2_adjustedIC_com4.mat','t','x','usol');
pcolor(t,x,usol); shading interp, axis tight, colormap(jet);