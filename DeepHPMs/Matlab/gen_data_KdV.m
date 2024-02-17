%% Korteweg-de Vries equation (two-soliton solution)
nn = 512;
steps = 500;
v = -1.0;
L = 20;
dom = [-L L]; x = chebfun('x',dom); tspan = linspace(0,40,steps+1);

S = spinop(dom, tspan);
S.lin = @(u) -1.0*diff(u,3);
S.nonlin = @(u) +v*0.5*diff(u.^2); % spin cannot parse "u.*diff(u)"
S.init = -sin(pi*x/20);

u = spin(S,nn,1e-5,'plot','off');
usol = zeros(nn,steps+1);
for i = 1:steps+1
    usol(:,i) = u{i}.values;
end

x = linspace(-L,L,nn+1);
x = x(1:end-1);
t = tspan;
pcolor(t,x,usol); shading interp, axis tight, colormap(jet);
save('/Users/pongpisit/Desktop/research/PDE-FIND/Datasets/KdV_sine_rep_acc_big.mat','t','x','usol');