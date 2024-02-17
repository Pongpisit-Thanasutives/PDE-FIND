%% Korteweg-de Vries equation (two-soliton solution)
nn = 512;
steps = 400;
v = -6.0;
L = 30;
dom = [-L L]; x = chebfun('x',dom); tspan = linspace(0,20,steps+1);

S = spinop(dom, tspan);
S.lin = @(u) -1.0*diff(u,3);
S.nonlin = @(u) +v*0.5*diff(u.^2); % spin cannot parse "u.*diff(u)"
S.init = kdv_exact(x, 1, -20)+kdv_exact(x, 0.5, 5);

u = spin(S,nn,1e-5,'plot','off');
usol = zeros(nn,steps+1);
for i = 1:steps+1
    usol(:,i) = u{i}.values;
end

x = linspace(-L,L,nn+1);
x = x(1:end-1);
t = tspan;
pcolor(t,x,usol); shading interp, axis tight, colormap(jet);
save('/Users/pongpisit/Desktop/research/PDE-FIND/Datasets/KdV_rudy_acc_big.mat','t','x','usol');