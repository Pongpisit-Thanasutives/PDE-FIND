%% Burgers equation and chaos
nx = 256;
L = 1;
dx = 2*L/(nx-1);
Lx = -L;
Rx = L+dx;
dom = [Lx Rx];

steps = 100;
tspan = linspace(0,1,steps+1);

x = chebfun('x',dom);
S = spinop(dom, tspan);

% c1 = 0.01/pi;
% c2 = 1;
% c1 = 0.00374609;
% c2 = 0.92514145;
% S.lin = @(u) + c1*diff(u,2);
% S.nonlin = @(u) - (c2/2)*diff(u.^2);
% com=3 From WF
% S.lin = @(u) + 0.34931585*u;
% S.nonlin = @(u) - (0.91928852/2)*diff(u.^2) + 0.00636898*((u.^2).*diff(u,2));

% com=2 From CWF
% S.lin = @(u) + 0.00288483*diff(u,2);
% S.nonlin = @(u) - 0.44707268*diff(u.^2);
% com=3 From CWF
S.lin = @(u) - 0.10564249*u + 0.0024246*diff(u,2);
S.nonlin = @(u) - 0.4259986*diff(u.^2);

% S.init = -sin(pi*x);
S.init = -sin(3.1763942*x);

tic, u = spin(S,nx,1e-5,'plot','off'); % 1e-4 for WF, 1e-5 for CWF
usol = zeros(nx,steps);
for i = 1:steps
    usol(:,i) = u{i}.values;
end

x = linspace(Lx, Rx, nx+1);
x = x(1:end-1);
t = tspan(1:end-1);

pcolor(x,t,transpose(real(usol))); shading interp, axis tight, colormap(jet);
% save('./Burgers_shock_sim/burgers_shock_reproduce.mat','t','x','usol');
% save('./Burgers_shock_sim/burgers_shock_sim_com2.mat','t','x','usol');
% save('./Burgers_shock_sim/burgers_shock_sim_com3.mat','t','x','usol');
% save('./Burgers_shock_sim/burgers_shock_sim_CWF_com2.mat','t','x','usol');
save('./Burgers_shock_sim/burgers_shock_sim_CWF_com3.mat','t','x','usol');