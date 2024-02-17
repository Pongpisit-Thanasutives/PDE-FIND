%% Burgers equation and chaos
nn = 256;
steps = 100;

dom = [-8 8]; x = chebfun('x',dom); tspan = linspace(0,10,steps+1);
S = spinop(dom, tspan);

% % GROUND
% S.lin = @(u) + 0.1*diff(u,2);
% S.nonlin = @(u) - (1/2)*diff(u.^2); % spin cannot parse "u.*diff(u)"

% % none | com=2 | 0.10333732, -1.01090519 | u_11+u*u_1
% % kalman | com=2 | 0.10599135, -0.99765862 | u_11+u*u_1
S.lin = @(u) + 0.10333732*diff(u,2);
S.nonlin = @(u) - (1.01090519/2)*diff(u.^2);
% S.lin = @(u) + 0.10599135*diff(u,2);
% S.nonlin = @(u) - (0.99765862/2)*diff(u.^2);

% % none, found | com=3 | 0.10181422, -1.01283594,  0.00382179 | u_11+u*u_1+u*u_11
% % none, brute | com=3 | -0.01370488, 0.10196836, -0.97529619 | u_1+u_11+u*u_1
% % kalman, both | com=3 | 0.10483621, -1.08522265,  0.18429072 | u_11+u*u_1+u*u*u_1
% S.lin = @(u) - 0.01370488*diff(u,1) + 0.10196836*diff(u,2);
% S.nonlin = @(u) - (0.97529619/2)*diff(u.^2);
% S.lin = @(u) + 0.10483621*diff(u,2);
% S.nonlin = @(u) - (1.08522265/2)*diff(u.^2) + (0.18429072/3)*diff(u.^3);

% S.init = -sin(pi*x/8);
% S.init = exp(-(x+2).^2);
S.init = exp(-1.0023574*(x+2.0152135).^2); % none
% S.init = exp(-1.0219623*(x+1.9614727).^2); % kalman

tic, u = spin(S,nn,1e-4,'plot','off'); % acc = 4-6 enough
usol = zeros(nn,steps+1);
for i = 1:steps+1
    usol(:,i) = u{i}.values;
end

x = linspace(-8,8,nn+1);
x = x(1:end-1);
t = tspan;

pcolor(x,t,transpose(real(usol))); shading interp, axis tight, colormap(jet);
% save('burgers_kalman_est_com3.mat','t','x','usol')