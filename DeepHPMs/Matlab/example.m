% Define the domain and initial condition
dom = [-8 8];
x = chebfun('x', dom);
u0 = exp(-(x+2).^2);

% Define the differential operator and boundary conditions
pdefun = @(x, t, u, ux) [diff(u, 2) + u.*diff(u, 1) + u.*diff(u, 2); 0; 0];
bcfun = @(ya, yb) [ya(1); yb(1)];

% Define the times at which the solution should be returned
tspan = linspace(0, 10, 101);

% Solve the PDE using pdepe
sol = pdepe(1, pdefun, u0, bcfun, tspan);

% Plot the solution
u = sol.';
pcolor(tspan, x, u);
shading interp;
colormap(jet);
