function usol = sim_burgers_func(com, coeff)
    nn = 256;
    steps = 100;

    dom = [-8 8]; x = chebfun('x',dom); tspan = linspace(0,10,steps+1);
    S = spinop(dom, tspan);

    if com ==2
        S.lin = @(u) + coeff(1)*diff(u,2);
        S.nonlin = @(u) + (coeff(2)/2)*diff(u.^2);
    elseif com == 3
        S.lin = @(u) + coeff(1)*diff(u,1) + coeff(2)*diff(u,2);
        S.nonlin = @(u) + (coeff(2)/2)*diff(u.^2);
    end

    % S.init = exp(-(x+2).^2);
    rps = [-4.15259791, -4.97142427e-01, 2.95012304e-03, 1.00797403];
    S.init = rps(3)+rps(4)*exp(rps(1)*((rps(2)*x-1).^2));

    tic, u = spin(S,nn,1e-3,'plot','off'); % acc = 3-6 enough
    usol = zeros(nn,steps+1);
    for i = 1:steps+1
        usol(:,i) = u{i}.values;
    end

    x = linspace(-8,8,nn+1);
    x = x(1:end-1);
    t = tspan;
    usol = real(usol);

    pcolor(x,t,transpose(usol)); shading interp, axis tight, colormap(jet);
end
