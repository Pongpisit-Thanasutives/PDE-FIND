ode_params=[];
tspan=0:0.01:1000;
x0=[];
tol_ode=10^-13;

[weights,x,t,ode_name,ode_params,rhs] = gen_data(ode_num,ode_params,tspan,x0,tol_ode);

U_obs = mat2cell(x',ones(size(x,2),1),size(x,1));
xs_obs = {0,t};
lhs = [eye(length(U_obs)) zeros(length(U_obs),1) ones(length(U_obs),1)];
pde_name = '';

true_nz_weights = weights;
for i=1:length(weights)
	true_nz_weights{i} = [true_nz_weights{i}(:,1:length(weights)) zeros(size(true_nz_weights{i},1),2) true_nz_weights{i}(:,end)];
end

function [weights,x,t,ode_name,ode_params,rhs] = gen_data(ode_num,ode_params,tspan,x0,tol_ode)

    if isempty(tol_ode)
        tol_ode=10^-12;
    end

% Simulate ODE
    ode_names = {'Linear','Logistic_Growth','Van_der_Pol','Duffing','Lotka_Volterra','Lorenz','Rossler','rational','Oregonator','Hindmarsh-Rose','Pendulum','custom'};
    ode_name = ode_names{ode_num};
    if strcmp(ode_name, 'Linear')%1
        % ll-dim
        ll = 2;
        x = zeros(ll,1);y= zeros(ll,1); x([2 end])=[1 -1]; y([2 end]) = [-1 1];
        A = toeplitz(y,x');
        A(logical(eye(size(A,1))))= -0.2;
        if isempty(ode_params)
            ode_params = {A};
        end
        if isempty(x0)
            x0 = zeros(ll,1); 
            x0(1) = 10; 
        end
        if isempty(tspan)
            tspan = 0:0.025:40;
        end
        % 2-dim
    %     ode_params = {[[-0.5 2];[-2 -0.5]]};x0 = [2;0]; tspan = 0:0.01:15;
    %     ode_params = {[[-0.1 2];[-2 -0.1]]};x0 = [2;0]; tspan = 0:0.01:20;

    elseif strcmp(ode_name, 'Logistic_Growth')%2
        if isempty(ode_params)
            ode_params = {2};
        end
        if isempty(x0)
            x0=0.01;
        end
        if isempty(tspan)
            tspan =  0:0.01:10;
        end


    elseif strcmp(ode_name, 'Van_der_Pol')%3
        if isempty(ode_params)
            ode_params = {4};
        end
        if isempty(x0)
            x0=[0;1];
        end
        if isempty(tspan)
            tspan = 0:0.01:30;    
        end
        %ode_params = {0.5}; x0 = [-2;1]; tspan = 0:0.01:10;    

    elseif strcmp(ode_name, 'Duffing')%4
    %        mu = 0.2; ode_params = {0.1, 0.2,1}; x0 = [-2;-2]; tspan = 0:0.01:25;
    %        ode_params = {0.2, 0.05,1}; x0 = [0;2]; tspan = 0:0.005:25;
        if isempty(ode_params)
            mu = 0.2; ode_params = {mu, mu, 1}; 
        end
        if isempty(x0)
            x0=[0;2];
        end
        if isempty(tspan)
            tspan = 0:0.01:30;    
        end
        
        
    elseif strcmp(ode_name, 'Lotka_Volterra')%5
    %    alpha = 1; beta = 10; ode_params = {alpha,beta,beta,2*alpha}; x0 = [1;1]; tspan = 0:0.01:25;
    %    alpha = 2/3; beta = 4/3; ode_params = {alpha,beta,1,1}; x0 = [10;10]; tspan = 0:0.02:125;
    %    alpha = 1; beta = 0.5; ode_params = {alpha,beta,beta,2*alpha}; x0 = [1;2]; tspan = 0:0.01:10;
        if isempty(ode_params)
            alpha = 3; beta = 1; ode_params = {alpha,beta,beta,2*alpha};     
        end
        if isempty(x0)
            x0=[1;1];
        end
        if isempty(tspan)
            tspan = 0:0.01:10;
        end
        
        
    elseif strcmp(ode_name, 'Lorenz')%6
        %x0 = [5 5 25]';%x0 = [rand(2,1)*30-15;rand*30+10]
        if isempty(ode_params)
            ode_params = {10, 8/3, 28};     
        end
        if isempty(x0)
            x0=[-8 7 27]';
        end
        if isempty(tspan)
            tspan = 0:0.001:10;
        end
        

    elseif strcmp(ode_name, 'Rossler')%7
    %        ode_params = {0.2, 0.2, 5.7}; x0 = [3 5 0]'; tspan = 0:.005:25;
        if isempty(ode_params)
            ode_params = {0.2, 0.2, 5.7}; 
        end
        if isempty(x0)
            x0 = [3 5 0]';
        end
        if isempty(tspan)
            tspan = 0:.01:25;
        end
        

    elseif strcmp(ode_name, 'rational')%8
        if isempty(ode_params)
            ode_params = {8,1}; 
        end
        if isempty(x0)
            x0 = -5;
        end
        if isempty(tspan)
            tspan = 0:.005:20;
        end
        

    elseif strcmp(ode_name, 'Oregonator')%9
        if isempty(ode_params)
            ode_params = {0.5, 5, 5, 0.5, 0.5, 0.95};
        end
        if isempty(x0)
            x0 = [1 0.2 50]';    
        end
        if isempty(tspan)
            tspan = 0:.005:50;
        end
        
        
    elseif strcmp(ode_name, 'Hindmarsh-Rose')%10
        if isempty(ode_params)
            a = 1; b = 3; c = 1; 
            d = 5; r = 10^-3; s = 4;
            xR = -3.19/4; I = 0;
            ts = 10;
            ode_params = {a,b,c,d,r,s,xR,I,ts};
        end
        if isempty(x0)
            x0 = [-1.3095   -7.5901   -0.2020]';
        end
        if isempty(tspan)
            tspan = 0:.001:100;
        end
        

    elseif strcmp(ode_name, 'Pendulum')%11
        if isempty(ode_params)
            ode_params = {1};
        end
        if isempty(x0)
            x0 = [pi-1*pi/16 0]';
        end
        if isempty(tspan)
            tspan = 0:0.1:50;
        end
        

    elseif strcmp(ode_name, 'SIRD')%12
    elseif strcmp(ode_name, 'custom')%13
        % nonlinear Schrodinger
    %     ode_func = @(x) [0.1*x(2).*(x(1).^2+x(2).^2);-x(1).*(x(1).^2+x(2).^2)];
    %     weights = {[[2 1 0.1];[0 3 0.1]], [[3 0 -1];[1 2 -1]]};
        % cubic oscillator
        if isempty(ode_params)
            ode_func = @(x) [-0.1*x(1).^3+2*x(2).^3;-2*x(1).^3-0.1*x(2).^3];
            weights = {[[3 0 -0.1];[0 3 2]], [[3 0 -2];[0 3 -0.1]]};
            ode_params = {ode_func,weights};
        end
        if isempty(x0)
            x0 = [1;0];
        end
        if isempty(tspan)
            tspan = 0:.01:20;
        end
        

    end
    [weights,t,x,rhs] = sim_ode(x0,tspan,tol_ode,ode_name,ode_params,'o');

end

function [weights,t,x,rhs] = sim_ode(x0,tspan,tol_ode,ode_name,params,odemethod)

if strcmp(ode_name,'Logistic_Growth')
    pow = params{1};
    rhs = @(x) x-x.^pow; 
    weights = {[[1 1];[pow -1]]};

elseif strcmp(ode_name,'Linear')
    A = params{1};
    rhs = @(x) A*x; 
    weights = {};
    for dim = 1:size(A,1)
        weights{dim} = [eye(size(A,1)) A(dim,:)'];
    end
    
elseif strcmp(ode_name,'Duffing')
    mu = params{1}; alpha= params{2}; beta=params{3};
    rhs = @(x) duff(x,mu,alpha,beta);
    weights = {[0 1 1], [[1 0 -alpha];[0 1 -mu];[3 0 -beta]]};

elseif strcmp(ode_name,'Lotka_Volterra')
    alpha = params{1}; beta = params{2}; delta = params{3}; gamma = params{4};
    rhs = @(x) LoVo(x,alpha,beta,delta, gamma);
    weights = {[[1 0 alpha];[1 1 -beta]],[[0 1 -gamma];[1 1 delta]]};

elseif strcmp(ode_name,'Van_der_Pol')
    mu = params{1};
    rhs = @(x) vanderpol(x,mu);
    weights = {[0 1 1], [[1 0 -1];[0 1 mu];[2 1 -mu]]};
    
elseif strcmp(ode_name,'Lorenz')
    sigma = params{1}; beta = params{2}; rho = params{3};
    rhs = @(x) lorenz(x,sigma,beta,rho);    
    weights = {[[0 1 0 sigma];[1 0 0 -sigma]],...
                      [[1 0 0 rho];[1 0 1 -1];[0 1 0 -1]],...
                      [[1 1 0 1];[0 0 1 -beta]]};                  
elseif strcmp(ode_name,'Rossler')
    [a,b,c] = params{:};
    rhs = @(x) rossler(x,a,b,c);    
    weights = {[[0 1 0 -1];[0 0 1 -1]],...
                      [[1 0 0 1];[0 1 0 a]],...
                      [[0 0 0 b];[1 0 1 1];[0 0 1 -c]]};                  

elseif strcmp(ode_name,'Oregonator')
    k1 = params{1}; k2 = params{2}; k3 = params{3}; k4 = params{4}; k5 = params{5}; f = params{6};
    rhs = @(x) oregonator(x,k1,k2,k3,k4,k5,f);
    weights = {[[0 1 0 k1];...
                [1 1 0 -k2];...
                [1 0 0 k3];...
                [2 0 0 -2*k4]],...
               [[0 1 0 -k1];...
                [1 1 0 -k2];...
                [0 0 1 0.5*f*k5]],...
               [[1 0 0 2*k3];...
                [0 0 1 -k5]]};

elseif strcmp(ode_name,'Hindmarsh-Rose')
    a = params{1}; b = params{2}; c = params{3}; 
    d = params{4}; r = params{5}; s = params{6};
    xR = params{7}; I = params{8};
    ts = params{9};
    rhs = @(x) HindmarshRose(x,a,b,c,d,r,s,xR,I,ts);
    weights = {[[0 1 0 ts];...
                [3 0 0 -a*ts];...
                [2 0 0 b*ts];...
                [0 0 1 -1*ts];...
                [0 0 0 I*ts]],...
               [[0 0 0 c*ts];...
                [2 0 0 -d*ts];...
                [0 1 0 -1*ts]],...
               [[1 0 0 r*s*ts];...
                [0 0 0 -r*s*xR*ts];...
                [0 0 1 -r*ts]]};
            
elseif strcmp(ode_name,'Pendulum')
    g = params{1};
    rhs = @(x) pendulum(x,g);
    weights = {[0 1 1], [-1i 0 -g]};

elseif strcmp(ode_name,'rational')
    rhs = @(x) (params{2}-x)./(1+x.^2);    
    v = repmat([1;-1],params{1}/2+1,1);
    weights = {[(0:2:params{1})' v(1:params{1}/2+1)]};

elseif strcmp(ode_name,'custom')
    rhs = @(x) params{1}(x);    
    weights = params{2};
end

options = odeset('RelTol',tol_ode,'AbsTol',tol_ode*ones(1,length(x0)));
if odemethod=='s'
    [t,x]=ode15s(@(t,x)rhs(x),tspan,x0,options);  
else
    [t,x]=ode45(@(t,x)rhs(x),tspan,x0,options);
end
end


function dx = duff(x,mu,alpha,beta)
    dx(1) = x(2);
    dx(2) = -mu*x(2)-alpha*x(1)-beta*x(1).^3;
    dx = dx';
end

function dx = LoVo(x,alpha,beta,delta, gamma)
    dx(1) = alpha*x(1)-beta*x(1).*x(2);
    dx(2) = delta*x(1).*x(2)- gamma*x(2);
    dx = dx';
end

function dx = lorenz(x,sigma,beta,rho)
    dx(1) = sigma*(x(2)-x(1));
    dx(2) = x(1)*(rho-x(3))-x(2);
    dx(3) = x(1)*x(2)-beta*x(3);
    dx=dx';
end

function dx = rossler(x,a,b,c)
    dx(1) = -x(2)-x(3);
    dx(2) = x(1)+a*x(2);
    dx(3) = b+x(1).*x(3)-c*x(3);
    dx=dx';
end

function dx = vanderpol(x,mu)
    dx(1) = x(2);
    dx(2) = mu*x(2)-mu*x(1).^2.*x(2)-x(1);
    dx=dx';
end

function dx = oregonator(x,k1,k2,k3,k4,k5,f)
    dx(1) =  k1*x(2) - k2*x(1).*x(2) + k3*x(1) - 2*k4*x(1).^2;
    dx(2) = -k1*x(2) - k2*x(1).*x(2) + 0.5*f*k5*x(3);
    dx(3) =  2*k3*x(1) - k5*x(3);
    dx=dx';
end

function dx = HindmarshRose(x,a,b,c,d,r,s,xR,I,ts)
    dx(1) =  ts*(x(2) - a*x(1).^3+b*x(1).^2 -x(3)+I);
    dx(2) = ts*(c-d*x(1).^2 - x(2));
    dx(3) =  ts*(r*(s*(x(1)-xR)-x(3)));
    dx=dx';
end

function dx = pendulum(x,g)
    dx(1) = x(2);
    dx(2) = -g*sin(x(1));
    dx=dx';
end
