
function [Cfs_x,Cfs_t,p_x,p_t,sub_inds] = get_testfcn_weights(dims,max_dx,max_dt,m_x,m_t,s_x,s_t,tols,phi_class)

    dim = length(dims);

    if or(2*m_x+1>min(dims(1:end-1)),2*m_t+1>min(dims(end)))
        disp('ERROR: Test function not compactly supported')
        return
    end

    sub_inds = cell(1,dim);
    mm = [repmat(m_x,1,dim-1) m_t];
    ss = [repmat(s_x,1,dim-1) s_t];
    for j=1:dim
        N = dims(j);
        m = mm(j);
        s = ss(j);
        sub_inds{j} = 1:s:N-2*m;
    end

    if length(phi_class)==1
        phi_class = repmat(phi_class,1,2);
    end
    
    if dim>1
        [Cfs_x,p_x] = phi_int_weights(m_x,max_dx,tols(1),phi_class{1});
    else
        Cfs_x = 1; p_x = 0;
    end
    [Cfs_t,p_t] = phi_int_weights(m_t,max_dt,tols(2),phi_class{2});

end

function [Cfs,p] = phi_int_weights(m,d,tol,phi_class)
    if isequal(class(phi_class),'double')
        if phi_class == 1
            if tol<0
                p = -tol;
            else
                p = ceil(max(log(tol)/log((2*m-1)/m^2),d+1));     % choose p so that the penultimate grid-point has value eps (for decay)
            end
            t = (0:m)/m;
            t_L = zeros(d+1,m+1);                             % store (1+t)^q, (1-t)^q
            t_R = zeros(d+1,m+1);                                  
            for j=1:m
                t_L(:,j)  = (1+t(j)).^(fliplr(p-d:p))';         
                t_R(:,j)  = (1-t(j)).^(fliplr(p-d:p))';
            end

            ps = ones(d+1,1);                                  % derivative coefficients
            for q=1:d
                ps(q+1) = (p-q+1)*ps(q);
            end
            t_L = ps.*t_L;
            t_R = ((-1).^(0:d)'.*ps).*t_R;

            Cfs = zeros(d+1,2*m+1);                            % Values of derivatives at grid points
            Cfs(1,:) = [fliplr(t_L(1,:).*t_R(1,:)) t_L(1,2:end).*t_R(1,2:end)];
            P = fliplr(pascal(d+1));    
            for k=1:d
                binoms = diag(P,d-k);
                Cfs_temp = zeros(1,m+1);
                for j=1:k+1
                    Cfs_temp = Cfs_temp + binoms(j)*t_L(k+2-j,:).*t_R(j,:);
                end
                Cfs(k+1,:) = [(-1)^k*fliplr(Cfs_temp) Cfs_temp(2:end)];
            end
        elseif phi_class == 2
            if and(tol<0,m>0)
                p = -tol;
                a = m*p;
            elseif and(tol > 0, m>0)
                a = sqrt(-2*log(tol));
                p = (1-1/m)/a;
            elseif and(tol > 0, m <= 0)
                p = -m;
                a = sqrt(-2*log(tol));
                m = ceil(1+a*p);
                p = p / m;
            end
    %         x = linspace(-a,a,2*m-1);
    %         dx = x(2)-x(1);
    %         x=[x(1)-dx x x(end)+dx];
            x = linspace(-a,a,2*m+1);

            Cfs = ones(d+1,2*m+1);
            Cfs(2,:) = x;

            for k=3:d+1
                Hnp = (k-2)*Cfs(k-2,:);
                Cfs(k,:) = x.*Cfs(k-1,:)-Hnp;
            end

            e = exp(-x.^2/2);
            s = (-p*m).^((0:d)');

            Cfs = Cfs.*(s*e);
        end
    elseif isequal(class(phi_class),'function_handle')
        x = linspace(-1,1,2*m+1);
        Cfs = zeros(d+1,2*m+1);
        syms y;
        f = @(y)phi_class(y);
        for j=1:d+1
            Df = matlabFunction(diff(f(y),j-1));
            Cfs(j,:) = fillmissing(Df(x),'constant',Df(eps));
            inds = find(isinf(abs(Cfs(j,:))));
            for k=1:length(inds)
                Cfs(j,inds(k)) = Df(x(inds(k))+eps);
            end
        end
        p = 0;
    end
    Cfs = Cfs./norm(Cfs(1,:));
end