function [U,noise,snr,sigma] = add_noise(U_exact,stdv,sigma_NR,noise_dist,noise_alg)

    dims = size(U_exact);
    if noise_dist == 0 % white noise
        sigma = sigma_NR*sqrt(stdv);
        noise = normrnd(0,sigma,dims);
    elseif noise_dist == 1 % uniform noise
        sigma = (3*sigma_NR^2*stdv)^(1/2);
        noise = 2*sigma*rand(dims)-sigma;
    end    
    if noise_alg == 0 % additive
        U = U_exact + noise;
    elseif noise_alg == 1 % multiplicative
        U = U_exact.*(1 + noise);
    end
    snr = norm(U(:)-U_exact(:))/norm(U_exact(:));
    
end
