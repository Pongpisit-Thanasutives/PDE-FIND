import pymc as pm
import numpy as np

def horseshoe_prior(X, tau0=1e-2, v=4, s=2):
    '''
    Regularizing horseshoe prior as introduced by Piironen & Vehtari
    https://arxiv.org/pdf/1707.01694.pdf
    
    name: variable name
    X: X (2-d array)
    y: y (for setting pseudo-variance)
    m: expected number of relevant features (must be < total N)
    v: regularizing student-t df
    s: regularizing student-t sd
    '''
    
    half_v = v/2
    n = X.shape[0]
    M = X.shape[1]
    
    tau_t = pm.HalfCauchy(f"tauT", beta = 1)
    tau = tau0*tau_t

    c2_t = pm.InverseGamma(f"c2", half_v, half_v)
    c2 = np.power(s,2) * c2_t

    lambda_m = pm.HalfCauchy(f"lambdaM", beta = 1, shape = M)
    lambda_t = (pm.math.sqrt(c2)*lambda_m) / pm.math.sqrt(c2 + pm.math.sqr(tau) * pm.math.sqr(lambda_m))
    

    beta_t = pm.Normal(f"coeffs", mu=0, sigma = 1, shape= M)
    beta = tau * lambda_t * beta_t
    
    return beta
