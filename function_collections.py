import numpy as np

def ssr2llf(ssr, nobs):
    nobs2 = nobs / 2.0
    llf = -nobs2 * np.log(2 * np.pi) - nobs2 * np.log(ssr / nobs) - nobs2
    return llf

def log_like_value(prediction, ground):
    nobs = float(ground.shape[0])
    ssr = np.sum(np.abs(ground - prediction)**2)
    return ssr2llf(ssr, nobs)

def BIC_AIC(prediction, ground, nparams, reg_func = lambda x: x):
    nparams = reg_func(nparams)
    llf = log_like_value(prediction, ground)
    return -2*llf + np.log(ground.shape[0])*nparams, -2*llf + 2*nparams

def AIC_Loss(A,b,x,epsilon=1e-5):
    N = A.shape[0]
    k = np.count_nonzero(x)
    # Rudy et al., 2019
    rss = ((b-A.dot(x))**2).sum()
    aic = N*np.log(rss/N+epsilon) + 2*k + (2*k**2+2*k)/(N-k-1)
    return aic

def PDE_FIND_Loss(As,bs,x,epsilon=1e-5,const_coeff=True,cv=0,ic_type="bic",version=0):
    # D: Number of candidates | m: either len(t) or len(x) (temporal or spatial group)
    D,m = x.shape
    # n: Number of horizon
    n,_ = As[0].shape
    N = n*m

    rss = [np.linalg.norm(bs[j] - As[j].dot(x[:,j].reshape(D,1)))**2 for j in range(m)]
    rss = np.sum(rss)
    llf = ssr2llf(rss, N)

    k = np.count_nonzero(x)
    if const_coeff:
        k = np.count_nonzero(x)/m
    k = k + k*cv

    if ic_type == "aic":
        aic1 = N*np.log(rss/N+epsilon) + 2*k + (2*k**2+2*k)/(N-k-1)
        aic2 = -2*llf + 2*k + (2*k**2+2*k)/(N-k-1)
        if version > 0:
            return aic2
        return aic1
    elif ic_type == "bic":
        bic1 = N*np.log(rss/N+epsilon) + np.log(N)*k
        bic2 = -2*llf + np.log(N)*k
        if version > 0:
            return bic2
        return bic1

def percent_coeff(pred, ground):
    return 100*np.abs(pred-ground)/np.abs(ground)

def msemse(AA, BB):
    return ((AA-BB)**2).mean()

def colvec(vec):
    return vec.reshape(-1, 1)
