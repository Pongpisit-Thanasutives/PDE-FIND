import numpy as np
import statsmodels.api as sm

def find_modes(arr):
    assert len(arr) > 0
    uniques, counts = np.unique(arr, return_counts=True)
    max_count = np.max(counts)
    return uniques[counts == max_count], list(counts).index(max_count)

def decreasing_values_indices(arr):
    mini = max(arr)+1; out = []
    for i, e in enumerate(arr):
        if e < mini:
            mini = e
            out.append(i)
    return np.array(out)

def log_like_value(prediction, ground): 
    nobs = float(ground.shape[0])
    ssr = np.sum(np.abs(ground - prediction)**2)
    def ssr2llf(ssr, nobs):
        # keep the constant term nobs2
        nobs2 = nobs / 2.0
        llf = -nobs2 * np.log(2 * np.pi) - nobs2 * np.log(ssr / nobs) - nobs2
        return llf
    return ssr2llf(ssr, nobs)

def BIC_AIC(prediction, ground, nparams, reg_func=lambda x:x):
    nparams = reg_func(nparams)
    llf = log_like_value(prediction, ground)
    return -2*llf + np.log(ground.shape[0])*nparams, -2*llf + 2*nparams

def baye_uncertainties(best_subsets, dataset, u_type='var', take_sqrt=True, ridge_lambda=0, unbiased=False):
    # if you want u_type='std', then call u_type='var' and take_sqrt=True
    XX, yy = dataset
    assert u_type == 'var' or 'cv' in u_type
    assert len(XX) == len(yy)
    yy = yy.reshape(-1, 1)

    post_means = np.zeros((XX.shape[-1], len(best_subsets)))
    bics = []
    uns = []
    for k, efi in enumerate(best_subsets):
        com = len(efi)
        Phi = XX[:, list(efi)]
        w = np.linalg.lstsq(Phi, yy, rcond=None)[0]
        err = yy-Phi@w
        # By MLE, we have variance_y written as follows:
        variance_y = np.mean(err**2)
        if unbiased: variance_y = variance_y*len(err)/(len(err)-com)
        w = w[np.abs(w)>0].reshape((com, 1))

        # prior_mean = np.zeros((com, 1))
        prior_mean = w
        prior_cov = np.identity(com)
        if ridge_lambda > 0: prior_cov = (variance_y/ridge_lambda)*prior_cov
        prior_cov_inv = np.linalg.inv(prior_cov)

        posterior_cov = variance_y*np.linalg.inv(variance_y * prior_cov_inv + Phi.T@Phi)
        posterior_mean = posterior_cov@(prior_cov_inv@prior_mean + (Phi.T@yy)/variance_y)
        post_means[:, k:k+1][list(efi)] = posterior_mean

        # collecting bics
        bics.append(BIC_AIC(Phi@posterior_mean, yy, com)[0])
        # collecting uns
        posterior_variance = np.diag(posterior_cov)
        if take_sqrt:
            posterior_variance = np.sqrt(posterior_variance)
        mm = posterior_mean
        ss = posterior_variance
        if u_type == 'var':
            uns.append(ss.sum())
        elif 'cv' in u_type:
            code = u_type.replace('cv', '')
            if len(code) == 0: order = 1
            else: order = int(u_type.replace('cv', ''))
            mm = np.linalg.norm(mm[:, 0], ord=order)
            ss = np.linalg.norm(ss, ord=order)
            uns.append(ss/mm)

    uns = np.array(uns)
    uns = uns/min(uns)
    return post_means, bics, uns

def baye_uncertainties_aic(best_subsets, dataset, u_type='var', take_sqrt=True, ridge_lambda=0):
    # if you want u_type='std', then call u_type='var' and take_sqrt=True
    XX, yy = dataset
    assert u_type == 'var' or 'cv' in u_type
    assert len(XX) == len(yy)
    yy = yy.reshape(-1, 1)

    post_means = np.zeros((XX.shape[-1], len(best_subsets)))
    bics = []
    uns = []
    for k, efi in enumerate(best_subsets):
        com = len(efi)
        Phi = XX[:, list(efi)]
        w = np.linalg.lstsq(Phi, yy, rcond=None)[0]
        err = yy-Phi@w
        # By MLE, we have variance_y written as follows:
        variance_y = np.mean(err**2)
        w = w[np.abs(w)>0].reshape((com, 1))

        # prior_mean = np.zeros((com, 1))
        prior_mean = w
        prior_cov = np.identity(com)
        if ridge_lambda > 0: prior_cov = (variance_y/ridge_lambda)*prior_cov
        prior_cov_inv = np.linalg.inv(prior_cov)

        posterior_cov = variance_y*np.linalg.inv(variance_y * prior_cov_inv + Phi.T@Phi)
        posterior_mean = posterior_cov@(prior_cov_inv@prior_mean + (Phi.T@yy)/variance_y)
        post_means[:, k:k+1][list(efi)] = posterior_mean

        # collecting bics
        bics.append(BIC_AIC(Phi@posterior_mean, yy, com)[1])
        # collecting uns
        posterior_variance = np.diag(posterior_cov)
        if take_sqrt:
            posterior_variance = np.sqrt(posterior_variance)
        mm = posterior_mean
        ss = posterior_variance
        if u_type == 'var':
            uns.append(ss.sum())
        elif 'cv' in u_type:
            code = u_type.replace('cv', '')
            if len(code) == 0: order = 1
            else: order = int(u_type.replace('cv', ''))
            mm = np.linalg.norm(mm[:, 0], ord=order)
            ss = np.linalg.norm(ss, ord=order)
            uns.append(ss/mm)

    uns = np.array(uns)
    uns = uns/min(uns)
    return post_means, bics, uns

def BICs(best_subsets, dataset, u_type='var', take_sqrt=True):
    assert u_type == 'var' or 'cv' in u_type
    XX, yy = dataset
    bics = []
    uncertainties = []
    for efi in best_subsets:
        fit_res = sm.OLS(yy, XX[:, efi]).fit()
        bics.append(fit_res.bic)
        mm = fit_res.params
        ss = fit_res.bse
        if not take_sqrt:
            ss = ss**2
        if u_type == 'var':
            uncertainties.append(ss.sum())
        elif 'cv' in u_type:
            code = u_type.replace('cv', '')
            if len(code) == 0: order = 1
            else: order = int(u_type.replace('cv', ''))
            mm = np.linalg.norm(mm, ord=order)
            ss = np.linalg.norm(ss, ord=order)
            uncertainties.append(ss/mm)

    bics = np.array(bics)
    uncertainties = np.array(uncertainties)
    uncertainties = uncertainties/min(uncertainties)
    return bics, uncertainties

def UBIC(BICs, uncertainties, n_samples, hyp=1, scale=None):
    assert len(BICs) == len(uncertainties)
    if scale is None:
        scale = np.log(n_samples)
    return BICs + hyp*scale*uncertainties

def UBICs(best_subsets, dataset, u_type='var', take_sqrt=True, use_baye=False, ridge_lambda=0, delta=1, n_lams=3, max_lam=11):
    assert u_type == 'var' or 'cv' in u_type
    assert n_lams > 1 and len(dataset) == 2
    print(f"n_lams = {n_lams}") # Use ics[-n_lams]
    delta = float(delta)
    n_samples = dataset[0].shape[0]

    if use_baye:
        print("Using baye_uncertainties")
        _, bics, uncertainties = baye_uncertainties(best_subsets, dataset, u_type, take_sqrt, ridge_lambda)
    else:
        print("Using OLS's uncertainties")
        bics, uncertainties = BICs(best_subsets, dataset, u_type, take_sqrt)
    print(uncertainties)

    lam = 0
    ics = []
    bcs = []
    bc2lam = {}

    uniq = []
    while lam <= max_lam:
        hyp = 10**lam
        ic = UBIC(bics, uncertainties, n_samples, hyp)
        bc = np.argmin(ic)

        if bc > 0:
            ics.append(ic)
            bcs.append(bc)
            if bc not in bc2lam:
                bc2lam[bc] = lam
            print(lam, '--->', bc+1)
        else:
            lam -= delta
            break

        bcl = bcs[-n_lams:]
        if len(bcl) == n_lams:
            uniq = np.unique(bcl)
            if len(uniq) == 1:
                break

        lam += delta

    if len(uniq) == 0:
        uniq, idx = find_modes(bcs)
        uniq = sorted(uniq)
    uniq = min(uniq)
    lam = bc2lam[uniq]
    print(UBIC(bics, uncertainties, n_samples, 10**lam))

    # Checking improvement
    bad_condition = True
    if bics[uniq] < bics[uniq-1]:
        percent_improve = abs((bics[uniq]-bics[uniq-1])/bics[uniq-1])
        # print(percent_improve)
        if percent_improve > 0.08:
            bad_condition = False
    if bad_condition:
        print(f"{uniq} to {uniq+1} complexity may not be worthy with percent_improve = {percent_improve}.")
        print(f"Staying at {uniq} complexity...")
        uniq = uniq-1
        if uniq < 0:
            print("Warning: consider decreasing delta to get more sensible results!")
    print(f"The optimal complexity is currently at the support sizes of {uniq+1}.")

    return ics, uniq, lam

# Algorithm 1 in the UBIC paper
def ubic_auto_select(best_subsets, dataset, scale=None, thres=0.02, percent=None, tau=3, verbose=False):
    if scale is None:
        scale = np.log(len(dataset[-1]))

    post_means, b_bics, b_uns = baye_uncertainties(best_subsets, (dataset[0], dataset[-1]), u_type='cv1', take_sqrt=True)

    b_bics = np.array(b_bics)
    complexities = np.arange(len(b_bics))+1
    d_complexities = complexities[decreasing_values_indices(b_bics)]
    d_bics = b_bics[decreasing_values_indices(b_bics)]

    if percent is not None:
        thres = np.percentile(np.abs(np.diff(d_bics)/(np.diff(d_complexities)*d_bics[:-1])), percent)

    predictions = X_pre@post_means
    lower_bounds = []
    for k, efi in enumerate(best_subsets):
        assert len(efi) == np.count_nonzero(post_means[:, k:k+1])
        com = len(efi)
        # lower_bound = 2*log_like_value(predictions[:, k:k+1], dataset[-1])/np.log(len(dataset[-1]))-com
        lower_bound = 2*log_like_value(predictions[:, k:k+1], dataset[-1])-np.log(len(dataset[-1]))*com
        lower_bounds.append(lower_bound)
    # last_lam = np.log10(max(max(lower_bounds/b_uns)), 0)
    last_lam = np.log10(max(max(lower_bounds/(b_uns*scale)), 0))
    if verbose:
        print("max_lam:", last_lam)

    delta = last_lam/tau
    now_lam = last_lam-delta
    last_ubic = UBIC(b_bics, b_uns, len(y_pre), hyp=10**last_lam)
    last_bc = np.argmin(last_ubic)
    while now_lam > 0:
        now_ubic = UBIC(b_bics, b_uns, len(y_pre), hyp=10**now_lam)
        now_bc = np.argmin(now_ubic)

        diff_com = now_bc-last_bc
        diff_bic = b_bics[now_bc]-b_bics[last_bc]

        imp = np.nan
        if diff_com > 0:
            imp = abs(diff_bic/(b_bics[last_bc]*diff_com))

        if verbose:
            print(min(last_bc, now_bc), '<--->', max(last_bc, now_bc),
                  np.nan_to_num(imp, nan=np.inf))

        if (diff_com > 0 and (diff_bic > 0 or imp < thres)) or \
            (diff_com < 0 and diff_bic > 0 and imp > thres):
            break

        last_lam = now_lam
        now_lam = last_lam-delta
        last_ubic = now_ubic
        last_bc = now_bc

    best_bc = last_bc
    if abs((b_bics[last_bc]-b_bics[last_bc-1])/b_bics[last_bc-1]) < thres:
        best_bc = best_bc - 1
    last_lam = round(last_lam, 10)

    return last_lam, last_ubic, last_bc, best_bc
