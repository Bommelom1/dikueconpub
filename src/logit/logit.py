import numpy as np
import scipy 
from scipy import optimize 
import copy
import matplotlib.pyplot as plt
import time 

def compute_utils(x,z,beta,MAX_RESCALE=True): 
    """
        x: N*Kx matrix of characteristics for individuals 
        z: Kz*J characteristics of the discrete choices 
        beta: Kx*Kz matrix of coefficients 
        MAX_RESCALE: switch for whether to rescale utilities linearly with the max
        across utilities (for numerical stability). 
    """ 
    
    # 1. linear algebra
    U = np.dot(np.dot(x,beta), z) # N*J matrix 
    
    # 2. subtract max for numerical stability 
    #    (probabilities are unaffected by adding the same scalar 
    #     to all utilities)
    if MAX_RESCALE: 
        Umax = np.max(U,axis=1,keepdims=True) # N-vector 
        U -= Umax # elementwise, broadcasts Umax t- N*J
    
    return U

def ccps(x,z,beta): 

    # 1. unpack 
    Kx,Kz = unpack_dimensions(x,z)
    N = x.shape[0]
    if Kz>1: 
        J = z.shape[1]
    else: 
        J = len(z)
    
    # something about subtracting a vector from a matrix fails when N==J
    assert N is not J , 'Broadcasting of denominator in division will fail when N == J.'

    EU =  np.exp(compute_utils(x,z,beta,True)) # N*J
    prob = EU / np.sum(EU,axis=1,keepdims=True)
    
    return prob

def simulate_choices(ccp): 
    """
        returns an N-vector of simulated discrete choices (integers, d=0,...,J-1)

        ccp: N*J matrix of choice probabilities
        (i.e. np.sum(ccp, axis=1) == 1.0 for all rows)
    """

    N,J = ccp.shape

    # pre-allocate
    d = np.zeros((N,),dtype=int)
    
    ccp_cum = np.cumsum(ccp, axis=1) # output is also N*J 
    u = np.random.uniform(0,1,size=(N,))

    # iteratively update
    # alternatively, we could find the highest column where u[i]>ccp_cum[i,:]
    for j in range(J):
        I = u>ccp_cum[:,j] # N-vector of booleans 
        d[I] = j
        
    return d 

def loglikelihood(x,z,d,beta): 
    """
        OUTPUT: 
            N-vector of log-likelihoods (float)

        INPUTS: 
            x: N*Kx, individual characteristics 
            z: Kz*J, choice characteristics 
            d: N vector of choices (int) 
            beta: Kx*Kz, coefficients
    """ 

    N = x.shape[0]
    U = compute_utils(x,z,beta,MAX_RESCALE=True)
    
    logsum = np.log(np.sum(np.exp(U), axis=1))
    U_d = U[range(N), d]
    return U_d - logsum 

def unpack_dimensions(x,z): 
    """
        finds Kx,Kz (ints) based on dimensions of x,z
    """
    if len(x.shape) == 1: 
        N = x.shape
        Kx = 1
    else: 
        N,Kx = x.shape 
    if len(z.shape) == 1: 
        J = z.shape
        Kz = 1
    else: 
        Kz,J = z.shape
    return Kx,Kz

def neg_ll(beta,x,z,d): 
    """
        Negative log likelihood function 
        Accepts beta to be in vector form and rescales to the matrix-form 
    """
    Kx,Kz = unpack_dimensions(x,z)
    
    B = np.reshape(beta, (Kx,Kz))
    return - np.mean(loglikelihood(x,z,d,B))

def monte_carlo(x,z,beta_true,R,method='BFGS',PERTURB_STARTING_VALUE=True): 
    """
        Run R iterations where in each step, we simulate choices (based
        on the choice probabilities calculated from the true parameters)
        and estimate parameters on that data. 

        OUTPUT: 
            betas: K*R matrix of estimated coefficients, where K = Kz*Kx

        INPUTS: 
            method: passed to scipy.optimize.minimize()
            PERTURB_STARTING_VALUE: bool, whether to add a little noise 
            to the starting values so we don't start directly in the true
            parameters. It's just a uniformly distributed scalar that 
            scales all parameters.    
    """

    Kx,Kz = unpack_dimensions(x,z)
    
    # compute true choice probabilities 
    ccp = ccps(x,z,beta_true)
    
    # pre-allocate space for output 
    betas = np.empty((Kx*Kz, R), dtype=float)

    # convenient output 
    nits = np.empty((R,), dtype=int)
    success = np.empty((R,), dtype=bool)
    times = np.empty((R,), dtype=float)

    for r in range(R): 

        t = time.time()

        # 1. resimulate data 
        d = simulate_choices(ccp)

        # 2. estimate 
        if PERTURB_STARTING_VALUE: 
            beta0 = beta_true*np.random.normal(loc=1.0, scale=0.05) 
        else: 
            beta0 = beta_true
        betahat,res = estimate(x,z,d,beta0,method=method)

        # 3. store 
        betas[:,r] = np.reshape(betahat, Kx*Kz)
        nits[r] = res.nit
        success[r] = res.success
        times[r] = time.time() - t
    
    print(f'MC done. Mean nit = {np.mean(nits)}; mean time = {np.mean(times): 5.2f} sec/est. Number of convergence failures = {np.sum(success == False)}')

    return betas 

def simulate_demographics(N,J,Kz,Kx): 
    x = np.ones((N,Kx))
    x[:,1:] = 1.*np.random.normal(0,1,size=(N,Kx-1))

    z = 1.*np.random.normal(0,1,size=(Kz,J))
    
    return x,z

def estimate(x,z,d,beta0,method='BFGS'): 
    """
        Calls scipy.optimize.minimize() to estimate parameters using maximum likelihood. 
    """
    Kx,Kz = unpack_dimensions(x,z)
    res = scipy.optimize.minimize(fun=neg_ll, x0=beta0, args=(x,z,d), method=method)
    betahat = np.reshape(res.x, (Kx,Kz))
    return betahat,res 

def plot_negll_profile(x,z,d,beta,ix=0,iz=0,h=1e-1,num_points=100): 
    """
        INPTS: 
            ix,iz: indexes the parameter we will be plotting the criterion over (beta[ix,iz])
            h: relative step around the truth 
            num_points: number of points on which to evaluate the criterion 
    """
    bb = np.linspace(1.-h, 1.+h, num_points) * beta[0, 0]
    yy = np.empty(bb.shape)
    for i,b_ in enumerate(bb): 
        bet = copy.deepcopy(beta)
        bet[0,0] = b_ 
        yy[i] = neg_ll(bet,x,z,d)

    fig,ax = plt.subplots()
    ax.plot(bb, yy, '-x', label='Log-likelihood')
    ax.plot(beta[0,0], neg_ll(beta,x,z,d), '-or', label='Truth'); 
    ax.legend(); 

def analyze_monte_carlo_estimates(beta_true, beta_mc, ix=0, iz=0): 
    R = beta_mc.shape[1]
    Kx,Kz = beta_true.shape
    i_mc = iz + ix*Kz # beta_mc has shape K*R, where K = Kx*Kz

    print(f'--- Monte carlo results over R={R} simulations ---')
    _b = beta_mc[i_mc,:]
    print(f'beta_00: truth = {beta_true[ix,iz]: 8.4f}, MC mean = {np.mean(_b) : 8.4f} (P10 = {np.percentile(_b, 10): 8.4f}, P90 = {np.percentile(_b, 90): 8.4f})')

    fig,ax = plt.subplots(); 
    ax.hist(beta_mc[i_mc,:], 50);
    ax.axvline(beta_true[ix,iz], linestyle='dashed',color='red'); 
    ax.set_xlabel(f'beta[ix={ix},iz={iz}]'); 
    print(f'Horizontal red line indicates the truth (beta_00 = {beta_true[ix,iz] : 8.4f})')

