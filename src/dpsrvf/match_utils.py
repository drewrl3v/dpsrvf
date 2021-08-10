from dpsrvf import dpmatch
import numpy as np
from scipy.interpolate import interp1d

def group_action_by_gamma(q, gamma):
    '''
    Computes composition of q and gamma and normalizes by gradient
    Inputs:
    -q: An (n,T) matrix representing a Square-Root Velocity Function
        n is the dimension of the space, i.e. R^n
        T is the number of points that make up the curve q
    -gamma: A (T,) dimensional vector representing the warp to apply to q
    '''
    n, T = q.shape
    gamma_t = np.gradient(gamma, 1/(T-1))
    f = interp1d(np.linspace(0, 1, T, True), q, kind = 'linear', fill_value = 'extrapolate')
    q_composed_gamma = f(gamma)

    sqrt_gamma_t = np.tile(np.sqrt(gamma_t), (n,1))
    qn = np.multiply(q_composed_gamma, sqrt_gamma_t)

    return qn


def match(q1, q2):
    '''
    Input: Two Square-Root Velocity Functions, represented as (n,T) matrices
    Output: A warped SRVF q2n matching q1 and the warping function gamma
    '''
    gamma = dpmatch().match(q1,q2)
    gamma = gamma/np.max(gamma)
    q2n = group_action_by_gamma(q2, gamma)

    return q2n, gamma

def fisher_rao(q1,q2):
    '''
    Input: Two Square-Root Velocity Functions, represented as (n,T) matrices
    Output: The Fisher-Rao distance between the two SRVF's
    '''

    q2n, gamma = match(q1, q2)
    gamma_dot = np.gradient(gamma, 1 / (len(gamma) - 1))
    gamma_dot_sqrt = np.sqrt(gamma_dot)
    q2n_tilde = np.zeros_like(q2n)

    n, _ = q2n.shape
    for i in range(n):
        q2n_tilde[i,:] = np.multiply(gamma_dot_sqrt, q2n[i,:])
    
    return induced_norm_L2(q1 - q2n_tilde)

