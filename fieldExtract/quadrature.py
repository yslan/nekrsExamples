import numpy as np
import scipy.special as sps

def zwgl(n):
    lege_nodes, lege_wts, _ = sps.legendre(n).weights.T
    return lege_nodes, lege_wts

def zwgll(n):
    # L_N: N-th order Legendre polyn.
    # GLL: zeros of L'_N and {-1, 1}
    
    N = n-1
    tol = 1e-15
    
    z = np.zeros(n)
    w = np.zeros(n)
    
    w0 = 2.0 / (N*n)
    z[0], z[N] = -1.0, 1.0
    w[0], w[N] = w0, w0
    
    # initial guess = Chebyshev nodes
    cheby = np.array([np.cos( (2*k+1)*np.pi/(2*n) ) for k in range(n)])
    cheby = cheby[::-1]

    for i in range(1,N):
        x = cheby[i]
        dx = np.inf
        it = 0
        while np.abs(dx)>tol and it < 10:
            it += 1
        
            # compute p = L(x) and pp = L'(x), ppp=L''(x)
            pk1, pk0 = 1.0, x    # p_{k}, p_{k-1}
            ppk1, ppk0 = 0.0, 1.0
            pppk1, pppk0 = 0.0, 0.0
            for k in range(1,n-1):
                k1, k2 = k+1, 2*k+1
                p = (k2*x*pk0 - k*pk1) / k1
                pp = (k2*(x*ppk0+pk0) - k*ppk1) / k1
                ppp = (k2*(x*pppk0 + 2*ppk0) - k*pppk1) / k1
        
                pk1,ppk1,pppk1 = pk0,ppk0,pppk0
                pk0,ppk0,pppk0 = p,pp,ppp
        
            # newton
            dx = - pp/ppp
            x = x + dx
        
        z[i] = x
        w[i] = 2.0 / ((n-1)*n * p**2)
    return z, w

def dhat(x):
    # bary centric formulation
    assert x.ndim == 1
    
    n = x.size
    a = np.ones((n,), dtype=np.float64)
    
    for i in range(n):
        for j in range(i):
            a[i] = a[i] * (x[i] - x[j])
        for j in range(i + 1, n):
            a[i] = a[i] * (x[i] - x[j])
    a = 1.0 / a
    D = np.zeros((n, n))

    for i in range(n):
        D[i, :] = a[i] * (x[i] - x)
        D[i, i] = 1.0
    for j in range(n):
        D[:, j] = D[:, j] / a[j]
    D = 1.0 / D

    for i in range(n):
        D[i, i] = 0
        D[i, i] = -np.sum(D[i, :])

    return D


def semhat(N, basis='GLL'):
    basis_map = {
        'GL': zwgl,
        'GLL': zwgll,
    }

    try:
        z, w = basis_map[basis](N + 1)
    except KeyError:
        raise ValueError(f"Unknown basis '{basis}'. Options: {list(basis_map.keys())}")

    D = dhat(z)
    return z, w, D

def gradr2(Dh, u):
    assert len(u.shape)==2, 'u must be 2D arr'
    nh = Dh.shape[0]
    nxyz, E = u.shape
    dim = 2 
    assert nh**dim==nxyz, 'size mismatched'

    ur = np.zeros_like(u)
    us = np.zeros_like(u)
    for e in range(E):
        ue = u[:,e].reshape((nh,nh), order='F')
        ure = Dh @ ue
        use = ue @ Dh.T
        ur[:,e] = ure.reshape(-1, order='F')
        us[:,e] = use.reshape(-1, order='F')
    return ur, us

def gradrs2(Dr, Ds, u):

    E = u.shape[-1]
    nr = Dr.shape[0]
    ns = Ds.shape[0]
    if len(u.shape)==2:
        assert u.shape[0] == nr*ns
        u = u.reshape((nr,ns,E), order='F')
    else:
        nx, ny, E = u.shape
        assert nr == nx, 'size mismatched'
        assert ns == ny, 'size mismatched'

    dim = 2 

    ur = np.zeros((nr, ns, E), order='F')
    us = np.zeros((nr, ns, E), order='F')
    for e in range(E):
        ue = u[:,:,e]
        ur[:,:,e] = Dr @ ue
        us[:,:,e] = ue @ Ds.T
    return ur, us

