# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 10:33:01 2018

@author: wangronin
"""


import numpy as np

from numpy import exp
from numpy.random import rand, randn, zeros, ceil
from .utils import boundary_handling
from .base import Solution

class self_adaptive_ES(object):
    def __init__(self, dim, obj_func, bounds=None, max_eval=np.inf,
                 minimize=True, mu_=4, lambda_=10, sigma0=None):
        
        self.dim = dim
        self.mu_ = mu_
        self.lambda_ = lambda_
        self.eval_count = 0
        self.iter_count = 0
        self.max_eval = max_eval
        self.plus_selection = False
        
        self.pop = rand(self.mu_, self.dim) * (self.bounds[1, self.id_r] - self.bounds[0, self.id_r]) + \
            self.bounds[0, self.id_r]
        
        sigma0 = 0.05 * (self.bounds[1, self.id_r] - self.bounds[0, self.id_r]) \
            if sigma0 is None else sigma0
        self.sigma = np.tile(sigma0, (self.mu_, 1))
            
        self.fopt = min(self.f_mu)
        self.xopt = self.pop_mu[self.f_mu == self.fopt]
        
        self.pop_lambda = np.zeros((self.lambda_, self.dim))
        self._set_hyperparameter()

        # stop criteria
        self.tolfun = 1e-5
        self.nbin = int(3 + ceil(30. * self.dim / self.lambda_))
        self.histfunval = zeros(self.nbin)

    def optimize(self):
        pass
    
    def recombine(self, id1, id2):
        x1 = self.pop[id1]
        s1 = self.sigma[id1]
        if id1 != id2:
            x2 = self.pop[id2]
            s2 = self.sigma[id2]
            # intermediate recombination for the mutation strengths
            s1 = (np.array(p1[self._id_hyperpar]) + \
                np.array(p2[self._id_hyperpar])) / 2
            # dominant recombination for solution parameters
            mask = randn(self.dim) > 0.5
            x1[mask] = x2[mask]
            
        return p1
        
    def mutate(self, x, sigma):
        if len(self._id_sigma) == 1:
            sigma = sigma * exp(self.tau_r * randn())
        else:
            sigma = sigma * exp(self.tau_r * randn() + self.tau_p_r * randn(self.N_r))
        
        r = randn(self.N_r)
        new = x + sigma * r
        if 11 < 2:
            new_in_bound = boundary_handling(new.copy(), self.bounds[0, self.id_r], \
                self.bounds[1, self.id_r])
        else:
            new_in_bound = new
        
        if 11 < 2:
            sigma_new = np.abs((new_in_bound - x) / r)
        else:
            sigma_new = sigma

        return new_in_bound, sigma_new
        
    def _set_hyperparameter(self):
        # hyperparameters: mutation strength adaptation
        self.tau_r = 1. / np.sqrt(2 * self.N_r)
        self.tau_p_r = 1. / np.sqrt(2 * np.sqrt(self.N_r))

        