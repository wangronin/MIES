# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 11:10:18 2017

@author: wangronin
"""
from __future__ import print_function
import pdb

import numpy as np
from numpy import exp, nonzero, argsort, ceil, zeros, mod
from numpy.random import randint, rand, randn, geometric

from utils import boundary_handling
from base import Solution
from SearchSpace import ContinuousSpace, OrdinalSpace, NominalSpace

class evolution_strategy(object):
    
    def __init__(self, search_space, obj_func, x0=None, ftarget=None, max_eval=np.inf,
                 minimize=True, mu_=4, lambda_=10, sigma0=None, eta0=None, P0=None,
                 verbose=False):
                     
        self.verbose = verbose
        self.mu_ = mu_
        self.lambda_ = lambda_
        self.eval_count = 0
        self.iter_count = 0
        self.max_eval = max_eval
        self.plus_selection = False
        
    def evaluate(self, pop):
        if not hasattr(pop[0], '__iter__'):
            pop = [pop]
        N = len(pop)
        f = np.zeros(N)
        for i, individual in enumerate(pop):
            var = individual[self._id_var]
            f[i] = np.sum(self.obj_func(var)) # in case a 1-length array is returned
            individual.fitness = f[i]
            self.eval_count += 1
        return f
    
    def select(self):
        pop = self.pop_mu + self.pop_lambda if self.plus_selection else self.pop_lambda
        fitness = np.r_[self.f_mu, self.f_lambda] if self.plus_selection else self.f_lambda
        
        fitness_rank = argsort(fitness)
        if not self.minimize:
            fitness_rank = fitness_rank[::-1]
        
        sel_id = fitness_rank[:self.mu_]
        self.pop_mu = pop[sel_id]
        self.f_mu = fitness[sel_id]
        

class self_adaptive_es(evolution_strategy):

    def mutate(self, individual):
        sigma = np.array(individual[self._id_sigma])
        if len(self._id_sigma) == 1:
            sigma = sigma * exp(self.tau_r * randn())
        else:
            sigma = sigma * exp(self.tau_r * randn() + self.tau_p_r * randn(self.N_r))
        
        # Gaussian mutation
        R = randn(self.N_r)
        x = np.array(individual[self.id_r])
        x_ = x + sigma * R
        
        # Interval Bounds Treatment
        x_ = boundary_handling(x_, self.bounds_r[:, 0], self.bounds_r[:, 1])
        
        # Repair the step-size if x_ is out of bounds
        individual[self._id_sigma] = np.abs((x_ - x) / R)
        individual[self.id_r] = x_
        
        return individual


# TODO: improve efficiency, e.g. compile it with cython
# TODO: try to use advanced python default parameter 
class mies(object):
    def __init__(self, search_space, obj_func, x0=None, ftarget=None, max_eval=np.inf,
                 minimize=True, mu_=4, lambda_=10, sigma0=None, eta0=None, P0=None,
                 verbose=False):

        self.mu_ = mu_
        self.lambda_ = lambda_
        self.eval_count = 0
        self.iter_count = 0
        self.minimize = minimize
        self.obj_func = obj_func
        self.stop_dict = {}
        self.verbose = verbose
        self.max_eval = max_eval
        self.ftarget = ftarget
        self.plus_selection = False
        
        self._space = search_space
        self.var_names = self._space.var_name.tolist()
        self.param_type = self._space.var_type

        # index of each type of variables in the dataframe
        self.id_r = self._space.id_C       # index of continuous variable
        self.id_i = self._space.id_O       # index of integer variable
        self.id_d = self._space.id_N       # index of categorical variable

        # the number of variables per each type
        self.N_r = len(self.id_r)
        self.N_i = len(self.id_i)
        self.N_d = len(self.id_d)
        self.dim = self.N_r + self.N_i + self.N_d

        # by default, we use individual step sizes for continuous and integer variables
        # and global strength for the nominal variables
        self.N_p = min(self.N_d, int(1))
        self._len = self.dim + self.N_r + self.N_i + self.N_p  # total length of the solution vector
        
        # unpack interval bounds
        self.bounds_r = np.asarray([self._space.bounds[_] for _ in self.id_r])
        self.bounds_i = np.asarray([self._space.bounds[_] for _ in self.id_i])
        self.bounds_d = np.asarray([self._space.bounds[_] for _ in self.id_d])   # actually levels...
        self._check_bounds(self.bounds_r)
        self._check_bounds(self.bounds_i)
        
        # step default step-sizes/mutation strength
        par_name = []
        if sigma0 is None and self.N_r:
            sigma0 = 0.05 * (self.bounds_r[:, 1] - self.bounds_r[:, 0])
            par_name += ['sigma' + str(_) for _ in range(self.N_r)]
        if eta0 is None and self.N_i:
            eta0 = 0.05 * (self.bounds_i[:, 1] - self.bounds_i[:, 0]) 
            par_name += ['eta' + str(_) for _ in range(self.N_i)]
        if P0 is None and self.N_d:
            P0 = 1. / self.N_d
            par_name += ['P' + str(_) for _ in range(self.N_d)]

        # column indices: used for slicing
        self._id_var = np.arange(self.dim)                    
        self._id_sigma = np.arange(self.N_r) + len(self._id_var)
        self._id_eta = np.arange(self.N_i) + len(self._id_var) + len(self._id_sigma) 
        self._id_p = np.arange(self.N_p) + len(self._id_var) + len(self._id_sigma) + len(self._id_eta)
        self._id_hyperpar = np.arange(self.dim, self._len)
        
         # initialize the populations
        if x0 is not None:                         # given x0
            self.pop_mu = Solution(np.tile(np.r_[x0, sigma0, eta0, [P0] * self.N_p], (self.mu_, 1)),
                                   var_name=self.var_names + par_name, verbose=self.verbose)
            fitness0 = self.evaluate(self.pop_mu[0])
            self.f_mu = np.repeat(fitness0, self.mu_)
            self.xopt = x0
            self.fopt = sum(fitness0)
        else:                                     # uniform sampling
            x = np.asarray(self._space.sampling(self.mu_), dtype='object')  
            
            if sigma0 is not None:
                x = np.c_[x, np.tile(sigma0, (self.mu_, 1))]
            if eta0 is not None:
                x = np.c_[x, np.tile(eta0, (self.mu_, 1))]
            if P0 is not None:
                x = np.c_[x, np.tile([P0] * self.N_p, (self.mu_, 1))]
            
            self.pop_mu = Solution(x, var_name=self.var_names + par_name, verbose=self.verbose)
            self.f_mu = self.evaluate(self.pop_mu)
            self.fopt = min(self.f_mu) if self.minimize else max(self.f_mu)
            self.xopt = self.pop_mu[self.fopt == self.f_mu, self._id_var]
            
        self.pop_lambda = self.pop_mu[0] * self.lambda_
        self.f_lambda = [self.f_mu[0]] * self.lambda_
        self._set_hyperparameter()

        # stopping criteria
        self.tolfun = 1e-5
        self.nbin = int(3 + ceil(30. * self.dim / self.lambda_))
        self.histfunval = zeros(self.nbin)
        
    def _check_bounds(self, bounds):
        if len(bounds) == 0:
            return
        if any(bounds[:, 0] >= bounds[:, 1]):
            raise ValueError('lower bounds must be smaller than upper bounds')

    def _set_hyperparameter(self):
        # hyperparameters: mutation strength adaptation
        if self.N_r:
            self.tau_r = 1 / np.sqrt(2 * self.N_r)
            self.tau_p_r = 1 / np.sqrt(2 * np.sqrt(self.N_r))

        if self.N_i:
            self.tau_i = 1 / np.sqrt(2 * self.N_i)
            self.tau_p_i = 1 / np.sqrt(2 * np.sqrt(self.N_i))

        if self.N_d:
            self.tau_d = 1 / np.sqrt(2 * self.N_d)
            self.tau_p_d = 1 / np.sqrt(2 * np.sqrt(self.N_d))

    def recombine(self, id1, id2):
        p1 = self.pop_mu[id1].copy()  # IMPORTANT: make a copy
        if id1 != id2:
            p2 = self.pop_mu[id2]
            # intermediate recombination for the mutation strengths
            p1[self._id_hyperpar] = (np.array(p1[self._id_hyperpar]) + \
                np.array(p2[self._id_hyperpar])) / 2

            # dominant recombination for solution parameters
            mask = np.nonzero(randn(self.dim) > 0.5)[0]
            p1[mask] = p2[mask]
        return p1

    def select(self):
        pop = self.pop_mu + self.pop_lambda if self.plus_selection else self.pop_lambda
        fitness = np.r_[self.f_mu, self.f_lambda] if self.plus_selection else self.f_lambda
        
        fitness_rank = argsort(fitness)
        if not self.minimize:
            fitness_rank = fitness_rank[::-1]
        
        sel_id = fitness_rank[:self.mu_]
        self.pop_mu = pop[sel_id]
        self.f_mu = fitness[sel_id]

    def evaluate(self, pop):
        if len(pop.shape) == 1:
            f = np.asanyarray(self.obj_func(pop[self._id_var]))
        else:
            f = np.asarray(map(self.obj_func, pop[:, self._id_var]))
        self.eval_count += pop.N
        pop.fitness = f
        return f

    def mutate(self, individual):
        if self.N_r:
            self._mutate_r(individual)
        if self.N_i:
            self._mutate_i(individual)
        if self.N_d:
            self._mutate_d(individual)
        return individual

    def _mutate_r(self, individual):
        sigma = np.asarray(individual[self._id_sigma], dtype='float')
        # mutate step-sizes
        if len(self._id_sigma) == 1:
            sigma = sigma * exp(self.tau_r * randn())
        else:
            sigma = sigma * exp(self.tau_r * randn() + self.tau_p_r * randn(self.N_r))
        
        # Gaussian mutation
        R = randn(self.N_r)
        x = np.asarray(individual[self.id_r], dtype='float')
        x_ = x + sigma * R
        
        # Interval Bounds Treatment
        x_ = boundary_handling(x_, self.bounds_r[:, 0], self.bounds_r[:, 1])
        
        # the constraint handling method will (by chance) turn really bad cadidates
        # (the one with huge sigmas) to good ones and hence making the step size explode
        # Repair the step-size if x_ is out of bounds
        if 1 < 2:
            individual[self._id_sigma] = np.abs((x_ - x) / R)
        else:
            individual[self._id_sigma] = sigma
        individual[self.id_r] = x_
        
    def _mutate_i(self, individual):
        eta = np.asarray(individual[self._id_eta], dtype='float')
        x = np.asarray(individual[self.id_i], dtype='int')
        if len(self._id_eta) == 1:
            eta = max(1, eta * exp(self.tau_i * randn()))
            p = 1 - (eta / self.N_i) / (1 + np.sqrt(1 + (eta / self.N_i) ** 2))
            x_ = x + geometric(p, self.N_i) - geometric(p, self.N_i)
        else:
            eta = eta * exp(self.tau_i * randn() + self.tau_p_i * randn(self.N_i))
            eta[eta > 1] = 1
            p = 1 - (eta / self.N_i) / (1 + np.sqrt(1 + (eta / self.N_i) ** 2))
            x_ = x + np.array([geometric(p_) - geometric(p_) for p_ in p])

        # Interval Bounds Treatment
        x_ = boundary_handling(x_, self.bounds_i[:, 0], self.bounds_i[:, 1])
        individual[self.id_i] = x_
        individual[self._id_eta] = eta

    def _mutate_d(self, individual):
        P = np.asarray(individual[self._id_p], dtype='float')
        P = 1 / (1 + (1 - P) / P * exp(-self.tau_d * randn()))
        individual[self._id_p] = boundary_handling(P, 1 / (3. * self.N_d), 0.5)

        idx = np.nonzero(rand(self.N_d) < P)[0]
        for i in idx:
            level = self.bounds_d[i]
            individual[self.id_d[i]] = level[randint(0, len(level))]

    def stop(self):
        if self.eval_count > self.max_eval:
            self.stop_dict['max_eval'] = True

        if self.eval_count != 0 and self.iter_count != 0:
            fitness = self.f_lambda
            # sigma = np.atleast_2d([__[self._id_sigma] for __ in self.pop_mu]) 
            # sigma_mean = np.mean(sigma, axis=0)
            
            # tolerance on fitness in history
            self.histfunval[int(mod(self.eval_count / self.lambda_ - 1, self.nbin))] = fitness[0]
            if mod(self.eval_count / self.lambda_, self.nbin) == 0 and \
                (max(self.histfunval) - min(self.histfunval)) < self.tolfun:
                    self.stop_dict['tolfun'] = True
            
            # flat fitness within the population
            if fitness[0] == fitness[int(min(ceil(.1 + self.lambda_ / 4.), self.mu_ - 1))]:
                self.stop_dict['flatfitness'] = True
            
            # TODO: implement more stop criteria
            # if any(sigma_mean < 1e-10) or any(sigma_mean > 1e10):
            #     self.stop_dict['sigma'] = True

            # if cond(self.C) > 1e14:
            #     if is_stop_on_warning:
            #         self.stop_dict['conditioncov'] = True
            #     else:
            #         self.flg_warning = True

            # # TolUPX
            # if any(sigma*sqrt(diagC)) > self.tolupx:
            #     if is_stop_on_warning:
            #         self.stop_dict['TolUPX'] = True
            #     else:
            #         self.flg_warning = True
            
        return any(self.stop_dict.values())

    def _better(self, f1, f2):
        if self.minimize:
            return f1 < f2
        else:
            return f1 > f2

    def optimize(self):
        while not self.stop():
            for i in range(self.lambda_):
                p1, p2 = randint(0, self.mu_), randint(0, self.mu_)
                individual = self.recombine(p1, p2)
                self.pop_lambda[i] = self.mutate(individual)
            
            self.f_lambda = self.evaluate(self.pop_lambda)
            self.select()

            curr_best = self.pop_mu[0]
            xopt_, fopt_ = curr_best[self._id_var], self.f_mu[0]
            self.iter_count += 1

            if self._better(fopt_, self.fopt):
                self.xopt, self.fopt = xopt_, fopt_

            if self.verbose:
                print('iteration ', self.iter_count + 1)
                print(self.pop_mu)

        self.stop_dict['funcalls'] = self.eval_count
        return self.xopt, self.fopt, self.stop_dict


if __name__ == '__main__':

#    np.random.seed(1)
    
    if 11 < 2:
        def fitness(x):
            x_r, x_i, x_d = np.array(x[:2]), x[2], x[3]
            if x_d == 'OK':
                tmp = 0
            else:
                tmp = 1
            return np.sum(x_r ** 2) + abs(x_i - 10) / 123. + tmp * 2
    
        x0 = [2, 1, 80, 'B']
        bounds = [[-5, -5, -100], [5, 5, 100]]
        levels = [['OK', 'A', 'B', 'C', 'D', 'E', 'F', 'G']]
    
        opt = mies(x0, fitness, bounds, levels, ['C', 'C', 'O', 'N'], 5e5, verbose=True)
    
    else:
        def fitness(x):
            x = np.asarray(x, dtype='float')
            return np.sum(x ** 2.)
        
        dim = 2
        space = ContinuousSpace([-5, 5]) * dim
        if 11 < 2:
            # test for continous maximization problem
            opt = mies(space, fitness, max_eval=500, minimize=True, verbose=True)
            xopt, fopt, stop_dict = opt.optimize()
            print(stop_dict)
            
        else:
            N = int(500)
            fopt = np.zeros((1, N))
            for i in range(N):
                opt = mies(space, fitness, max_eval=500, verbose=False)
                xopt, fopt[0, i], stop_dict = opt.optimize()
            
            np.savetxt('mies.csv', fopt, delimiter=',')
