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

        # unpack the search space
        self._space = search_space
        self.var_names = self._space.var_name.tolist()
        self.param_type = self._space.var_type

        assert hasattr(self.obj_func, '__call__')

    def _better(self, f1, f2):
        return f1 < f2 if self.minimize else f1 > f2
    
    def recombine(self, id1, id2):
        p1 = self.pop[id1].copy()  # IMPORTANT: make a copy
        if id1 != id2:
            p2 = self.pop[id2]
            # intermediate recombination for the mutation strengths
            p1[self._id_hyperpar] = (np.array(p1[self._id_hyperpar]) + \
                np.array(p2[self._id_hyperpar])) / 2

            # dominant recombination for solution parameters
            _, = np.nonzero(randn(self.dim) > 0.5)
            p1[_] = p2[_]
        return p1

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
        pop = self.pop + self.offspring if self.plus_selection else self.offspring
        fitness = np.r_[self.fitness, self.f_offspring] if self.plus_selection else self.f_offspring
        
        fitness_rank = argsort(fitness)
        if not self.minimize:
            fitness_rank = fitness_rank[::-1]
        
        sel_id = fitness_rank[:self.mu_]
        self.pop = pop[sel_id]
        self.fitness = fitness[sel_id]

    def optimize(self):
        hist_sigma = []
        while not self.stop():
            # hist_sigma += np.mean(np.asarray(self.pop[:, self._id_sigma], dtype='float'), axis=0).tolist()
            hist_sigma += [np.mean(self.fitness)]
            # TODO: vectorize this part
            for i in range(self.lambda_):
                p1, p2 = randint(0, self.mu_), randint(0, self.mu_)
                individual = self.recombine(p1, p2)
                self.offspring[i] = self.mutate(individual)
            
            self.f_offspring[:] = self.evaluate(self.offspring)
            self.select()

            curr_best = self.pop[0]
            xopt_, fopt_ = curr_best[self._id_var], self.fitness[0]

            if self._better(fopt_, self.fopt):
                self.xopt, self.fopt = xopt_, fopt_
            self.iter_count += 1

            if self.verbose:
                print('iteration {}, fopt: {}'.format(self.iter_count + 1, self.fopt))
                print(self.pop)

        self.stop_dict['funcalls'] = self.eval_count
        return self.xopt, self.fopt, self.stop_dict, hist_sigma
    
    def stop(self):
        if self.eval_count > self.max_eval:
            self.stop_dict['max_eval'] = True

        if self.eval_count != 0 and self.iter_count != 0:
            fitness = self.f_offspring
        
        return any(self.stop_dict.values())
    


class self_adaptive_es(evolution_strategy):
    def __init__(self, search_space, obj_func, x0=None, ftarget=None, max_eval=np.inf,
                 minimize=True, mu_=4, lambda_=10, sigma0=None, eta0=None, P0=None,
                 verbose=False):
        super(self_adaptive_es, self).__init__(search_space, obj_func, x0, ftarget, max_eval,
                 minimize, mu_, lambda_, sigma0, eta0, P0,
                 verbose)
        # index of each type of variables in the dataframe
        self.id_r = self._space.id_C       # index of continuous variable

        # the number of variables per each type
        self.dim = self.N_r = len(self.id_r)
        self._len = self.dim + self.N_r
        
        # unpack interval bounds
        self.bounds_r = np.asarray([self._space.bounds[_] for _ in self.id_r])
        # self._check_bounds(self.bounds_r)
        
        # step default step-sizes/mutation strength
        par_name = []
        if sigma0 is None and self.N_r:
            sigma0 = 0.05 * (self.bounds_r[:, 1] - self.bounds_r[:, 0])
            par_name += ['sigma' + str(_) for _ in range(self.N_r)]

        # column indices: used for slicing
        self._id_var = np.arange(self.dim)                    
        self._id_sigma = np.arange(self.N_r) + len(self._id_var)
        self._id_hyperpar = np.arange(self.dim, self._len)
        
        # initialize the populations
        if x0 is not None:                         # given x0
            self.pop = Solution(np.tile(np.r_[x0, sigma0], (self.mu_, 1)),
                                var_name=self.var_names + par_name, verbose=self.verbose)
            fitness0 = self.evaluate(self.pop[0])
            self.fitness = np.repeat(fitness0, self.mu_)
            self.xopt = x0
            self.fopt = sum(fitness0)
        else:                                     # uniform sampling
            x = np.asarray(self._space.sampling(self.mu_), dtype='float')  
            x = np.c_[x, np.tile(sigma0, (self.mu_, 1))]
            
            self.pop = Solution(x, var_name=self.var_names + par_name, verbose=self.verbose)
            self.fitness = self.evaluate(self.pop)
            self.fopt = min(self.fitness) if self.minimize else max(self.fitness)
            self.xopt = self.pop[self.fopt == self.fitness, self._id_var]
            
        self.offspring = self.pop[0] * self.lambda_
        self.f_offspring = np.asarray([self.fitness[0]] * self.lambda_)
        self._set_hyperparameter()

        # stopping criteria
        self.tolfun = 1e-5
        self.nbin = int(3 + ceil(30. * self.dim / self.lambda_))
        self.histfunval = zeros(self.nbin)  

    def _set_hyperparameter(self):
        self.tau_r = 1 / np.sqrt(2 * self.N_r)
        self.tau_p_r = 1 / np.sqrt(2 * np.sqrt(self.N_r))

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


class integer_es(evolution_strategy):
    def mutate(self, individual):
        eta = np.asarray(individual[self._id_eta].tolist(), dtype='float')
        x = np.asarray(individual[self.id_i], dtype='int')
        if len(self._id_eta) == 1:
            eta = eta * exp(self.tau_i * randn())
        else:
            eta = eta * exp(self.tau_i * randn() + self.tau_p_i * randn(self.N_i))
        eta[eta > 1] = 1

        p = 1 - (eta / self.N_i) / (1 + np.sqrt(1 + (eta / self.N_i) ** 2.))
        g1, u1 = zip(*[rgeom(_) for _ in p])
        g2, u2 = zip(*[rgeom(_) for _ in p])
        g1, g2 = map(lambda x: np.asarray(x, dtype='int'), (g1, g2))
        u1, u2 = map(lambda x: np.asarray(x, dtype='float'), (u1, u2))
        x_ = x + g1 - g2

        # Interval Bounds Treatment
        x__ = boundary_handling(x_, self.bounds_i[:, 0], self.bounds_i[:, 1])

        if 11 < 2:
            idx, = np.nonzero(x_ != x__)
            p_ = 1. - (u1 / u2) ** (1. / abs(x__ - x))
            if any(np.isinf(p_)):
                pdb.set_trace()
            eta[idx] = 2. * (1 - p_[idx]) / (p_[idx] * (2. - p_[idx]))

        individual[self.id_i] = x__
        individual[self._id_eta] = eta