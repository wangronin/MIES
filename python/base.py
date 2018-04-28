# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 17:16:39 2018

@author: Hao Wang
@email: wangronin@gmail.com

"""

import pdb
import numpy as np
from tabulate import tabulate

# TODO: test the performance against Pandas dataframe
# TODO: maybe set var_name as a constant attribute
# TODO: register the objective function the class and implement an eval function  

class Solution(np.ndarray):
    """
    Subclassing numpy array to represent set of solutions in the optimization
    Goal to achieve: 
        1) heterogenous data types, like pandas
        2) easy indexing as np.ndarray
        3) extra attributes (e.g., fitness) sliced together with the solution
    """
    # ger rid of self.__dict__ for speed concern
    __slots__ = 'N', 'dim', 'var_name', 'index', 'fitness', 'n_eval', 'verbose'
    def __new__(cls, x, fitness=None, n_eval=0, index=None, var_name=None, verbose=False):
        """
        Parameters
        ----------
            x : array-like,
            fitness : float (array),
                objective value of solutions
            n_eval : int (array),
                number of evaluations per each solution           
            index : int (array),
                indices of solutions
            var_name : str (array),
                column names
            verbose : bool,
                controls if additional information are printed when calling __str__
                and to_dict
        Note
        ----
            If attributes index, fitness, n_eval are modified in a slice of Solution, 
            the corresponding attributes in the original object are also modified. 
            var_name is not affected by this behavior
        """
        obj = np.asarray(x, dtype='object').view(cls)

        if len(obj.shape) > 2:
            raise Exception('More than 2D is not supported')
        
        obj.N = 1 if len(obj.shape) == 1 else obj.shape[0]
        obj.dim = obj.shape[0] if len(obj.shape) == 1 else obj.shape[1]
        
        if not hasattr(fitness, '__iter__'):
            fitness = [fitness] * obj.N
        if not hasattr(n_eval, '__iter__'):
            n_eval = [n_eval] * obj.N
            
        if index is None:
            index = list(range(obj.N))
        assert len(index) == obj.N
        
        if var_name is None:
            var_name = ['x' + str(i) for i in range(obj.dim)]
        assert len(var_name) == obj.dim
        
        index = index[0] if len(index) == 1 else index
        var_name = var_name[0] if len(var_name) == 1 else var_name
        
        # np.ndarray is set for those attributes because slicing it returns references
        # avoid calling self.__setattr__ for attributes fitness, n_eval and index
        super(Solution, obj).__setattr__('fitness', np.asarray(fitness, dtype='float'))
        super(Solution, obj).__setattr__('n_eval', np.asarray(n_eval, dtype='float'))
        super(Solution, obj).__setattr__('index', np.asarray(index, dtype='float'))
        obj.var_name = np.asarray(var_name)
        obj.verbose = verbose
        
        return obj      

    def _check_attr(self):
        # TODO: check for the compatibility of fitness and n_eval
        pass

    def __add__(self, other):
        """
        Concatenate two sets of solutions 
        """
        assert isinstance(other, Solution)
        assert self.dim == other.dim
        _ = [self.tolist()] if len(self.shape) == 1 else self.tolist()
        __ = [other.tolist()] if len(other.shape) == 1 else other.tolist()
        return Solution(_ + __,  self.fitness.tolist() + other.fitness.tolist(),
                        self.n_eval.tolist() + other.n_eval.tolist(), 
                        var_name=self.var_name, verbose=self.verbose)

    def __mul__(self, N):
        """
        repeat a solution N times
        """
        assert isinstance(N, int)
        if self.N > 1:
            raise Exception('Replication is not supported for 2D')
        return Solution([self.tolist()] * N, self.fitness.tolist() * N, 
                         self.n_eval.tolist() * N, var_name=self.var_name, verbose=self.verbose)

    def __rmul__(self, N):
        return self.__mul__(N)
    
    def __setattr__(self, name, value):
        attr = getattr(self, name, None)
        if name in ['fitness', 'n_eval', 'index'] and hasattr(attr, '__iter__'):
            attr[:] = value  # IMPORTANT: copy the value (not reference) to the attribute
        else:
            super(Solution, self).__setattr__(name, value)

    def __setitem__(self, index, value):
        # TODO: maybe add variable type checker here...
        super(Solution, self).__setitem__(index, value)
    
    def __getitem__(self, index):
        _, __ = index, slice(None, None)

        if isinstance(index, tuple):
            index = list(index)
            _ = index[0]
            if len(index) == 2 and isinstance(index[1], int):
                index[1] = slice(index[1], index[1] + 1)
            if len(index) == 2:
                __ = index[1]
        
        subarr = super(Solution, self).__getitem__(index)

        # sub-slicing the attributes
        if isinstance(subarr, Solution):
            # make sure an array is always returned after slicing attributes fitness, n_eval
            _ = slice(_, _ + 1) if isinstance(_, int) else _
            if len(self.shape) == 1:   # only one solution
                subarr.var_name = subarr.var_name[_]
            else:
                super(Solution, subarr).__setattr__('fitness', subarr.fitness[_])
                super(Solution, subarr).__setattr__('n_eval', subarr.n_eval[_])
                super(Solution, subarr).__setattr__('index', subarr.index[_])

                # multiple solutions and the column is indexed
                subarr.var_name = subarr.var_name[__]

            subarr.N = 1 if len(subarr.shape) == 1 else subarr.shape[0]
            subarr.dim = subarr.shape[0] if len(subarr.shape) == 1 else subarr.shape[1]

        return subarr

    def __array_finalize__(self, obj):
        """
        __array_finalize__ is called after new 'Solution' instance is created: from calling
        1) __new__, 2) view casting (ndarray.view()) or 3) slicing (__getitem__) 
        """
        if obj is None: return
        # Needed for array slicing (__getitem__)
        super(Solution, self).__setattr__('fitness', getattr(obj, 'fitness', None))
        super(Solution, self).__setattr__('n_eval', getattr(obj, 'n_eval', None))
        super(Solution, self).__setattr__('index', getattr(obj, 'index', None))
        self.var_name = getattr(obj, 'var_name', None)
        self.verbose = getattr(obj, 'verbose', None)
        self.dim = getattr(obj, 'dim', None)
        self.N = getattr(obj, 'N', None)
    
    def to_dict(self):
        # avoid calling self.__getitem__
        res = {k : super(Solution, self).__getitem__((slice(None, None), i)).tolist() \
            for i, k in enumerate(self.var_name)} 

        if self.verbose:
            res['index'] = self.index.tolist()
            res['fitness'] = self.fitness.tolist()
            res['n_eval'] = self.n_eval.tolist()
        return res
    
    def __str__(self):
        var_name = self.var_name.tolist()
        headers = var_name + ['n_eval', 'fitness'] if self.verbose else var_name
        if len(self.shape) == 1:
            t = [self.tolist() + self.n_eval.tolist() + self.fitness.tolist()] if self.verbose else \
                [self.tolist()]
        else:
            t = np.c_[self, self.n_eval, self.fitness] if self.verbose else self
        
        return tabulate(t, headers=headers, 
                        showindex=self.index.tolist(), tablefmt='fancy_grid').encode('utf-8')
                        
    def to_csv(self, fname, delimiter=',', header=True, index=True, verbose=True):
        if header:
            _header = self.var_name.tolist()
            if index:
                _header = [''] + _header
            if verbose:
                _header += ['n_eval', 'fitness']
            _header = ','.join(_header) + '\n'
        
        data = self.reshape(1, -1) if len(self.shape) == 1 else self
        if index:
            data = np.c_[self.index, data]
        if verbose:
            data = np.c_[data, self.n_eval, self.fitness] 

        out = [','.join(map(str, row)) + '\n' for row in data.tolist()]
        with open(fname, 'w') as f:
            if header:
                f.writelines(_header)
            f.writelines(out)

if __name__ == '__main__':
    # test for 2D solution
    s = Solution(np.c_[np.random.randn(5, 3), ['simida', 'niubia', 'bang', 'GG', 'blyat']], 
                       verbose=True)
    print(s.to_dict())
    print(s)
    print(s[0])
    print(s[0:1])
    print(s[0, 0:3])
    print(s[:, 0])
    print(s[0:2][0, 0:2])
    print(s[0][0:2])
    
    
    s[:, 0] = np.asarray(['wa'] * 5).reshape(-1, 1)
    print(s)

    a = s[0]
    a.fitness = 3
    print(s)

    s[0].fitness = 2
    print(s)
    
    s[0:2].fitness = [1, 2]
    s[0:2].n_eval += 1
    print(s)

    s.to_csv('test.csv', header=False, verbose=False, index=False)

    # test for 1D solution
    ss = Solution(np.random.randn(10))

    print(ss * 2)
    print(ss[0])
    print(ss)

    print(s[0] + s[3:4])