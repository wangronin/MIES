# -*- coding: utf-8 -*-
"""
Created on Wed Apr 25 16:12:26 2018

@author: wangronin
"""

import numpy as np
import subprocess

N = int(500)
fopt = np.zeros((1, N))
for i in range(N):
    r = np.random.randint(0, 1000)
    a = subprocess.check_output(['./MIES', '-d', '2', '-t', './sphere', '--mu', '4', '--lambda', '10',
                             '--kappa', '1', '-b', '-5', '5', '-5', '5',
                             '-m', '500', '-r', str(r), '-V', '0'])
    print(a)
    fopt[0, i] = eval(a)
    
np.savetxt('mies.csv', fopt, delimiter=',')