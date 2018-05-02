#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 16:28:39 2018

@author: wangronin
"""

import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt

def pmd_geometric(k, p):
    return p / (2. - p) * (1. - p) ** abs(k)

def pdf_discrete_Gaussian(x, scale=1):
    _ = np.asarray(map(round, x))
    return norm.cdf(_ + 0.5, scale=scale) - norm.cdf(_ - 0.5, scale=scale)

def entropy_geometric(p):
    a = 1. - p
    res = a * np.log(a) / (p * (2. - p)) + (p / (2. - p)) * np.log(p / (2. - p)) * (1. / p + 1)
    return -res * 2.

def entropy_discrete_Gaussian(scale=1):
    k = np.arange(1, 201)
    p = pdf_discrete_Gaussian(k)
    res = 0
    for pp in p:
        if not np.isinf(np.log(pp)):
            res += np.sum(pp * np.log(pp)) * 2.
    res += pdf_discrete_Gaussian([0]) * np.log(pdf_discrete_Gaussian([0]))
    return -res

def rgeom(p, size=1):
    u = np.random.rand(size)
    g = np.floor(np.log(u) / np.log(1. - p))
    return g

def MLE(sample):
    Z = sum(np.abs(sample))
    n = len(sample)
    return (n + Z - np.sqrt(n ** 2. + Z ** 2.)) / Z

sample = rgeom(0.3, int(1e2)) - rgeom(0.3, int(1e2))
print(MLE(sample))

k = np.arange(-10, 11)
x = np.arange(-10, 11, 0.01)

p = np.sqrt(3) - 1
scale = 1

y1 = pmd_geometric(k, p)
y2 = pdf_discrete_Gaussian(x, scale)
y3 = norm.pdf(x, scale=scale)
            
with plt.style.context('ggplot'):
    fig0, ax = plt.subplots(1, 1, dpi=100)
    ax.hold(True)
    
    ax.plot(k, y1, 'r+', ms=10)
    ax.plot(x, y2, 'k-', alpha=0.7)
    ax.plot(x, y3, 'k--', alpha=0.7)
    
    ax.set_xticks(k)
    ax.set_xlabel('$X$')
    
    plt.legend(['differential geometric', 'discrete Gaussian', 'Gaussian'])
    plt.show()