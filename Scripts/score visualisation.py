# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 18:09:23 2020

@author: yann-stanislas.barral@roche.com

"""

import numpy as np
import matplotlib.pyplot as plt

def score(x):
    f = np.power(abs(x), 1/3)
    return f

def weight(x):
    f = np.array(x)
    for i in range(0, len(x)):
        if x[i] < 0 :
            f[i] = 2.5 * np.exp((x[i]) / 5) + 0.5
        elif x[i] < 20:
            f[i] = 3
        elif x[i] > 20:
            f[i] = 2.5 * np.exp((20 - x[i]) / 20) + 0.5
    return f

x = np.linspace(-100, 100, 201)
score = score(x)

plt.figure(figsize = (15, 6))
plt.subplot(1, 2, 1)
plt.plot(x, score, LineWidth = 5)
plt.xlabel('$|APD_{90, exp} - APD_{90, sim}|$ (in ms)', Fontsize = 20)
plt.ylabel('Score', Fontsize = 20)


x = np.linspace(-50, 100, 201)
weight_to_plot = weight(x)

plt.subplot(1, 2, 2)
plt.plot(x, weight_to_plot, LineWidth = 5)
plt.xlabel('$|APD_{90, drug} - APD_{90, baseline}|$ (in ms)', Fontsize = 20)
plt.ylabel('Weight', Fontsize = 20)

