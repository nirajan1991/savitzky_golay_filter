
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 27 17:38:30 2021

@author: Nirajan
The function is based on the function available at
#https://dsp.stackexchange.com/questions/1676/savitzky-golay-smoothing-filter-for-not-equally-spaced-data
#The weights were added to the function follwing similar function for R available at
#https://rdrr.io/github/ranghetti/sen2rts/src/R/w_savgol.R

"""
import numpy as np
#import functools
def non_uniform_savgol(x, y, weight, window=5, polynom=2):
    """
    Applies a Savitzky-Golay filter to y with non-uniform spacing
    as defined in x

    This is based on https://dsp.stackexchange.com/questions/1676/savitzky-golay-smoothing-filter-for-not-equally-spaced-data
    The borders are interpolated like scipy.signal.savgol_filter would do

    Parameters
    ----------
    x : array_like
        List of floats representing the x values of the data
    y : array_like
        List of floats representing the y values. Must have same length
        as x
    weight: array_like 
        List of floats between 0 and 1 representing y values. Must have 
        same length as x and y
    window : int (odd)
        Window length of datapoints. Must be odd and smaller than x
    polynom : int
        The order of polynom used. Must be smaller than the window size

    Returns
    -------
    np.array of float
        The smoothed y values
    """
    if len(x) != len(y) or len(x) != len(weight) or len(y) != len(weight):
        raise ValueError('"x", "y" and "weight" must be of the same size')

    if len(x) < window:
        raise ValueError('The data size must be larger than the window size')

    if type(window) is not int:
        raise TypeError('"window" must be an integer')

    if window % 2 == 0:
        raise ValueError('The "window" must be an odd integer')

    if type(polynom) is not int:
        raise TypeError('"polynom" must be an integer')

    if polynom >= window:
        raise ValueError('"polynom" must be less than "window"')

    half_window = window // 2
    polynom += 1

    # Initialize variables
    A = np.empty((window, polynom))     # Matrix
    tA = np.empty((polynom, window))    # Transposed matrix
    t = np.empty(window)                # Local x variables
    y_smoothed = np.full(len(y), np.nan)

    # Start smoothing
    for i in range(half_window, len(x) - half_window, 1):
        # Center a window of x values on x[i]
        for j in range(0, window, 1):
            t[j] = x[i + j - half_window] - x[i]
        
        # Diagonal matrix of weights
        w = weight[range(i-half_window-1, i+half_window)] # vector of weights
        w = w * len(w) / sum(w) # normalise
        W = np.diag(w) # diagonal matrix
        
        # Create the initial matrix A and its transposed form tA
        for j in range(0, window, 1):
            r = 1.0
            for k in range(0, polynom, 1):
                A[j, k] = r
                tA[k, j] = r
                r *= t[j]

        # Multiply the three matrices
        #AA = functools.reduce(np.matmul,tA,W, A)
        AA = tA @ W @ A

        # Invert the product of the matrices
        tAA = np.linalg.inv(AA)

        # Calculate the pseudoinverse of the design matrix
        #coeffs = functools.reduce(np.matmul,tAA, tA,  W)
        coeffs = tAA @ tA @ W

        # Calculate c0 which is also the y value for y[i]
        y_smoothed[i] = 0
        for j in range(0, window, 1):
            y_smoothed[i] += coeffs[0, j] * y[i + j - half_window]

        # If at the end or beginning, store all coefficients for the polynom
        if i == half_window:
            first_coeffs = np.zeros(polynom)
            for j in range(0, window, 1):
                for k in range(polynom):
                    first_coeffs[k] += coeffs[k, j] * y[j]
        elif i == len(x) - half_window - 1:
            last_coeffs = np.zeros(polynom)
            for j in range(0, window, 1):
                for k in range(polynom):
                    last_coeffs[k] += coeffs[k, j] * y[len(y) - window + j]

    # Interpolate the result at the left border
    for i in range(0, half_window, 1):
        y_smoothed[i] = 0
        x_i = 1
        for j in range(0, polynom, 1):
            y_smoothed[i] += first_coeffs[j] * x_i
            x_i *= x[i] - x[half_window]

    # Interpolate the result at the right border
    for i in range(len(x) - half_window, len(x), 1):
        y_smoothed[i] = 0
        x_i = 1
        for j in range(0, polynom, 1):
            y_smoothed[i] += last_coeffs[j] * x_i
            x_i *= x[i] - x[-half_window - 1]

    return y_smoothed
#%%
'''
# Test the function
y = np.array([1,2,4,6,8,9,4,5,1,6,8,5,2,5])
x = np.array([1,2,4,5,6,7,8,9,11,12,13,14,15, 18])
q = np.random.rand(len(y))

ysm = non_uniform_savgol(x,y,q)
#%%
#Plot the unsmoothed as well as smoothed data
from matplotlib import pyplot as plt
plt.plot(y)
plt.plot(ysm)
'''
#%%