from typing import Callable, Any

import numpy as np


def root(f: Callable[[float], float], x0: float, df: Callable[[float], float] = None,
         epsilon: float = 10e-9, delta: float = 10e-6, stop: int = 10e3, xn1: float = None):
    """
    find a root of the function, using the newton's method if derivative is available, or the secant method if not.

    :param f: univariate function
    :param df: derivative of f(x)
    :param x0: initial guess for the solution of f(x) = 0
    :param epsilon: largest permissible value of abs(f(x)) when solution is found
    :param delta: largest permissible distance between two consecutive approximations when solution is found
    :param stop: maximum iterations
    :param xn1: one of the first initial guesses needed to initialize secant method

    :return: a root
    """
    if df is None:
        if xn1 is None:
            xn1 = x0 + 0.1
        new, old = x0, xn1
        while (abs(f(new) > epsilon) or abs(new - old) > delta) and stop > 0:
            oldest = old
            old = new
            new = old - f(old) * (old - oldest) / (f(old) - f(oldest))
            stop -= 1
    else:
        new, old = x0, x0 + 1
        while (abs(f(new) > epsilon) or abs(new - old) > delta) and stop > 0:
            old = new
            new = old - f(old) / df(old)
            stop -= 1
    return new


def integrate(f: Callable[[float], float], a: float, b: float, epsilon: float = 10e-8, stop: int = 10e3):
    """
    calculate the definite integral of a function
    :param f: function
    :param a: left endpoint
    :param b: right endpoint
    :param epsilon: maximum tolerance

    :return: approximated definite integral with error less than epsilon
    """
    old = simpsons(f, a, b, 4)
    new = simpsons(f, a, b, 8)
    n = 8
    while abs(old - new) > epsilon:
        old = new
        n *= 2
        new = simpsons(f, a, b, n)
    return new


def simpsons(f: Callable[[Any], Any], a: float, b: float, n: int):
    """
    approximate the definite integral of a function using simpson's method
    :param f: function
    :param a: left endpoint
    :param b: right endpoint
    :param n: partitions

    :return: approximated definite integral
    """
    x = np.linspace(a, b, n + 1)
    y = f(x)
    s = np.sum(y[0:-1:2] + y[2::2] + 4 * y[1::2])
    return s * (b - a) / (3 * n)
