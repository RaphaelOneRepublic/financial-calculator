from typing import Callable
import warnings

import numpy as np


def integral_with_simpsons(func: Callable,
                           begin: float,
                           end: float,
                           tolerance: float = 10e-6):
    """
    calculate the definite integral of function within a certain interval with error less than tolerance
    :param func: function of which definite integral is to be calculated
    :param begin: starting point of integral interval
    :param end: terminating point of integral interval
    :param tolerance: maximum error value
    :return: calculated value of definite integral
    """
    k = 256
    i = _integral_with_simpsons_with_steps(func, begin, end, 128)
    j = _integral_with_simpsons_with_steps(func, begin, end, 256)
    while abs(i - j) > tolerance:
        k *= 2
        i = j
        j = _integral_with_simpsons_with_steps(func, begin, end, k)
    return j


def _integral_with_simpsons_with_steps(func: Callable,
                                       begin: float,
                                       end: float,
                                       steps: int):
    """
    integrate the function using Simpson's method with the interval partitioned into <steps> sub-intervals
    :param func: function of which definite integral is to be calculated
    :param begin: starting point of integral interval
    :param end: terminating point of integral interval
    :param steps: number of sub-intervals to divide into
    :return:
    """
    x = np.linspace(begin, end, steps + 1)
    y = func(x)
    s = np.sum(y[0:-1:2] + y[2::2] + 4 * y[1::2])
    return s * (end - begin) / (3 * steps)


def newtons(func: Callable,
            partial: Callable,
            initial: float,
            tolerance: float = 10e-8,
            error: float = 10e-6):
    """
    use the Newton's method to find a zero near an initial guess of a certain differentiable function
    :param func: function
    :param partial: derivative with respect to the variable in interest
    :param tolerance: maximum error in the underlying variable
    :param error: maximum error in the function value
    :param initial: initial guess
    :return: a zero
    """
    while True:
        temp = initial - func(initial) / partial(initial)
        if abs(temp - initial) < error and abs(func(initial) - func(temp)) < tolerance:
            return temp
        initial = temp


def bisect(func: Callable, begin, end, tolerance=10e-8, error=10e-6):
    """
    use bisection method to find a zero of the function
    :param func: function
    :param begin: starting point of an interval
    :param end:: terminating point of an interval
    :param tolerance: maximum error in the underlying variable
    :param error: maximum error in the function value
    :return:
    """
    warnings.warn("the bisection method converges slowly")
    if func(begin) == 0: return begin
    if func(end) == 0: return end
    if func(begin) > 0 and func(end) > 0 or func(begin) < 0 and func(end) < 0:
        raise ValueError("same sigh at both interval endpoints")
    if func(begin) < 0:
        def increase(x):
            return func(x)
    else:
        def increase(x):
            return -func(x)
    mid = 0
    while abs(func(end) - func(begin)) > tolerance or abs(end - begin) > error:
        mid = (begin + end) / 2
        mv = increase(mid)
        if mv == 0:
            return mid
        elif mv < 0:
            begin = mid
        elif mv > 0:
            end = mid
    return mid


def secant(func: Callable, primary: float, secondary: float, tolerance=10e-8, error=10e-6):
    """
    use secant method to find a zero of the function
    :param func: function
    :param primary: initial guess
    :param secondary: initial support
    :param tolerance: maximum error in the underlying variable
    :param error: maximum error in the function value
    :return:
    """
    while abs(func(primary) - func(secondary)) > tolerance or abs(primary - secondary) > error:
        temp = primary + func(primary) * (primary - secondary) / (func(secondary) - func(primary))
        secondary = primary
        primary = temp
    return primary
