from typing import Sequence

import numpy as np


def pv(y: float, C: Sequence[float], due=False):
    """
    Compute the present value of a set of cash flows

    :param y: yield
    :param C: cash received
    :param due: True if the first flow occurs at the present time, or the next period otherwise
    :return: present value of cash flow
    """
    disc = np.power(1 + y, -np.arange(len(C)) - (0 if due else 1))
    return np.sum(disc * C)


def fv(y: float, C: Sequence[float]):
    """
    Compute the future value of a set of cash flows
    :param y:
    :param C:
    :return:
    """
    disc = np.power(1 + y, -np.arange(0, len(C))[::-1])
    return np.sum(disc * C)


def compound2(r: float, m: int, n: int = 0):
    """
    Convert an interest rate compounded m times per annum
    into the equivalent rate compounded n times per annum

    :param r: interest rate
    :param m: compound periods per annum of input, 0 for continuous compounding
    :param n: compound periods per annum of output, 0 for continuous compounding
    :return:
    """
    if m == n:
        return r
    c = r if m == 0 else m * np.log(1 + r / m)
    if n == 0:
        return c
    return n * (np.exp(c / n) - 1)


def amortize(PV: float, n: int, r: float, m: int = 12):
    """
    Compute the level payment of fully amortized traditional mortgages.
    Note that r is bond-equivalent if m = 2, or mortgage-equivalent if m = 12

    <! WARNING: outdated comments>
    Newton's method is used to find the solution.
    However, it is guaranteed to find the solution in less than one iteration
    due to the constant derivative of present value with respect to level payment.
    <! outdated comments>

    The current computation uses the fact that
    the first derivative of present value with respect to level payment is constant

    :param PV: present value of original balance
    :param n: maturity in years
    :param r: annual interest rate, default is 12
    :param m: compound period
    :return: level payment of fully amortized mortgage
    """
    y = 1 + r / m
    return PV * r / m / (1 - y ** (-m * n))


if __name__ == '__main__':
    print(amortize(250000, 15, 0.08))
