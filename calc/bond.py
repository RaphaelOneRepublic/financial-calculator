import logging
from typing import Sequence

import numpy as np

from calc.optimize import root


class Bond(object):
    """
    Represents a coupon paying bond.
    Upon creation, the time to maturity, coupon periods per year, coupon rate must be provided.
    If yield to maturity is provided, bond value would be ignored.
    If yield to maturity is provided, bond value would be used to compute the implied yield to maturity.
    Face value is assumed to be 100 if not provided.
    """

    def __init__(self, T: float, R: float, m: int = 2, y: float = None, F: float = 100, B: float = None):
        """
        construct a coupon paying bond

        :param T: time to maturity in years
        :param m: coupon payments per year
        :param R: quoted annual coupon rate
        :param y: (implied) yield to maturity
        :param (optional) F: face value
        :param B: traded bond price
        """
        self._T = T
        self._m = m
        self._R = R
        self._F = F

        if y is not None:
            self._y = y
            self.__refresh_value_cache__()
        elif B is not None:
            self.B = B
        else:
            raise ValueError("one of yield to maturity or bond price must be provided")

    def __refresh_value_cache__(self):
        """
        recompute cached bond properties.

        :return:
        """
        self.__refresh_primary_cache__()
        self._d2Bdy2 = np.sum(self._ts * self._ts * self._dcs)
        self._duration = -self._dBdy / self._B
        self._convexity = self._d2Bdy2 / self._B

    def __refresh_primary_cache__(self):
        """
        recompute frequently accessed bond properties except for duration, convexity and second order derivative.

        :return:
        """
        self._ts = np.arange(self._T, 0, -1 / self._m)[::-1]
        self._cs = [self._R * self._F / 100 / self._m for _ in range(len(self._ts))]
        self._cs[-1] += self._F
        self._dcs = np.exp(-self._y * self._ts) * self._cs

        self._B = np.sum(self._dcs)
        self._dBdy = float(np.sum(-self._ts * self._dcs))

    @property
    def T(self):
        """
        time to maturity

        :return:
        """
        return self._T

    @T.setter
    def T(self, value):
        self._T = value
        self.__refresh_value_cache__()

    @property
    def m(self):
        """
        coupon payments per year

        :return:
        """
        return self._m

    @m.setter
    def m(self, value):
        self._m = value
        self.__refresh_value_cache__()

    @property
    def R(self):
        """
        coupon rate

        :return:
        """
        return self._R

    @R.setter
    def R(self, value):
        self._R = value
        self.__refresh_value_cache__()

    @property
    def y(self):
        """
        yield to maturity

        :return:
        """
        return self._y

    @property
    def ytm(self):
        """
        yield to maturity

        :return:
        """
        return self._y

    @property
    def current(self):
        """
        the current yield of the bond
        = annual interest payment / bond price

        :return:
        """
        return self._R / 100 * self._F / self._B

    @property
    def bankeq(self):
        """
        the bank equivalent yield of the bond
        = (par - value) / par * 360 / days to maturity

        :return:
        """
        assert self._R == 0
        return (self._F - self._B) / self._F * 360 / (self._T * 365)

    @property
    def cdeq(self):
        """
        the money market equivalent yield of the bond
        = (par - value) / value * 360 / days to maturity

        :return:
        """
        assert self._R == 0
        return (self._F - self._B) / self._B * 360 / (self._T * 365)

    @y.setter
    def y(self, value):
        self._y = value
        self.__refresh_value_cache__()

    @property
    def B(self):
        """
        bond value

        :return:
        """
        return self._B

    @B.setter
    def B(self, value):
        def f(x):
            self._y = x
            self.__refresh_primary_cache__()
            return self._B - value

        def df(x: float):
            return self._dBdy

        try:
            # compute implied yield to maturity with initial guess = 0.1
            self.y = root(f, 0.1, df, epsilon=10e-9, delta=10e-9)
        except RuntimeError:
            logging.error("invalid bond value")

    @property
    def F(self):
        """
        face value

        :return:
        """
        return self._F

    @F.setter
    def F(self, value):
        self._F = value
        self.__refresh_value_cache__()

    @property
    def duration(self):
        """
        modified duration of the bond

        :return:
        """
        return self._duration

    @property
    def convexity(self):
        """
        :convexity of the bond

        :return:
        """
        return self._convexity


def find_curve(bond, known: np.array, epsilon: float = 10e-10):
    t = np.arange(bond.T, 0, - 1. / bond.m)[::-1]
    c = np.array([bond.R / bond.m] * len(t))
    c[-1] += bond.F

    def f(x: float) -> float:
        r = np.linspace(x, known[-1], len(t) + 1 - len(known), endpoint=False)[::-1]
        rr = np.concatenate([known[1:], r])
        return float(np.sum(c * np.exp(-rr * t))) - bond.B

    def df(x: float) -> float:
        r = np.linspace(x, known[-1], len(t) + 1 - len(known), endpoint=False)[::-1]
        cc = c[len(known) - 1:]
        tt = t[len(known) - 1:]
        return float(np.sum(-cc * tt * np.exp(-r * tt) * np.arange(1, len(tt) + 1) / len(tt)))

    x = root(f, 0.05, df=df)
    r = np.linspace(x, known[-1], len(t) + 1 - len(known), endpoint=False)[::-1]
    rr = np.concatenate([known[:], r])
    return rr


def bootstrap(bonds: Sequence[Bond], overnight: float, epsilon: float = 10e-10):
    """
    Bootstrap a zero rate curve from the given bonds and bond values.
    Note that the bonds must have equal coupon payment periods (equal <m>s).
    Zero rates at times for which we do not have a bond are calculated
    by a linear line connecting the two nearest rates at times for which we do have a bond.

    :param overnight:
    :param epsilon:
    :param bonds:
    :return:
    """
    bonds = sorted(bonds, key=lambda x: x.T)
    known = [overnight]
    for bond in bonds:
        known = find_curve(bond, known)
    return known
