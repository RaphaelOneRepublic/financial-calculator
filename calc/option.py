import logging

import numpy as np
from scipy.stats import norm

from calc.optimize import root


class Vanilla(object):
    """
    Represents a plain vanilla European option, either e a call or a put.
    Upon creation, the underlying spot price, the strike price, the risk-free interest rate ust be provided.
    A put can be specified by setting the indicating boolean value to True,
    The volatility can be provided as <sigma>, which will shadow the argument <price>
    If <sigma> is not provided, <price> must be provided in its place, where the volatility is implied
    """

    def __init__(self, S: float, K: float, T: float, r: float, sigma: float = None, q: float = 0,
                 put: bool = False, price: float = None):
        """
        construct a plain vanilla European option
        :param S: underlying spot price
        :param K: strike price
        :param T: time to maturity
        :param r: continuously compounded risk-free interest rate
        :param sigma: (implied) black-scholes volatility
        :param q: continuously distributed dividend rate
        :param put: whether the option is a put
        :param price: traded price
        """
        self._S = S
        self._K = K
        self._T = T
        self._r = r
        self._q = q
        self._put = put

        if sigma is not None:  # use sigma to compute premium
            self._sigma = sigma
            self.__refresh_value_cache__()
        elif price is not None:
            self.premium = price  # use premium to compute implied volatility
        else:
            raise ValueError("one of implied volatility or traded price must be present")

    def __refresh_value_cache__(self):
        """
        recompute frequently accessed cache values
        :return:
        """
        # black scholes coefficients
        self._rootT = np.sqrt(self._T)
        self._d1 = (np.log(self.S / self.K) + (self._r - self._q + 0.5 * self._sigma ** 2) * self._T) \
                   / (self._sigma * self._rootT)
        self._d2 = self._d1 - self._sigma * self._rootT

        # discount factor of dividend and risk-free rate respectively
        self._expnqt = np.exp(-self._q * self._T)
        self._expnrt = np.exp(-self._r * self._T)

        # present value of S & K, discounted at dividend and risk-free rate respectively
        self._pvs = self._S * self._expnqt
        self._pvk = self._K * self._expnrt

        # the following three sections are related to the standard normal distribution
        # probability density of ds
        self._nd1 = norm.pdf(self._d1)
        self._nd2 = norm.pdf(self._d2)

        # cumulative density of ds
        self._Nd1 = norm.cdf(self._d1)
        self._Nd2 = norm.cdf(self._d2)

        # cumulative density of negative ds
        self._Nnd1 = 1 - self._Nd1
        self._Nnd2 = 1 - self._Nd2

    @property
    def S(self):
        return self._S

    @S.setter
    def S(self, value):
        self._S = value
        self.__refresh_value_cache__()

    @property
    def K(self):
        return self._K

    @K.setter
    def K(self, value):
        self._K = value
        self.__refresh_value_cache__()

    @property
    def T(self):
        return self._T

    @T.setter
    def T(self, value):
        self._T = value
        self.__refresh_value_cache__()

    @property
    def r(self):
        return self._r

    @r.setter
    def r(self, value):
        self._r = value
        self.__refresh_value_cache__()

    @property
    def sigma(self):
        return self._sigma

    @property
    def implied(self):
        return self._sigma

    @sigma.setter
    def sigma(self, value):
        self._sigma = value
        self.__refresh_value_cache__()

    @property
    def q(self):
        return self._q

    @q.setter
    def q(self, value):
        self._q = value
        self.__refresh_value_cache__()

    @property
    def put(self):
        return self._put

    @put.setter
    def put(self, value):
        self._put = value

    @property
    def premium(self):
        if not self._put:
            return self._pvs * self._Nd1 \
                   - self._pvk * self._Nd2
        else:
            return - self._pvs * self._Nnd1 \
                   + self._pvk * self._Nnd2

    @premium.setter
    def premium(self, value):
        def f(x):
            self._sigma = x
            self.__refresh_value_cache__()
            return self.premium - value

        def df(x):
            return self.vega

        try:
            self.sigma = root(f, 0.1, df)
        except RuntimeError:
            logging.error("invalid option price ")

    @property
    def delta(self):
        """
        first order derivative of option value with respect to underlying spot price
        :return:
        """
        if not self._put:
            return self._expnqt * self._Nd1
        else:
            return -self._expnqt * self._Nnd1

    @property
    def vega(self):
        """
        first order derivative of option value with respect to implied volatility
        :return:
        """
        return self._pvs * self._rootT * self._nd1

    @property
    def gamma(self):
        """
        second order derivative of option value twice with respect to underlying spot price
        :return:
        """
        return self._expnqt / self._S / self._sigma / self._rootT * self._nd1

    @property
    def theta(self):
        """
        first order derivative of option value with respect to time to maturity
        :return:
        """
        if not self._put:
            return - self._sigma * self._pvs / 2 / self._rootT * self._nd1 \
                   + self._q * self._pvs * self._Nd1 \
                   - self._r * self._pvk * self._Nd2
        else:
            return - self._sigma * self._pvs / 2 / self._rootT * self._nd1 \
                   - self._q * self._pvs * self._Nnd1 \
                   + self._r * self._pvk * self._Nnd2

    @property
    def rho(self):
        """
        first order derivative of option value with respect to risk-free interest rate
        :return:
        """
        if not self._put:
            return self._T * self._pvk * self._Nd2
        else:
            return -self._T * self._pvk * self._Nnd2

    # TODO Other Greeks
