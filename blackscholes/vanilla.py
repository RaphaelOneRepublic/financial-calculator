import numpy as np
from scipy.stats import norm


class PlainVanilla(object):
    """
    Class to represent a plain vanilla option
    """
    S: float  # spot price of the underlying

    K: float  # strike price of the option

    r: float  # prevailing interest rate

    q: float  # continuously distributed dividend rate

    T: float  # time to maturity

    sigma: float  # price volatility

    direction: str  # option direction 'C' | 'P'

    def __init__(self, S, K, r, T, sigma, direction='C', q=0):
        self.q = q
        self.direction = direction
        self.sigma = sigma
        self.T = T
        self.r = r
        self.K = K
        self.S = S

    def d1(self):
        return (np.log(self.S / self.K) + (self.r - self.q + 0.5 * self.sigma * self.sigma) * self.T) / \
               (self.sigma * np.sqrt(self.T))

    def d2(self):
        return (np.log(self.S / self.K) + (self.r - self.q + 0.5 * self.sigma * self.sigma) * self.T) / \
               (self.sigma * np.sqrt(self.T))

    def price(self):
        if self.direction == 'C' or self.direction == 'c' or self.direction == 'call':
            return self.S * np.exp(-self.q * self.T) * norm.cdf(self.d1()) \
                   - self.K * np.exp(-self.r * self.T) * norm.cdf(self.d2())
        elif self.direction == 'P' or self.direction == 'p' or self.direction == 'put':
            return - self.S * np.exp(-self.q * self.T) * norm.cdf(-self.d1()) \
                   + self.K * np.exp(-self.r * self.T) * norm.cdf(-self.d2())

    def delta(self):
        if self.direction == 'C' or self.direction == 'c' or self.direction == 'call':
            return norm.cdf(self.d1())
        elif self.direction == 'P' or self.direction == 'p' or self.direction == 'put':
            return norm.cdf(self.d2())

    def gamma(self):
        return self.K * np.exp(-self.r * self.T) * norm.cdf(self.d2()) / \
               (self.S * self.S * self.sigma * np.sqrt(self.T))

    def vega(self):
        return self.K * np.exp(-self.r * self.T) * norm.cdf(self.d2()) * \
               np.sqrt(self.T)
