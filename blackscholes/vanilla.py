import numpy as np
from scipy.stats import norm

from blackscholes import numerical


class Vanilla(object):
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

    def __init__(self, S, K, r, T, sigma, direction='C', q=0.0):
        """
        construct a plain vanilla option object


        :param S: currently prevailing spot price of underlying
        :param K: strike price
        :param r: risk-free interest rate
        :param T: time to maturity
        :param sigma: black-scholes price volatility
        :param direction: one of 'C' or 'P'
        :param q: continuously distributed dividend rate
        """
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
        return (np.log(self.S / self.K) + (self.r - self.q - 0.5 * self.sigma * self.sigma) * self.T) / \
               (self.sigma * np.sqrt(self.T))

    def price(self):
        """
        compute the fair price of the option


        :return:
        """
        if self.direction == 'C' or self.direction == 'c' or self.direction == 'call':
            return self.S * np.exp(-self.q * self.T) * norm.cdf(self.d1()) \
                   - self.K * np.exp(-self.r * self.T) * norm.cdf(self.d2())
        elif self.direction == 'P' or self.direction == 'p' or self.direction == 'put':
            return - self.S * np.exp(-self.q * self.T) * norm.cdf(-self.d1()) \
                   + self.K * np.exp(-self.r * self.T) * norm.cdf(-self.d2())

    def delta(self):
        """
        compute the sensitivity of option price with regard to the underlying price


        :return:
        """
        if self.direction == 'C' or self.direction == 'c' or self.direction == 'call':
            return np.exp(-self.q * self.T) * norm.cdf(self.d1())
        elif self.direction == 'P' or self.direction == 'p' or self.direction == 'put':
            return - np.exp(-self.q * self.T) * norm.cdf(-self.d1())

    def gamma(self):
        """
        compute the sensitivity of option delta with regard to the underlying price,
        or the second derivative of option price with regard to underlying price


        :return:
        """
        return self.K * np.exp(-self.r * self.T) * norm.cdf(self.d2()) / \
               (self.S * self.S * self.sigma * np.sqrt(self.T))

    def vega(self):
        """
        compute the sensitivity of option price with regard to volatility


        :return:
        """
        return self.K * np.exp(-self.r * self.T) * norm.pdf(self.d2()) * np.sqrt(self.T)

    def theta(self):
        """
        the time decay
        compute the sensitivity of option price with regard to time to maturity


        :return:
        """
        if self.direction == 'C' or self.direction == 'c' or self.direction == 'call':
            return - np.exp(-self.q * self.T) * (self.S * norm.pdf(self.d1()) * self.sigma) / (2 * np.sqrt(self.T)) \
                   - self.r * self.K * np.exp(-self.r * self.T) * norm.cdf(self.d2()) \
                   + self.q * self.S * np.exp(-self.q * self.T) * norm.cdf(self.d1())
        elif self.direction == 'P' or self.direction == 'p' or self.direction == 'put':
            return - np.exp(-self.q * self.T) * (self.S * norm.pdf(-self.d1()) * self.sigma) / (2 * np.sqrt(self.T)) \
                   + self.r * self.K * np.exp(-self.r * self.T) * norm.cdf(-self.d2()) \
                   - self.q * self.S * np.exp(-self.q * self.T) * norm.cdf(-self.d1())

    def rho(self):
        """
        compute the sensitivity of option price with regard to risk-free interest rate


        :return:
        """
        if self.direction == 'C' or self.direction == 'c' or self.direction == 'call':
            return self.K * self.T * np.exp(-self.r * self.T) * norm.cdf(self.d2())
        elif self.direction == 'P' or self.direction == 'p' or self.direction == 'put':
            return - self.K * self.T * np.exp(-self.r * self.T) * norm.cdf(-self.d2())

    def omega(self):
        """
        also called lambda / gearing
        compute the percentage change in option value per percentage change in underlying price


        :return:
        """
        return self.delta() * self.S / self.price()

    def vanna(self):
        """
        compute the sensitivity of option delta with regard to sigma
        also the second order derivative of option price with regard to spot price and sigma


        :return:
        """

        return self.vega() * (1 - self.d1() / (self.sigma * np.sqrt(self.T))) / self.S

    def charm(self):
        raise NotImplementedError()

    def implied_volatility(self, market: float):
        """
        compute the implied volatility at the provided market price
        This will put the computed volatility in place of the originally stored value


        :param market: trending market price of option
        :return:
        """

        def target(x: float):
            self.sigma = x
            return self.price() - market

        def derivative(x: float):
            self.sigma = x
            return self.vega()

        return numerical.newtons(target, derivative, 0.1)


if __name__ == '__main__':
    print(Vanilla(100, 100, 0.05, 0.5, 0, 'P', q=0.1).implied_volatility(market=9.396990676896607))
