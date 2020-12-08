import warnings

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
        self.sigma = sigma
        self.T = T
        self.r = r
        self.K = K
        self.S = S
        direction = direction.strip().capitalize()
        if direction == 'C' or direction == "CALL":
            self.direction = 'C'
        elif direction == 'P' or direction == "PUT":
            self.direction = 'P'

    def d1(self):
        return (np.log(self.S / self.K) + (self.r - self.q + 0.5 * self.sigma * self.sigma) * self.T) / \
               (self.sigma * np.sqrt(self.T))

    def d2(self):
        return (np.log(self.S / self.K) + (self.r - self.q - 0.5 * self.sigma * self.sigma) * self.T) / \
               (self.sigma * np.sqrt(self.T))

    def premium(self):
        """
        compute the fair price of the option


        :return:
        """
        if self.direction == 'C':
            return self.S * np.exp(-self.q * self.T) * norm.cdf(self.d1()) \
                   - self.K * np.exp(-self.r * self.T) * norm.cdf(self.d2())
        elif self.direction == 'P':
            return - self.S * np.exp(-self.q * self.T) * norm.cdf(-self.d1()) \
                   + self.K * np.exp(-self.r * self.T) * norm.cdf(-self.d2())

    def delta(self):
        """
        compute the sensitivity of option price with regard to the underlying price


        :return:
        """
        if self.direction == 'C':
            return np.exp(-self.q * self.T) * norm.cdf(self.d1())
        elif self.direction == 'P':
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
        if self.direction == 'C':
            return - np.exp(-self.q * self.T) * (self.S * norm.pdf(self.d1()) * self.sigma) / (2 * np.sqrt(self.T)) \
                   - self.r * self.K * np.exp(-self.r * self.T) * norm.cdf(self.d2()) \
                   + self.q * self.S * np.exp(-self.q * self.T) * norm.cdf(self.d1())
        elif self.direction == 'P':
            return - np.exp(-self.q * self.T) * (self.S * norm.pdf(-self.d1()) * self.sigma) / (2 * np.sqrt(self.T)) \
                   + self.r * self.K * np.exp(-self.r * self.T) * norm.cdf(-self.d2()) \
                   - self.q * self.S * np.exp(-self.q * self.T) * norm.cdf(-self.d1())

    def rho(self):
        """
        compute the sensitivity of option price with regard to risk-free interest rate


        :return:
        """
        if self.direction == 'C':
            return self.K * self.T * np.exp(-self.r * self.T) * norm.cdf(self.d2())
        elif self.direction == 'P':
            return - self.K * self.T * np.exp(-self.r * self.T) * norm.cdf(-self.d2())

    def omega(self):
        """
        also called lambda / gearing
        compute the percentage change in option value per percentage change in underlying price


        :return:
        """
        return self.delta() * self.S / self.premium()

    def vanna(self):
        """
        compute the sensitivity of option delta with regard to sigma
        also the second order derivative of option price with regard to spot price and volatility


        :return:
        """

        return self.vega() * (1 - self.d1() / (self.sigma * np.sqrt(self.T))) / self.S

    def charm(self):
        """
        delta decay

        compute the sensitivity of option delta with regard to time to maturity
        also the second order derivative of option price with regard to spot price and time to maturity


        :return:
        """
        residual = - np.exp(-self.q * self.T) * norm.pdf(self.d1()) \
                   * (2 * (self.r - self.q) * self.T - self.d2() * self.sigma * np.sqrt(self.T)) / \
                   (2 * self.T * self.sigma * np.sqrt(self.T))
        if self.direction == 'C':
            return self.q * np.exp(-self.q * self.T) * norm.cdf(self.d1()) + residual
        elif self.direction == 'C':
            return -self.q * np.exp(-self.q * self.T) * norm.cdf(-self.d1()) + residual

    def vomma(self):
        """
        vega convexity

        compute the sensitivity of vega with regard to volatility
        also the second order derivative of option price with respect to volatility twice


        :return:
        """
        return self.vega() * self.d1() * self.d2() / self.sigma

    def veta(self):
        """
        compute the sensitivity of vega with regard to time to maturity
        also the second order derivative of option price with respect to volatility and time to maturity


        :return:
        """
        return - self.vega() * (self.q
                                + ((self.r - self.q) * self.d1()) / (self.sigma * np.sqrt(self.T))
                                - (1 + self.d1() * self.d2()) / (2 * self.T))

    def speed(self):
        """
        compute sensitivity of gamma with respect to the underlying price
        also the third order derivative of option price with respect to the underlying price thrice


        :return:
        """
        return - self.gamma() / self.S * (self.d1() / (self.sigma * np.sqrt(self.T)) + 1)

    def zomma(self):
        """
        compute the sensitivity of gamma with regard to volatility
        also the third order derivative of option price with respect once to volatility and twice to underlying price


       :return:
        """
        return self.gamma() * (self.d1() * self.d2() - 1) / self.sigma

    def color(self):
        """
        gamma decay

        compute the sensitivity of gamma with regard to time to maturity
        also the third order derivative of option price with respect once to time and twice to underlying price


        :return:
        """
        return - np.exp(-self.q * self.T) * norm.pdf(self.d1()) / (2 * self.S * self.T * self.sigma * np.sqrt(self.T)) \
               * (2 * self.q * self.T + 1
                  + (2 * (self.r - self.q) * self.T - self.d2() * self.sigma * np.sqrt(self.T))
                  / (self.sigma * np.sqrt(self.T)) * self.d1())

    def ultima(self):
        """
        compute the sensitivity of vomma with regard to volatility
        also the third order derivative of option price with respect tp volatility thrice


        :return:
        """
        return - self.vega() / (self.sigma ** 2) * \
               (self.d1() * self.d2() * (1 - self.d1() * self.d2()) + self.d1() ** 2 + self.d2() ** 2)

    def implied_volatility(self, market: float, method: str = 'newton', initial: float = 0.1):
        """
        compute the implied volatility at the provided market price
        This will put the computed volatility in place of the originally stored value

        :param market: trending market price of option
        :param method: numerical method used to find the implied volatility
        :param initial: hint for answer
        :return:
        """

        def target(x: float):
            self.sigma = x
            return self.premium() - market

        def derivative(x: float):
            self.sigma = x
            return self.vega()

        method = method.strip().lower()
        if method == 'newton':
            return numerical.newtons(target, derivative, initial)
        elif method == 'secant':
            return numerical.secant(target, initial, initial + 0.01)
        elif method == 'bisect':
            warnings.warn("the bisection method converges slowly")
            return numerical.bisect(target, 0.0001, initial)
        else:
            raise ValueError(f"unsupported numerical method {method}")


if __name__ == '__main__':
    print(Vanilla(100, 100, 0.05, 0.5, 0.3, q=0.01).implied_volatility(2.75, method='newton'))
