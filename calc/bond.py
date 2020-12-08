import numpy as np

from calc.cashflow import Cashflow


class Bond(object):
    """
    a coupon paying bond
    """
    T: float  # time to maturity
    period: float  # coupon payment period in years
    coupon: float  # coupon rate in percentage
    F: float  # face value
    B: float  # bond value

    def __init__(self, T, period, coupon, B, F=100):
        """
        a coupon bond

        :param T: time to maturity
        :param period: coupon paying periods in years
        :param coupon: annual coupon rate in percentage
        :param (optional) F: face value, default to 100
        """
        self.T = T
        self.period = period
        self.coupon = coupon
        self.F = F
        self.B = B

        self.times = np.arange(self.T, 0, -self.period)[::-1]
        self.cash = [self.coupon * self.F / 100 * self.period for _ in range(len(self.times))]
        self.cash[-1] += self.F

    def yield_to_maturity(self):
        """
        compute the yield to maturity of the bond

        :return:
        """
        return Cashflow(flow=self.cash, time=self.times, rate=0).internal_rate_of_return(pv=self.B)

    def modified_duration(self):
        """
        compute the modified duration of the bond
        this is the same as the Macaulay duration assuming continuously compounded interest


        :return:
        """
        return np.sum(- self.times * self.cash
                      * np.exp(- self.times * np.array([self.yield_to_maturity()] * len(self.cash)))) / (-self.B)

    def macaulay_duration(self):
        """
        compute the macaulay duration of the bond
        this is the same as the modified duration assuming continuously compounded interest


        :return:
        """
        return np.sum(self.cash * np.exp(- self.times * self.yield_to_maturity()) / self.B * self.times)

    def convexity(self):
        """
        compute the convexity of the bond


        :return:
        """
        return np.sum(self.times * self.times * self.cash * np.exp(
            - self.times * np.array([self.yield_to_maturity()] * len(self.cash)))) / self.B
