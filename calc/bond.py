import numpy as np
from calc.numerical import newtons


class Bond(object):
    """
    Represents a coupon paying bond
    Upon creation, the time to maturity, coupon periods per year, coupon rate must be provided
    If yield to maturity is provided, bond value would be ignored.
    If yield to maturity is provided, bond value would be used to compute the implied yield to maturity
    Face value is assumed to be 100 if not provided
    """

    def __init__(self, T: float, m: int, R: float, y: float = None, F: float = 100, B: float = None):
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
            raise ValueError("one of yield to maturity of bond price must be provided")

    def __refresh_value_cache__(self):
        """
        recompute cached bond properties
        :return:
        """
        self.__refresh_primary_cache__()
        self._d2Bdy2 = np.sum(self._ts * self._ts * self._dcs)
        self._duration = -self._dBdy / self._B
        self._convexity = self._d2Bdy2 / self._B

    def __refresh_primary_cache__(self):
        """
        recompute frequently accessed bond properties except for duration, convexity and second order derivative
        :return:
        """
        self._ts = np.arange(self._T, 0, -1 / self._m)[::-1]
        self._cs = [self._R * self._F / 100 / self._m for _ in range(len(self._ts))]
        self._cs[-1] += self._F
        self._dcs = np.exp(-self._y * self._ts) * self._cs
        self._B = np.sum(self._dcs)
        self._dBdy = np.sum(-self._ts * self._dcs)

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
        def func(x):
            self._y = x
            self.__refresh_primary_cache__()
            return self._B - value

        def derivative(x):
            return self._dBdy

        self.y = newtons(func, derivative, 0.1)

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


if __name__ == '__main__':
    bond = Bond(3, 2, 4, B=101)
    print(bond.ytm)
    print(bond.duration)
    print(bond.convexity)
    Bond()
