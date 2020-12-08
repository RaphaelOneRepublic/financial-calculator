import numpy as np
from typing import Union, List
from calc import numerical


class Cashflow(object):
    flows: np.ndarray  # cash flows

    times: np.ndarray  # times at which the flows are received

    rates: np.ndarray  # discount factor

    def __init__(self, flow, time, rate: Union[float, List[float]] = 0.0):
        if isinstance(rate, float) or rate == 0:
            rate = [rate] * len(flow)
        self.flows = np.array(flow, dtype=float)
        self.times = np.array(time, dtype=float)
        self.rates = np.array(rate, dtype=float)

        if (len(self.flows) != len(self.times)) or \
                (len(self.flows) != len(self.rates)) or \
                (len(self.rates) != len(self.times)):
            raise ValueError("lengths of input vectors do not match")

    def present_value(self):
        """
        compute the present value of the cash flow
        :return:
        """
        return np.sum(self.flows * np.exp(- self.times * self.rates))

    def internal_rate_of_return(self, pv: float = None):
        """
        compute the internal rate of rate for the cash flow

        :param pv: present value of cash flow, if provided, onw copy of rates would be ignored
        :return:
        """
        if pv is None:
            pv = self.present_value()

        def func(x):
            rates = np.array([x] * len(self.flows))
            return np.sum(self.flows * np.exp(- self.times * rates)) - pv

        def partial(x):
            rates = np.array([x] * len(self.flows))
            return np.sum(- self.times * self.flows * np.exp(- self.times * rates))

        return numerical.newtons(func, partial, 0.1)
