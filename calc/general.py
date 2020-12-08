import numpy as np


def continuous(rate: float, n: int = 1):
    """
    compute the continuously compounded interest rate from discretely compounded rate

    :param rate: discretely compounded interest rate
    :param n: compound periods
    :return:
    """
    return n * np.log(1 + 1 + rate / n)
