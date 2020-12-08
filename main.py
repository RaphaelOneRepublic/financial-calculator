from calc.vanilla import Vanilla
from calc.bond import Bond

S = 40
K = 40
r = 0.025
T = 5 / 12
sigma = 0.1
q = 0.01

if __name__ == '__main__':
    print(Bond(3, 0.5, 4, 101).modified_duration())
    print(Bond(3, 0.5, 4, 101).macaulay_duration())
    print(Bond(3, 0.5, 4, 101).convexity())
