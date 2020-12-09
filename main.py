from calc.bond import Bond
from calc.option import Vanilla

if __name__ == '__main__':
    bond = Bond(3, 2, 4, B=101)
    print(bond.ytm)
    print(bond.duration)
    print(bond.convexity)

    option = Vanilla(100, 100, 1, 0.05, price=8)
    print(option.implied)
