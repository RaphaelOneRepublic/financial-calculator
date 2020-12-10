from calc.option import Vanilla

if __name__ == '__main__':
    option = Vanilla(100, 100, 0.5, 0.05, q=0.01, price=0.1)
    print(option.implied)
