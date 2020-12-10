from calc.bond import Bond, bootstrap

if __name__ == '__main__':
    bonds = [
        Bond(0.5, 0, B=99),
        Bond(1, 4, B=102),
        Bond(2, 4, B=103.5),
        Bond(5, 4, B=109)
    ]
    print(bootstrap(bonds))
