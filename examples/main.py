from calc.bond import Bond

if __name__ == '__main__':
    bond = Bond(34, 8 + 5 / 8, F=1000, B=1065)
    print(bond.current)
