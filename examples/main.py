from calc.bond import Bond, bootstrap

if __name__ == '__main__':
    # # 0.2
    # # 0.3034547400773509
    # # 0.3032954578149987
    # # 0.3033023937192005
    # # 0.30330239401059217
    # option = Vanilla(50, 45, 9 / 12, 0.02, q=0.01, price=8)
    # print(option.implied)

    # [0.015      0.01613126 0.01726252 0.02027111 0.0183568  0.01644248
    #  0.01452817]
    #
    # 0.012593315743269974
    # 0.014522623446027394
    # 0.01452817173607332
    # 0.014528171781783183

    bonds = [
        Bond(1, 3, B=101.25),
        Bond(1.5, 2, B=99.95),
        Bond(3, 5, B=110.3),
    ]
    curves = bootstrap(bonds, 0.015)
    print(curves)

    # T = 28 / 12
    # r = 4
    # m = 2
    #
    # ts = np.arange(T, 0, -1 / m)[::-1]
    # cs = [r / m for _ in range(len(ts))]
    # cs[-1] += 100
    # rs = 0.015 + (1 + 2 * ts * ts) / (100 + 100 * ts * ts)
    #
    # dcs = np.exp(-rs * ts) * cs
    # B = sum(dcs)
    #
    # print(ts)
    # print(cs)
    # print(rs)
    # print(np.exp(-rs * ts))
    #
    # bond = Bond(T, 4, B=B)
    # print(B)
    # print(bond.ytm)
    # print(bond.duration * bond.B)
    # print(bond.convexity)

    # sigma = 0.25
    # S = 50
    # q = 0.01
    # r = 0.03
    # T = 0.5
    #
    # option = Vanilla(S, 50, T, r, sigma, q=q, put=True)
    #
    #
    # def f(K):
    #     option.K = K
    #     return option.premium - K + S
    #
    # def df(K):
    #     return option.dVdK - 1
    #
    # k = root(f, 50, df, delta=10e-6, epsilon=10e-6,progress=True)
    # print(k)
