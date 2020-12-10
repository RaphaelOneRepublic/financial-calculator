# Numerical Finance Stuff

Inspired by the Advanced Calculus for Financial Engineering Applications seminar at Baruch College but more.

## Numerical Implementations of Several Dull Routines
### Integration
1. Simpson's for Numerical Integration
### Optimization
1. Newton's method
2. Secant method

## Encapsulation of Plain Vanilla European Options
Using Black-Scholes model to evaluate the following for plain vanilla European options.
The Greeks are not implemented in the most efficient way, but in the most straightforward manner. 
Performance optimization is possible by taking into account the repeated computation of parameters d1 and d2, and the cumulative probability thereof.

|   d/dx*  | underlying price | implied volatility | time to maturity | risk-free interest rate |
| ------- | ---------------- | ------------------ | ---------------- | ----------------------- |
| premium | delta            | vega               | theta            | rho                     |
| delta   | gamma            | vanna              | charm            |                         |
| vega    | vanna            | vomma              | veta             |                         |
| gamma   | speed            | zomma              | color            |                         |
| vomma   |                  | ultima             | totto**          |                         |
| rho     |                  | vera**             |                  |                         |

\* each cell is obtained by taking the derivative of the row with respect to the column header

\** not implemented

For call options
```python
# N(x)
#      = norm.cdf(x)
# d1 = (ln(S/K) + (r - q + 1 / 2 * sigma ** 2) * T) / (sigma * sqrt(T))
# d2 = (ln(S/K) + (r - q - 1 / 2 * sigma ** 2) * T) / (sigma * sqrt(T))


# European call options using Black-Scholes formula
# premium = S * e ** (-q * T) * N(d1) - K * e ** (-r * T) * N(d2) 
```
### Example
```python
from calc.option import Vanilla
# compute the implied volatility of a plain vanilla ATM put option
# struck at $100, expiring in 6 months, with risk-free rate being 5% and traded at $5
option = Vanilla(100, 100, 0.5, 0.05,put=True, price=5)
print(option.implied)
# and compute its theta and gamma
print(option.theta, option.gamma)
```
## Bond Math
### Yield to Maturity and Bond Value
Yield to maturity can be calculated from traded price, and vice versa.
### Modified Duration, Convexity
In the continuous case, Macaulay duration is the same as modified duration and therefore not implemented.
Same for convexity
### Yields for Discount Bond
Current yield, bank-equivalent yield, and CD-equivalent yield
### Bootstrapping for Finding Zero Rate Curves
Linear interpolation assumed.

*TODO*: Smooth the curve (equalizing derivatives of curves at each point)

### Example
```python
from calc.bond import Bond, bootstrap
# compute the value of a semiannual coupon bond expiring in 4 years 
# with coupon rate 4% with yield to maturity 0.4
bond = Bond(4, 4, y=0.4)
print(bond.B)
# compute the yield to maturity of a semiannual coupon bond expiring in 4 years 
# with coupon rate 4% with value 101
bond.B = 101
print(bond.ytm)
# compute the above bond's duration and convexity
print(bond.duration, bond.convexity)

# bootstrap zero rate curves using four semiannual coupon bonds
# with 4% coupon rate and expiring in 6 months, 1, 2, and 5 years respectively
bonds = [
    Bond(0.5, 0, B=99),
    Bond(1, 4, B=102),
    Bond(2, 4, B=103.5),
    Bond(5, 4, B=109)
]
print(bootstrap(bonds))
```

## TIme Value of Money
Compute the present value / future value of a set of cash flows
Convert between different compound periods
Level-payment of fully amortized mortgage

