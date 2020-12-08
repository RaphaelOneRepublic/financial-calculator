# financial-calculator

Inspired by the Advanced Calculus for Financial Engineering Applications seminar at Baruch College but more.

## Numerical Implementations of Several Dull Routines
### Integration
1. Simpson's for Numerical Integration
### Optimization
1. Newton's method
2. Secant method
3. Bisection method (converges slowly)

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

