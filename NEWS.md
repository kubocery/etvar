# etvar 0.1.0

* Initial CRAN submission.
* `etvar()` estimates an event-targeted VAR in one call: reduced-form VAR,
  target-shock optimisation, structural impulse responses, historical
  decomposition with an intervention path, and plotting helpers.
* The target shock is parameterised by `n - 1` direction angles; the event
  objective is precomputed once, so each optimiser evaluation is cheap. With
  two variables the single angle is found by a global grid search refined
  with Brent's method.
* The estimated direction is signed so that the target shock is positive on
  the first event date, and impulse responses are responses to a positive
  one-standard-deviation shock.
* Single event dates and multi-date target windows are supported, as is
  restricting the VAR estimation sample with `start_date` and `end_date`.
* The historical decomposition separates the event-day target shock from
  target shocks before and after it.
* `compute_fevd()` gives a structural forecast error variance decomposition
  (target shock versus all other shocks), and `compute_irf_confint()` gives
  residual-bootstrap confidence bands for impulse responses; normalised
  impulse responses get bands computed from normalised draws.
* New `etvar_example` dataset of daily US market data (S&P 500 and VIX from
  Yahoo Finance; Treasury yields, breakeven inflation, EUR/USD and WTI oil
  from FRED).
