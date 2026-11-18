# etvar

Event-Targeted Vector Autoregression in R.

etvar estimates structural VAR models identified by targeting one or more
event dates. Starting from a reduced-form VAR, it looks for the rotation of
the Cholesky factor whose single structural shock best explains how all
variables moved on the event date(s). From the fitted model you get impulse
responses (with bootstrap bands), a structural forecast error variance
decomposition, and a historical decomposition that separates the event-day
shock from target shocks that arrive after it.

The method follows Pinter, Smets and Üslü (2025),
[*Market whiplash after the 2025 tariff shock: an event-targeted VAR approach*](https://www.bis.org/publications/working-paper-1282-market-whiplash-after-2025-tariff-shock-event-targeted-var-approach),
BIS Working Papers No. 1282.

## Installation

From CRAN:

```r
install.packages("etvar")
```

Development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("kubocery/etvar")
```

## Example

```r
library(etvar)

# Daily S&P 500, VIX, Treasury yields, WTI oil, breakeven inflation and EUR/USD
data("etvar_example")

fit <- etvar(
  data        = as.matrix(etvar_example[, -1]),
  dates       = etvar_example$date,
  target_date = "2025-04-03",
  lags        = 8,
  shock_index = 1,
  n_starts    = 4,
  seed        = 1
)

summary(fit)

plot(fit, type = "irfs")            # structural impulse responses
plot(fit, type = "fevd")            # target shock vs. all other shocks
plot(fit, type = "hd")              # historical decomposition
plot(fit, type = "counterfactual")  # realised vs. target-shock-only paths
```

Passing a vector of dates to `target_date` targets a window of events jointly.
`start_date` and `end_date` restrict the sample used to estimate the
reduced-form VAR.

## Data

The bundled `etvar_example` dataset combines the S&P 500 and VIX from Yahoo
Finance with public-domain series from FRED (Federal Reserve Board H.15 and
H.10 releases; US Energy Information Administration). See `?etvar_example`
for definitions and sources.

WTI oil and the five-year breakeven inflation rate replace the copper price
and one-year inflation swap rate used in the paper, which are not freely
available. The substitutes behave differently around the tariff
announcements; in particular, the five-year breakeven rate fell rather than
rose, so the example does not reproduce the paper's finding on short-term
inflation expectations.

## Citation

Run `citation("etvar")` for how to cite the method and the package.

## License

MIT © Jakub Cery
