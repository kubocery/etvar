#' Daily US market data around the 2025 tariff announcements
#'
#' Seven daily financial-market series from 4 January 2021 to 15 May 2025,
#' suitable for estimating an event-targeted VAR around the US tariff
#' announcement of 2 April 2025 (first trading-day reaction on 3 April 2025).
#' The variables mirror those used by Pinter, Smets and \enc{Üslü}{Uslu} (2025);
#' the oil price and five-year breakeven inflation rate stand in for the
#' copper price and one-year inflation swap rate used in the paper.
#'
#' These substitutes are not equivalent to the paper's series, so results based
#' on this dataset will differ from the published ones. In particular,
#' `Breakeven5Y` is not a proxy for one-year inflation expectations: it
#' averages expected inflation over five years, moves closely with oil prices,
#' and includes a TIPS liquidity premium that widens in risk-off episodes. On
#' 3 April 2025 the one-year inflation swap rate rose sharply while the
#' five-year breakeven rate was little changed and fell over the following
#' days, so the identified tariff shock lowers `Breakeven5Y` here, whereas the
#' paper finds that it raises short-term inflation expectations.
#'
#' Observations are kept only on days on which every series is available, so
#' days on which US equity, bond, foreign-exchange or oil markets were closed
#' are dropped.
#'
#' @format A data frame with 1086 rows and 8 columns:
#' \describe{
#'   \item{date}{Trading date (`Date`).}
#'   \item{SP500}{Log of the S&P 500 index close.}
#'   \item{VIX}{Log of the CBOE Volatility Index (VIX) close.}
#'   \item{Yield1Y}{One-year Treasury constant-maturity yield, percent
#'     (FRED series `DGS1`).}
#'   \item{TermSpread}{Ten-year minus one-year Treasury constant-maturity
#'     yield, percentage points (`DGS10 - DGS1`).}
#'   \item{Oil}{Log of the WTI crude oil spot price, US dollars per barrel
#'     (`DCOILWTICO`).}
#'   \item{Breakeven5Y}{Five-year breakeven inflation rate: five-year nominal
#'     minus five-year TIPS constant-maturity yield, percentage points
#'     (`DGS5 - DFII5`).}
#'   \item{USDEUR}{Log of US dollars per euro, noon buying rate in New York
#'     (`DEXUSEU`).}
#' }
#'
#' @source
#' * `SP500`, `VIX`: Yahoo Finance, <https://finance.yahoo.com>, symbols
#'   `^GSPC` and `^VIX`, daily closing values.
#' * `Yield1Y`, `TermSpread`, `Breakeven5Y`: Board of Governors of the Federal
#'   Reserve System, H.15 Selected Interest Rates, retrieved from FRED, Federal
#'   Reserve Bank of St. Louis, <https://fred.stlouisfed.org>.
#' * `USDEUR`: Board of Governors of the Federal Reserve System, H.10 Foreign
#'   Exchange Rates, retrieved from FRED.
#' * `Oil`: US Energy Information Administration, Spot Prices, retrieved from
#'   FRED.
#'
#' Retrieved on 17 November 2026.
#'
#' @references
#' Pinter, G., Smets, F. and \enc{Üslü}{Uslu}, S. (2025). Market whiplash after
#' the 2025 tariff shock: an event-targeted VAR approach. *BIS Working Papers*,
#' No. 1282. Bank for International Settlements.
#' <https://www.bis.org/publications/working-paper-1282-market-whiplash-after-2025-tariff-shock-event-targeted-var-approach>
#'
#' @examples
#' data("etvar_example")
#' str(etvar_example)
#' \donttest{
#' fit <- etvar(as.matrix(etvar_example[, -1]), etvar_example$date,
#'              target_date = "2025-04-03", lags = 8)
#' fit
#' }
"etvar_example"
