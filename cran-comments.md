## Submission

This is the first submission of etvar.

The package implements the event-targeted VAR method of Pinter, Smets and
Üslü (2025), BIS Working Papers No. 1282, which is cited in the DESCRIPTION,
the package documentation and inst/CITATION. The authors of the paper were
contacted and asked to be cited rather than listed as contributors.

## Test environments

* local Windows 11, R 4.6.1
* win-builder, R-devel and R-release
* macOS builder (mac.r-project.org), R-release

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Data

The example dataset `etvar_example` (data/etvar_example.rda, about 30 KB)
contains daily closing values of the S&P 500 and VIX indices obtained from
Yahoo Finance, and public-domain series from FRED (Federal Reserve Board H.15
and H.10 releases; US Energy Information Administration). All sources are
documented in ?etvar_example.
