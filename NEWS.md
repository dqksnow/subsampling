# subsampling 0.4.0

* Expanded `ssp.glm.rF()` with a more comprehensive rare-feature subsampling implementation, including additional criteria, automatic rare-feature detection, response balancing options, and union-based combination.
* Improved summaries for rare-feature GLMs and rare-event logistic regression, including clearer sample-composition information.
* Added broader tests for rare-feature GLM scenarios and updated rare-event logistic regression summary tests.
* Revised package vignettes, including the rare-feature GLM vignette and wording cleanup across method vignettes.

# subsampling 0.3.0

* Added `ssp.glm.rF()`. It provides balanced subsampling methods to deal with rare features in GLMs.

# subsampling 0.2.0

## Changes to `ssp.glm()`
* Unified the representation of weights in the square score and information matrices  to improve readability.
* Added support for `family = "gaussian"`.
* Updated the default value of `alpha` from `0` to `0.01`.
* Added a new output element `comp.time` that reports computation time.

# subsampling 0.1.1

# subsampling 0.1.0

* Initial CRAN submission.

* First Version submitted to CRAN.

* After this submission, all updates will be logged here.
