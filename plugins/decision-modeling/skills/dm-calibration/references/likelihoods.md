# Choosing and building the calibration likelihood

Read this when the targets are anything other than a simple proportion with a standard
error, when a target is overdispersed relative to the obvious sampling distribution, or
when calibrating a cancer natural history model.

The governing principle: **the likelihood should describe how the target was generated,
not how convenient it is to code.** The normal-with-known-SE default is an approximation
that happens to be right often enough that people stop noticing when it is not.

## Quick selection

| Target | Distribution | R |
|---|---|---|
| Proportion with numerator and denominator | Binomial | `dbinom(x = n_event, size = n_total, prob = mu, log = TRUE)` |
| Proportion, overdispersed (clustered, multi-registry) | Beta-binomial | see below |
| Counts, or rates with person-time | Poisson | `dpois(x = n_obs, lambda = mu, log = TRUE)` |
| Counts, variance > mean | Negative binomial | `dnbinom(x = n_obs, mu = mu, size = k, log = TRUE)` |
| Distribution across >2 categories | Multinomial | `dmultinom(x = v_n, prob = v_p, log = TRUE)` |
| Only a mean and an SE are available | Normal | `dnorm(x = value, mean = mu, sd = se, log = TRUE)` |
| Strictly positive, right-skewed summary | Lognormal | `dlnorm(...)` on the log scale |

## Counts and rates

Age-specific cancer incidence is the common case. Incidence is a **rate**, so the model's
predicted rate has to be multiplied by the person-time in the target stratum to become a
count that the Poisson can score:

```r
# v_rate_model: model-predicted incidence rate per person-year, by age group
# df_target$n_cases, df_target$pyr: observed cases and person-years per age group
v_mu   <- v_rate_model * df_target$pyr        # expected count = rate x exposure
n_llik <- sum(dpois(x = df_target$n_cases, lambda = v_mu, log = TRUE))
```

Forgetting the person-time offset is a silent error: the calibration still runs and simply
fits the wrong thing, because a rate compared against a count is off by a factor that
varies across strata.

**When the Poisson is too tight.** The Poisson fixes variance equal to the mean. Registry
incidence data are almost always more variable than that — geographic heterogeneity,
reporting artifacts, secular trends the model does not represent. Fitting a Poisson to
overdispersed data produces credible intervals that are far too narrow and a posterior
that is confidently wrong. Use the negative binomial, in the `mu`/`size` parameterization:

```r
# size (k) is the dispersion parameter: Var = mu + mu^2 / k
n_llik <- sum(dnbinom(x = df_target$n_cases, mu = v_mu, size = n_k, log = TRUE))
```

`size` controls how far from Poisson you are, and the direction is easy to get backwards:
**large `size` means less overdispersion.** As `size` goes to infinity the negative
binomial converges to the Poisson, which is a useful check:

```r
sum(dnbinom(v_obs, mu = v_exp, size = 1e8, log = TRUE)) ==
  sum(dpois(v_obs, lambda = v_exp, log = TRUE))          # TRUE to numerical precision
```

For reference, the PRE-CISE paper used `size = 5` for the SIR testbed and `size = 1` for
the COVID-19 case study; `size = 1` is strong overdispersion.

## Proportions across more than two categories

Stage distribution at diagnosis, or the distribution of a cohort across health states in a
cSTM trace, is multinomial rather than a set of independent binomials. Treating the
categories as independent binomials double-counts the information, because the proportions
are constrained to sum to one.

```r
# v_n_stage: observed counts by stage at one time point
# v_p_stage: model-predicted stage probabilities, same order, summing to 1
n_llik <- dmultinom(x = v_n_stage, prob = v_p_stage, log = TRUE)
```

**`dmultinom` silently renormalizes `prob` to sum to one.** If the model's predicted
probabilities do not sum to one — you conditioned on the wrong denominator, or dropped a
category — you get a finite log-likelihood and no warning whatsoever. Assert it:

```r
stopifnot(abs(sum(v_p_stage) - 1) < 1e-8)
```

`dmultinom` is not vectorized over time points. Loop or `mapply` across them, and remember
that each time point contributes one term, not one per category.

For overdispersed category counts (pooling registries with different case mixes), the
Dirichlet-multinomial is the multi-category analogue of the beta-binomial.

## Overdispersion as a calibrated parameter

The dispersion parameter is usually unknown. Two options, and the choice matters:

**Fix it** at a defensible value and report a sensitivity analysis over it. Defensible when
you have external information — a published variance estimate, or a dispersion estimated by
fitting a GLM to the target data alone, outside the model.

**Calibrate it** as an extra parameter alongside the natural history parameters. This is
the honest choice when nothing external pins it down. Practically:

```r
v_param_names <- c("p_onset", "p_progress", "hr_mort", "n_k_inc")  # dispersion included
v_lb <- c(..., n_k_inc = 0.1)
v_ub <- c(..., n_k_inc = 50)
# Sample the dispersion on the log scale: it spans orders of magnitude and a
# uniform prior on the natural scale puts almost all mass on "nearly Poisson".
```

Three things to know before adding it:

- **It is weakly identified** unless there are many target points. With a handful of
  strata there is little information to separate "the model fits badly" from "the data are
  overdispersed," and the sampler will happily inflate dispersion to absorb misfit.
- **Inflating dispersion is a way of ignoring your targets.** Watch for a posterior
  pushing toward maximum overdispersion; that usually means structural misfit, not noisy
  data. Check the posterior predictive before accepting it.
- **Include it in the identifiability analysis** (see `precise-workflow.md`) along with
  everything else. It is a calibrated parameter and can be collinear with the rest.

Beta-binomial, with no extra package, parameterized by the mean `prob` and an
intra-cluster correlation `rho` in (0, 1):

```r
dbetabinom_log <- function(x, size, prob, rho) {
  a <- prob * (1 - rho) / rho
  b <- (1 - prob) * (1 - rho) / rho
  lchoose(size, x) + lbeta(x + a, size - x + b) - lbeta(a, b)
}
# rho -> 0 recovers the binomial exactly; verify this when you write it
```

## Model discrepancy

Everything above assumes that with the right parameters the model is *correct*, and that
all remaining difference from the targets is sampling noise. That is false for every
natural history model. Systematic model-target difference gets absorbed into the
parameters, biasing them, and the posterior becomes overconfident because the likelihood
has no way to express "the model is a bit wrong here."

The Kennedy–O'Hagan formulation makes the discrepancy explicit:

    y_observed(x) = f(x, theta) + delta(x) + epsilon

where `delta(x)` is a systematic discrepancy term with its own parameters, estimated
jointly with `theta`.

Two pragmatic versions, in increasing order of ambition:

1. **A variance inflation term.** Add a calibrated parameter that widens the target
   variance beyond sampling error — the negative binomial dispersion above is exactly this
   for counts. Cheap, and it stops the posterior claiming precision the model does not
   have. It does not remove bias.
2. **A structured discrepancy** — for example an age-varying multiplicative bias on
   incidence, smooth in age, with its own small number of parameters.

**The caveat that matters, and it is severe.** Calibration parameters and a flexible
discrepancy term are generically confounded: the discrepancy can absorb almost any amount
of parameter misspecification, so you trade bias for non-identifiability. Brynjarsdóttir
and O'Hagan (2014) showed that useful inference requires genuinely informative priors that
constrain the discrepancy — on its smoothness, its magnitude, or its sign. A vague prior
on `delta` produces wide, honest-looking, useless posteriors.

So: run the collinearity analysis *with the discrepancy parameters included*. If adding
them pushes the index above threshold, the discrepancy is not identified from these
targets and needs either a tighter prior or removal.

Do not introduce a discrepancy term as a first response to poor fit. Check first that the
misfit is not a target alignment bug (`SKILL.md` Step 1), a wrong likelihood family, or a
structural error in the model. Discrepancy is the right tool for known, irreducible
model limitation, not for an unexplained bad fit.

## Cancer natural history models: a worked target set

The recurring combination in cancer natural history modeling — a natural history model calibrated to registry
and screening-study data. Each target type with its likelihood:

| Target | Typical data | Likelihood |
|---|---|---|
| Age-specific cancer incidence | Cases and person-years by age group | Poisson, or negative binomial if overdispersed; rate x person-years offset |
| Stage distribution at diagnosis | Counts by stage, possibly by age | Multinomial per stratum; assert probabilities sum to 1 |
| Prevalence of precursor / precancerous lesions | Positives out of examined, by age | Binomial; beta-binomial if pooled across studies |
| Prevalence of preclinical (screen-detectable) cancer | Screen detections out of screened | Binomial; note this is prevalence at a screening round, not incidence |
| Cancer-specific mortality | Deaths and person-years | Poisson or negative binomial |

Points specific to this setup:

- **Stage distribution and incidence are not independent.** They come from the same
  registry cases. Summing their log-likelihoods double-counts those individuals. This is
  standard practice and usually tolerated, but it should be stated, and it is a reason to
  be skeptical of very narrow posteriors.
- **Preclinical prevalence is the target that identifies sojourn time**, and it is usually
  the scarcest. Expect the collinearity analysis to show that dropping it makes the
  preclinical duration parameters non-identifiable — that is the single most valuable
  thing the pre-calibration analysis will tell you about this model class.
- **Incidence targets span orders of magnitude across age.** Under a Poisson or negative
  binomial this is handled correctly and automatically, because the variance scales with
  the mean. Under a normal likelihood with a constant SE it is not, and the oldest age
  groups dominate the fit. This alone is a sufficient reason to prefer the count
  likelihood.
- **Weighting across target types** is unavoidable here: 18 age-specific incidence points
  against 4 stage categories against 2 prevalence points. State the implicit weighting and
  decide deliberately (`SKILL.md` Step 2).

## References

Kennedy MC, O'Hagan A. Bayesian calibration of computer models. *J R Stat Soc B.*
2001;63(3):425–464.

Brynjarsdóttir J, O'Hagan A. Learning about physical parameters: the importance of model
discrepancy. *Inverse Problems.* 2014;30(11):114007.
