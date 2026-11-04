---
name: dm-calibration
description: "Model calibration in health decision science: estimating parameters no study reports by making the model reproduce data that does exist. Use in R to set a calibration up from scratch, extend one to new targets or parameters, choose the likelihood each target deserves (binomial; Poisson or negative binomial for counts and incidence rates with person-time; multinomial for stage or state distributions; beta-binomial, overdispersion, model discrepancy, target weighting), run pre-calibration coverage, sensitivity and collinearity checks before burning compute (the PRE-CISE workflow, via FME::sensFun and FME::collin), select or abandon a search algorithm (Latin hypercube random search, Nelder-Mead, IMIS, BayesianTools DEzs MCMC), judge whether parallelizing a likelihood is worth it, or debug output that looks wrong: collapsed effective sample size, one dominant importance weight, chains not converging, credible intervals that are huge or pinned to prior bounds, non-identifiable parameters. Targets are typically registry or SEER incidence, prevalence, survival, mortality, or stage at diagnosis; calibrated parameters are transition probabilities or natural-history rates that cannot be observed directly. Also for teaching labs comparing calibration methods. Do NOT use to review existing model code for bugs, speed or reproducibility (use model-code-review); nor for PSA, cost-effectiveness or value-of-information analysis on an already-calibrated model; nor for risk-prediction calibration plots, slopes or Hosmer-Lemeshow; nor for laboratory instrument calibration."
---

# Setting up a model calibration

Calibration estimates model parameters that cannot be observed directly by making the model reproduce quantities that were. The thing that makes it dangerous is that **a calibration almost never fails loudly.** A broken model throws an error or draws an obviously wrong trace. A broken calibration returns a tidy table of estimates with credible intervals, and nothing about it looks wrong.

So the order below is not about efficiency. It front-loads the checks that catch silent failure.

## Step 0: establish the three things before writing any code

Get these explicitly. Guessing any of them produces a calibration that runs and is wrong.

1. **Which parameters are being calibrated, and why they can't be observed.** If a parameter has a published estimate, calibrating it is usually a mistake — you are discarding evidence and inviting non-identifiability. Push back and ask what the direct estimate is.
2. **What each target is, and where its uncertainty comes from.** A target needs a value, a time or stratum it refers to, and a measure of precision. If the precision is unknown, say so rather than inventing a standard error; the choice silently sets the weight the target receives.
3. **Whether the parameters are identifiable from these targets.** Run Step 3 *before* the full calibration, not after.

State what you assumed. If item 1 or 3 is unavailable, ask one question rather than proceeding on a guess.

## Step 1: the model function contract

Everything downstream depends on the model being callable in one specific shape. Establish or refactor to this before anything else:

```r
run_<model>(v_params) -> named list, one element per calibration target
```

- **Input** is a single named vector. The function checks that every parameter it needs is present and stops if not. Without that check a misspelled name falls through to a same-named object in the global environment and the model runs with the wrong value — silently, because `with(as.list(v_params), ...)` resolves names in the enclosing scope.
- **Output** is a named list whose elements align **element-by-element with the rows of the corresponding target**, and only those rows. Return the time points the targets cover, not the model's full output — the cycles of a cohort trace, the solver time points of a compartmental model, or the summary across individuals of a microsimulation, as applicable.
- **Infeasible input raises an error.** Not a warning, not a silently clamped probability. The search loop catches the error and assigns the worst possible fit. Clamping hides the infeasible region and biases the estimates toward its boundary.
- **No side effects and no global reads** beyond the targets. A model that reads a global `n_t` cannot be parallelized or unit-tested.

**Verify the alignment numerically rather than by eye.** This is the most common silent calibration bug, and it is invisible in the output: the calibration simply fits the wrong time points to each other.

```r
stopifnot(identical(lengths(run_model(v_params_test))[v_target_names],
                    sapply(l_targets, nrow)[v_target_names]))
```

Length agreement is necessary but not sufficient — also confirm the *times* correspond. A trace row indexed from cycle 0 and a target indexed from year 1 will have matching lengths and be off by one cycle throughout.

**If the model is stochastic** (a microsimulation, or any model with sampling inside it), two further requirements apply, and ignoring them is why stochastic calibrations misbehave in ways that look like sampler problems:

- **Fix the random seed inside the model call**, so the likelihood is a deterministic function of the parameters. Otherwise the same parameter set scores differently on re-evaluation, Nelder-Mead chases noise and never converges, and MCMC acceptance ratios are corrupted — the chain will wander without ever settling.
- **Check that the Monte Carlo noise is small relative to the target standard errors.** Evaluate the likelihood several times at one parameter set with different seeds; if the spread is comparable to the differences the calibration needs to resolve, increase the simulated population before calibrating. No sampler can recover a signal smaller than its own noise floor.

## Step 2: build the likelihood in log space

**Choose the distribution from how the target was measured, not from habit.** The
normal-with-known-SE default is an approximation that is right often enough that people
stop noticing when it is not. Pick from the data-generating process:

| Target | Distribution |
|---|---|
| Proportion with numerator and denominator | Binomial |
| Proportion, overdispersed across sources | Beta-binomial |
| Counts, or incidence rates with person-time | Poisson |
| Counts with variance exceeding the mean | Negative binomial |
| Distribution over >2 categories (stage at diagnosis, cSTM state occupancy) | Multinomial |
| Only a mean and an SE available | Normal |

```r
# proportion with a numerator and denominator
sum(dbinom(x = df_target$n_event, size = df_target$n_total, prob = v_mu, log = TRUE))

# incidence: model gives a RATE, so multiply by person-time to get an expected count
sum(dpois(x = df_target$n_cases, lambda = v_rate_model * df_target$pyr, log = TRUE))

# overdispersed counts; large `size` means LESS overdispersion
sum(dnbinom(x = df_target$n_cases, mu = v_mu, size = n_k, log = TRUE))

# category distribution at one time point
dmultinom(x = v_n_stage, prob = v_p_stage, log = TRUE)
```

If the targets carry counts and the code uses `dnorm` on the proportion, say so and show
both — it is a real improvement, not a stylistic one. Two traps worth knowing by heart:
omitting the person-time offset for a rate target fits the wrong quantity silently, and
`dmultinom` **renormalizes `prob` without warning**, so assert `sum(v_p) == 1` yourself.

**`references/likelihoods.md` covers the rest**: the full selection table, how to fit
overdispersion as a calibrated parameter (and why it is weakly identified), beta-binomial
without an extra package, model discrepancy and its severe confounding caveat, and a
worked target set for cancer natural history models — age-specific incidence, stage
distribution, precursor lesion prevalence, preclinical cancer prevalence. Read it whenever
targets go beyond a simple proportion with an SE.

**Stay in log space everywhere.** `exp(log_lik)` overflows to `Inf` once the summed
log-likelihood exceeds about 709, which happens easily with a few dozen targets and small
standard errors. Samplers that require a non-log likelihood (the archived `IMIS` package
was one; its successor `imisr` takes the log-likelihood directly) need the log-likelihood
computed first and exponentiated once; check the maximum before trusting the run:

```r
max(v_llik)  # approaching 709 means subtract a constant or rescale
```

**Return `-Inf`, not `NA`, for parameter sets the model rejects.** `-Inf` propagates
correctly through sums and exponentiates to zero. `NA` silently poisons every downstream
summary.

**Two assumptions worth stating out loud in the write-up**, because reviewers ask and
because they are usually false:

- Summing log-likelihoods across time points assumes the targets are *independent*.
  Survival or prevalence measured on one cohort over time is strongly autocorrelated, so
  this overstates the information and produces credible intervals that are too narrow. It
  is standard practice, but it is an approximation; name it. The same applies across
  targets drawn from the same individuals — registry incidence and stage distribution
  share cases.
- Summing across targets with different numbers of data points assigns implicit weights.
  Thirty survival points against three proportion points means survival carries ten times
  the influence. Weight deliberately and justify it, rather than letting the count decide.

**On parallelizing the likelihood:** measure before assuming it helps. For a fast cohort
model it is routinely *slower* than serial — a measured 0.30 s against 0.077 s on 15
cores. See `references/parallel.md` for backend selection by platform, worker counts that
respect cluster allocations, chunking, cluster lifetime, and the RNG stream handling that
otherwise silently duplicates random draws across workers.

## Step 3: run the pre-calibration analysis before you calibrate

Non-identifiability produces the most confident wrong answers in this field, and almost
all of it is detectable *before* spending compute on a calibration. Run the PRE-CISE
workflow — coverage, local sensitivity, collinearity — as a separate script that precedes
calibration and outputs the refined prior bounds.

**PRE-CISE**: Gracia V, Goldhaber-Fiebert JD, Alarid-Escudero F. *PRE-CISE: A
PRE-calibration Coverage, Identifiability, and SEnsitivity analysis workflow to streamline
model calibration.* medRxiv 2026.02.27.26346591. Code: https://github.com/vgraciaol/PRE-CISE

1. **Coverage** — draw from the prior, run the model, and check the 95% interquantile
   range of predicted outputs actually spans each target. A prior that cannot generate the
   targets makes everything downstream meaningless, and this is the cheapest check there
   is.
2. **Local sensitivity** — one-at-a-time elasticities around the prior midpoint, giving a
   sensitivity matrix whose column norms rank parameter influence and whose signs say
   which direction to move a bound. This is the tool that *repairs* coverage failures.
   Widen generously, narrow cautiously: a truncated posterior is far more damaging and far
   less visible than an inefficient search.
3. **Collinearity** — the index `gamma = 1/sqrt(min(eigenvalue(S'S)))` over subsets of
   parameters and targets. Above roughly 15 means practical non-identifiability.

**Do not hand-roll any of this.** It is in the `FME` package, which PRE-CISE uses:

```r
library(FME)
# func must return a data frame: first column the independent variable, rest the outputs
sens <- sensFun(func = f_model_fme, parms = v_params_mid)
summary(sens)                              # L1 and L2 norms rank parameter influence
collin(sens)                               # index for every parameter subset
collin(sens, which = c("Surv", "Prev"))    # ... for a chosen set of targets
```

Comparing `collin()` across target subsets is the highest-value single output: it tells you
how many parameters your targets can support and which additional target would fix a
shortfall. On the Sick-Sicker model the three-parameter index is **1137** with survival
alone and **5.3** once prevalence is added — established in seconds, before any
calibration.

`references/precise-workflow.md` has the full workflow with verified code and output,
the sensitivity-to-bound-adjustment formula, and what each threshold means.

**Also run the structural symmetry test**, which collinearity can miss because it is a
local linear diagnostic and exact symmetry is a global property. If two parameters enter
so that predicted targets are invariant to swapping them, the likelihood has exactly equal
mirror-image modes:

```r
all.equal(run_model(c(a = 0.05, b = 0.12))$Surv,
          run_model(c(a = 0.12, b = 0.05))$Surv)   # TRUE => not identifiable
```

**After calibration, use the diagnostics that match the method you ran:**

| Symptom | Diagnostic | What it means |
|---|---|---|
| Multi-start Nelder-Mead reaches the same fit at *different* parameter sets | Sort the results, look at the top 10, not just the top 1 | Multiple modes; "the" best set is an artifact of the starting point |
| Posterior mean fits the targets far worse than the MAP | Compare `calc_log_lik(post_mean)` with `calc_log_lik(MAP)` | Multimodal posterior — the mean sits in a valley between modes |
| Posterior mass pressed against a prior bound | Fraction of draws within a few percent of `v_lb` or `v_ub` | The prior, not the data, is determining the estimate; widen and rerun |
| Ridge in the pairwise posterior, correlation near ±1 | `pairs.panels()` or `GGally::ggpairs()` on the draws | Parameters trade off; only their combination is identified |

The posterior-mean check is worth writing into the script permanently, because it is one line and it catches the case where someone reports a mean that fits nothing:

```r
if (calc_log_lik(v_post_mean) < calc_log_lik(v_map) - 10) {
  warning("Posterior mean fits far worse than the MAP: posterior is likely multimodal. ",
          "Do not report the mean, marginal CrIs or correlations as they stand.")
}
```

**What to do when a model is non-identifiable** — in order of preference:

1. **Add a target that distinguishes the parameters.** The real fix, and the collinearity analysis tells you which target does it.
2. **Increase the temporal resolution** of a target you already have. Daily rather than weekly aggregates changed identifiability in PRE-CISE's COVID-19 case study.
3. **Fix or tightly constrain the low-influence parameters**, which the L2 norms identify — a parameter with an elasticity norm near zero is not being estimated by these data anyway.
4. **Reparameterize** to calibrate the combination that *is* identified rather than its components (in an SIR during exponential growth, only `r = beta - phi` is determined).
5. **Impose an identifying constraint** (an ordering, say) in *both* the prior sampler and the prior density. Applying it to only one breaks the sampler's validity — the excluded region gets non-zero density and the second mode returns. Renormalize: an ordering constraint on two exchangeable parameters multiplies the density by 2.
6. **Report it honestly** — summarize each mode separately, and carry the full posterior (not a summary) into the downstream analysis.

Do not resolve non-identifiability by narrowing the prior until the problem disappears. That hides it rather than fixing it, and the estimate becomes a restatement of the prior.

See `references/identifiability.md` for a worked example with real numbers.

## Step 4: choose the search method deliberately

| Method | Gives you | Use when |
|---|---|---|
| Random search (LHS) | A picture of the GOF surface; crude best sets | Always, first. Cheap, and it reveals multimodality and ridges before you commit |
| Nelder-Mead, multi-start | Point estimate; local CI from the Hessian | You need a single best-fitting set, and ~100 starts to detect multiple modes |
| IMIS — `imisr::imis()` | Full joint posterior, and the marginal likelihood | Low-to-moderate dimension (say under ~10 parameters), fast model. The DARTH default, and well matched to teaching material. Use `imisr`, not the archived `IMIS` package |
| MCMC — `BayesianTools::runMCMC(sampler = "DEzs")` | Full joint posterior, with convergence diagnostics | Higher dimension, correlated or multimodal posteriors, or when IMIS degenerates (effective sample size collapsing, one dominant weight) |

Start with random search even when the destination is Bayesian. It costs minutes and tells
you whether the search space is sane before you spend hours.

**Choosing between IMIS and MCMC.** IMIS is importance sampling: it gives independent
draws, needs no burn-in or convergence assessment, and is excellent in low dimension. It
degrades as dimension grows — the importance weights concentrate, and the effective sample
size collapses even though the algorithm reports completion. Watch the `ESS`, `unique` and
`max_w` columns `imis()` prints at each stage, and the warning it gives when the stopping
rule is not met: a max weight that stays high or an ESS far below `n_resample` means the
posterior is being represented by a handful of points, and the answer is not trustworthy
no matter how tidy the output looks.

```r
library(imisr)  # successor of the archived IMIS: remotes::install_github("FerAlLab-git/imisr")
fit_imis <- imis(log_likelihood = calc_log_lik,  # matrix in, one value per row out
                 log_prior      = calc_log_prior,
                 sample_prior   = sample_prior,
                 n_resample     = 1000,
                 n_optim        = 0)  # 5-10 when the posterior may be multimodal (IMIS-opt)
fit_imis                  # stages, ESS, unique points, log marginal likelihood
plot(fit_imis)            # weight diagnostics by stage
m_post <- fit_imis$resample
```

`imis()` passes the likelihood a matrix with one parameter set per row, which is the shape
`calc_log_lik()` already takes. A likelihood written for a single parameter vector works
with `vectorized = FALSE`, and `catch_errors = TRUE` turns model errors into `-Inf` and
reports how many there were. Mapping from old `IMIS()` calls: `B` is `n_per_component`,
`B.re` is `n_resample`, `number_k` is `max_iter + 1`, `D` is `n_optim`, and the initial
sample of `10 * B` is `n_initial`.

Differential-evolution MCMC (`DEzs`) handles correlated and multimodal posteriors better
and scales to more parameters, at the cost of needing burn-in, thinning, and convergence
checks. It is the right default above roughly ten parameters, and the right escape hatch
when IMIS degenerates.

```r
library(BayesianTools)
n_iter <- 3e5
setup <- createBayesianSetup(
  likelihood = function(v_p) calc_log_lik(setNames(v_p, v_param_names)),
  prior      = createUniformPrior(lower = v_lb, upper = v_ub),
  names      = v_param_names,
  parallel   = FALSE          # let the sampler parallelize, or your likelihood, never both
)
l_settings <- list(iterations = n_iter, nrChains = 3, message = FALSE)
fit_mcmc <- runMCMC(bayesianSetup = setup, sampler = "DEzs", settings = l_settings)

gelmanDiagnostics(fit_mcmc)   # psrf should be < 1.05 for every parameter
tracePlot(fit_mcmc)
# `start` is indexed PER INTERNAL CHAIN -- see the warning below
m_post <- getSample(fit_mcmc, start = floor(n_iter / 3 / 3), thin = 10)
MAP(fit_mcmc)$parametersMAP
```

Three things that will bite:

- **The likelihood receives an *unnamed* vector.** Restore the names before passing it to
  the model, or the Step 1 parameter check fires — and if you skipped that check, the
  model reads whatever happens to be in scope instead.
- **`getSample(start = )` counts rows within one internal chain, not total iterations.**
  DEzs splits `iterations` across three internal chains, so with `iterations = 3e5` each
  chain holds 1e5 rows and the legal range for `start` is 1 to 1e5. Passing `start = 1e5`
  to discard "the first hundred thousand" leaves **nine draws**, and `start = 1e5 + 1`
  errors with a bare `subscript out of bounds`. Neither failure mentions burn-in, so this
  reads as the sampler having broken. Compute it as a fraction of `iterations / 3`.
- **`burnin` in `settings` is for the diagnostics, not for extraction.** The discarding
  that matters happens in `getSample(start = )`. Setting only `burnin` and then calling
  `getSample()` bare returns the full chain including warm-up.

`nrChains = 3` runs three independent replicates, which is what makes Gelman-Rubin
meaningful — the three *internal* DEzs chains share a Z-matrix and are not independent.
Convergence assessment is not optional here in the way it is simply absent for IMIS:
report the psrf.

**Method-specific traps:**

- **Nelder-Mead is unconstrained.** `optim()` will propose parameters outside the bounds and outside the feasible region. The objective must return `-Inf` for those. Separately, `optim()` refuses to *start* where the objective is not finite ("function cannot be evaluated at initial parameters") — so filter or resample infeasible starting points before the loop.
- **The Hessian from `optim()` is the Hessian of the objective itself, not the objective divided by `fnscale`.** With the log-likelihood as the objective and `fnscale = -1`, it is negative definite at the maximum and the observed information is `-hessian`, so `solve(-hessian)` is the right covariance. This looks like a sign error and gets "corrected" into a real one; leave a comment.
- **Hessian-based intervals are local to one mode** and asymptotic. Compare them against a Bayesian posterior before reporting them as uncertainty.
- **Bayesian samplers need the prior density and the prior sampler to agree.** Any constraint has to be in both.
- **Parallelizing the likelihood usually is not worth it.** If a single model run takes well under ~10 ms, cluster setup dominates and the parallel version is measurably slower. Time the serial version first, parallelize the outer loop over parameter sets and never the inner time loop, and chunk the sets across workers. Full treatment in `references/parallel.md`, including platform detection, worker counts that respect cluster allocations, and RNG streams.
- **Do not parallelize in two places at once.** If the sampler parallelizes internally (`createBayesianSetup(parallel = TRUE)`), the likelihood must not. Nested parallelism oversubscribes the machine and is usually slower than either alone.

## Step 5: the script skeleton

Follow the section order below. It matches the DARTH teaching scripts, so material written this way stays consistent with published DARTH tutorials and course materials.

Split it across numbered scripts rather than one long file. The refined prior bounds are
an *input* to calibration and deserve their own provenance — reviewers increasingly ask
how the ranges were chosen, and the coverage and collinearity plots are the answer.

```
01_precise.R    Coverage -> local sensitivity -> collinearity (Step 3)
                Writes: refined v_lb / v_ub, coverage plot, sensitivity table,
                        collinearity table across target subsets

02_calibrate.R  Header: model, parameters, targets, method, GOF measure
                Seed, declared once as n_seed and used everywhere
                Packages
                Load targets          -> inspect and plot before anything else
                Load model function   -> test-run; verify output/target alignment
                Read refined bounds from 01
                Specify: v_param_names, n_param, v_target_names, n_target
                Calibration functions -> sample_prior, calc_log_prior,
                                         calc_log_lik, calc_log_post
                Calibrate
                Export draws          -> for 03 and for downstream PSA/CEA

03_diagnostics.R  Convergence (MCMC) or ESS and max weight (IMIS)
                  Identifiability, bound proximity, posterior-mean-vs-MAP check
                  Explore posterior -> pairs, marginals, prior vs posterior
                  Posterior predictive -> model output vs targets, credible band
```

For a teaching script, collapse these into one file with the same section order — students
need to read it top to bottom. For research code, keep them separate: `02` is the
expensive step and you will rerun `03` many times against saved draws.

**Index parameters by name, never by position.** `v_lb[[s_par]]`, not `v_lb[i]`. Positional indexing across `v_param_names`, `v_lb`, and `v_ub` gives no error when it goes wrong — it silently assigns a hazard-ratio bound to a probability. Loop over the names:

```r
for (s_par in v_param_names) {
  m_param_samp[, s_par] <- qunif(m_lhs_unit[, match(s_par, v_param_names)],
                                 min = v_lb[[s_par]], max = v_ub[[s_par]])
}
```

**Pass only parameter columns to the model.** Result matrices usually carry fit and posterior-probability columns; `m_calib_res[1, v_param_names]`, not `m_calib_res[1, ]`.

Naming follows the DARTH convention as set out in the `lab-conventions` skill: `<x>_<y>_<var_name>`, with data type prefixes `v_` vector, `m_` matrix, `a_` array, `l_` list, `df_` data frame, `dtb_` data table (scalars take none), and variable type prefixes such as `n_` number, `p_` probability, `r_` rate, `hr_` hazard ratio, `or_` odds ratio. `n_` means "number", not "scalar". Watch for collisions: `lb`/`ub` are bound vectors *and* conventional column names for target confidence limits. Name the bounds `v_lb`/`v_ub` and leave the columns alone.

## Step 6: hand off to the downstream analysis

The reason to run a Bayesian calibration is to carry the joint posterior into the PSA — including its correlations, which a set of independent marginal distributions would destroy. Export the draws themselves, not a summary:

- Write the posterior draws to a file, with the seed and the target file recorded alongside.
- The PSA samples rows from that matrix. It does not re-sample from fitted marginals.
- If the analysis will use a point estimate, use the MAP and say why, and note what uncertainty was discarded.

## Output

When setting up a new calibration, produce: the script following Step 5, a short statement of the assumptions from Steps 0 and 2 (likelihood choice, independence, target weighting), and the identifiability finding from Step 3 stated plainly — including "no evidence of non-identifiability" when that is the answer.

When extending or fixing an existing calibration, give targeted changes rather than a rewritten script, and lead with anything from Steps 1–3 that is currently silent and wrong.

Flag honestly when a calibration should not be run at all: when a parameter has a perfectly good published estimate, when the targets cannot identify the parameters and no additional data is proposed, or when the targets are internally inconsistent and no parameter set can reproduce them.
