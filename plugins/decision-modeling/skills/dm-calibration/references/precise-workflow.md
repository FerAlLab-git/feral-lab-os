# PRE-CISE: pre-calibration coverage, identifiability, and sensitivity

Read this before running any expensive calibration. The three analyses below run in
minutes and routinely prevent days of wasted computation on a search space that cannot
work.

PRE-CISE is Gracia V, Goldhaber-Fiebert JD, Alarid-Escudero F. *PRE-CISE: A PRE-calibration
Coverage, Identifiability, and SEnsitivity analysis workflow to streamline model
calibration.* medRxiv 2026.02.27.26346591.
Reference implementation: https://github.com/vgraciaol/PRE-CISE

The key practical point: **do not hand-roll the sensitivity and collinearity machinery.**
It is already in the `FME` package (Flexible Modelling Environment), which PRE-CISE uses.
`FME::sensFun()` builds the elasticity matrix; `summary()` gives the norms; `FME::collin()`
gives the collinearity index over parameter subsets.

## The order, and why it is this order

1. **Coverage** — can the prior generate outputs that span the targets at all?
2. **Local sensitivity** — which parameters move which outputs, and in which direction?
   Use this to *fix* coverage failures.
3. **Collinearity** — can these targets separate these parameters, even in principle?

Coverage comes first because a prior that cannot produce the targets makes the other two
meaningless. Sensitivity comes second because it is the tool that repairs coverage.
Collinearity comes last because it needs the sensitivity matrix, and because it is only
worth asking once the search space is capable of fitting.

## Step 1: coverage analysis

Draw N parameter sets from the prior, run the model on each, and compare an interquantile
range of the predicted outputs (95% is typical) against the calibration targets. This is
not a fit check — it asks only whether the targets lie inside the envelope the prior can
generate.

```r
run_coverage <- function(m_params, f_model, v_times) {
  m_out <- t(apply(m_params, 1, function(v_p) {
    l_out <- tryCatch(f_model(v_p), error = function(e) NULL)
    if (is.null(l_out)) rep(NA_real_, length(v_times)) else l_out
  }))
  data.frame(time = v_times,
             lb   = matrixStats::colQuantiles(m_out, probs = 0.025, na.rm = TRUE),
             med  = matrixStats::colQuantiles(m_out, probs = 0.500, na.rm = TRUE),
             ub   = matrixStats::colQuantiles(m_out, probs = 0.975, na.rm = TRUE))
}
```

Read the result per target, not in aggregate:

- **Target above the envelope** — the model systematically underpredicts it across the
  whole prior. Shift bounds, do not merely widen them.
- **Target below the envelope** — systematic overprediction; same remedy, other direction.
- **Target inside but near an edge** — marginal coverage. Only a sliver of the prior can
  produce it, so the sampler will spend nearly all its effort in useless regions.
- **Envelope enormously wider than the targets** — the prior is too vague. Not wrong, but
  it makes calibration slow, and IMIS in particular will burn its early iterations.

Plot it. A coverage plot with the target error bars overlaid on the prior envelope
communicates the problem instantly and belongs in the supplement.

## Step 2: local sensitivity analysis, and using it to fix coverage

One-at-a-time perturbation around a baseline (the prior mean or midpoint), scaled into
**elasticities** — proportional change in output per proportional change in parameter — so
that parameters on different scales (a probability and a hazard ratio) are comparable:

    s_ij = (dy_i / dtheta_j) * (w_theta_j / w_y_i)

`FME::sensFun` computes this. The function passed to it must return a data frame whose
first column is the independent variable (time, age) and whose remaining columns are the
outputs. Leaving `varscale` and `parscale` at their defaults gives elasticities scaled by
the local output and parameter values, which is what PRE-CISE uses.

Verified against the Sick-Sicker model from the DARTH teaching materials:

```r
library(FME)
source("SickSicker_MarkovModel_Function.R")

f_model_fme <- function(v_params, v_times = c(10, 20, 30)) {
  l_out <- run_sick_sicker_markov(v_params)
  data.frame(time     = v_times,
             Surv     = l_out$Surv[v_times],
             Prev     = l_out$Prev[v_times],
             PropSick = l_out$PropSick)
}

v_params_mid <- c(p_S1S2 = 0.255, hr_S1 = 2.75, hr_S2 = 10)
sens <- sensFun(func = f_model_fme, parms = v_params_mid)
summary(sens)
```

which gives:

```
       value scale    L1    L2     Mean    Min    Max N
p_S1S2  0.26  0.26 0.536 0.691 -3.0e-01 -1.354 0.3639 9
hr_S1   2.75  2.75 0.017 0.022 -1.6e-02 -0.045 0.0029 9
hr_S2  10.00 10.00 0.254 0.307  1.9e-09 -0.407 0.6254 9
```

Read it as: `p_S1S2` dominates (L2 = 0.69), `hr_S2` matters (0.31), `hr_S1` is nearly
inert (0.02) at this baseline. The `Min`/`Max` columns carry the **sign**, which is what
tells you which direction to move a bound. A parameter whose elasticity changes sign
across outputs cannot be adjusted to help both at once — prioritize.

**Turning sensitivities into bound adjustments.** For output `i` and parameter `j`, the
proportional parameter change needed to close the gap is the proportional output gap
divided by the elasticity:

```r
# proportional change required in parameter j to move output i onto its target
n_pct_change <- ((y_target - y_model) / abs(y_model)) / s_ij
```

Compute this for every output the parameter influences, then combine across outputs with a
single stated criterion — minimum, maximum, or mean — depending on how conservative you
want to be. PRE-CISE's guidance, which is sound: **widen generously, narrow cautiously.**
A prior narrowed too far produces a truncated posterior, which is a much more damaging and
much less visible error than an inefficient search. Coverage iterations are cheap; the
calibration is not.

Re-run coverage after each adjustment and iterate until the envelope covers every target.

Local sensitivities are local — they describe the neighbourhood of the baseline point, not
the whole prior. After a large bound change, recompute at the new midpoint rather than
reusing the old matrix.

## Step 3: collinearity analysis

The collinearity index asks whether the columns of the (normalized) sensitivity matrix are
nearly linearly dependent. If they are, different parameter combinations produce
indistinguishable output changes and the parameters cannot be separately recovered — the
structural, pre-calibration version of what a bimodal posterior shows you afterwards.

    gamma = 1 / sqrt(min(eigenvalue(S_hat^T S_hat)))

```r
collin(sens)                              # every parameter subset, all targets
collin(sens, which = "Surv")              # one target only
collin(sens, which = c("Surv", "Prev"))   # a pair of targets
```

**Thresholds.** Roughly 15 indicates practical non-identifiability in the inverse-modelling
literature; PRE-CISE reports against 10 and 20 as well. Treat these as soft. An index of
14 is not meaningfully better than 16, but 5 versus 200 is decisive.

**Evaluate at more than one point — the index is local and the point matters enormously.**
The sensitivity matrix is computed at one parameter vector, so the collinearity index
describes identifiability *there*, not across the prior. On the CRS model, with survival as
the only target:

```
symmetric midpoint (0.086, 0.086)   collinearity = 94,906,266
at a mode          (0.048, 0.124)   collinearity =          5.7
```

Same model, same target, same two parameters — and the diagnostic swings between
"catastrophically non-identifiable" and "comfortably below threshold." At the symmetric
midpoint the two elasticities are exactly opposed, so the columns are perfectly dependent;
at a mode they are locally distinguishable and the index looks fine, even though the
global mirror-image mode is still there.

PRE-CISE's prescription to evaluate at the prior midpoint happens to catch this case,
because the midpoint of a symmetric prior *is* the symmetric point. Do not rely on that
coincidence. Evaluate `collin()` at the midpoint **and** at several draws from the prior,
and report the range. A single low value is not evidence of identifiability; a single high
value may be an artifact of landing on a symmetry. Pair it with the structural swap test in
`SKILL.md` Step 3, which is global and cheap.

**Verified output for the Sick-Sicker model**, which reproduces the paper's headline
result. Survival alone:

```
  p_S1S2 hr_S1 hr_S2 N collinearity
1      1     1     0 2          6.3
2      1     0     1 2        233.2
3      0     1     1 2          6.1
4      1     1     1 3       1136.7
```

Survival and prevalence together:

```
  p_S1S2 hr_S1 hr_S2 N collinearity
1      1     1     0 2          1.5
2      1     0     1 2          1.1
3      0     1     1 2          2.9
4      1     1     1 3          5.3
```

All three parameters are hopelessly non-identifiable from survival alone (1137), and the
pair (`p_S1S2`, `hr_S2`) is already unrecoverable (233). Adding prevalence drops every
index below 10. This is the paper's conclusion — recovering all three parameters requires
at least two targets — and it takes seconds to establish, before any calibration runs.

**What the analysis is actually for**, beyond a pass/fail:

- **How many parameters these targets can support.** Read down the `N` column. A single
  target may sustain two parameters and not three.
- **Which targets to collect or add.** Compare indices across target subsets; the subset
  that collapses the index is the data worth having.
- **Whether resolution matters.** PRE-CISE's COVID-19 case study found daily incidence
  targets gave indices below threshold where weekly aggregates did not. Worth testing
  before aggregating data for convenience.
- **When identifiability appears.** In the SIR testbed, collinearity fell below threshold
  only once the targets spanned the epidemic peak — during exponential growth only the net
  growth rate `r = beta - phi` is determined, so any (beta, phi) preserving it gives near
  identical trajectories. The index flags non-identifiability for the correct structural
  reason and resolves it exactly when the data first contain the information.

That last point is the one to internalize: a high index is not a defect to be suppressed.
It is a statement about what the data can support.

## What to do when the index is over threshold

In order of preference, matching `SKILL.md` Step 3:

1. **Add or refine a target** — the index tells you which one helps.
2. **Increase temporal resolution** of an existing target.
3. **Fix or tightly constrain the low-influence parameters** — the L2 norms identify them.
   `hr_S1` above, at L2 = 0.02, is a candidate to fix rather than calibrate.
4. **Reparameterize** to calibrate the identified combination instead of its components
   (the SIR case: calibrate `r = beta - phi` and one component, not both separately).
5. **Report the non-identifiability** and carry a posterior that reflects it.

## Fitting PRE-CISE into the script

Run it as a separate script that precedes the calibration and writes out the refined
bounds, rather than as a preamble inside the calibration script. The bounds are an input
to calibration and deserve their own provenance — the coverage plots, the sensitivity
table, and the collinearity table are supplement figures, and reviewers increasingly ask
how the prior ranges were chosen.

```
01_precise_coverage_sensitivity.R  -> refined v_lb, v_ub, saved to file + plots
02_calibrate.R                     -> reads those bounds
03_posterior_diagnostics.R
```
