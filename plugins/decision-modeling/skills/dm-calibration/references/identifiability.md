# Non-identifiability: a worked example

Read this when a calibration shows any of the Step 3 symptoms, or when explaining
non-identifiability to someone who has not met it before. The numbers below are real,
from the 3-state cancer relative survival (CRS) model used in DARTH teaching materials.

## The setup

A 3-state model: NED → Mets → Death, governed by two transition probabilities,
`p_Mets` and `p_DieMets`. One calibration target: 59 annual survival proportions.
Prior: uniform on [0.04, 0.16] for both. Calibrated with IMIS, 1,000 posterior draws.

## What the output looked like

Entirely unremarkable:

| Quantity | Value |
|---|---|
| Posterior mean | `p_Mets` 0.088, `p_DieMets` 0.084 |
| 95% CrI, `p_Mets` | 0.045 – 0.146 |
| 95% CrI, `p_DieMets` | 0.045 – 0.144 |
| Posterior correlation | −0.97 |

A reader would conclude the two probabilities are around 0.085, imprecisely estimated,
and strongly negatively correlated. Every part of that is wrong.

## What was actually happening

Time to death in this model is the time spent in NED plus the time spent in Mets. That
sum is symmetric in the two transition probabilities, so the predicted survival curve is
**identical** when they are swapped — not approximately, but to machine precision:

```r
a <- run_crs_markov(c(p_Mets = 0.048, p_DieMets = 0.124))$Surv
b <- run_crs_markov(c(p_Mets = 0.124, p_DieMets = 0.048))$Surv
max(abs(a - b))   # 2.2e-16
```

The likelihood therefore has two exactly equal modes, at (0.048, 0.124) and
(0.124, 0.048). The posterior is exactly bimodal, with about half the draws at 
each mode — 47.4% had `p_Mets` < `p_DieMets`.

Every summary above is an artifact of averaging across the two:

- The **mean** (0.088, 0.084) is the midpoint of the two modes, which is a valley. Its
  log-likelihood is **−53.6**, against **+156.0** at the MAP. It is not a compromise
  estimate; it is a parameter set that fits nothing.
- The **credible intervals** span both modes, so they describe the distance between the
  modes rather than the uncertainty within either.
- The **−0.97 correlation** is the anti-diagonal line joining the two clusters. It is not
  a posterior correlation.

Nelder-Mead showed the same thing from the other direction: 100 random starts all reached
log-likelihood 156.03, but landed on *two* distinct parameter sets depending on where they
started. Reading only the single best row hides that completely; reading the top ten makes
it obvious.

## The diagnostics that catch it

Cheapest first:

1. **Swap test**, before calibrating at all — the `all.equal` call above. Catches exact
   structural symmetry for free.
2. **Top-10 from a multi-start search**, not the top 1. Equal fit at different parameter
   sets is the signature.
3. **Posterior mean versus MAP fit.** A gap of more than ~10 log-likelihood units means
   the mean is not in a mode. This is the general diagnostic: it catches multimodality
   that has no tidy symmetry behind it, which the swap test would miss.
4. **Pairwise plot of the draws.** Two clouds, or a ridge, rather than one blob.

## Resolving it, in this case

An ordering constraint, `p_Mets < p_DieMets`, applied in **both** the prior sampler and
the prior density. Sorting each sampled row enforces it in the sampler; the density needs
`-Inf` outside the ordered region and a `log(2)` renormalization inside it, since the
constraint halves the support of two exchangeable parameters.

Result after the constraint:

| Quantity | Before | After |
|---|---|---|
| Posterior mean fit (log-lik) | −53.6 | 155.9 |
| MAP fit (log-lik) | 156.0 | 156.0 |
| 95% CrI, `p_Mets` | 0.045 – 0.146 | 0.045 – 0.054 |
| 95% CrI, `p_DieMets` | 0.045 – 0.144 | 0.101 – 0.147 |

The marginal likelihood barely moved — 151.08 to 151.11 — which is the check that the
renormalization was right. The constraint relabelled the posterior; it did not change how
well the model fits.

Note what the constraint does and does not buy. The intervals are now interpretable, but
the model still cannot tell you which physical quantity is the larger one. The constraint
is a labelling decision, and it has to be defended on external grounds — here, knowing
that progression is slower than death after metastasis. Without such an argument, the
honest move is option 4 from Step 3: report both modes.

The residual within-mode correlation of about −0.93 is real, and worth reporting: the two
probabilities trade off against each other in preserving the survival curve even within a
single mode.

## Why this matters beyond the estimates

The two modes are indistinguishable in survival and give **identical life expectancy** —
27.21 life-years either way, agreeing to 1e-14. But they distribute that time completely
differently between the states:

| | Cycles in NED | Cycles in Mets | Total life-years |
|---|---|---|---|
| Mode A (`p_Mets` = 0.048) | 19.80 | 7.41 | 27.21 |
| Mode B (`p_Mets` = 0.124) | 8.06 | 19.15 | 27.21 |

**2.58 times more time in the metastatic state**, with no difference in survival and none
in life expectancy. Metastatic disease carries the costs and the utility decrement, so the
two modes that the survival target cannot distinguish produce materially different total
costs and QALYs — and therefore potentially different ICERs and different decisions.

That is the whole argument in one table. A calibration that reports a single tidy
parameter estimate here is not making a small technical error; it is picking one of two
answers arbitrarily and hiding a decision-relevant difference. The posterior feeds the
PSA: propagated honestly, a bimodal posterior gives a bimodal distribution of incremental
costs and effects, which may straddle a decision threshold. Collapsed to a mean that fits
nothing, it gives a single confident and meaningless answer.

Non-identifiability is a decision problem, not just an estimation problem.

The canonical treatment is Alarid-Escudero F, MacLehose RF, Peralta Y, Kuntz KM, Enns EA.
Nonidentifiability in model calibration and implications for medical decision making.
*Med Decis Making.* 2018;38(7):810–821.
