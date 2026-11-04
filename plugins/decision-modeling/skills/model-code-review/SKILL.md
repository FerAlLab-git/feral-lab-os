---
name: model-code-review
description: "Review R code and simulation model implementations for correctness, numerical soundness, computational efficiency, and reproducibility. Use whenever the user asks to review, critique, debug, speed up, or check someone's code, including their own, a trainee's, or code accompanying a manuscript. Trigger on 'review this code', 'why is this slow', 'check my model implementation', 'is this right', 'can you look at this script', and on any mention of microsimulation, Markov model, transition probabilities, calibration code, PSA loops, Rcpp, or a model repository. Prioritize correctness first and computational efficiency second, since decision models are often large and long-running. Do NOT use for reviewing manuscript prose. Do NOT use for setting up, writing, or designing a NEW calibration, choosing a calibration method, or diagnosing non-identifiability (use dm-calibration); this skill is for reviewing calibration code that already exists."
---

# Model code review

Review order is not negotiable: **correctness, then reproducibility, then speed.** A fast wrong model is worse than a slow wrong model, because nobody re-examines something that finishes quickly.

## Step 0: establish what this code is supposed to do

Before reading a line, get or infer: what the model computes, what the units of every returned quantity are, what the intended cycle length is, and whether a reference implementation or expected output exists. Reviewing simulation code without knowing the expected output is guessing.

State up front what you assumed. If a critical assumption is unavailable, ask one question rather than reviewing on a guess.

## Step 0b: establish the audience, because it reorders the review

Ask or infer whether this is didactic or research code, and state which you 
assumed.

**Didactic** (tutorial appendices, course and workshop scripts, anything a 
trainee is meant to read and learn from): clarity is a correctness requirement 
here, not a style preference. An explicit loop over cycles or over individuals 
is not a defect. Do not flag readable-but-slow constructions as problems, and 
do not run the Step 3 list as a defect inventory. Instead, close the review with 
one short section that shows the faster implementation beside the didactic one, 
states the expected speedup, and gives one line on why the slower version is the 
right teaching choice. Never propose replacing the didactic version with the 
fast one.

**Research** (the user's own models, CISNET and applied project code, anything 
long-running): efficiency is in scope at full weight, in its normal position 
after correctness and reproducibility.

Correctness and reproducibility are reviewed identically for both. A teaching example that converts a rate incorrectly teaches the error to everyone who reads it, which is worse than the same bug in a private model.

If the user has not said the code is for teaching or illustration, assume it is research code.

## Step 1: correctness

The bug classes below cause published errors in this field. Check every one explicitly rather than reading for general smell.

**Rate and probability conversion.** The single most common serious bug. A probability over one cycle is not a rate, and `p_new = p_old * (new_length / old_length)` is wrong for anything but very small probabilities. Check every rescaling of a probability to a different cycle length, every conversion from a published annual rate, and every place a hazard ratio is applied to a probability rather than a rate. The correct route is rate `r = -log(1 - p) / t`, rescale in rate space, convert back with `p = 1 - exp(-r * t)`.

**Transition matrices.** Rows sum to one, within tolerance, at every cycle and for every individual or subgroup, not just at initialization. Add the assertion if it is absent. Check that no probability is negative or above one after any covariate or hazard-ratio adjustment; clamping silently is a bug, not a fix.

**Time-dependence.** Distinguish state-dependence, time-in-state (tunnel), and calendar-time dependence. Check that age indexing into life tables is off by neither one cycle nor one year, which is the second most common bug here. Check the first and last cycle explicitly; loops are almost always right in the middle and wrong at the edges.

**Half-cycle correction and accumulation.** Applied consistently to costs and effects, or deliberately absent. Check whether discounting is applied at the right point in the accumulation, and whether the discount vector is indexed from zero or one consistently with the cycle indexing.

**Random number handling.** Is the seed set, and set in a place that makes the whole pipeline reproducible rather than just one function? Under parallelism, are streams handled properly (L'Ecuyer-CMRG or an equivalent), or do workers share or duplicate streams? Duplicated streams across workers silently destroy the independence the PSA assumes. Are common random numbers used across strategies where variance reduction is intended, and if so is the pairing actually preserved?

**Sampling and indexing.** `sample()` with a length-one first argument silently switches to `1:n`. `apply()` over a data frame coerces to a matrix and can turn numerics into characters. `seq_len(0)` versus `1:0`. Integer overflow on large population counts. `which.max` ties. Off-by-one in cumulative-probability inversion sampling: check the boundary condition explicitly.

**Numerical.** Sums of many small probabilities, `log(0)`, division by a denominator that can be zero when a state empties, and comparisons of floating-point values with `==`.

**Silent recycling.** Vector recycling in an arithmetic operation between mismatched lengths is R's most dangerous convenience. Flag every operation where the lengths are not obviously equal.

## Step 2: verify, do not just inspect

Reading code finds some bugs. Running it finds more. Where the environment allows, propose or write:

- **A null test.** Set the intervention effect to zero and assert the output equals the comparator arm exactly.
- **Conservation checks.** Total population is conserved across states each cycle; nobody leaves an absorbing state.
- **Extreme-value tests.** Perfect sensitivity, zero cost, 100 percent adherence; assert the expected bounds.
- **A regression test.** Pin current output for a fixed seed so future refactoring is detectable. This is the highest-value single artifact to add to a model repository that lacks one.
- **Seed stability.** Rerun with several seeds and report Monte Carlo error on the headline outputs. If the between-seed spread is comparable to the difference between strategies, say so loudly.

Write these as `testthat` tests (`tests/testthat/test-<source file>.R`) where the project already uses them, plain scripts otherwise.

## Step 3: efficiency

Only after correctness. And **profile before optimizing**: propose `profvis` or 
`Rprof` output rather than guessing at the hot path. Guessed optimizations in 
simulation code are usually wrong about where the time goes.

For didactic code, do not work through the list below as a set of findings. 
Collapse it into the single side-by-side comparison described in Step 0b: take 
the two or three items with the largest payoff, show the faster code, state the 
magnitude, and stop. The ordered list is for research code.

Ordered by typical payoff in this kind of model:

1. **Growing objects in a loop.** `rbind`, `c(x, new)`, and `df[i, ] <- ` inside a loop are quadratic. Preallocate a matrix of the final size.
2. **Data frames in hot loops.** Data frame subsetting is orders of magnitude slower than matrix subsetting. Use matrices with dimnames inside the simulation and convert to a data frame once, at the end.
3. **Vectorize across individuals, loop over time.** In a microsimulation the time loop generally cannot be vectorized, but everything inside it should operate on the whole population at once. A nested loop over individuals inside the time loop is the classic 100x slowdown.
4. **Sample once, not per iteration.** Draw all random numbers for a cycle in a single call. `runif(n)` once beats `runif(1)` n times by a large factor.
5. **Avoid repeated `apply` over rows.** Row-wise `apply` is a loop with overhead. Use `rowSums`, `colSums`, `matrixStats`, or a matrix multiply.
6. **Inversion sampling via `findInterval` or a matrix multiply** rather than per-individual `sample()`.
7. **Move the invariant out.** Anything not depending on the loop index gets computed once above the loop. In PSA code, check what is being recomputed for every parameter set that could be computed once.
8. **Parallelize the outer loop only** (PSA sets, calibration candidates), never the inner time loop. Watch the cost of copying large objects to workers; `future` with appropriate globals handling, or forking on Linux, beats naive socket clusters for big shared state.
9. **Rcpp only when profiling justifies it,** and only for a tight, well-tested inner function with a regression test already in place. Crossing the R/C++ boundary repeatedly in a loop can be slower than staying in R. Never introduce Rcpp into code that lacks tests.
10. **Memory.** Check for unnecessary copies of large matrices, and for storing full individual-level trajectories when summary statistics per cycle would do. Memory pressure often masquerades as a CPU problem.

State the expected speedup magnitude for each suggestion, and be honest when it is small. A change that risks introducing a bug for a 5 percent gain should be labeled as not worth doing.

## Step 4: reproducibility and structure

- Is the environment pinned (`renv`, a lockfile, a container)? Without it, the results are not reproducible in any meaningful sense.
- Are paths relative and project-rooted (`here`) rather than absolute machine paths?
- Are parameters in a data file rather than hard-coded in the analysis script? Hard-coded magic numbers scattered across scripts are the reason models cannot be audited.
- Does the repository run end to end from a clean clone, and is there a single entry point that says so?
- Are outputs written to a versioned location with the seed and parameter file recorded alongside them?
- Is there a `README` that states how to reproduce the paper's headline number specifically?
- Function length and single responsibility; anything over roughly 100 lines that both computes and plots should be split.
- Naming follows the DARTH convention as set out in the `lab-conventions` skill: 
`<x>_<y>_<var_name>`. Data type prefixes: `v_` vectors, `m_` matrices, `a_` 
arrays, `l_` lists, `df_` data frames, `dtb_` data tables, and none for 
scalars. Variable type prefixes include `n_` number, `p_` probability, `r_` 
rate, `hr_` hazard ratio, `or_` odds ratio, `u_` utility, `c_` cost, `d_` 
discount rate. `n_` means "number", not "scalar": `n_samp` is right, and calling 
`n_` a scalar marker is a mistake. Parameter lists as `l_params_all`, state 
names as `v_names_states`. `lst_` is a legacy list prefix in older course 
code and `imisr`; flag it `[MINOR]` in new code only.
- Dot notation (`m.P`, `p.HS1`, `v.dwc`) predates that convention and is legacy. 
In new code, flag it `[MINOR]` and give the underscore equivalent. Do not 
propose converting a file that is already consistently dot-notated unless the 
whole file is being migrated. A file mixing both conventions is `[MAJOR]`, 
because the prefix stops carrying information.
- A prefix that lies about the object's type is worse than no prefix. Check that 
`m_` objects are actually matrices after every subsetting operation, since 
`m_x[i, ]` drops to a vector by default.
- Comments explain *why* rather than restating the code.

## Output format

Anchor every finding so the recipient can find it with search. The `ANCHOR` line is code
copied verbatim from the file, long enough to be unique in it; file, function and line number
are helpful redundancy, but line numbers drift and the quote does not. Never paraphrase an
anchor. Emit findings in file order, and before returning, re-check every anchor against the
source and drop any that does not match exactly.

```
[BLOCKER] R/simulate.R, function `run_cycle`, around line 88
ANCHOR: p_die <- rate_annual / 12
ISSUE: annual rate divided to a monthly probability without conversion.
FIX:
  p_die <- 1 - exp(-rate_annual / 12)
IMPACT: overstates monthly mortality; the error grows with the rate and
        propagates to every downstream cost and effect.
```

Structure the review as: **Verdict** (does this code do what it claims, yes or no), **Correctness blockers**, **Correctness majors**, **Verification tests to add**, **Efficiency** (with profiling first), **Reproducibility**, **Style and structure** (short), **What is good and should be protected**, **What I could not evaluate** (name the files you did not see).

## Two things not to do

Do not rewrite the whole script and present it as the review; the recipient cannot diff a rewrite against their intent. Give targeted changes.

Do not suggest a refactor for elegance in code that has no tests. Add the regression test first, then refactor.
