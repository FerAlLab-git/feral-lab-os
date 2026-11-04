# Parallelizing a model run over parameter sets

The advice applies wherever a decision model is evaluated many times over different
parameter sets — a calibration likelihood, a probabilistic sensitivity analysis, a scenario
sweep, an EVPPI inner loop — not only to calibration.

Read it before parallelizing anything, and especially before assuming it will help.

## First: measure, because it often makes things slower

Parallelism trades communication overhead for compute. It wins only when a **single model
evaluation** is slow enough to dominate the cost of shipping work to a worker.

A measured example from the DARTH teaching CRS Markov model — 3 states, 60 cycles, evaluated on
1,000 parameter sets:

| | Elapsed |
|---|---|
| Serial | 0.077 s |
| Parallel, 15 workers, one task per parameter set | 0.296 s |

Nearly **four times slower**, on 15 cores. Cluster setup costs a few hundred milliseconds
and the model runs in microseconds, so every worker spent longer receiving work than
doing it.

Time the serial version first:

```r
system.time(v_llik <- calc_log_lik(m_param_samp))
```

Rough guidance: if a single model run is under ~10 ms, parallelism will almost certainly
lose. Between 10 ms and 100 ms it depends on the number of sets and the size of the data
shipped to workers. Above ~100 ms per run — microsimulation, a natural history model with
fine age and time steps, an ODE system solved numerically — it pays, and pays well.

When the model is fast and the job is still slow, the bottleneck is usually the *number* of
evaluations, not their cost — and that is almost always the cheaper thing to fix. In a
calibration, shrink the search space (`references/precise-workflow.md`): a
coverage-corrected prior can cut required iterations by more than parallelism would. In a
PSA, check convergence before adding draws; many analyses run ten thousand when two
thousand had already stabilised the decision.

## Choosing the backend by operating system

Forking is faster to start and shares memory copy-on-write, so workers inherit the whole
session without explicit export — but it does not exist on Windows, and it is unsafe
inside RStudio and alongside multithreaded BLAS or GPU libraries. Socket clusters work
everywhere and require every object to be exported explicitly.

Detect capability rather than the operating-system name — `parallelly::supportsMulticore()`
accounts for RStudio and for the `parallelly.fork.enable` option, which a bare
`Sys.info()[["sysname"]] == "Darwin"` test does not:

```r
make_cluster_for_platform <- function(n_cores) {
  if (parallelly::supportsMulticore()) {
    parallel::makeForkCluster(n_cores)          # macOS, Linux, outside RStudio
  } else {
    parallel::makeCluster(n_cores, type = "PSOCK")  # Windows, or RStudio
  }
}
```

If `parallelly` is unavailable, the fallback test is
`.Platform$OS.type == "unix" && !identical(Sys.getenv("RSTUDIO"), "1")`.

## How many workers

`parallel::detectCores()` reports the physical machine and ignores cgroup limits, Slurm
allocations, and container quotas. On a cluster node it will happily report 128 cores when
the job was granted 4, and oversubscription makes the job slower, not faster. Prefer:

```r
n_cores <- parallelly::availableCores(omit = 1)   # respects Slurm, cgroups, env vars
n_cores <- max(1L, n_cores)                        # never request zero
```

`detectCores() - 1` also returns 0 on a single-core machine and `NA` when detection fails,
either of which errors out at cluster creation. If you must use it:
`max(1L, parallel::detectCores() - 1, na.rm = TRUE)`.

On a shared HPC node, take the scheduler's allocation, never the machine's core count.

## Chunk the work

One task per parameter set means one scheduling round trip per few microseconds of work.
Split into one chunk per worker:

```r
l_chunks <- split(seq_len(n_samp), cut(seq_len(n_samp), n_cores, labels = FALSE))

v_llik <- foreach::foreach(v_idx = l_chunks, .combine = c) %dopar% {
  calc_log_lik(m_params[v_idx, , drop = FALSE])
}
```

`drop = FALSE` matters: a single-row chunk silently becomes a vector without it, and the
likelihood function then treats one parameter set as several parameters.

## Create the cluster once

A Bayesian sampler calls the likelihood once per iteration. Creating and destroying a
cluster inside the likelihood pays the startup cost every time. Create it outside, pass it
in, and close it at the end:

```r
cl <- make_cluster_for_platform(n_cores)
on.exit(parallel::stopCluster(cl), add = TRUE)
doParallel::registerDoParallel(cl)
```

`on.exit` is not optional. Without it, any error leaves orphaned R processes holding
memory, and across a long interactive session they accumulate until the machine swaps.
Check with `ps aux | grep -c "[R] --slave"` if a session starts feeling slow.

## Exporting for socket clusters

Fork workers inherit everything. Socket workers start empty and need the model function,
the targets, and any helpers:

```r
parallel::clusterExport(cl, varlist = c("run_model", "l_targets", "v_target_names"))
parallel::clusterEvalQ(cl, library(matrixStats))
```

`.export = ls(globalenv())` in a `foreach` call is a blunt instrument — it works, but it
ships everything in the workspace to every worker on every call, which for a large
posterior or trace object is a real cost. Name what is needed.

The deeper fix is the model function contract: a function that takes all its inputs as
arguments and reads nothing from the global environment needs almost nothing exported, and
is also the only kind that can be unit-tested. See `SKILL.md` Step 1.

## Random numbers

This is the part that silently invalidates results. Workers started from the same state
produce **identical** random draws, so a stochastic model evaluated in parallel gets
duplicated noise across parameter sets — the independence a PSA assumes is destroyed, and
a calibration scores correlated noise as if it were signal. Nothing errors.

```r
parallel::clusterSetRNGStream(cl, iseed = n_seed)
```

This sets independent L'Ecuyer-CMRG streams per worker, and it is reproducible: rerunning
with the same `iseed` gives the same draws, and the draws differ across workers. Verified:

```r
cl <- parallel::makeCluster(2)
parallel::clusterSetRNGStream(cl, iseed = 42)
a <- parallel::parSapply(cl, 1:4, function(i) runif(1))
parallel::clusterSetRNGStream(cl, iseed = 42)
b <- parallel::parSapply(cl, 1:4, function(i) runif(1))
identical(a, b)              # TRUE  - reproducible
length(unique(a)) == 4       # TRUE  - streams are independent
```

Note that results become reproducible **for a fixed number of workers**. Changing
`n_cores` changes the stream assignment and therefore the answers. Record the worker count
alongside the seed, or the run is not reproducible.

For a deterministic model this does not apply — but set it anyway, because models acquire
stochastic components later and nobody revisits the parallel block.

## What to parallelize

**The outer loop over parameter sets.** Never the inner loop over cycles: those are
sequentially dependent, and the overhead per cycle is larger than the work.

For Bayesian samplers, check whether the sampler parallelizes for you before building your
own. `BayesianTools::createBayesianSetup(..., parallel = TRUE)` handles it internally, and
doing both gives nested parallelism that oversubscribes the machine badly.
`imisr::imis(..., parallel = TRUE)` does the same: it splits each batch of draws into one
block per worker through the future framework (pick the backend with `future::plan()`), so
a likelihood that parallelizes itself must be passed with `parallel = FALSE`. The archived
`IMIS` package had no built-in support, which is why older scripts parallelize inside the
likelihood.

Nested parallelism is worth checking for generally: if the model itself uses a
multithreaded BLAS, set `RhpcBLASctl::blas_set_num_threads(1)` on the workers so that
`n_cores` workers each spawning `n_cores` BLAS threads does not thrash.

## Reporting

Report elapsed time in units matched to the duration. `round(difftime(..., units = "hours"), 2)`
prints `0 hrs.` for anything under 18 seconds, which is the least useful possible
diagnostic. Use seconds and let the reader scale, or pick the unit from the magnitude.
