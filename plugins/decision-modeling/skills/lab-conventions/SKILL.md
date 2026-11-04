---
name: lab-conventions
description: "FerAl Lab coding conventions for decision-analytic models in R, based on the DARTH coding framework. Use whenever writing, naming, organizing or refactoring R code for a decision model, cost-effectiveness analysis, microsimulation, Markov or state-transition model: naming variables, files, functions or test files; choosing prefixes such as v_, m_, a_, l_, df_, p_, r_, hr_; splitting an analysis into the five DARTH components (inputs, model, calibration, validation, analysis); writing the model as a function; choosing supporting packages (roxygen2, testthat, dampack, bookdown, pkgdown); or starting a new model repository. Also use to answer 'what does this prefix mean' or 'is this name right'. Do NOT use to review existing code for bugs or speed (use model-code-review) or to set up a calibration (use dm-calibration), though both follow these conventions."
---

# FerAl Lab coding conventions

These conventions follow the DARTH coding framework as closely as possible, taken from the
published paper rather than from memory. Where the lab deliberately departs from the paper,
this skill says so.

- Alarid-Escudero F, Krijkamp E, Pechlivanoglou P, Jalal H, Kao SY, Yang A, Enns EA.
  "A need for change! A coding framework for improving transparency in decision modeling."
  *PharmacoEconomics* 2019;37(11):1329-1339. doi:10.1007/s40273-019-00837-x
- Reference implementation: https://github.com/DARTH-git/darthpack

If a change you are about to propose (an optimization, a refactor, a rename) would break one
of these conventions, say so explicitly and let the author decide. Do not silently rename
objects or restructure a project.

## Project structure

Do not restate or improvise a folder layout. The lab's standard layout, starter files and
`.gitignore` live in one place, the template repository
[FerAlLab-git/project-template](https://github.com/FerAlLab-git/project-template). Start new
model repositories from it ("Use this template" on GitHub), and point to it when asked where
a file belongs. That layout is for decision-analysis model repositories, not for lab
operations or other kinds of projects.

## Code priorities

R is the primary language. Models are large and long-running, so when writing or improving
code, optimize for computational efficiency first: vectorize, preallocate, never grow objects
inside loops, and prefer matrix operations over `apply` where it matters. Rcpp is in scope
once profiling justifies it. If the efficient version is hard to follow, say so and offer the
simpler version alongside it.

Teaching code is the exception: there, clarity comes first (see `model-code-review`, Step 0b).

## Naming

### Variables: `<x>_<y>_<var_name>`

`x` is the data type prefix, `y` is the variable type prefix, and `var_name` is a brief
descriptor in lowercase, separated by underscores. Either prefix is omitted when it does not
apply.

| Prefix | Data type |
|---|---|
| (none) | scalar |
| `v` | vector |
| `m` | matrix |
| `a` | array |
| `df` | data frame |
| `dtb` | data table |
| `l` | list |

| Prefix | Variable type |
|---|---|
| `n` | number |
| `p` | probability |
| `r` | rate |
| `u` | utility |
| `c` | cost |
| `hr` | hazard ratio |
| `rr` | relative risk |
| `ly` | life years |
| `q` | QALYs |
| `se` | standard error |
| `or` | odds ratio (lab addition, not in the paper) |
| `d` | discount rate (lab addition, not in the paper) |

Examples from the paper: `n_samp`, `hr_S1D`, `v_r_mort_by_age`, `a_M`, `l_params_all`,
`df_out_ce`.

**The most common mistake: `n_` is not a scalar marker.** Scalars carry no data type prefix.
`n` is a *variable type* meaning "number", so `n_samp` is a scalar number of samples, and a
vector of counts by age is `v_n_cases_by_age`. Do not describe `n_` as "scalar" in code,
comments or reviews.

**A prefix must tell the truth.** `m_x[i, ]` drops to a vector unless `drop = FALSE`; an
object named `m_` must still be a matrix after subsetting. A lying prefix is worse than none.

### Legacy forms

Recognize these in existing code; do not use them in new code.

- **`lst_` for lists** appears in existing DARTH course materials and in `imisr`. Read it as
  `l_`. Do not rename it in existing code unless that code is being migrated as a whole.
- **Dot notation** (`m.P`, `p.HS1`, `v.dwc`) predates the underscore convention. In new code,
  give the underscore equivalent. Do not convert a file that is consistently dot-notated
  unless the whole file is being migrated. A file mixing both styles is a real problem,
  because the prefix stops carrying information.

### Files

`<component number>_<description>.<ext>` within their directory, for example
`analysis/01_model_inputs.R` and `R/02_simulation_model_functions.R`. The component numbers
are the five components below.

### Functions

`<action>_<description>`, starting with a verb, for example `generate_init_params()` and
`generate_psa_params()`.

### Test files

`test-<source file name>.R`, the form `usethis::use_test()` creates, for example
`tests/testthat/test-02_simulation_model_functions.R`. The paper suggests `test_`; the lab
uses the hyphen to match `usethis` and `imisr`.

## The five components

A decision analysis is split into:

1. Model inputs
2. Decision model implementation
3. Model calibration
4. Model validation
5. Analysis, with subcomponents 5a probabilistic analysis, 5b scenario and deterministic
   sensitivity analysis, and 5c value of information

Components 1 to 4 are the model development phase. The probabilistic analysis (5a) is the
primary analysis, not a deterministic base case.

**Inputs** are categorized as external, estimated, or calibrated. Calibrated parameters get
valid placeholder values in component 1 and are replaced once component 3 has run.

**The decision model is a function** mapping inputs to outputs, general enough to serve
calibration, validation and every analysis without modification. Interventions that open
different pathways are passed as an explicit categorical variable rather than by changing
parameter values. Do not mix the two approaches in one model.

## Supporting tools

- **roxygen2** for function and dataset documentation.
- **testthat** for unit tests, written alongside every new function and whenever a bug is
  found, one test file per source file (named as above).
- **dampack** for CEA output, PSA handling, CEACs, the CEAF and expected loss curves.
- **R Markdown with bookdown** for reports; **pkgdown** for the project website.
- Repositories are **structured as R packages**, so the model loads with
  `devtools::load_all(".")`, with an RStudio `.Rproj` file at the root. The template
  repository already has this structure.
