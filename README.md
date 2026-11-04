# feral-lab-os

Claude Code skills for health decision science, from FerAl Lab. The repository is a
[Claude Code plugin marketplace](https://code.claude.com/docs/en/plugins/create-marketplace.md)
that currently offers one plugin, `decision-modeling`.

## Skills in `decision-modeling`

| Skill | Use it for |
|---|---|
| `dm-calibration` | Setting up, extending or debugging a model calibration in R: choosing the likelihood for each target (binomial, Poisson, negative binomial, multinomial, overdispersion, model discrepancy), running the PRE-CISE pre-calibration checks (coverage, sensitivity, collinearity) before spending compute, choosing a search method (random search, Nelder-Mead, IMIS, DEzs MCMC), diagnosing non-identifiability, and deciding whether parallelizing is worth it |
| `model-code-review` | Reviewing existing R and simulation-model code: correctness first (rate and probability conversion, transition matrices, time indexing, random numbers), then reproducibility, then computational efficiency |
| `lab-conventions` | Coding and naming conventions for decision models, based on the DARTH coding framework. For the project folder layout it points to [FerAlLab-git/project-template](https://github.com/FerAlLab-git/project-template) |

Claude uses each skill on its own when a request matches it. You can also call one directly,
for example `/decision-modeling:dm-calibration`. `dm-calibration` is for building a new
calibration and `model-code-review` is for reviewing code that already exists.

## Install

You need Claude Code and `git`. In a Claude Code session:

```
/plugin marketplace add FerAlLab-git/feral-lab-os
/plugin install decision-modeling@feral-lab-os
```

Or from your shell:

```bash
claude plugin marketplace add FerAlLab-git/feral-lab-os
claude plugin install decision-modeling@feral-lab-os
```

Check that it loaded:

```bash
claude plugin list                        # decision-modeling@feral-lab-os, Status: enabled
claude plugin details decision-modeling   # Component inventory lists the three skills
```

If your machine has no GitHub SSH key, set `CLAUDE_CODE_PLUGIN_PREFER_HTTPS=1` in your
environment so Claude Code clones over HTTPS. The repository is public, so no GitHub
account is needed.

## Update

A new release reaches you only when the plugin's `version` changes. To update by hand, run
this in a session:

```
/plugin marketplace update feral-lab-os
```

or this in your shell:

```bash
claude plugin update decision-modeling@feral-lab-os
```

Then run `/reload-plugins` in an open session, or start a new one.

### Turn on auto-update

Auto-update is off by default for marketplaces outside Anthropic's own. To turn it on, run
`/plugin`, open **Marketplaces**, select `feral-lab-os`, and choose **Enable auto-update**.

With auto-update on, Claude Code checks for a new version shortly after the first message of
a session (after a random delay of up to ten minutes) and installs it. You will see
`Plugin updated: decision-modeling · Run /reload-plugins to apply`; the new version loads on
your next launch either way.

Setting `DISABLE_AUTOUPDATER=1`, `DISABLE_UPDATES=1` or
`CLAUDE_CODE_DISABLE_NONESSENTIAL_TRAFFIC=1` turns plugin auto-update off, unless you also
set `FORCE_AUTOUPDATE_PLUGINS=1`.

## Get notified of new releases

Each version is published as a GitHub release, with notes taken from
[CHANGELOG.md](CHANGELOG.md). To get notified, open this repository on GitHub, click
**Watch**, choose **Custom**, tick **Releases**, and click **Apply**.

## Contributing

`main` is protected, so changes come in as pull requests. A pull request that changes a skill
should also:

1. Raise `version` in `plugins/decision-modeling/.claude-plugin/plugin.json`. Without this,
   installed copies don't update.
2. Add an entry to `CHANGELOG.md`.
3. Pass `claude plugin validate --strict .` from the repository root.

## Credits

The worked examples (the Sick-Sicker and cancer relative survival models) come from the
DARTH workgroup's teaching materials and from HRP 263 at Stanford University. The coding
conventions follow Alarid-Escudero F, et al. "A need for change! A coding framework for
improving transparency in decision modeling." *PharmacoEconomics* 2019;37(11):1329-1339,
[doi:10.1007/s40273-019-00837-x](https://doi.org/10.1007/s40273-019-00837-x).

## License

[MIT](LICENSE)
