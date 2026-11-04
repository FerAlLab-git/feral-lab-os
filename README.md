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

You need Claude Code and `git`. In a Claude Code session in a terminal:

```
/plugin marketplace add FerAlLab-git/feral-lab-os
/plugin install decision-modeling@feral-lab-os
```

**In Positron, VS Code or the desktop app**, the chat panel does not run `/plugin`.
Install from Terminal instead:

```bash
claude plugin marketplace add FerAlLab-git/feral-lab-os
claude plugin install decision-modeling@feral-lab-os
```

The same two commands work from any shell. The plugin is installed for your user, so it is
then available in every Claude Code surface on that machine; start a new session to load
it.

Check that it loaded:

```bash
claude plugin list                        # decision-modeling@feral-lab-os, Status: enabled
claude plugin details decision-modeling   # Component inventory lists the three skills
```

If your machine has no GitHub SSH key, set `CLAUDE_CODE_PLUGIN_PREFER_HTTPS=1` in your
environment so Claude Code clones over HTTPS. The repository is public, so no GitHub
account is needed.

## Releasing

Maintainers follow these steps for every change to a skill.

1. **Make the change on a branch** and open a pull request. `main` is protected.
2. **Bump the version** in `plugins/decision-modeling/.claude-plugin/plugin.json` using semantic versioning:
   - patch (0.1.0 to 0.1.1): wording fixes, clarifications
   - minor (0.1.0 to 0.2.0): new guidance, a new skill, changed behavior
   - major (0.x to 1.0.0): a renamed or removed skill
3. **Validate** before merging:
   ```bash
   claude plugin validate .
   ```
4. **Merge, then tag the release** from `main`:
   ```bash
   git checkout main && git pull
   git tag decision-modeling--vX.Y.Z
   git push origin decision-modeling--vX.Y.Z
   ```
   Add a short GitHub release note saying what changed and whether users need to do anything.
5. **Re-upload any skill you also keep on claude.ai.** Skills uploaded to claude.ai do not update from this repo. If you changed one, zip it and upload it again in claude.ai skill settings, replacing the old copy:
   ```bash
   cd plugins/decision-modeling/skills
   zip -r ~/Desktop/<skill-name>.zip <skill-name>
   ```
   An out-of-date copy on claude.ai will also sync into Claude Code and run alongside the plugin version.

### Getting updates (lab members)

In Terminal (not Positron or the desktop app):

```bash
claude plugin marketplace update feral-lab-os
claude plugin update decision-modeling@feral-lab-os
```

Then start a new Claude Code session.

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
