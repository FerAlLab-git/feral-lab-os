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
claude plugin list
```

You should see `decision-modeling@feral-lab-os` with `Status: enabled`. Then start
`claude`, type `/`, and look for `dm-calibration`, `lab-conventions` and
`model-code-review`.

If your machine has no GitHub SSH key, set `CLAUDE_CODE_PLUGIN_PREFER_HTTPS=1` in your
environment so Claude Code clones over HTTPS. The repository is public, so no GitHub
account is needed.

## Update

A new release reaches you only when the plugin's `version` changes. To update, run these in
Terminal (not in the chat panel of Positron, VS Code or the desktop app):

```bash
claude plugin marketplace update feral-lab-os
claude plugin update decision-modeling@feral-lab-os
```

Then start a new Claude Code session.

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

### Get notified of new releases

Each version is published as a GitHub release, with notes taken from
[CHANGELOG.md](CHANGELOG.md). To get notified, open this repository on GitHub, click
**Watch**, choose **Custom**, tick **Releases**, and click **Apply**.

## Contributing

Anyone can propose a change to a skill. `main` is protected, so every change comes in as a
pull request.

If you are not a member of FerAlLab-git, fork the repository first and open the pull request from your fork.

1. **Create a branch:**
   ```bash
   git checkout -b <short-description>
   ```
2. **Edit the skill** under `plugins/decision-modeling/skills/<skill-name>/`.
3. **Add an entry to [CHANGELOG.md](CHANGELOG.md)** under a `## [Unreleased]` heading at the
   top (add the heading if it is not there). Say what changed and whether users need to do
   anything.
4. **Validate** from the repository root:
   ```bash
   claude plugin validate .
   ```
5. **Push the branch and open a pull request.**

Leave the version in `plugin.json` alone; a maintainer sets it before merging.

## Releasing (maintainers)

For each Contributing PR you accept:

1. **Bump the version on the pull request's branch, before merging.** Raise `version` in
   `plugins/decision-modeling/.claude-plugin/plugin.json` using semantic versioning:
   - patch (0.1.0 to 0.1.1): wording fixes, clarifications
   - minor (0.1.0 to 0.2.0): new guidance, a new skill, changed behavior
   - major (0.x to 1.0.0): a renamed or removed skill

   In the same commit, rename `## [Unreleased]` in `CHANGELOG.md` to
   `## [X.Y.Z] - YYYY-MM-DD`. Then merge the pull request.

   If you can't push to the contributor's branch (for example, a fork that does not allow
   edits by maintainers), merge the pull request as it is and make these two changes in a
   small release pull request instead.
2. **Tag the release** from `main` once the version bump is merged:
   ```bash
   git checkout main && git pull
   claude plugin tag --push plugins/decision-modeling
   ```
   `claude plugin tag` reads the version from `plugin.json`, checks that it agrees with the
   marketplace entry, tags the current commit `decision-modeling--vX.Y.Z`, and with
   `--push` pushes that tag to `origin`. It stops if the working tree has uncommitted
   changes or the tag already exists.
3. **Publish a GitHub release** for the new tag. On GitHub, open **Releases**, choose
   **Draft a new release**, select the tag, and paste that version's entry from
   `CHANGELOG.md` as the notes.
4. **Re-upload any skill you also keep on claude.ai.** Skills uploaded to claude.ai do not
   update from this repo. If you changed one, zip it and upload it again in claude.ai skill
   settings, replacing the old copy:
   ```bash
   cd plugins/decision-modeling/skills
   zip -r ~/Desktop/<skill-name>.zip <skill-name>
   ```
   An out-of-date copy on claude.ai will also sync into Claude Code and run alongside the
   plugin version.

## Credits

The worked examples (the Sick-Sicker and cancer relative survival models) come from the
DARTH workgroup's teaching materials and from HRP 263 at Stanford University. The coding
conventions follow Alarid-Escudero F, et al. "A need for change! A coding framework for
improving transparency in decision modeling." *PharmacoEconomics* 2019;37(11):1329-1339,
[doi:10.1007/s40273-019-00837-x](https://doi.org/10.1007/s40273-019-00837-x).

## License

[MIT](LICENSE)
