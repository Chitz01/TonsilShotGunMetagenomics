# Instructions for coding agents (Codex etc.) working in this repository

- Never run `sbatch`, `srun`, `submit.sh` (without `--dry-run`), downloads,
  `prefetch`, `fasterq-dump`, `curl`/`wget` of sequencing data, or database
  installs. The user runs jobs on LUNARC COSMOS themselves.
- Never delete, move or rename the user's existing files, papers or PDFs.
- All configuration lives in `config.sh`; scripts must `source` it and must not
  hard-code paths, accounts or e-mails. Unknown values stay `__SET_ME__`.
- Outputs go under `ANALYSIS_ROOT` only. `DATA_ROOT` is read-only.
- Explain intended changes before editing anything outside the requested scope.
- Check shell scripts with `bash -n` (and `shellcheck` if available).
