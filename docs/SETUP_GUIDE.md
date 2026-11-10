# Setup guide: from zip to pinned GitHub repo

Works on Linux, macOS (Intel or Apple Silicon) and Windows via WSL2. On Windows, do everything inside the WSL Ubuntu terminal, not PowerShell; bwa and samtools do not run natively on Windows.

## 1. Install micromamba (once)

```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
# accept defaults, then open a new terminal
micromamba --version
```

## 2. Create the environment

```bash
cd ngs-vv-demo
# Linux / WSL: exact versions used for the verification report
micromamba create -y -n ngs-vv-demo -f envs/conda-linux-64.lock
# macOS (the lock file is linux-only):
micromamba create -y -f envs/environment.yml
micromamba activate ngs-vv-demo
```

Apple Silicon: if bioconda packages are not found, create the env as x86_64:
`CONDA_SUBDIR=osx-64 micromamba create -y -f envs/environment.yml`.

## 3. Run it locally

```bash
make test            # or: pytest --cov=vvdemo
make trace           # rebuild docs/traceability_matrix.md and the README table
make demo-change     # CR-001: regression suite catches the harmful change
```

Expected: `71 passed`, `release gate PASSED`, and `RESULT: regression suite CAUGHT the harmful change`.

## 4. Put it on GitHub

1. Create an empty public repository named `ngs-vv-demo` on github.com (no README, no licence; the repo has both).
2. Replace `YOUR_GITHUB_USERNAME` in `README.md` (badge line) with your username. Put your name in `LICENSE`.
3. Push:

```bash
git init -b main
git add .
git commit -m "ngs-vv-demo 1.0.0: verified mini variant-calling pipeline"
git remote add origin https://github.com/YOUR_GITHUB_USERNAME/ngs-vv-demo.git
git push -u origin main
```

GitHub asks for a password: use a personal access token (Settings → Developer settings → Tokens), or install the GitHub CLI and run `gh auth login`.

## 5. Watch CI go green

Open the **Actions** tab. The first run takes ~4–6 minutes (environment download); later runs use the cache. Both jobs should be green:
- *Verify pipeline*
- *Regression suite catches harmful change (CR-001)*

If a job is red, open it, read the failing step, fix, push again. Common causes: an action version GitHub has retired (update `actions/checkout` or `actions/upload-artifact` to the version GitHub suggests), or a conda solve problem (switch `environment-file` to `envs/conda-linux-64.lock`).

## 6. Take the README screenshot

1. Open the green run. The **Summary** page shows both green jobs and, below them, the traceability matrix the workflow writes to the job summary.
2. Screenshot the top of that page (job graph + start of the matrix). Crop to ~1400 px wide.
3. Save as `docs/img/ci-green.png`, overwriting the placeholder.
4. `git add docs/img/ci-green.png && git commit -m "CI screenshot" && git push`.

The README also has a live CI badge, so the screenshot is for first impression and the badge proves it is current.

## 7. Pin it on your profile

Your GitHub profile → **Customize your pins** → tick `ngs-vv-demo`. Add a one-line description in the repo's About box, e.g. "Variant-calling pipeline with requirement-traced tests, CI release gate and change-impact demo". Add topics: `bioinformatics`, `ngs`, `software-testing`, `pytest`, `snakemake`, `verification`.

## 8. Optional: GIAB run

See [giab_mode.md](giab_mode.md). Add your GIAB metrics to the verification report; real numbers are more convincing than synthetic 1.0s.
