# frankensim-vvref — external FEM deck runner (dev-only)

The external half of the Level-B thermal cross-code corpus
(`data/vv-corpus/thermal-level-b/`, bead
`frankensim-extreal-program-f85xj.4.3`).

This is deliberately OUTSIDE the FrankenSim Cargo workspace (the
`tools/oracle` precedent): scikit-fem/numpy/scipy are external
development oracles and must never enter the runtime dependency graph
checked by `xtask check-deps`.

* `solve_skfem.py` — solves committed case decks with scikit-fem and
  freezes probe temperatures into the reference TSV manifest. Runs
  known-answer self-checks before any freeze.
* `pyproject.toml` + `uv.lock` — the pinned environment
  (scikit-fem 12.0.2, numpy 2.5.3, scipy 1.18.1, blake3 1.0.9,
  CPython 3.12).

IMPORTANT: materialize the virtualenv OUTSIDE the repository. `xtask`'s
Rust-source inventory walks `tools/` directly and refuses any symlink,
and a Python venv always contains them, so an in-tree `tools/vvref/.venv`
breaks `check-all` and `generate-source-manifest` for everyone. Run:

```bash
UV_PROJECT_ENVIRONMENT="${TMPDIR:-/tmp}/frankensim-vvref-venv" \
  uv run --project tools/vvref python tools/vvref/solve_skfem.py ...
```

Why scikit-fem: no CalculiX/Elmer/OpenFOAM-class binary was available
on the reference host (dropped from Homebrew; no container runtime up),
and a pip-pinned pure-Python FEM stack gives an exactly re-derivable
lock with an independent assembly + SuperLU solve path. The corpus
schema and manifest format are code-agnostic: references from
CalculiX-class solvers can be added later under new case/corpus ids.
