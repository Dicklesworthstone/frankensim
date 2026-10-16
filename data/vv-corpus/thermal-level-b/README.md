# Level-B thermal cross-code frozen references

Bead: `frankensim-extreal-program-f85xj.4.3`.

Each `<case-id>.case.json` deck plus the pinned external environment in
`tools/vvref/` (pyproject + `uv.lock`) exactly re-derives every row of
`thermal-level-b-references-v1.tsv`:

```bash
UV_PROJECT_ENVIRONMENT="${TMPDIR:-/tmp}/frankensim-vvref-venv" \
uv run --project tools/vvref python tools/vvref/solve_skfem.py \
  --cases data/vv-corpus/thermal-level-b/*.case.json \
  --out   data/vv-corpus/thermal-level-b/thermal-level-b-references-v1.tsv
```

(The venv must live outside the repo: `xtask`'s source inventory refuses
symlinks under `tools/`; see `tools/vvref/README.md`.)

Regeneration is byte-stable on the recorded interpreter/platform; the
deck runner refuses to freeze anything unless its known-answer
self-checks (Dirichlet–Dirichlet and Dirichlet–Robin slabs against
closed forms) pass first.

## What these references are

Independent-implementation checks: the external code (scikit-fem on
numpy/scipy, SuperLU linear solves — no code shared with the FrankenSim
workspace) solves the SAME discrete system fs-conduction assembles
(same Kuhn/Freudenthal mesh witnessed by a BLAKE3 mesh-identity hash,
P1 elements, element-mean `k(T)`, consistent source/Robin mass rules).
Agreement therefore checks assembly, boundary-condition, and solver
implementations against each other at tight envelopes — it does NOT
measure discretization error and it is NOT physical validation. The
corpus registers every row as an Estimated-colour `CrossCode` (Level-B)
reference: two codes agreeing is not truth.

## Change discipline

A mismatch beyond a declared envelope opens an investigation bead —
never a silent envelope widening (golden-bump discipline). Any change
to a deck, the runner, or the pinned environment that changes the TSV
is a corpus version bump with a recorded reason.

The external side of these references is dev-only tooling under
`tools/vvref/`; nothing here enters the workspace runtime dependency
graph (`xtask check-deps` is unaffected).
