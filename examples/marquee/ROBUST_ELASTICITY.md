# Robust multi-load elasticity marquee

`fs-marquee-elasticity-robust` runs the level-set/CutFEM topology optimizer with
independent structural load scenarios in the actual shape descent, then
independently re-solves every returned candidate under the same scenarios before
publication.

```bash
cargo run -p fs-marquee --bin fs-marquee-elasticity-robust -- \
  /tmp/frankensim-robust \
  examples/marquee/robust-elasticity-loads.csv \
  4 12 0.45 5 worst
```

The CSV format is one scenario per non-comment line:

```text
edge,start,end,fx,fy,weight
right,0.375,0.625,0,-1,0.4
right,0.375,0.625,0,1,0.4
right,0.20,0.35,1,0,0.2
```

`edge` is `bottom`, `right`, or `top`; the left edge is the canonical clamp and
cannot also carry a robust external traction. `start`/`end` are normalized edge
coordinates in `[0,1]`. Tractions are signed normalized `(fx, fy)` values.
Weights are finite and nonnegative, with at least one positive case. The initial
geometry must contain each declared loaded segment; after the first admitted
solve the optimizer retains every load pad as non-design material.

`AGGREGATE` is `sum` for the unnormalized weighted sum of independently solved
compliances, or `worst` for the largest weighted scenario compliance. Loads are
never summed before the elasticity solve, so equal and opposite scenarios cannot
cancel.

Candidate generation uses the same independent load cases directly in the
level-set shape field. Weighted-sum mode combines each case's self-adjoint strain
energy and topological derivative after separate equilibrium solves. Worst mode
uses the active weighted scenario as an explicit active-branch subgradient; an
exact tie selects the lowest input index. Each bounded candidate is then
independently replayed under every CSV case again, and only that final replay and
cut-quadrature area decide publication.

Outputs are create-only:

* `summary.json` — baseline and accepted robust metrics and explicit claim limits.
* `candidates.jsonl` — every bounded candidate with all case compliances/refusals.
* `level-set.csv` — accepted geometry only.
* `trajectory.jsonl` — simultaneous multi-load evolution evidence for the accepted case.

The normalized workflow makes no 3-D, physical-validation, KKT, global-optimum,
or global smooth-max-gradient claim.
