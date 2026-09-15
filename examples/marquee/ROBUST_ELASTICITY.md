# Robust multi-load elasticity marquee

`fs-marquee-elasticity-robust` runs the existing level-set/CutFEM topology
optimizer transactionally, then independently re-solves every returned candidate
under an explicit set of structural load scenarios before publication.

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
top,0.75,1.0,1,0,0.2
```

`edge` is `bottom`, `right`, or `top`; the left edge is the canonical clamp and
cannot also carry a robust external traction. `start`/`end` are normalized edge
coordinates in `[0,1]`. Tractions are signed normalized `(fx, fy)` values.
Weights are finite and nonnegative, with at least one positive case.

`AGGREGATE` is `sum` for the unnormalized weighted sum of independently solved
compliances, or `worst` for the largest weighted scenario compliance. Loads are
never summed before the elasticity solve, so equal and opposite scenarios cannot
cancel.

Candidate generation deliberately remains the existing canonical downward
mid-edge cantilever descent. Each bounded move-size candidate is subsequently
replayed under every CSV case; only its robust replay and cut-quadrature area
decide publication. This is a bounded robust selection layer, not yet a
simultaneous multi-load shape-gradient optimization.

Outputs are create-only:

* `summary.json` — baseline and accepted robust metrics and explicit claim limits.
* `candidates.jsonl` — every bounded candidate with all case compliances/refusals.
* `level-set.csv` — accepted geometry only.
* `trajectory.jsonl` — nominal-driver evolution evidence for the accepted case.

The normalized workflow makes no 3-D, physical-validation, KKT, global-optimum,
or simultaneous multi-load-gradient claim.
