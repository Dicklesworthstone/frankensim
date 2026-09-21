# Authored stress-constrained cantilever studies

`fs-marquee-elasticity-stress --projected` accepts a new initial geometry and
mechanics without giving up same-area acceptance, sampled-stress admission,
cooperative wall budgets or exact accepted-state continuation:

```sh
cargo run --release -p fs-marquee --bin fs-marquee-elasticity-stress -- \
  --projected ./authored-study 1e12 3 30 0.6 8 0 300 \
  --initial-field ./initial.csv --youngs 2 --poisson 0.3 \
  --load 2 --load-band 0.125 --pause-after 1
```

The stress limit above is deliberately loose for a software exercise; it is not
an engineering allowable. Geometry is still a normalized unit square, mechanics
are 2-D plane strain, the left edge is clamped, and downward traction is applied
on the right edge over `0.5 - half_width ..= 0.5 + half_width`. `--load-band` is
that HALF width, not its total width. Loads and moduli share the normalized
model's units, not implicitly pascals or newtons. The canonical material law
refuses unsupported near-incompressible parameters.

Defaults preserve the original beam, load 1, half width 0.125, Young's modulus 1,
and Poisson ratio 0.3. Each optional flag may appear only once. The field CSV uses
`x_normalized,y_normalized,phi_normalized` followed by exactly `(2^level+1)^2`
finite nodal rows in x-fastest order. Existing elasticity exports are accepted
without resampling. Coordinates must match the requested lattice; missing,
extra, shifted or reordered rows refuse. Reads are bounded to 8 MiB.

The COMPLETE traction band must lie strictly inside the authored material.
Checking endpoints and midpoint is insufficient: the level-set enclosure checks
every crossed interval. Both left/right nodal traces are then fixed through
projection and all accepted updates. A disappeared load is never counted as a
compliance improvement. Other invalid geometry or infeasible area/stress still
refuses under the existing numerical gates.

`input-level-set.csv` retains the actual source field. `baseline-level-set.csv`
is the separately projected, solved feasible baseline against which progress is
measured. Import is a NEW study, not a resumed iteration. Summaries retain the
actual material and load, and the existing checkpoint retains them with exact
node bits, multiplier and ordinal. `--resume` needs no original input CSV and
refuses all problem overrides. See `STRESS_RESUME.md` for continuation and exit
semantics. These features do not add physical validation, stress-adjoint
optimality, a continuous-stress certificate or 3-D support.
