# Nozzle interference against deformed deposition geometry

`fs_solid::bridge::contact::nozzle` checks a supplied finite nozzle envelope
against the actual accepted beam top surface. `Bridge::nozzle_clearance`
checks one affine motion; `nozzle_path_clearance` checks consecutive timed
positions, including vertical moves, holds and reversals. These are read-only
queries: they neither command a printer nor advance mechanics or heat.

## Geometry and results

Supply `Profile` with flat-tip radius, outer radius and nonnegative conical
flank slope. At radial coordinate r the underside rises by
`slope * max(r-flat_radius, 0)` above the tip, up to the outer radius.
The body above that lower envelope is treated as occupied. This is an explicit
approximation to a nozzle, not a measured machine default. Its axis stays on
the beam's transverse centreline; lateral offsets and rotations are not admitted.

A `Sweep` specifies start/end [world x, underside-tip z] positions in metres
and a required clearance margin. The entire finite footprint and motion are
checked continuously in the represented polynomial geometry, not by sampling
only the tool tip or the two motion endpoints. Each beam cell's world-x map
includes the accepted axial strain and top-fibre rotation. Its top height is
quadratic. On each nozzle flank, clipping the position/time square gives a
polygon on which the gap is quadratic in position and affine in time; edge
endpoints and stationary points determine its minimum.

The report returns signed minimum gap, source cell, reference/world position
and motion fraction. Negative gap means penetration; gap at or below the
requested margin sets `interferes`. `minimum=None` means this surface was not
overlapped, NOT that other objects or the whole scene are clear. The witness
locates the minimum gap, not necessarily the first instant of interference.

`PathPoint` adds absolute seconds to a tip position. `PathConfig` bounds both
segment count and TOTAL cell visits across the path. All times must be finite,
nonnegative and strictly increasing. `PathReport` retains the worst witness
with its interpolated time, the first interfering segment, and complete work
counts. It does not retain every field or stop at the first detected collision.
Late invalid inputs, exhausted budgets or cancellation return no partial result.

## Physical example

```sh
cargo test -p fs-solid bridge::contact::nozzle
cargo test -p fs-couple --example sagging_deposition
cargo run -p fs-couple --example sagging_deposition
```

The existing example cools and deforms real thermoviscoelastic beam states,
with and without a solved midpoint support. It then checks a deliberately low
planned nozzle path against both accepted geometries. That path clears the
sagged road but intersects the higher supported road. The existing contact-area
and reheating comparisons remain intact. Coefficients and nozzle geometry are
synthetic. The snapshot is frozen at 5 s throughout the labelled 5–5.02 s path.

Seven focused unit regressions cover flat and penetrated surfaces, clearance
margins, intervening swept collisions, flank strikes despite a clear tip,
interior warp, thermal axial motion, reversal, complete-path accounting and
refusals. The coupled-example regression also checks both clearance outcomes.
Native Rust tests and the example have NOT run in the authoring environment:
DSR, its documented fallback, RCH and Cargo were absent. Python/SciPy references
check independent equations, not Rust execution or physical calibration.

This is nominal floating-point TOP-surface geometry, not an interval clearance
certificate, side-wall/overhang check, scene-wide detector, or impact model.
It assumes the accepted geometry is frozen during the query. Changes in
material/shape during actual motion need an independently resolved timeline.
No contact force, damage, material removal, controller safety guarantee or
calibrated printer prediction follows from this query alone.
