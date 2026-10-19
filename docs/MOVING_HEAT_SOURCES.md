# Moving absorbed-power sources and scanned phase change

```sh
cargo run -p fs-conduction --example scanned_melting -- 2 100 > scan.csv
cargo run -p fs-conduction --example scanned_melting -- 2 400 > higher-latent.csv
```

The example moves a volumetric heating footprint through a fixed 3D slab,
turns it off at one second, and evolves the original enthalpy solver to 1.2 s.
CSV records include actual node positions, source position, temperature,
specific enthalpy, liquid fraction and the original per-step energy balance.
All material parameters and geometry are synthetic examples, not calibration.

`bc::history::march::moving::MovingHeatSource` accepts three independently
dimensioned world-position histories, positive axis-aligned half widths, an
absorbed-power history in watts, a reference time, and `MovingHeatBudget`.
`at` samples right-continuous endpoint values; `before` uses exact left limits
at held-table switches. Both return a nodal `ScalarField` directly consumable
by existing conduction, enthalpy and radiative solvers.

At each node, `r² = sum_i [((x_i-center_i)/a_i)²]` and the footprint weight
is `max(0, 1-r²)²`. The returned P1 interpolant is normalized using the
actual tetrahedral nodal volumes. Its separately integrated watts must meet
the caller's relative power tolerance. A nonzero source with insufficient
active non-orphan nodes refuses; it never silently deposits zero or selects a
nearest node. Zero-power travel outside the domain is permitted.

`enthalpy::scanning::ScanningEnthalpyRun` attaches a bounded set of these
sources to an unadvanced `EnthalpyHistoryRun`. Their fields add to the original
background generation, including overlapping sources. Original phase charts,
reference masses, nonlinear conductivity, phase-conductivity multipliers and
fixed thermal contacts remain unchanged. `ScanRadiation` optionally retains
sourced emissivities and independent reservoirs through every solve and retry.
This interface does not attach temperature-dependent contact iteration; that
separate contact-history runner is not replaced or silently downgraded.

`advance_to` takes an explicit endpoint. `advance_event_step` stops at all
source-position, power, boundary, background and reservoir knots, under the
requested time-step cap. `advance_adaptive_step` applies the existing
step-doubling test in temperature, enthalpy and liquid fraction. Both accepted
half steps retain full source and physical reports and consume two lifetime
slots; coarse/rejected trials consume none. A final cancellation check precedes
publication of the complete pair. Failed calls leave accepted state intact.

## Resolve motion independently of thermal endpoint agreement

Before taking a physical step, attach an explicit geometric limit:

```rust
use fs_conduction::bc::history::march::moving::motion::ScanMotionLimit;
run.set_motion_limit(cx, ScanMotionLimit {
    max_normalized_travel: 0.25,
    max_chebyshev_coefficients: 256,
})?;
```

The limit uses the sum of each axis's speed divided by that source's ellipsoid
half width. Every source constrains the next event/adaptive interval, including
off-source travel. Ramps and linear tables use their local slopes; a Chebyshev
axis uses all coefficients and the derivative inequality `|T_n'| <= n²`, not a
chord between sampled positions. The example selects a 0.05 half-width allowance.

This matters even with phase-aware step doubling: a smooth returning path can
visit other locations while its start, midpoint and end coincide. Both candidate
thermal trajectories then deposit heat at the same wrong location and agree
exactly. A motion limit restricts the proposed interval independently of that
agreement. It is a conservative path-speed calculation in exact arithmetic,
not an interval-certified bound or a temperature/phase error guarantee.

`motion_limit()` reports the retained policy. Selection or replacement is only
permitted before the first accepted step; cancellation or invalid selection
leaves it unchanged. Explicit `advance_to` refuses an interval exceeding this
limit or crossing a path knot rather than silently changing the requested time.
It keeps its existing right-continuous endpoint rule. Event/adaptive methods
keep their exact left-limit rule at jumps. A motion-limited adaptive interval
below its configured minimum must still pass the normal error test; no small
interval is automatically accepted. With no policy, existing behavior is unchanged.

Focused regressions exercise anisotropic travel, returning-polynomial motion,
held jumps, domain clamping, cancellation, manual endpoint equivalence, multiple
sources and retained lifetime budgets. Native execution must be reported
separately from independent equation tests.

The power is **absorbed in the represented domain**, not incident beam power.
Clipping the footprint redistributes that declared power; no optical absorption,
reflection, penetration, escaped beam, or calibrated Gaussian/Goldak law is
inferred. P1 support includes elements incident to active nodes. Conserving
watts and counting active nodes do not bound source-interpolation error.
Continuous travel still needs temporal/spatial refinement: neither step doubling
nor a travel allowance certifies field accuracy. There is no moving geometry, melt flow,
recoil pressure, nucleation, element birth or manufacturing-fidelity claim.
