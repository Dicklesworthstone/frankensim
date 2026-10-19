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

The power is **absorbed in the represented domain**, not incident beam power.
Clipping the footprint redistributes that declared power; no optical absorption,
reflection, penetration, escaped beam, or calibrated Gaussian/Goldak law is
inferred. P1 support includes elements incident to active nodes. Conserving
watts and counting active nodes do not bound source-interpolation error.
Continuous travel still needs temporal/spatial resolution: even step doubling
cannot certify unsampled excursions. There is no moving geometry, melt flow,
recoil pressure, nucleation, element birth or manufacturing-fidelity claim.
