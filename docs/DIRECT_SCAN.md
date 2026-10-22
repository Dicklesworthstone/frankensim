# Deterministic global exploration of thermal source powers

`fs_dfo::direct::DirectState` implements the original DIRECT rectangle search.
The `simplex_scan_fit --direct` consumer uses it through the same value-only
P1 enthalpy simulation as the existing simplex/population modes. It needs no
objective derivatives or random initialization. The specimen, material charts,
moving footprint, normalization of absorbed watts, event-left sampling,
spatial fraction observations and physical energy checks are unchanged.

```sh
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example simplex_scan_fit -- --direct --initial-enthalpy 80
```

The explicitly synthetic cube starts at 80 J/kg (308 K), below its 310 K
solidus. The default targets are 0.1 and 0.3 volume-averaged nodal liquid mass
fraction after one and two seconds. The corresponding reference heat balance
is 30 W followed by 20 W. DIRECT's first evaluation is the midpoint of the
finite power box, not the all-off template. Its search covers this declared
box; it does not silently enlarge physically allowable heater powers.

```sh
# Explicit targets, maximum watts and lifetime physical-evaluation allowance:
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example simplex_scan_fit -- --direct --initial-enthalpy 80 0.1 0.3 50 3000

# Intentionally insufficient heater power. The output retains the mismatch,
# and the executable exits unsuccessfully rather than announcing recovery:
cargo run --manifest-path examples/scan_power_fit/Cargo.toml \
  --example simplex_scan_fit -- --direct --initial-enthalpy 80 0.1 0.3 1 65
```

One physical evaluation is reserved for a final replay. The replayed objective
must match the selected score bit-for-bit, and source export preserves the
fixed exact-zero shutdown. The executable initially admits eight search
attempts, then extends the same search and physical allowance without restarting
or refunding work. Its normal success criterion is actual replayed loss at or
below the declared target, not a small rectangle or a completed work allowance.

## Engine and resource behavior

The engine normalizes the box, selects rectangles satisfying the original
potential-optimality inequalities for some positive slope, samples both sides
of every longest axis, then trisects in order of the best axis samples. Equal
size/value ties are retained. A flat objective is not called converged: it may
be an unheated solid regime before the observed phase fraction changes.

`DirectConfig` explicitly bounds callback attempts, retained rectangle count,
rectangle-coordinate entries, and pair comparisons in one selection round.
`advance` has a second per-call attempt cap. Failed callbacks consume work and
retain their original candidate for retry; they are never converted into a
penalty. Successful callbacks are cached before cancellation is observed.
Partial samples do not publish a partial partition. Both evaluation and
partition allowances can be extended explicitly while retaining the same state.

The current pairwise selection is quadratic in the retained rectangle count.
This implementation targets bounded, low-dimensional expensive objectives; no
performance advantage over population search is claimed. The thermal CLI
admits at most `evaluations+4` rectangles, twice that many coordinate entries,
and the square of that rectangle cap in selection comparisons. These are
shape/work limits, not total-memory limits on the physical solver.

Stops distinguish target reached, evaluation exhaustion, a requested pause,
rectangle/storage exhaustion, selection-work exhaustion, and unresolvable
floating-point child centers. None except an actually met objective target is
reported as target recovery. No stop supplies a global lower-bound certificate,
KKT certificate, or proof that an unobserved better design does not exist.

`--direct` is mutually exclusive with `--population`, `--seed` and `--workers`.
DIRECT is serial and supports in-process continuation only. Existing disk
checkpoint/resume flags remain population-only and refuse in DIRECT mode.
No nonlinear constraints are silently applied: use the existing constrained
population or SQP study when physical feasibility requirements are needed.

## Verification scope

Seven engine regressions exercise selection, partition volume/nonoverlap,
atomic split publication, pause/evaluation boundaries, typed failure retry,
cancellation, finite resolution, Branin and a latent-balance objective. Two
consumer regressions exercise cold-start physical fitting, replay, bounded
unreachable targets, pause equivalence and command admission.

Independent Python DIRECT equations and SciPy's original DIRECT implementation
were compared on Branin. A separate NumPy P1/SciPy nonlinear endpoint reference
recovered approximately 29.99999059/20.00000941 W in 459 callbacks with loss
4.43e-13. The underpowered reference retained loss 5. Across these two reference
runs and their replays, 1,050 physical endpoints were solved. These are isolated
reference equations, not execution of the Rust implementation. Native Rust
compilation and tests were unavailable in the authoring environment; no compiler,
runtime, concurrency, benchmark or measured-speedup pass is claimed.

Algorithm: Jones, Perttunen and Stuckman (1993), *Lipschitzian optimization without
the Lipschitz constant*. Definition and partition discussion: Jones and Martins
(2021), *The DIRECT algorithm: 25 years Later*, doi:10.1007/s10898-020-00952-6.
