# Material memory in curved cymbal bending

`--shell-relaxation material.fssr` supplies proportional generalized-Maxwell
bending memory to `splash`, `splash-wav`, `splash-mic`, `hihat`, `hihat-wav`, and
`hihat-mic`. The existing shell reduction, material law, and joint impact time
owner produce its restoring force and dissipation. This changes physical shell
motion, not a microphone gain or an output decay envelope.

```sh
cargo run --release -p fs-couple --example percussion -- \
  splash 4096 --strike-position-m 0.06 0.01 --strike-speed-m-s 0.8 \
  --shell-relaxation crates/fs-couple/examples/percussion/illustrative-shell-relaxation.fssr \
  --analytic-newton --impact-substeps 8 511
```

The included spectrum is an illustrative numerical input, not identified
bronze, manufacturer data, or a calibrated cymbal. The focused regressions were
authored and reviewed in a session without a working Rust toolchain/dependency
environment; no native test result or full render is claimed.

## Explicit loss ownership and bounded input

```text
frankensim-shell-relaxation-v1
intrinsic_loss,replace
initial,relaxed
band_hz,0,5000
shell,single
branch,single,0.015,0.002
branch,single,0.005,0.02
```

`branch,target,ratio,tau_seconds` specifies excess flexural stiffness divided
by the actual equilibrium bending stiffness, with a positive relaxation time.
The ratio is dimensionless. It proportionally scales the existing per-facet
operator, including heterogeneous sections; no homogeneous Young modulus is
inferred. The material owner requires nonnegative ratios and positive times.
Each selected shell allows at most eight branches. Zero ratios create no state.
A declared shell with no branches is explicitly elastic.

For paired cymbals use `upper` and/or `lower` in place of `single`, declaring
each selected shell before its branches. Unselected shells retain their
supplied intrinsic damping. All selected shells share the file's initial state
and material band. No shell selection, band, or loss decision is defaulted.
Unknown, repeated singleton, malformed, nonfinite, and oversized inputs refuse;
the file limit is 64 KiB.

The required `intrinsic_loss,replace` record explicitly replaces the selected
shell's intrinsic modal-loss approximation. On a single splash, this removes
the host's documented damping ratio of 0.001 at construction. On a paired
hi-hat, each selected shell MUST already have zero intrinsic damping in the
supplied `.fshh` input's `damping_ratio,upper,lower` row; a conflicting nonzero
value refuses. Stand felt, contact loss, viscous mufflers, mallet/shaft loss,
pedal drag, squeeze-film/gas dissipation, and radiation reaction retain their
own physical meanings. They are not removed by material selection.

`initial,relaxed` equilibrates the new branches at the supplied initial shape.
`initial,unrelaxed` starts with zero viscous strain. Both are cold physical
initial conditions. The input does not replace a material during a performance.

## Complete bending energy and composition

The core projects the actual physical DKT bending form into the retained shell
basis, preserving all cross-mode terms. Factoring its nonzero energy subspace
provides strain-energy rows for the existing Maxwell storage and resistance
owner. An exact free translation has identically zero bending rows and gains
no spring, loss, or artificial memory. Nonzero singular subspaces refuse;
neither a pseudoinverse nor a rank cutoff makes them fit.

The original mechanical/force/source coordinates stay fixed. Material memories
follow existing felt/Kelvin and optional gap-gas states, before CLI radiation
loading and numerical preparation. Both selected shells attach together under
the existing 256-state material limit; complete spectra are never truncated.
Actual and requested frequencies must fit the declared band. The core also
checks instantaneous stiffness, the original Nyquist guard, and `dt/tau <= 0.25`.

CSV adds `shell_memory_energy_j`, already included in total storage, and
`shell_relaxation_power_w`, an endpoint power diagnostic. Accepted interval loss
remains in the ordinary shared energy/work ledger. Analytic preparation,
substeps, physical mallets, flexible shafts, and radiation feedback retain the
same state and rollback ownership.

This image relaxes linear bending only. Nonlinear membrane stretching and
numerical drilling remain elastic; it does not identify metallurgical loss,
viscous large-strain deformation, changing temperature, or measured decay.
The retained basis and fixed material band still need refinement for a claim
about impact transients or full-band sound.

Focused tests cover strict input and six command forms, free-translation
exclusion, changed contact-driven motion, selected/unselected loss ownership,
paired gas/history/source ordering, and accepted energy with rollback.

```sh
cargo test --release -p fs-couple --example percussion shell_relaxation -- --test-threads=1
cargo test --release -p fs-couple --example percussion hihat::relaxation_tests -- --test-threads=1
```
