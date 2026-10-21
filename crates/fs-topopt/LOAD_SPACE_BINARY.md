# Binary material under a simultaneous-load requirement

```sh
cargo run --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 -- --load-space --binary-material 1500 200 250000 0.8
```

Arguments are the authored physical compliance cap, accepted gray optimizer
updates (default 200), total Krylov iterations (default 250000), and an optional
binary DILATED volume limit (default 1, admitted in [0,1]). Values shown are
illustrative requirements in this dimensionless example, not calibrated ratings.
The authored load set is the same Y and correlated Y/Z basis used by
`--load-space`: `f = F q`, `||q||_2 <= 1`. There is no probability model.

This command reuses the existing minimum-material restoration and optimization,
freshly rechecks the retained gray design, then selects a solid/void endpoint.
The original `--load-space --minimum-material` route is unchanged. Passing gray
material or passing each basis load separately does not authorize binary success.
All phases share ONE physical solve budget; the cap, basis, robustness and
spectral smoothing policy stay fixed. A failed optimization or gray recheck
returns before starting binary selection.

## Selection API and physics

`scan_binary_load_space_cap` admits the complete authored threshold schedule
before work. It filters/projects the source raw controls once, then uses the
existing binary evaluator to solve every eroded basis load at every threshold.
The original load-space spectrum owner recomputes all cross compliances, the
worst combined force/displacement and its physical work. There is no frozen
gray direction, force summation, skipped load, refiltered binary mask, or new
elasticity/eigensolver implementation.

`BinaryLoadSpaceScanOptions3` shares endpoint-policy validation with the finite
case scanner, not its feasibility decision. Each selected entry must pass the
full combined-load requirement, the smooth constraint and the separate dilated
material-volume limit. The selected entry has the least dilated volume among
completed feasible candidates. Exact ties keep the first authored entry. No
monotonicity or global optimum over thresholds/material layouts is claimed.

The command uses thresholds `[0.2, 0.35, 0.5, 0.65, 0.8]` and explicit binary
void stiffness `1e-3`. It prints every candidate's physical and spectral results,
then exports the selected three masks, actual eroded stiffness multipliers,
basis diagnostics, maximizing coefficients and combined force/displacement.
The prior gray export is labeled separately. Only current/best displacement
families and completed candidate summaries are retained by the scan.

A stopped scan preserves its completed feasible endpoint, when present, but is
not reported as complete. The command exports that endpoint before returning
an error. No qualifying binary candidate means `best = None` and a nonzero
command result, not a fallback to gray feasibility. All returned outcomes
preserve the caller's nominal operator state and charge spent physical work.

This is a fixed-grid ersatz-material result. Positive-stiffness void cells still
carry their authored loads/supports and cut integration. Passing can depend on
that declared void stiffness; there is no zero-void, extracted CAD, continuum,
connectivity, feature-size, stress/buckling or manufacturing certificate.

## Focused native checks

```sh
cargo test --release -p fs-topopt --features cutfem-marquee \
  --lib sdf3::robust::envelope::load_space::design::limits::binary
cargo test --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 load_space::material::binary::tests
```
