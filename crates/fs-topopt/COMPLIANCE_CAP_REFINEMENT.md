# Compliance requirements across a refined design grid

```sh
cargo run --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 -- --refine-material 1000 1500 200 250000
```

The first two numbers are authored compliance caps for the independent Y and Z
body-load laws, in the example's dimensionless units, not calibrated ratings.
They are never inferred from a baseline or loosened after refinement. The last
two numbers bound accepted optimizer updates **per grid** and cumulative Krylov
iterations. One quadrature control separately bounds all geometry builds. The
original worst-load and `--minimum-material` routes remain unchanged.

## What the command does

It restores a cap-feasible coarse seed, minimizes dilated material while every
load is solved on eroded material, then checks the retained design on an
enriched octree. Every physical compliance is compared directly to its cap on
both grids, alongside the two smooth envelopes. Each cap-normalized change
must be at most 0.02; that explicit numerical tolerance is printed with the
other policy values. A favorable average cannot waive an individual cap.

A passing coarse/fine check ends the command without unnecessary refitting.
Otherwise the command uses cap-normalized absolute DWR indicators to mark up
to two cells at a target marking fraction of 0.5. These marks refine the SOURCE
background, not the already refined diagnostic background. The candidate
inherits raw densities and protected regions, rebuilds its filter, reintegrates
the same reference loads, restores feasibility, and minimizes material under
the unchanged requirements. One selective refit is permitted.

A fresh solve checks the candidate incumbent, and another enriched-grid probe
assesses that NEW material. Success requires this current two-grid cap and
change check, not just an optimizer iteration cap or a probe of the old design.
Unresolved checks, missing marks, restoration/optimizer limits, cancellation,
and exhausted budgets return nonzero after exporting any solved prefix. That
prefix may still be useful; its kind and grid are explicitly labeled. An
interrupted seed is never described as optimized material.

## Library boundaries

`estimate_compliance_cap_enrichment` assesses a raw design at explicit material
parameters. It reuses the existing eroded reference-load DWR path, rejects an
unchanged active space, retains individual values/caps and normalized changes,
and scales absolute cell indicators by each case's requirement. Fine physical
stiffness is inherited, not obtained by filtering raw controls on the fine grid.
Success and returned errors preserve the source nominal model and scales.

`refit_compliance_caps` solves a DIFFERENT, finer-grid design problem. It uses
the checked parent map for raw controls and physical solid/void labels, preserves
the filter radius and raw floor, then reconstructs the fine design map. The
returned candidate includes restoration, accepted history, a separate feasible
incumbent and a final recheck. The source is immutable and is not automatically
replaced. All solves share the supplied physical control; the caller owns the
geometry budget and unchanged SDF, material, support and load-law identity.

The diagnostic probe and fine-grid reoptimization must not be conflated. These
are discrete numerical results, not continuum load guarantees, global optima,
stress/buckling checks, manufactured length scales, or estimates of dilated
volume discretization error. No cross-grid optimization descent is promised.

## Focused native checks

```sh
cargo test --release -p fs-topopt --features cutfem-marquee \
  --lib sdf3::robust::envelope::limits
cargo test --release -p fs-topopt --features cutfem-marquee \
  --example worst_load_sdf3 refinement::tests
```

The editing environment lacks DSR, RCH and Cargo/rustc; native execution is
unverified. Independent 3-D body-fitted P1 FEM/SciPy checks exercise the
formulation, not production CutFEM or these Rust implementations.
