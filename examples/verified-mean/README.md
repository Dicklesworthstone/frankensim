# Worked example: a Verified mean temperature on a real material card

This project asks the default `frankensim` binary for the one QoI it reports
as **Verified**: `temperature-volume-mean-bound`, a numerical interval for the
volume-mean temperature of the `temperature-max` requirement's region. It is
the tracked reference project (`data/reference-project/cooling-reference.fsim`)
moved onto a 1 m cube at schema 17, with its own metadata and source
names, the cube's source hash, its part and region renamed (`cube`, `block`)
and one added output row:

```lisp
(outputs
  (qoi :name "temperature-max" :kind "scalar")
  (qoi :name "temperature-volume-mean-bound" :kind "report"))
```

The mesh is the tracked `examples/contact-pair/cold-body.stl` and the card pack
the tracked `data/reference-project/aa6061.fsmcdpk`; no bytes are copied. The
fan system and vent are the reference project's fixture flow network, and the
conduction boundary is its declared convection coefficient (h = 10 W/(m² K) to
293.15 K on every face), not one derived from the network. The block
dissipates 5 W.

```bash
WORK=$(mktemp -d)
cargo run -p fs-cli --bin frankensim -- --json validate \
  examples/verified-mean/verified-mean.fsim
cargo run -p fs-cli --bin frankensim -- --json import \
  examples/verified-mean/verified-mean.fsim \
  examples/contact-pair/cold-body.stl "$WORK/mean.db" --unit m --max-hole-edges 0
cargo run -p fs-cli --bin frankensim -- --json solve \
  examples/verified-mean/verified-mean.fsim "$WORK/mean.db" \
  --materials data/reference-project/aa6061.fsmcdpk
# <run> is the "run" field of the solve output.
cargo run -p fs-cli --bin frankensim -- --json report <run> "$WORK/mean.db"
cargo run -p fs-cli --bin frankensim -- --json package <run> "$WORK/mean.db"
```

`the_tracked_verified_mean_example_bounds_the_real_card_end_to_end` in
`crates/fs-cli/tests/bounded_native_mean.rs` runs import, solve, report and
package on the tracked files and checks the properties below, not the measured
floats. A build with `--no-default-features` refuses the row at conduction by
name (`cli-solve-mean-bound-feature`).

## What the bound says

MEASURED 2026-11-26 on `8b7ee8987a` plus the change that made
`thermal-verification` a default feature of fs-cli (bead `frankensim-k0org`),
default build. The conduction receipt's `volume_mean_bound`:

| Field | Value |
|-------|-------|
| `authority` | `Verified` |
| `value_k` | 293.2338999324541 K |
| `enclosure_k` | [293.22893386035275, 293.2388659993437] K |
| `error_upper_k` | 0.0049660721039686 K |
| `region`, `selected_cells` | `block`, 22 |
| `region_volume_m3` | [0.9999999999999969, 1.0000000000000042] |
| `material_validity.method` | `continuum-slab-comparison`, two-sided along x |
| `continuum_temperature_range_k` | [293.15, 293.4037425149701] K |
| `retained_material_domain_k` | [233.15, 398.15] K |

The report shows the row as 293.2339 K, color Verified, with no uncertainty
budget. The evidence package (authority `structural-integrity-only`) carries
the same interval as an Estimated claim (`retained-nominal-mean-enclosure`):
its standalone checker does not reverify the continuum producer. That checker
passes.

An energy balance checks the scale without the solver: in steady state the
5 W leave through the cube's 6 m² at h = 10 W/(m² K), so the surface-mean
temperature is 293.15 + 5/60 = 293.2333 K. Heating from inside puts the volume
mean slightly above that, as `value_k` is.

## Why a cube

The AA6061 card's conductivity claim is valid only on [233.15, 398.15] K. A
finite card is admitted only when the whole-continuum slab comparison places
the PDE's entire temperature range inside that domain, and the comparison
needs an exterior on bounding-box planes (`docs/FINITE_MATERIAL_MEAN_BOUNDS.md`).
The cube meets that. The reference project with the same output row on its
tetrahedron `plate.stl` (4 facets) refuses at conduction: "finite spans
require a slab barrier inside the original material temperature validity".

## What it does not claim

The interval bounds the volume mean of the nominal linear conduction PDE on the
published mesh, with the declared boundary data and constant conductivity. It
is not a bound on the maximum, on material or boundary uncertainty, or on the
as-built geometry. The `temperature-max` QoI (293.234878 K) stays Estimated and
the run's verdict stays indeterminate. The receipt's `scope` string states the
full scope.
