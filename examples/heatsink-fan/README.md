# Worked example: a finned heatsink with a declared fan

This is the deepest walkthrough the product honestly supports today. It
declares the full cooling contract — fan system, vent area, airflow-leakage
bypass — so the solve prefix executes all four durable producer stages,
including the interval-certified flow-network operating point, and then
refuses at the conduction stage naming the bead that owns that gap. The
refusal is the product working: nothing downstream of a typed gap is
invented.

```bash
cargo run -p fs-cli --bin frankensim -- --json validate examples/heatsink-fan/heatsink-fan.fsim
```

Expected: exit 0, `"status":"ok"`, `"finding_count":0`. The G0 battery
(`crates/fs-cli/tests/cli.rs`,
`g0_the_heatsink_fan_example_reaches_the_flow_network_operating_point`)
executes validate, import, and solve continuously against exactly these
bytes, so this walkthrough cannot silently rot.

## The geometry

`heatsink.stl` (ASCII, 60 facets) is five axis-aligned boxes in meters:

| Part | Extents (x, y, z) |
|------|-------------------|
| Base plate | 80 × 60 × 5 mm at the origin |
| Fin k = 0..3 | 6 × 60 × 20 mm, x-offsets 8 + 18k mm |

Each box is closed; the soup is five disconnected shells by construction.
Import admits it through the quarantine path with `--max-hole-edges 0`
(zero boundary edges — checkable arithmetic, not a claim of manufacturability).

## The cooling declaration

- **Fan system**: one series bank whose curve is a two-point fixture line
  (0 Pa at 0.002 m³/s, 120 Pa at zero flow). The `source` string says
  plainly that this is NOT manufacturer data — an example may teach the
  schema with a declared synthetic curve; a product study may not.
- **Vent**: 0.0004 m² — the sum of the three inter-fin channels
  (3 × 12 mm × 20 mm).
- **Airflow leakage**: 0.0001 m² declared bypass (an imperfect shroud).
  Absence here would make the flow-network stage refuse rather than invent
  a path; declaring a small honest bypass is what lets the operating point
  execute.
- **Power**: 3 W at duty 1.0 dissipated in the metal region — a volumetric
  simplification of a chip heating the base, documented as such.

## What solve does with this

```
import-verify → assign → material-resolve → flow-network ✓
conduction ✗ (typed gap: bead frankensim-s93ej)
```

The flow-network stage computes the fan/vent/leakage network's operating
point with interval certification — the only stage in the pipeline whose
output currently carries outward-rounded arithmetic authority end to end.
Conduction then refuses with code `cli-solve-stage-gap`, naming its owner;
report and package refuse the same way. See
`examples/cooling-enclosure/README.md` for receipt anatomy and
`examples/heated-plate/README.md` for the minimal schema tour.
