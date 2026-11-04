# Specimen-specific Bond ball-mill tests

`joanna-hosco-sgs-bond-2011` retains five composite test rows from the
[Joanna Gold Project technical report](https://www.sec.gov/Archives/edgar/data/913955/000113717112000140/report.htm),
effective 2011-12-31, section 13.1.3, Table 13.2, page 90. Each row supplies
work index, feed F80, product P80 and closing screen. The source explicitly
uses metric tonnes. Section 13 identifies this as historical PFS testwork.

The six required source coordinates retain the row identity, unknown test
temperature, pressure and preparation details, historical-reference status,
and lack of operating-ore qualification. Integer selectors 1 through 5 choose
distinct specimens; no interpolation or average work index is provided.
The material source uses `kW*h/Mg`: a megagram is exactly 1000 kg, preserving
the printed kWh/metric-tonne value while the compiler records its SI transform.
Uncertainty remains unstated; displayed precision is not a confidence interval.

## Compile and use

```bash
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/joanna-hosco-sgs-bond-2011/manifest.tsv \
  --out /tmp/joanna-hosco-sgs-bond-2011.fsmatpk
cargo run -p fs-material --example bond_grinding_reference -- \
  --material /tmp/joanna-hosco-sgs-bond-2011.fsmatpk
cargo test -p xtask --test mining_bond_sources
```

The example evaluates the existing Bond relation at each reported F80/P80
pair. Its 1000 kW net-pinion power and 90% utilization are explicitly authored
inputs. Calculated energy and dry throughput are estimates, not observations.

`comminution::source::resolve_bond_ball_mill` requires an explicit
`BondSizeDomain` in addition to the source query. A laboratory test point does
not define an operational size range. The resolver retains all four original
claims, their observation, uncertainty and receipts inside `BondGrindingCard`;
the existing blend, inverse and grinding-circuit consumers keep that card.
It refuses a short-ton interpretation of this metric source, mixed specimen
fields, unsupported conditions and missing measurements. It does not infer a
deposit-wide hardness, grinding kinetics, liberation, recovery or wear law.
