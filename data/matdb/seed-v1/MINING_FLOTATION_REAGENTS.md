# Flotation frother reference: MIBC

`mibc-arkema-flotation-reference` adds five source claims for the actual Arkema
MIBC product, CAS 108-11-2. The [manufacturer's technical data](https://www.arkema.com/global/en/products/product-finder/product/thiochemicals/solvents/methylisobutylcarbinol/)
identify its use as a flotation frother and provide these reference quantities:

| Property | Retained value | Basis |
| --- | ---: | --- |
| Liquid density endpoints | 806 and 808 kg/m³ | 293.15 K; no midpoint or probability law |
| Minimum product purity | 0.985 | Mass fraction; not actual lot composition |
| MIBC in water-rich saturated phase | 0.017 | Mass fraction at 293.15 K |
| Water in MIBC-rich saturated phase | 0.058 | Mass fraction at 293.15 K |

The two solubilities have different phase and component meanings. They cannot
be interpreted as g/L, exchanged, or used as frother dosing recommendations.
The page's boiling-pressure entry has an apparent unit inconsistency, so it
is excluded. Pressure, qualified lot, temperature response, surface tension,
bubble/coalescence behavior and mineral-specific flotation kinetics are absent.

Use the ordinary compiler or catalog; no special reagent registry is needed:

```bash
cargo run -p xtask -- matdb-catalog \
  --pack mibc-arkema-flotation-reference \
  --property density-reported-range-minimum \
  --property density-reported-range-maximum \
  --at temperature=293.15 --at source-arkema-mibc-product=1 \
  --at engineering-reference-model=1 --at source-lot-qualified=0 \
  --at source-test-pressure-known=0
```

For a declared product volume flow, mass flow is volume flow times the selected
density endpoint. The focused test uses an authored 1 L/min feed at 60 dry
tonnes/h to obtain 0.806–0.808 kg product/t. That is an accounting example;
it is not an operating dose or evidence of flotation performance. Product purity
remains separate from delivered product mass.

`cargo test -p xtask --test mining_flotation_reagent_sources` compiles the actual
manifest, verifies round-trip/source receipts and catalog discovery, checks
both phase compositions, and refuses unsupported temperature, lot, units and
recovery properties. Every claim retains unstated statistical uncertainty.
