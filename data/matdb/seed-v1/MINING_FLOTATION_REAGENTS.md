# Flotation product references: frother, collectors and copper sulfate

## MIBC frother

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

## Xanthate collectors and copper sulfate pentahydrate

Four additional product-specific packs supply eight numerical claims for
collector-feed accounting and copper-sulfate reference identity. They are
discoverable by the same compiler and catalog as the MIBC pack.

| Pack | Source quantity retained | Meaning |
| --- | --- | --- |
| `pax-ixom-solid-reference` | PAX 0.90; KOH 0.015 | Reported mass fractions in the commercial solid; neither is a certified lot assay or a bound |
| `pax22-ixom-solution-reference` | PAX 0.22 | Nominal mass fraction on the liquid product label, explicitly W/W |
| `sibx-ixom-solid-reference` | SIBX ≥0.90; NaOH ≤0.01 | Separate lower and upper component mass-fraction bounds |
| `sibx-ixom-solid-reference` | Approximately 500 g/L at 293.15 K | Reported water solubility; solvent-volume versus solution-volume denominator unspecified |
| `copper-sulfate-pentahydrate-ixom-reference` | CuSO₄·5H₂O ≥0.98 | Minimum pentahydrate component mass fraction; five waters per formula unit |
| `copper-sulfate-pentahydrate-ixom-reference` | Relative density 2.284 | Dimensionless ratio; source temperature, pressure and reference density unspecified |

The original supplier sources are:

- [IXOM solid PAX SDS, 000033326701, revision 9, 1 November 2023](https://sds.ixom.com/pdf/shess-en-mtr-UGHS-000033326701.pdf),
  section 3, and the associated [solid product page](https://www.ixom.com/p/PAX-900_AU).
- [IXOM PAX22-BULK_AU liquid product page](https://www.ixom.com/p/PAX22-BULK_AU).
  Its 22% W/W product label is retained as nominal, with no invented tolerance.
- [IXOM SIBX SDS, 000033039301, revision 8, 24 November 2023](https://sds.ixom.com/pdf/shess-en-mtr-UGHS-000033039301.pdf),
  section 3; the same revision's [New Zealand SDS](https://sds.ixom.com/pdf/shess-en-mtr-ZGHS-000033039301.pdf),
  section 9, supplies the approximate water-solubility entry.
- [IXOM copper sulfate pentahydrate SDS, 000031007601, revision 8, 4 August 2022](https://sds.ixom.com/pdf/shess-en-mtr-UGHS-000031007601.pdf),
  sections 3 and 9.

Supplier PDF links can redirect to the current SDS catalogue; the product
codes and revision dates identify the cited disclosure. These selected facts
retain unstated statistical uncertainty. The repository does not redistribute
the source documents or upgrade an SDS disclosure into a shipment certificate.

### Mass accounting and product substitution

For an authored 20 kg product/h feed and 100 dry tonnes/h ore throughput,
the disclosed solid PAX fraction gives 0.18 kg PAX/t. The same product mass
flow of nominal 22% solution gives 0.044 kg PAX/t. Matching the first active
mass flow would require about 4.091 times the solution mass. These arithmetic
examples compare declared product inputs; they do not recommend a collector
dose or establish equivalent flotation performance.

Likewise, the SIBX lower bound and NaOH upper bound remain separate. Their
endpoints cannot reconstruct an exact batch composition. The copper-sulfate
lower bound describes **pentahydrate mass**; converting that mass to copper
or anhydrous CuSO₄ requires a separate formula/molar-mass calculation. The
unlisted mass of every product remains unidentified.

```bash
cargo run -p xtask -- matdb-catalog \
  --pack pax22-ixom-solution-reference \
  --property pax-solution-component-mass-fraction-nominal \
  --at source-ixom-pax22-solution-product=1 \
  --at source-product-label-nominal=1 \
  --at engineering-reference-model=1 --at source-lot-qualified=0

cargo run -p xtask -- matdb-catalog \
  --pack copper-sulfate-pentahydrate-ixom-reference \
  --property copper-sulfate-pentahydrate-product-component-mass-fraction-minimum \
  --at source-ixom-copper-sulfate-pentahydrate-product=1 \
  --at source-hydration-water-moles-per-formula-unit=5 \
  --at engineering-reference-model=1 --at source-lot-qualified=0
```

The collector packs supply no density for volume-feed calculations. Copper
sulfate's relative-density value remains dimensionless and cannot substitute
for bulk-bed or solution density. SIBX's approximate solubility requires its
293.15 K query, an acknowledged unknown volume basis and an acknowledged
approximate value; it cannot silently become a saturation mass fraction or
temperature-dependent law.

The integration test also compiles and round-trips each new pack, verifies
provenance receipts, repeats deterministic compilation, exercises ordinary
catalog queries and rejects other product identities, invented lot
certification, wrong hydrate states and unsupported thermochemical or
flotation-rate properties. Product availability and these component facts
are prerequisites for flowsheet reagent inventories. Mineral-specific
adsorption, activation, kinetics and recovery still require their own
experimental evidence and process models.
