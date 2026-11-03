# Measured gold adsorption on a specific Calgon carbon fraction

The `calgon-dg11-53um-gold-cyanide-martinez2023-reference` profile preserves
published mass-basis fit parameters and four paired observations. It does not
supply a dry-carbon capacity or a CIL rate coefficient.

## Source and specimen

[Martínez-Peñuñuri et al., Materials 16(14), 4961 (2023)](https://doi.org/10.3390/ma16144961),
Tables 1–3 and sections 2.1, 2.4 and 2.7.1, reports mortar-ground, sieved
53 µm Calgon DG-11 coconut-shell carbon, lot HB184AN-4, in synthetic
1000 mg/L NaCN solution at 25 °C and pH 10.5–11. The nominal initial bath
contains 10 mg/L Au in 0.5 L. The laboratory carbon masses are not explicitly
qualified as dry mass.

The source's 53 µm Langmuir fit gives `qmax = 39.5257 mg Au/g carbon`,
`KL = 3.9531 L/mg Au`, and linearized `R² = 0.9737`. The article generally
prefers Freundlich. Its separate observed pairs are:

| Reported carbon mass (g) | Equilibrium Au (mg/L) | Reported Au loading (mg/g carbon) |
|---:|---:|---:|
| 0.010 | 9.258 | 37.100 |
| 0.050 | 5.766 | 42.340 |
| 0.125 | 1.837 | 32.652 |
| 0.250 | 0.244 | 19.512 |

One observation exceeds the fitted asymptote. Retaining that disagreement
prevents treating a convenient fit as a measured hard site capacity.

## Quantities and boundaries

The pack contains 12 claims: three published fit quantities, the explicitly
derived half-saturation concentration, and the eight observed values above.
All uncertainties remain `Unstated`. Each observation is selected at its exact
reported carbon mass; no measured-curve interpolation is provided.

| Property | Normalized unit |
|---|---|
| `langmuir-saturation-mass-loading` | kg Au/kg reported carbon |
| `langmuir-affinity-mass-basis` | m³/kg Au |
| `langmuir-half-saturation-mass-concentration` | kg Au/m³ solution |
| `langmuir-linearized-fit-r-squared` | 1 |
| `observed-equilibrium-gold-mass-concentration` | kg Au/m³ solution |
| `observed-gold-mass-loading` | kg Au/kg reported carbon |

From `q = qmax KL C/(1 + KL C)`, the derived half-saturation mass
concentration is `1/KL = 0.2529660266626192 mg Au/L`.
This reciprocal is an algebraic conversion, not a new experiment.
The reported equilibrium concentration span is 0.244–9.258 mg/L; this does
not establish a validated extrapolation to zero concentration or to arbitrary
ore liquor.

Every query must identify the specimen, 53 µm sieve fraction, synthetic
cyanide medium, preparation NaCN concentration, temperature, reported pH
limits, engineering reference use, and `source-dry-carbon-mass-basis-established=0`.
Use `source-published-langmuir-fit=1` for fit quantities and zero for
the paired observations. Their additional `source-carbon-mass` coordinate
is in kg after normalization. The exact per-case pH is unreported:
`source-reported-pH-lower=10.5`, `source-reported-pH-upper=11` and
`source-exact-pH-known=0` retain that context. These are exact source-context
coordinates, not a constant-property model over a pH interval.

The preparation NaCN concentration is not a measurement of free CN⁻ at
equilibrium. The profile neither supplies metal competition nor qualifies
other particle sizes, industrial 6×12 carbon, natural preg-robbing material,
attrition losses, regeneration, desorption kinetics, or finite-time capture.

## Connecting to the existing equilibrium owner

`AqueousLangmuir` requires `saturation-loading` in mol Au/kg **dry**
adsorbent and `half-saturation-concentration` in mol Au/m³. Those exact
properties are deliberately absent here.

A future adapter needs an independently supplied dry-mass fraction
`f_dry = dry mass / reported mass`, an explicit Au molar mass `M_Au`,
and a separately qualified model domain. The unit conversion would be:

- `qmax_molar_dry = qmax_mass / (f_dry M_Au)`.
- `c_half_molar = c_half_mass / M_Au`.

No default dry fraction is justified by this profile. The conversion would
still preserve the published fit's disagreement with observations. Its
parameters must not be presented as a validated physical capacity.

The paper's pseudo-second-order kinetic fit is a separate law. It does not
supply the local linear-driving-force coefficient used by the CIL owner.

## Compile and verify

```sh
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/calgon-dg11-53um-gold-cyanide-martinez2023-reference/manifest.tsv \
  --out /tmp/calgon-dg11-53um.fsmatpk

cargo test -p xtask --test mining_calgon_gold_sources
```

The focused checks verify literal normalization and receipts, the four
nominal-bath gold mass accounts, reciprocal affinity, the retained fit
discrepancy, and refusals for unsupported source conditions or invented
dry-molar and kinetic properties.
