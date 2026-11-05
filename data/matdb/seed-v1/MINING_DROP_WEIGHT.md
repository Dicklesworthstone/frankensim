# JK drop-weight ore breakage references

`bueno-jkdwt-2021` provides **five ore-specific JKDWT calibrations and ten
claims**, connected to `fs_material::comminution::drop_weight`. This fills the
impact-energy response named in FrankenMining's comminution plan alongside the
existing Bond ball-mill models.

## Source and quantities

Bueno, Torvela, Chandramohan, Chavez Matus, Liedes and Powell (2021),
[*The double wheel breakage test*, Minerals Engineering 168, 106905](https://doi.org/10.1016/j.mineng.2021.106905),
Tables 5–7, supplies the original JKDWT specimen identities and independently
reported A and b fits. The article is CC-BY-4.0; this profile transcribes numerical
facts with original unit and scope annotations. It excludes the distinct SMC
and double-wheel fits.

The source identifies gold, polymetallic, nickel-zinc, polymetallic and
phosphate ores as A–E. Mine identities, temperature, pressure and coefficient
uncertainty remain unreported. No mineral name is substituted for these specimens.

The model follows the [JKTech / SMC Testing SCSE report, page 1](https://www.smctesting.com/documents/Reporting_JK_Drop_Weight_and_SMC_Test_Results.pdf):

`fraction passing one tenth of parent size = A * (1 - exp(-b * specific_energy))`

The compiler converts A from percent to a fraction and b from inverse
kWh/metric-tonne to kg/J. Thus `1 kWh/metric-tonne = 3600 J/kg`, and a reported
`b = 0.55/(kWh/metric-tonne)` becomes `0.55/3600 kg/J`. The inverse unit is written
`Mg/kW/h` in the source TSV so ordinary dimensional parsing performs this conversion.

Table 6's JKDWT settings span 0.1–2.5 kWh/t and five particle-size fractions
between 13.2 and 63 mm. Different sizes have different tested energies; their
overall envelope is **not a fully tested size-by-energy rectangle**. No
particle-size law or operational transfer is fitted here. The paper's Table 8
contains calculated apparent JKDWT t10 values at rounded energies, separately
from double-wheel measurements; the integration test compares the calculated
column with a 0.2 percentage-point allowance for the printed precision.

## Running the source-backed comparison

```bash
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/bueno-jkdwt-2021/manifest.tsv \
  --out /tmp/bueno-jkdwt-2021.fsmatpk
cargo run -p fs-material --example drop_weight_reference -- \
  --material /tmp/bueno-jkdwt-2021.fsmatpk
```

The example compares t10 at 1 kWh/metric-tonne and calculates the impact energy
for 20% t10 for each original specimen. It explicitly selects a 0.1–2.5
kWh/metric-tonne calculation interval. These outputs are model estimates.

`resolve_jk_drop_weight` selects both scalar parameters using the same typed
specimen query. Each card retains the original pack hash, claims, uncertainty,
observation and verified receipts. Missing parameters, mixed source support,
unknown specimens, unsupported operating conditions and claims of plant
qualification refuse. A and b cannot be recovered from their product alone.

`response` preserves both fractions around the one-tenth-size threshold, including
a small exponential retained tail when A is one. `energy_for_t10` inverts a
reachable target, checks the original declared energy interval and refuses targets
at or above A. The inverse becomes sensitive to rounded t10 near saturation.

A single t10 statistic is not a full product-size distribution. Applying it to
mill power or throughput still requires circuit configuration and a suitable
model; applying it to the existing population balance requires independently
supported daughter-size and selection-rate information. No ta, recovery,
liberation, liner wear or grinding-media consumption is inferred.
