# Oxygen partitioning and transfer in mineral leaching

The matched NBS gas and dissolved-oxygen profiles now feed
`fs_thermochem::stoichiometry::thermodynamics::gas_liquid::IdealGasLiquidPartition`.
It evaluates an explicitly ideal gas/aqueous equilibrium and a finite sealed
headspace partition. Oxygen-transfer rates remain a separate constitutive model.

## Original sources and reference basis

The data are the original O2 gas and neutral aqueous `ao` entries in
Wagman et al., *The NBS Tables of Chemical Thermodynamic Properties*, Table 1:O,
printed page 2-37, as distributed in the
[NIST digitization with 1989 errata](https://doi.org/10.18434/M32124).
See the [corrected numerical CSV](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv)
and [original publication](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf).

| Profile | Species | Formation H, J/mol | Direct formation G, J/mol | Standard S, J/(mol K) |
| --- | --- | ---: | ---: | ---: |
| `oxygen-gas-nbs1989-thermochemistry` | O2(g) | 0 | 0 | 205.138 |
| `oxygen-aqueous-nbs1989-thermochemistry` | O2(aq) | -11700 | 16400 | 110.9 |

Both retain molar mass 0.0319988 kg/mol, 298.15 K, 100000 Pa and
`source-elemental-reference-nbs1989`. Every selected datum retains its original
receipt and `Unstated` uncertainty. The direct-G route uses the tabulated Gibbs
column; it does not replace it with formation H minus a species entropy term.

The exactly conserved equation is `O2(g) -> O2(aq)`. Its direct standard reaction
G is 16400 J/mol. Therefore `ln K = -16400/(R*298.15) = -6.615685469396972`
and `K = 0.0013391965031285026` under the declared source conventions.

## What the ideal calculation means

Let `n_ref` be the **amount at unit aqueous activity**. The caller independently
declares it, with `a_aq=n_aq/n_ref`. For the retained ideal molality convention,
`n_ref` is solvent mass in kilograms times 1 mol/kg. It is not the solution
volume, and unit activity is not 1 mol/m3 or 1 mol/L.

With the additional ideal-gas assumption, gas activity is `p_O2/p_ref`, where
`p_ref=100000 Pa`. Equilibrium gives

```text
n_aq = n_ref * K * p_O2 / p_ref
```

For example, an explicitly declared oxygen partial pressure of 21000 Pa and
`n_ref=1 mol` give 0.00028123126565698554 mol dissolved O2. This is an ideal
calculation from the original thermochemical reference, not a measured saline
slurry saturation concentration. A liquid volume is needed separately if the
consumer wants concentration for a rate law.

For a closed, rigid ideal headspace of volume `V_g`, the finite oxygen pool is
partitioned by two pressure capacities:

```text
C_gas = V_g/(R*T)
C_aq  = n_ref*K/p_ref
p_O2  = n_total/(C_gas+C_aq)
n_gas = C_gas*p_O2
n_aq  = C_aq*p_O2
```

The result retains both amounts, partial pressure, the original model,
component-balance residual and logarithmic equilibrium residual. Empty total
oxygen gives zero stocks and pressure, with an undefined logarithmic residual.
Positive trace stocks remain separately represented or arithmetic is refused.
Numerical residuals are not uncertainty bounds. Gas composition, water-vapor
pressure, nonideal fugacity/activity, salinity and temperature dependence are
not inferred.

## Equilibrium and transfer rate are separate inputs

These source profiles supply the equilibrium partition, **not** a mass-transfer
coefficient, gas-flow setting or ore-specific oxygen demand. A local transfer
law such as `dn_aq/dt = kLa*(n_aq_equilibrium-n_aq)` needs its own independently
declared inverse-second coefficient and source applicability. The equilibrium
amount must use the actual gas partial pressure at that rate evaluation.
An oxygen-saturated solution must not silently imply an infinite oxygen supply;
finite gas consumption and supplied gas must have explicit accounts.

## Compile the sources and run the focused integration

```bash
oxygen_pack_dir="$(mktemp -d)"
for pack in oxygen-gas-nbs1989-thermochemistry oxygen-aqueous-nbs1989-thermochemistry; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$oxygen_pack_dir/$pack.fsmatpk"
done
cargo test -p xtask --test mining_oxygen_partition
```

The test compiles these actual TSVs, verifies source values and receipts,
checks the independently calculated equilibrium constant, reconstructs the
finite-headspace mass balance, and exercises solvent-basis scaling, empty/trace
stocks and source/domain refusals. These checks establish the declared model's
arithmetic and source linkage; they do not calibrate leach kinetics or forecast
plant recovery.

## Couple finite oxygen exchange to gold cyanidation

`fs_material::cyanide_leach::gas_transfer::AeratedGoldCyanidation` attaches
this partition and an independently declared `OxygenTransfer` card to an
original `GoldCyanidationBatch`. At every numerical stage, the original
chemical owner recalculates free CN, pH and the gold dissolution rate `r`.
For cumulative dissolved gold `x`, liquid oxygen `D`, gas oxygen `G` and
signed transfer `J`, the coupled accounts are

```text
J = kLa*(D_equilibrium(G)-D)
dx/dt = r
dD/dt = J-r/4
dG/dt = -J
```

`advance` can disable gas contact while existing liquid oxygen continues to
support dissolution. `dose_oxygen` meters an explicit amount into the finite
headspace between advances; it changes neither liquid inventories nor elapsed
time. This is a fixed-volume headspace with explicit oxygen additions, not
continuous air flow or an infinite oxygen reservoir. Temperature, liquid
volume and solvent-activity basis remain fixed. Original source support is
checked at every evaluated stage. Numerical error estimates remain distinct
from source/model uncertainty.

The `aerated_gold_leach` example uses seven actual NBS packs. Compile them and run:

```bash
aerated_pack_dir="$(mktemp -d)"
for pack in hydrogen-cyanide-aqueous-nbs1989 hydrogen-ion-aqueous-nbs1989 cyanide-ion-aqueous-nbs1989 hydroxide-ion-aqueous-nbs1989 water-liquid-nbs1989-thermochemistry oxygen-gas-nbs1989-thermochemistry oxygen-aqueous-nbs1989-thermochemistry; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$aerated_pack_dir/$pack.fsmatpk"
done
cargo run -p fs-material --example aerated_gold_leach -- \
  --materials "$aerated_pack_dir"
```

The authored example uses `k=0.5 m6/(mol2 s)` and `kLa=0.01/s`, with 0.001 mol
accessible Au, 0.000001 mol initial dissolved O2 and 0.0001 mol headspace O2.
After 120 s without contact, dissolved Au is approximately 2.924e-6 mol and
remains below the initial-DO-only ceiling of 4e-6 mol. With finite headspace
exchange, the corresponding amount is approximately 1.969e-4 mol. Contact
then stops for 60 s: exchange is zero while reaction continues. A metered
0.0002 mol headspace dose followed by 240 s of resumed contact produces
approximately 7.369e-4 mol dissolved Au. These results describe this synthetic
rate/inventory scenario, not a calibrated ore or plant forecast.

Every stage prints gas/liquid O2, oxygen partial pressure, free CN, pH, Au
inventories, independently integrated transfer and reaction consumption, and
Au/CN/oxygen/charge residuals. The final gas plus liquid plus consumed oxygen
matches the original plus metered 0.000301 mol oxygen account. The example
does not infer air purity, bubble hydrodynamics, heat feedback, competing
cyanide complexes or carbon adsorption chemistry.

## Use measured oxygen solubility in the coupled leach

`oxygen-water-henry-sander2023-rettich2000` supplies the directly measured
oxygen-in-water reference **Hcp = 1.3e-5 mol/(m3 Pa) at 298.15 K**, using the
Rettich, Battino and Wilhelm (2000) row in
[Sander's normalized oxygen compilation](https://henrys-law.org/henry/casrn/7782-44-7).
The original measurement paper is
[J. Chem. Thermodynamics 32, 1145–1156](https://doi.org/10.1006/jcht.1999.0581);
the compilation is [Sander (2023)](https://doi.org/10.5194/acp-23-10901-2023).
The imported datum is one rounded reference point, with `Unstated` uncertainty.
It does not supply kLa, a finite-pressure validity range, a temperature curve,
or a saline/slurry correction.

`gas_transfer::OxygenHenrySource::from_pack` resolves the original typed
`henry_molar_solubility` claim and retains its complete pack, hash, typed query,
answer and replayable receipt. `OxygenTransfer::from_henry_source` connects it
to the same finite gas/liquid/reaction solver and a separately supplied kLa
card. This route does not need oxygen formation Gibbs profiles or invent an
elemental reference convention. The existing `from_card` Gibbs route remains
available; `equilibrium()` returns `Some` for that route and `None` for a
measured Henry source, which is accessible through `henry_source()`.

The IUPAC amount-concentration/partial-pressure convention gives

```text
p_O2 = n_gas*R*T/V_gas
c_O2_equilibrium = Hcp*p_O2
n_O2_equilibrium = Hcp*p_O2*V_liquid
```

The multiplier is the physical solution volume, independent of the solution's
amount-at-unit-activity basis. The selected claim must own the molecular-O2,
pure-water, infinite-dilution, engineering-reference and unqualified-specimen
flags. The consumer admits only that declared reference context plus typed
temperature. A claim that additionally depends on pressure, concentration,
ionic strength or another unsupported coordinate is refused. Every numerical
stage re-queries the same original claim with the actual solution temperature
and all five actual declared flags, preserving the source's semantic kelvin
descriptor. Each accepted state carries its resulting `henry_answer` receipt.

Pure-water infinite-dilution solubility applied to a cyanide-bearing solution is
an **explicit engineering reference assumption**, not electrolyte calibration.
The example prints that assumption and keeps its ore rate, kLa and inventory
choices visibly synthetic. Compile the new profile and select this route:

```bash
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/oxygen-water-henry-sander2023-rettich2000/manifest.tsv \
  --out "$aerated_pack_dir/oxygen-water-henry-sander2023-rettich2000.fsmatpk"
cargo run -p fs-material --example aerated_gold_leach -- \
  --materials "$aerated_pack_dir" \
  --henry-pack "$aerated_pack_dir/oxygen-water-henry-sander2023-rettich2000.fsmatpk"
```

This option loads only the five NBS cyanide/acid/base profiles from the example's
list, followed by the independent measured Henry pack. It exercises the same
finite-headspace, contact-stop, metered-dose and restart sequence, printing
the actual Hcp receipt alongside the chemical and oxygen-accounting results.
Consumer tests independently check finite exchange in both directions with
unequal solution-volume/activity bases, coupled reaction, retained clock and
dose accounts, receipt replay, and original source/context refusals.
