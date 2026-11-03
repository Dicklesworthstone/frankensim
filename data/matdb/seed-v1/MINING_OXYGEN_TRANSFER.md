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
