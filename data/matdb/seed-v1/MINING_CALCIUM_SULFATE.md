# Calcium sulfate precipitation, hydration and lime demand

The compatible NBS calcium, sulfate, water and lime sources now connect to
phase-specific anhydrite and gypsum records. Four conserved reactions support
reagent and hydrate-water accounting, standard reaction heat, and an explicit
finite-inventory equilibrium along a selected free-ion precipitation coordinate.

These capabilities support calculations of dissolved calcium/sulfate stocks and
solid formation. They do not assign a fixed sulfate-removal percentage or infer
a commercial scale-deposition rate.

## Original source identities

The new records are literal rows from the
[NIST digitization of the corrected NBS tables](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv),
Table 94:Ca, Calcium (Prepared 1968), at **298.15 K and 100000 Pa**:

| Profile | Original row | M, g/mol | Formation H, kJ/mol | Formation G, kJ/mol | S, J/(mol K) |
|---|---|---:|---:|---:|---:|
| `anhydrite-solid-nbs1989-thermochemistry` | CaSO4, cr, insoluble anhydrite | 136.1416 | -1434.11 | -1321.79 | 106.7 |
| `gypsum-selenite-solid-nbs1989-thermochemistry` | CaSO4·2H2O, cr, selenite | 172.1724 | -2022.63 | -1797.28 | 194.1 |

Selenite is a crystalline variety of gypsum
([USGS Bulletin 1105, p. 2](https://pubs.usgs.gov/bul/1105/report.pdf)).
The original qualification is retained with both `source-phase-gypsum=1` and
`source-variety-selenite=1`. Anhydrite requires
`source-phase-anhydrite=1`. The adjacent NBS soluble alpha/beta CaSO4 and
hemihydrate rows are separate materials. The source description “insoluble”
does not imply a zero thermodynamic solubility.

Each pack has four claims with `Unstated` uncertainty, the original
`source-elemental-reference-nbs1989` axis, pure-solid reference flags, and
no specimen qualification. These records remain separate from the existing
USGS gypsum and anhydrite evaluation. Neither a caller-supplied flag nor a
species name changes a source's reference identity.

Compatible existing participants are:

- `calcium-ion-aqueous-nbs1989`: individual Ca2+, not total dissolved calcium.
- `sulfate-ion-aqueous-nbs1989`: free SO4(2-), not HSO4- or paired sulfate.
- `water-liquid-nbs1989-thermochemistry`: pure liquid water.
- `hydrogen-ion-aqueous-nbs1989`: conventional individual H+.
- `calcium-hydroxide-solid-nbs1989-thermochemistry`: the original crystalline
  hydrated-lime reference.

## Conserved endpoints and source-derived energies

`MineralReaction` provides these four separately named endpoints:

| Variant | Written equation |
|---|---|
| `GypsumPrecipitation` | Ca2+(aq) + SO4(2-)(aq) + 2 H2O(l) → CaSO4·2H2O(s) |
| `AnhydritePrecipitation` | Ca2+(aq) + SO4(2-)(aq) → CaSO4(s) |
| `AnhydriteHydration` | CaSO4(s) + 2 H2O(l) → CaSO4·2H2O(s) |
| `GypsumPrecipitationWithHydratedLime` | Ca(OH)2(s) + 2 H+(aq) + SO4(2-)(aq) → CaSO4·2H2O(s) |

The exact reaction owner checks atoms and charge, and the hydrate waters are
explicit participants. Shared definitions also compose with
`HydratedLimeNeutralization` through the existing `ReactionNetwork`.

At the declared standard state, source sums give:

| Endpoint | ΔH°, J/mol written extent | Direct ΔG°, J/mol written extent | ΔS°, J/(mol K) |
|---|---:|---:|---:|
| Gypsum precipitation | 1130 | -24912 | 87.28 |
| Anhydrite precipitation | 17990 | -23680 | 139.7 |
| Anhydrite hydration | -16860 | -1232 | -52.42 |
| Gypsum from hydrated lime and sulfate | -127270 | -154260 | 90.61 |

Direct ΔG uses the independent formation-Gibbs column. The separate caloric
route evaluates the conserved reaction sum ΔH - TΔS and retains its own source
receipts; small differences between the two routes are not overwritten.
Standard reaction heat is not a temperature-rise or process heat-duty model.

For example, a declared **0.00370474 kg of pure referenced Ca(OH)2 reacted**
is 0.05 mol written neutralization extent. It consumes 0.10 mol H+ and
0.05 mol free sulfate and produces 0.05 mol gypsum, with standard reaction
heat **-6363.5 J**. With initial inventories of 0.10 mol lime, 0.12 mol H+
and 0.10 mol sulfate, the inventory ceiling is 0.06 mol extent. That ceiling
does not assert equilibrium or actual completion.

## Finite selected-phase equilibrium

The existing `ConservedReaction::equilibrate_ideal` method now consumes these
source sets. The caller declares every participant's stock and activity law:

- Ca2+ and SO4(2-) use `IdealDilute`, with a common positive amount at unit
  activity. For the NBS unit-molality convention, 1 mol corresponds to the
  nominal 1 kg solvent basis.
- The selected solid and, for gypsum, liquid water use `PureCondensed`.
  Their activities are one, but their available amounts are finite.
- Positive precipitation extent removes calcium and sulfate in a 1:1 ratio.
  Negative extent dissolves the supplied solid. The feasible interval includes
  all hydrate-water and phase stocks.

For the gypsum dissolution direction, the independent source arithmetic is

```text
ΔG° = -553580 -744530 + 2(-237129) -(-1797280) = 24912 J/mol
Ksp = exp(-24912 / (R × 298.15)) ≈ 4.32122047e-5
```

For a 1 mol activity-reference amount and equal dissolved free-ion stocks, a
present gypsum phase gives approximately **0.006573599 mol each** Ca2+ and
SO4(2-). An initial 0.020 mol gypsum charge retains about 0.013426401 mol
solid; a smaller 0.001 mol charge dissolves completely and remains
undersaturated. Starting with 0.050 mol of each free ion instead precipitates
approximately 0.043426401 mol gypsum.

The common-ion example supplies 0.050 mol sulfate and explicitly retains
0.100 mol fixed Na+ spectator. Its dissolved calcium follows the independent
quadratic `s(s+0.05)=Ksp`, giving approximately **0.000849801 mol**.
The selected anhydrite equation has its own dissolution ΔG° of 23680 J/mol;
a separate anhydrite calculation gives about 0.008427946 mol of each ion.
These separate runs do not constitute competitive phase selection.

All figures above are authored ideal-model scenarios. The result retains the
original formation-source receipts, signed extent, final inventory, numerical
bracket, residual and interior/phase-exhaustion status. A phase-exhaustion
inequality is reported as such rather than forcing a saturation equality.

## Run the actual-source example

From the repository root:

```bash
sulfate_pack_dir=$(mktemp -d)
for pack in gypsum-selenite-solid-nbs1989-thermochemistry anhydrite-solid-nbs1989-thermochemistry calcium-ion-aqueous-nbs1989 sulfate-ion-aqueous-nbs1989 water-liquid-nbs1989-thermochemistry hydrogen-ion-aqueous-nbs1989 calcium-hydroxide-solid-nbs1989-thermochemistry; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$sulfate_pack_dir/$pack.fsmatpk"
done

cargo run -p fs-thermochem --example calcium_sulfate -- \
  --materials "$sulfate_pack_dir"

# Explicit common-ion scenario with its fixed counterion inventory:
cargo run -p fs-thermochem --example calcium_sulfate -- \
  --materials "$sulfate_pack_dir" --sulfate-mol 0.05 --sodium-mol 0.10

cargo test -p xtask --test mining_calcium_sulfate_sources
```

The default example prints all four source-derived reaction energies and
receipts, five finite-phase scenarios, and the declared lime-reaction account.
Custom options choose gypsum or anhydrite, free-ion stocks, solid stock,
finite water, activity-reference amount, and fixed Na+/Cl- spectator stocks.
The initial free ions and explicit spectators must balance charge. The
unchanged spectators contribute no reaction-energy term.

Results include calcium, sulfate, hydrate-water, charge and source-mass
residuals. The source integration tests independently verify the literal
eight new claims; all four reaction H/G/S values and source-mass closure;
reagent limits; saturation and the common-ion quadratic; finite phase/water
exhaustion; changed activity basis; phase/reference refusal; and solve-budget
refusal.

## Physical scope

This is a fixed-temperature, selected-phase, **one-reaction ideal free-ion
model**. Its activity-reference amount stays fixed even while the explicit
hydrate-water stock changes. It does not infer a solvent-volume law or correct
the molality reference for that change.

The existing [sulfate acid/base model](MINING_SULFATE_SPECIATION.md) and
[lime/pH model](MINING_LIME_EQUILIBRIUM.md) remain separate owners. Running
them sequentially and freezing their pH or sulfate fraction would not provide
a simultaneous gypsum/HSO4-/water equilibrium. Calcium-sulfate ion pairs,
nonideal activities, multiple competing solids, nucleation, transport,
deposition kinetics and measured plant removal require additional models.
