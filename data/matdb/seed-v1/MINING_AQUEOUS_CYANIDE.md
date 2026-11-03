# Aqueous cyanide, water dissociation and dissolved oxygen

Three profiles connect the existing six NBS gold-cyanidation profiles to
explicit aqueous reactions. They add **12 source claims**: molar mass, standard
formation enthalpy, direct formation Gibbs energy and standard entropy for each
new species. All nine profiles retain the same
`source-elemental-reference-nbs1989` evaluation.

| New profile | Source formula/state | M, g/mol | Formation H, kJ/mol | Formation G, kJ/mol | S, J/(mol K) |
| --- | --- | ---: | ---: | ---: | ---: |
| `hydrogen-cyanide-aqueous-nbs1989` | HCN, ao | 27.0259 | 107.1 | 119.7 | 124.7 |
| `hydrogen-ion-aqueous-nbs1989` | H+, ao | 1.0080 | 0 | 0 | 0 |
| `oxygen-aqueous-nbs1989-thermochemistry` | O2, ao | 31.9988 | -11.7 | 16.4 | 110.9 |

## Source and chemical identity

These are literal rows from the
[NIST NBS tables with 1989 errata](https://doi.org/10.18434/M32124),
using its [official corrected CSV](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv).
The [original NBS tables](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf)
place dissolved O2 on printed page 2-37, H+ on 2-38, and HCN on 2-89.
Sections 5.5-5.8 define the reference states; section 9.5 illustrates
acid-ionization and gas-solubility equilibrium constants.
The [NIST digitization report](https://doi.org/10.6028/jres.125.007)
describes preservation of the source tables and errata.

All claims are point values at **298.15 K and 100000 Pa**. The original
oxygen/hydrogen tables were prepared in 1963 and the carbon table in 1966.
The NBS/1989 name identifies the distributed evaluation; it does not date a
new measurement. Each claim keeps `Unstated` uncertainty because no per-row
probability distribution or cross-species covariance is supplied.

**HCN uses the undissociated `ao` row.** The adjacent `ai` row represents
the ionized H+ plus CN- reference component; its H/G/S values are
150.6 kJ/mol, 172.4 kJ/mol and 94.1 J/(mol K). Using that row as molecular HCN
would incorrectly erase the dissociation energy. The new profile is also
distinct from HCN gas, pure liquid HCN, total cyanide and metal-bound cyanide.

**Dissolved O2 is an aqueous neutral species.** Its nonzero formation H/G
must remain distinct from the zero formation H/G of reference O2 gas.
An aqueous amount is a finite dissolved inventory; it does not include a
gas bubble or imply replenishment from the atmosphere.

The neutral `ao` standard state is a hypothetical ideal solution at unit
molality, with activity referred to infinite dilution. It is neither a
measured 1 mol/kg solution nor a 1 mol/L concentration standard state.
HCN and O2 use `standard-molar-entropy`; H+ uses
`conventional-standard-molar-entropy` and the same H+-zero convention as the
existing CN- and OH- packs. H+ zeros define the ionic reference, not zero
absolute hydration properties or a measured pH.

The existing USGS H+ profile remains a separate evaluation. Its coincident
zero values are not permission to relabel a mixed NBS/USGS reaction.

## Conserved reactions and standard results

Three `MineralReaction` variants pass through the existing exact atom and
charge checker:

```text
HydrogenCyanideDissociation:
    HCN(aq) -> H+(aq) + CN[-](aq)

WaterDissociation:
    H2O(l) -> H+(aq) + OH[-](aq)

GoldCyanidationDissolvedOxygen:
    4 Au(s) + 8 CN[-](aq) + O2(aq) + 2 H2O(l)
        -> 4 Au(CN)2[-](aq) + 4 OH[-](aq)
```

The original gas-O2 `GoldCyanidation` variant remains distinct.
Reverse dissociation is represented by a negative equilibrium extent.

| Reaction | Standard H, J/mol extent | Standard S, J/(mol extent K) | Direct G, J/mol extent | H - T S, J/mol extent |
| --- | ---: | ---: | ---: | ---: |
| HCN dissociation | 43500 | -30.6 | 52700 | 52623.39 |
| Water dissociation | 55836 | -80.66 | 79885 | 79884.779 |
| Dissolved-O2 Elsner | -572216 | -548.12 | -407118 | -408794.022 |

The direct formation-G and caloric columns remain separate source routes.
Their differences are exposed and are not attributed solely to rounding,
removed by offsets, or used to alter a source datum.
The new aqueous gold endpoint retains the historical NBS gold-complex row;
the distinct later USBM evaluation discussed in
[MINING_GOLD_CYANIDATION.md](MINING_GOLD_CYANIDATION.md) is not substituted.

For the direct-G route, ln K = -G/(RT). At the declared source state:

- HCN: Ka = 5.85274292644451e-10 and pKa = 9.232640551245206.
- Water: Kw = 1.0110054626220328e-14 and pKw = 13.995246497841048,
  with pure liquid water assigned unit activity.
- O2(g) -> O2(aq): standard G = 16400 J/mol and
  K = 0.0013391965031285008.

The oxygen constant follows by subtracting the dissolved-O2 Elsner equation
from the gas-O2 equation. It relates aqueous activity to gas fugacity divided
by the source reference pressure. It is not directly a mol/m3 saturation
concentration or a gas-transfer coefficient.

## Consumers and meaningful checks

`ConservedReaction::equilibrate_ideal` already accepts these aqueous and
pure condensed phases for a **single reaction**. Each aqueous participant
requires an explicit `IdealActivityModel::IdealDilute` amount at unit
activity, giving a = n/n_ref. The caller must declare the matching fixed
solvent basis; neither tank volume nor molarity is silently substituted for
the NBS molal standard state. Pure water uses `PureCondensed`, with unit
activity and a finite stock.

The source integration target compiles all nine real packs through
`xtask matdb-pack`, then uses the original public owners:

```bash
aqueous_pack_dir=$(mktemp -d)
for pack in hydrogen-cyanide-aqueous-nbs1989 hydrogen-ion-aqueous-nbs1989 oxygen-aqueous-nbs1989-thermochemistry cyanide-ion-aqueous-nbs1989 water-liquid-nbs1989-thermochemistry hydroxide-ion-aqueous-nbs1989 gold-solid-nbs1989-thermochemistry dicyanoaurate-ion-aqueous-nbs1989 oxygen-gas-nbs1989-thermochemistry; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$aqueous_pack_dir/$pack.fsmatpk"
done

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction hydrogen-cyanide-dissociation --gibbs --equilibrium \
  "$aqueous_pack_dir/hydrogen-cyanide-aqueous-nbs1989.fsmatpk" \
  "$aqueous_pack_dir/hydrogen-ion-aqueous-nbs1989.fsmatpk" \
  "$aqueous_pack_dir/cyanide-ion-aqueous-nbs1989.fsmatpk" 0.0270259

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction water-dissociation --gibbs --equilibrium \
  "$aqueous_pack_dir/water-liquid-nbs1989-thermochemistry.fsmatpk" \
  "$aqueous_pack_dir/hydrogen-ion-aqueous-nbs1989.fsmatpk" \
  "$aqueous_pack_dir/hydroxide-ion-aqueous-nbs1989.fsmatpk" 0.0180154

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction gold-cyanidation-dissolved-oxygen --gibbs --equilibrium \
  "$aqueous_pack_dir/gold-solid-nbs1989-thermochemistry.fsmatpk" \
  "$aqueous_pack_dir/cyanide-ion-aqueous-nbs1989.fsmatpk" \
  "$aqueous_pack_dir/oxygen-aqueous-nbs1989-thermochemistry.fsmatpk" \
  "$aqueous_pack_dir/water-liquid-nbs1989-thermochemistry.fsmatpk" \
  "$aqueous_pack_dir/dicyanoaurate-ion-aqueous-nbs1989.fsmatpk" \
  "$aqueous_pack_dir/hydroxide-ion-aqueous-nbs1989.fsmatpk" 0.787868

cargo test -p xtask --test mining_aqueous_cyanide_sources
```

The three mass arguments declare one mole of HCN, one mole of water and four
moles of gold actually reacted, respectively. They are thermodynamic accounting
bases, not predicted equilibrium conversions. The CLI's `--equilibrium` flag
reports the standard caloric constant; the finite-stock equilibrium calculations
below exercise `equilibrate_ideal` directly.

Its four checks cover:

1. Literal source rows, typed SI conversion, original uncertainty and receipts,
   neutral versus ionic entropy, and refusal of wrong phase, temperature or
   reference convention.
2. HCN dissociation against an independent finite-stock quadratic. Added
   proton suppresses free cyanide; scaling both stocks and the solvent basis
   preserves activities and scales the equilibrium amount.
3. Water dissociation against sqrt(Kw), with charge and solvent accounting,
   solvent-basis scaling, and an inventory-blocked result without water.
4. Dissolved-O2 gold H/G/S, its distinct gas-reference difference, and a finite
   dissolved-oxygen material ceiling. Starting with 4 mol Au, 8 mol CN-,
   0.00025 mol O2(aq) and 55.5 mol water, at most 0.00025 mol extent is available:
   only 0.001 mol Au can enter the complex. A request beyond that ceiling refuses.

The scalar ideal-equilibrium owner does not solve several reactions jointly.
Sequentially equilibrating HCN and water once is not a coupled speciation
solution. A joint calculation needs shared inventories and simultaneous
reaction closure; counterions or a declared background-charge inventory are
also required when initial dissolved species do not balance charge.
Activity coefficients, copper-cyanide competition, HCN volatilization,
gas-transfer kinetics and finite leach rates remain separate physical models.

## Coupled pH and free cyanide

`stoichiometry::thermodynamics::acid_base::IdealAcidBaseSystem` now solves
water and one or more conserved protonation families simultaneously. Each
one-proton dissociation retains its original equation, selected Gibbs claims
and reference conditions. All families share one proton activity; a monotone
charge-balance solve includes the explicit spectator-ion inventories. The
HCN/CN- analytical total is the **uncomplexed pool**, not total cyanide,
WAD cyanide or cyanide already in metal complexes.

After compiling the source packs above, run:

```bash
cargo run -p fs-thermochem --example mineral_speciation -- \
  --materials "$aqueous_pack_dir"
```

The example uses an explicit 10 kg solvent basis, so the amount at unit
molal activity is 10 mol. It compares pure water, 0.1 mol molecular HCN,
the same analytical pool with 0.1 mol Na+, additional strong-base equivalents,
and acid equivalents represented by explicit Cl-. These are equilibrium
scenarios on declared stocks, not salt-dissolution or operating-dose advice.
With the retained NBS constants, the matching cyanide-salt case gives
pH approximately 10.6049 and 0.0959293 mol free CN-. Explicit acid equivalents
change both pH and free cyanide while the HCN+CN- analytical total remains
0.1 mol.

Every aqueous activity uses the same `n / n_ref` approximation with activity
coefficient one. Water activity and the solvent basis remain fixed; this
owner does not claim a closed finite-solvent H/O balance. It refuses conflicting
shared sources, overlapping pools, double-counted spectators, unsupported pH
brackets, numerical range loss or exhausted work. The reported bracket and
charge/mass-action residuals are numerical diagnostics, not uncertainty bounds
or a validated concentrated-electrolyte model.

Eight focused analytical test groups cover simultaneous mono/polyprotic
families, water/strong-ion quadratic solutions, counterions, dilution and
extensive scaling, tiny positive members, source/topology admission and bounded
atomic refusals. `xtask/tests/mining_speciation_sources.rs` adds three tests
through the actual five-pack compiler and consumer, including independent
HCN/cyanide-salt bounds retaining the water contribution.
