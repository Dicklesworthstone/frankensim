# Ferric hydrolysis and lime treatment with named solid phases

The goethite and portlandite profiles complete a compatible source set for
precipitating ferric iron and accounting for the accompanying acid and lime.
They use the normal `matdb-pack` compiler and typed `NormalizedPack` queries.
The source values come from [USGS Bulletin 2131](https://doi.org/10.3133/b2131),
printed page 16, verified in the original report through its
[government-document mirror](https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/22/).

## Source values and phase identity

All values apply at **298.15 K and 100000 Pa**, with
`source-elemental-reference-usgs-b2131=1`. H is standard molar formation
enthalpy; G is directly tabulated standard molar formation Gibbs energy;
S is standard molar entropy of the named pure solid.

| Profile | Formula | Molar mass, kg/mol | H, J/mol | G, J/mol | S, J/(mol K) | Reported H magnitude, J/mol | Reported S magnitude, J/(mol K) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `goethite-solid-usgs-b2131-thermochemistry` | FeO(OH) | 0.088854 | -562600 | -491800 | 60.4 | 2100 | 0.6 |
| `portlandite-solid-usgs-b2131-thermochemistry` | Ca(OH)2 | 0.074093 | -986100 | -898000 | 83.4 | 1300 | 0.4 |

Each pack has six claims: molar mass, formation enthalpy, formation Gibbs
energy, entropy and the two reported H/S uncertainty magnitudes. Every claim
requires the named phase selector
(`source-phase-goethite=1` or `source-phase-portlandite=1`),
`source-phase-solid=1`, `source-pure-species=1`, `source-standard-state=1`,
`engineering-reference-model=1` and `source-specimen-qualified=0`, together
with typed temperature and pressure. These are reference phases, not measured
commercial reagents or sludge compositions.

Goethite FeO(OH) remains distinct from lepidocrocite, ferrihydrite and generic
Fe(OH)3. The reaction species `FeO(OH)(s)` explicitly denotes the goethite
endpoint in these helpers. The existing generic ferric-hydroxide pyrite
equation is unchanged. The USGS portlandite profile retains its own reference
identity alongside `calcium-hydroxide-solid-nist-mining`; source formulas do
not authorize relabeling one evaluation as another.

Bulletin 2131 printed page 3 describes a calorimetric uncertainty convention
of twice the standard deviation of the mean, without a universal coverage
probability or cross-species covariance. All claims retain `Unstated`
uncertainty semantics. The dimensional magnitude claims preserve source
information without treating it as independent Gaussian error or a certified
reaction-error interval.
The separately printed Gibbs-energy uncertainty magnitudes, 2100 J/mol for
goethite and 1300 J/mol for portlandite, remain in their Gibbs observations
with the same `Unstated` semantics. The direct G values retain the original
printed Gibbs column; they are not reconstructed from a species' formation
enthalpy and absolute entropy or adjusted to agree with a different compilation.

The compatible participants already in the database are
`iron-iii-ion-aqueous-usgs-b2131`, `hydrogen-ion-aqueous-usgs-b2131`,
`calcium-ion-aqueous-usgs-b2131` and `water-liquid-usgs-b2131`.
Their formation enthalpies are respectively -49900, 0, -543000 and
-285800 J/mol; their direct formation Gibbs energies are -16700, 0, -553600
and -237100 J/mol. The ions use the source's unit-activity, infinite-dilution
aqueous reference and explicit H+(aq)-zero convention; their entropy key is
`conventional-standard-molar-entropy`. Water uses the pure liquid phase.

## Conserved reaction endpoints, compatible heat and direct Gibbs energy

The hydrolysis equation is printed in
[USGS Open-File Report 2005-1148](https://pubs.usgs.gov/of/2005/1148/of2005-1148.pdf)
and is also the reverse of the goethite dissolution reaction in the
[USGS PHREEQC database](https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/html/final-90.html).
The equation supplies the reaction identity; all heats and direct Gibbs sums
below use the selected Bulletin 2131 rows, without mixing in PHREEQC's separate
numeric evaluation.

```text
FerricHydrolysisGoethite:
Fe3+(aq) + 2 H2O(l) -> FeO(OH)(goethite) + 3 H+(aq)
Delta H = +58900 J/mol of the written equation
Direct Delta G = -900 J/mol of the written equation

Existing HydratedLimeNeutralization, now with compatible USGS portlandite:
Ca(OH)2(portlandite) + 2 H+(aq) -> Ca2+(aq) + 2 H2O(l)
Delta H = -128500 J/mol of the written equation
Direct Delta G = -129800 J/mol of the written equation

FerricGoethitePrecipitationWithHydratedLime:
2 Fe3+(aq) + 3 Ca(OH)2(portlandite)
    -> 2 FeO(OH)(goethite) + 3 Ca2+(aq) + 2 H2O(l)
Delta H = -267700 J/mol of the written equation
Direct Delta G = -391200 J/mol of the written equation
```

Each variant builds through the exact atom/charge conservation checker. Two
hydrolysis extents plus three lime-neutralization extents cancel the intermediate
protons and four of the produced water molecules, giving the combined equation.
The source heat obeys the same Hess relation:
`2 * 58900 + 3 * (-128500) = -267700 J/mol`.
The independently tabulated Gibbs values obey
`2 * (-900) + 3 * (-129800) = -391200 J/mol`.
Original participant queries and verified receipts remain available in both
reaction results and in activity evaluations that borrow the standard Gibbs
result.

## Finite inventories and scope

For **10 mol Fe3+ and 6 mol portlandite**, lime limits the combined extent to
2 mol of the written equation. Advancing that extent leaves 6 mol Fe3+ and
no portlandite, and forms 4 mol goethite, 6 mol Ca2+ and 4 mol water. The
standard enthalpy change is -535400 J. A larger requested extent is refused.

The same bookkeeping can be expressed as explicit sequential hydrolysis and
neutralization stages. Starting with 2 mol Fe3+, 3 mol portlandite and 4 mol
water, hydrolysis first supplies the six protons consumed by the following
neutralization stage. The final water inventory is 6 mol, including the
4 mol initially present. Requesting both gross stage demands from the original
proton-free inventory is refused; a future reaction's products do not fund
the current stage. The combined endpoint remains available when its net
stoichiometry is the intended model.

Extents are caller inputs. These endpoints supply material and standard-state
heat bookkeeping; they do not determine pH, saturation, stable precipitate,
activity coefficients, solubility, reaction rates, sludge separation, reagent
purity or plant dose. Standard entropy claims can be used with the existing
`standard_equilibrium` helper, retaining the same source and ionic conventions.
The complete direct formation-Gibbs source set also supports
`standard_gibbs_change` and `at_activities`. Every active species needs a
strictly positive dimensionless activity, including water and the named solids.
These diagnostics evaluate the local reaction driving force at the supplied
state; they do not select an equilibrium composition or a finite conversion.

## Compile and check the actual sources

```bash
precipitation_pack_dir="$(mktemp -d)"
for pack in goethite-solid-usgs-b2131-thermochemistry portlandite-solid-usgs-b2131-thermochemistry iron-iii-ion-aqueous-usgs-b2131 water-liquid-usgs-b2131 hydrogen-ion-aqueous-usgs-b2131 calcium-ion-aqueous-usgs-b2131; do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "$precipitation_pack_dir/$pack.fsmatpk"
done
cargo test -p xtask --test mining_iron_precipitation
```

Run each source-backed endpoint with those compiled packs:

```bash
cargo run -p fs-thermochem --example mineral_processing -- \
  --equilibrium --gibbs --reaction ferric-hydrolysis-goethite \
  "$precipitation_pack_dir/iron-iii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$precipitation_pack_dir/water-liquid-usgs-b2131.fsmatpk" \
  "$precipitation_pack_dir/goethite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$precipitation_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" 55.847

cargo run -p fs-thermochem --example mineral_processing -- \
  --equilibrium --gibbs --reaction ferric-goethite-precipitation-with-hydrated-lime \
  "$precipitation_pack_dir/iron-iii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$precipitation_pack_dir/portlandite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$precipitation_pack_dir/goethite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$precipitation_pack_dir/calcium-ion-aqueous-usgs-b2131.fsmatpk" \
  "$precipitation_pack_dir/water-liquid-usgs-b2131.fsmatpk" 55.847

cargo run -p fs-thermochem --example mineral_processing -- \
  --equilibrium --gibbs --reaction hydrated-lime-neutralization \
  "$precipitation_pack_dir/portlandite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$precipitation_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" \
  "$precipitation_pack_dir/calcium-ion-aqueous-usgs-b2131.fsmatpk" \
  "$precipitation_pack_dir/water-liquid-usgs-b2131.fsmatpk" 74.093
```

The last argument is kilograms of the **first declared species actually
consumed**. Thus 55.847 kg in the iron cases means 1000 mol of Fe3+, not
55.847 kg of solution. The lime case uses 1000 mol of pure portlandite.
These commands produce the following reference results:

| Endpoint | Extent, mol | Total standard enthalpy change, MJ | Reaction entropy, J/(mol K) | H/S-derived reaction G, J/mol | Direct reaction G, J/mol |
| --- | ---: | ---: | ---: | ---: | ---: |
| Ferric hydrolysis | 1000 | +58.9 | 200.4 | -849.26 | -900 |
| Combined lime precipitation | 500 | -133.85 | 402.0 | -387556.3 | -391200 |
| Lime neutralization | 1000 | -128.5 | 0.4 | -128619.26 | -129800 |

The combined command consumes 1500 mol portlandite and produces 1000 mol
goethite, 1500 mol Ca2+ and 1000 mol water. Hydrolysis instead produces
3000 mol H+ while consuming 2000 mol water. The neutralization command
consumes 2000 mol H+ and produces 2000 mol water.

`--equilibrium` derives Gibbs quantities from the compatible H/S claims.
`--gibbs` sums the independent directly tabulated formation-Gibbs claims.
The two routes preserve their source values, including their differences;
neither silently replaces the other. Omit both flags to request amounts and
standard heat alone.

For an explicit hydrolysis activity scenario, declare all four activities:

```bash
cargo run -p fs-thermochem --example mineral_processing -- \
  --gibbs --reaction ferric-hydrolysis-goethite \
  --activity 'Fe[3+](aq)=0.001' --activity 'H2O(l)=1' \
  --activity 'FeO(OH)(s)=1' --activity 'H+(aq)=0.001' \
  "$precipitation_pack_dir/iron-iii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$precipitation_pack_dir/water-liquid-usgs-b2131.fsmatpk" \
  "$precipitation_pack_dir/goethite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$precipitation_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" 55.847
```

Here `Q = a(goethite) * a(H+)^3 / (a(Fe3+) * a(water)^2) = 1e-6`,
giving a local reaction Gibbs energy of approximately -35.148 kJ/mol.
Changing only the declared H+ activity to 1 gives Q = 1000 and approximately
+16.224 kJ/mol, reversing the thermodynamic driving direction. These are
supplied relative activities under the retained source conventions, not
solution concentrations or a calculated pH.

The four focused tests compile the real TSVs, check the literal source values
and named-phase admission, verify compatible source heat and direct Gibbs
energy, preserve Hess cancellation at shared activities, and exercise finite
combined and staged inventories. Missing participants or activities and wrong
solid-phase selectors are refused.
