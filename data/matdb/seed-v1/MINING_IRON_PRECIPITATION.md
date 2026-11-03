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
enthalpy; S is standard molar entropy of the named pure solid.

| Profile | Formula | Molar mass, kg/mol | H, J/mol | S, J/(mol K) | Reported H magnitude, J/mol | Reported S magnitude, J/(mol K) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `goethite-solid-usgs-b2131-thermochemistry` | FeO(OH) | 0.088854 | -562600 | 60.4 | 2100 | 0.6 |
| `portlandite-solid-usgs-b2131-thermochemistry` | Ca(OH)2 | 0.074093 | -986100 | 83.4 | 1300 | 0.4 |

Each pack has five claims: molar mass, formation enthalpy, entropy and the two
reported uncertainty magnitudes. Every claim requires the named phase selector
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

The compatible participants already in the database are
`iron-iii-ion-aqueous-usgs-b2131`, `hydrogen-ion-aqueous-usgs-b2131`,
`calcium-ion-aqueous-usgs-b2131` and `water-liquid-usgs-b2131`.
Their formation enthalpies are respectively -49900, 0, -543000 and
-285800 J/mol. The ions use the source's unit-activity, infinite-dilution
aqueous reference and explicit H+(aq)-zero convention; their entropy key is
`conventional-standard-molar-entropy`. Water uses the pure liquid phase.

## Conserved reaction endpoints and compatible heat

The hydrolysis equation is printed in
[USGS Open-File Report 2005-1148](https://pubs.usgs.gov/of/2005/1148/of2005-1148.pdf)
and is also the reverse of the goethite dissolution reaction in the
[USGS PHREEQC database](https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/html/final-90.html).
The equation supplies the reaction identity; all heats below use the selected
Bulletin 2131 rows, without mixing in PHREEQC's separate numeric evaluation.

```text
FerricHydrolysisGoethite:
Fe3+(aq) + 2 H2O(l) -> FeO(OH)(goethite) + 3 H+(aq)
Delta H = +58900 J/mol of the written equation

Existing HydratedLimeNeutralization, now with compatible USGS portlandite:
Ca(OH)2(portlandite) + 2 H+(aq) -> Ca2+(aq) + 2 H2O(l)
Delta H = -128500 J/mol of the written equation

FerricGoethitePrecipitationWithHydratedLime:
2 Fe3+(aq) + 3 Ca(OH)2(portlandite)
    -> 2 FeO(OH)(goethite) + 3 Ca2+(aq) + 2 H2O(l)
Delta H = -267700 J/mol of the written equation
```

Each variant builds through the exact atom/charge conservation checker. Two
hydrolysis extents plus three lime-neutralization extents cancel the intermediate
protons and four of the produced water molecules, giving the combined equation.
The source heat obeys the same Hess relation:
`2 * 58900 + 3 * (-128500) = -267700 J/mol`.
Original participant queries and verified receipts remain available in the
reaction-enthalpy result.

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
  --equilibrium --reaction ferric-hydrolysis-goethite \
  "$precipitation_pack_dir/iron-iii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$precipitation_pack_dir/water-liquid-usgs-b2131.fsmatpk" \
  "$precipitation_pack_dir/goethite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$precipitation_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" 55.847

cargo run -p fs-thermochem --example mineral_processing -- \
  --equilibrium --reaction ferric-goethite-precipitation-with-hydrated-lime \
  "$precipitation_pack_dir/iron-iii-ion-aqueous-usgs-b2131.fsmatpk" \
  "$precipitation_pack_dir/portlandite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$precipitation_pack_dir/goethite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$precipitation_pack_dir/calcium-ion-aqueous-usgs-b2131.fsmatpk" \
  "$precipitation_pack_dir/water-liquid-usgs-b2131.fsmatpk" 55.847

cargo run -p fs-thermochem --example mineral_processing -- \
  --equilibrium --reaction hydrated-lime-neutralization \
  "$precipitation_pack_dir/portlandite-solid-usgs-b2131-thermochemistry.fsmatpk" \
  "$precipitation_pack_dir/hydrogen-ion-aqueous-usgs-b2131.fsmatpk" \
  "$precipitation_pack_dir/calcium-ion-aqueous-usgs-b2131.fsmatpk" \
  "$precipitation_pack_dir/water-liquid-usgs-b2131.fsmatpk" 74.093
```

The last argument is kilograms of the **first declared species actually
consumed**. Thus 55.847 kg in the iron cases means 1000 mol of Fe3+, not
55.847 kg of solution. The lime case uses 1000 mol of pure portlandite.
These commands produce the following reference results:

| Endpoint | Extent, mol | Total standard enthalpy change, MJ | Reaction entropy, J/(mol K) | Reaction Gibbs energy, J/mol |
| --- | ---: | ---: | ---: | ---: |
| Ferric hydrolysis | 1000 | +58.9 | 200.4 | -849.26 |
| Combined lime precipitation | 500 | -133.85 | 402.0 | -387556.3 |
| Lime neutralization | 1000 | -128.5 | 0.4 | -128619.26 |

The combined command consumes 1500 mol portlandite and produces 1000 mol
goethite, 1500 mol Ca2+ and 1000 mol water. Hydrolysis instead produces
3000 mol H+ while consuming 2000 mol water. The neutralization command
consumes 2000 mol H+ and produces 2000 mol water.

`--equilibrium` derives Gibbs quantities from the compatible H/S claims.
Omit that flag to request amounts and standard heat alone. The new solid
profiles contain no separate formation-Gibbs claims, so the independent
`--gibbs` route is unavailable with this particular source set.

The three focused tests compile the real TSVs, check the literal source
values and named-phase admission, verify compatible source heat and exact
Hess cancellation, and exercise finite combined and staged inventories.
