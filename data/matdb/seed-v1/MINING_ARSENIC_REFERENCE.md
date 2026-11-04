# Arsenic acid/base and oxidation reference data

Six matched aqueous NBS profiles add 24 claims for arsenate acid/base
speciation, the first arsenite proton step, and balanced As(III)-to-As(V)
oxidation accounts. Every claim is restricted to 298.15 K and 100000 Pa.

## Literal source rows

| Profile | Species, state `ao` | Molar mass g/mol | Formation H kJ/mol | Formation G kJ/mol | S J/(mol K) |
|---|---|---:|---:|---:|---:|
| `arsenic-acid-aqueous-nbs1989` | H3AsO4 | 141.9432 | -902.5 | -766.0 | 184. |
| `dihydrogen-arsenate-ion-aqueous-nbs1989` | H2AsO4- | 140.9352 | -909.56 | -753.17 | 117. |
| `hydrogen-arsenate-ion-aqueous-nbs1989` | HAsO4(2-) | 139.9272 | -906.34 | -714.60 | -1.7 |
| `arsenate-ion-aqueous-nbs1989` | AsO4(3-) | 138.9192 | -888.14 | -648.41 | -162.8 |
| `arsenous-acid-aqueous-nbs1989` | H3AsO3 | 125.9438 | -742.2 | -639.80 | 195.0 |
| `dihydrogen-arsenite-ion-aqueous-nbs1989` | H2AsO3- | 124.9358 | -714.79 | -587.13 | 110.5 |

All six rows are from Wagman et al., *The NBS Tables of Chemical
Thermodynamic Properties* (1982), Table 20:As, prepared 1964, as digitized by
NIST with the 1989 errata. The evaluation year is preserved; `nbs1989` names
the source convention, not a new measurement date.

- Dataset: <https://doi.org/10.18434/M32124>.
- Corrected CSV: <https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv>.
- Original tables and standard-state conventions, sections 5.5-5.8:
  <https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf>.
- Digitization methods: <https://doi.org/10.6028/jres.125.007>.

The corrected CSV checked on 2026-11-04 has SHA256
`cb8a8b8b94a287e6e8c82eed8bad2cf9b611d574d1cfd95c8528d318dfef4b44`.
Each manifest records that hash. Only attributed numerical facts and original
annotations are redistributed; the source document and table layout are not.

The `ao` standard is hypothetical unit molality with activity coefficients
referenced to infinite dilution. Neutral species use standard molar entropy;
ions use conventional standard molar entropy and the source H+(ao)=0 reference.
Negative ionic entropies are preserved. Heat capacities are unavailable for
these rows, and all uncertainties are explicitly Unstated.

## Identity and component conventions

The four As(V) profiles require `source-arsenic-oxidation-state-v=1`.
The two As(III) profiles require `source-arsenic-oxidation-state-iii=1` and
`source-arsenite-component-plus-water=1` on every claim. The latter retains
the NBS source definitions:

- H3AsO3 is equivalent to HAsO2(ao) + H2O(l).
- H2AsO3- is equivalent to AsO2-(ao) + H2O(l).

The complete H3AsO3 and H2AsO3- formulas must be used consistently in
conservation equations. Replacing either with an unhydrated component changes
the water stoichiometry. This convention is not a new structural or hydration
measurement. The neighboring H3AsO4 `cr` and `aq` rows, and the H3AsO3 `aq`
row, are not the selected standard-state rows.

## Derived reaction checks

With the matched NBS H+ reference, products-minus-reactants sums of direct
formation G give the following dissociation checks. pKa is
`delta_G / (R * 298.15 * ln(10))`, using R = 8.31446261815324 J/(mol K).

| Dissociation | Direct reaction G kJ/mol | Reaction H kJ/mol | Derived pKa |
|---|---:|---:|---:|
| H3AsO4 -> H+ + H2AsO4- | 12.830 | -7.060 | 2.247718752798405 |
| H2AsO4- -> H+ + HAsO4(2-) | 38.570 | 3.220 | 6.757171652021395 |
| HAsO4(2-) -> H+ + AsO4(3-) | 66.190 | 18.200 | 11.595986301459584 |
| H3AsO3 -> H+ + H2AsO3- | 52.670 | 27.410 | 9.227384778635386 |

The displayed pKa digits describe reproducible arithmetic from the tabulated
formation values, not additional experimental precision. As(V) has all three
proton steps here; As(III) has only its first step. A two-state As(III) account
does not establish strongly alkaline speciation with further deprotonation.

The balanced endpoint `H3AsO3 + H2O2(aq) -> H3AsO4 + H2O(l)` combines these
profiles with existing NBS aqueous peroxide and liquid-water profiles. Its
direct standard reaction G is -229.299 kJ/mol and reaction H is -254.960 kJ/mol.
These are standard-state sums; the sign of G does not select conversion.
Finite reagent inventories and an explicit extent are still required.

These profiles supply no reaction-rate law, peroxide selectivity or
decomposition model, nonideal electrolyte correction, temperature extrapolation,
arsenic sorption, scorodite/ferric-arsenate precipitation, fixation yield or
discharge prediction. Compatible solid-phase data and an explicit precipitation
owner remain necessary for fixation claims.

## Run the existing chemistry examples

From the repository root, compile the six new profiles and the four existing
reaction partners into a caller-chosen directory:

```bash
arsenic_pack_dir=target/arsenic-reference-packs
mkdir -p "$arsenic_pack_dir"
for profile in \
  arsenic-acid-aqueous-nbs1989 \
  dihydrogen-arsenate-ion-aqueous-nbs1989 \
  hydrogen-arsenate-ion-aqueous-nbs1989 \
  arsenate-ion-aqueous-nbs1989 \
  arsenous-acid-aqueous-nbs1989 \
  dihydrogen-arsenite-ion-aqueous-nbs1989 \
  hydrogen-ion-aqueous-nbs1989 \
  hydroxide-ion-aqueous-nbs1989 \
  water-liquid-nbs1989-thermochemistry \
  hydrogen-peroxide-aqueous-nbs1989
do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$profile/manifest.tsv" \
    --out "$arsenic_pack_dir/$profile.fsmatpk" || break
done

cargo run -p fs-thermochem --example mineral_speciation -- \
  --family arsenate --materials "$arsenic_pack_dir"

cargo run -p fs-thermochem --example mineral_processing -- \
  --reaction arsenous-acid-peroxide-oxidation --gibbs --equilibrium \
  "$arsenic_pack_dir/arsenous-acid-aqueous-nbs1989.fsmatpk" \
  "$arsenic_pack_dir/hydrogen-peroxide-aqueous-nbs1989.fsmatpk" \
  "$arsenic_pack_dir/arsenic-acid-aqueous-nbs1989.fsmatpk" \
  "$arsenic_pack_dir/water-liquid-nbs1989-thermochemistry.fsmatpk" \
  0.1259438
```

The speciation example prints all four forms for six declared-stock scenarios
on a fixed 10 kg solvent basis, with charge, analytical-pool and mass-action
residuals. The oxidation invocation specifies one mole of reacted H3AsO3 and
prints its stoichiometric requirements; it does not assume that a supplied
feed reaches that extent. Direct formation-Gibbs and caloric H-T*S routes
remain separately visible because rounded source tables need not close
identically between those routes.
