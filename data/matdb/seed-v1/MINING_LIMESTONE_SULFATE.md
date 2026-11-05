# Finite limestone dissolution and gypsum formation

`fs-thermochem` can now equilibrate selected salts that compete for one finite
cation inventory. The working calcite/gypsum consumer connects limestone acid
consumption to sulfate precipitation: dissolving calcite changes free calcium,
carbonate speciation and pH, which changes how much gypsum can form. Both solids
use the same final aqueous state and the same original calcium stock.

The reusable owner is
`stoichiometry::thermodynamics::acid_base::precipitation::salt::assemblage::CommonCationSaltSystem`.
The executable is
[`limestone_sulfate_equilibrium.rs`](../../../crates/fs-thermochem/examples/limestone_sulfate_equilibrium.rs).
The source compiler consumer test is
[`mining_limestone_sulfate_sources.rs`](../../../xtask/tests/mining_limestone_sulfate_sources.rs).

## Matched source data

The new `calcite-solid-nbs1989-thermochemistry` profile supplies the missing
calcite standard state from the same NBS evaluation as the existing calcium,
carbonate, sulfate, gypsum and water profiles. The independently evaluated
`calcite-solid-usgs-b2131-thermochemistry` profile remains a different source;
its formation values cannot be substituted into this model by changing a label.

The original numerical facts are the `CaCO3`, `calcite`, `cr` row in Table 94:Ca
(prepared 1968) of the
[NIST corrected NBS CSV](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv),
distributed with the 1989 errata. Its adjacent aragonite row is a different phase.
The original tables and conventions are in
[Wagman et al., 1982](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf);
[NIST's digitization methods](https://doi.org/10.6028/jres.125.007) explain the
column identities and application of errata. The corrected CSV retrieved for
this profile had SHA-256
`cb8a8b8b94a287e6e8c82eed8bad2cf9b611d574d1cfd95c8528d318dfef4b44`.

| Calcite property | Original value |
| --- | ---: |
| Molar mass | 100.0894 g/mol |
| Standard molar formation enthalpy | -1206.92 kJ/mol |
| Standard molar formation Gibbs energy | -1128.79 kJ/mol |
| Standard molar entropy | 92.9 J/(mol K) |
| Standard molar isobaric heat capacity | 81.88 J/(mol K) |

Every claim is restricted to 298.15 K, 100000 Pa, the NBS elemental reference,
pure crystalline calcite and the engineering reference model. Uncertainty is
`Unstated`, as in the source. The heat capacity is a single reference-state
value; it supplies no temperature law. The profile supplies neither limestone
purity nor ore texture, dissolution kinetics or a specimen qualification.

The example loads 11 immutable normalized packs, verifies their contents and
identities, and queries the exact phase/reference axes. Individual ions retain
the NBS `H+(ao)=0` convention. `H2CO3(aq)` retains the original formal carbonic
component convention (`CO2(aq) + H2O`); it is not a newly measured molecular
carbonic-acid concentration. The original receipts stay attached to the water,
acid dissociation and mineral reaction owners.

## Equations and finite inventories

The selected precipitation reactions are exactly conserved:

```
Ca[2+](aq) + CO3[2-](aq)                  -> CaCO3(s)
Ca[2+](aq) + SO4[2-](aq) + 2 H2O(l)      -> CaSO4.2H2O(s)
```

The same aqueous owner solves water autoprotolysis, the two carbonic-component
dissociations and hydrogen-sulfate dissociation with full electroneutrality.
Source-derived precipitation Gibbs changes are -47400 J/mol for calcite and
-24912 J/mol for gypsum. For the dissolution activity products, the selected
ideal model uses `Ksp = exp(delta_G_precipitation / (R*T))`.

The calcite reaction's independent enthalpy and entropy columns give
`delta_H = +13050 J/mol` and `delta_S = +202.9 J/(mol K)`, hence
`delta_H - T*delta_S = -47444.635 J/mol`. The direct Gibbs value remains
-47400 J/mol. The implementation does not reconcile this source discrepancy
or treat the numerical stopping tolerance as source uncertainty.

Inputs are analytical calcium, inorganic-carbon and sulfate amounts **including
their contents in any supplied solids**, fixed spectator-ion amounts and a
common amount at unit aqueous activity. Supply each component once. The
selected ideal approximation uses `a_i = n_i / n_ref`, with unit pure-solid and
solvent activities. For an ideal molality interpretation, 1 mol of activity
basis corresponds to 1 kg of solvent at the 1 mol/kg standard activity scale;
an amount basis alone is not an inferred slurry volume or solvent inventory.

The final independent component accounts are:

```
Ca_total = free_Ca + calcite + gypsum
C_total  = H2CO3 + HCO3[-] + CO3[2-] + calcite
S_total  = HSO4[-] + SO4[2-] + gypsum
```

Amounts in these equations are moles of the named component or formula unit.
Water in the hydrate has its original exact stoichiometry but fixed activity;
this model does not close a finite solvent H/O account.

For every selected phase, a positive final solid amount requires its dissolution
`ln(IAP/K)` to be zero within the declared tolerance. An absent phase requires
`ln(IAP/K) <= 0`; an empty original component has no defined saturation ratio.
This finite-stock distinction follows the usual pure-phase equilibrium
semantics described in the
[USGS PHREEQC equilibrium-phase documentation](https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-13.htm).
No PHREEQC constants, solver or activity model are embedded here.

## Reusing the existing owners

`CommonCationSaltSystem::new` accepts 1 through 8 original `SaltPrecipitation`
owners. Every phase must share the same cation, charge, complete aqueous
reaction set and exact selected formation-Gibbs receipts. Each selected solid
must consume a distinct acid family. Overlapping families and polymorph
competition are refused because their phase amounts need another phase-selection
model. Different integer cation/anion stoichiometries retain their original
formula-unit counts through the existing single-salt owner.

The system uses a canonical solid-identity order. Each coordinate invokes the
existing bounded single-salt solve while holding all other solid stocks fixed.
Only the dissolved cation plus that coordinate's current solid stock are
available to it. The original pH and every acid family are solved again for
each trial. Dissolved amounts are carried directly, preserving positive traces
even when a solid stock rounds to the original total.

At the end of a complete pass, the system recomputes **all** saturation residuals
at the one final aqueous state. It includes the inner pH resolution contribution
in every phase decision. Previously converged coordinate residuals are not
accepted as a joint solution. Independently summed cation and complete family
accounts must also meet their declared absolute and relative allowances.

`SaltAssemblageOptions` exposes the original phase/aqueous controls plus:

| Control | Default | Hard admission bound |
| --- | ---: | ---: |
| Complete phase sweeps | 128 | 1 through 256 |
| Total aqueous solves, including nested trials | 16384 | 1 through 65536 |
| Absolute component-account tolerance | 1e-14 mol | Finite, positive |
| Relative component-account tolerance | 1e-11 | Finite, nonnegative |

Before every coordinate call, the remaining global aqueous budget limits that
call's possible endpoint and iteration solves. Inner pH iterations retain the
existing owner bound. Exhausted work, unresolvable phase boundaries, cancellation
and numerical range loss return no result; the caller's stocks remain unchanged.
The returned `SaltAssemblageSpeciation` records all component accounts, shared
aqueous state, phase diagnostics, original owners, numerical options and actual
work counts. It is a numerical result under the selected model, not an interval
certificate or a propagated source-uncertainty envelope.

## Run the source-backed example

From the repository root, compile the 11 source profiles into one directory:

```bash
mkdir -p /tmp/limestone-sulfate-packs
for name in \
  calcite-solid-nbs1989-thermochemistry \
  gypsum-selenite-solid-nbs1989-thermochemistry \
  calcium-ion-aqueous-nbs1989 \
  carbonic-component-aqueous-nbs1989 \
  bicarbonate-ion-aqueous-nbs1989 \
  carbonate-ion-aqueous-nbs1989 \
  hydrogen-sulfate-ion-aqueous-nbs1989 \
  sulfate-ion-aqueous-nbs1989 \
  water-liquid-nbs1989-thermochemistry \
  hydrogen-ion-aqueous-nbs1989 \
  hydroxide-ion-aqueous-nbs1989
do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$name/manifest.tsv" \
    --out "/tmp/limestone-sulfate-packs/$name.fsmatpk"
done
cargo run -p fs-thermochem --example limestone_sulfate_equilibrium -- \
  /tmp/limestone-sulfate-packs
```

These are the normal workspace commands; they require the repository's declared
toolchain and sibling dependencies. The executable also accepts six explicit
amounts after the pack directory, in this order: calcium total, carbon total,
sulfate total, sodium, chloride and the unit-activity amount, all in mol.

The default demonstration uses authored inputs, not a plant calibration:
0.05 mol calcium, 0.03 mol inorganic carbon, 0.05 mol sulfate, 0.10 mol sodium,
0.06 mol chloride and a 1 mol activity basis. One possible bookkeeping origin is
0.03 mol calcite, 0.02 mol calcium chloride, 0.05 mol sodium sulfate and 0.02 mol
HCl. The second case adds 0.005 mol HCl, represented by the additional chloride
stock and a newly solved proton balance. Both requested cases must solve before
the executable emits a successful report.

| Ideal-model output | Initial | After additional 0.005 mol HCl |
| --- | ---: | ---: |
| pH | 6.6436764 | 6.5367336 |
| Dissolved free calcium, mol | 0.00246011315 | 0.00287659861 |
| Calcite, mol | 0.0151054091 | 0.0121458109 |
| Gypsum, mol | 0.0324344778 | 0.0349775905 |
| Dissolved inorganic carbon, mol | 0.0148945909 | 0.0178541891 |
| Dissolved sulfate family, mol | 0.0175655222 | 0.0150224095 |

Both minerals are present in these two modeled equilibria. The additional
calcite dissolution supplies calcium for both the increase in gypsum and the
increase in free calcium; it cannot supply those destinations twice. The
executable prints every aqueous form, all three component-account residuals,
charge residual, final phase saturation diagnostics and work counts.

## Checks and practical limits

The direct owner tests include a closed-form two-phase pH/stock solution,
bitwise reduction to the original single-salt owner, phase absence and empty
stocks, extensive scaling, canonical ordering, preservation of a positive
dissolved trace, incompatible receipts and ambiguous phase-family rejection,
and atomic refusal under insufficient work or cancellation.

The compiler-backed integration tests build the original 11 profiles, verify
the five literal calcite facts and receipts, reject wrong source domains, and
check the direct formation-property sums. An independent known answer chooses
pH 8, free-calcium activity 0.005, 0.02 mol calcite and 0.03 mol gypsum, then uses
literal NBS energy differences and acid-speciation polynomials to construct the
complete input stocks. The owner must recover that state. A separate acid-addition
test independently sums all Ca/C/S and charge accounts and final activity
products using the actual consumer.

This supplies a thermodynamic material connection for reagent accounting and
scaling scenarios. It does not predict how quickly limestone reacts or gypsum
coats its surface. The selected model excludes nonideal activity coefficients,
calcium ion pairing, metal complexes, other carbonate/sulfate solids, atmospheric
CO2 exchange, gas loss, redox, thermal evolution and metastable phase persistence.
Other dissolved species, ore purity, particle area, residence time, mixing,
passivation, recovery and operating cost remain explicit inputs or future
physical owners. Do not infer a plant acid dose or treatment guarantee from
these reference-state results.
