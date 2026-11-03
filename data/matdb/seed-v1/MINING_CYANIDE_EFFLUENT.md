# Cyanide oxidation and retained effluent nitrogen

Five matched NBS profiles supply the peroxide, cyanate and ammonia chemistry
needed to keep reagent demand and downstream nitrogen in an explicit mineral
processing account. They connect directly to the existing cyanide, water,
hydrogen-ion and carbonate profiles. The data are immutable source claims;
`ConservedReaction` and `IdealAcidBaseSystem` execute the calculations.

## Literal source rows

All rows are aqueous `ao` at 298.15 K and 100000 Pa. Values below retain
the printed signs and precision. The source labels the cyanate formula `OCN-`.

| Profile | Species | Molar mass g/mol | Formation H kJ/mol | Formation G kJ/mol | S J/(mol K) | Cp J/(mol K) |
|---|---|---:|---:|---:|---:|---:|
| `hydrogen-peroxide-aqueous-nbs1989` | H2O2 | 34.0148 | -191.17 | -134.03 | 143.9 | unavailable |
| `hydroperoxide-ion-aqueous-nbs1989` | HO2- | 33.0068 | -160.33 | -67.3 | 23.8 | unavailable |
| `cyanate-ion-aqueous-nbs1989` | OCN- | 42.0173 | -146.0 | -97.4 | 106.7 | unavailable |
| `ammonia-aqueous-nbs1989` | NH3 | 17.0307 | -80.29 | -26.50 | 111.3 | unavailable |
| `ammonium-ion-aqueous-nbs1989` | NH4+ | 18.0387 | -132.51 | -79.31 | 113.4 | 79.9 |

The 21 claims come from Wagman et al., *The NBS Tables of Chemical
Thermodynamic Properties* (1982), as digitized by NIST with the 1989 errata:

- Table 2:H, prepared 1963: H2O2 and HO2-.
- Table 18:N, prepared 1964: NH3 and NH4+.
- Table 23:C, prepared 1966: OCN-, source name `cyanate`.
- Dataset: <https://doi.org/10.18434/M32124>.
- Corrected CSV: <https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv>.
- Original tables and conventions: <https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf>,
  sections 5.5–5.8.
- Digitization methods: <https://doi.org/10.6028/jres.125.007>.

The checked corrected CSV has SHA256
`cb8a8b8b94a287e6e8c82eed8bad2cf9b611d574d1cfd95c8528d318dfef4b44`.
Only attributed numerical facts and original annotations are redistributed.

The aqueous standard is hypothetical unit molality with activity coefficients
referenced to infinite dilution. Individual-ion formation properties, entropy
and heat capacity retain the source H+(ao)=0 convention. The neutral H2O2 and
NH3 entropies are molecular standard entropies. Every query retains phase,
reference, temperature and pressure; uncertainty is explicitly Unstated.

## Executable conserved pathways

The integration test `xtask/tests/mining_cyanide_effluent_sources.rs` builds
these equations from explicit C/H/N/O atom counts, integer coefficients and
charges through the public `ConservedReaction` API. Species names do not infer
formula or charge.

| Conserved equation | Direct standard reaction G kJ/mol | Standard reaction H kJ/mol |
|---|---:|---:|
| CN- + H2O2(aq) -> OCN- + H2O(l) | -372.899 | -391.260 |
| OCN- + 2 H2O(l) -> NH3(aq) + HCO3- | -41.612 | -54.620 |
| OCN- + 2 H2O(l) -> NH4+ + CO3(2-) | -35.462 | -91.990 |

These numbers are products-minus-reactants sums of the directly tabulated
NBS formation claims. They describe standard reaction properties, not the
finite bath's free energy or temperature change. The two hydrolysis product
choices differ by the explicitly balanced proton transfer
`NH4+ + CO3(2-) -> NH3 + HCO3-`, whose direct standard G is -6.150 kJ/mol.
They must not be combined as two independent sinks consuming the same cyanate.

`limiting_extent` and `advance_inventory` enforce finite reagents, explicit
zero products and the original stoichiometry. An illustrative account starting
with 1 mol free CN-, 0.25 mol aqueous H2O2 and 55.5 mol water can oxidize at most
0.25 mol cyanide along the specified equation. Advancing the hydrolysis equation
by the same chosen extent leaves 0.75 mol CN-, produces 0.25 mol NH3 and
0.25 mol bicarbonate, and retains 55.25 mol water. Carbon, nitrogen, hydrogen,
oxygen, charge and source molar-mass totals close. Nitrogen has moved to ammonia;
it has not disappeared. The standard enthalpy sum at that extent is -111.470 kJ,
excluding mixing, dilution, heat transfer and secondary chemistry.

No reaction extent is selected from the sign of standard G alone. The test
uses an explicitly chosen extent within the finite-inventory limit. Actual
conversion requires a supported kinetic or equilibrium owner and the
appropriate activities.

## Ammonia and peroxide speciation

The same profiles construct `AcidDissociation::from_sources` and
`AcidFamily` owners without hard-coded textbook pKa overrides:

| Dissociation | pKa from direct NBS G at 298.15 K |
|---|---:|
| H2O2(aq) -> H+ + HO2- | 11.69059020843629 |
| NH4+ -> H+ + NH3(aq) | 9.2519117174812 |

The check solves both families and water at a shared pH and explicit spectator
charge. Added sodium countercharge raises the solved pH and increases both
the free-ammonia fraction and the hydroperoxide fraction. Analytical ammonia
is NH3+NH4+; aqueous NH3 alone is not total ammonia. In this acid/base owner the
water activity is fixed at one and the caller supplies the common solvent
amount at unit activity; that is a different assumption from the finite-water
reaction account above.

## Scope and verification

This addition supplies source-backed thermodynamics, exact reaction definitions,
finite reagent accounting and ideal acid/base integration. It supplies no
peroxide reaction-rate fit, catalytic performance, competing peroxide
decomposition model, copper-cyanide destruction law, ammonia stripping model,
nonideal electrolyte model or treatment/discharge prediction. The lone NH4+
heat-capacity value is constrained to 298.15 K and does not authorize temperature
extrapolation. Missing heat capacities remain unavailable.

The source/integration target compiles all five new packs and their matched
existing participants, verifies source receipts and original values, checks
finite reagent and nitrogen balances, and refuses missing sources, unsupported
reference conditions and unbalanced equations:

```sh
cargo test -p xtask --test mining_cyanide_effluent_sources
```

To compile a single new profile into a runtime artifact:

```sh
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/cyanate-ion-aqueous-nbs1989/manifest.tsv \
  --out /tmp/cyanate-ion-aqueous-nbs1989.fsmatpk
```

No additional data store or parsing path is introduced; the existing immutable
pack loader, property queries, conserved reactions and acid/base solver consume
the new source rows.
