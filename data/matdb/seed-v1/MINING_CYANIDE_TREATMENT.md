# Peroxide treatment, cyanate hydrolysis and ammonia speciation

Four matched NBS profiles complete a finite free-cyanide treatment path and
the resulting ammonia acid/base family. The existing reaction-network owner
closes prescribed extents against available stocks; the existing aqueous
owner solves ammonia, carbonate and uncomplexed cyanide at one common ideal pH.

## Source facts and identities

The treatment consumes four M/H/G/S claims from each profile at **298.15 K,
100000 Pa**. This table retains the source's displayed units; compilation
converts mass to kg/mol and energy to J/mol.

| Profile | Species | M, g/mol | Hf°, kJ/mol | Gf°, kJ/mol | S°, J/(mol K) |
| --- | --- | ---: | ---: | ---: | ---: |
| `hydrogen-peroxide-aqueous-nbs1989` | H2O2(aq) | 34.0148 | -191.17 | -134.03 | 143.9 |
| `cyanate-ion-aqueous-nbs1989` | OCN[-](aq) | 42.0173 | -146.0 | -97.4 | 106.7 |
| `ammonia-aqueous-nbs1989` | NH3(aq) | 17.0307 | -80.29 | -26.50 | 111.3 |
| `ammonium-ion-aqueous-nbs1989` | NH4[+](aq) | 18.0387 | -132.51 | -79.31 | 113.4 |

All sixteen scalar strings and the `ao` solution states agree with both the
[original NBS tables](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf), printed
pages 2-38, 2-89, 2-64 and 2-65 respectively, and the
[NIST distribution incorporating 1989 errata](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv).
Original table images resolve an OCR error in cyanate's phase marker.
`OCN[-](aq)` identifies cyanate explicitly; a bare CNO formula must not silently
identify the fulminate isomer. Aqueous peroxide and ammonia remain distinct
from pure liquid peroxide and ammonia gas.

Neutral species use `standard-molar-entropy`; ions use
`conventional-standard-molar-entropy` and H+(aq)=0. Aqueous queries require
the NBS unit-activity/infinite-dilution selector. Uncertainty is unstated.
The existing ammonium profile also retains conventional ionic Cp = 79.9
J/(mol K), independently confirmed in the corrected CSV. This point value is
not a temperature law. These data supply no rate, catalyst demand or
extrapolation range.

## Conserved endpoints and thermochemistry

The new `MineralReaction` variants are:

| Variant | Complete integer equation |
| --- | --- |
| `CyanidePeroxideOxidation` | CN- + H2O2(aq) -> OCN- + H2O(l) |
| `CyanateHydrolysisAmmoniumCarbonate` | OCN- + 2 H2O(l) -> NH4+ + CO3(2-) |
| `AmmoniumDissociation` | NH4+ -> NH3(aq) + H+ |
| `HydrogenPeroxideDecomposition` | 2 H2O2(aq) -> 2 H2O(l) + O2(g) |
| `HydrogenPeroxideDecompositionDissolvedOxygen` | 2 H2O2(aq) -> 2 H2O(l) + O2(aq) |

Oxidation and hydrolysis follow the
[EPA treatment review, EPA530-R-94-037, section 3.6](https://archive.epa.gov/epawaste/nonhaz/industrial/special/web/pdf/cyanide.pdf).
The review describes copper-catalyzed treatment; the equations supply no
catalyst loading or rate. Peroxide decomposition follows
[EPA/600/2-90/006, equation 1](https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=20007M9S.TXT).
All equations pass the original exact atom/charge constructor. The gas and
dissolved oxygen endpoints are separate alternatives, not a transfer model.

Actual compatible NBS packs give these values per mole of complete equation:

| Endpoint | ΔH°, J/mol reaction | Direct ΔG°, J/mol reaction | ΔS°, J/(mol reaction K) |
| --- | ---: | ---: | ---: |
| Cyanide oxidation | -391260 | -372899 | -61.39 |
| Cyanate hydrolysis | -91990 | -35462 | -190.02 |
| Ammonium dissociation | 52220 | 52810 | -2.1 |
| Peroxide decomposition, O2(g) | -189320 | -206198 | 57.158 |
| Peroxide decomposition, O2(aq) | -201020 | -189798 | -37.08 |

Direct formation-Gibbs sums and H−TS remain independent. Rounded historical
source differences are reported, not reconciled by changing data.

## Compile the source packs

The two runnable modes together use thirteen profiles:

```bash
mkdir -p /tmp/frankensim-cyanide-treatment-packs
for pack in \
  hydrogen-peroxide-aqueous-nbs1989 \
  cyanate-ion-aqueous-nbs1989 \
  ammonia-aqueous-nbs1989 \
  ammonium-ion-aqueous-nbs1989 \
  cyanide-ion-aqueous-nbs1989 \
  hydrogen-cyanide-aqueous-nbs1989 \
  hydrogen-ion-aqueous-nbs1989 \
  hydroxide-ion-aqueous-nbs1989 \
  carbonate-ion-aqueous-nbs1989 \
  bicarbonate-ion-aqueous-nbs1989 \
  carbonic-component-aqueous-nbs1989 \
  oxygen-gas-nbs1989-thermochemistry \
  water-liquid-nbs1989-thermochemistry
do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "/tmp/frankensim-cyanide-treatment-packs/$pack.fsmatpk"
done
```

The dissolved-oxygen decomposition API additionally uses the existing
`oxygen-aqueous-nbs1989-thermochemistry` profile. Integration tests exercise
that endpoint as well.

## Finite treatment with shared peroxide

```bash
cargo run -p fs-thermochem --example cyanide_treatment -- \
  treat /tmp/frankensim-cyanide-treatment-packs 10 8 100 6 1 6
```

Arguments after the directory are initial CN-, aqueous H2O2 and liquid-water
amounts, then oxidation, decomposition and hydrolysis extents, all in mol.
All initial products are zero. Stage one requests oxidation extent 6 and
decomposition extent 1 against the **same 8 mol peroxide stock**. Only after
that stage succeeds can its cyanate and water fund hydrolysis of extent 6.

| Species | Initial, mol | After oxidation/decomposition, mol | After hydrolysis, mol |
| --- | ---: | ---: | ---: |
| CN- | 10 | 4 | 4 |
| H2O2(aq) | 8 | 0 | 0 |
| H2O(l) | 100 | 108 | 96 |
| OCN- | 0 | 6 | 0 |
| NH4+ | 0 | 0 | 6 |
| CO3(2-) | 0 | 0 | 6 |
| O2(g) | 0 | 1 | 1 |

Standard chemical enthalpy change is **−3,088,820 J** and the source-mass
residual is zero in this example. Output includes stage extents, both Gibbs
paths, final species masses and original pack identities/hashes. These are
standard-state chemical accounts, not operating cooling duty.

Oxidation extent 7 plus decomposition extent 1 would consume 9 mol peroxide
and refuses. Excess cyanate or water demand in stage two also refuses. Both
stages and all diagnostics must succeed before accepted output is printed.
The network cannot borrow another reaction's products within the same stage.

Nitrogen remains in cyanide, cyanate or ammonium. Prescribed extents establish
no treatment time, efficiency, WAD/total-cyanide removal, nitrogen removal or
water-quality limit. Complexed cyanide is not relabelled as free CN-.

## All three aqueous families at one pH

```bash
cargo run -p fs-thermochem --example cyanide_treatment -- \
  speciate /tmp/frankensim-cyanide-treatment-packs 6 6 4 1000 10 0
```

Arguments are analytical ammonia, inorganic carbon and uncomplexed cyanide
totals, common amount at unit activity `n_ref`, and explicit sodium/chloride
spectator stocks. Counterions are not inferred. Original source-backed
dissociations and water solve all families and electroneutrality together.

This supplied ideal composition gives pH **9.9164202612**, NH4+ =
1.0678853321 mol, NH3(aq) = 4.9321146679 mol, HCN = 0.6863255008 mol and
CN- = 3.3136744992 mol. Ammonia's constant comes from selected formation-Gibbs
sources, not a hard-coded pKa.

This separate analytical-pool calculation uses `a=n/n_ref` and unit water
activity. Under the NBS unit-molality convention, `n_ref` is solvent kg times
1 mol/kg, **not solution volume**. It is not silently derived from the treatment
example's finite water. Formal H2CO3 retains the NBS CO2(aq)+H2O(l) convention
and its explicit selector. Cyanate/peroxide are outside these acid families.
Finite-solvent H/O balance, nonideal activities, ammonia/HCN gas loss,
metal complexes, precipitation and kinetics require their respective owners.

`xtask/tests/mining_cyanide_treatment_sources.rs` exercises actual packs,
literal claims/receipts, five H/G/S sums, finite reagent and stage limits,
source-mass closure and coupled family/charge/mass-action balances.
