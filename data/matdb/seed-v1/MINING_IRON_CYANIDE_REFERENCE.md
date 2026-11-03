# Iron–cyanide source profiles and finite reference speciation

Four profiles add the original NBS aqueous Fe(II), Fe(III), ferrocyanide and
ferricyanide rows. The runnable `iron_cyanide_reference` consumer joins these
profiles to the existing HCN/CN⁻ and water system. Declared Fe(II) and Fe(III)
pools compete for the **same finite cyanide inventory**, with six cyanides per
complex and independent iron, carbon, nitrogen and charge balances.

**This is an explicitly selected legacy reference calculation. The NBS
iron–cyanide equilibrium constants disagree materially with later critical
literature; the profiles are not a validated capture calibration.** No reference
is selected automatically. Preserve this distinction when using the records in
FrankenMining.

## Literal source records

All four rows are in Table 41:Fe, *Iron (Prepared 1966)*, of Wagman et al.,
*The NBS Tables of Chemical Thermodynamic Properties*, JPCRD 11 Supplement 2
(1982), as digitized by NIST with the 1989 errata. State **ao** denotes an
individual aqueous species; these are not fully ionized salts (state ai),
hydrated solids or total dissolved iron.

| Profile | Species | Original mass (g/mol) | ΔfH° (kJ/mol) | ΔfG° (kJ/mol) | Conventional S° (J/mol K) |
|---|---|---:|---:|---:|---:|
| `iron-ii-ion-aqueous-nbs1989` | Fe²⁺ | 55.8470 | −89.1 | −78.90 | −137.7 |
| `iron-iii-ion-aqueous-nbs1989` | Fe³⁺ | 55.8470 | −48.5 | −4.7 | −315.9 |
| `hexacyanoferrate-ii-ion-aqueous-nbs1989` | Fe(CN)₆⁴⁻ | 211.9544 | 455.6 | 695.08 | 95.0 |
| `hexacyanoferrate-iii-ion-aqueous-nbs1989` | Fe(CN)₆³⁻ | 211.9544 | 561.9 | 729.4 | 270.3 |

Each profile has four claims. Missing heat capacities remain unavailable;
uncertainties remain `Unstated`. Negative ionic entropies retain the source's
conventional H⁺(aq) zero reference. The pre-existing USGS iron profiles retain
their own values and reference identities.

Claims require 298.15 K, 100000 Pa, the NBS elemental/ionic convention and the
hypothetical ideal unit-molality aqueous state referred to infinite dilution.
The model uses a common `n_ref` amount with `a_i = n_i/n_ref` and unit activity
coefficients. A 1 kg solvent basis at the 1 mol/kg standard concentration gives
`n_ref = 1 mol`; process liquid volume is not a substitute for this basis.

Source links:

- [NIST dataset, DOI 10.18434/M32124](https://doi.org/10.18434/M32124)
- [Corrected NIST CSV](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv)
- [Original NBS publication and standard-state definitions](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf)
- [NIST digitization methods](https://doi.org/10.6028/jres.125.007)

The corrected CSV used for transcription has SHA256
`cb8a8b8b94a287e6e8c82eed8bad2cf9b611d574d1cfd95c8528d318dfef4b44`.

## Material source disagreement

For the explicitly selected, matched NBS records, exact conserved reactions give:

| Association | ΔrG° (kJ/mol) | ΔrH° (kJ/mol) | NBS-derived log10 β |
|---|---:|---:|---:|
| Fe²⁺ + 6CN⁻ → Fe(CN)₆⁴⁻ | −260.42 | −358.9 | 45.62361010161814 |
| Fe³⁺ + 6CN⁻ → Fe(CN)₆³⁻ | −300.3 | −293.2 | 52.61028382426822 |

These numbers follow directly from the tabulated formation Gibbs energies and
`ln β = −ΔrG°/(RT)`; they are not new equilibrium measurements. The DOE report
[PNL-5441, Table 2](https://www.osti.gov/servlets/purl/6345492) independently
reproduces the NBS complex rows, including both charge states.

[Keefe and Miller (1996), *Was ferrocyanide a prebiotic reagent?*, pp. 112–114](https://ntrs.nasa.gov/api/citations/19980201273/downloads/19980201273.pdf)
uses thermodynamic Fe(II) β = 2.5 × 10³⁵ at 25 °C and separately calculates
apparent β = 7.1 × 10³⁶ at ionic strength 0.72. Its thermodynamic log10 β
(35.397940) differs from the NBS result by more than ten orders. The distinction
between thermodynamic and apparent constants does not remove this discrepancy.

The consumer prints the discrepancy before reporting its results. It does not
replace any literal Gibbs energy with a formation constant from another
evaluation, attach an unverified uncertainty to reconcile the sources, or claim
that either dataset validates plant capture. A future calibrated model needs
its own internally compatible evaluation and declared activity treatment.

## Run the finite inventory consumer

Compile these nine profiles to `PACK_DIR/PACK_ID.fsmatpk` with the existing
material compiler:

```text
hydrogen-cyanide-aqueous-nbs1989
hydrogen-ion-aqueous-nbs1989
cyanide-ion-aqueous-nbs1989
hydroxide-ion-aqueous-nbs1989
water-liquid-nbs1989-thermochemistry
iron-ii-ion-aqueous-nbs1989
iron-iii-ion-aqueous-nbs1989
hexacyanoferrate-ii-ion-aqueous-nbs1989
hexacyanoferrate-iii-ion-aqueous-nbs1989
```

For each profile:

```sh
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/PACK_ID/manifest.tsv \
  --out PACK_DIR/PACK_ID.fsmatpk
```

Then supply the six explicit input amounts, in order: total cyanide, Fe(II),
Fe(III), Na⁺, Cl⁻, and the amount corresponding to unit activity. All are mol.
The following deliberately declared example has excess alkaline spectator
charge and a 1 kg solvent basis:

```sh
cargo run -p fs-thermochem --example iron_cyanide_reference -- \
  --reference nbs1989-legacy --materials PACK_DIR \
  --amounts 0.001 0.0001 0.00005 0.003 0.00035 1
```

The selected NBS reference gives pH ≈ 11.296462 and free CN⁻ ≈
9.914405686 × 10⁻⁵ mol. Reducing the cyanide stock to 0.0006 mol with Na⁺ =
0.0006 mol, the same Fe stocks and Cl⁻ = 0.00035 mol leaves appreciable
uncomplexed Fe(II); the solver cannot bind 0.0009 mol cyanide from a
0.0006 mol input. These outputs illustrate the selected mathematical
assemblage, not observations or proof of its completeness.

Na⁺ and Cl⁻ are independent caller declarations. The consumer does not
manufacture balancing salts, solve iron mineral accessibility, exchange Fe(II)
with Fe(III), or infer redox rates. Water has fixed unit activity, so this is
not a finite solvent H/O ledger. Hydrolysis, precipitation, protonated
hexacyanoferrates, photochemistry and nonideal activities are absent; omitted
iron hydroxides can matter strongly at alkaline pH. The analytical cyanide
account includes HCN, CN⁻ and the two declared complexes, and is neither a
WAD-CN assay model nor a regulatory total-cyanide method.

The original solver retains the selected source receipts, rejects unsupported
source conditions and preserves each oxidation-state pool separately. The
focused integration test exercises both ligand-rich and ligand-limited
competition against an independent high-precision scalar balance calculation:

```sh
cargo test -p xtask --test mining_iron_cyanide_sources
```
