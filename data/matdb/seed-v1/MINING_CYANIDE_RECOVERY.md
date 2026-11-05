# Copper sulfide precipitation and cyanide recovery

Three matched NBS source profiles complete a working selected-phase calculation
for copper-bearing cyanide liquor. `CopperSulfidePrecipitation` couples one
Cu2S phase to the existing `IdealLigandSystem`: finite copper, sulfide and
cyanide inventories, every declared copper complex, HCN/CN-, H2S/HS- and charge
balance are solved together. Acidification can release cyanide from copper
complexes; the runnable consumer explicitly separates the calculated solid
before adding base to the remaining liquor.

This is the chemical connection used by the SART recovery route described in
[Ford, Henderson and Fleming, SGS/Kinross Technical Paper 2008-51, pp. 4-5](https://www.sgs.com/-/media/sgscorp/documents/corporate/brochures/sgs-min-tp-2008-51-application-of-sart-to-heap-leaching-of-gold-copper.cdn.en-CA.pdf).
The implementation supplies a selected ideal equilibrium model, not the paper's
pilot-plant performance, rates, separation efficiency or operating calibration.

## Source identity and conditions

The new profiles retain the original NBS rows at **298.15 K and 100000 Pa**:

| Profile | Source identity | M (g/mol) | Delta-f H (kJ/mol) | Delta-f G (kJ/mol) | S (J/mol K) | Cp (J/mol K) |
|---|---|---:|---:|---:|---:|---:|
| `hydrogen-sulfide-aqueous-nbs1989` | Table 14:S, H2S, ao | 34.0800 | -39.7 | -27.83 | 121. | missing |
| `hydrosulfide-ion-aqueous-nbs1989` | Table 14:S, HS-, ao | 33.0720 | -17.6 | 12.08 | 62.8 | missing |
| `copper-i-sulfide-alpha-solid-nbs1989` | Table 36:Cu, Cu2S, alpha, cr | 159.1440 | -79.5 | -86.2 | 120.9 | 76.32 |

These are **13 source facts**. All uncertainties remain `Unstated`; absent heat
capacities remain unavailable. HS- entropy retains the individual-ion,
H+(aq)-zero convention through `conventional-standard-molar-entropy`.
The neutral H2S row is dissolved molecular H2S, not the neighboring gas row.
The Cu2S formula unit contains two Cu atoms and one S atom. The source's alpha
label is retained as `Cu2S(s-nbs-alpha)` and
`source-phase-copper-i-sulfide-nbs-alpha`.

The existing USGS chalcocite profile has a different evaluation and source
convention; it is not silently substituted into this NBS calculation. Neither
record represents a mixed commercial precipitate or its particle properties.

Primary numerical sources:

- [NIST dataset DOI 10.18434/M32124](https://doi.org/10.18434/M32124).
- [NIST corrected CSV](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv).
- [Wagman et al., NBS Tables, JPCRD 11 Supplement 2 (1982)](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf), Tables 14:S and 36:Cu; source-state conventions in sections 5.5-5.8.
- [NIST digitization and errata methods](https://doi.org/10.6028/jres.125.007).

The transcribed corrected CSV has SHA256
`cb8a8b8b94a287e6e8c82eed8bad2cf9b611d574d1cfd95c8528d318dfef4b44`.

## Exact equations and source-derived constants

`MineralReaction` adds these exact atom/charge-conserved equations:

```text
HydrogenSulfideDissociation:
  H2S(aq) -> H+(aq) + HS-(aq)

CopperSulfideHydrosulfidePrecipitation:
  2 Cu+(aq) + HS-(aq) -> Cu2S(s-nbs-alpha) + H+(aq)

CopperCyanideHydrosulfideRecovery:
  2 Cu(CN)3(2-)(aq) + HS-(aq) + 5 H+(aq)
    -> Cu2S(s-nbs-alpha) + 6 HCN(aq)
```

The last equation is the hydrosulfide form of the SGS sulfide endpoint. It
follows from the existing copper/cyanide associations, HCN dissociation and
the new phase equation; it does not assume that all dissolved copper starts as
the tricyano complex or that reaction proceeds completely.

| Equation | Direct source Delta-r G (J/mol equation) | log10 K |
|---|---:|---:|
| H2S first dissociation | 39910 | -6.99192949526 |
| Cu2S precipitation from Cu+ and HS- | -198240 | 34.73014540567 |
| Tricyanocuprate net recovery | -187680 | 32.88011344702 |

The phase's formation-enthalpy sum is -205240 J/mol equation. **The net SART
heat remains unavailable** because the selected NBS copper-cyanide complex
rows lack formation enthalpies. Direct Gibbs sums and enthalpy-derived
quantities are not substituted for one another.

For the precipitation-direction constant, the solved dissolution saturation
residual is

```text
r = 2 ln(aCu+) + ln(aHS-) - ln(aH+) + ln(K_precip).
```

A positive solid requires `r = 0` within the declared numerical tolerance.
An absent solid requires resolved `r <= 0`. Every aqueous activity uses the
same caller-supplied `n_ref`, with `a_i = n_i/n_ref`. For the source's ideal
unit-molality reference, 10 kg solvent gives `n_ref = 10 mol`; liquid volume
is not an interchangeable basis. Solid and solvent activities are one.

## Finite partition and admitted assemblage

The public owner is
`fs_thermochem::stoichiometry::thermodynamics::acid_base::complexation::precipitation::CopperSulfidePrecipitation`.
Its constructor binds the original ligand model, named copper and sulfide
families, exact phase equation and source selections. Shared species must
retain identical original formation-Gibbs receipts. Inputs supply every
original metal/acid total once and explicitly declare spectator counterions.

The Cu and sulfide totals include any initial solid at two Cu and one S per
formula unit. Cyanide stays in the aqueous/complex inventory; Cu2S contains no
cyanide. Each outer trial invokes the original ligand and acid/base solver,
so copper release cannot reuse a stale free-CN amount or pH.

The phase coordinate is the remaining limiting dissolved formula-equivalent
amount, `min(Cu_dissolved/2,S_dissolved)`, with any stoichiometric excess kept
separately. Small dissolved traces survive even when subtracting twice the
solid from total copper would round to zero. Results retain independent
Cu/S/CN balances, source owners, the phase bracket and the inner numerical
spread. These brackets describe solver resolution, not experimental
uncertainty or rigorous interval certificates. A caller-declared trace floor,
work budget or cancellation can refuse a calculation; no partial state is
returned as success.

**The entire final pH bracket must lie in [2,10].** This is a conservative
selected-model domain, not a measured accuracy interval. The sulfide family
contains only H2S and HS-. The historical NBS S(2-) row would imply a second
pKa near 12.91, materially different from later investigations; adding it to
fit a salt API would introduce a poorly supported high-pH prediction. See
[Migdisov et al. (2002)](https://doi.org/10.1016/S0016-7037(01)00896-1)
and [Schoonen and Barnes (1988)](https://doi.org/10.1016/0016-7037(88)90326-2).

The selected assemblage excludes S(2-), polysulfides, copper hydrosulfide
complexes, CuCN and competing solids, Cu(I)/Cu(II) redox, sulfide oxidation,
gas escape and nonideal activity coefficients. It does not establish
nucleation, residence time, particle settling or filtration. Water remains a
fixed solvent reference, so no complete H/O solvent inventory is claimed.

## Compile and run the three-stage consumer

Compile these twelve original profiles to `PACK_DIR/PACK_ID.fsmatpk`:

```text
hydrogen-sulfide-aqueous-nbs1989
hydrosulfide-ion-aqueous-nbs1989
copper-i-sulfide-alpha-solid-nbs1989
hydrogen-cyanide-aqueous-nbs1989
hydrogen-ion-aqueous-nbs1989
cyanide-ion-aqueous-nbs1989
hydroxide-ion-aqueous-nbs1989
water-liquid-nbs1989-thermochemistry
copper-i-ion-aqueous-nbs1989
dicyanocuprate-i-ion-aqueous-nbs1989
tricyanocuprate-i-ion-aqueous-nbs1989
tetracyanocuprate-i-ion-aqueous-nbs1989
```

For each profile, use the existing compiler:

```sh
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/PACK_ID/manifest.tsv \
  --out PACK_DIR/PACK_ID.fsmatpk
```

The eight following arguments are total CN, dissolved Cu(I), feed Na+, feed
Cl-, NaHS addition, HCl addition, NaOH addition and `n_ref`, all in mol:

```sh
cargo run -p fs-thermochem --example copper_cyanide_recovery -- \
  PACK_DIR 0.010 0.002 0.009 0.002 0.0011 0.0082 0.008 10
```

These are explicitly authored demonstration amounts, not pilot measurements.
The selected model gives:

| Stage | pH | Dissolved Cu (mol) | Free CN- (mol) | Copper-bound CN (mol) |
|---|---:|---:|---:|---:|
| Feed | 9.54451481733 | 0.002 | 0.00277690696214 | 0.00586887977378 |
| After NaHS and HCl | 4.98721781719 | 0.00000141916125625 | 0.000000568105709596 | 0.00000284171488609 |
| After solid removal and NaOH | 9.67043951985 | 0.00000141916125625 | 0.00732329957238 | 0.00000426124000235 |

The explicitly removed solid is 0.000999290419372 mol Cu2S. The dissolved
sulfide remainder is 0.000100709580628 mol. Residual copper binds somewhat
more cyanide after re-alkalization; the final calculation retains that demand
instead of assuming all formerly bound cyanide becomes free CN-.

The first stage computes the feed with no sulfide addition. The second adds
the declared NaHS and HCl ion equivalents. The third explicitly removes the
calculated primary solid and adds the declared NaOH equivalents to the
remaining liquor. The additions use a fixed solvent basis; they do not model
commercial reagent assay, dissolution heat, mixing volume or gypsum.

Primary solid removal is an explicit ideal separation action, not a predicted
filter efficiency. Residual dissolved Cu and S both enter the final phase
calculation, and any secondary solid remains in its reported inventory.
Whole-sequence balances include the removed primary solid. All stages must
succeed before the executable prints any successful result, including the
final pH-domain check on re-alkalization.

Focused verification uses the original source compiler and the same runnable
consumer:

```sh
cargo test -p xtask --test mining_cyanide_recovery_sources
cargo test -p fs-thermochem \
  stoichiometry::thermodynamics::acid_base::complexation::precipitation::tests
```
