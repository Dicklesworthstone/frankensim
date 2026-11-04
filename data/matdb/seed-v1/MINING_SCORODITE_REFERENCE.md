# Matched arsenate and scorodite Gibbs reference

Nine profiles supply a coherent 298.15 K, 100000 Pa reference for arsenate
acid/base speciation and a selected crystalline scorodite equilibrium. There
are nine formation-Gibbs claims and eight separate source-reported uncertainty
magnitudes. Every claim requires `source-elemental-reference-nmk2014=1`.

## Source values

The first six values are the selection of Nordstrom, Majzlan and Königsberger
(2014), abbreviated NMK2014. The remaining three are compatible CODATA
auxiliaries. All energies and printed plus/minus magnitudes below are kJ/mol.

| Profile | Species | Formation G | Source-reported +/- magnitude |
|---|---|---:|---:|
| `arsenic-acid-aqueous-nmk2014` | H3AsO4(aq) | -766.420 | 0.35 |
| `dihydrogen-arsenate-ion-aqueous-nmk2014` | H2AsO4-(aq) | -753.575 | 0.35 |
| `hydrogen-arsenate-ion-aqueous-nmk2014` | HAsO4(2-)(aq) | -713.727 | 0.35 |
| `arsenate-ion-aqueous-nmk2014` | AsO4(3-)(aq) | -647.618 | 1.5 |
| `iron-iii-ion-aqueous-nmk2014` | Fe3+(aq) | -16.28 | 1.1 |
| `scorodite-solid-nmk2014` | FeAsO4.2H2O(cr) | -1284.8 | 2.9 |
| `water-liquid-nmk2014` | H2O(l) | -237.140 | 0.041 |
| `hydrogen-ion-aqueous-nmk2014` | H+(aq) | 0.000 | Defined zero; no error assigned |
| `hydroxide-ion-aqueous-nmk2014` | OH-(aq) | -157.220 | 0.072 |

The arsenate rows are NMK2014 Table 8, p. 227; Fe3+ and scorodite are Table 9,
p. 229. NMK2014 explicitly makes its Fe-As reference consistent with Parker
and Khodakovsky (1995) for Fe3+, Nordstrom and Archer (2003) for elemental
arsenic, and Cox et al. (1989) CODATA for water and elemental hydrogen/oxygen.
This selection is distinct from the older `nbs1989` arsenic profiles.

The water, H+ and OH- Gibbs values are printed in OECD NEA *Chemical
Thermodynamics of Uranium*, Table IV.1, p. 64. That table identifies the
water/OH- Gibbs values as internally calculated, rounded to three decimals.
They agree with the CODATA formation enthalpies and elemental/individual-ion
entropy conventions. The auxiliary table specifies 298.15 K, 0.1 MPa and
aqueous infinite dilution. All aqueous profiles retain unit-activity standard
states referenced to infinite dilution; individual ions retain H+(aq)=0.

The property `source-reported-formation-gibbs-uncertainty` preserves each
printed plus/minus magnitude as a dimensional scalar. This scalar does not
assert its own confidence or a distribution. The main Gibbs uncertainty is
Unstated for the six NMK2014 rows because no common coverage or covariance is
established here. For water and OH-, the NEA table describes total
uncertainties corresponding **in principle** to 95% confidence; those two
main claims retain `HalfWidth` values of 41 and 72 J/mol at confidence 0.95.
H+ is a defining zero with no experimental error claim. No independent-error,
Gaussian or propagated reaction-uncertainty assumption is made.

## Source documents and redistribution

- Nordstrom, D. K., Majzlan, J. and Königsberger, E. (2014), *Thermodynamic
  Properties for Arsenic Minerals and Aqueous Species*, Reviews in Mineralogy
  and Geochemistry 79, 217-255: <https://doi.org/10.2138/rmg.2014.79.4>.
  [Author-deposited full text](https://www.researchgate.net/publication/267751439_Thermodynamic_Properties_for_Arsenic_Minerals_and_Aqueous_Species).
  The inspected PDF has SHA256
  `1a63faf42405186c1cbacf10a2bf9a989022eb8efbf8c4702b5ebedc5f41c2c0`.
- Parker, V. B. and Khodakovsky, I. L. (1995), *Thermodynamic Properties of
  the Aqueous Ions (2+ and 3+) of Iron and the Key Compounds of Iron*, JPCRD
  24, 1699-1745: <https://srd.nist.gov/jpcrdreprint/1.555964.pdf>.
  The inspected PDF has SHA256
  `ff7f8aee30aac01e8f42a89c579ddd1c6df5ec39fc82b91e2d8992e435335bef`.
- Cox, J. D., Wagman, D. D. and Medvedev, V. A. (1989), *CODATA Key Values
  for Thermodynamics*, official selected table:
  <https://codata.info/resources/databases/key1.html>.
  The inspected HTML has SHA256
  `88049de2c7700bf44d233dacee019ac70c9b46ad280596270ebb3d5d0235d460`.
- Grenthe et al., *Chemical Thermodynamics of Uranium*, OECD Nuclear Energy
  Agency, Table IV.1 p. 64:
  <https://cms.oecd-nea.org/upload/docs/application/pdf/2019-12/uranium.pdf>.

Sources were checked on 2026-11-04. These packs redistribute selected
attributed numerical facts and original notes. They grant no license to the
source articles, books, figures or table layouts, none of which are included.

## Independent numerical checks

Use the products-minus-reactants sum of the matched formation-Gibbs values.
At 298.15 K with R = 8.31446261815324 J/(mol K):

| Reaction | Reaction G, kJ/mol | Derived -log10(K) |
|---|---:|---:|
| H3AsO4 -> H+ + H2AsO4- | 12.845 | 2.2503466391 |
| H2AsO4- -> H+ + HAsO4(2-) | 39.848 | 6.9810675652 |
| HAsO4(2-) -> H+ + AsO4(3-) | 66.109 | 11.5817957154 |
| H2O(l) -> H+ + OH- | 79.920 | 14.0013782326 |
| FeAsO4.2H2O(cr) -> Fe3+ + AsO4(3-) + 2 H2O(l) | 146.622 | 25.6870630532 |
| FeAsO4.2H2O(cr) + 3 H+ -> Fe3+ + H3AsO4 + 2 H2O(l) | 27.820 | 4.8738531335 |

The extra digits document arithmetic, not experimental precision. NMK2014
Table 10 gives log10(Ksp) = -25.68 and an acid-form dissolution constant
log10(*Ks0) = -4.86. The literal rounded formation values give the nearby
numbers above; they are not adjusted to force exact agreement. The older
log10(Ksp) = -25.83 uses a different arsenate selection and must not be
substituted into this family.

## Physical scope

The solid is pure crystalline bulk scorodite with exactly two waters of
crystallization. Its dissolution releases one Fe(III) and one As(V) per mole.
The conventional Ksp uses **free** Fe3+ and AsO4(3-) activities. Analytical
dissolved arsenic must be partitioned among all four arsenate protonation
states before applying this product. Water activity must be declared; setting
it to one is the pure-solvent approximation.

An ideal calculation that selects only this solid and free ferric ion is a
restricted equilibrium calculation. NMK2014 describes aqueous Fe-As
complexes, Fe(III) hydrolysis and competing ferric oxyhydroxides. Competition
can cause incongruent scorodite dissolution above roughly pH 2-4, while
preparation, grain size and aging affect observed behavior. Those mechanisms
are not represented by these nine profiles. They do not establish total
mine-water solubility, a precipitation rate, permanent fixation or plant
removal performance. The packs supply only point Gibbs data; no caloric
properties, heat capacity or temperature extrapolation is inferred.
