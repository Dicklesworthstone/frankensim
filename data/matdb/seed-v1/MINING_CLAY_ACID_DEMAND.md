# Clay acid demand and staged aluminum precipitation

Kaolinite dissolution can consume the same acid inventory as a copper-oxide
leach. Subsequent precipitation of the dissolved aluminum can release acid,
but that later production cannot fund an earlier leaching step. The existing
`ReactionNetwork` now accepts both endpoints with a matched set of source
thermochemistry, and `mineral_processing` runs the finite sequence.

Reaction extents and the selected precipitate are caller inputs. These
calculations do not establish how much kaolinite dissolves at room temperature
or how quickly it reacts. A density/hardness record for a kaolinite specimen
does not supply its acid-accessible fraction.

## Five compatible source packs

The new profiles transcribe the 298.15 K columns of the
[NBS tables with 1989 errata](https://doi.org/10.18434/M32124), using the
[NIST corrected CSV](https://data.nist.gov/od/ds/mds2-2124/NBS_Tables_with_ErratumValues.csv).
The selected CSV has SHA-256
`cb8a8b8b94a287e6e8c82eed8bad2cf9b611d574d1cfd95c8528d318dfef4b44`.
The [original tables](https://srd.nist.gov/JPCRD/jpcrdS2Vol11.pdf) define the
standard states in sections 5.5–5.8; the
[NIST digitization paper](https://doi.org/10.6028/jres.125.007) describes the
columns and state notation. These are attributed numerical facts and original
notes, not redistributed source documents or table images.

| Pack | Source row | M, g/mol | Formation H, kJ/mol | Formation G, kJ/mol | Entropy, J/(mol K) | Cp, J/(mol K) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `kaolinite-solid-nbs1989-thermochemistry` | Table 29:Al, Al2Si2O7·2H2O, kaolinite, cr | 258.1616 | −4119.6 | −3799.7 | 205.0 | 246.14 |
| `aluminum-iii-ion-aqueous-nbs1989` | Table 29:Al, Al3+, ao | 26.9815 | −531. | −485. | −321.7, conventional ionic | Unavailable |
| `silicic-acid-component-aqueous-nbs1989` | Table 24:Si, H4SiO4, ao | 96.1156 | −1468.6 | −1316.6 | 180. | Unavailable |
| `gibbsite-solid-nbs1989-thermochemistry` | Table 29:Al, Al2O3·3H2O, gibbsite, cr | 156.0074 | −2586.67 | −2310.21 | 136.90 | 183.47 |
| `copper-ii-oxide-solid-nbs1989-thermochemistry` | Table 36:Cu, CuO, cr | 79.5394 | −157.3 | −129.7 | 42.63 | 42.30 |

All **23 claims** are points at **298.15 K and 100000 Pa**, with
`source-elemental-reference-nbs1989=1`, `source-standard-state=1`,
`engineering-reference-model=1` and `source-specimen-qualified=0`. The Cp
points do not license temperature extrapolation. Row-level uncertainty and
cross-species covariance are not supplied, so uncertainty remains `Unstated`.
The separate JANAF and USGS CuO evaluations remain distinct source families.

Three identities are essential to using these records correctly:

- **Kaolinite:** the source writes Al2Si2O7·2H2O, which has the same atom counts
  as Al2Si2O5(OH)4. The reaction ID `Al2Si2O5(OH)4(s-kaolinite)` retains the
  crystalline phase. It does not represent metakaolin or a generic clay mixture.
- **Silicic component:** NBS explicitly defines the aqueous H4SiO4 row as
  **H2SiO3(ao) + H2O(l)**. Its water formula unit remains in mass, element and
  energy accounting. The required selector is
  `source-silicic-component-h2sio3-plus-water=1`. The adjacent crystalline
  H4SiO4 row, unhydrated component, total dissolved silica, and a silica gel
  are not interchangeable with it.
- **Gibbsite:** NBS reports **Al2O3·3H2O**, equivalent to **two Al(OH)3 units**.
  Every molar quantity and reaction extent retains that doubled basis through
  `Al2O3.3H2O(s-gibbsite)` and `source-formula-aluminum-atoms=2`. It is not an
  amorphous hydroxide or another aluminum oxyhydroxide phase.

The aqueous Al3+ record uses the hypothetical ideal unit-molality standard
referred to infinite dilution. Its negative entropy is conventional under
the H+(aq)-zero reference. It represents free Al3+, not total dissolved aluminum
or hydrolyzed/sulfate-complexed aluminum.

## Conserved endpoints and source-driven thermochemistry

The kaolinite and gibbsite equations follow the
[USGS PHREEQC phase listing](https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/html/final-90.html);
the [USGS reaction-path example](https://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc/phreeqc3-html/phreeqc3-61.htm)
also retains kaolinite's six-proton, H4SiO4 endpoint. These references establish
the equation forms. Their equilibrium constants are **not** substituted for
the separately tabulated NBS formation data. The study
[*The effect of pH and temperature on kaolinite dissolution rate under acidic conditions*](https://doi.org/10.1016/S0016-7037(02)00966-3)
uses the same net dissolution equation while treating kinetics separately.

| `MineralReaction` | Written equation | Standard ΔH, J/mol extent | Direct standard ΔG, J/mol extent |
| --- | --- | ---: | ---: |
| `KaoliniteAcidLeach` | Al2Si2O5(OH)4 + 6 H+ → 2 Al3+ + 2 H4SiO4 + H2O | −165430 | −40629 |
| `AluminumHydrolysisGibbsite` | 2 Al3+ + 6 H2O → Al2O3·3H2O + 6 H+ | +190310 | +82564 |
| Existing `CupricOxideAcidLeach`, new NBS source route | CuO + 2 H+ → Cu2+ + H2O | −63760 | −41939 |

These sums use the three existing compatible packs
`hydrogen-ion-aqueous-nbs1989`, `water-liquid-nbs1989-thermochemistry` and
`copper-ii-ion-aqueous-nbs1989`. All species, including standard-state zeros,
need a source. Gibbs energy comes directly from the formation-G table; it is
not inferred by subtracting a species entropy from formation enthalpy.

Combining one clay extent with one hydrolysis extent cancels six protons and
two aluminum ions, leaving kaolinite + 5 H2O → one doubled gibbsite unit +
2 H4SiO4. That aggregate identity does not make acid available before the
hydrolysis stage, establish a stable solid, or predict that either extent occurs.

## Run the finite shared-acid calculation

Compile the eight participants into a directory using the existing compiler:

```bash
mkdir -p /tmp/clay-acid-packs
for pack in \
  kaolinite-solid-nbs1989-thermochemistry \
  aluminum-iii-ion-aqueous-nbs1989 \
  silicic-acid-component-aqueous-nbs1989 \
  gibbsite-solid-nbs1989-thermochemistry \
  copper-ii-oxide-solid-nbs1989-thermochemistry \
  copper-ii-ion-aqueous-nbs1989 \
  hydrogen-ion-aqueous-nbs1989 \
  water-liquid-nbs1989-thermochemistry
do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$pack/manifest.tsv" \
    --out "/tmp/clay-acid-packs/$pack.fsmatpk"
done

cargo run -p fs-thermochem --example mineral_processing -- \
  --clay-acid-inventory /tmp/clay-acid-packs \
  2 1 5 100 1 0.5 0.5 1
```

The first four numbers are initial CuO, kaolinite, proton and water amounts
in moles. Products start at zero. The next four numbers prescribe the first
CuO extent, clay extent, gibbsite extent and later CuO extent, each in moles of
its complete written equation. The constructed example has this result:

| State | CuO, mol | Kaolinite, mol | H+, mol | Al3+, mol | Cu2+, mol | Gibbsite source units, mol | H4SiO4, mol | Water, mol |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Initial | 2 | 1 | 5 | 0 | 0 | 0 | 0 | 100 |
| After shared-acid leach | 1 | 0.5 | 0 | 1 | 1 | 0 | 1 | 101.5 |
| After selected hydrolysis | 1 | 0.5 | 3 | 0 | 1 | 0.5 | 1 | 98.5 |
| After later CuO leach | 0 | 0.5 | 1 | 0 | 2 | 0.5 | 1 | 99.5 |

The standard chemical enthalpy changes are **−146475 J**, **+95155 J** and
**−63760 J**, totaling **−115080 J**. Initial and final source-molar-mass
inventories both sum to **2.2238204 kg** within floating-point arithmetic.
No acid or metal disappears into a balancing stream.

Requesting two first-stage CuO extents with the same 0.5 clay extent demands
seven moles of H+ from five and refuses. After the accepted first stage,
requesting hydrolysis and later CuO leaching simultaneously also refuses:
the initial proton stock for that stage is zero. The caller must first accept
the hydrolysis stage, then pass its returned inventory to the later leach.
Aluminum and water availability independently bound the hydrolysis extent.

## Individual source modes

`mineral_processing --help` lists the exact pack order for:

- `--reaction kaolinite-acid-leach`: kaolinite, H+, Al3+, H4SiO4 component, water.
- `--reaction aluminum-hydrolysis-gibbsite`: Al3+, water, doubled-unit gibbsite, H+.
- `--reaction cupric-oxide-acid-leach-nbs1989`: CuO, H+, Cu2+, water.

Each takes compiled pack paths followed by kilograms of the **first species
actually reacted**. For example, 0.0269815 kg of Al3+ corresponds to one mole
of aluminum and **0.5 mol** of the written gibbsite-hydrolysis equation.
`--gibbs` resolves the direct source formation-G route. Existing `--activity`
diagnostics require all relative activities; they do not solve concentrations.

The finite example is net ionic accounting. Unchanged acid counterions lie
outside these active reactions. Proton equivalents are not pH, commercial acid
mass, an electroneutral solution specification, or a plant dosing instruction.
Room-temperature dissolution, reaction time, nonideal activity coefficients,
aluminum complexes, precipitation completion and heat of mixing remain separate
models. No operating temperature or heat duty is inferred from the reference
enthalpy balance.

The focused compiler/consumer regression is
`cargo test -p xtask --test mining_clay_acid_sources`. It checks all 23 facts,
source domains, phase/component/formula-unit refusals, three independent H/G
sums, finite acid competition, staged regeneration, and retained mass/Al/Si.
