# Solution recovery: solvent extraction and carbon adsorption

FrankenMining's comprehensive plan §17.3 calls for source-specific extraction
and stripping isotherms between copper leaching and electrowinning. The
`fs_material::solvent_extraction` owner now supplies a finite ideal contact
that can be composed into those explicit process stages.

## Original isotherm and complete metal balance

`SolventExtractionIsotherm::from_card` reads the aqueous/organic concentration
knots from the original `ConstitutiveModelCard`. Each pair is expressed in
mol of the same analytical element per m³ of its own phase. The ordinate is
always organic concentration, including for stripping. There are 2–256
ordered pairs; aqueous concentration increases strictly and organic loading
does not decrease. No origin, smoothing or extrapolated segment is invented.

The card retains its original law version, source hashes, provenance and typed
validity. Every actual contact supplies both phase volumes, both initial
metal inventories, temperature and all additional source conditions. The
owner supplies the solved concentrations and actual organic/aqueous volume
ratio; callers cannot override those coordinates.

For aqueous concentration `x` and source organic isotherm `y(x)`, the contact
solves `Vaq*x + Vorg*y(x) = naq_initial + norg_initial`. The piecewise-linear
intersection is algebraic. Results retain the original model, initial/final
phases, signed metal transfer, source conditions and a separately recomputed
conservation residual. A loaded organic phase may release metal under the
same axis convention. An unsupported equilibrium refuses without changing
the inputs. Initial nonequilibrium concentrations need not lie on the curve.

The concentration convention follows [IUPAC D01817](https://doi.org/10.1351/goldbook.D01817).
Explicit phase-volume balances are also described in
[Metals 12(6), 1053, equations 4–7](https://doi.org/10.3390/met12061053).
These references define the model form; they do not calibrate the authored
example for a commercial extractant or a particular ore.

## Runnable extraction → stripping → electrowinning

```bash
cargo run -p fs-material --example mineral_solvent_extraction
cargo test -p fs-material solvent_extraction::tests
```

The example has two explicitly synthetic source tables and a 100 mol Cu feed.
Its extraction stage loads 66.666667 mol Cu into organic. The following
stripping contact transfers 53.333333 mol into 0.1 m³ of strip solution.
The example explicitly assigns that analytical Cu to dissolved Cu(II) and
supplies water and proton inventories before calling the existing Faraday
and exact reaction-inventory owner. At 3,600,000 C and 95% current efficiency,
the modeled productive pathway deposits 1.126219474 kg Cu. Copper remaining
in raffinate, organic and electrolyte remains in the reported metal account.

Optional arguments are `FEED_CU_MOL PASSED_CHARGE_C`. An unsupported feed
equilibrium or charge demand exceeding strip-liquor inventory is refused.
Zero passed charge leaves the strip copper available for subsequent work.

This operator tracks one analytical element. Aqueous species/charge balance,
proton and extractant consumption, heat, contact kinetics, phase entrainment,
stage efficiency, volume changes and multistage countercurrent convergence
need their own admitted models. Binary64 balance residuals are numerical
diagnostics, not experimental uncertainty bounds or certificates.

## Finite solution/carbon adsorption and desorption

`fs_material::aqueous_adsorption::AqueousLangmuir` supplies the missing aqueous
counterpart to the separate gas/solid adsorption owner. Its original card
declares saturation loading in mol metal/kg dry adsorbent, half-saturation
concentration in mol metal/m³ solution, source hashes, temperature support
and every additional chemical condition. The isotherm is
`q = q_max*c/(c_half+c)`, following
[USGS PHREEQC example 19, equation 35](https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqc3-html/phreeqc3-81.htm).

Each contact solves that relation together with
`V*c + m*q = initial_dissolved + initial_adsorbed`. The result retains the
original model, both phase inventories, final concentration/loading, independently
evaluated vacant capacity, signed uptake or release, and the metal-balance
residual. A stable scaled quadratic avoids subtracting nearly equal positive
roots. Unrepresentable positive inventories refuse instead of becoming zeros.
Source concentration support is checked at the solved equilibrium. Initially
loaded carbon can desorb; its stock must not exceed the supplied dry-carbon
capacity. Zero carbon preserves the dissolved stock and leaves loading undefined.

```bash
cargo run -p fs-material --example mineral_carbon_adsorption
cargo run -p fs-material --example mineral_carbon_adsorption -- 0
cargo test -p fs-material aqueous_adsorption::tests
```

The explicitly synthetic example starts with 10 mol analytical gold in 100 m³
solution and 100 kg dry carbon. With `q_max=0.2 mol/kg` and
`c_half=0.05 mol/m³`, equilibrium leaves 2.807764064 mol in solution and
7.192235936 mol on carbon. Separating the original liquid and contacting that
loaded carbon with 100 m³ of clean liquid releases 1.831121225 mol into the
new bath under the same isotherm. These values demonstrate the finite balance;
they are not measured gold/carbon performance or an elution-process prediction.

The owner conserves one analytical element at fixed solution volume and dry
adsorbent mass. It does not infer cyanide consumption, charge/speciation,
competitive adsorption, preg-robbing, kinetic contact time or CIL/CIP recovery.
Measured chemical support and separate kinetic/transport models remain needed
for those process claims. Six focused tests check independent equilibrium
answers and a bisection oracle, source admission, scaling, loaded-carbon release,
zero inventories and retained small vacant-site amounts.
