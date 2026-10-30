# Violin varnish: source materials, section mechanics, and sound

The material database distinguishes historical coating observations from measured modern varnish specimens. The modern films have executable, frequency-dependent storage and loss moduli. They can be resolved with their original source receipts, applied as layers, and traced into added mass, membrane/bending stiffness, coupling, and material damping.

The runnable [violin_varnish example](../crates/fs-material/examples/violin_varnish.rs) performs that calculation with the existing violin spruce source. It answers a specific question: **how does an explicitly chosen outer film change the local mechanical section that a violin shell solver consumes?** Its directional frequency and decay estimates hold a spatial curvature pattern fixed. They are not a whole-instrument modal solution or a prediction of perceived quality.

## Historical evidence and executable surrogates

The [Stradivari stratigraphy pack](../data/matdb/seed-v1/varnish-stradivari-echard2010-stratigraphy/manifest.tsv) records observations from Echard and colleagues' five-instrument study: drying oil in the lower region and oil with Pinaceae resin in the upper coating, with instrument-dependent pigments. The reported upper coating was 10–30 micrometres; penetration into the sampled spruce and maple differed. These are observations of particular sampled locations, not a complete historical recipe. Chemical identification does not supply a density, modulus, loss spectrum, or proof of acoustic superiority. The historical pack therefore does not manufacture a mechanical law. [Echard et al., DOI 10.1002/anie.200905131](https://doi.org/10.1002/anie.200905131).

The six executable film packs retain two fitted specimens for each modern formulation in Lämmlein et al.: oil/copal varnish, an alcohol/shellac formulation, and an oil varnish containing pumice used as ground. Their fitted reference temperature is 20 °C; film conditioning includes at least three months at 20 °C/35% RH. Frequency dependence comes from dynamic mechanical measurements and time-temperature superposition, not direct measurements at every audible frequency. The committed use band is 20 Hz–20 kHz. These are modern comparison materials, not samples of Stradivari's varnish. [Lämmlein et al., DOI 10.1016/j.rinma.2020.100137](https://doi.org/10.1016/j.rinma.2020.100137).

| Modern family | Source specimen packs | Group density, kg/m³ |
| --- | --- | --- |
| Oil/copal | `varnish-oil-copal-lammlein2021-s1`, `-s2` | 1097; SD 26; n = 5 |
| Alcohol/shellac | `varnish-spirit-shellac-lammlein2021-s1`, `-s2` | 984; SD 30; n = 7 |
| Pumice-containing ground | `varnish-mineral-ground-lammlein2021-s1`, `-s2` | 1155; SD 30; n = 4 |

Density measurements are group statistics and are not matched to fitted specimen 1 or 2. Standard deviation is not a minimum/maximum range or a propagated uncertainty on a simulated violin. Film Poisson ratio and coating geometry remain authored inputs. See the individual packs for composition, conditioning, coefficients, and source restrictions. [Lämmlein et al., Tables 1–2](https://doi.org/10.1016/j.rinma.2020.100137).

## Run the quantitative section study

From the repository root, compile the existing spruce pack and all six films:

```bash
packs=target/violin-varnish-packs
mkdir -p "$packs"
for slug in \
  violin-spruce-norway-viala2020-shell \
  varnish-oil-copal-lammlein2021-s1 \
  varnish-oil-copal-lammlein2021-s2 \
  varnish-spirit-shellac-lammlein2021-s1 \
  varnish-spirit-shellac-lammlein2021-s2 \
  varnish-mineral-ground-lammlein2021-s1 \
  varnish-mineral-ground-lammlein2021-s2
do
  cargo run -p xtask -- matdb-pack \
    --manifest "data/matdb/seed-v1/$slug/manifest.tsv" \
    --out "$packs/$slug.fsmatpk"
done

cargo run -p fs-material --example violin_varnish -- "$packs" \
  > target/violin-varnish.tsv
```

Optional positional arguments are wood thickness in metres, external-film thickness in metres, and film Poisson ratio. Defaults are **0.003 m, 0.000020 m, and 0.30**. Compare 40 micrometres and the zero-film control using the same source cards:

```bash
cargo run -p fs-material --example violin_varnish -- "$packs" 0.003 0.000040 0.30
cargo run -p fs-material --example violin_varnish -- "$packs" 0.003 0 0.30
```

Zero thickness omits the film layer and recovers the bare section. Invalid inputs, unsupported source conditions, or a failed fit refuse. The successful report is emitted only after all source and section calculations succeed.

### Inputs and scope

The substrate uses the committed [violin spruce shell pack](../data/matdb/seed-v1/violin-spruce-norway-viala2020-shell/manifest.tsv), not a new table of invented wood constants. Its axes are longitudinal grain L, radial R in the plate, and tangential T through thickness. It retains the Viala tonewood-quarter means: density 460 kg/m³, E_L = 13.8 GPa, E_R = 1.2 GPa, G_LR = 0.8 GPa, G_LT = 0.85 GPa, and G_RT = 47 MPa, with longitudinal/radial losses near 0.0075/0.0175. Those specimens were modern Norway spruce quarters, not a particular historical violin. [Viala et al., DOI 10.1016/j.culher.2019.09.004](https://doi.org/10.1016/j.culher.2019.09.004).

The shell pack already declares nu_LR = 0.372 as a **Sitka spruce surrogate** from the USDA Wood Handbook because the Viala study did not identify it. That choice remains visible; no third Young modulus or additional Poisson ratios are invented for the plane-stress wood layer. [Wood Handbook, FPL-GTR-282, Chapter 5](https://research.fs.usda.gov/treesearch/62200).

Additional choices are explicit:

- A flat, aligned, perfectly bonded section, with one homogeneous film on its outer face. The reference surface stays at the original wood midsurface. Applied thickness is independent of the DMA specimen thickness; the study assumes the fitted bulk film response transfers to the authored coat.
- A real, frequency-independent isotropic film Poisson ratio of 0.30 and an authored small-strain envelope of 0.0003. The uniaxial film data do not identify that three-dimensional symmetry or Poisson ratio.
- Fixed source states. Wood and film measurements were conditioned separately; assembling their cards is an engineering comparison, not validation of an equilibrated coated specimen. The wood pack's temperature support includes an author-declared range. No moisture, thermal, curing, or aging correction is supplied.
- The source wood pack's existing **constant-loss** interpretation: `Q''wood = omega_reference * Vwood`. Multiplying its tabulated viscosity by each evaluation frequency would instead impose a different Kelvin–Voigt law.
- Ground specimens use the same authored external-film geometry as the other materials. This comparison does not establish where or how thick an actual ground should be. Penetration into wood needs a separately characterized affected wood zone, not another external film representing the same depth.

The wood source is limited to about 100–5000 Hz. At **20 Hz and 20 kHz**, the report supplies the valid film response and marks section estimates unavailable. At nominal 5000 Hz, the existing wood pack's printed upper bound is 31415.9265 rad/s, or 4999.999994286667 Hz. The report labels and evaluates that actual source edge and prints both requested and evaluated frequencies. Its printed lower bound, 628.318531 rad/s, is slightly above 2 pi 100; the body-response commands use 101 Hz to stay inside that retained bound rather than silently extrapolate or round it away.

## Physics represented by the example

### 1. The measured-fit film law

`ModifiedFractionalZener` implements the source's two-exponent law:

$$
E^*(\omega)=E_0+(E_{\mathrm{fit}}-E_0)
\frac{(i\omega\tau)^\alpha}{1+(i\omega\tau)^\beta}.
$$

Storage and loss are its real and imaginary parts; `eta = E''/E'`. The admitted region is `0 < beta <= alpha <= 1`. With unequal exponents the source symbol `Einf` is a **fit scale**, not a finite infinite-frequency modulus. Coefficients and frequency support travel through `ResolvedMaterialStatePoint` into `ResolvedCoating`. [Lämmlein et al., equations 1–3](https://doi.org/10.1016/j.rinma.2020.100137).

### 2. Mass and layer position

The laminate owner integrates each layer about the same reference:

$$
I_0=\int\rho\,dz,\qquad I_1=\int\rho z\,dz,\qquad I_2=\int\rho z^2\,dz.
$$

`I0` is areal mass, `I1` retains translational/rotational coupling, and `I2` retains rotary inertia. Material at the surface has a larger bending lever arm than material near the centre. Default oil-film added mass is `1097 * 20e-6 = 0.02194 kg/m²`, against bare wood's `460 * 0.003 = 1.38 kg/m²`: approximately a 1.590% increase.

### 3. Membrane, coupling, and bending stiffness

For oriented reduced plane-stress stiffness `Q*`, the section computes

$$
A^*=\int Q^*dz,\qquad B^*=\int Q^*z\,dz,\qquad D^*=\int Q^*z^2\,dz.
$$

A one-sided film generally makes `B` nonzero. The report includes all three complex 3×3 matrices in engineering coordinates `[xx, yy, xy]`. For **zero membrane resultant**, equilibrium gives

$$
\epsilon_0=-(A^*)^{-1}B^*\kappa,\qquad
D_{\mathrm{free}}^*=D^*-B^*(A^*)^{-1}B^*.
$$

The inverse action uses the existing complex-LU solver. This condensation applies to the stated local free-membrane comparison. A constrained or curved violin shell must retain its coupled membrane and bending degrees of freedom; it cannot universally replace its law by `Dfree`.

### 4. Directional frequency estimates

For unchanged spatial mode shape, fixed sampled material spectrum, and translational inertia, the local Rayleigh estimate is

$$
\frac{f_c}{f_b}\approx
\sqrt{\frac{\kappa^T\operatorname{Re}D_{\mathrm{free},c}^*\kappa}
{\kappa^T\operatorname{Re}D_{\mathrm{free},b}^*\kappa}
\frac{I_{0,b}}{I_{0,c}}},\qquad
\frac{\Delta f}{f}\approx\frac12\left(\frac{\Delta K}{K}-\frac{\Delta M}{M}\right).
$$

Pure longitudinal and radial curvatures isolate mechanisms; they are not named violin modes. An actual mode mixes curvature, membrane strain and shear across a curved, varying-thickness body, and can change shape after coating. A frequency-dependent eigenproblem must also evaluate the constitutive response at the resulting modal frequency.

### 5. Damping and decay from strain energy

The film's loss factor is not assigned to the whole plate. Each layer uses the relaxed harmonic strain `e(z) = epsilon0 + z*kappa`:

$$
S_i=\int e^\dagger Q_i'e\,dz,\qquad
L_i=\int e^\dagger Q_i''e\,dz,\qquad
\eta_{\mathrm{section}}=\frac{\sum_i L_i}{\sum_i S_i}.
$$

Summed layer weights are checked against the condensed bending response. The report separates `Lwood / sum(S)` and `Lfilm / sum(S)`. The film contribution differs from the net increase in damping because its presence also changes the wood's energy share. At each stated evaluation frequency, light-damping material-only estimates are

$$
\tau_{\mathrm{amplitude}}\approx\frac{1}{\pi f\eta},\qquad
T_{60}\approx\ln(1000)\tau_{\mathrm{amplitude}}.
$$

These include the retained wood and film losses. They exclude radiation, cavity losses, joints, contacts, bow/string interaction, and room reverberation. `T60` describes a conditional modal amplitude decay, not an audible recording or the decay of the entire violin.

## Reproduced default results at 1 kHz

These are **calculated section results** from the executable and compiled source packs, with the assumptions above. Frequency shifts use the square-root ratio. Total losses are dimensionless; decay includes the wood loss.

| Film specimen | L frequency shift | R frequency shift | Total eta, L | Total eta, R | Material T60, L | Material T60, R |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Bare spruce | 0% | 0% | 0.007500 | 0.017500 | 0.29317 s | 0.12565 s |
| Oil/copal 1 | −0.66284% | +0.60597% | 0.007710 | 0.019545 | 0.28519 s | 0.11250 s |
| Oil/copal 2 | −0.66174% | +0.61826% | 0.007733 | 0.019795 | 0.28434 s | 0.11108 s |
| Alcohol/shellac 1 | −0.55639% | +0.98009% | 0.007702 | 0.019382 | 0.28549 s | 0.11344 s |
| Alcohol/shellac 2 | −0.55290% | +1.01885% | 0.007720 | 0.019572 | 0.28481 s | 0.11235 s |
| Pumice ground 1 | −0.65608% | +1.09612% | 0.007684 | 0.019127 | 0.28616 s | 0.11496 s |
| Pumice ground 2 | −0.64832% | +1.18186% | 0.007739 | 0.019704 | 0.28411 s | 0.11159 s |

Spruce is much stiffer along the grain than across it. In the longitudinal comparison, film mass outweighs its small stiffness contribution; radially, its relative stiffness contribution is larger. For oil/copal specimen 1, the film stores about 0.2466% of the longitudinal strain-energy weight and 2.7215% of the radial weight, making its damping contribution much stronger radially.

Measured varnished-wood studies also report direction-dependent changes, including radial frequency increases and longitudinal decreases. Agreement of signs is useful context, not calibration of this section comparison or validation of a whole violin. [Lämmlein et al., Scientific Reports 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6901594/).

## Frequency slices, material memory, and the violin solver

`ResolvedCoating::at_frequency(h, omega)` returns a layer with source storage stiffness and viscosity `E''/omega`. At that frequency its complex response is the fitted law. The section example regenerates the slice at **each** frequency. Reusing one slice across a broad spectrum would keep storage constant and make loss proportional to frequency; that additional approximation must stay visible.

`ResolvedCoating::lower_to_tensor_prony` instead produces passive material memory using the existing tensor-Prony owner. The section example also performs this lowering for each source and prints its certificate. Sixteen candidate relaxation times yield thirteen positive scalar arms for these fits. Over the 512-point verification grid, worst relative storage error is about `8.314e-5` and worst relative loss error about `4.533e-4`, below independently authored `0.001` gates. These measure approximation to the fitted equation, not experimental uncertainty; a grid is not a continuous analytic bound. The section table uses exact frequency slices, not the Prony approximation.

### Two explicit coating choices

The default `ResolvedShearShell::with_coatings` route freezes storage at an authored reference and checks each retained mode against an explicit relative storage-drift limit. Each film's loss reads `E''(omega)/omega` at that mode's frequency; positive diagonal congruence preserves a passive modal damping matrix with geometric-mean cross terms. Wood retains its own source loss law. This is a modal approximation, not hereditary film memory or a self-consistent dispersive eigenproblem.

The additional `with_coating_relaxation` selection uses the source's **equilibrium modulus E0** in the original laminate/MITC eigensolve. Film mass, centroid and rotary inertia remain in I0/I1/I2; signed face offsets and membrane–bending coupling remain in A/B/D. The same equilibrium modal/bubble transformation maps every film's complete unit-modulus stiffness form H, including off-diagonal terms. Each positive Prony increment supplies a restoring branch proportional to H. Wood/glue damping and assembly-interface losses retain their own operators; the frozen film damping is not added a second time.

Equal relaxation times are combined across the entire body before the existing positive factorization. For the current O1 fit, 13 positive scalar times and 30 body modes require at most 390 histories, even though the varnish covers 60 regions in 45 components. A separate history bank per face would make the same physical transfer needlessly large. The aggregation preserves relaxed or uniformly unrelaxed initialization; it does not recover arbitrary independent initial histories of individual films.

The existing coupled dynamics owner advances these internal strains together with strings, body modes, cavity, bow/contact and external loads. For one factored branch, with r = cᵀq − z,

```text
v' = -K0 q - C v - Σ c r + applied_force
z' = r / tau
H  = (vᵀv + qᵀK0q + Σ r²) / 2
H' = applied_forceᵀv - vᵀCv - Σ r²/tau
```

Thus the same material branch changes the body restoring force and deposits nonnegative heat. The prepared backend now includes those histories in its existing exponential and work/heat integrals, with no new time stepper. Its original mode, byte, multiply-work and augmented-dimension limits remain binding.

The [body profile grammar](../crates/fs-couple/examples/material_strings/BODY_PROFILES.md#hereditary-varnish-and-a-structural-frequency-response) exposes the explicit policy:

```text
coating-memory prony-continuation 16 0.001 0.001 512
```

`prony-continuation` means that transient calculation accepts the globally passive fit outside its finite measured frequency band. The original source band and separate sampled storage/loss certificates remain visible; this is not a claim of measured DC or high-frequency behavior. Every retained body mode must still have source support. The fixed Galerkin/bubble basis also remains an approximation: it is not frequency-dependent local condensation.

The [hereditary il Cannone configuration](../crates/fs-couple/examples/material_strings/cannone-h8-varnish-memory.assembly) applies the same 20 µm modern O1 film to the exterior top, back and ribs as the [frozen-reference configuration](../crates/fs-couple/examples/material_strings/cannone-h8-varnished.assembly). The [violin recipe](../data/instruments/descriptions/violin/README.md#apply-a-measured-varnish-to-the-existing-violin) supplies build and render commands. The current recipe's conservative prepared-memory admission is 3,954,407,664 bytes, within its existing 4 GiB declaration. Its independent 512-history cap is not a promise that every case at that cap fits the byte budget. This resource calculation does not establish that a full violin render has completed.

Library users explicitly call `ReciprocalStringBank::with_material_shell_coating_relaxation` before playback; the frontend performs that binding when the policy record is present. Shell constructors alone return the equilibrium image and retained projection, without silently selecting a transient continuation.

### Trace bridge force to structural motion and loss

Both single-shell and assembled-shell images now expose `coating_harmonic_response`. It evaluates each film's **original fitted law at the common drive frequency**, rather than at separate natural frequencies. With the exp(-iωt) phasor convention and mass-normalized coordinates, the dynamic matrix is

```text
Z(omega) = K0 - omega² I - i omega (Cwood/glue + Cinterface)
           + Σfilm [(E'(omega) - E0) - i E''(omega)] Hfilm.
```

The existing complex LU solves Zq = F and checks the residual against the original equations. Returned complex bridge mobilities and signed exterior/cavity volume velocities show how changed stiffness and damping affect mechanical motion. Separate wood/glue, interface and per-film dissipation, supplied power and power residual make the loss accounting inspectable. Source bands are enforced at every frequency. Bare images are supported; frozen-film images refuse this exact-source calculation.

The renderer's `--body-response NEW.tsv MIN_HZ MAX_HZ POINTS BRIDGE_INDEX` exports these quantities on a logarithmic grid alongside an ordinary WAV render. For structural analysis alone, `--body-response-only` uses the same arguments but omits the WAV positional argument. It resolves and projects the original body geometry and materials without constructing the string bank, cavity dynamics, bow performer or prepared time propagator. The performance file still supplies body dimensions, material state, mode count and input limits. Only required body/interface source packs are loaded.

Both selections validate all samples before creating a new report, refuse overwrite, and limit total solve work and output size. The transfer is normalized to a 1 N peak bridge force; that linear normalization does not certify the material's finite-amplitude strain envelope. The [body-profile guide](../crates/fs-couple/examples/material_strings/BODY_PROFILES.md#hereditary-varnish-and-a-structural-frequency-response) shows a staged-run command.

This response is the **unloaded retained structural body**, including its fixed reduced wood/glue viscosity. It does not include string impedance, cavity pressure feedback or a full exterior wave solution. The ordinary WAV path separately couples the modified body to the strings, bow and cavity, then uses the existing acoustic observer. A changed bridge response is a mechanism for changed sound, not a perceptual quality score.

### Trace the loaded violin to pressure and coating heat

`ReciprocalStringBank::material_shell_bridge_harmonic_response` and the equivalent method on the actual `MultiBowStringBank` apply a physical bridge force through the retained shell force map into the **complete connected linear performer**. They call the original `CoupledModalSystem::harmonic_response`; they do not introduce another dynamics solver. All retained strings, including sympathetic strings, bridge connections, selected cavity coordinates, material histories and admitted exterior radiation poles remain in that solve. Bow/contact traction is explicitly released. This is a linear released transfer, not a linearization of nonlinear bowed playing.

The returned native pressure is the original coupled time owner's pressure observation evaluated on those same displacement and velocity phasors. Acceleration weights contribute `-omega² * weight * q`; existing real modal pressure states retain `Re(T) * v + Im(T) * omega_mode * q`. Replacing the latter by `T * v` at an arbitrary forcing frequency would change the configured observer. No amplitude normalization or acoustic gain is added. This native pressure retains the configured observer's physical scope; a compact volume observer does not become a converged exterior-wave receiver merely by querying it in frequency space.

An explicitly configured `shell_observer` or `one_way_observer` now also exposes the selected microphone through `ShellHarmonicReceiver`. It borrows the actual connected bank and performs the mechanics query internally, so unrelated same-sized modal results cannot be substituted. The closed-shell route uses the original reciprocal P1/P0 observation owner tied to the same complete BEM load, gas, basis and body/mouth traces. The acoustic-skin route uses the original ruled Burton–Miller producer, every retained body mode and the actual aperture-flow shares. Both reuse playback's source-coordinate maps, evaluate the existing producer at the forcing frequency, and retain complex pressure and its three spatial gradients. Raw propagation phase is preserved without a playback FIR's delay or taper, and observation adds no mechanical reaction or extra loss debit. The report names the operator/loading/skin identities and distinguishes loaded from one-way body modes; it does not certify omitted feedback as negligible.

Unlike the unloaded exact-source film calculation, this loaded response uses the **actual positive Prony realization bound for playback**. For each original Maxwell branch, with peak mismatch phasor `r = cᵀq - z`, the existing harmonic owner now also returns ordered branch diagnostics:

```text
mean branch storage = |r|² / 4
mean branch heat    = |r|² / (2 tau)
```

The bank sums the retained coating branch range to report coating heat and stored energy; other material heat is summed separately. Equal-time factorization can mix films, so the loaded report does not invent a per-film allocation. Component viscosity, connection damping, material branches and radiation retain distinct loss channels. Supplied power and the original solve/power residuals expose accounting errors rather than hiding them in a sound-pressure multiplier.

`--loaded-body-response NEW.tsv MIN_HZ MAX_HZ POINTS BRIDGE_INDEX`, with the normal performance and new-WAV positional arguments, exports all bridge mobilities, signed shell volume velocities, native and selected receiver pressure in Pa/N, the loss channels and stored energies before continuing the WAV render. The v2 TSV adds `drive_bridge`; it retains native pressure separately and includes receiver gradients when an exterior receiver is selected. Without an exterior declaration, the selected receiver pressure equals native pressure. It interrogates the assembled performer without advancing accepted states or histories. This option still constructs the complete playback bank and prepared propagator. The original `prepared-linear` multiply-work declaration bounds the **entire sweep**, divided across its frequencies; small nonprepared runs retain their existing work allowance. No numerical cap is silently raised.

Set the final bridge argument to `all` to export the full bridge transfer matrix, with one 1 N peak input per column and rows ordered by frequency then declared bridge. `CoupledModalSystem::harmonic_responses` assembles and factors the original dynamic stiffness once per frequency and uses the existing cancellable complex LU solve for every right-hand side. Its batch work allowance is `d² [d + G + r(C + M + 3)]`, for `d` mechanical-plus-air coordinates, `G` constant-loss links, `C` connections, `M` material branches and `r` force vectors; at `r = 1` it equals the scalar allowance. The existing unknown cap remains unchanged; batches additionally admit at most 64 right-hand sides and 8 MiB of retained response storage. The selected microphone shares one BEM frequency sample across the matrix. Its boundary, quadrature and projection work is reserved from the same per-frequency allowance before mechanics.

The banks' `material_shell_bridge_forces_harmonic_response` and prepared receiver's `bridge_forces_harmonic_response` accept one complex peak force per physical bridge. They combine generalized forces **before** solving and calculating energy, so coherent cross-power and cross-heat remain present. Linear pressure or displacement can be combined from unit-force matrix columns; isolated column powers and coating losses cannot simply be added for coherent excitation. Focused tests compare batch and scalar responses bit-for-bit and verify phase rotation, exact cancellation and these cross terms. The three-input work fixture requires 198 work units with shared factorization versus 270 for three scalar solves; this is an admission/operation-count example, not an end-to-end violin timing.

For frequency analysis alone, `--loaded-body-response-only` omits the WAV positional argument and constructs the same connected bank with `ModalPreparation::HarmonicOnly`. This explicit original-owner policy retains all physical coefficients and the same constant-loss damping directions, while constructing **no coupled time propagator or held-step bow compliance**. The analysis still resolves source materials, reduces the shell and initializes the original physical components. Material, air, connection and internal-component rebuilds preserve the construction purpose. Time advancement and held-step compliance requests refuse before changing accepted states.

The frontend shares one physical bank-construction sequence between analysis and playback, then schedules gestures only for playback. Loaded-only analysis requires the existing explicit `prepared-linear` resource record; its mode, byte and work allowances now bound harmonic construction and the complete sweep, without requesting time exponentials. Admitted geometry-derived shell BEM and skin loads remain available in this mode; the existing prepared-playback restriction on non-skin computed loads applies only to actual playback preparation. Source, geometry, frequency-band and acoustic-load admission remain binding. The [violin recipe](../data/instruments/descriptions/violin/README.md#inspect-the-loaded-bridge-to-pressure-response) gives the no-WAV command.

The original harmonic owner also assembles each connection or material rank form on its exact nonzero support. For 390 film branches confined to 30 body coordinates inside a 445-mode bank, this reduces material outer-product updates from 77,229,750 to 351,000 per frequency. This is an operation-count comparison, not a measured end-to-end speedup. Every retained mode, acoustic pole and material history remains present; dense LU and the original-equation residual remain unchanged. No coefficient threshold, mode truncation or new material approximation is introduced.

Every sample must lie within the intersection of retained string, wood, coating, interface and selected exterior-fit support, and any bound cylinder-flow law must admit the frequency. An exterior receiver also retains its declared pass band and frequency-sample allowance. Frozen coating images, unbound coating memory, unsupported nonlinear/adhesive/contact or local-strain attachments, and residual-flexibility extensions without complete direct bridge-force columns refuse. A frontend plate microphone remains unsupported by this shell-specific report. Both report forms retain the 2–1024 point and 8 MiB limits and create a new file only after every sample passes.

For a controlled varnish comparison, stage bare and hereditary-coated versions with identical geometry, mode count, strings, gas, connections and observation. Compare pressure magnitude **and phase**, bridge mobility and coating heat at common supported frequencies. The loaded comparison includes the shift in string/body/cavity response caused by the coating; the separate unloaded report isolates the body's exact fitted film law. Neither report alone establishes an audible improvement or recording-equivalent violin sound.

### Verification and remaining qualification

Seven focused native shell checks cover the physical section and projection path, equilibrium stiffness counted once, unchanged coating mass, full off-diagonal reconstruction, unique-time aggregation, and transactional refusals. Four original prepared-owner tests cover state, work and heat against an independent full fs-phs realization, step subdivision, histories beyond the previous small-state ceiling, and stiff relaxation. Four unloaded harmonic checks cover force/velocity and power identities, source support and bare/assembly response.

Four further native checks exercise the production coupled harmonic assembly, complex LU, pressure observation and branch diagnostics against independently eliminated loaded equations, including an air-pole antiresonance, phase/energy formulas and state-preserving refusals. Seven structural/parser/frontend checks use the source pack and original report writer, including the analysis-only route; one loaded-writer check carries an actual production coupled solution through a small transport facade and validates the physical TSV. The focused native containers do not compile the entire bank dependency graph. Bank integration and source-to-projection equivalence regressions are included in the normal Cargo test seams but were not executed here.

The harmonic-only extension has two additional focused native core checks using the original coupled, prepared, coefficient and harmonic owners, complete `fs-phs` and original linear-algebra kernels. Analysis and prepared playback produce bit-for-bit identical responses at four frequencies with constant-loss connections, signed Maxwell branches and an air pole, including after coefficient updates. A 32 KiB budget admits that small analysis fixture while refusing time preparation; this is not a memory measurement of the complete violin. Time/prediction/compliance refusals preserve state, history and clock. Two exact-support tests compare the rank assembly against independent dense outer products; the two production CLI parser tests also cover loaded-only operation without a WAV argument. Full-bank, closed-shell BEM and friction integration regressions remain authored but unrun.

All four new production batch tests pass in the focused native construction harness, covering scalar parity, coherent cross-power and branch heat, reciprocity, phase cancellation, bounded work/output and unchanged accepted state on refusal. The complete production harmonic test module reports 13 passed and one pre-existing failure: `harmonic_two_body_response_matches_independent_dynamic_stiffness_and_power` refuses an unresolved negative physical-heat-resistance eigenvalue during its initial time step, before calling a harmonic method. The identical setup failure was reproduced using the pre-batch source at `c72b13d1881ace6718aefc3c4937ef2a5358b25f`. The four independent loaded-equation/pressure checks still pass. Three production parser checks pass, including `all` selection. Three focused receiver-helper checks pass for complex phase, cancellation, shared neck-flow projection, coordinate/refusal checks, combined budgets and source batches beyond 256. These helper checks do not execute a BEM boundary solve. Closed-shell microphone, acoustic-skin microphone and full-bank matrix regressions are authored in the normal frontend test seams and remain unrun.

Three current v2 writer checks also pass, retaining the real production coupled and released-pressure calculations behind small bank/receiver transport facades. They verify single/all bridge row ordering, actual native pressure/mobility/coating heat, distinct selected-receiver pressure, all six pressure-gradient columns and their units, one receiver preparation per sweep, create-new preservation, and refusal without partial output. The selected microphone transfer in this writer seam is a transport fixture, not a BEM validation. Together with the three exact production parser tests, these checks verify the new report path without claiming a complete CLI build.

These checks exercise the original numerical owners through focused native builds. The complete Cargo workspace and full violin WAV were not executed in that environment because required sibling constellation repositories were absent. No mesh/mode/time convergence, evolving coating-strain certificate, listening comparison or recording-level validation is claimed. Those are still required before attributing a particular audible difference to varnish.

## Remaining physical limits

Moisture protection is another pathway: altered transport can later change wood mass, stiffness, swelling, prestress and damping. Experiments show that complete coatings retard uptake and worn regions protect less effectively; uncoated interiors and f-holes still permit exchange. A coating is not a hermetic boundary. [Lämmlein et al., Scientific Reports 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6901594/).

These packs do not identify a complete moisture permeability/sorption law, cure kinetics, aging law, adhesion failure model, thermal transport law, or penetrated-wood constitutive tensor for a particular historical varnish. Those effects need additional measurements and source-bound models. Pigment or resin identification, historical thickness observations, and modern dynamic-modulus fits cannot fill those missing properties by themselves.
