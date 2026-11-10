# CONTRACT: fs-couple

> **Layout (restored 2026-11-09).** This contract is an index plus part files. Once this file grew
> past 1 MiB (1,048,576 bytes), some commits replaced the whole contract with their one new
> section; no edit below that size ever did. The cause has not been identified: the cutting commits
> are ordinary multi-file commits, so a tool that cannot read past 1 MiB and then writes the file
> whole is the likely kind of cause.
> This file was cut to a single section 18 times between 2026-10-30 and 2026-11-09; the text written before 2026-11-05 19:45 (2.3 MB) was absent from then until this restore.
>
> - All earlier contract text is in `contracts/part-NN.md`, unchanged. A higher part number is newer;
>   inside a part, sections are newest first, as they were in the single file.
> - Add a new section under "New sections" below, newest first.
> - There is no size limit on this file or on the parts (owner ruling, 2026-11-09). If a tool shows
>   you this file empty or cut short, do not write it back whole: add your section with a patch, and
>   check that `git diff --stat` shows no large deletion here before you commit.
> - The ten required sections are in the parts named under their headings below.

## New sections (newest first)

### Reed opening tap

`ReedBoreVoice::record_reed_opening` and `ReedPerformance::record_reed_opening` give a massive reed an opt-in tap: at the end of every solver sample one text line `y v` (the blade's opening in metres and velocity in m/s, as the closure stop and the lay see them) goes to a caller-supplied writer. It reads state after the step and changes no audio sample (test `reed_trace_tests`: the pressure bits with and without the tap are equal, and with the stop on `y` never goes below zero, with `closure_stop off` it can). The `music_render wind` command routes the existing `--aperture-flows FILE` to it on a `frankensim-reed-performance-v1` file (`--labium` stays flue-only); on any other file type the flag keeps its earlier meaning and refusals.

### The reed-noise row against the print's power (measured)

A held, passive reed (damping ratio 30, it cannot oscillate) on an 8 m cylinder (radius 6.935 mm; the characteristic-line ladder refuses this card at 9, 10, 15 and 25 m) with the blow held at 1800 Pa was run with `constriction_noise` on and off, and the outgoing wave p+ read in the window of samples 800 to 2200 (the first reflection returns at 2239), where the bore carries only the wave the slit launches. The print gives, per direction, W = K(St)^2 A dp^2 / (4 rho c) (Oldham and Ukpoho as printed in Guedel, AIVC CR02 eq. 3), an outgoing wave of rms K dp / 2 with dp the static drop across the slit (1368 Pa here). Measured per third octave from 630 to 8000 Hz (12 bands) the row's wave is below that by 6.0, 13.5, 8.2, 9.1, 7.4, 10.9, 7.9, 9.0, 11.4, 9.3, 8.9 and 9.7 dB, band mean 9.29 dB (6.3 dB if the printed power is the total over both sides of the duct, which is not settled). The likely reason is that the row is a series pressure source on a slit whose small-signal impedance is a few times the bore's, so the launched wave is only part of the source (not measured separately). The test that measured it (`noise_power_tests`, band mean within 1.5 dB) fails by that amount and is not part of this change; it returns when the print's definition (one- or two-sided power, a dipole force or a monopole volume source) is settled. No row, value or engine path changes.

### Reed residual flow as a card row

Almeida, Vergez and Causse 2007 (JASA 121(1) 536-546), Sec. III.1 and Fig. 9, measure a natural cane oboe reed on an
artificial mouth: the flow "first stabilises at a certain minimum value" instead of vanishing at the beating pressure
("very hard to completely close the reed"; "narrow channels remaining between the two blades ... behaving like rigid
capillary ducts"), and Fig. 10 shows a 1/2 power dependence of that residual flow. The engine's reed flow went to exactly zero
at zero opening. The reed card now accepts `residual_flow AREA_M2` right after `closure_stop` (or `contraction`,
`staple_recovery`, `slit_length`, the reed row), before `tongue_contact`, for a massive reed on the inviscid law (refused with
`slit_length` and with `staple_recovery`). The jet area becomes `C_c w max(h, 0) + AREA` (`C_c` the `contraction` value, else 1): the
voice hands the junction the opening `max(h, 0) + AREA / (C_c w)` through the existing flow-opening argument of
`step_with_junction`, so no kernel and no function signature changes; the blade's mechanics, lay contact and closure stop keep
the blade opening `h` alone. The row has no default and carries a measured quantity: one reed's Fig. 9 read at 450 dpi gives
0.094 to 0.100 mm^2 (mean 0.097e-6 m^2, 22 to 27 ml/s at 32 to 47 kPa), about 3 % of that reed's fitted rest area, and it is
one reed's. A card without the row renders as before, bit for bit (the flow-opening argument is unchanged when the voice has no
residual). `ReedBoreVoice::set_residual_flow_area` and `residual_flow_area_m2` are the voice-level switch. Tests (G0 to G2):
`render::schedule::reed::residual_flow_tests`.

`AREA` is the contracted JET area, with discharge coefficient 1: a shut reed passes `AREA sqrt(2 dp / rho)` and the junction sees the floor opening `AREA / (C_c w)`, 27.3 um for 9.33e-8 m^2 on a 5.6 mm reed with `C_c` 0.611 (`AREA / w` would be 16.7 um). A shut reed with no row takes the junction's exact closed-valve branch (`characteristic_junction_with_recovery`, `a == 0`); any positive `AREA`, however small, takes the general branch, which reproduces it to rounding (about 1e-16 relative). So a vanishing area agrees with no row to rounding and not bit for bit: `g1` pins it at 1e-20 (the pressure within 1e-9 of the peak), and on the Bb4 oboe card `residual_flow 1e-300` moves the first closure by 1.5e-21 m and leaves the clip level equal to five digits (3185.7 LSB, rms of 1 to 2 s). Measured on that card (rest jet area 1.222 mm^2, blow 7.4 kPa, rho 1.2 assumed): 9.33e-8 m^2 is 7.64 % of the rest jet area, 11 to 12 % of the mean flow, and cancels 30 % of the reed's quasi-static negative slope (reed -2.317e-9, leak +7.0e-10 m^3/s/Pa). It raises the oscillation threshold. Clip level (rms of 1 to 2 s) in dB against AREA 0: heavy-reed Bb4 at 7.4 kPa 0, -0.25, -0.76, -1.57, -2.45 for AREA 0, 1e-8, 3e-8, 6e-8, 9.33e-8; F4 heavy reed at 6.865 kPa 0, -0.20, -0.56 and silent at 9.33e-8; light reed (7.69e-6 kg, damping ratio 0.32) Bb4 at 7.4 kPa 0, -0.21, -0.66, -1.03 (4.5e-8) and silent at 6e-8, where the same card sounds at blows of 8.5, 9.5 and 10.5 kPa (-0.36, -0.14, -0.97 dB): the leak raises the blowing pressure the reed needs to start, and a card that holds a leak the print names says so.

### Reed closure stop as a card row

Commit a857c6f17 made the reed step's rigid closure stop at zero opening unconditional (a reed that the free implicit
midpoint solve would carry below zero opening ends the step at exactly zero with end velocity 0; no restitution is printed
for a reed on its lay). The reed card now accepts `closure_stop on|off` right after `contraction` (or `staple_recovery`,
`slit_length`, the reed row), before `tongue_contact`, for a massive reed. `on` is the default and a card without the row
renders as before, bit for bit (`step_with_contact_response` still passes the stop). `off` selects
`OpeningPropagator::MidpointFree`, the same implicit midpoint step without the stop, so a reed that goes below zero keeps that
opening and the lay contact alone closes it; a step whose free solution stays at or above zero is the same step either way.
The row carries no value: it is an A/B selection, not a physical claim. `ReedBoreVoice::set_closure_stop` and `closure_stop`
are the voice-level switch (massive reed, voice not advanced). The loaded finite-gas path keeps passing `false` as before.
Tests (G0, G1): `render::schedule::reed::closure_stop_tests`.

### A prescribed-flow step of the DG port load

`DgPortLoad::step(volume_flow_m3_s)` advances the load alone by one sample
under a prescribed bore volume flow: the crate-private trial of
`preview_step`, accepted at once. It adds no numerics and no state; a refusal
publishes nothing, as for the trial. It is the public counterpart of
`step_with_mouth_flow` for a load with no driven mouth, and it exists so that
a test outside the crate can march the load as its own reference:
`tests/dg_finite_gas/molecular.rs` called the crate-private trial and the
`dg_finite_gas` target did not compile from 2026-10-29 until 2026-11-10. No
claim: unlike `step_with_mouth_flow` it takes no cancellation gate, and it is
not a coupling path (nonlinear participants keep using the trial).

### Thick isotropic wall: Lame above thickness/radius 0.05

`LocalMaterialWall::base_patch` (graph, DG and acoustic-cell walls all pass through it) keeps the
membrane expressions at thickness/radius <= 0.05, bit for bit. Above that an ISOTROPIC wall takes
Lame's thick cylinder (Timoshenko and Goodier, thick-walled cylinder under internal pressure; free
outer face; plane stress for a stress-free axial restraint, plane strain for a strain-fixed one):
the stiffness per inner area is `p / u(a)` and the surface density is the kinetic energy of the same
radial field `u = A r + B / r` per inner area (a Rayleigh estimate, below the through-thickness
resonances). The Poisson ratio is the wall's own resolved material's; the Kelvin-Voigt resistance and
Maxwell arms scale with `stiffness / E` as before. An ORTHOTROPIC wall over the limit is still
refused. At the gate for nu = 1/3 the two laws differ by 4.06 % in stiffness and 0.68 % in mass
(stress-free), 4.82 % and 0.12 % (strain-fixed); the tests carry those amounts in their names.
This replaces the refusal `local membrane wall requires positive geometry and thickness/radius <= 0.05`
for isotropic walls (test `g3_the_same_input_over_a_twentieth_of_the_radius_now_takes_the_thick_wall_law`).

### Adaptive laminar enthalpy and external hydraulic continuation

`LaminarEnthalpy::advance_windkessel_adaptive` reuses the original coupled
BDF1 coarse/two-fine controller. Coarse and first-half solves borrow the same
accepted flow/compliance history; the second half borrows the actual first-half
pair. Rejected comparisons advance no external state. The existing ordinary,
Boussinesq and driven entries retain their original semantics and leave the
new optional hydraulic observations empty. No second time controller or
extrapolated physical field is introduced.

An explicit hydraulic norm joins the original caloric, velocity and work
criteria. It checks each port's compliance and added pressure, actual integrated
native Q, separate C/L storage, and signed plus resolved-magnitude resistor,
reservoir and port work. Signed and magnitude Q*dt*depth throughput have their
own cubic-metre tolerance, so equal endpoint flow cannot conceal different
transported volume. Pressure scaling removes the held reservoir offset and
uses the start compliance/Ohmic state; it invents no initial acceleration.
Storage changes and errors use difference of squares before dimensional
multiplication. All spans and local corrections come from the original
first-order controller, including exact target clipping and one-ULP cap tails.

Only accepted fine receipts contribute compensated per-port work and signed/
absolute transported volume totals. Initial Q is independently obtained from
the original accepted BDM edge transport; initial/final C*x²/2 and L*Q²/2
changes use the paired physical states and declared depth. Each port and their
aggregate must satisfy port work minus reservoir work, both resistor losses,
both storage changes and the separate native BDF defect within the explicit
window allowance. This external-circuit identity does not add heat, replace
the independent native caloric balance, or certify whole-system passivity.

All coarse/rejected work consumes the shared native-step and thermal-endpoint
ceilings. Source, material, momentum, stream and thermal refusals remain
terminal. Cancellation or a failed final circuit account rejects the complete
window, including any internally accepted prefix. The final step's retained
compliance coordinates continue with the returned original velocity and clock.
Circuit parameters/references and physical boundary data remain held; this
entry adds no scheduled circuit events or Boussinesq/load combination.

Native regression definitions cover symmetric compliance discharge with zero
fluid/thermal change, hydraulic-only refinement and an independent accepted BE
recurrence; signed equilibrium all-port plug transport and exact work/storage;
nonzero resistor losses excluded from heat; decimal target completion; and
late source, work-budget and cancellation rollback after an accepted prefix.
All seven changed Rust files passed the pinned parser, formatter and
`rustfmt --check`. Bounded independent source review found no paired-history
or accounting blocker. Full native FEM/thermal compilation and these native
integration test executions remain unverified in the current partial restored
environment; no runtime pass is inferred from those checks.

## Parts

| Part | Bytes | `##` sections | First (newest) section | Last (oldest) section |
|---|---|---|---|---|
| [`contracts/part-08.md`](contracts/part-08.md) | 449466 | 145 | Driven full-curve Hertz friction with original thermal sources | CST/DKT hereditary material memory on the original modal basis |
| [`contracts/part-07.md`](contracts/part-07.md) | 442691 | 141 | Thermal porous follower pressure through enthalpy feedback | Electrical contact heat in thermally live actuator supports |
| [`contracts/part-06.md`](contracts/part-06.md) | 446940 | 130 | Complete material and acoustic memory in exponential impact | Shared material cavitation, films and bubble acoustic history |
| [`contracts/part-05.md`](contracts/part-05.md) | 420894 | 27 | Finite solvent heat during gas absorption and exsolution | Conformance tests |
| [`contracts/part-04.md`](contracts/part-04.md) | 449410 | 65 | No-claim boundaries | Live gas changes through the shared acoustic port |
| [`contracts/part-03.md`](contracts/part-03.md) | 444523 | 95 | Finite thermochemical feedback through the complete acoustic load | Nonlinear spatial current with complete material heat exchange |
| [`contracts/part-02.md`](contracts/part-02.md) | 425077 | 113 | Source-resolved Hertz islands on shared finite thermal cells — v0 | Complete sliding-film layer heat on the impact clock |
| [`contracts/part-01.md`](contracts/part-01.md) | 285661 | 57 | Generalized-Newtonian film resistance in impact mechanics | Native joint friction with finite-cell thermal feedback |

## Purpose and layer

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Public types and semantics

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Invariants

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Error model

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Determinism class

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Cancellation behavior

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Unsafe boundary

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Feature flags

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## Conformance tests

Full text: [`contracts/part-05.md`](contracts/part-05.md), under the same heading.

## No-claim boundaries

Full text: [`contracts/part-04.md`](contracts/part-04.md), under the same heading.
