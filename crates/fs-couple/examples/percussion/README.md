# Physical percussion reference: sourced anchors, editable estimates

Run from the workspace root:

```sh
cargo run -p fs-couple --example percussion -- splash 4096 > splash.csv
cargo run -p fs-couple --example percussion -- drum 4096 > drum.csv
```

The commands above are **mechanical reference experiments**. The additional
`splash-wav` and `drum-wav` commands in [AUDIO.md](AUDIO.md) connect these mechanics
to actual two-sided/closed-boundary BEM transfers and fixed-receiver pressure.
Neither path is a finished sampled instrument,
not measured digital twins and not a real-time-qualified plugin. The executable
uses the original fs-phs nonlinear stepper, which allocates and builds a dense
finite-difference Newton matrix. Native compilation and execution were not
available during implementation. Numerical Python references do not establish
that this Rust example builds, completes its modal solve, or runs in real time.

The generated CSV contains physical motion, *internal* cavity pressure and actual
energy. It deliberately does not label displacement or interior pressure as a
microphone waveform. A cymbal is unbaffled: its signed integrated surface motion
cannot simply be substituted into a compact baffled-monopole sound formula.

## Research anchors (checked 19 September 2026)

| Reference object | Published information used | Still an estimate or missing |
|---|---|---|
| Zildjian splash studied by Kaselouris et al., Acoustics 2023, 5, 165–176 | 203.2 mm outer diameter; 78 mm bell diameter; 12.3 mm mounting hole; 0.5 mm edge thickness. The study adopts B20 density 8607 kg/m³, E=112.6 GPa, nu=0.342 from literature. | Full meridian height and thickness map, hammer coordinates/depths, lathe grooves, specimen mass, residual forming stress, local hardness, measured damping. The executable's interior profile is labelled estimated, not extracted from a figure. |
| Zildjian Z5A | 406.4 mm length, 0.560 inch diameter (14.224 mm), hickory, medium taper, oval wooden tip. | The example's radius stations, 800 kg/m³ density, grip pivot, tip curvature and transverse contact modulus. The paper's 394 mm stick is a DIFFERENT stick; it is not substituted for Z5A. |
| Pearl Masters Maple MM6 | Selected nominal 14 × 6.5 inch size, six-ply 7.5 mm maple shell, 45-degree bearing edge. | The present drum example uses a rigid cylindrical enclosure, not elastic maple plies, actual bearing-edge geometry, hoops or lug compliance. Clear vibrating radius is estimated from nominal radius minus shell thickness. |
| Remo Ambassador Clear and Ambassador Hazy Snare Side | Single-ply 10 mil clear film and 3 mil snare-side film: 0.254 mm and 0.0762 mm. | Polymer constants E=4 GPa, nu=.38, density=1390 kg/m³, and installed tensions 3000/1500 N/m are declared research inputs, not manufacturer measurements. |
| Stand felts | Wool-felt compression/hysteresis is represented by the existing fs-material law. | Zildjian stand-felt dimensions and a dynamic force–compression curve were not found in the reviewed manufacturer information. OD30/ID13 mm, 6 mm thickness, preload and all crush/recovery constants are estimates. No piano-hammer parameters are claimed as measured stand-felt data. |

Primary sources:

* Kaselouris et al., **FEM-BEM Vibroacoustic Simulations of Motion Driven Cymbal-Drumstick Interactions**, doi:10.3390/acoustics5010010, https://www.mdpi.com/2624-599X/5/1/10
* Zildjian, **5A Drumsticks (Z5A)**, https://zildjian.com/products/5a-drumsticks
* Pearl, **Masters Maple**, https://pearldrum.com/global/products/snares/drum-set-series-snare-drums/masters-maple
* Remo, **Ambassador Clear**, https://remo.com/partnumber/ba-0312-00 ; **14-inch Ambassador Hazy Snare Side**, https://remo.com/partnumber/sa-0114-00
* Copper Development Association interview with Zildjian's Paul Francis, **The Legacy of Zildjian Cymbals Signature Sound Lives On** (80% copper/20% tin description; manufacturing remains proprietary), https://copper.org/consumers/arts/2012/november/Legacy_Zildjian_Cymbals_Signature_Sound_Lives_On.php
* Daniel Russell, Penn State, **The piano hammer as a nonlinear spring**, https://www.acs.psu.edu/drussell/Piano/NonlinearHammer.html — qualitative felt hysteresis and why static curves do not identify dynamic impact. Not a stand-felt calibration.

## What actually runs

The splash path constructs a free curved, tapered shell, reduces its numerical
modes, and projects **its own** Green–Lagrange membrane strain. No rectangular
plate coefficients are borrowed. The first example retains vertical rigid
translation plus elastic modes in an explicit 50–1200 Hz window. It omits the
other rigid coordinates and is not a full rocking/sliding stand model. No
continuum accuracy or sufficient crash bandwidth is implied by that window.

Six opposing felt patches reuse WoolFelt's loading, unloading and crush state.
Their stored recovery energy and rate-dependent Kelvin deformation participate
in the same implicit solve as the shell and striker. Irreversible conditioning
loss is recorded on history acceptance; a refused sample changes neither body
motion nor felt history. The initial symmetric preload is a declared initial
condition, not a solved gravity/wingnut equilibrium.

The stick's mass, center of mass and pivot inertia come from integrating its
editable axial radius profile and declared density. The contact-coordinate
mass is I/lever². It is then a free effective mass, not an oscillator tethered
by a hidden spring. The elastic Hertz tip approximation feeds fs-dcontact;
impacts and rebound arise from contact and initial motion, not a force envelope.
Wood flexure, anisotropic tip compliance, hand forces and contact damping are
not inferred. Nonzero Hunt–Crossley loss is rejected by this particular host
rather than silently ignored.

The drum path creates two tensioned films through the existing DKT/prestress
pencils and couples their signed swept volumes through one sealed-air
compliance. Its initially resting second head is driven only by that coupling.
It is **not yet a snare**: no wire bed/rattle or vent is represented, and the
heads remain linear within their retained basis. Timpani additionally need
bowl/air eigenmodes and radiation; a volume spring alone does not reproduce them.

All reductions and contacts are generic. A gong changes the supplied shell
profile/material/supports; a membrane changes actual radius, film and tension.
Measured hammer and lathe detail can enter `profile::revolve` as explicit local
geometry/thickness relief; unresolved features are reported. No random hammer
pattern is passed off as an exact Zildjian manufacturing map.

## The remaining path to realistic real-time audio

The immediate accuracy gate is a measured meridian/thickness/hammer survey,
impact and felt load–unload/rate tests, plus multi-position impulse responses.
Those identify uncertainties that alloy labels alone cannot determine. Retain
in-plane and high-frequency modes until nonlinear energy transfer and radiated
pressure converge. Match damping and contact only within held-out force/response
tests, not by matching one attractive sound.

The runtime gate is a bounded, workspace-based fast realization of the existing
discrete-gradient law, benchmarked at actual sample rates, mode counts, contact
counts and polyphony. It must reproduce this reference before being advertised
as real time. A first one-way two-sided/closed-boundary BEM and fixed-receiver
filter path now exists (AUDIO.md), but native execution, boundary refinement,
full-band accuracy and radiation loading remain open. A moving receiver should
use the existing broadband directional facility rather than this fixed bank. A room response is downstream, not a substitute for source
physics. Finite-amplitude membrane/shell validity and aliasing must be checked
separately. None of those missing gates is promoted by the tests here.

## Supply the actual shell specimen

`--shell-profile PATH` replaces the hard-coded splash geometry, isotropic material
and elastic frequency window before shell assembly. It works with `splash`,
`splash-wav` and `splash-mic`, including `--prepared-nonlinear`. The existing shell
FEM, nonlinear reduction, two-sided acoustic surface and BEM microphone path all
consume that same geometry and thickness; nothing is converted into an output EQ.

```sh
(set -C; cargo run --release -p fs-couple --example percussion -- \
  splash 4096 --prepared-nonlinear \
  --shell-profile crates/fs-couple/examples/percussion/estimated-splash.profile \
  --strike-position-m 0.067 0 > supplied-shell.csv)
```

The distributed file reproduces the **estimated** reference meridian and material,
not a new measurement. Replace its records with specimen data. The strict text
format begins with `frankensim-shell-profile-v1`, accepts `#` comments, and uses
comma-separated numeric records. `material,E_Pa,nu,rho_kg_m3`, `azimuths,N` and
`band_hz,lower,upper` are each mandatory exactly once. At least two
`station,radius_m,height_m,thickness_m` rows describe the meridian in strictly
increasing radius. Zero inner radius makes a disk; nonzero radius makes a hole.
All heights are reference midsurface coordinates, all thicknesses are physical.
There is no implicit unit conversion, random detail, material-name lookup or
fallback to the example when input is invalid.

Optional `ring,radius_m,half_width_m,height_delta_m,thickness_delta_m` records use
the existing annular-relief law. Optional
`dent,center_x_m,center_y_m,radius_m,height_delta_m,thickness_delta_m` records use
its compact C1 indentation law. Both alter actual mechanical geometry and section
stiffness/mass, not a rendering normal map. They describe a stress-free reference
shape, **not** simulated manufacturing, work hardening or residual forming stress.
The existing mesh/work budgets remain in force. Nonintersecting detail and details
flagged underresolved by the geometry owner refuse; add meridian stations and
azimuths rather than treating invisible detail as simulated. Passing this coarse
resolution screen is not convergence certification.

The shell can have different diameter, bell, taper, relief and material, but this
command still supplies the same **estimated stand, stick and Hertz-contact model**.
Its pad centres are at 12 mm radius and are located on the supplied surface, never
snapped across the hole. Geometry that does not support those pads refuses.
Without an explicit strike position, supplied profiles use the point two-thirds
of the way from inner to outer radius on positive x. The no-file reference retains
its original contact/pad sampling. Thus identical profile geometry alone does not
assert bit-identical trajectories between the two hardware sampling choices.
A differently suspended gong still needs its own support model; this importer
does not silently make the stand appropriate for every shell.

The window remains limited to 31 elastic modes plus vertical translation and to
the mechanical clock's Nyquist guard. Increasing its upper edge is not proof of
an adequate in-plane basis or nonlinear bandwidth. The fixed-receiver acoustic
bake remains 40–1640 Hz. Large scans, sharp features or broader frequency windows
can exceed the existing budgets and should refuse, not truncate silently. An
input file carries no automatic measured/calibrated status. Native execution and
measured response comparisons remain required for instrument-fidelity claims.
