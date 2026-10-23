# Two cymbals, physical closing and reopening

`hihat`, `hihat-wav` and `hihat-mic` construct two independently reduced curved
shells and one force-driven axial carriage in the existing nonlinear impact
owner. Optional bounded cymbal rocking adds real angular inertia. Distributed
inter-cymbal contacts, both shell potentials, stand felts,
carriage inertia/return spring, and stick contacts are solved together. This
is a new physical instrument composition, not a crossfade between open/closed
recordings or two separately rendered cymbals mixed after impact.

```sh
# Illustrative eight-inch pair, NOT a measured commercial hi-hat.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  hihat crates/fs-couple/examples/percussion/estimated-hihat.fshh 12000 \
  --analytic-newton --impact-substeps 8 511 --strike-speed-m-s 0 \
  > pedal-cymbals.csv)

# Both shells in one reference BEM scene, two independent receivers.
# This can require expensive offline preparation and can refuse a poor fit.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  hihat-mic crates/fs-couple/examples/percussion/estimated-hihat.fshh 4800 20 \
  --analytic-newton --impact-substeps 8 511 --strike-speed-m-s 0.8 \
  --second-stick-position-m -0.06 0.01 --second-stick-speed-m-s 0.6 \
  --microphone-right -0.08,0.05,0.35 > paired-cymbals.wav)
```

Command comes first. Positional forms are `hihat INPUT [steps]`,
`hihat-wav INPUT [frames] [full_scale_Pa]`, and
`hihat-mic INPUT [frames] [full_scale_Pa] [x_m y_m z_m]`.
The input's entire pedal force program must fit in the requested duration.
The supplied program ends at 0.020 s; the 12000-step example covers 0.024 s.
These are supported invocations, not a claimed completed render or calibration.

## Complete input, no inferred instrument parameters

The 16 KiB UTF-8 format starts with `frankensim-hihat-v1`. Blank lines and `#`
comments are allowed. The included `.fshh` supplies every required record:

- `shell,upper|lower,profile|mesh,path` selects each source independently; paths
  are relative to the `.fshh` file. The existing profile or explicit 3D `.fss`
  reader retains all geometry, nodal thickness and material data. Supplying the
  same path twice explicitly creates two physical bodies with separate states.
- `separation_m,D` places the upper local origin at world z=+D/2. The lower is
  rotated 180 degrees about x and placed at -D/2. This is a proper rotation,
  not an inverted mesh. Both local +z axes point outward from the pair.
- `carriage,mass_kg,return_N_m,drag_Ns_m` defines one true translating mass with
  a grounded linear return spring and passive drag. Its equilibrium travel is
  zero. Positive travel/force is downward, closing the upper cymbal.
- `damping_ratio,upper,lower` supplies physical modal loss on each nonzero
  structural mode. No loss is manufactured from the material's elastic label.
- `contact,K_N_per_m_alpha,alpha,chi_s_m` binds directly to the existing
  distributed power-law contact and Hunt-Crossley loss owner.
- `site,x_m,y_m,fraction` supplies 1..64 distinct reference-XY rim-contact
  stations. Fractions are positive and sum to one; they partition the total
  effective contact stiffness K, NOT independent full-strength duplicates.
- `ambient,<gas record>` declares the air around both shells; a shell file that
  declares its own `ambient` row must agree with it.
- `strike_m,x,y` supplies the first stick's default world-XY station. Existing explicit
  strike position and speed controls override it; first launch defaults to
  0.8 m/s unless zero is selected. A second stick requires its own position.
- `mount,upper|lower,radius_m,face_area_m2,thickness_m,precompression_m,
  f_ref_Pa,eps_ref,p,q,crush_fraction,eps_densify,prior_strain,K_N_m,eta_Ns_m`
  supplies each opposed felt mount. Each face has three equal-area sites at
  the supplied radius. Its Kelvin K and eta are whole-face values, divided by
  three together. The original WoolFelt law and separate patch histories apply.
- `pedal,time_s,force_N` is a complete signed external carriage-force program.
  It uses the existing force player's linear interpolation and exact tick
  integration: strictly increasing nonnegative times, first/last force zero,
  and zero force outside the program. Negative force is active retraction;
  zero force leaves carriage motion and the passive return spring intact.

All singleton records occur once (each side once for shells/mounts). Missing,
unknown, duplicate, nonfinite, invalid or over-budget records refuse. The
estimated mounts/contact/forces are editable research inputs, not measured
pedal hardware, bronze restitution, or manufacturer dimensions.

## Contact geometry and reciprocal work

Contacts use the actual thickness-offset NEGATIVE skin of each shell. The
same area-weighted directors create the two-sided acoustic skin. Local point
velocity includes translation plus `theta cross director_offset`; its
transpose applies both force and moment to the shell. Barycentrics are located
on the inner skin's actual XY chart, including mounting holes and ambiguity
checks. Stick stations use the selected shell's positive skin's actual XY chart;
an inner-skin barycentric location is not reused at a different physical face.

With local outward surface displacements u and l, the signed point gap is
`D + z_upper_skin + z_lower_skin + u + l`. The owner receives closure row
`[-b_upper,-b_lower]`, so its forces are equal and opposite in the world frame.
Neither sticks nor pedal receive direct inter-cymbal reaction coefficients.
Pedal force reaches the upper cymbal only through the opposed felt mount,
which supplies reciprocal carriage/shell reactions. The lower mount is grounded.
All original mechanical source addresses, contact losses, Kelvin histories,
energy gates, and accepted-step rollback remain in the same time owner.

The shared shell preparation retains each complete requested elastic slice
plus its true geometry-derived vertical translation. Optional rocking adds two
geometry-derived angular coordinates per shell. Original mode, panel,
material-validity and solver limits are unchanged; overflow never drops modes.
CSV reports upper/downward and lower/downward point motion, pedal travel/speed,
minimum sampled gap, active sites and the accepted total work/energy/loss ledger.
Gap/activity are endpoint observations, not averaged contact reactions.

## Existing controls and limits

Both existing stick-force files and two-stick position/speed controls work.
`--prepared-nonlinear`, `--analytic-newton`, `--impact-substeps`,
`--microphone-right` (finite microphone only), and `--radiation-spec` retain
their existing meanings. A cancelled/refused tick consumes no pedal or stick
input. `--compliant-mute` admits an explicitly selected outer hi-hat surface
as described below. Other drum/single-cymbal options retain their admission rules.

### Cymbal rocking against the stand felts

`--cymbal-rocking MAXIMUM_ANGLE_RAD` enables small-angle rotation about the two
horizontal axes on **both** shells. The maximum rotation magnitude must be in
`(0,0.1]` radians; `--cymbal-rocking 0.03` is an example declared bound, not a
measured hi-hat limit. The default continues to retain vertical translation and
the complete elastic slice without these additional angular coordinates.

Each shell's own geometry and shell mass matrix determine its angular inertia.
Off-center strikes, the three-site opposed stand washers, inter-cymbal contacts,
and optional thin-gap air apply reciprocal tilting moments at their physical
stations. All twelve stand pads and their individual histories remain; the
carriage continues to move axially through its supplied spring and pedal force.

The option composes with independent targets, flexible sticks, felt mallets,
a moving choke, shell material memory and the existing acoustic image. The pair
gains four coordinates; all elastic modes must still fit the unchanged budget.
CSV adds `upper_rigid_rotation_{x,y,z}_rad` and
`lower_rigid_rotation_{x,y,z}_rad`, expressed as rotation vectors in the world
frame. A step outside either shell's declared angular bound refuses without
consuming the player's input or retaining a partial state.

See [ROCKING.md](ROCKING.md) for the constant-inertia, first-order pose contract
and its fixed contact/acoustic geometry limits. The paired regression exercises
both angular responses with gap air, independent hand/pedal forces, preserved
histories and rejected-step rollback; native execution remains required.

### Independent upper and lower exterior strikes

Each hand chooses its own outer surface with `--strike-target upper|lower` and
`--second-stick-target upper|lower`. Both default to `upper`, preserving existing
performances. A second target requires `--second-stick-position-m X Y`; it does
not create another hand implicitly. These flags apply to all three hi-hat output
forms. They select physical contact geometry and force direction, not a sound
layer or an acoustic gain.

```sh
cargo run --release -p fs-couple --example percussion -- \
  hihat crates/fs-couple/examples/percussion/estimated-hihat.fshh 12000 \
  --strike-target upper --strike-position-m 0.06 0.01 --strike-speed-m-s 0.5 \
  --second-stick-target lower --second-stick-position-m -0.06 0.01 \
  --second-stick-speed-m-s 0.6 --analytic-newton --impact-substeps 8 511 \
  > opposed-strikes.csv
```

All supplied XY positions are in the pair's world frame. Lower stations map to
local `(x,-y)` through the same proper rotation used by shell collision and
radiation. An upper exterior striker travels downward; a lower exterior striker
travels upward. For either hand, positive launch speed, displacement and force
mean **inward toward its selected shell**, and negative force retracts that hand.
The scalar contact penetration is `stick_inward + shell_outward`; the reaction
uses that same row's transpose. Contacts retain the selected outer skin's
thickness-offset moments, without changing either shell's acoustic source sign.

Wood tips, flexible shafts, and the existing v1/v2 felt mallets all follow the
selected target. The two hand-force files remain independent. Lower felt faces
use their complete outer-skin footprint and geometric clearances; the shaft
azimuth in a v2 card remains a **world-XY** direction, with each site's moment arm
properly transformed. All twelve stand contacts, the pedal, optional gas/material
state, moving choke, original source addresses and existing mode/pad limits stay
in the joint solve. Two felt mallets, or one felt mallet plus a choke, still exceed
the sixteen-pad limit and refuse.

The `upper_down_m` and `lower_down_m` probes remain world-downward. The upper
probe follows the first hand's station when that hand targets the upper shell;
when the first hand targets the lower shell, the upper probe uses the first
inter-cymbal contact station. The lower probe remains at that first contact
station. Lower-target flexible-hand displacement columns say `inward_m`; their
speed columns share the inward convention. Default-upper CSV names are unchanged.
Only exterior axial impacts are modeled: inner-face strikes, striker-to-striker
collision, and interaction with stand hardware require additional geometry and
contacts. Focused regressions cover reflected finite-thickness virtual work,
loaded-mallet moments and clearances, a lower wood impact with independent hand
forces, and a lower felt strike with gap air and rejected-tick rollback. Native
execution and acoustic renders remain required before sound or performance claims.

### Moving felt choke on either outer shell

`--compliant-mute file.fsm` adds one translating felt pad on `surface,upper`
with `jaw,above,...`, or on `surface,lower` with `jaw,below,...`. These targets
are explicit: an existing `surface,shell` card does not implicitly choose a
member of the pair. Only an exterior jaw is admitted; a pad between the closing
shells would need additional jaw/shell collision and displaced-air mechanics.

```sh
cargo run --release -p fs-couple --example percussion -- \
  hihat crates/fs-couple/examples/percussion/estimated-hihat.fshh 12000 \
  --strike-speed-m-s 0.5 --analytic-newton --impact-substeps 8 511 \
  --compliant-mute crates/fs-couple/examples/percussion/estimated-hihat-mute.fsm \
  > choked-pair.csv
```

The force program drives the pad's actual mass inward and later retracts it.
Each positive-area site contacts the selected outer skin through its physical
displacement and rotation port. World XY stations on the lower shell are
transformed through the same proper rotation used by the pair geometry. Felt
compression, release, creep and conditioning evolve in the joint time owner;
the resulting shell motion supplies the existing acoustic sources. See
[COMPLIANT_MUTE.md](COMPLIANT_MUTE.md) for the card and force conventions.

All twelve mount sites, the pedal, both shells and their collision remain.
The moving jaw precedes any flexible shaft coordinates, while source and pedal
addresses remain unchanged. One to four choke sites fit the existing sixteen-pad
limit. A four-site mallet plus a choke exceeds that limit and refuses before
shell preparation; no stand or contact site is removed. Optional shell material
memory, squeeze-film/gas storage, acoustic feedback and both flexible hands
continue through the same shared state and force clock. The three hi-hat output
forms accept this option; CSV also reports inward pad travel, velocity and
endpoint contact force.

The authored pad is an axial, conforming contact quadrature. It does not model
hand anatomy, tangential friction, pad acoustic scattering, or jaw contact with
stand/striker hardware. Site area is not a certification of a resolved footprint
boundary. Focused tests cover both exterior signs, actual thickness moments,
preserved source/material addresses, and a driven lower pad with retraction,
gap air, shaft flexure, pedal work and rejected-tick rollback. These tests still
require native execution; no calibrated choke duration or sound is claimed.

### Felt mallets and loaded shafts

Either hand can replace its wood tip with `--mallet-spec` or
`--second-mallet-spec` on all three hi-hat output forms. The selected hand
requires an explicit strike position. A v1 card supplies its effective moving
mass; a v2 card requires the same hand's flexible shaft and loads that shaft
with the head's actual mass and rotary inertia before reduction. See
[MALLETS.md](MALLETS.md) and [LOADED_MALLETS.md](LOADED_MALLETS.md).

```sh
cargo run --release -p fs-couple --example percussion -- \
  hihat crates/fs-couple/examples/percussion/estimated-hihat.fshh 12000 \
  --strike-target lower --strike-position-m 0.06 0.01 --strike-speed-m-s 0.5 \
  --mallet-spec crates/fs-couple/examples/percussion/estimated-loaded-mallet.fsmallet \
  --flexible-stick crates/fs-couple/examples/percussion/estimated-flexible-stick.fst \
  --second-stick-position-m -0.06 0.01 --second-stick-speed-m-s 0.3 \
  --analytic-newton --impact-substeps 8 511
```

The four-site felt face uses the selected shell's actual positive skin and the
same footprint-derived clearances as a single cymbal. It replaces that hand's
Hertz contact. A loaded shaft also retains each site's reciprocal bending
moment, with physical hand forcing at the supplied hand station. The other
hand may remain a wood tip with or without shaft flexure.

All twelve opposed mount sites remain first in the material-history layout;
the four mallet sites follow them. **One felt mallet fits the existing
sixteen-pad limit. Two felt mallets refuse before shell preparation.** Neither
stand sites nor footprint samples are dropped to make a selection fit.
Inter-cymbal contact, pedal work, optional squeeze-film/gas storage, radiation
feedback, and both acoustic shell sources retain their existing composition.
Mallet or shaft coordinates are not additional acoustic sources. Failed ticks
preserve all material histories and the shared pedal/hand clock.

Focused regressions cover either-hand admission, the complete second-mallet
layout with compressible gap air, a loaded felt strike with two flexible hands
and pedal forcing, and absence of a parallel hard tip before physical contact.
Native execution and full renders are still required; these invocations do
not establish calibrated sound or real-time performance.

Both skins are placed in ONE stationary BEM solve, retaining cross-body
scattering at the separated reference geometry. Source velocities are not
sign-flipped during the proper lower-body rotation. Pedal/sticks do not gain
invented radiation channels; the actual resulting shell motion is observed.
The complete reference skins must pass a conservative separating-plane check.

This is a fixed-reference, small-slope image with axial carriage motion and
optional bounded cymbal tilt. It does NOT include motion of the stand body,
tangential friction, changing contact normals/pairings, pedal
lever geometry, gravity, or acoustics recomputed as
the gap closes. Sampling stations do not certify absence of interpenetration
between stations. Stationary one-way BEM during a closing gesture is explicitly
an approximation, not a gap-dependent scattering or full-band hi-hat claim.
Material, spatial, modal and temporal convergence remain separate obligations.

Optional `--squeeze-film INPUT.fsf` adds reciprocal, geometry-derived thin-gap
air resistance inside the same mechanical solve. See [SQUEEZE_FILM.md](SQUEEZE_FILM.md)
for its required inputs and incompressible applicability limits. Trapped pockets
require a different storage model and refuse; the flag does not make the BEM
scene move with the closing gap.

Focused native checks:

```sh
cargo test --release -p fs-plate --lib shell::reduction::radiation
cargo test --release -p fs-couple --example percussion hihat -- --test-threads=1
cargo test --release -p fs-couple --example percussion paired_scene -- --test-threads=1
```

The time fixture retains the complete requested elastic slices, true translation
masses, and washer histories. Construction and virtual-work regressions also
exercise different shell materials and both skin projections. Authored tests require native execution; no acoustic listening or
real-time performance result follows from construction alone.

### Native acoustic bake limit at 0.8 mm separation (`frankensim-nhf66`)

Native execution of `hihat-mic crates/fs-couple/examples/percussion/estimated-hihat.fshh 24000 50 0.08 0.05 0.35 --microphone-right -0.08,0.05,0.35 --analytic-newton --impact-substeps 8 511 --strike-speed-m-s 0.8` selects `PlainCbie` across `40..1640 Hz` (`2048` panels, `panels_per_wavelength = 16.2824` at `1640 Hz`, `condition_lower_bound = 3596.634`) and refuses at `1640 Hz, input 18` (upper-shell mode 18 at `892.43 Hz`):

- `surface_power = (-5.226533e-10, -2.645980e-7) W` (`power_interval = (-5.226532984933255e-10, -5.226532984929689e-10)`)
- `positive_real_power = +2.003103e-9 W`, `negative_real_power = -2.525756e-9 W`
- `apparent_power = 6.344901e-7 VA`, `real / apparent = -8.237375e-4` (`-0.0824%`)
- `plane_wave_reference_power = 3.459663e-7 W`, dense-solve roundoff floor `n * eps * cond_lb * apparent = 1.037747e-15 W`

Because the `0.8 mm` inter-cymbal rim separation is 30x smaller than the `24.7 mm` outer facet edges and the proper lower-shell rotation `(x, y, z) -> (x, -y, -z)` crosses quad triangulation diagonals across the gap, P0 centroid collocation (`PlainCbie`) incurs a `-8.24e-4` non-symmetric cross-gap resistance error on sub-coincident multipoles even with analytical static and `O(k^2 r^2)` double-layer regularization. This exceeds the LU roundoff floor by `5e5x` and is rightly refused by the passivity gate; the mechanical CSV (`hihat`) runs unaffected.

