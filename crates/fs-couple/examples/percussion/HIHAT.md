# Two cymbals, physical closing and reopening

`hihat`, `hihat-wav` and `hihat-mic` construct two independently reduced curved
shells and one force-driven axial carriage in the existing nonlinear impact
owner. Distributed inter-cymbal contacts, both shell potentials, stand felts,
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
- `strike_m,x,y` supplies the first stick's default station. Existing explicit
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
checks. Stick stations use that same material chart on the upper positive face.

With local outward surface displacements u and l, the signed point gap is
`D + z_upper_skin + z_lower_skin + u + l`. The owner receives closure row
`[-b_upper,-b_lower]`, so its forces are equal and opposite in the world frame.
Neither sticks nor pedal receive direct inter-cymbal reaction coefficients.
Pedal force reaches the upper cymbal only through the opposed felt mount,
which supplies reciprocal carriage/shell reactions. The lower mount is grounded.
All original mechanical source addresses, contact losses, Kelvin histories,
energy gates, and accepted-step rollback remain in the same time owner.

The shared shell preparation retains each complete requested elastic slice
plus its true geometry-derived vertical translation. Original mode, panel,
material-validity and solver limits are unchanged; overflow never drops modes.
CSV reports upper/downward and lower/downward point motion, pedal travel/speed,
minimum sampled gap, active sites and the accepted total work/energy/loss ledger.
Gap/activity are endpoint observations, not averaged contact reactions.

## Existing controls and limits

Both existing stick-force files and two-stick position/speed controls work.
`--prepared-nonlinear`, `--analytic-newton`, `--impact-substeps`,
`--microphone-right` (finite microphone only), and `--radiation-spec` retain
their existing meanings. A cancelled/refused tick consumes no pedal or stick
input. Other drum/single-cymbal options, including mallets and compliant mutes,
are not admitted here; they are not silently ignored or converted.

Both skins are placed in ONE stationary BEM solve, retaining cross-body
scattering at the separated reference geometry. Source velocities are not
sign-flipped during the proper lower-body rotation. Pedal/sticks do not gain
invented radiation channels; the actual resulting shell motion is observed.
The complete reference skins must pass a conservative separating-plane check.

This is an axial, fixed-reference, small-slope image. It does NOT include full
stand rocking, tangential friction, changing contact normals/pairings, pedal
lever geometry, gravity, air squeeze-film loading, or acoustics recomputed as
the gap closes. Sampling stations do not certify absence of interpenetration
between stations. Stationary one-way BEM during a closing gesture is explicitly
an approximation, not a gap-dependent scattering or full-band hi-hat claim.
Material, spatial, modal and temporal convergence remain separate obligations.

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
