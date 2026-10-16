# Finite-inertia felt mounting on the actual shell skins

`--mount-spec input.fsmount` replaces the six fixed legacy stand patches with
file-supplied upper/lower felt supports moving on **one shared carrier**. The
carrier has its own mass and explicit passive ground spring/resistance. It
receives the instrument's reaction in the same nonlinear solve: it is not an
externally prescribed displacement or an audio filter. Each site retains its
own WoolFelt conditioning and optional series Kelvin recovery state.

```sh
# Mechanics, including physical support motion and total energy/loss columns.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  splash 4096 --analytic-newton --impact-substeps 8 511 \
  --cymbal-rocking 0.03 --strike-position-m 0.068 0.01 \
  --shell-profile crates/fs-couple/examples/percussion/estimated-splash.profile \
  --mount-spec crates/fs-couple/examples/percussion/estimated-splash-mount.fsmount \
  > mounted-shell.csv)

# Same mechanics, observed through the existing physical microphone/BEM path.
# Acoustic fitting, passivity and mechanical admission still may refuse.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  splash-mic 48000 20 --analytic-newton --impact-substeps 8 511 \
  --cymbal-rocking 0.03 --strike-position-m 0.068 0.01 \
  --mount-spec crates/fs-couple/examples/percussion/estimated-splash-mount.fsmount \
  > mounted-shell.wav)
```

These are invocation examples, not a report that the native renderer has run or
that its BEM fit, sound quality or real-time deadline has passed. The supplied
card is **illustrative**. Its 0.4 kg carrier, ground coefficients, backing planes,
loaded-site locations and felt data are not measured Zildjian hardware. No
manufacturer-specific material or geometry identification is asserted.

## One geometry and one power-conjugate contact row

Each specified XY station is located independently on the actual positive and
negative thickness-offset skins. The upper/lower designation comes from their
physical heights, not a guessed triangle orientation. Skin displacement uses
the existing `u + theta cross offset` row, including the contact moment. BEM
radiation uses this same thickness/director definition. A dimple, thickness map
or change of bell shape therefore changes the actual gap and force coupling,
not just a render normal or a post-processing gain.

With carrier displacement `X` positive in the shell's reference Z direction,
reference skin heights `z_top,z_bottom`, backing planes `b_top,b_bottom`, and
felt thicknesses `h_top,h_bottom`, the geometric compressions are

```text
upper: h_top    - (b_top    - z_top)    + u_top    - X
lower: h_bottom - (z_bottom - b_bottom) - u_bottom + X
```

The material owner subtracts each local series creep deformation before
computing stress. Multiplying a patch force by the transpose of that same row
returns equal/opposite shell/carrier work. Both backings travel together; they
are not two separately moving mute jaws, and the carrier mass is counted once.
The original stick, shell, second-stick, mute and flexible-shaft coordinates
stay in place. The carrier is appended after them, before private material
states. Contact rows, force vectors and observations are extended with zero;
existing source-to-microphone addresses are unchanged. The carrier is not
invented as an additional acoustic radiator. It changes sound through its
reaction on the already-radiating shell.

## Explicit SI input

The first noncomment line is `frankensim-felt-mount-v1`. Records are comma
separated, with `#` comments permitted. `carrier`, `upper`, `lower` and `felt`
are each required exactly once. There are no inferred material defaults.

```text
carrier,mass_kg,ground_stiffness_N_per_m,ground_resistance_Ns_per_m,initial_Z_m,initial_velocity_m_per_s
upper,backing_Z_m,uncompressed_felt_thickness_m,prior_maximum_compression_strain
lower,backing_Z_m,uncompressed_felt_thickness_m,prior_maximum_compression_strain
felt,reference_stress_Pa,reference_strain,loading_exponent,unloading_exponent,residual_crush_fraction,densification_strain
creep,whole_washer_stiffness_N_per_m,whole_washer_viscosity_Ns_per_m
site,x_m,y_m,loaded_projected_area_m2
```

Both washers use the explicit supplied law but independent thickness and prior
conditioning. One to eight distinct positive-area sites are required; each
creates an upper and a lower patch. Areas are **loaded area projected normal to
Z**, not the shell's curved surface area or an annulus bounding-box area. Each
`creep` row describes one series branch for a whole washer. Its stiffness and
viscosity are apportioned by each site's fraction of the total area, so refining
site quadrature does not accidentally multiply either coefficient. There may
be zero to four creep rows. Zero explicitly means no rate-dependent branch;
conditioning/crush remains governed by WoolFelt.

The reference file preserves the legacy illustrative washer's total annular
area and whole-washer Kelvin coefficients but does **not** reproduce its fixed
boundary or constant manually assigned preload. At each site, skin height and
backing position determine preload. The initial carrier position is relative
to the ground spring's rest position and moves both backing planes equally.
Initial gaps, densification and the final composed owner's energy/frequency
limits are checked; there is no silent clamp, clipping or snapping over a hole.

## Scope and limits

This is a fixed-reference, fixed-axis translational support. Optional shell
rocking remains the existing bounded small-angle image. It does not simulate a
stand's complete beam network, lateral sliding, sleeve/hole collision, bolt or
wingnut mechanics, finite rotations, gravity, or a solved static equilibrium.
Preload can excite a startup transient, especially with initially relaxed
Kelvin branches. Supply identified data and settled initial conditions before
claiming a measured instrument response; an alloy name does not identify them.

Sites are explicit **area quadrature**, not proof of a continuous footprint or
a converged contact-pressure field. The importer checks the actual sampled
locations, not the unseen extent of each site's assigned area. Survey the
loaded footprint and refine its quadrature against response measurements; do
not assign area across mounting holes or outside the shell. The three-site
example is deliberately coarse, not an exact washer survey.

The `felt::mount::MountedFelt` library composition accepts arbitrary original
modal rows and signed carrier projections, independent of instrument type. This
CLI adapter currently admits only single-shell `splash`, `splash-wav` and
`splash-mic` commands, including supplied profile or explicit shell meshes.
Different profiles and site/backing geometry can reuse the same owners; a
suspended gong or paired hi-hat still needs its actual support kinematics, not
an automatically repurposed stand. No new time integrator or dependency is
introduced. This extension does not establish full-band modal convergence,
calibrated damping/contact, successful BEM-to-WAV rendering, or real-time speed.
