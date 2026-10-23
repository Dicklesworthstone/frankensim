# Profile-defined chambers in the existing percussion player

`--chamber-profile input.fscp` selects the actual piecewise-linear inner and
outer wall profile. It implies `--cavity-modes`. Both drumhead pencils, physical
stick impacts, felt mallets, shaft flexure, material history and the selected
time solver are unchanged. It works on drum/snare CSV, WAV and microphone
commands, including the existing reciprocal radiation-feedback route.

```sh
cargo run --release -p fs-couple --example percussion -- --drum-spec crates/fs-couple/examples/percussion/estimated_drum.fsd \
  drum 128 --analytic-newton --strike-speed-m-s 4 \
  --chamber-profile crates/fs-couple/examples/percussion/estimated-waisted-chamber.fscp
```

The shipped file is an illustrative waist, **not measured hardware**. With no
profile, the original cylindrical construction remains unchanged. Supply a
matching `--drum-spec` for different end radii, depth, materials or head tensions.

The file has a `frankensim-chamber-profile-v1` header, then mandatory records:

- `mesh,radial_intervals,maximum_azimuthal_order`
- `band_hz,maximum_pressure_frequency`
- `budget,maximum_nodes,maximum_matrix_entries,maximum_modes,residual_tolerance`
- Ordered `station,z_m,inner_radius_m,outer_radius_m` records, starting at z=0.

Z increases downward from the batter head. Both end inner radii equal the clear
head radius; both end outer radii equal the bearing-edge outer radius. The final
z equals the drum depth. No snapping or inferred thickness repairs mismatches.
At most eight pressure modes are admitted by this player; incomplete Fourier
pairs or an overfull frequency window refuse, rather than discard resonances.

The shared meridional P1/Fourier owner derives pressure modes from the inner
profile. Its exact constant member uses the profile's frustum volume, replacing
rather than supplementing the old bulk spring. The same head quadrature supplies
signed reciprocal pressure work. Extra sticks, wires, pads and shafts receive no
direct pressure coupling. Nonuniform momentum drag remains explicit.

The outer profile and the actual head meshes form one closed exterior BEM
boundary. Wall stations are retained, including geometric corners, and rigid
wall source weights are zero. The two CSV interior probes move to the axis at
one-quarter and three-quarters of the actual chamber depth. They are pressure
observations, not additional forces or exterior audio channels.

Sidewall vents and elastic barrels currently use cylindrical geometric adapters;
combining either with this profile explicitly refuses. Wall flexibility, measured
manufacturing stress, arbitrary non-axisymmetric cavities and realistic porous
losses are not inferred. Algebraic residual bounds do not certify mesh or angular
convergence; refine stations, radial/angular resolution, retained modes and BEM
panels independently. Existing acoustic fitting/clock/work limits remain in force.
No native callback deadline, full-band accuracy or calibrated specimen is claimed.
