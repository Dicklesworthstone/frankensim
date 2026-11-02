# Resolved flow, finite breath and a complete acoustic field

`resolved-thermoviscous-breath.performance` combines the existing manufactured
NSF grid and six-tetrahedron acoustic field from `resolved-thermoviscous.performance`
with the separate supply/escape openings and finite-source control law from
`resolved-dg-actuated-mixture.performance`. These are numerical fixtures, not
measurements of a flute, player or windway.

```bash
cargo run --release -p fs-couple --bin music_render -- resolved \
  crates/fs-couple/examples/resolved-thermoviscous-breath.performance \
  /tmp/resolved-thermoviscous-breath.wav --block 37 --decimate --self-contained
```

The finite supply retains its actual mass, energy, volume and gas composition.
Its declared volume motion and heat change pressure and temperature using the
original material law. The same accepted mass and total-energy flux from the
resolved flow updates that inventory. The nonlinear flow couples to the full
three-dimensional thermoviscous field and then the original finite DG duct.
The added field does not reset or reconstruct the supply from pressure scalars.

The microphone includes the duct outlets and the independently declared escape
tiles. Each surface retains its original signed flow, physical travel distance,
gas attenuation and phase. Cold observer histories move with those unchanged
surfaces and their output clock when the complete field is inserted. Adding a
microphone path does not change the source trajectory.

The embedded `gas` and `dg_exterior` records accept the same named mixtures and
caller-supplied species, NASA9 and transport packs as ordinary `resolved` input.
The existing `--gas-component`, `--gas-water`, `--gas-phase` and source attachment
options also apply to this schema. Original `dg_wall` material-card records use
the same admitted material loader. For example, after declaring a supplied
thermochemical species in the embedded gas records and its `energy_solve`:

```text
--gas-component SPECIES.pack NASA9.pack TRANSPORT.pack --gas-water none --gas-phase frozen
```

In `--self-contained` mode every referenced pack must be included through the
existing attachment mechanism. Missing or incompatible material data refuses;
it never selects dry air or a replacement wall. Changing composition preserves
the physical mesh and output clock. The flow owner derives its own physical
time step from that gas, and the existing conservative observation clock handles
the difference without retiming the physics.

The sidecar now includes the original embedded source information, including
its gas and material records, finite inventory and observed escape placement,
alongside the complete-field geometry. It is valid JSON like the other native
routes.

## Limits of this composition

The inserted field is a linear thermoviscous acoustic region at the original
reference gas state. It does not transport mean breath mass or enthalpy into
the bore. The finite supply uses the existing held-state, first-order partition;
physical and observation clocks still require convergence checks.

The escape in this fixture retains its prescribed ambient boundary and a
one-way physical microphone. This native handoff still refuses an already
attached upstream acoustic load or reciprocal mouth-radiation load, because
their direct acoustic port basis must be explicitly rebound to the complete
field. The typed complete-field runtime supports such constructions separately.
No reduced jet card, empirical displacement gain, synthesized noise or measured
instrument qualification is introduced by this example.
