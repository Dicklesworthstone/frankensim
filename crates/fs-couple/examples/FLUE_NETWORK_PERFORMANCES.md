# Branched flue resonators and regional atmospheres from data

The `music_render wind` command accepts `flue-network.performance` through the
same `frankensim-flue-performance-v1` source format as the single-pipe example:

```sh
cargo run --release -p fs-couple --bin music_render -- wind \
  crates/fs-couple/examples/flue-network.performance flue-network.wav
```

This is a synthetic geometry and performance, not a characterized flute. One
jet drives the complete connected acoustic graph. Branches are not separate
voices. Every section derives travel time, characteristic impedance and any
selected viscothermal loss from its own coherent material gas state. Every
selected opening contributes its actual terminal volume flow to the existing
multirate Rayleigh receiver. The source never substitutes internal pressure for
exterior microphone pressure.

## Geometry and material gas

Replace the original `tube LENGTH RADIUS LENGTH_ERROR REFLECTION` row with:

```text
network NODE_COUNT SECTION_COUNT
duct_node KIND ...
... one duct_node per node ...
duct_section NODE_A NODE_B LENGTH_M RADIUS_M LENGTH_ERROR_M [gas GAS_RECORD]
... one duct_section per section ...
```

Node indices and section indices are zero-based in source order. The shared
network owner admits one degree-one inlet, connected sections, valid junction
and load degrees, and passive loads. The file limits are 64 nodes and 128
sections. All dimensions and numerical error allowances are explicit.

| Node declaration | Physical meaning |
| --- | --- |
| `duct_node inlet` | The unique port driven by the jet |
| `duct_node junction` | Common-pressure junction with no added storage |
| `duct_node reflection R` | Explicit memoryless pressure reflection in [-1,1] |
| `duct_node impedance R L C` | Degree-one series resistance, inertance and compliance |
| `duct_node series R L C` | Two-section inline pressure-drop load |
| `duct_node shunt R L C` | Load from an interior common-pressure junction to its reference |

For R-L-C records, units are Pa s/m³, Pa s²/m³ and m³/Pa. Zero `C` explicitly
omits compliance. Positive compliance and nonnegative resistance/inertance are
checked by the existing passive impedance owner; the file parser contains no
replacement impedance or propagation equations. Supplied coefficients are not
promoted to measured radiation, tissue or structural material properties.

A section without `gas` uses the top-level gas record. A supplied record uses
the shared material gas parser, including explicit species mixtures and loaded
thermochemical sources. All regions must have the same static pressure. The
jet and a selected windway use the **actual inlet section's** density,
viscosity, sound speed and other gas properties, even when that section is not
first in the file. `exterior same` uses the top-level default gas; an explicit
exterior record describes the receiver environment independently.

The example includes a warmer N₂/O₂/CO₂ inlet, an ambient main tube and a
separate branch temperature. These are prescribed frozen states. They do not
claim to solve breath advection, thermal exchange, gas mixing or condensation.

## Distributed loss and exterior openings

A network uses `viscothermal none`, or explicitly addresses its original
physical sections:

```text
viscothermal sections 2
loss_section 0 CELLS MINIMUM_FREQUENCY_HZ MAXIMUM_FREQUENCY_HZ
loss_section 2 CELLS MINIMUM_FREQUENCY_HZ MAXIMUM_FREQUENCY_HZ
```

The original section indices survive the shared loss subdivision. No gas is
averaged across branches. Existing gas, wide-tube, fractional phase, memory and
frequency-domain admission checks remain in their physics owners. A plain
single tube retains its existing short viscothermal record.

After `receiver`, an explicit network must select its openings:

```text
outlets PRIMARY_NODE ADDITIONAL_COUNT
outlet OTHER_NODE CENTER_X_M CENTER_Y_M
... one outlet per additional opening ...
```

The primary opening defines the origin of the common baffle plane. Additional
centers and receiver coordinates use that same frame. The shared receiver
checks terminal topology, opening radii, overlaps, common pressure and its
quadrature limits. The example's two outlets contribute to one pressure
history, preserving propagation timing and interference. A single-tube source
may omit this block and retain its original sole terminal receiver.

Prescribed terminal loads retain their one-way infinite-baffle receiver. The
explicit MIMO selection below instead supplies shared exterior feedback. These
remain plane-wave cylindrical sections and ideal junctions; the source file
does not establish measured flute fidelity, a resolved jet or a validated
three-dimensional radiation model.

## Explicit mutual radiation from a supplied artifact

`flue-mimo.performance` and its adjacent `flue-mutual.mimo-radiation` exercise
the existing common terminal operator through the native loader:

```sh
cargo run --release -p fs-couple --bin music_render -- wind \
  crates/fs-couple/examples/flue-mimo.performance flue-mimo.wav
```

The artifact is an explicitly **synthetic passive two-port example**. Its
coefficients demonstrate the complete source-to-feedback-to-microphone route;
they are not a measured flute or a validated exterior-radiation bake.

After optional wall and moving-boundary records, select:

```text
radiation mimo flue-mutual.mimo-radiation 32 1048576 257
```

The fields are the relative artifact path, maximum rational states, maximum
runtime memory bytes and sampled-passivity grid points. The common MIMO decoder
also limits artifact input to 8 MiB. The original runtime enforces the declared
state/memory limits and checks both continuous and actual sampled responses.
The sampled result remains `GridOnly`; it is not a discrete storage certificate.

Artifact ports follow the exact `outlets`/`outlet` order. Their circle centers
and radii must match those openings. The artifact separately declares its
exterior gas using the same gas-record grammar as the performance, and its
frequency band must cover the entire mechanical Nyquist band, even when the
receiver emits a narrower audio band. Existing reactive terminal loads cannot
be replaced. No matrix is inferred from the instrument name, fitted at load
time or rescaled to conceal a stale gas declaration.

Each accepted terminal flow drives the one shared rational state, whose full
matrix pressure returns to the same duct and jet solve. The same openings and
exterior gas drive the multirate microphone. Changing off-diagonal coefficients
therefore changes the jet's returning wave as well as emitted pressure. Material
walls, regional section gases, windway dynamics, mouth contraction, moving-key
mechanics and the gesture clock remain attached to their original owners.

An atmosphere change that affects the baked operator requires a matching
artifact. Changing only the performance gas while retaining an incompatible
artifact refuses before rendering. The library entry point
`FlueFilePerformance::from_bytes_with_all_sources` accepts the existing material
loader plus a radiation-byte loader; the CLI resolves both relative to the
performance directory. `mimo_radiation_record()` retains the exact artifact
hash, declared identities, gas record, geometry, band and limits. The existing
constructors and `prescribed`/`compact-baffled` selections retain their behavior.

## Pressure and jet direction in one performance

`flue-embouchure.performance` adds a prescribed jet-angle phrase and the
existing separated-flow mouth loss to the branched, regional-gas example:

```sh
cargo run --release -p fs-couple --bin music_render -- wind \
  crates/fs-couple/examples/flue-embouchure.performance flue-embouchure.wav
```

The source is a manufactured fixed-window demonstration. Its geometry, player
controls and loss coefficient are authored inputs; none was inferred by
matching an instrument recording.

The canonical gesture suffix requires one `blowing-pressure` track and accepts
one additional `jet-angle` track, whose values are `angle-rad`. Track IDs are
arbitrary and distinct. Both tracks use the original gesture sampler, including
interrupted ramps. A tick applies at
`ceil(tick * mechanics_rate / control_rate)`, before that mechanical sample.
Changing callback size does not alter the physical control clock or discard
partial microphone frames. The `controls` limit bounds the number of tracks
times control ticks before held-value coalescing.

The direction enters the existing geometric cut,
`labium_offset + cut_up * tan(jet_angle)`. It changes the jet volume entering
the mouth, and the same source-pressure derivative and implicit mouth solve
propagate that change through every duct and selected outlet. Angles must be
finite and strictly inside `(-pi/2, pi/2)`, with a representable transverse
offset. The complete phrase is checked before any source sample.

This is the existing instantaneous tilted-cut approximation at fixed mouth
area, cut-up, exit width and exit height. It does not resolve transient jet
turning or lip movement, change lip coverage or remesh convection. The
windway velocity field, jet feedback history, wave propagation, material-wall
memory and radiation history persist. A pressure release still leaves stored
fluid and acoustic motion intact. Dynamic cut-up or mouth-area changes need a
moving-domain and boundary-work model and are not exposed as simple controls.

## Mouth separation from the actual geometry and gas

After the optional `windway` row and before `jet_card`, select:

```text
mouth_loss contraction 0.85 0.3 authored-fixed-window-demonstrator
```

The fields are contracted jet area / physical opening area, maximum contracted
Mach number, and one explicit source token. The existing `ContractionSpec`
owner validates them and binds the loss to the physical mouth window
`flue_width * cut_up` and the **actual inlet section's complete gas state**.
There is no independently entered density or sound speed. The example's source
label explicitly identifies a manufactured coefficient, not measured edge data.

The mouth solve includes the signed pressure drop
`rho * Q * abs(Q) / (2 * (Cc * area)^2)` at the accepted flow, together with
mouth inertance, the supplied linear mouth resistance and the same acoustic
port. Its power is recorded as separation loss; it is distinct from radiation
and windway viscosity. The solved contracted Mach limit remains active during
rendering. `Cc = 1` still supplies this loss law. Omit the row or use
`mouth_loss none` for the original linear mouth.

The CLI sidecar includes the supplied coefficient/source, bound effective area
and actual gas density/sound speed. This remains a quasi-steady separation
closure with no vortex memory, resolved shedding or measured flute accuracy.
It does not qualify the full Boehm assembly or the plan's register experiment.
