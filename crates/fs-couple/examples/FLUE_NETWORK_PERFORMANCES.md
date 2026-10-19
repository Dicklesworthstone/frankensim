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

This remains a plane-wave graph of cylindrical sections, ideal junctions and
supplied passive loads with a one-way infinite-baffle receiver. It does not
provide moving finger pads, branch end corrections, matched mutual radiation,
a resolved jet or measured flute fidelity. These limitations are unchanged by
making the existing physical construction accessible from data.
