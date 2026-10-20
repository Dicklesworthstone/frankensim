# Source-backed thermal joints in electrothermal scenes

Thermal resistance can be assigned to individual voxel faces, including a joint
between two pieces of the same material. No pseudo-material IDs, averaged joint
coefficient, new heat solver or electrical-to-thermal conversion is introduced.
The existing steady command, prescribed-fluid steady command and all-solid
transient command share this path.

## Scene inputs

Add an independent `thermal_interfaces` list to the existing
`frankensim.electrothermal.v1` scene:

```json
{
  "thermal_interfaces": [{
    "pack": "joint.fsintpk",
    "faces": [[0, 1]],
    "temperature_axes": ["Ta", "Tb"],
    "temperature_kind": "legacy",
    "state": [{
      "name": "normal_pressure",
      "value": "100000 Pa",
      "kind": "legacy"
    }]
  }]
}
```

This is a scene fragment, not a complete runnable scene. `faces` contains ordered
linear cell indices; x varies fastest. The two cells must share a face. Binding a
pack assigns its ordered surface/finish descriptors there; its material states
must exactly match the scene's bulk assignments. Temperature and fixed-coordinate
kinds must match the actual source schema. The example's untyped `legacy` axes
are not a recommendation to strip semantic source metadata. The source must
support both supplied cell-temperature coordinates, not an inferred average or
an assumed microscopic contact temperature.

The pack must contain exactly one claim named
`area-specific-thermal-contact-resistance`, in m² K/W. Resolution reuses
`fs_conduction::interface::InterfaceResistance::from_card_pinned`, including
its strictly positive, finite resistance and source-support rules. The selected
claim is pinned for the entire run: ending its support is a refusal, not a switch
to another curve. A missing electrical property is irrelevant to a thermal-only
joint. A bulk material declared `conducting: false` also needs no resistivity.

Supply an explicit positive, less-than-one convergence tolerance in `numerics`:

```json
{"thermal_contact_relative_tolerance": 1e-10}
```

The existing thermal iteration budget and temperature/conductivity tolerances
remain in force. Old scenes without thermal interface bindings need no additional
input. The existing `contacts` material-pair shorthand remains available, but a
face cannot also have a source-backed thermal law: competing owners are refused,
not added together or resolved by last-write-wins.

Electrical and thermal interface bindings are independent. The same FSINTPK
pack and ordered face may appear in both lists when the card actually supplies
both responses. Electrical area resistance is in ohm m²; it is NOT a thermal
resistance. Thermal joints have no `heat_fraction_to_a`: they transfer heat,
whereas the electrical joint's heat partition deposits its generated Joule heat.

## One heat-transport calculation

For cell widths w, molecular/effective conductivities k, face area A and the
resolved contact resistance R'', the existing face flux uses:

```text
G = A / (w_a/(2*k_a) + R'' + w_b/(2*k_b))
Q_a_to_b = G * (T_a - T_b)
```

The same coefficient participates in steady and implicit transient assembly,
energy accounting and solid-to-fluid heat observation. Internal transfer enters
neighboring balances with opposite signs; it is neither a new heat source nor a
reservoir. Conductivity, electrical heating and stored enthalpy remain owned by
their existing implementations.

Every thermal iteration resolves the contact at its current cell temperatures,
then re-evaluates it at the proposed solution. Acceptance requires both the usual
temperature/conductivity/storage checks and the new resistance-consistency check.
The electrical feedback loop sees the resulting temperatures and re-evaluates
its own material and interface laws. Neither contact heating nor storage is
counted twice.

Both JSON result formats include `thermal_interfaces`, with ordered cells,
final source-resolved resistance, resistance actually used in the last linear
solve, pack/card identities and the final query-receipt identity. The top-level
`thermal_contact_relative_defect` reports their maximum relative discrepancy.
They agree within the declared tolerance; finite iteration is not exact equality.

The low-level `ThermalSetup::face_contacts` / `FaceContactResistance` interface
also accepts explicit zero resistance. Source-backed CLI input keeps the existing
thermal source owner's stricter positive-resistance contract. Face contacts
refuse duplicate faces, overlapping material-pair rules, compact-component cells,
unrepresentable conductance and any through-flow. A solid/fluid wall is supported;
a zero-thickness stationary contact is not a porous-flow or advective membrane law.

## Focused checks and limits

Six low-level tests exercise graded same-material joints, reversed ordering,
solid/fluid heat reporting, transient energy exchange, material-pair equivalence,
invalid geometry, conflicting laws and cancellation. Five command tests exercise
actual file reads, canonical packs, source queries and numerical owners: thermal-
only swaps, nonlinear resistance, shared electrical/thermal faces, transient
redistribution, source support and input refusals.

Independent NumPy reference calculations were executed. In the synthetic
three-cell, unit-conductivity fixture with a 4 W/m² input and 300 K outlet,
changing only the joint resistance from 1 to 3 m² K/W changes the first cell from
314 K to 322 K; the other cells remain 306 K and 302 K. For
R'' = 1 + 0.02*(T_a - 300), the reference first cell is 315.2173913043 K.
An insulated eight-step heating/redistribution reference retains 0.8 J for either
joint data set (maximum reference per-step closure below 6e-17 J).
These are synthetic equation checks, not Rust execution or material measurements.

Native Rust compilation, formatting and tests were not executable in the
implementation environment because Cargo/rustc/rustfmt were unavailable.
Run the focused tests under the repository's normal execution policy:

```sh
cargo test -p fs-lbm --test thermal_face_contacts
cargo test -p fs-cli --bin frankensim electrothermal_command::thermal_interface_tests
```

The model remains fixed-geometry with nominal source values. This does not solve
contact pressure, microscopic contact area, interfacial heat capacity, thermal
expansion, debonding, phase change, uncertainty propagation or model calibration.
Prescribed-fluid steady transport keeps its existing frozen reference assumptions;
the command's fluid-transient restriction is unchanged.
