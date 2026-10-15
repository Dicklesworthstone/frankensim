# Source-bound electrical contacts

The existing `electrothermal` command accepts optional `electrical_interfaces`
in its scene. The same bindings work with the solid-only `--transient` mode;
there is no separate contact solver, command, material catalog or heating loop.

```sh
cargo run -p fs-cli --bin frankensim -- --json electrothermal scene.json
cargo run -p fs-cli --bin frankensim -- --json electrothermal scene.json --transient time.json
```

## Scene binding

Add this optional block to a `frankensim.electrothermal.v1` scene:

```json
"electrical_interfaces": [{
  "pack": "joint.fsintpk",
  "faces": [[1, 2]],
  "temperature_axes": ["cell_temperature_a", "cell_temperature_b"],
  "temperature_kind": "legacy",
  "state": [{"name": "normal_pressure", "value": "100000 Pa", "kind": "legacy"}],
  "heat_fraction_to_a": 0.25
}]
```

This is a scene fragment, not a complete runnable material dataset. The pack is
an existing canonical `NormalizedInterfacePack` (`FSINTPK`), containing an
ordered `InterfaceSystemCard` and its source claims. Paths are relative to the
scene. No measured joint properties are invented by this example.

Each face lists two **linear voxel cell indices**, in surface A/B order. They
must share a face. The interface card's named material states must match the
actual material assignments in that order. Reversed/duplicate/nonadjacent
faces refuse. One pack may bind many faces; the overall limit is 64 packs and
65,536 faces. Pack bytes share the existing aggregate material-input budget.

The existing `fs_material::electrical_interface` owner resolves the canonical
`electrical_contact_resistance` property. Its dimension-only schema is
`[L,M,t,T,I,N] = [4,1,-3,0,-2,0]`, **ohm m²**: area-specific resistance, not bulk
resistivity in ohm m or a lumped resistance in ohms. The CLI does not define a
second property resolver or remove source quantity semantics. Values must be
finite and non-negative; an explicit zero is a perfect contact. A missing,
ambiguous, wrong-unit or unsupported claim is an error, not a fallback zero.

The two named temperature coordinates receive the adjacent **cell-state
samples**, in kelvin. They are not inferred microscopic contact-surface
temperatures. The selected source must be appropriate for this convention
and spatial resolution. `temperature_kind` may be `legacy`, `dimensional`, or
`absolute-temperature`; it must match the source's coordinate descriptors.
Temperature differences are not absolute temperatures. Fixed coordinates such
as pressure use the same explicit quantity parser as bulk material queries.
Both temperature axes must be distinct and absent from the fixed `state`.

`heat_fraction_to_a` is an explicit caller-owned, zero-thickness heat-deposition
model: the fraction of **interface-generated** heat assigned to cell A. The
remainder goes to B. It is not an inferred equal split, a measured partition,
a thermal-effusivity calculation, or a function of current direction. A more
detailed interface thermal model would need its own resolved state and storage.

## One conservative electrical and thermal path

For face area A, bulk half-widths w_A/2 and w_B/2, and area resistance R_c,

```text
s_A = rho_e,A * w_A/2
s_B = rho_e,B * w_B/2
G_face = A / (s_A + R_c + s_B)
Q_face = G_face * (phi_A - phi_B)^2
Q_contact = Q_face * R_c/(s_A + R_c + s_B)
```

An ideal electrode has zero bulk half resistance. The shared spatial current
balance uses this G; the shared Joule calculation uses the same G. Bulk heat
stays in its own half-cell. Contact heat goes to the two thermal cells according
to the declared fraction, including electrode cells. Reversing current leaves
this partition unchanged. Direct contact between two distinct ideal electrodes
remains refused; a resistive face does not make an insulating cell conduct.

The returned `joule_w` contains **bulk plus contact heat exactly once**. The
existing thermal coupling adds it to authored sources. Do not add another
I²R_contact heat source. Thermal conductivity, authored thermal contact
resistances, boundary conditions and storage continue through their existing
owners. The scene's separate `contacts` block is still its independently
declared thermal contact law; electrical resistance does not infer it.

At every proposed thermal endpoint the command re-queries all electrical
interface claims at both current cell temperatures and the complete fixed
state. Source-domain refusal aborts publication. The existing coupling loop
checks temperature and per-cell Joule-source consistency. Transient mode uses
the same path, preserving accepted temperatures and fixed cell masses across
scheduled drive changes; contact heat enters stored enthalpy without a second
deposition pass.

Both JSON result schemas include `electrical_interfaces`: bound ordered cells,
final area resistance, declared heat fraction, canonical pack/card identities,
and the final resolved interface-state identity. Bulk material identities and
the existing power/energy defects remain alongside them.

The existing prescribed-fluid steady command is preserved. Its source-resolved
reference transport and flow field pass unchanged to the heat owner, while
per-cell electrical assignments and contact faces use the same potential
solver. This does not add a momentum solver or variable-density advection.
Fluid scenes remain refused by the solid-only transient storage model.

## Focused verification

Five electrical-kernel regressions cover voltage/current drive signs, electrode
heating, zero resistance, reversed surface roles, conducting-fluid parity and
invalid faces. Five file-driven command regressions exercise canonical
interface decoding, material-role matching, data-only swaps, nonlinear joint
feedback, source support, replay and transient heat conservation.

Independent NumPy calculations checked 1,000 resistor networks: maximum power
discrepancy was 4.43e-14 W. For the four-cell unit-area/unit-conductivity command
fixture with two 2-ohm bulk cells, a 2-ohm contact and 3 V, the reference current
is 0.5 A and deposited powers are 0.625 W and 0.875 W. The steady temperatures
are 300.359375, 301.078125, 301.171875 and 300.390625 K. With
R_c=2*(1+0.02*(T_A-300)), an independent nonlinear reference gives total power
1.489377453809321 W. These are synthetic equation references, not real metal
measurements or outputs from an executed Rust binary.

The insulated transient reference receives 0.6 J over four 0.1-s powered steps
and retains it during four unpowered redistribution steps; independent maximum
per-step closure was below 7e-17 J. Rust compilation/tests were **not executed**
in the implementation environment: Cargo/rustc/rustfmt are unavailable.
Run the focused tests under the repository's normal execution policy:

```sh
cargo test -p fs-lbm --lib conjugate::electric::contact::tests
cargo test -p fs-cli --bin frankensim electrothermal_command::interface_tests
```

This is fixed-geometry ohmic conduction with explicitly bound resistive faces.
It does not determine true contact area, pressure evolution, oxide rupture,
microscopic surface temperature, tunnelling, dielectric breakdown, melting or
AC behavior. Source nominal values and selected heat partitions are not
experimental qualification or propagated uncertainty.
