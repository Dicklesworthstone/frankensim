# Shared-state material gas networks

`fs_material::thermochemical_gas::energy::parcel::transfer::snapshot`
provides multi-inlet/outlet conservative transfer from one accepted snapshot.
`transfer::nozzle::network::snapshot` supplies the pressure-driven consumer.
The previous `nozzle::network::step_candidate` remains a separate, unchanged
Lie-split method; no caller silently changes time-discretization semantics.

## All connections see the same physical state

A nozzle's pressure, composition, stagnation enthalpy, throat sound speed and
choking condition are resolved from its ORIGINAL upstream parcel. After all
rates are known, each volume combines every incoming/outgoing extensive flux
and reconstructs its material state once. Renaming link IDs or rearranging the
input vector cannot select a different physical transfer sequence. Identical
receivers connected through identical restrictions receive identical updates.

For a directed integrated mass `dm`, the unchanged donor composition `Y` and
specific enthalpy `h` contribute:

- `delta(m) = +/- dm`
- `delta(m_species) = +/- dm * Y_species`
- `delta(U) = +/- dm * h`

An intermediate volume's outgoing gas is its old gas, not a mixture containing
material that just arrived during this interval. That new mixture can leave on
the next interval. This is a first-order explicit finite-volume discretization,
not an implicit simultaneous pressure solve or a finite propagation-time model.

The generic transfer owner accepts prescribed fluxes, including simultaneous
counterflow. The nozzle adapter instead selects direction from the original
pressures and reuses the existing source-backed subsonic/choked nozzle law.
No constant-air gamma, duplicate NASA coefficients, new species registry,
caloric solver or transport rule is introduced. The existing source join
retains matching species laws, transport cards, explicit water identity and
pack bindings; incompatible shared-species sources still refuse.

## API and acceptance

The generic `snapshot::transfer_candidate` accepts an immutable
`BTreeMap<u64, ClosedGasParcel>`, integrated `GasSnapshotTransfer` rows,
exactly one `GasEnergySolve` per parcel and `GasSnapshotSolve` controls.
Parallel rows are combined in canonical numeric order. Zero-flow and isolated
volumes remain bitwise unchanged. Hard work limits are 256 volumes and 4096
submitted transfers; the species bound remains the existing material owner's.

The nozzle `snapshot::step_candidate` takes `SnapshotNozzleLink` rows with
explicit reporting id, two parcel ids, `GasNozzle` geometry and original throat
solve controls, plus physical seconds and the same endpoint/control maps. Its
proposal retains every source-bearing throat result and the conservative
multiport account. `into_parts()` moves all endpoints together. Do not deposit
the reported enthalpy again as heat or resolved jet kinetic energy.

Gross withdrawals are limited against ORIGINAL inventory before incoming gas
can conceal depletion. Local total mass, every local species, local energy,
global mass, every global species and global energy are checked from actual
endpoints with absolute kg/J tolerances. There is no automatic clipping.

After the joint update, initially active nozzle paths must retain their
pressure ordering. Two individually acceptable incoming pair steps can jointly
overshoot equilibrium; the whole graph proposal then refuses. Reduce/refine
the physical interval. No event time or unconditional stability is inferred.
Source failures, unbracketed caloric states, exhausted work and cancellation
return no partial graph, including at the final publication checkpoint.

## Run the four-volume material example

From the repository root:

```sh
cargo run --release -p fs-material --example material_gas_network -- \
  crates/fs-material/examples/material_gas_spring/mixture-a.tsv \
  crates/fs-material/examples/material_gas_spring/mixture-b.tsv \
  200 0.0002 0.000001 > gas-network.csv

cargo test -p fs-material --lib thermochemical_gas::energy::parcel::transfer::
cargo test -p fs-material --example material_gas_network
```

The fixed example topology is `0 -> 2`, `1 -> 2`, `2 -> 3`, with direction free
to change on later intervals. Each volume is 0.001 m3. Initial states are
(500 K, 300 kPa), (450 K, 220 kPa), (400 K, 140 kPa), (350 K, 100 kPa).
Node 0 uses the first file; the others use the second. Reverse only the two
files to change material composition without rebuilding or changing physics.
Both files share one explicit synthetic source set through the original loader.
They are not experimentally validated real gases or a production input schema.

CSV contains one row per node/time with mass, temperature, pressure, sound
speed, impedance, viscosity, conductivity, incoming/outgoing mass, global energy
drift, species residual and the number of choked connections. The original
pair example remains unchanged.

## Verification and physical limits

Seven core, six nozzle-graph and one executable Rust regressions are authored.
Native execution was unavailable during authoring: Cargo and rustc are absent.
Independent Python equations checked 256 random networks and 256 reordered
copies (maximum global energy discrepancy 9.20e-11 J and species discrepancy
4.34e-18 kg). A separate nozzle/trajectory reference executed 624 steps and
compared pressure against DOP853; first-order refinement ratio was 2.0025 and
maximum energy drift 1.71e-11 J. These are equation checks, not Rust passes.

Fixed-volume, well-mixed, frozen-composition ideal-gas transport only. No reaction,
resolved jet momentum, vacuum filling, wall motion, thermal boundary solve,
spatial acoustic wave or live renderer reconfiguration is implied. Existing
heat/work/acoustic consumers may use the returned common states; composing
other operators over the same duration requires explicit time/error accounting.
Conservation does not prove entropy production, uniqueness or timestep accuracy.

Physical balance reference: Cantera, Control Volume Reactor and Reactor
Interactions, https://cantera.org/stable/reference/reactors/interactions.html.
