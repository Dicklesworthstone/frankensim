# A driven gas cylinder and finite receiver

This example connects the existing Otto slider-crank geometry to the public
`fs_couple::prescribed_gas` control-volume solver. A prescribed drive compresses
a source-loaded gas while a restriction carries mass and upstream total enthalpy
into a finite receiver. The restriction then closes and the drive returns through
an expansion stroke. Both gas inventories retain their evolving temperature,
pressure, composition and source material models.

The CSV contains the actual pressure–volume trajectory, both gas temperatures
and masses, and signed work and torque at the drive. Positive drive work leaves
the gas; negative work must be supplied by the drive. Average torque is the
accepted interval's net pressure work divided by its crank-angle change, so the
displayed torque integrates to exactly the reported mechanical work.

The same upstream enthalpy removed from the cylinder enters the receiver once.
It already includes flow work and the jet energy thermalized in the receiver.
There is no second jet-heat term. Every step checks local and global energy;
the example also independently reconstructs cumulative gas storage plus gross
boundary work. CSV is published after the complete trajectory succeeds.

## Run

From the repository root, compile the existing nitrogen source packs:

```bash
run_dir="$(mktemp -d)"
demo=crates/fs-couple/examples/driven_gas_cycle
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nitrogen/manifest.tsv \
  --out "$run_dir/nitrogen.species"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nitrogen-nasa9-mcbride2002/manifest.tsv \
  --out "$run_dir/nitrogen.nasa9"
cargo run -p xtask -- matdb-pack \
  --manifest data/matdb/seed-v1/nitrogen-transport-svehla1995/manifest.tsv \
  --out "$run_dir/nitrogen.transport"
cargo run --release -p fs-couple --example driven_gas_cycle -- \
  "$run_dir/nitrogen.species" "$run_dir/nitrogen.nasa9" \
  "$run_dir/nitrogen.transport" "$demo/cylinder.gas" "$demo/receiver.gas" 256 \
  > "$run_dir/cycle.csv"
```

The optional final argument is the number of intervals **per stroke**, from 8
to 4,096. Repeat with 512 intervals to assess time refinement. Source-domain,
withdrawal or pressure-crossing failures refuse the trajectory rather than
clipping its mass, heat or pressure.

## Explicit physical model

The supplied engineering geometry is a 50 mm bore, 1 mm crank radius, 4 mm
connecting rod and 80 cm³ clearance, driven at 120 rpm. The receiver holds
2 litres. The restriction area is `5e-8 m²`, with an explicitly selected
discharge coefficient of 0.8. Exterior pressure is held at 80 kPa. These are
demonstration inputs, not dimensions or measurements of Otto's historical engine.

The short stroke deliberately stays inside the existing transport source cards'
authored **75–101.325 kPa** application rectangle. The library supports larger
volume changes when the caller supplies material sources valid for that path;
the constant-capacity and nonlinear-caloric tests exercise substantial
compression and an air-standard Otto cycle. This example does not widen a
source's physical validity to make a high-pressure result appear admissible.

The model is spatially mixed, ideal-gas and adiabatic, with frozen composition.
The valve is explicitly open for the compression stroke and closed for the
return; it is not an independently solved check-valve actuator. This is one
compression/discharge and return journey, not a periodic intake/combustion engine.
The finite receiver changes state as it receives gas. No infinite pressure
reservoir, invented shaft inertia or external heat source is substituted.

Flow is evaluated from one accepted snapshot and then the original nonlinear
midpoint material law advances pressure work. The combined open-system method
is first order. Its full energy closure does not establish exact finite-step
isentropy, experimental compressor efficiency, or time-error certification.
