# One material table: electrical heating, thermal history and mechanism response

`fs-mbd-material::thermal::heating::HeatedArticulated` connects complete
material-defined mechanisms to one `fs-thermal` enthalpy timeline. The optional
`electrothermal-coupling` feature adds `advance_joule` to it and to `HeatedBody`.
This calls the existing LBM spatial current/enthalpy driver, not a new circuit,
thermal or mechanical solver.

## Run a data-only comparison

From the repository root, after the pinned dependency bootstrap:

```sh
cargo run -p fs-mbd-material --features electrothermal-coupling \
  --example electrically_heated_rotor -- \
  crates/fs-mbd-material/examples/electrically_heated_rotor/cp10.txt 20 2 4 0.25

cargo run -p fs-mbd-material --features electrothermal-coupling \
  --example electrically_heated_rotor -- \
  crates/fs-mbd-material/examples/electrically_heated_rotor/cp20.txt 20 2 4 0.25
```

Arguments are material file, applied voltage [V], voltage-shutoff time [s],
final time [s], and maximum thermal/output step [s]. The example is deterministic
and has no random forcing. Use successively smaller steps to check temporal
accuracy. There is no assertion that the requested maximum step resolves an
arbitrary nonlinear source. Both times and the step must be finite and positive;
shutoff precedes the end. The example bounds requested outputs and actual thermal
substeps instead of allowing an unbounded trajectory.

The two supplied tables differ **only in Cp**. Their coefficients are synthetic,
not measurements of a metal. Geometry, reference density, expansion, electrical
resistivity, conductivity, electrode drive and mechanical probe are unchanged.
No expected-output curve or precomputed acceleration is read by the program.

The table has a version header followed by six whitespace-separated SI columns:

```text
frankensim-electrothermal-material-v1
# T[K] density[kg/m3] Cp[J/(kg K)] alpha[1/K] resistivity[ohm m] k[W/(m K)]
300 2 10 0.002 5 2
500 2 10 0.002 15 6
```

Two through 256 increasing temperature rows are accepted, with a 64 KiB input
limit. Numbers must be finite; all columns except alpha must be positive. The
example compiles these data into one immutable `MaterialCard` using ordinary
source claims. `fs-matdb` owns interpolation and state-domain refusal;
`fs-material` owns the Cp/enthalpy chart and integrated alpha(T). No extrapolation
or material-name inference is added. The first row defines the reference
absolute temperature and zero specific-enthalpy datum. Pressure is explicitly
fixed at 100 kPa in this example's source context.

## What the output represents

Three one-cubic-metre reference cells form an idealized two-link assembly.
The outer cells are ideal electrical electrodes and belong to the fixed mount;
the central conducting cell belongs to the rotor. Ordinary thermal links join
adjacent cells. Their conductances use the actual adjacent conductivity samples,
one-half-metre series lengths and a perfect thermal interface. No mechanical
joint silently implies thermal conductance.

Each actual thermal substep resolves current and bulk/contact Joule power,
advances the existing heat timeline, then configures all link regions through
the original material expansion owner. Each region keeps its exact reference
mass. Density at the initial temperature defines that mass; heating does not
replace it with a new rho(T) times an old volume. Current geometric density is
mass divided by expanded volume.

CSV rows contain the accepted time, voltage applied over that interval, three
cell temperatures, endpoint current under that drive, stored thermal energy,
cumulative electrode work, rotor volume and mass, and the **existing articulated
solver's acceleration under a one-newton-metre test torque** at q=qdot=0.
This is a probe of the configured mechanism, not an integrated rotor trajectory.
At the shutoff endpoint the row still describes the just-finished powered
interval; subsequent rows use zero-voltage (shorted) electrodes. Zero voltage is
not an open-circuit model. Stored heat continues to redistribute after shutoff,
so geometry and test acceleration can continue changing without new source work.
Earlier complete CSV rows remain valid if a later advance or output write fails.

## Shared ownership and limits

The library accepts heterogeneous per-link/per-region material descriptions,
explicit global cell identities and contact schedules. Sources and thermal
conductances can be composed through the existing observer interface. Missing
link/region assignments, competing transport ownership, source-domain failures,
work-budget exhaustion and cancellation refuse the complete candidate. The
original electrical driver's final source evaluation and energy-transfer checks
run before the thermal state and complete mechanism are returned together.

`advance_joule` retains generic concrete bulk/interface source payloads. Supply
resolvers for the original bound cards and full state contexts; `regions()`
exposes those descriptions. The example binds all five properties to the same
card. The thermal audit **already includes** deposited Joule heat: adding the
report's `deposited_joule_j` again would double count it.

This remains uniform-region, solid, isotropic free-expansion mechanics on fixed
electrical/thermal reference geometry. Joint frames are prescribed. It does not
solve constrained thermal stress, evolving contact geometry, rotating electrical
connections, induction, melt flow, q/qdot integration or mechanical work during
expansion. Its thermal energy check is not a whole-system thermomechanical
conservation certificate or empirical material-accuracy claim.

Focused native checks:

```sh
cargo test -p fs-mbd-material --test heated_mechanism
cargo test -p fs-mbd-material --features electrothermal-coupling --test electrothermal_mechanism
cargo test -p fs-mbd-material --features electrothermal-coupling --example electrically_heated_rotor
```
