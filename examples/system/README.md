# Devices described as data

Each file here is a `frankensim.system.v0` scene: elements from
`fs_system::elements` or `fs_system::constitutive` and the connections between
their ports. Run one with

    frankensim system examples/system/quarter_car.json
    frankensim --json system examples/system/quarter_car.json

`frankensim system --help` lists the elements, their ports and the
connections. Nothing in the CLI integrates anything: the scene is lowered to
a `fs_system::System` and stepped by the composition kernel, and the report
is the kernel's state and its energy and entropy ledger.

| scene | what it shows | checked against (crates/fs-cli/tests/system_command.rs) |
|---|---|---|
| `driven_oscillator.json` | a mass on a spring and damper, driven by a sine force; the damper heats a housing that loses heat to a fixed ambient | the monolithic midpoint rule for the same oscillator |
| `adaptive_oscillator.json` | a loaded oscillator chooses shorter intervals to meet explicit coordinate and source-work tolerances | the analytic oscillation phase, tolerance refinement and work from the accepted displacement |
| `heated_body.json` | a 30 W heater in a body that loses heat through a wall to 293.15 K | the closed-form discrete approach to 308.15 K |
| `heated_plane.json` | 30 W generated on a plane of no heat capacity (`heat_node`) between a 50 J/K pad at 400 K and a 120 J/K disc at 300 K, joined to the plane by 2 and 6 W/K: no fraction is given, the heat divides by where the two temperatures put the plane | the midpoint rule for the two bodies and the plane's balance, solved independently each step: both temperatures and the plane's to 1e-9, every joule in the two bodies; the same device with the heat rate prescribed on the plane |
| `friction_belt.json` | a block follows a slowly accelerating belt by static friction, breaks away when the acceleration exceeds the static limit, and catches it again as the lining warms | independent constant-acceleration segments: 3 m/s final speed, 2.875 m travel, 3 J of frictional heat and 12 J supplied by the belt |
| `dc_motor.json` | a DC motor across three domains: a 12 V supply, winding inductance and resistance in a current loop, a gyrator to a rotor with a viscous bearing, and both heats into one case | the monolithic midpoint rule for the two coupled equations, the supply's energy against stored energy plus heat |
| `parallel_rlc.json` | a current source drives a capacitor, inductor and resistor in parallel; all Joule heat warms the case | independent coupled voltage/current equations, source work, electrical storage and heat |
| `nonlinear_rlc.json` | the same topology with voltage-dependent capacitance and flux-dependent inductance, crossing both curves' knees | independent charge/flux differential equations, original constitutive energies and Joule heat |
| `geared_oscillator.json` | an ideal 2.5:1 transformer couples a mass to an anchored spring | reflected stiffness `ratio² × stiffness`, the midpoint oscillation phase and conserved energy |
| `projected_suspension.json` | a real 3-D rigid body suspended by a scalar spring along an oblique direction through a rectangular power map | the projected oscillator, free transverse motion, signed/scaled transmissions and energy closure |
| `anisotropic_damping.json` | a 3-D rigid body slowed by a full damping tensor; off-diagonal terms induce transverse motion and every lost joule heats the case | independently diagonalized decay modes, position, last-step forces and kinetic-to-thermal energy transfer |
| `rack_and_pinion.json` | one pinion inertia carrying a motor torque and a bearing directly, and a rack's spring, friction and load through the pinion radius: three members of one join, each behind its own ratio | the midpoint rule for the pinion's one equation with the reflected stiffness and friction; motor work, the load's received work and the heat |
| `pneumatic_cylinder.json` | a 2 l tank of air at 6 bar feeds a cylinder through a valve; the piston pushes a spring and the outside pressure out by 104 mm; the tank cools to 281 K as it empties, the chamber warms to 325 K as it fills and gives heat to the barrel, which also takes the seal's friction | an independent RK4 integration of the cylinder's seven equations (4.9e-9 at 3 s, second order in mid-stroke), conservation of the gas, and the force balance at rest |
| `damped_pair.json` | two rigid bodies in three dimensions closing at 2 m/s on a viscous connector, stepped by the fs-mbd network owner (`rigid_network`); a quarter of the connector's heat warms a shell, half an oil bath, the rest goes to the air | the closed form: each momentum decays as `exp(-0.7 t)`, the heat is `2 (1 - exp(-1.4 t))` J, split 25 / 50 / 25 |
| `viscoelastic_mount.json` | two rigid bodies at rest on a mount with viscous damping and one Maxwell branch stretched by 0.5 m (`viscoelastic_network`, the fs-mbd relaxation owner): the stored 0.21 J sets the bodies moving and heats the casing | the pair's own two equations integrated independently: relative speed, branch energy and heat |
| `self_heating_damper.json` | the damped pair with all of its heat in 1 J/K of damper fluid, whose temperature thins the damper by an Arrhenius shift: the damping falls to 0.81 of its value as the fluid warms by 0.93 K | the exact recursion of the scheme (damping held over each step at the temperature the step began with) to 1e-11, and first-order convergence to the coupled differential equations |
| `quarter_car.json` | a car body and wheel on a suspension spring and shock absorber, a tyre spring driven over a 30 mm bump; 70 % of the shock's heat warms it, the rest goes to the air | the ledger: road work = mechanical energy + heat, and the 70/30 split |

## Dry friction and prescribed motion

`friction_belt.json` puts a 2 kg block on a belt under 20 N of normal
force. The static coefficient is 0.5 and the kinetic coefficient is 0.2.
At first the belt accelerates at 4 m/s², requiring 8 N, so the block stays
stuck. The second ramp needs 16 N, exceeding the 10 N static limit; the
block slips with 4 N of friction. Once the belt holds 3 m/s, the block
catches it at 1.25 s. At 1.5 s the block has travelled 2.875 m, the lining
has received 3 J and warmed by 0.06 K, and the belt has supplied 12 J.
These are outputs of the coupled equations, not inputs in the scene.

The `slider` has power ports `a` and `b` and a thermal port `heat`.
`slip_m_s`, when needed initially, is velocity at `a` minus velocity at
`b`; it must agree with the connected bodies' initial velocities. Zero
starts stuck. The report gives the same physical slip, a numeric `stuck`
flag (1 or 0), force on side `a`, and the last step's heat rate.

An accelerating prescribed end needs both its step mean and endpoints.
`velocity_source` publishes those through `p` and the signal `span`.
The slider declares `"prescribed_ends": ["b"]`, creating `span_b`, and
the scene connects the belt to `friction.b` and its signal to `span_b`.
A direct boundary can instead prescribe `velocity_m_s` on the end and
the same waveform on its span; the span takes endpoint values. A direct
boundary publishes no event schedule, so its table knots must align with
steps; `velocity_source` lets the kernel land on those knots. An untold
prescribed end is held during each step. Error-controlled stepping does
not cover sliders; use a fixed `step_s`. Coefficients and normal force
are constant authored inputs. This is ideal dry friction without contact
geometry, impact or material-card claims, and switching steps are first
order.

## Rigid bodies from a physics owner

`damped_pair.json` is the first scene whose element is not a lumped law of
`fs_system::elements` but a physics owner behind the Component contract: a
`rigid_network` lowers to `fs_mbd::component::network::ViscousNetwork`, which
steps rigid bodies in three dimensions and their viscous connectors with the
fs-mbd owner's own splitting. Its ports carry three components: `force<i>`
takes the world force on body `i` and `torque<i>` (for the bodies named in
`torque_ports`) the torque in that body's axes, each prescribed with an array
of three waveforms; `heat<e>` is routed like a damper's heat. The report
gives every body's position, velocity, orientation and angular velocity and
every connector's mean heat rate over the last step. Error-controlled
stepping does not cover a `rigid_network`; a scene with one runs at a fixed
`step_s`.

`viscoelastic_mount.json` uses the second fs-mbd owner,
`fs_mbd::component::relaxation::RelaxationNetwork`, through a
`viscoelastic_network`. Its connectors may carry Maxwell `branches`: a
stiffness on each axis, a `relaxation_time_s` (without one the branch never
relaxes) and a starting `extension_m` or `twist_rad`. A branch is a rate law
in its first body's axes: it stores the integral of the relative velocity
it has seen, less what has relaxed, and is not a spring with a rest
geometry. The report adds each branch's stored energy,
`connector<e>.branch<b>.energy_j`. A connector with `"drive": "body"` or
`"world"` has an actuator: ports `drive<e>` and `couple<e>` take the force
and couple it holds on the connector's second attachment, three waveforms
each, in the first body's axes or in world axes, and the work it does is
that boundary's account in the ledger.

## Heat generated between bodies

`heated_plane.json` puts a heat source on a `heat_node`: a thermal node that
stores nothing, joined to each of its bodies by a conductance
(`conductances_w_k`, one per body; ports `body0`, `body1`, ...), with a port
`source` for the heat put into the node. Its temperature is fixed by its
balance, `T = (P + sum G_k T_k) / sum G_k`, and each body receives
`G_k (T - T_k)`. A hot body therefore gives heat up through the node while a
cold one takes the source's and its neighbour's, with no fraction written
anywhere. This is the plane a sliding interface generates its heat on, a
thin heater between two parts, or a glue line.

- Join each `body<k>` to a `heat_capacity` (or prescribe its
  `temperature_k`), and join a `heater` to `source` or prescribe
  `heat_rate_w` there. With `plane.body1` prescribed at 300 K in place of
  the disc, the pad settles at `300 K + P / G1 = 305 K`, where it exchanges
  nothing, and goes on stepping there.
- `temperature_k` on the element is optional and is not a state. With a
  heater on `source` the node's temperature is solved in a loop, and this is
  where that solve starts. Every part of this scene's loop is linear, and
  its report is the same to the byte from 250, 301, 325, 600 and 5000 K; in
  a loop with a part that is not linear the result would depend on the
  start within the solve's tolerance. Without it such a scene is admitted
  and its first step refused (exit 4), naming the heater's port.
- The node reports `temperature_k`: what its balance gave in the last step.

## Temperature feedback

`self_heating_damper.json` closes a loop between two owners with nothing
written for the pair: the damper's heat warms a `heat_capacity`, and that
temperature sets the damper's viscosity.

- Every `heat_capacity` has a signal output `temperature`: its temperature at
  the start of the step. `{"signal": "fluid.temperature", "to": [ports]}`
  copies it to signal inputs.
- A network with a `shift` has one signal input `temperature<e>` per
  connector. The shift is `a_T = tau(T) / tau_ref` from the fs-material
  owner, as `{"arrhenius": {activation_energy_j_mol, reference_k, valid_k}}`
  or `{"wlf": {c1, c2_k, reference_k, valid_k}}`; a temperature outside
  `valid_k` refuses the step. In a `rigid_network` it multiplies each
  connector's damping; in a `viscoelastic_network` it multiplies the damping
  and every relaxation time, and leaves the stiffnesses alone.
- A `temperature<e>` input can also be held from outside with
  `{"prescribe": port, "temperature_k": T}`, for an isothermal run away from
  the reference temperature.

The damping is held over each step at the temperature the step began with,
so the feedback is first order in the step: halving `step_s` halves the
distance to the solution of the coupled equations. The heat itself is still
accounted exactly by the ledger at any step.

## Error-controlled stepping

`adaptive_oscillator.json` starts with a maximum interval of 0.5 s and asks
the kernel to choose the interval from explicit tolerances. Its 1 kg mass and
4 N/m spring oscillate about the 0.1 m equilibrium set by a 0.4 N force. The
frequency is therefore 2 rad/s. Tightening the tolerances improves the phase
relative to this independent analytic solution while the energy ledger still
closes.

The optional top-level `adaptive` object supplies `order`, `min_step_s`,
`max_attempts` and a `tolerances` object keyed by `element.state_field`:

```json
"adaptive": {
  "order": 2,
  "min_step_s": 1e-12,
  "max_attempts": 32,
  "tolerances": {
    "block.position_m": {"absolute": 1e-5, "relative": 0},
    "block.velocity_m_s": {"absolute": 1e-5, "relative": 0},
    "spring.extension_m": {"absolute": 1e-5, "relative": 0},
    "drive.energy_in_j": {"absolute": 1e-6, "relative": 0}
  }
}
```

Every dynamic coordinate needs its own tolerance: both coordinates of a mass,
inertia or linear inductor; spring extension or capacitor charge; nonlinear
inductor flux; thermal-body temperature; and each source's `energy_in_j`.
Source tolerances apply to work during the trial interval, so earlier supplied
energy does not loosen a later interval's budget. Derived voltage/current
probes and stateless last-step means cannot substitute for those coordinates.
Absolute tolerances are positive in the named field's SI units; relative
tolerances are finite and nonnegative.

The driver compares a full step with two half steps and retains the audited
two-half-step solution when its normalized estimate is at most one. Each
coordinate uses `absolute + relative * max(|initial|, |refined|)` as its
scale; the maximum coarse/refined difference is divided by `2^order - 1`.
The default order is 1; choose order 2 only for a smooth second-order model,
as in this oscillator. Estimates remain conditional at curve knees and
discontinuous inputs; a lower order assumption does not establish an error
bound there. These are local estimates, not global trajectory bounds.

The `step_s` ceiling still applies. Only an externally capped final interval
may be smaller than `min_step_s`; its own accuracy retries cannot reduce that
interval further. Rejected trials contribute no accepted state, source work
or ledger entries. JSON reports an `adaptive` object with the order, total
attempts and trial steps, last/maximum accepted normalized error, and the next
suggested interval. The ordinary step count includes the two accepted half
steps per adaptive interval. Error control covers the declared coordinates
and source-work increments; heat and interface balances retain their separate
ledger audits. Omitting `adaptive` preserves fixed-step behavior and output.

## Parallel electrical networks

`parallel_rlc.json` uses a `node_capacitor` as the common-voltage node. An
explicit `common_effort` connection joins its `p` port to the inductor,
`resistor_branch` and `current_source`. The resistor branch's other terminal
is grounded by prescribing `voltage_v: 0`. Its `heat` port must also be
connected; the example routes every joule to the case.

`node_capacitor.charge_c` is the charge on the node, so its endpoint voltage is
`charge_c / capacitance_f`. Current is positive into each component, including
the source: the negative source waveform in the example supplies
`0.4 + 0.2 sin(2π × 3t + 0.3)` amperes to the node. The existing `capacitor`
and `resistor` elements remain the forms used in series current loops.

## Nonlinear electrical storage

`nonlinear_rlc.json` supplies a `nonlinear_capacitor` with a
`charge_voltage_curve` of `[charge_c, voltage_v]` pairs, and a
`nonlinear_inductor` with a `flux_current_curve` of `[flux_wb, current_a]`
pairs. These are illustrative authored laws, with no measured-material claim.
The capacitor's differential capacitance falls from 0.01 F to 0.005 F beyond
0.02 C in either polarity. The inductor's differential inductance falls from
0.05 H to 0.025 H beyond 0.05 Wb. The initial state excites both regions, and
the resistor's heat flows into the case while the current source supplies
0.05 A.

Each curve needs 3–4096 knots, with both coordinates strictly increasing and
an interior `(0,0)` knot. The existing `fs-phs` constitutive owner retains the
knots and exact piecewise integrals. The capacitor stores `integral V dq` and
the inductor `integral I d(flux)`. A local differential capacitance or
inductance multiplied by a squared endpoint value is not their energy.
States outside the supplied curve support refuse by name; the curves never
extrapolate or clamp. They do not model hysteresis, dielectric loss or
temperature dependence.

Both elements have a single `p` port. The capacitor receives entering current
and stores `charge_c`; the inductor receives voltage and stores `flux_wb`.
Their reported `voltage_v` and `current_a` are endpoint constitutive values.
The kernel uses integral mean efforts for work during a finite step, retaining
energy when a step crosses knots. The original constant `capacitor`,
`node_capacitor` and `inductor` scene forms keep their existing conventions.

## Gas in volumes

`pneumatic_cylinder.json` has two `gas_volume` elements and a `restriction`
between them. A volume holds a perfect gas (a gas constant and a heat
capacity ratio, here air) and is the node its streams meet at:

```json
{"join": ["tank.stream", "valve.a"]},
{"join": ["chamber.stream", "valve.b"]}
```

The valve passes `conductance × (p_a − p_b)` kg/s. What crosses carries the
enthalpy of the side it leaves, so the gas arriving in the chamber brings
the tank's `c_p T` with it and the chamber warms as it fills (to 325 K in
mid-stroke) while the tank cools as it expands (to 281 K). Nothing in the
scene says so. Several restrictions may meet at one volume; two volumes do
not join directly.

The chamber has `"wall": true`, which gives it a port `p`: the pressure on a
moving wall. It joins the piston as a member with the face area as its ratio,
beside the outside pressure (a constant force), the load spring and the seal:

```json
{"join": ["piston.p", {"port": "chamber.p", "ratio": 2.0e-3},
          "outside.p", "load.a", "seal.a"]}
```

A volume's `heat` port is a node like a heat capacity's: here a film
conductance joins the chamber's gas to the barrel, which also takes the
seal's friction heat. The tank is insulated (`heat_rate_w: 0`).

A volume may carry other quantities with its gas (`"carried": [0.2, 0.0]`,
mass fractions, and `"carried": 2` on each restriction): they move with the
flow by the same rule. From outside, `{"prescribe": "vessel.stream",
"stream": [mass rate, enthalpy rate, carried rates...]}` feeds a volume, and
`{"prescribe": "inlet.a", "contents": [pressure, specific enthalpy,
fractions...]}` holds a restriction's far end at a supply's state.

All joined volumes must hold the same gas; this is not checked. The gas has
constant heat capacities. A very small volume between restrictions forces
small steps: no more than about twice its mass can pass through it in one.

## Ideal transformers

`geared_oscillator.json` connects the mass's effort-input port to the spring's
flow-input port. A ratio `n` enforces `flow_b = -n × flow_a` and
`effort_a = n × effort_b`, conserving interface power. It reflects the
spring's stiffness by `n²`. The initial spring extension is `-n` times the
mass position, and reversing the signed ratio reverses that extension.
This ideal transformer has no losses, backlash or contact dynamics.

A transformer joins two ports. When one node carries several loads at
different ratios, write the ratio on the member of the join instead:

```json
{"join": ["pinion.p", "motor.p", "bearing.a",
          {"port": "rack_spring.a", "ratio": 0.02},
          {"port": "rack_friction.a", "ratio": 0.02},
          {"port": "load.p", "ratio": 0.02}]}
```

`rack_and_pinion.json` is this join. The pinion's inertia sets the angular
velocity `w`. A member with ratio `r` receives the flow `-r w` (here the
rack's velocity for a radius of 0.02 m) and its effort returns to the node
multiplied by `r`, so a 4000 N/m spring is felt as `r² k = 1.6` N m/rad
and a -12 N load as -0.24 N m. Members without a ratio share the node's
flow unchanged. Members may be in other power domains than the node; heat
ports take no ratio, and neither does the node port. The pinion settles at
`(0.6 - 0.24) / 1.6 = 0.225` rad, 4.5 mm of rack travel. The same form
works in a `common_effort` list.

## Rectangular power maps

`projected_suspension.json` connects a 2 kg rigid body to an 18 N/m scalar
spring along the fixed world direction `[0.6, -0.8, 0]`:

```json
{"power_map": ["carriage.force0", "spring.a"], "matrix": [[0.6, -0.8, 0.0]]}
```

Rows follow the second port and columns the first. The spring receives minus
the projected velocity; its force returns through the transpose, preserving
virtual work. Its initial extension is minus the projected position. The
body oscillates along that direction and moves freely in the perpendicular
plane. Doubling the map reflects four times the stiffness; entries are not
normalized. Matrices must have the exact port widths and finite entries.
Heat, stream and signal ports refuse. A fixed matrix neither tracks a
rotating attachment direction nor infers geometry.

## Matrix resistance and damping

`anisotropic_damping.json` starts a 2 kg body moving in the x-z plane and
applies the damping matrix `[[4, 1, 0], [1, 4, 0], [0, 0, 2]]` N s/m through
an identity power map. The off-diagonal terms produce a y velocity despite
zero initial y motion. Its independent decay rates are 2.5/s and 1.5/s along
the two diagonals in the x-y plane, and 1/s along z. Every joule removed from
motion warms the 10 J/K case, starting at 300 K.

The `resistance_matrix` element declares `domain` as `translation`,
`rotation` or `electrical` and a square `matrix` of 1 to 64 coordinates.
Coefficients are respectively in N s/m, N m s/rad or ohm. The existing
resistance owner requires finite entries, exact symmetry and positive
definiteness. Its `u` port receives the vector flow and returns the vector
effort; its `heat` port must be routed or prescribed. Last-step probes
`u0.force_n`, `u1.force_n`, ... (or `torque_n_m` / `voltage_v` for the other
domains) and `heat_w` report the accepted mean response. The element stores
no energy and needs no adaptive state tolerance. These are authored constant
coefficients, with no material or temperature law inferred.
