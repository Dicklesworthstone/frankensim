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

## Pause and resume a run

Save an accepted prefix, including sampled history and adaptive step memory:

    frankensim --json system examples/system/adaptive_oscillator.json --until 0.5 --checkpoint prefix.fssnap

Keep the `checkpoint.resume_root` from that report independently of the saved
file. Continue with the same scene contents and executable:

    frankensim --json system examples/system/adaptive_oscillator.json --resume prefix.fssnap --resume-root <saved-root>

`--until` is an absolute simulated time within the scene's original duration;
zero saves the admitted initial state. A resumed run can stop and save another
prefix by adding a later `--until` and a new `--checkpoint` destination. Existing
files are never overwritten. Stopping inside a proposed step adds a numerical
step boundary; resuming preserves the accepted state and its continuation.
Checkpoint files and restored sample histories are each bounded to 64 MiB.
An element whose owner does not provide a native state codec refuses saving
instead of omitting its state.
This includes support for authored or material-card-backed `enthalpy_body`
and `layer_conductor` laws; keep referenced material packs unchanged when
resuming. A `carreau_pipe` saves and resumes as well; its law is part of the
scene. `rigid_network` saves every body's exact pose and momenta;
`viscoelastic_network` also saves every Maxwell branch's full elastic memory.
Their Arrhenius and WLF calibrations are bound in full, so self-heating
feedback continues under the original temperature law. For example, the same
save/resume commands work with `self_heating_damper.json` and
`viscoelastic_mount.json` at `--until 0.125`.

| scene | what it shows | checked against (crates/fs-cli/tests/system_command.rs) |
|---|---|---|
| `driven_oscillator.json` | a mass on a spring and damper, driven by a sine force; the damper heats a housing that loses heat to a fixed ambient | the monolithic midpoint rule for the same oscillator |
| `adaptive_oscillator.json` | a loaded oscillator chooses shorter intervals to meet explicit coordinate and source-work tolerances | the analytic oscillation phase, tolerance refinement and work from the accepted displacement |
| `adaptive_suspension.json` | a three-dimensional rigid carriage coupled to a projected spring, with explicit tolerances for its motion and rotation | the analytic projected oscillator and free transverse motion, with improved phase at tighter tolerances |
| `adaptive_viscoelastic_mount.json` | a preloaded Maxwell mount with adaptive control of both rigid bodies and all six memory amplitudes | the independent two-mode motion and memory solution, routed heat, and durable continuation |
| `heated_body.json` | a 30 W heater in a body that loses heat through a wall to 293.15 K | the closed-form discrete approach to 308.15 K |
| `heated_plane.json` | 30 W generated on a plane of no heat capacity (`heat_node`) between a 50 J/K pad at 400 K and a 120 J/K disc at 300 K, joined to the plane by 2 and 6 W/K: no fraction is given, the heat divides by where the two temperatures put the plane | the midpoint rule for the two bodies and the plane's balance, solved independently each step: both temperatures and the plane's to 1e-9, every joule in the two bodies; the same device with the heat rate prescribed on the plane |
| `friction_belt.json` | a block follows a slowly accelerating belt by static friction, breaks away when the acceleration exceeds the static limit, and catches it again as the lining warms | independent constant-acceleration segments: 3 m/s final speed, 2.875 m travel, 3 J of frictional heat and 12 J supplied by the belt |
| `stribeck_block.json` | a 2 kg block thrown at 3 m/s on a Stribeck lining (`friction_interface`, a dry law of fs-tribo): the coefficient is 0.35 at 3 m/s and grows toward 0.5 as the block slows; all 9 J become the lining's heat | the stopping distance from Simpson's rule on the law's own integral, the error falling about fourfold when the step halves; heat equal to the initial kinetic energy |
| `dragged_lining.json` | a Stribeck lining (`sliding_interface`) dragged at a prescribed 0.5 m/s under 40 N; its heat warms a 20 J/K pad that loses heat through a wall to 293.15 K | the law's traction mu(0.5) N along the slip and its work; the pad against the midpoint rule's closed recursion to 1e-10 K; the pad's heat plus the wall's equal to the work |
| `spun_puck.json` | a 0.5 kg puck (`rigid_network`) sliding at 2 m/s and spinning at 100 rad/s on a 1 cm `hertz_patch` under 5 N, mapped from its force and torque ports; its heat warms the ice | the last step's force and moment against a quadrature of local Coulomb friction over the Hertz pressure to 1e-3; both below their pure-slide and pure-spin values; the lost kinetic energy equal to the ice's heat |
| `bowed_mode.json` | one undamped 2 Hz mode of a string (`modal_structure`) rubbed by a bow (the interface's prescribed end) at 0.1 m/s through Coulomb friction (`friction_interface` reading the end of the string's flow): kinetic friction pushes the mode until it catches the bow, then it rides stuck | the closed forms: the catch at `sin w t1 = V w / Fk` and the heat of the slide; stuck, the position as the held effort gives it, `2 V tan(w s / 2) / w` per step of length `s` (2.2e-7 m past the continuous ride at 0.2 s, the held force's work on that gap booked as the interface's numerical loss); a 0.5 m/s bow, never caught, against the sliding solution |
| `dc_motor.json` | a DC motor across three domains: a 12 V supply, winding inductance and resistance in a current loop, a gyrator to a rotor with a viscous bearing, and both heats into one case | the monolithic midpoint rule for the two coupled equations, the supply's energy against stored energy plus heat |
| `parallel_rlc.json` | a current source drives a capacitor, inductor and resistor in parallel; all Joule heat warms the case | independent coupled voltage/current equations, source work, electrical storage and heat |
| `nonlinear_rlc.json` | the same topology with voltage-dependent capacitance and flux-dependent inductance, crossing both curves' knees | independent charge/flux differential equations, original constitutive energies and Joule heat |
| `carreau_pipe.json` | a shear-thinning liquid (an illustrative Carreau-Yasuda law, no measured material) driven through a 5 mm pipe by 10 kPa; all of the viscous heat warms the liquid's 0.2 J/K, whose Arrhenius shift thins it, so the flow rises as it warms | the Rabinowitsch-Mooney integral of the same law at the temperature the pipe read, to 1e-7; supply work = pipe heat = the liquid's stored heat |
| `geared_oscillator.json` | an ideal 2.5:1 transformer couples a mass to an anchored spring | reflected stiffness `ratio² × stiffness`, the midpoint oscillation phase and conserved energy |
| `projected_suspension.json` | a real 3-D rigid body suspended by a scalar spring along an oblique direction through a rectangular power map | the projected oscillator, free transverse motion, signed/scaled transmissions and energy closure |
| `anisotropic_damping.json` | a 3-D rigid body slowed by a full damping tensor; off-diagonal terms induce transverse motion and every lost joule heats the case | independently diagonalized decay modes, position, last-step forces and kinetic-to-thermal energy transfer |
| `port_hamiltonian_oscillator.json` | the driven oscillator written as port-Hamiltonian matrices and stepped by the exact held-input owner of fs-system; its damping work warms a case | the closed form of `m q'' + c q' + k q = F` to 1e-12 of scale after 500 steps (2e-16 measured), where the element-built oscillator is about 3e-6 m off; supplied work = stored energy + heat to 1e-12 J |
| `rack_and_pinion.json` | one pinion inertia carrying a motor torque and a bearing directly, and a rack's spring, friction and load through the pinion radius: three members of one join, each behind its own ratio | the midpoint rule for the pinion's one equation with the reflected stiffness and friction; motor work, the load's received work and the heat |
| `pneumatic_cylinder.json` | a 2 l tank of air at 6 bar feeds a cylinder through a valve; the piston pushes a spring and the outside pressure out by 104 mm; the tank cools to 281 K as it empties, the chamber warms to 325 K as it fills and gives heat to the barrel, which also takes the seal's friction | an independent RK4 integration of the cylinder's seven equations (4.9e-9 at 3 s, second order in mid-stroke), conservation of the gas, and the force balance at rest |
| `n2_gas_spring.json` | a 1 kg piston rings on 0.1 L of insulated N2 (`gas_cavity`, the gas of the seed McBride NASA-9 and species packs in `n2/`) held at 90 kPa by 90 N | the adiabatic period from the NASA-9 γ at 300 K with its amplitude and step corrections, to 1e-6; the ledger closes |
| `damped_pair.json` | two rigid bodies in three dimensions closing at 2 m/s on a viscous connector, stepped by the fs-mbd network owner (`rigid_network`); a quarter of the connector's heat warms a shell, half an oil bath, the rest goes to the air | the closed form: each momentum decays as `exp(-0.7 t)`, the heat is `2 (1 - exp(-1.4 t))` J, split 25 / 50 / 25 |
| `damped_pair_lining.json` | the damped pair with all of its connector heat routed into the first of three steel cells of 5 J/K (`thermal_timeline`, the fs-thermal enthalpy timeline), joined in a row by 1 W/K: the heat enters the first cell and spreads along the lining | the closed form: the lining stores the pair's heat `2 (1 - exp(-1.4 t))` J to 1e-9 J, the cells warmest first |
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

`stribeck_block.json` uses `friction_interface` instead, which steps a dry
law of fs-tribo: `coulomb`, `velocity_dependent`, `stribeck`, `hyperbolic`,
`shifted_power`, or the measured rosin fit by its name,
`"rosin_woodhouse_galluzzo_2025_eq13"`. Its ports and `slip_m_s` read as
the slider's. The law's coefficient here is 0.2, plus 0.3 decaying as a
Gaussian of the speed over 0.5 m/s, plus 0.05 per m/s. The block stops where
`(m / N) ∫ v / μ(v) dv` from 0 to 3 m/s puts it: each sliding step takes
the law's traction at the step's mean slip, and the step the stop falls in
is shortened to end at it. A moving surface under a `friction_interface`,
such as a belt, is a prescribed end (`{"prescribe": "friction.b",
"velocity_m_s": 1.0}`) or a `velocity_source` joined to the end. A
`friction_interface` has no span; either end is held over each step, so
give an accelerating end short steps. Fixed `step_s` only.

`bowed_mode.json` rubs a `modal_structure` with one. A modal structure
steps its modes exactly under the effort held over a step, so the end of
its flow does not follow the midpoint relation; each attachment tells that
end, and the interface reads it on the side it names in `reads_ends`
(`["a"]` here). The kernel wires the end input and refuses an interface
beside such a neighbour that does not read it. The string's modes,
damping and mode-shape weights are authored; nothing is solved from its
geometry.

A `modal_structure` can instead give a lumped `network`: a stiffness
and a mass over its degrees of freedom, and optional Rayleigh damping
`[a, b]`. fs-nlmodal solves its modes, and the attachments and the
`start` are written per degree of freedom:

```json
"chain": {"kind": "modal_structure", "sample_rate_hz": 1000,
          "network": {"stiffness": [[600.0, -200.0], [-200.0, 200.0]],
                      "mass": [[1.0, 0.0], [0.0, 0.5]],
                      "rayleigh": [0.5, 0.001]},
          "start": {"displacement": [0.0, 0.01]},
          "attachments": {"tip": {"domain": "translation",
                                  "weights": [[0.0, 1.0]]}}}
```

`dragged_lining.json` uses `sliding_interface`: the same laws, with the
motion prescribed instead of solved. A lining is dragged at a constant
`slip_rate_m_s` of [0.3, 0.4] m/s under 40 N, and nothing in the scene
moves. The interface turns that motion's work into heat on its `heat`
port, joined to the pad as a heater would be, and reports the traction
from its signal output `traction`. The work is the interface's own
account, so the ledger closes with the pad's heat and the wall's.

`spun_puck.json` rubs a rigid body on a `hertz_patch`: one circular patch
whose slip and spin friction fs-tribo solves together. Its `wrench` port
takes `[v_x, v_y, radius_m w]`, the patch's slip and its radius times its
spin. A `power_map` of members builds that from the puck's `force0` (its
velocity) and `torque0` (its spin), each through its own matrix, negated
because the map hands the patch minus the mapped flow:

```json
{"power_map": [[{"port": "puck.force0",
                 "matrix": [[-1, 0, 0], [0, -1, 0], [0, 0, 0]]},
                {"port": "puck.torque0",
                 "matrix": [[0, 0, 0], [0, 0, 0], [0, 0, -0.01]]}],
               "patch.wrench"]}
```

The patch sits at the puck's centre of mass, so its force does not tip
the puck; a contact arm `r` adds `r x t1` and `r x t2` to the torque
port's rows. A spinning patch slides more easily, and a sliding one spins
more easily, than either alone. The run stops before the puck does: a
step in which the slip or the spin would reverse is refused.

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
every connector's mean heat rate over the last step. Both network types can
use error-controlled stepping with the complete native accounts described
below; `adaptive_suspension.json` demonstrates a rigid carriage and spring.

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

`damped_pair_lining.json` gives the pair's heat to a second physics owner:
a `thermal_timeline` lowers to `fs_thermal::component::ThermalTimeline`,
the enthalpy timeline of fs-thermal, whose cells, links and boundaries are
stepped by that owner's explicit rule. Each cell named in `ports` has a heat
port `cell<id>`, joined or routed to like a heat capacity's `heat`. A
boundary that exchanges heat (a film, a reservoir or a source) needs the
port `environment`, whose temperature is prescribed and not read; its
ledger row is the heat that crossed the boundaries. `events` make and break
links and replace boundaries at later times; every cell is born at the
start. Fixed `step_s` only; the kernel also cuts steps to the owner's
stability bound and to event times. No checkpoint codec.

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

For a `rigid_network` or `viscoelastic_network`, declare these accounts for
every body, prefixed by its element name:

| Account | Error controlled | Absolute unit |
|---|---|---|
| `body<i>.position_m` | difference of the complete world-position vectors | m |
| `body<i>.velocity_m_s` | difference of the complete world-velocity vectors | m/s |
| `body<i>.angular_velocity_rad_s` | difference of body-frame angular-velocity vectors | rad/s |
| `body<i>.orientation_rad` | shortest rotation between coarse and refined orientations | rad |
| `connector<i>.branch<j>.amplitude_sqrt_j` | difference of all six signed amplitudes of each viscoelastic branch | square root of J |

For example, the mount's first branch needs
`mount.connector0.branch0.amplitude_sqrt_j`. Its energy alone is insufficient:
opposite memory amplitudes have the same energy and different force responses.
All bodies and all branches need their accounts, even if some coordinates
start at zero. `adaptive_viscoelastic_mount.json` provides a complete example.
Vector errors use vector lengths for relative scaling. Rotation uses angular
motion during the trial interval, so its tolerance does not depend on the
world frame or a quaternion's sign. State quaternions are never renormalized by
the error controller. Temperature feedback still requires the order of the
complete coupled method; the start-of-step feedback examples use order 1.

The driver compares a full step with two half steps and retains the audited
two-half-step solution when its normalized estimate is at most one. Each
scalar coordinate uses `absolute + relative * max(|initial|, |refined|)` as
its scale; vector and rotation accounts use the scales above. The maximum
normalized coarse/refined difference is divided by `2^order - 1`.
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

## Liquid circuits

The parallel forms carry a liquid as well, with volume, flow and pressure in
place of charge, current and voltage: an `accumulator`
(`capacitance_m3_pa`, `volume_m3`) sets the pressure of a `common_effort`
node, a `hydraulic_resistor` (`resistance_pa_s_m3`) passes
`(p_a - p_b) / R` and sends its heat out of `heat`, a `flow_source`
(`flow_m3_s`, positive into the source) supplies or draws flow, and a
`pressure` source (`pressure_pa`) holds a node's pressure. A far end is
held with `{"prescribe": port, "pressure_pa": 0}`. These are the same kernel
elements as the electrical ones; a scene in either set of names gives the
same numbers.

`carreau_pipe.json` runs the fs-material laminar pipe of a Carreau-Yasuda
liquid. Its viscosity at the shift's reference temperature is
`mu_inf + (mu_0 - mu_inf) (1 + (time_constant_s rate)^a)^((n - 1) / a)`; the
owner solves the stress over the radius for the flow, so the pipe is a
nonlinear resistance between `a` and `b`. Its `heat` port reads the liquid's
temperature and returns all of the viscous heat; the `shift` (Arrhenius or
WLF, as for a network) multiplies the viscosity, and with
`"shift_mode": "vertical-and-rate"` also the rate the curve is read at.
Joined to a `heat_capacity`, the heat warms the liquid and the warmer liquid
flows faster, solved within each step. The law in the example is
illustrative, not a measured liquid. The pipe holds no liquid and has no
inertia, entrance loss or turbulence; a wall rate above
`maximum_wall_rate_per_s` or outside `rate_support_per_s`, or a temperature
outside the shift's `valid_k`, refuses the step.

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

`n2_gas_spring.json` uses a `gas_cavity` instead: 0.1 L of N2 at 300 K and
90 kPa whose heat capacities, pressure and the states it may reach come
from cited sources. These are the seed species, McBride (2002) NASA-9 and
Lemmon (2004) transport cards, compiled by `xtask matdb-pack` into the packs
in `n2/`. The cavity is insulated, and 90 N on a 1 kg piston of 10 cm²
holds it at rest; the piston starts at 5 cm/s and rings at the adiabatic
frequency, with γ from the NASA-9 heat capacity at 300 K. A mass moving the
positive way enlarges the cavity, and a prescribed `flow_m3_s` on `p` is the
compression flow (positive shrinks it). The cavity is closed and frozen
(no reactions, condensation or leakage). Every step resolves the packs at
its end state, so a run that leaves a card's validity (the transport card's
75 to 101.325 kPa here) is refused at that step. The caloric bracket must
stay inside one NASA-9 chart (below 1000 K for N2).

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

## Port-Hamiltonian matrices

`port_hamiltonian_oscillator.json` writes the driven oscillator as matrices
over `x = [momentum, displacement]`: `energy` is `[[1/m, 0], [0, k]]`,
`structure` `[[0, -1], [1, 0]]`, `dissipation` puts `c` on the momentum row,
and `input` takes the force into the momentum. A held 2 N force drives it
from rest, and the dissipated work warms a 10 J/K case.

The element is the exact held-input owner of fs-system: for an input held
over a step, the step is exact. After 500 steps of 1 ms its state agrees with
the closed form of `m q'' + c q' + k q = F` to about 2e-16 (the test asserts
1e-12 of scale), while the same oscillator built from `mass`, `spring` and
`damper` (the midpoint rule) is about 3e-6 m off. The supplied work, the stored energy and the case's heat
close to 1e-12 J. Inputs that come from other elements are the kernel's step
means, so a coupled scene is not exact.

`energy`, `structure` and `dissipation` are n by n with n from 1 to 64:
`energy` and `dissipation` symmetric positive semidefinite, `structure` skew.
`input` has one column per port. Each port names a `domain` (`translation`,
`rotation` or `electrical`) and a `role`: an `effort_in` port takes an effort
and returns the step mean of its output as a flow; a `flow_in` port takes a
flow and returns an effort. A one-state `flow_in` element with `energy`
`[[k]]` and zero `structure` and `dissipation` follows a `spring` under the
same end motion. Reports give `x0`, `x1`, ..., each port's returned mean
(`p.velocity_m_s` here) and `heat_w`; adaptive runs name `oscillator.x0`,
`oscillator.x1`, ... as state tolerances. The storage is quadratic only, and
the element has no checkpoint codec.
