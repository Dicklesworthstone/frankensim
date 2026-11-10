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
| `dc_motor.json` | a DC motor across three domains: a 12 V supply, winding inductance and resistance in a current loop, a gyrator to a rotor with a viscous bearing, and both heats into one case | the monolithic midpoint rule for the two coupled equations, the supply's energy against stored energy plus heat |
| `parallel_rlc.json` | a current source drives a capacitor, inductor and resistor in parallel; all Joule heat warms the case | independent coupled voltage/current equations, source work, electrical storage and heat |
| `nonlinear_rlc.json` | the same topology with voltage-dependent capacitance and flux-dependent inductance, crossing both curves' knees | independent charge/flux differential equations, original constitutive energies and Joule heat |
| `geared_oscillator.json` | an ideal 2.5:1 transformer couples a mass to an anchored spring | reflected stiffness `ratio² × stiffness`, the midpoint oscillation phase and conserved energy |
| `quarter_car.json` | a car body and wheel on a suspension spring and shock absorber, a tyre spring driven over a 30 mm bump; 70 % of the shock's heat warms it, the rest goes to the air | the ledger: road work = mechanical energy + heat, and the 70/30 split |

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

## Ideal transformers

`geared_oscillator.json` connects the mass's effort-input port to the spring's
flow-input port. A ratio `n` enforces `flow_b = -n × flow_a` and
`effort_a = n × effort_b`, conserving interface power. It reflects the
spring's stiffness by `n²`. The initial spring extension is `-n` times the
mass position, and reversing the signed ratio reverses that extension.
This ideal transformer has no losses, backlash or contact dynamics.
