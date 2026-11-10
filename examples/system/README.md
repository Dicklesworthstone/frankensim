# Devices described as data

Each file here is a `frankensim.system.v0` scene: elements from
`fs_system::elements` and the connections between their ports. Run one with

    frankensim system examples/system/quarter_car.json
    frankensim --json system examples/system/quarter_car.json

`frankensim system --help` lists the elements, their ports and the
connections. Nothing in the CLI integrates anything: the scene is lowered to
a `fs_system::System` and stepped by the composition kernel, and the report
is the kernel's state and its energy and entropy ledger.

| scene | what it shows | checked against (crates/fs-cli/tests/system_command.rs) |
|---|---|---|
| `driven_oscillator.json` | a mass on a spring and damper, driven by a sine force; the damper heats a housing that loses heat to a fixed ambient | the monolithic midpoint rule for the same oscillator |
| `heated_body.json` | a 30 W heater in a body that loses heat through a wall to 293.15 K | the closed-form discrete approach to 308.15 K |
| `dc_motor.json` | a DC motor across three domains: a 12 V supply, winding inductance and resistance in a current loop, a gyrator to a rotor with a viscous bearing, and both heats into one case | the monolithic midpoint rule for the two coupled equations, the supply's energy against stored energy plus heat |
| `parallel_rlc.json` | a current source drives a capacitor, inductor and resistor in parallel; all Joule heat warms the case | independent coupled voltage/current equations, source work, electrical storage and heat |
| `geared_oscillator.json` | an ideal 2.5:1 transformer couples a mass to an anchored spring | reflected stiffness `ratio² × stiffness`, the midpoint oscillation phase and conserved energy |
| `quarter_car.json` | a car body and wheel on a suspension spring and shock absorber, a tyre spring driven over a 30 mm bump; 70 % of the shock's heat warms it, the rest goes to the air | the ledger: road work = mechanical energy + heat, and the 70/30 split |

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

## Ideal transformers

`geared_oscillator.json` connects the mass's effort-input port to the spring's
flow-input port. A ratio `n` enforces `flow_b = -n × flow_a` and
`effort_a = n × effort_b`, conserving interface power. It reflects the
spring's stiffness by `n²`. The initial spring extension is `-n` times the
mass position, and reversing the signed ratio reverses that extension.
This ideal transformer has no losses, backlash or contact dynamics.
