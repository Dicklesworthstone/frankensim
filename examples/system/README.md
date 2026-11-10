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
| `quarter_car.json` | a car body and wheel on a suspension spring and shock absorber, a tyre spring driven over a 15 mm bump; 70 % of the shock's heat warms it, the rest goes to the air | the ledger: road work = mechanical energy + heat, and the 70/30 split |
