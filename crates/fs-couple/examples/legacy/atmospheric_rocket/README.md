# Source gas through atmospheric rigid flight

This example connects the actual source-bearing gas, variable-caloric nozzle,
finite-reservoir blowdown, variable-mass rigid dynamics, rotating-atmosphere
wind, drag and central-gravity owners. The valve closes at an explicit step
boundary, after which the identical gas reservoir remains aboard a coasting
vehicle. There is no independent thrust or propellant-consumption curve.

From the repository root:

```sh
cargo run -p fs-couple --example atmospheric_rocket -- \
  crates/fs-material/examples/material_gas_spring/mixture-a.tsv \
  1000 1000000 > atmospheric_rocket.csv
```

The optional trailing arguments are total steps, physical timestep in seconds,
and powered steps. The defaults are 200, 0.01 and 100: one second of propulsion
followed by one second with the valve closed. For a temporal-refinement check
of the same physical journey, double both step counts and halve the timestep:

```sh
cargo run -p fs-couple --example atmospheric_rocket -- \
  crates/fs-material/examples/material_gas_spring/mixture-a.tsv \
  1000 1000000 400 0.005 200 > atmospheric_rocket_refined.csv
```

## Explicit demonstration model

The input is the existing manufactured two-species NASA-9/transport fixture,
not measured rocket propellant. Its original source identities, temperature
and pressure domains govern the entire tank/nozzle path. The example starts
with a 2 L tank, a 0.1 kg dry body, a 1 mm² throat and an exit/throat area ratio
of two. Dry principal moments are explicit physically admissible inputs. The
gas-specific moments are derived from the same 2 L volume and a declared
0.2 m uniform cylindrical tank. The circular outlet lies on the body
principal x axis; its co-rotating exhaust moment is declared separately from
the tank inertia and effective thrust application point.

The rigid owner uses a fixed common COM, fixed principal axes and declared
diagonal angular exhaust flux. This reduced model does not infer internal
relative-flow or transverse outlet linear momentum, slosh, shifting COM or a
general outlet inertia tensor. Those effects need a more general open-system
flow/momentum model; adding another axial lever torque would count the existing
thrust moment twice.

The vehicle starts airborne at 100 m geometric altitude with 30 m/s outward
radial speed in addition to the local rigid atmosphere-rotation velocity. It
has a 0.02 rad initial tilt and 0.02 rad/s body spin. There is no launchpad,
hold-down, ground collision or automatic launch event.

Air uses the lower U.S. Standard Atmosphere 1976 profile through 80 km geometric
altitude. This demonstration chooses its spherical zero-altitude radius from
the profile's effective radius, sets central mu to reference g times radius
squared, declares a 24-hour rotation, and adds a 5 m/s inertial background wind.
These are explicit demonstration choices, not a measured Earth gravitational
parameter, an ellipsoid, terrain, weather or a rotating hydrostatic solution.
Dynamics runs in an inertial frame; fictitious rotating-frame forces are not
added. The bulk air velocity includes the declared rotation at each actual
sample point.

The supplied aerodynamic law is explicitly **Cd = 0.5**, with 0.003 m² reference
area and 0.2 m reference length. Its force acts 0.04 m behind the common COM.
This is a manufactured drag law, not calibrated rocket aerodynamics. A real
consumer supplies a repeatable coefficient law valid at the passed Mach,
Reynolds number, body-relative flow direction and vehicle configuration. The
API imposes no hidden coefficient table or extrapolation. Aerodynamic heating,
lift/stability derivatives, rarefied flow and structural flexibility need their
own models.

## What is actually coupled

Each outer interval predicts the environment at a tentative half step. The
accepted full step then starts from the original gas/vehicle state with those
midpoint loads. The prediction cannot consume accepted propellant. The gas
owner retains its own bounded midpoint blowdown solve, and the rigid owner
retains its bounded attitude substeps and original momentum algebra. This is a
second-order composition for smooth conditions, so changing environments and
shock regimes still require temporal refinement.

Atmosphere is evaluated at the actual drag point and nozzle ambient point.
The drag-point velocity includes body rotation and is measured relative to the
local moving air. Its force supplies the matching body torque. Gravity and
specific potential are evaluated at the vehicle COM. The same accepted nozzle
mass and axial impulse drive mechanics, with pressure thrust included once.
Scalar axial impulse remains distinct from the world-vector impulse as the
body turns.

The final atmosphere and coefficient-law domains are checked before either
physical owner is published. Unsupported source states, a work-budget failure
or cancellation return no partial endpoint. Nozzle failure never closes the
valve implicitly. The caller explicitly splits an interval at commanded valve
closure; this example does so using the powered-step count. The flight layer
does not locate ground impact or crossings of a source-domain boundary.

## Output

CSV rows identify endpoint time and midpoint load time separately. They report
propulsion mode and nozzle regime, altitude, local pressure, Mach, Reynolds
number, dynamic pressure, tank state, radial velocity, native quaternion and
angular speed, thrust, cumulative axial impulse and expelled mass. The final
columns retain actual mass, axial-impulse and angular-momentum residuals.
During closed-valve flight, engine thrust and mass flow are zero while gravity,
drag and attitude evolution continue.

The nozzle supports ideal subsonic, internal-normal-shock, attached external-
compression and matched/underexpanded branches. It does not predict viscous
separation, side loads, oblique plume structure or transient shock motion. Gas
composition is frozen in this example; the separate equilibrium chamber API
can prepare source-preserving frozen product parcels before a flight.
