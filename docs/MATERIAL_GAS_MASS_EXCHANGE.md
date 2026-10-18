# Material-driven gas mass exchange and pressure-driven flow

`ClosedGasParcel::transfer_mass_to` connects two existing material gas volumes
with a prescribed mass transfer. `through_nozzle_to` determines that mass from
an ideal converging restriction, the actual upstream material, downstream
pressure and explicit physical duration. Both operations return paired tentative
endpoints. Publish both together; never deposit the reported enthalpy again.

## One species inventory and energy model

A mass `dm` transports the donor's mass fractions and upstream specific
**enthalpy**, including flow work. The donor loses `dm*h_upstream` and the receiver
gains exactly that energy. Using internal energy instead would incorrectly leave
the withdrawing gas isothermal. Each receiving species mass is its old inventory
plus the incoming mass of that same canonical species. The existing quantity
owner converts the resulting mass fractions to mole fractions, and the existing
NASA-9/EOS/transport owners reconstruct the complete endpoint.

The actual donor, receiver and pair energy residuals are checked in joules.
Total and individual species masses are checked in kilograms. No tolerance is
normalized by a large formation energy. Upstream enthalpy may be negative under
a valid source reference: its sign is retained, not clamped. The donor must retain
positive represented mass, and transfers smaller than the supplied absolute mass
resolution refuse instead of disappearing. Zero transfer preserves both states.

Shared species must have identical original caloric models, transport cards and
selection policies; conflicting sources refuse rather than being chosen by
material name. The model also retains canonical pack receipts and the explicit
water-species label. Each volume keeps its own non-temperature/pressure query
axes. Unsupported incoming species at the receiver's local conditions refuse.
The constructor's 128-species bound remains in force. Composition conversion
corrects only a final floating-point normalization residue, then checks the
resulting individual species inventories; positive trace species cannot vanish.

## Pressure-driven restriction and choking

`GasNozzle { area_m2, discharge_coefficient }` declares flow area and a coefficient
in `[0,1]`. Zero area or coefficient shuts the path. Upstream and downstream are
explicit; swap the arguments for reverse pressure direction.

The flow solve uses the original material functions to satisfy

- `h_upstream = h_throat + velocity^2/2`,
- `s_throat = s_upstream`,
- `mass_rate = coefficient * area * density_throat * velocity`.

Back pressure sets the subsonic exit until velocity reaches the **same material
sound speed** used by acoustic consumers. Beyond that limit the flow is choked;
further lowering downstream pressure does not increase the mass rate. The
bounded temperature solve uses enthalpy and entropy from NASA-9, not a gas-name
lookup, an incompressible Bernoulli rate, or a gamma frozen at the upstream state.
For a small pressure difference only the actual subsonic interval must be
supported; there is no unnecessary demand for data at a hypothetical sonic state.

The final throat is fully material-resolved, retaining its transport/pack data
and independently checked entropy and sonic residuals. A physical step consumes
`mass_rate * duration` via the original paired transfer. A donor-fraction cap and
an endpoint pressure-order check refuse oversized physical steps. They do not
clip transferred mass or prove a timestep-error bound: shorten and refine time.

The receiving uniform tank thermalizes the jet. The transferred upstream total
enthalpy already accounts for its kinetic energy; adding that energy as a second
heat source is incorrect. These conventions follow the NASA Glenn isentropic
flow and choking derivations, and the ideal-mixed-volume balance convention:

- https://www.grc.nasa.gov/www/BGH/isndrv.html
- https://www.grc.nasa.gov/www/BGH/mflchk.html
- https://doc.modelica.org/om/Modelica.Fluid.Interfaces.PartialLumpedVolume.html

## Run the data-only two-volume example

From the repository root:

```sh
cargo run --release -p fs-material --example material_gas_exchange -- \
  crates/fs-material/examples/material_gas_spring/mixture-a.tsv \
  crates/fs-material/examples/material_gas_spring/mixture-b.tsv \
  200 0.0002 0.000001 > gas-exchange.csv

# Reverse only the two mixture files to change the material response.
cargo test -p fs-material --lib thermochemical_gas::energy::parcel::transfer::
cargo test -p fs-material --example material_gas_exchange
```

The example reuses the existing synthetic gas-spring input loader and its two
files; no new material parser or production schema is introduced. Both input
paths are recorded together as one declared **synthetic** source set, allowing
shared-species laws to match when only fractions differ. Numerical/source-law
conflicts still refuse. This paired example provenance is not measured-gas
qualification or a substitute for the production source-pack constructor.

The two rigid volumes are each `0.001 m3`, initially `500 K / 300000 Pa` and
`400 K / 100000 Pa`. The fixed discharge coefficient is explicitly `0.8`.
The optional arguments set steps, duration per step and area; no heat or motion
is imposed. CSV reports both inventories, temperatures, pressures and sound
speeds; receiver gas constant, impedance and transport; signed mass rate,
choking, accumulated transfer and actual pair mass/energy errors. Flow is driven
by the current pressure difference, and every accepted receiver retains the new
composition for its next step. A closed path preserves both initial states.

## Checks and limits

Six transfer, six nozzle and two new example tests exercise the actual public
owners. They cover enthalpy-versus-internal-energy falsification, new species,
source conflicts, no-op identity, active/choked source domains, analytic
constant-capacity references, variable-capacity choking, area scaling, time
refinement, original heat/work/force composition and cancellation/refusal.
The included shared loader also contributes its existing malformed-input test.

Native Rust compilation/tests could not run in the authoring environment:
`cargo` and `rustc` are absent. Independent Python equations checked 512 mass
transfers (maximum energy error `2.99e-13 J`, species error `4.07e-19 kg`) and
400 nozzle cases (maximum relative reference disagreement `9.03e-12`). Five
flow histories totaling 512 steps gave pressure refinement ratio `2.001` and
maximum step/history energy discrepancy `1.36e-12 J`. These are mathematical
references, **not** executed Rust tests or experimental material validation.

This is first-order quasi-steady exchange between rigid, quiescent, well-mixed
ideal-gas volumes, with no species reaction during transit. The nozzle assumes
an ordinary continuous source-supported gas branch with one sonic crossing.
There is no resolved valve friction, actuator, pipe residence delay, shock,
supersonic diverging nozzle, phase equilibrium, jet momentum field, combustion,
or spatial acoustic renderer. The discharge coefficient is caller data, not
calibrated here. Energy checks do not prove nonnegative entropy production.
Existing heat/piston operations may consume the returned states, but splitting
transfer and moving-boundary updates is not a new simultaneous open-system
integrator. Refine the full composed trajectory separately.
