# Local material transport in evolving thermal flow

`MaterialSlab::step_with_material_transport` evaluates viscosity and thermal
conductivity at each fluid cell's actual accepted Kelvin temperature, then uses
those local properties in the existing 2D or 3D heat/flow step. Heating changes
subsequent material evaluations without resetting temperature or momentum.

The material callback returns any existing `FluidTransport` implementation. The
adapter extracts one snapshot from that value and retains the original concrete
object in `MaterialTransportStep::applied`. Source-backed objects keep their own
receipts and conditions; raw synthetic coefficients do not acquire measurement
authority. The same callback evaluates the candidate endpoint before publication,
with those results retained separately in `endpoint`. An unsupported endpoint,
failed conversion, invalid source, failed heat balance or observed cancellation
leaves the accepted slab unchanged. Externally captured callback state remains
caller-owned, including receipts from rejected attempts.

## Ownership and physical coordinates

The callback receives `(linear_cell_index, temperature_kelvin, retained_transport)`.
It may capture material-card assignments, pressure, composition, source-selection
policy and specimen conditions. Neither an instrument name nor a material label
selects a hidden model. The callback must consistently represent its stated law.
Only fluid cells are queried, in increasing index order; the two fixed-temperature
wall planes are separate reservoirs and do not require fluid material properties.

Power is supplied in W/m3 and external force density in N/m3. The existing
material-bound heating adapter owns their conversion and independently measured
source/storage/wall balance in joules. The 2D method requires an explicit extrusion
depth; the 3D method uses the retained cubic cell volume. Signed bottom-minus-top
temperature coordinates are preserved, so positive power still heats when the
top boundary is warmer. Buoyancy remains the existing Boussinesq force and is
composed with the supplied external loads.

`plan_scaling_from_material` remains the sole SI transport-to-lattice conversion
owner. D2Q5 and D3Q7 have different thermal quadrature coefficients; this adapter
uses the matching conversion, not the same numerical thermal relaxation in both
dimensions. The original collision/streaming kernels consume the resulting local
rates. No replacement current, flow, thermal or material-database implementation
is introduced.

## Storage and time integration boundaries

Density and specific heat must equal the retained material snapshot bit-for-bit.
This path evolves sensible heat with uniform, fixed volumetric storage. Changing
those quantities requires conservative mass/enthalpy evolution; it is not an
in-place coefficient edit. Viscosity and conductivity may vary in space and time.
Expansion/buoyancy parameters, geometry and wall-temperature coordinates remain
unchanged. There is no latent heat, melting, compressible work or mixture advection.

Coefficients are evaluated at the initial accepted temperature and held over one
step. Endpoint evaluation checks the material and lattice conversion at the
resulting state; it is not a nonlinear fixed-point solve or a proof of continuous
source support between the query points. Refinement in time remains necessary
for strongly temperature-dependent properties. The underlying BGK, low-Mach,
fixed-topology and advection-accuracy limitations still apply; a successful step
is not a stability, general transport-accuracy or experimental certificate.

Local coefficients are borrowed for one step. The background scalar rates are
restored without altering the evolved populations, so a later legacy call cannot
accidentally retain only the last local viscosity while reverting heat diffusion.
Call this method on each step to continue local feedback. Scalar `diffusivity`,
`rayleigh` and `nusselt` diagnostics still describe the background setup, not the
local coefficient field. The returned samples carry the actual applied rates.

Existing electrical/magnetic convenience methods retain their current behavior;
this addition does not silently enable local thermal-property feedback in those
methods. A caller can supply heat and force fields produced by those owners, but
owns their coupling schedule and work accounting. There is no new electrical,
magnetic, viscous-dissipation or chemical heat-source law here.

Cancellation is polled between material queries and around whole thermal steps,
not inside every lattice update. The callback must itself return. Original,
trial and thermal-owner tentative grids coexist, plus material samples and rate
arrays; `max_cells` bounds grid size, not total bytes. No throughput claim is made.

## Data-only runnable example

The example-only table format has a fixed `density specific_heat` first row,
followed by increasing `temperature conductivity dynamic_viscosity` rows in SI.
It explicitly uses piecewise-linear interpolation, rejects extrapolation and
limits the input to 64 KiB and 256 temperature knots. These are synthetic tables,
not sourced properties of a named substance or a new production material format.

```sh
cargo run -p fs-lbm --example material_transport_feedback -- \
  crates/fs-lbm/examples/material_transport_feedback/constant.txt 200
cargo run -p fs-lbm --example material_transport_feedback -- \
  crates/fs-lbm/examples/material_transport_feedback/temperature_dependent.txt 200
cargo test -p fs-lbm --test local_thermal_transport --test material_transport_feedback
cargo test -p fs-lbm --example material_transport_feedback
```

Both runs use the same force and heat schedule: heating is switched off halfway
through without resetting the simulation. CSV records temperature, velocity,
local conductivity range and physical heat balance. The two tables differ only
in transport data; their density and specific heat agree.

Four core and six material integration tests plus one example parser test are
authored. They cover original-path equivalence, spatial and temperature-dependent
conductivity against independent analytic profiles, real momentum collision,
concrete source retention, depth and reversed coordinates, endpoint failures,
cancellation and retry. Native Rust compilation/tests/example execution were not
available in the implementation environment (Cargo/rustc absent). Independent
NumPy stencil checks show approximately fourfold spatial-error reduction when
doubling resolution in the quiescent variable-conductivity cases. Separate
current-free flow/heat equations ran both 200-step example schedules; these are
reference equations, not execution of the Rust implementation or empirical data.
