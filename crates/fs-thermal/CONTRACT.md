# fs-thermal contract

Layer L3. Shared conservative thermal transport and solid/vapor exchange.
Runtime dependencies: std, fs-la, fs-math. Spatial discretizations and material
constitutive data remain owned by their existing crates.

## Existing thermal timeline

The enthalpy, advection, contact, observation and event APIs retain their
documented control-volume semantics. A cell stores total energy as mass times
specific enthalpy. Every material chart in one balance uses a compatible energy
datum. The shared timeline is re-exported by fs-conduction; consumers must not
implement independent copies of its storage update.

## Sublimation and deposition

`sublimation::SublimationLaw` supplies a bounded immutable saturation-pressure,
latent-enthalpy, density and molecular-mass chart. The material crate supplies
reference ice/CO2 correlations and an explicitly parameterized
Clausius–Clapeyron approximation. A numerical correlation is not a sourced
specimen certificate, and supplied density/accommodation/film coefficients
retain their caller-declared applicability.

Positive mass flux means solid becomes vapor; negative flux deposits vapor.
The input pressure is the volatile partial pressure. Outgoing molecules use
surface temperature and incoming molecules use gas temperature. One common
accommodation coefficient preserves detailed balance at equal temperature.
The optional dilute gas film is combined with interface kinetics by series
resistance. Zero accommodation or film coefficient seals the surface.

`surface_flux` returns signed mass and latent-energy fluxes, their analytic
surface-temperature derivatives, and outward-normal interface velocity.
`isothermal_mass_transfer` integrates an ideal, thermostated finite vapor
reservoir exponentially at fixed surface temperature. It uses the shared
deterministic expm1 implementation to preserve short timesteps and does not
overshoot its frozen-temperature equilibrium. Its caller supplies the thermal
work needed to hold the surface at that temperature.

`advance_surface` couples that transfer to the caller's original T(h) chart
with a bounded implicit endpoint solve. The endpoint vapor carries h_end + L.
The energy equation is:

    m_initial * (h_end - h_initial) + transferred_mass * L(T_end) = 0

This is algebraically the balance of remaining solid energy and transported
vapor enthalpy. The latent account is part of vapor energy, not another
independent energy debit. All mass limits are applied before deriving heat,
vapor inventory and recession. A shifted common caloric datum leaves the
physical temperature and mass solution unchanged.

The nonlinear solver requires a continuous nondecreasing caloric chart and a
bracketable endpoint in the declared constitutive domain. It refuses
unsupported temperatures and failed convergence rather than clamping them.
The default budget permits 96 evaluations; callers may select 1 through 4096.
Each evaluation polls the supplied checkpoint. Functions return candidate
values and mutate no caller state on failure.

## Numerical and physical limits

The coupled endpoint is a first-order reduced model. Area remains frozen over
one phase step; recession uses endpoint solid density. Finite gas temperature
is prescribed: its thermostat owns subsequent sensible-heat equilibration.
The dilute film does not model non-dilute Stefan flow, pore resistance, plumes,
curvature, nucleation, melt competition, radiation or unresolved gradients.
Those mechanisms require explicit compatible models, not hidden coefficients.

Returned zero solid mass denotes extinction. A spatial consumer must update
topology, contacts and conservative remapping or refuse that endpoint.
Reporting recession does not itself move a mesh, alter a level set or change a
flow/electrical geometry map. Uniform-cell use requires a justified internal
temperature approximation. Refining time cannot repair an invalid spatial
model.

The scalar flux and finite-reservoir kernels allocate no memory and perform
constant work. The implicit solver uses constant memory and at most its
declared evaluation budget. This is an algorithmic statement, not a measured
hardware throughput or cross-ISA reproducibility claim. Transcendentals use
the existing fs-math owner; current executed evidence is same-ISA.

## Focused verification

G0/G1/G3 tests cover signed deposition, equilibrium, unequal gas temperature,
analytic tangents, kinetic/film limits, finite-reservoir exact relaxation and
semigroup behavior, tiny timesteps, inventory limits, implicit cooling and
heating, mass/energy/recession consistency, caloric-datum shifts, first-order
refinement, cancellation and exhausted nonlinear budgets. Material tests
compare the actual correlations to published reference values and their
analytic derivatives to finite differences. These tests establish the
declared reduced model; they are not experimental validation of a complete
sublimating object or a general spatial Stefan solution.


## Phase exchange on the existing heat timeline

`enthalpy::sublimation` composes the surface operator with the original
`DepositionSimulation`. `sublimate_cell` atomically changes the active cell's
actual mass and specific enthalpy plus its vapor inventory, without advancing
the clock. `advance_sublimating_to` performs bounded first-order heat/phase
splitting with one total thermal work budget and one global mass/energy audit.
It honors birth and boundary/contact times exactly, uses the exposed area from
the start of each interval, and never sublimates a not-yet-born cell. The
existing caloric chart remains authoritative; rejected out-of-domain trial
points do not invalidate a successfully converged in-domain endpoint. All
failures and cancellation preserve the original timeline and reservoir.

The whole cell must remain a pure solid on its original caloric branch.
Temperature alone cannot establish solid fraction on a melting plateau;
mixed-phase inventories require a separately owned partition. The supplied
finite reservoir belongs to this one surface and remains thermostated. A
joint chamber shared by several surfaces needs a joint backpressure solve.
Reference geometry/contact areas remain fixed and recession is reported.
Complete cell exhaustion refuses until a consumer owns topology removal.
There is no implied moving mesh, gas-energy solve or material calibration.

Vapor energy includes the solid's original sensible-enthalpy datum plus latent
heat; the separately reported latent account is never an additional debit.
Ten new native integration regressions cover balances, frost growth, covered
surfaces, birth/contact timing, nonlinear chart trials, transactional refusal,
whole-advance budgets and first-order refinement. The combined exact-source
thermal/material suite passed 95 tests with one existing performance test
ignored. `cargo run -p fs-material --example sublimating_ice` demonstrates the
real ice law with explicit illustrative calorimetry in a maintained vacuum and
a finite chamber; both native example cases passed mass/energy assertions.
