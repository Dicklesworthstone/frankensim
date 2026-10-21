# Source-backed molecular gas mixing

This addition closes a different physical connection from nozzle flow: different
compositions can mix at equal total pressure. It does not replace the existing
pressure-driven network or alter its update order.

## Existing owners, connected operations

`ClosedGasParcel::exchange_species_with` accepts signed, integrated masses for
individual original `SpeciesId` values. Positive means left to right. Different
species may travel in opposite directions during the same interval. All fluxes
use the original states, and each volume is reconstructed only once.

Every species carries its own upstream ideal-gas partial specific enthalpy,
including formation reference and flow work. Using the mixture's enthalpy for a
selective flux would give the wrong energy. Source consistency, canonical
species/transport/pack joins, NASA-9 caloric inversion, EOS and the common medium
remain with the existing implementations. The complete new states expose their
updated density, heat capacity, viscosity, conductivity, sound speed and impedance.

`ClosedGasParcel::diffusion_flux_to` supplies an isothermal/isobaric ideal-gas
Maxwell-Stefan closure. At the central face it solves

```
sum(j != i) (x_j J_i - x_i J_j) / D_ij = -c grad(x_i)
sum(i) J_i = 0
c = p / (R T)
```

`J_i` is mol/(m2 s), relative to the **molar-average frame**. Zero total molar
flux does NOT imply zero total mass flux when species molar masses differ.
All species gradients and pair diffusivities participate in the coupled solve;
a species with zero own gradient can still diffuse in a multicomponent mixture.
The binary limit is the ordinary equimolar Fick law.

`diffuse_with` integrates these original-state rates over the supplied duration,
calls `exchange_species_with`, and checks the actual final pair entropy. Both
states are returned together or neither; no accepted input is changed.

## Explicit material data and geometry

Each unordered active species pair needs one `BinaryDiffusionMaterial` binding:
the existing `MaterialCard`, original `QueryPoint` context, explicit source
selection and two existing species identities. The canonical property is
`binary_diffusivity`, dimensions `[2,0,-1,0,0,0]`, in m2/s. Actual face temperature
and pressure replace those query coordinates; other source conditions survive.
Selected source records and evaluated coefficients are retained in the result.

No diffusivity is inferred from viscosity, conductivity, a species name or an
air preset. Missing/duplicate/foreign pairs and unsupported sources refuse.
Pair input order and the order of the two identities do not choose the answer.
Area [m2] and distance [m] come from `GasDiffusionGeometry`; zero area closes the
path. The local dense solve is bounded to 32 species and uses existing `fs-la`
LU, with cancellation checkpoints before and after that indivisible operation.
Every original unscaled equation is checked, including the equation replaced
by the molar-frame constraint.

The endpoint exchange checks gross withdrawal against each ORIGINAL species
inventory, not just total gas mass; incoming material cannot replenish that
allowance. Individual species, total mass and internal energy are checked both
locally and across the pair. No flux, trace species or temperature is clipped.
Prescribed selective transport alone does not prove a passive membrane or account
for separation/pumping work; the diffusion adapter is one explicitly scoped law.

## Run the synthetic example

From the repository root:

```sh
A=crates/fs-material/examples/material_gas_spring/mixture-a.tsv
B=crates/fs-material/examples/material_gas_spring/mixture-b.tsv
cargo run --release -p fs-material --example material_gas_diffusion -- \
  "$A" "$B" 0.00002 200 0.05 > diffusion.csv
cargo test -p fs-material --lib thermochemical_gas::energy::parcel::transfer::species::
```

The explicit diffusivity argument is synthetic constant material data, not a
measured coefficient for A/B. The original synthetic gas loader is reused.
Both volumes start at 400 K, 100000 Pa and 0.001 m3. The example uses 0.01 m2
area and 0.01 m separation. CSV reports composition relaxation, gas observations,
entropy change and energy balance. Changing input gas data or D requires no
solver change. The same resulting parcels can enter existing nozzle, thermal,
reaction and mechanical paths; this adapter does not automatically schedule those
operators or retune a live acoustic renderer.

## Model and verification limits

This is a two-point spatial approximation and an explicit first-order time step.
It neglects pressure diffusion, thermal diffusion/Soret and Dufour coupling,
external species forces, resolved bulk momentum, Knudsen effects, nonideal
mixtures, membranes and reactions DURING the diffusion substep. Caller-declared
T/p difference limits admit the local approximation; permissive limits do not
supply an error bound. Refine the interval and geometry separately. Positive
pair coefficients and a checked endpoint entropy change are not a proof of
unconditional finite-step stability or spatial accuracy.

Twelve authored pair tests cover partial enthalpy, whole-mixture parity,
new species, canonical ordering, Fick's limit, ternary cross diffusion, data-only
changes, repeated relaxation, original mechanical consumers and refusal/cancellation.
Eight network tests are described below: twenty library tests in total.
Native Rust compilation/tests were not run: Cargo and rustc were unavailable.
Independent Python checks executed 512 selective-transfer cases and 512
Maxwell-Stefan systems; the latter were compared with an independent symmetric
friction/velocity formulation. Binary relaxation gave refinement ratios
2.0144 and 2.0071. These are equation checks, not native or measured-gas results.

Formulation references (no foreign runtime dependency):

- COMSOL 6.3, *Multicomponent Gas Diffusion: Maxwell-Stefan Description*,
  equations 3-75--3-78:
  https://doc.comsol.com/6.3/doc/com.comsol.help.chem/chem_ug_chemsptrans.08.239.html
- Cantera, *Governing Equations for One-dimensional Flow*, species enthalpy
  transport and diffusive fluxes:
  https://cantera.org/stable/reference/onedim/governing-equations.html

These references use their stated velocity-frame conventions. The implementation
here explicitly selects the molar-average frame, not Cantera's mass-average frame.

## Shared-state diffusion networks

`transfer::species::network::exchange_candidate` accepts a bounded graph of
nonnegative individual-species masses with explicit donor and receiver IDs.
All outgoing material comes from each original inventory. Parallel records are
summed in canonical physical order, and every receiver is rebuilt once using
the existing source join and energy inverse. Gross per-species withdrawal is
checked before reconstruction; incoming material cannot conceal excessive
withdrawal. Idle nodes retain their exact accepted states.

`transfer::species::diffusion::network::step_candidate` connects those inventories
to the same Maxwell-Stefan face calculation. `DiffusionNetworkLink` supplies
physical endpoints, area/distance, actual pair cards, and `DiffusionFaceSolve`.
One endpoint-control map owns the energy inversions; `DiffusionNetworkSolve`
owns aggregate withdrawal, extensive balance and entropy limits. No separate
per-face endpoint reconstruction can spend the same node's inventory twice.
Face IDs determine report order only, not the physical result.

All faces see the same accepted composition, temperature and pressure. Their
individual species masses and partial enthalpies are combined before any gas
state is published. Newly received species cannot leave a second face within
this explicit interval. The whole-network entropy check comes from the actual
endpoints. The result retains source-bearing face samples and all original
local/global species, mass and energy accounts. All original isothermal/isobaric
and first-order limitations remain; pressure-driven flow, chemistry, thermal
feedback and wall motion are not automatically scheduled by this new adapter.

Eight authored network tests cover pair equivalence, input-order and report-ID
independence, symmetric fanout, no same-step downstream forwarding, aggregate
depletion despite incoming mass, a three-node diffusion eigenmode, closed/no-op
behavior and late failure/cancellation. They are unexecuted here. Independent
Python references passed 256 branching/cyclic network cases and three chain
refinements. Maximum energy discrepancy was 3.97e-13 J; maximum species-mass
discrepancy was 2.29e-19 kg; chain error-reduction ratios were 2.00778 and 2.00387.
These references do not run the Rust implementation.
