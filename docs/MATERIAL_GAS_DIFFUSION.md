# Source-backed molecular gas mixing

Molecular diffusion transports species through composition gradients and, with
an explicit mass-frame option, pressure gradients. The original equal-pressure
molar counterdiffusion APIs remain unchanged. These molecular substeps do not
replace the bulk-flow network or evolve mixture momentum.

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

## Explicit pressure diffusion in the mass-average frame

`GasDiffusionDriving::CompositionAndPressureMass` adds the ideal-gas pressure
force through the existing shared Maxwell–Stefan solver:

```
d_i = grad(x_i) + (x_i - Y_i) grad(ln p)
Y_i = x_i M_i / sum(k) x_k M_k
sum(i) M_i J_i = 0
```

The mode explicitly chooses the **mass-average frame** and requires its own
positive finite `mass_flux_tolerance_kg_m2_s`. The existing
`molar_flux_tolerance_mol_m2_s` applies to the default `IsobaricMolar` mode.
A pressure-driven mass-frame flux generally has a nonzero molar sum. Reports
retain the driving/frame choice, both signed flux sums, the log-pressure
gradient, original-equation residual and the shared kernel's nonnegative
pair-friction entropy rate. The sum for the other frame is diagnostic, not a
failed conservation constraint.

Use `diffusion_flux_to_with_driving` for the unchanged-state face calculation
or `diffuse_with_driving` for the conservative two-volume substep:

```rust
let driving = GasDiffusionDriving::CompositionAndPressureMass {
    mass_flux_tolerance_kg_m2_s: 1e-16,
};
let trial = left.diffuse_with_driving(
    &right, geometry, duration_s, &pairs, controls, driving, || false,
)?;
```

The supplied `maximum_pressure_difference_pa` and temperature-difference limit
still bound the local approximation. Pressure mode qualifies the actual left,
face and right T/p against the same selected diffusivity claim; an admissible
mean cannot hide an unsupported endpoint. The logarithmic pressure difference
uses the existing deterministic, cancellation-resistant ratio evaluation, so
adjacent representable pressures do not lose their driving force.

The mass frame matters thermodynamically. For positive face fractions,

```
sigma_drag = -R sum(i) J_i d_i/x_i
           = -R sum(i) J_i grad(ln(x_i p))  when sum(i) M_i J_i = 0.
```

Thus pressure driving contributes to the full chemical entropy balance. Using
pressure force with a zero-molar-sum constraint and treating that isolated
substep as a closed system would omit a mechanical contribution. This API uses
the mass frame and retains every existing finite-step entropy rejection.
Species still carry their original upstream partial enthalpies once, including
flow work. The local friction rate is separate from the actual endpoint entropy
change and does not certify an arbitrary finite step.

At equal pressure the new mode is the mass-frame transformation of the default
molar counterdiffusion, not generally the same species flux. Equal species molar
masses remove both this frame difference and the pressure-separation term.
Different molar masses allow separation even from uniform mole fractions.

`network::step_candidate_with_driving` applies the same explicit choice to all
faces in one original snapshot. `network::adaptive::advance_candidate_with_driving`
resamples the pressure force and source endpoints at every coarse/fine probe.
All aggregate withdrawal, species, energy, endpoint entropy and complete-window
audit gates remain in force. Neither API supplies bulk momentum, thermal
diffusion, or an automatic schedule for those other physical operators.

Six focused regressions cover an independent binary pressure-diffusion oracle,
the full chemical entropy force, reversal, equal-pressure frame conversion,
equal-mass cancellation, adjacent pressures, source/cancellation refusal, and
conservative pair, fixed-network and adaptive-window execution. Run them with
the existing module suite:

```sh
cargo test -p fs-material --lib thermochemical_gas::energy::parcel::transfer::species::diffusion
```

The pressure-force formulation follows COMSOL 6.4, *Multicomponent Gas
Diffusion: Maxwell–Stefan Description*, equation 5-44:
https://doc.comsol.com/6.4/doc/com.comsol.help.fce/fce_ug_chemsptrans.08.168.html

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
path. The local dense solve is bounded to 32 species and uses the shared
`fs-thermochem` scaled friction owner, with cancellation checkpoints during
assembly and elimination. Every original unscaled equation is checked,
including the equation replaced by the selected frame constraint.

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
The default mode neglects pressure diffusion; the explicit mass-frame mode
includes it. Both modes omit thermal diffusion/Soret and Dufour coupling,
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

These references use their stated velocity-frame conventions. The original
entry points select the molar-average frame; the explicit pressure-diffusion
mode selects the mass-average frame.

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
