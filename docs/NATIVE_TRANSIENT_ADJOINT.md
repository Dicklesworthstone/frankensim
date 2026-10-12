# Native transient temperature adjoints

A native `.fsim` project with a `conduction.transient` declaration can request
one of these additional report outputs alongside its existing scalar maximum:

```lisp
(outputs
  (qoi :name "temperature-max" :kind "scalar")
  (qoi :name "temperature-final-adjoint" :kind "report"))
```

Replace `temperature-final-adjoint` with `temperature-sampled-peak-adjoint`
to differentiate the earliest maximum over the accepted fine-grid endpoints
and the initial state. Both use the existing maximum requirement's volume
region. Do not attach an output `:region` (that field selects surfaces), and
do not request both derivative reports in the same project.

The ordinary `frankensim solve` path executes the request. The derivative
subreport is retained at `transient.adjoint` inside the conduction stage's
content-addressed JSON receipt. The HTML engineering report does not yet have
a dedicated transient derivative table; its conduction receipt hash remains
the link to these exact source data. A sampled-peak derivative does not change
the separate final-time scalar requirement or promote its uncertainty budget.

## Heat-capacity storage

Each reverse endpoint uses `C/dt + J(T_new)`, including the full smooth
nonlinear conductivity tangent and fixed matching contact. The old physical
field is never overwritten, including when the first prescribed-temperature
endpoint differs from the initial temperature. With prescribed exterior laws,
a reverse reconstruction must reproduce its retained forward temperatures bit
for bit before its derivative can be used. Natural-convection and
ambient-radiation endpoints instead bind the saved field directly under the
original forward physical residual target and energy tolerance. Their complete
transpose includes the changing convection/radiation coefficients and weighted
reservoir references. Both paths retain
their iteration, material-range and cancellation checks.

The subreport contains nodal and uniform initial-temperature derivatives,
a derivative with respect to a common multiplier on all regional heat
capacities, and reference-temperature/logarithmic-coefficient derivatives for
fixed Robin boundaries. It also contains absolute delivered-watt derivatives
for every volume region and every accepted fine interval. These controls use
the forward regional P1 projector, including shared-node mixing, and remain
available when the nominal regional power is zero.

For a declared constant-power workload segment, sum the fine-interval watt
derivatives whose intervals lie in that segment. The `regional_power` total
instead corresponds to adding the same watt increment throughout the entire
history. `dtemperature_dsource_multiplier_k` scales the interval's existing
volumetric source, not surface power. Duty factors are not reapplied to
absolute delivered-watt controls.

Controls after an earlier selected peak return zero. An initial-state peak
requires no reverse solves. Ties are reported through `active_ties`; the
selected earliest-state/smallest-vertex branch is deterministic but does not
become a unique derivative of a nonsmooth maximum. This is not a claim about
peaks between time samples.

### Ambient radiation

Driver 75 admits sourced gray ambient radiation with heat-capacity storage.
The accepted field, source history, prescribed-temperature reactions and
forward energy totals remain unchanged. The reverse pass reassembles the
actual radiation law at each consumed endpoint and checks the physical
tolerance retained from that endpoint's forward solve. It performs no
replacement primal solve.

The complete Jacobian includes the radiation secant's dependence on the
area-mean wall temperature and the combined convection/radiation reference.
The consistent pointwise Robin trace is retained, including nonuniform wall
temperatures. The same residual-checked multiplier drives initial temperature,
regional power, capacity, original convection and radiation controls; the
existing capacity-coordinate pullbacks do not receive another timestep factor.

`transient.adjoint.radiation.surfaces` preserves the declared surface order,
name, target, card, material state, emissivity usage receipt and query
temperature. Each row supplies:

| Control | Derivative field |
| --- | --- |
| Absolute hemispherical emissivity | `dtemperature_demissivity_k` |
| Black-reservoir temperature | `dtemperature_dreservoir_temperature_k_per_k` |

Each surface parameter remains constant throughout the selected prefix and
is independent even when two surfaces share a card. Emissivity directions
must remain inside (0, 1]; at one, only nonpositive increments are feasible.
The `robin` rows refer to original prescribed convective coefficients and
references; natural laws use the physical reservoir controls below. Radiation is never counted as heat delivered to an air stream.

`revalidated_solid_endpoints` counts consumed radiation endpoints.
`reconstructed_solid_endpoints` and reconstruction Krylov work are zero on
this path. The admitted workspace includes at most four nodal factor vectors
per Robin port, bounded parameter rows and one active endpoint workspace,
with at most 64 ports. The full physical targets are retained with the
temperature tape. An initial sampled peak has zero future boundary influence.

`data/reference-project/cooling-radiative-adjoint.fsim` requests a final
temperature adjoint for the existing synthetic startup-pulse model.
`tests/solve/transient_radiation.rs` compares all seven physical controls
against independent implicit scalar recurrences and fresh native project/card
perturbations for hotter and cooler radiation reservoirs.

### Natural convection with optional radiation

Driver 76 also admits the existing Churchill–Chu vertical-plate natural
convection law for heat-capacity trajectory adjoints. Heated and cooled walls
use the signed wall departure and film-density feedback from the same
correlation as the forward solver. Natural surfaces can radiate to separate
reservoirs in the same simulation.

The reverse pass evaluates each natural coefficient at the accepted Robin
wall mean. If radiation is present, both laws contribute to the combined
coefficient and weighted-reference derivative. It retains the **outer natural
solver's** physical acceptance target; the target from an inner frozen-natural
radiation trial is not substituted. The accepted temperature and energy
history remain unchanged.

`transient.adjoint.convection.natural_reservoirs` gives each natural
boundary's `target`, nominal `ambient_temperature_k`, and
`dtemperature_dambient_temperature_k_per_k`. This derivative includes both
the explicit reservoir reference and the natural coefficient's ambient and
film-density dependence. Pressure, geometry, characteristic length, transport
properties and correlation data remain fixed. Each natural reservoir is an
independent control, even when several surfaces use the same card.

Natural-law ports are omitted from the ordinary fixed-`robin` control rows.
Their effect on initial temperature, regional power, capacity and radiative
parameters still enters through the complete coupled transpose. An initial
sampled peak returns zero future reservoir influence without solving a
transpose. Equilibrium and unsupported correlation or validity-domain
derivatives retain their explicit refusals.

`data/reference-project/cooling-natural-radiative-adjoint.fsim` combines
synthetic heat storage, natural cooling and a distinct gray-radiation
reservoir. The native regressions in `tests/solve/transient_natural.rs`
compare hot and cold histories against independent physical recurrences and
fresh project/card perturbations. Distinct natural and radiative reservoirs
ensure the tests exercise both weighted-reference feedback terms.

## Total-enthalpy storage and coupled cooling

Projects declaring `transient.enthalpy` use the total-enthalpy reverse path.
It differentiates the accepted endpoint equations in specific-enthalpy
coordinates, including smooth temperature-dependent conductivity, declared
phase-conductivity feedback, fixed matching contact, and latent-plateau
interiors. Initial enthalpy controls hold the declared chart fixed; use
`:initial-temperature from-regional-enthalpy` when perturbing those controls
without maintaining a separate uniform-temperature assertion.

This path supports all of the following together:

- Prescribed Neumann and uniform Robin laws.
- In-domain Churchill–Chu vertical-plate natural convection on heated or
  cooled walls, including the signed wall difference and film-density feedback.
- Quasi-steady forced-air branches with each segment's dependence on every
  upstream wall on its branch.
- Sourced area-mean ambient radiation, including its complete temperature
  feedback and independent surface emissivity/reservoir controls.

Natural and forced-air laws must occupy distinct traces. Radiation can share
those traces while retaining its own reservoir heat account. Each reverse
endpoint rebuilds the actual laws at the saved field and checks its original
forward residual target and joule energy tolerance. The air check substitutes
actual stream enthalpy gain and keeps natural and radiation heat separate.
No accepted endpoint is re-solved or replaced.

The `fs-cli-native-enthalpy-adjoint-v1` report contains nodal and regional
initial specific-enthalpy derivatives, regional reference-density controls,
and absolute delivered-watt derivatives per fine interval and throughout the
trajectory. Radiation surfaces also carry absolute emissivity and reservoir-
temperature derivatives with their source identities. Every control uses the
same complete coupled transpose. Initial h carry remains energetic across a
latent plateau even when the instantaneous temperature slope is zero.
Driver 72 also returns natural-reservoir temperatures, independently supplied
air inlets, and absolute fan-bank speed ratios. These parameters remain
constant through the selected trajectory; they are accumulated using the
same complete multiplier at every consumed endpoint. These cooling coordinates hold constitutive parameters,
contact, geometry, pressure, the hydraulic density envelope and transport
properties fixed.

`transient.adjoint.convection` records the physical coupling method, natural
and air port counts, original-target rechecks, and the additional workspace.
At most 64 convection ports are admitted. Convection adds at most four nodal
factor vectors per port; radiation adds two per patch. Their factors operate
on the sparse enthalpy Jacobian without allocating a dense nodal matrix.
The fine-grid selection and tie rules above apply unchanged.

### Physical cooling controls

The following arrays appear inside `transient.adjoint.convection`:

| Array | Parameter identity | Derivative field |
| --- | --- | --- |
| `natural_reservoirs` | Boundary `target` and `ambient_temperature_k` | `dtemperature_dambient_temperature_k_per_k` |
| `air_inlets` | Air `branch` and `inlet_temperature_k` | `dtemperature_dinlet_temperature_k_per_k` |
| `fan_banks` | Fan `bank` and `speed_ratio` | `dtemperature_dspeed_ratio_k` when available |

A natural ambient derivative includes both its explicit Robin reference and
the coefficient's signed ambient/film-density partial. An inlet derivative
changes the branch supply temperature while retaining the hydraulic density
envelope, pressure and frozen transport properties. When editing a native
input for comparison, update every airflow declaration for that branch to the
same inlet temperature; do not change the envelope with it.

Fan derivatives include the hydraulic operating-point response, heat-capacity
rate, the convection card's Reynolds dependence, and the resulting downstream
air references. They differentiate one bank's absolute `speed_ratio` while
holding other banks' ratios fixed. They are not logarithmic-speed derivatives
or a common multiplier on all banks. The same retained hydraulic equations
supply series/parallel bank derivatives without another operating solve.

Each fan row reports `status: "available"` with its derivative, or
`status: "unsupported"` with an explicit `reason`. A missing physical
differential never becomes a partial derivative or an invented zero. A bank
that is unsupported at any consumed endpoint remains unsupported for the
trajectory. An initial-state selected peak has zero boundary influence,
because no future residual enters that observation.

`control_steps` equals the number of consumed endpoints in the selected
prefix. The aggregate arrays and bounded hydraulic scratch are charged once
in `control_workspace_bytes`; no control tape grows with the number of fine
steps. At most 64 fan banks are admitted, in addition to the convection-port
limit above.

### Material and contact controls

Driver 74 adds `transient.adjoint.material_controls` for the same selected
temperature observation. These controls remain constant through the consumed
trajectory and use the complete, residual-checked physical multiplier at each
endpoint.

| Array | Independent coordinate | Derivative field |
| --- | --- | --- |
| `conductivity_regions` | One scale on a volume region's entire conductivity tensor/curve, evaluated at one | `dtemperature_dconductivity_multiplier_k` |
| `contact_interfaces` | One scale on every face resistance of a named matching interface, evaluated at one | `dtemperature_dresistance_multiplier_k` |

A conductivity scale multiplies the existing `K(T)` and its declared phase
multiplier together. The retained element partial uses the actual accepted
temperature and mean nodal liquid mass fraction. Temperature and phase
feedback already belong to the complete state Jacobian; the parameter
contraction does not apply them a second time. Regional sums follow the
retained element labels and preserve anisotropic tensors and nonlinear curves.

A contact scale preserves each interface's resistance map and the independent
temperatures on its two traces. It uses the consistent P1 face mass matrix,
including nonuniform temperature and multiplier jumps. Separate interfaces
remain independent even when they share a material card. Contact geometry,
pairing, pressure, finish, and constitutive selection remain fixed.

The rows retain their material/interface identities. Derivatives are in
kelvin per dimensionless scale at one, which also equals the logarithmic-scale
derivative at zero. They are not absolute conductivity/resistance derivatives,
edits to material cards, or propagated material uncertainty. An initial
sampled peak returns zero for every material and contact control. Each array
admits at most 64 declared coordinates.

Only one set of four conductive-load entries per element, one element result
vector, a contact multiplier vector, and bounded aggregate rows are needed.
Their additional workspace is admitted before forward integration and does
not grow with the number of accepted steps. No additional state or transpose
solve is introduced.

## Scope and budgets

The capacity-storage path supports prescribed or in-domain Churchill–Chu
natural convection with optional ambient radiation, constant declared heat
capacities, static or scheduled volumetric sources, smooth heterogeneous
conductivity and fixed matching contact. Coupled airflow still requires the
total-enthalpy path described above. Prescribed temperatures are held fixed,
not returned as controls.
The capacity-storage path does not return conductivity or contact controls.
Neither path differentiates geometry, surface power, time-grid changes,
constitutive knots or latent-heat parameters.

The forward coarse/fine step cap remains a physical integration cap. The
optional derivative adds at most one reconstruction and one transpose per
fine endpoint up to the selected state, each with the existing linear
iteration limit. Their separate iteration totals and allowances are retained.
All reverse work shares the original wall deadline and cancellation context.
The additional fine-state tape, source bases, scratch and output estimate is
added to the forward memory admission before either grid executes; this is
not an allocator/RSS certificate. Without the output request there is no tape
or reverse solve, and forward receipt bytes remain unchanged.

Derivatives are Estimated discrete sensitivities, not continuum-error,
uncertainty, compliance or experimental-validation certificates.

## Focused verification

```sh
cargo test -p fs-conduction --test backward_euler_prescribed_adjoint
cargo test -p fs-conduction --test ambient_radiation accepted_radiative_endpoint_adjoint
cargo test -p fs-cli --test solve transient_radiation
cargo test -p fs-cli --test solve transient_natural
cargo test -p fs-cli --lib native_transient_adjoint
cargo test -p fs-conduction --test enthalpy_adjoint robin_feedback
cargo test -p fs-conduction --test enthalpy_radiation robin_feedback
cargo test -p fs-cli --test solve transient_enthalpy_convection_adjoint
cargo test -p fs-cli --test solve transient_enthalpy_material_adjoint
```

The core regressions compare linear/nonlinear prescribed launches and their
history/source/capacity derivatives against forward solves and finite
differences. The native-producer regressions parse real project declarations,
then exercise a manufactured mesh at the existing lowering seam: exact
insulated heating, zero-source controls, a pulse whose peak precedes the
final state, nonlinear conductivity, capacity and Robin derivatives,
initial-state peaks, repeated output, default-path equivalence, and invalid
or incomplete history refusal. They are not a full imported-geometry CLI test.
The coupled enthalpy cases add independent heated/cooled natural-reservoir
oracles, fixed-envelope inlet differences, absolute fan-speed differences at
a nonunit ratio, latent-history inlet propagation, initial-peak zero controls,
explicit unsupported bank reasons, and a spatial two-segment air/natural/
radiation control direction. The material-control regressions perturb fresh conductivity and contact cards,
including nonlinear conductivity through a latent-to-sensible history and
unequal-temperature bodies. Native execution of these new regressions remains
pending; no passing Rust result is asserted here.
