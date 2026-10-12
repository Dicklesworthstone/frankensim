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
endpoint differs from the initial temperature. A reverse reconstruction must
reproduce its retained forward temperatures bit for bit before its derivative
can be used; the reconstructed solve and transpose each keep their residual,
iteration and cancellation checks.

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
Natural reservoir temperatures, branch inlets, flow rates, constitutive
parameters, contact and geometry remain fixed in this report.

`transient.adjoint.convection` records the physical coupling method, natural
and air port counts, original-target rechecks, and the additional workspace.
At most 64 convection ports are admitted. Convection adds at most four nodal
factor vectors per port; radiation adds two per patch. Their factors operate
on the sparse enthalpy Jacobian without allocating a dense nodal matrix.
The fine-grid selection and tie rules above apply unchanged.

## Scope and budgets

The capacity-storage path supports fixed exterior laws, constant declared heat
capacities, static or scheduled volumetric sources, smooth heterogeneous
conductivity and fixed matching contact. Radiation, natural convection and
coupled airflow require the total-enthalpy path described above. Prescribed temperatures are held fixed, not returned as controls.
Geometry, conductivity parameters, contact resistance, surface power,
time-grid changes and latent heat are not differentiated by this path.

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
cargo test -p fs-cli --lib native_transient_adjoint
cargo test -p fs-cli --test solve transient_enthalpy_convection_adjoint
```

The core regressions compare linear/nonlinear prescribed launches and their
history/source/capacity derivatives against forward solves and finite
differences. The native-producer regressions parse real project declarations,
then exercise a manufactured mesh at the existing lowering seam: exact
insulated heating, zero-source controls, a pulse whose peak precedes the
final state, nonlinear conductivity, capacity and Robin derivatives,
initial-state peaks, repeated output, default-path equivalence, and invalid
or incomplete history refusal. They are not a full imported-geometry CLI test.
Native compilation and execution of these new regressions remain pending;
no passing Rust result is asserted here.
