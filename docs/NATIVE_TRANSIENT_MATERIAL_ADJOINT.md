# Native transient material sensitivities

Capacity-storage projects can explicitly add regional conductivity and named
contact-resistance controls to their complete transient temperature adjoint:

```lisp
(outputs
  (qoi :name "temperature-max" :kind "scalar")
  (qoi :name "temperature-final-material-adjoint" :kind "report"))
```

Use `temperature-sampled-peak-material-adjoint` instead to select the earliest
maximum over the accepted fine endpoints and the initial state. Select exactly
one adjoint report, without an output region; the existing maximum requirement
supplies its volume region. The final-time scalar requirement remains separate
from a sampled-peak sensitivity.

These are explicit extensions of the reports described in
[NATIVE_TRANSIENT_ADJOINT.md](NATIVE_TRANSIENT_ADJOINT.md). Existing capacity
`temperature-final-adjoint` and `temperature-sampled-peak-adjoint` requests keep
their original output, work allowance and scope. The new output name is bound
into the canonical project and run identity. Total-enthalpy reports already
include material controls; the new names also select the corresponding
enthalpy observation without changing its storage coordinates.

## Physical coordinates

Results are under `transient.adjoint.material_controls` in the retained
conduction receipt. Both multipliers have reference value one and apply
throughout the selected trajectory:

| Array | Independent coordinate | Derivative field |
| --- | --- | --- |
| `conductivity_regions` | A region's entire conductivity tensor or temperature curve, scaled together | `dtemperature_dconductivity_multiplier_k` |
| `contact_interfaces` | Every retained face resistance of one matching interface, scaled together | `dtemperature_dresistance_multiplier_k` |

Each derivative is kelvin per dimensionless multiplier at one, also equal to
the logarithmic-multiplier derivative at zero. A positive resistance derivative
means a local increase of that interface resistance raises the selected
observation. Separate regions remain separate even when their material cards
match. Interface rows retain card identity, mapping status and face-pair count.

Conductivity partials evaluate the original production element stiffness at
the accepted endpoint temperature, including anisotropy and nonlinear K(T).
The full state tangent already contains K-prime and every supported cooling
feedback term. The fixed-state parameter contraction must not apply K-prime
again. Matching-contact partials retain both independent temperature traces,
resistance maps and the consistent P1 face integration.

The same residual-checked physical multiplier drives material, contact, power,
capacity, inlet, fan, natural-reservoir and radiative controls. There is no
extra primal or transpose solve per material. Capacity adjoints already expose
watt-coordinate nodal load weights: neither material nor contact receives an
extra timestep factor. The shared enthalpy path retains its separate once-only
conversion from joule coordinates.

Only the selected prefix contributes. An initial sampled peak returns zero
material/contact influence while retaining all parameter identities. Ties
retain the existing deterministic branch and do not claim a unique derivative.
The accepted forward field, energy history and other gradients are unchanged.

The extra element ownership/partial scratch, contact multiplier and aggregate
rows enter memory admission before either forward grid. At most 64 material
regions and 64 matching interfaces are admitted. This is a bounded logical
workspace estimate, not a measurement of allocator overhead or peak RSS.
Material corners, unsupported contact models, failed residuals, insufficient
work or cancellation retain their existing refusal behavior; no incomplete
thermal field or gradient is published.

## Example: transient heating across a contact

The input `examples/contact-pair/capacity-material-adjoint.fsim` reuses the
existing two-block geometry and fixture cards. Both bodies start at 300 K with
explicitly synthetic capacities, the hot body receives 100 W for two seconds,
and all exterior faces are insulated. Its independent total-input reference
is **200 J**. Contact heat is internal and must not add another energy source.
The example asks for final-temperature conductivity and interface controls.

Run from the repository root with a writable work directory:

```sh
WORK="$(mktemp -d)"
PROJECT=examples/contact-pair/capacity-material-adjoint.fsim
cargo run -p fs-cli --bin frankensim -- --json validate "$PROJECT"
cargo run -p fs-cli --bin frankensim -- --json import "$PROJECT" \
  examples/contact-pair/cold-body.stl examples/contact-pair/hot-body.stl \
  "$WORK/pair.db" --unit m --max-hole-edges 0
cargo run -p fs-cli --bin frankensim -- --json solve "$PROJECT" "$WORK/pair.db" \
  --materials data/reference-project/aa6061.fsmcdpk \
  --interfaces examples/contact-pair/cold-hot.fsintpk
```

## Verification and scope

Focused native tests are in
`crates/fs-cli/tests/solve/transient_conjugate/material_controls.rs`:

```sh
cargo test -p fs-cli --test solve native_capacity_material
cargo test -p fs-cli --lib material_output_tests
```

They compare fresh immutable conductivity/contact-card perturbations with the
reported derivatives, include nonlinear conductivity with two-segment air,
natural convection and radiation, preserve the old report's gradients and
transpose work, distinguish shared-card regions, and check initial peaks and
ambiguous report refusal. These new Rust tests and the example have not been
executed in the implementation environment, which lacks the Rust toolchain.
An independent eight-node P1 equation reference agreed with fresh trajectory
perturbations within 2.55e-9; this is not execution of the repository solver.

These are local Estimated fixed-grid model sensitivities, not material-card
mutations, absolute tensor-entry derivatives, uncertainty propagation or
physical validation. They do not differentiate geometry, constitutive knots,
latent-heat parameters, prescribed temperatures or time-grid changes, and do
not bound continuous-time peaks or establish a global design optimum.
