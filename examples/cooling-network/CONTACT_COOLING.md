# Finite-resistance contact in coupled cooling

The experimental `cooling-network` workflow accepts explicitly paired thermal
contacts between separately owned tetrahedral solids. It runs the existing
matching-P1 contact operator, not a lumped connector or a perfect-contact
approximation. Native `.fsim` and ledger workflows are unchanged.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network examples/cooling-network/size-contact-slab.json
cargo test -p fs-conduction --lib adjoint::robin::contact_tests
cargo test -p fs-cli --bin frankensim network_command::contacts
cargo test -p fs-cli --test cooling_contact_controls
```

Remove the example's optional `design` object for a single baseline simulation.
Its geometry is the existing 50 x 100 x 100 mm slab split at x = 25 mm: the two
sides of the bond use separate vertex IDs at exactly matching coordinates.
**Do not weld those IDs.** They carry separate temperatures, including a finite
jump at the contact. The sample has R'' = 0.01 m² K/W over 0.01 m², hence a
lumped series resistance of 1 K/W. R'' is area-specific, not itself K/W.

## Request

Add `contacts` inside `solid`. Every field in each declaration is required:

```json
"contacts": [{
  "name": "bondline",
  "source": "Caller-declared resistance for this scenario",
  "side_a_material": "declared left solid",
  "side_b_material": "declared right solid",
  "resistance_m2_k_w": 0.01,
  "face_pairs": [
    {"side_a": [1,4,10], "side_b": [12,13,15]},
    {"side_a": [1,7,10], "side_b": [12,14,15]}
  ]
}]
```

Each side is an exterior triangle of the corresponding disconnected solid.
Coordinates must match exactly under the existing production contact rules,
with opposing outward normals. Triangle vertex order need not match; the
contact kernel establishes correspondence geometrically. A face cannot belong
to another contact or to an external cooling surface. Every coincident pair
must be bound; omission is not an insulated gap or an implicit perfect bond.
Zero, negative, unrepresentable and nonfinite resistance inputs refuse.

With `adiabatic_remainder=false`, every non-contact exterior face must have an
external cooling owner. Internally the contact faces remain untagged by external
boundary conditions because the contact operator owns them. Allowing those faces
does not silently insulate other undeclared faces.

Material labels and the source string are caller declarations. They are retained
in an inline `fs-matdb` interface card solely to use the existing contact API;
this creates no measured material authority, experimental observation or verified
license. Labels are not independently cross-checked against material chemistry.
Resistance is constant within one simulation. Its source uncertainty remains
unstated unless a separate probability law is explicitly supplied to UQ below.

## Coupled solves and design

Contact participates in every heterogeneous solid solve, and in the solid
matrix used by both forward and transposed Robin derivatives. Thus the existing
coupled inlet and effective-h gradients include the changed temperature field
and downstream heated-air feedback across the bond.

The same path is used during effective-h and fan-speed searches. Fan trials
continue to recompute hydraulics and correlation-derived convection as before;
they do not erase or perfect the contact. Component heating still uses nodal
P1 support: a duplicated contact vertex on side A does not inject directly into
side B's distinct node. Heat crosses between the solids through the declared
interface law.

The `contacts` result array records the declared resistance and source, area,
conductance, mean A-minus-B temperature jump, and signed A-to-B heat transfer.
Internal contact heat is not counted again as a source or as external air heat.
The ordinary whole-solid and solid/air energy checks remain in force. Reversing
A/B changes the reported signs, not the solved temperature field.

## Total contact-resistance gradients

A steady request with `objective.gradient=true` now also reports
`contact_sensitivities.rows`, sorted by contact name. Each row contains:

- `dobjective_dlog_resistance_k`: the temperature-objective derivative for a
  relative resistance change, `dT_objective / d ln(R'')`, in kelvin.
- `dobjective_dresistance_w_m2`: the derivative per absolute area-specific
  resistance, equal to the preceding derivative divided by R''. Its unit is
  K/(m² K/W), or W/m².

For a small relative perturbation `delta`, the local prediction is
`T_objective(R'' * exp(delta)) = T_objective(R'') + delta * derivative` to
first order. This is a local derivative, not a finite-change uncertainty bound
or a proof of global monotonicity. The existing peak-objective derivative
semantics, including active-vertex/tie handling, remain unchanged.

The calculation reuses the TOTAL coupled nodal-load adjoint, so air mixing and
upstream/downstream thermal feedback are included. For nonlinear conductivity,
the solid tangent still contains the existing K'(T) term. No extra perturbed
PDE solve is performed per interface. For contact block Kc and p=ln(R''),
`dKc/dp = -Kc`, so the contraction is `lambda^T Kc T`. The exact P1 triangle
mass entries retain nonuniform primal/adjoint jumps; multiplying face-average
jumps would be wrong. Geometry and the material/contact pairing remain fixed.

Without a steady adjoint, `contact_sensitivities` is null, not a vector of
invented zeros. Transient requests still do not provide a transient adjoint.
No contact-flux objective, shape, material-parameter, pressure-dependent
resistance or nonmatching-contact derivative is introduced here.

## Resistance uncertainty through actual solves

The UQ command accepts a named contact target in m² K/W:

```json
{
  "target": {"kind":"contact-resistance", "contact":"bondline"},
  "distribution": {"kind":"uniform", "lo":0.005, "hi":0.02}
}
```

Run the nonlinear transient pulse with the illustrative resistance law:

```bash
frankensim --json cooling-network-uq \
  examples/cooling-network/nonlinear-contact-pulse.json \
  examples/cooling-network/uq-contact-resistance.json \
  --checkpoint contact-study.uqcp
```

The output destination must be new. The target also works in steady UQ, where
`qoi.kind=steady-objective` is explicit or omitted. Transient UQ requires
`qoi.kind=transient-sampled-peak`, as described in `TRANSIENT_UQ.md`.

Every sample changes the actual `solid.contacts[*].resistance_m2_k_w` consumed
by the normal producer. It does not change convection, geometry, face ownership,
material labels, loads or schedule. One resistance is sampled per trajectory
and held through all of its steps and fixed repeated cycles; independent
per-step contact noise is not implied. Existing checkpoint/resume, exact sample
ordinal retry and sequential-compliance options apply unchanged.

Missing or duplicate contact names refuse. Zero/negative/nonfinite resistance
and unrepresentable reciprocal conductance refuse rather than clipping,
redrawing or dropping a sample. Uniform supports must be positive. Gaussian
samples outside the physical domain terminate the UQ execution. Dependence
among multiple random parameters must still be explicitly declared; a local
sensitivity is not used as a surrogate for their real cooling simulations.
The ordinary per-run `contacts[].uncertainty` remains null: the UQ result carries
the declared sampled law and outcome, not a fabricated material-card interval.

## Independent reference values

For the unpowered split-slab example at its baseline h = 80 W/(m² K), the
continuous series/NTU solution gives about 2.087463703 W across the bond, a
2.087463703 K interface jump and a first-wall mean of 325.530546290 K. The
otherwise identical welded model gives about 324.35144 K, which would pass a
325 K mean limit while the declared contact model does not. A separate NumPy
FEM/bisection reference reaches a passing 324.999994091 K with h approximately
135.187229859 W/(m² K). The baseline first-wall log-resistance derivative is
approximately 0.932982239 K. These are independent reference calculations,
not retained execution measurements from this Rust command.

Tests cover the continuous reference and its reversal, equal inlet temperatures,
heterogeneous/nonlinear component-powered hotspots, perturbed coupled-FEM
contact gradients, canonical per-contact identities, missing/nonmatching/
double-owned faces, real steady and transient UQ, and exact checkpoint replay.
Source availability does not assert that the Rust tests have been executed.
Results remain nominal discrete model estimates: no nonmatching/mortar contact,
pressure-dependent resistance, radiation, geometry evolution, continuum-error
bound, continuous-time maximum guarantee or physical uncertainty certification
is inferred. The example probability law is illustrative, not measured data.
