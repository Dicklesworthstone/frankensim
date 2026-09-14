# Finite-resistance contact in coupled cooling

The experimental `cooling-network` workflow accepts explicitly paired thermal
contacts between separately owned tetrahedral solids. It runs the existing
matching-P1 contact operator, not a lumped connector or a perfect-contact
approximation. Native `.fsim` and ledger workflows are unchanged.

```bash
cargo run -p fs-cli --bin frankensim -- --json cooling-network examples/cooling-network/size-contact-slab.json
cargo test -p fs-conduction --lib adjoint::robin::contact_tests
cargo test -p fs-cli --bin frankensim network_command::contacts::tests
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
Resistance is assumed constant; its uncertainty remains explicitly unstated.

## Coupled solves and design

Contact participates in every heterogeneous solid solve, and in the solid
matrix used by both forward and transposed Robin derivatives. Thus the existing
coupled inlet and effective-h gradients include the changed temperature field
and downstream heated-air feedback across the bond. Contact geometry and R''
are held fixed; no contact-resistance derivative is claimed.

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

## Independent reference values

For the unpowered split-slab example at its baseline h = 80 W/(m² K), the
continuous series/NTU solution gives about 2.087463703 W across the bond, a
2.087463703 K interface jump and a first-wall mean of 325.530546290 K. The
otherwise identical welded model gives about 324.35144 K, which would pass a
325 K mean limit while the declared contact model does not. A separate NumPy
FEM/bisection reference reaches a passing 324.999994091 K with h approximately
135.187229859 W/(m² K). These are independent reference calculations, not
retained execution measurements from this Rust command.

The added tests cover the continuous reference and its reversal, equal inlet
temperatures, heterogeneous component-powered hotspots, perturbed coupled-FEM
gradients, missing/nonmatching/double-owned faces and actual target sizing.
Compilation, formatting and Rust tests were not executed in the authoring
environment, which has no Rust/Cargo, DSR or RCH. Results remain nominal discrete
model estimates: no nonmatching/mortar contact, pressure-dependent resistance,
radiation, geometry evolution, continuum-error bound or physical uncertainty
certification is inferred.
