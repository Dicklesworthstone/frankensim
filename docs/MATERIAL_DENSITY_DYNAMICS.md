# Graded reference density in solid motion

The same total mass can have different centroids, rotational inertia, local
acceleration and contact mobility. For a P1 tetrahedron, density is interpolated
from four positive reference-vertex values and the consistent mass integrates
`rho_0 * N_a * N_b` exactly. Replacing these values by their element mean keeps
total mass but loses the distribution that controls those responses.

For example, on the unit tetrahedron, `rho_0 = 3 + 12*x` kg/m³ has mass 1 kg
and centroid `(0.325, 0.225, 0.225)` m. A uniform 6 kg/m³ tetrahedron has the
same mass but centroid `(0.25, 0.25, 0.25)` m. The consistent X mass at vertex
`(1,0,0)` is 0.15 kg for the graded specimen and 0.10 kg for the uniform one.
With that coordinate free and all others held, the same 0.3 N force therefore
produces initial acceleration 2 m/s² and 3 m/s², respectively.

## Source-backed finite motion

Resolve the original reference positions and density queries with
`fs_mbd_material::density_field::TetrahedralDensityBody::resolve`. Each vertex
retains its immutable material bundle, complete physical query and selected
density receipt. Supply that existing owner through an `Arc` to
`FiniteSolidBody::new_with_reference_density_field`, or use
`new_with_reference_density_field_and_pressure` for current-area pressure.

This is an explicit independent inertial-field choice. It replaces the scalar
inertial densities in the constitutive bindings; their elementwise elastic laws
remain authoritative for stress. No stiffness interpolation or equivalence
between different source states is inferred. Density sources must satisfy the
existing fixed-external-condition admission; an unowned evolving stress,
strain, rate or history axis cannot be silently frozen.

The reference vertices and tetrahedral partition must match the solid exactly.
The density owner's canonical cell order is used only to check this identity;
the original finite-element order continues to define material assignments,
pressure faces and contact-point indices. The full consistent mass is prepared
once and shared by time steps, point impulses, contact mobility and candidates.
`reference_density_field()` retains the original owner and receipts. Deformation
never changes the reference inventory; the existing density body's `deform`
method observes current matter with `rho = rho_0 / J`.

`reference_acceleration_load_n` projects a prescribed nodal acceleration through
the body's actual full mass. Pass its sorted nodal loads to the existing step
methods without adding the same gravity again as a force density. Supported
rows remain included. Position-dependent field resampling, potential and tangent
belong to the external field owner.

## Native solid APIs

| Consumer | Entry point | Retained behavior |
| --- | --- | --- |
| Finite hyperelasticity | `HyperTetProblem::reference_mass_matrix_with_nodal_density` | Exact reference inertia, including supported coordinates |
| Static acceleration histories | `TetLoadHistory3::at_with_nodal_reference_density` | Density-weighted force, torque and native reactions |
| Prescribed acceleration fields | `TetLoadGeometry3::assemble_acceleration_with_nodal_reference_density` | Exact `M*a` for uniform or P1 acceleration |
| Finite viscoelasticity | `FiniteViscoTetMotion::try_new_with_nodal_density` | The same mass plus complete accepted Maxwell memory |

The separate `loads3` transient constructors still use their existing material
mass. Supplying a graded static acceleration load does not replace that inertia.
The finite-motion and viscoelastic entries above do bind graded dynamic mass.

## Verification boundary

Focused G0/G1/G3 regressions cover analytical cubic integrals, centroids,
kinetic energy, scalar parity, shared vertices, gravity reactions, source-bound
motion and impulses, pressure initialization and viscoelastic memory. The
tests use authored numerical fixtures and do not establish experimental material
accuracy or a general mesh/time convergence guarantee.

During implementation, four native tests of the unchanged original mass-kernel
source passed against independent barycentric and physical-quadrature oracles.
The added package integration tests were not compiled or run: the checkout was
missing all seven pinned constellation siblings, and DSR/RCH were unavailable.
Once that environment is present, run:

```sh
cargo test --locked -p fs-solid --test hyper3_reference_density --test load_history3 --test hyper_visco_dynamics
cargo test --locked -p fs-couple --test finite_solid
```
