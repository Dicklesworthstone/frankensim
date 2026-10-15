# fs-mbd-material

L3 composition boundary from immutable material data and volume geometry to
`fs-mbd::articulated::SpatialInertia`. The numerical `fs-mbd` owner remains
independent of material databases; no replacement dynamics solver is introduced.

`resolve_body` requires a unique key, volume moments, material card, complete
state point and source selection for every region. Only density is required.
It retains the original resolved bundles and composes their mass, centre of
mass and full centroidal tensors, including products of inertia. Region keys
fix the reduction order. The parallel-axis merge works about the centroids,
not by subtracting origin tensors. Final spatial-matrix admission is performed
by the same downstream inertia owner used for momentum and kinetic energy.

`VolumeMoments` accepts moments from a geometry integration owner and supplies
analytic cuboid and straight-tetrahedron constructors. SE(3) placement uses
`fs-ga`; tensors remain centroidal after placement. Lengths are metres, volume
is m3, volume inertia is m5, density is kg/m3 and mass inertia is kg m2.

The caller supplies a non-overlapping rigid partition. Numeric admission does
not establish geometry coverage, mesh accuracy, rigidity or measured material
accuracy. Source uncertainty is retained, not propagated to a certified inertia
bound. Density at a changed state must be paired with the intended geometry:
compiling a new specimen is not a mass-conserving thermal-expansion update, a
live state transfer, or an implicit source of energy/momentum.

Compilation has a 4096-region hard cap and a caller-selected lower cap. Missing,
ambiguous, wrong-unit, nonpositive or unsupported density refuses. Cancellation
is checked before work, between fixed-size region operations, and before
publication. Sorting and each underlying material query are bounded cold-path
operations, not a real-time cancellation-latency guarantee.
