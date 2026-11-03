

## Prepared full tensor mechanics with exact shared relaxation history

`tensor_material::prepared::PreparedTensorBody` provides a fixed-clock cached
route for the original full three-dimensional tensor material, consistent inertia
and force/velocity coordinates. It consumes the existing `ProjectedTensorProny`
points, retaining every equilibrium tensor, Kelvin factor, Maxwell factor,
positive quadrature volume and exactly represented relaxation time. It does not
fit modal damping, lump mass, merge nearby poles, remove weak material directions,
or change the caller's displacement basis.

For each BIT-IDENTICAL positive source time tau, one complete displacement-memory
vector satisfies `y_tau' = M^-1 p - y_tau/tau`. Original point/branch springs are
`z_pb = A_pb y_tau`, and the mechanical force is
`p' = f - K_eq u - C M^-1 p - sum_tau K_tau y_tau`, where
`K_tau = sum_(pb at tau) A_pb^T A_pb`. Every source factor remains available for
original storage. Thus the dynamic state has `(2 + distinct_times)*coordinates`
entries instead of six extra states for every branch at every quadrature point.
This is an exact shared-coordinate realization for the admitted source history,
not a reduction of mechanical modes or a constitutive approximation.

The constructor binds the source and clock; borrowing histories explicitly choose
`InitialMemory::Relaxed` (all y=0) or `Instantaneous` (each y=u0). No prior work or
heat is invented for that initial state. Nonzero free eigenstrain refuses because
its initial offsets need not belong to the displacement-driven source image.
Arbitrary imported point memories, changing laws/geometry/supports, plane-stress
descriptors, finite rotation and thermal/source rebinding retain their separate
owners. Endpoint source-strain admission is enforced; a within-interval strain
or frequency-band enclosure is not claimed.

`PreparedQuadraticFlow` caches the SAME native coupled generator and propagation
kernel once. Each declared heat group uses `PreparedSourceLoss` on the ORIGINAL
positive Kelvin and Maxwell source squares. Diagonal similarity conditioning
changes only numerical coordinates. When source rows exceed the state dimension,
the existing native Householder QR supplies a full factor without rank truncation;
its normalized Gram matrix must match the original source Gram, including cross
terms. Group labels are explicit, contiguous and nonempty; grouping changes
observation cost and attribution only. One original element can collect all its
quadrature heat in one cache.

A warm interval performs cached contractions, original source-domain and energy
checks, and independently integrates held-force work `f dot (u1-u0)`. Equilibrium
and kinetic storage use original K_eq and M; Maxwell storage is the sum of the
squares of EVERY original A_pb*y_tau. Stable quadratic increments and compensated
work/heat histories check `W - D - Delta H` per interval and cumulatively before
publishing state, clock or ledgers. Heat is a positive source observation, never
an energy-deficit fill. Finite force jumps apply no impulse. Histories clone only
O(state) accepted data and scratch; immutable dense caches are borrowed.

Explicit capacities cover coordinates, points, complete states, heat groups,
ALL preparation/cache bytes and multiply terms, warm work and accepted steps.
Preliminary grouping is admitted before allocation. The existing native dense
limits remain: the base split, extended flow and doubled positive-loss state
must each fit their own 256-coordinate bound. No capacity escalation or hidden
truncation occurs. Preparation includes original Cholesky/eigensolve/QR and
bounded flow-integral calls; cancellation is checked around those existing
atomic native kernels and throughout source/cache work. No hard cold-path
cancellation-latency or real-time timing claim is made.

Six `tests/tensor_material_prepared.rs` tests pass natively: an independently
solved standard-linear-solid with exact roots -1,-2,-3 and analytical dashpot heat
for both prehistories; separate analytical heat for unequal elements with reversed
channel labels; coupled anisotropic, nonproportional Kelvin and consistent-mass
agreement with the original uncompressed owner (including every reconstructed
spring coordinate); exact versus adjacent-f64 source times; capacity/eigenstrain
admission; and bitwise rollback for malformed input, material-domain refusal,
cancellation, exhausted steps and work-budget refusal. The native source-closure
build uses the original material law/projection/history, full fs-phs and native
factor/eigen/GEMM/SIMD kernels. It is not a full workspace Cargo or DSR pass.
