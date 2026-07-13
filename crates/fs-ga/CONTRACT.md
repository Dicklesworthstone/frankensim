# fs-ga — CONTRACT

The geometric-algebra layer (plan §7.7, Bet 2): PGA Cl(3,0,1) as the
kinematics substrate — motors/screws replace the quaternion + matrix +
Plücker zoo and remove gimbal-class bugs by construction — plus CGA
Cl(4,1) for sphere/tangency-rich construction. Multiplication tables are
CONST-EVALUATED from the metric signatures; conventional Vec3/quaternion/
matrix façades sit at every API boundary.

Ambition tags: PGA motors/exp/log/incidence [S]; CGA rounds/tangency [S];
const-evaluated codegen + monomorphized kernels [S].

## Purpose and layer

Layer **L2** (MORPH). Runtime deps: `std`, `fs-math` (deterministic
`sin`/`cos`/`atan2`/`sqrt` so motor exp/log is bit-identical across
platforms). Consumers: fs-time SE(3) integrators, fs-solid rods/joints,
fs-scenario frames, the gradient stack's manifold variants.

## Public types and semantics

- `table`: `build_table` (const fn) computes the full Cayley table of any
  diagonal-metric Clifford algebra at COMPILE TIME; `PGA_TABLE` (16×16,
  metric (0,+,+,+)) and `CGA_TABLE` (32×32, metric (+,+,+,+,−)) are
  compile-time constants — no runtime blade bookkeeping anywhere.
- `Pga`/`Cga`: dense multivectors with `gp`, `wedge`, `lcontract`, `vee`
  (regressive via Poincaré right/left complements — correct in the
  degenerate PGA metric), `reverse`, `involute`, `grade_part`, `norm_sq`.
- `pga`: `Point` (trivector, ganja layout `e123 + x e032 + y e013 +
  z e021`), `Plane` (`a e1 + b e2 + c e3 + d e0`), `Line` (join of points
  / meet of planes), `Motor` (even versor) with `compose`,
  `transform_point`, `transform_plane`, `renormalize`, `slerp`;
  `exp_bivector`/`motor_log` via exact screw decomposition
  `B = θℓ + dℓ*`, `exp B = (cos θ + sin θ ℓ)(1 + d ℓ*)`.
- Monomorphized kernels (const-generated from the Cayley table): even⊗even
  (composition, 64 fused terms), even⊗point and odd⊗even→point (sandwich).
  `transform_point_dense` is the dense reference path, cross-checked in
  conformance.
- `cga`: `up`/`down` null embedding, `n_o`/`n_inf`, `sphere_through`
  (4 points), `circle_through`, `plane_through`, `dual_sphere`,
  `sphere_center_radius`, `incidence`, `tangency_residual`.
- `facade`: `Vec3` (+, −, dot, cross), `Quat` (Hamilton `*`, `rotate`,
  `from_axis_angle`), `Mat34` (motor lowered to a rigid-motion matrix),
  `Motor::from_parts`/`to_parts` — the no-formalism-tax boundary.

## Invariants

1. **Generated, not written**: every product sign comes from the const
   evaluator; the table itself is test-audited for blade-level
   associativity (exhaustive over all triples, both algebras), metric
   squares, and anticommutation.
2. **G0 identity battery**: associativity, distributivity, grade-partition,
   reverse anti-automorphism, wedge antisymmetry on randomized dense
   multivectors in both algebras.
3. **Motor sandwich = rigid motion**: agrees with the independent
   Rodrigues + translation matrix reference to < 1e−11 absolute over
   randomized screws (worst observed ~1e−13); the monomorphized kernel
   agrees with the dense path to < 1e−12.
4. **exp/log**: mutually inverse on the principal branch (< 1e−10 over
   randomized screws; pure translators exact through the θ→0 branch);
   `exp` emits unit motors (defect < 1e−12).
5. **Façade exactness**: quat ↔ rotor is a bitwise relabeling both ways;
   motor → (quat, translation) recovers the quat bitwise and the
   translation to ≲1e−13 relative.
6. **Versor hygiene**: `renormalize` divides out the full `M M̃ = a + bI`
   residue; a 20 000-product chain renormalized every 64 steps stays
   below 1e−11 drift (ledgered by the conformance run).

## Error model

`GaError`: `IdealPoint` (a trivector with zero e123 weight has no
Cartesian form), `ZeroWeight { context }` (conformal element with no
finite representative — flat blades, degenerate sphere construction).
Degenerate constructions return the honest degenerate blade or a
structured error; nothing silently normalizes garbage.

## Determinism class

**D0 on-target and cross-ISA**: products are fixed-order f64 chains over
compile-time integer tables; all transcendentals go through
`fs_math::det`. No platform intrinsics, no reassociation.

## Cancellation behavior

All operations are small, bounded, allocation-free (arrays on the stack)
— P7 satisfied by boundedness; there are no long-running loops.

## Unsafe boundary

Zero `unsafe`.

## Feature flags

None.

## Conformance tests

`tests/conformance.rs` (JSON verdicts, suite `fs-ga/conformance`):

- **ga-001** fixed identity battery (seed `0x6A_0001`, 200 rounds):
  associativity, distributivity, grade partition, reverse, and vector-wedge
  antisymmetry in PGA; associativity and reverse in CGA. Bead 4nh8 adds the
  missing complete generated coverage without changing this fixed pin:
  five 400-case PGA suites (`0x6A_1001..0x6A_1005`) and five 400-case CGA
  suites (`0x6A_2001..0x6A_2005`), all shrink-armed. Approximate residuals
  require every coefficient to be finite before applying the 1e−12 PGA or
  1e−11 CGA bound; grade partitions compare exactly.
- **ga-002** 1200 motor sandwiches vs Rodrigues+translation; kernel vs
  dense cross-check; `Mat34` lowering; tiny-rotation (1e−13 rad)
  stability near the screw-axis singularity.
- **ga-003** 300 exp/log round-trips (screws + exact translators); 32-step
  slerp through a 90°-pitch gimbal fixture: stays a unit motor with
  constant screw speed (spread < 1e−9), endpoints met.
- **ga-004** join/meet incidence; the plane incidence measure equals the
  implicit plane equation; meet lines lie in both parent planes.
- **ga-005** CGA: null embedding, the conformal distance law
  `up(a)·up(b) = −½‖a−b‖²`, sphere-through-4-points center/radius
  recovery, fifth-point incidence, circle/plane incidence, tangency
  residual (zero for tangent spheres, large for separated ones).
- **ga-006** façade exactness (bitwise quat relabel, ULP translation).
- **ga-007** versor drift statistics under the renormalization policy.

Unit tests audit the generated tables themselves (exhaustive blade-level
associativity, metric squares, complements, reverse/involution signs).

## Performance (ledgered by `ga_bench`)

Doctrine: **compose in motor land, apply in matrix land** (`Mat34`
lowered once per motor). Measured on the shared Linux worker (JSON output
of `cargo run --release -p fs-ga --bin ga_bench`): bulk point transform
via lowered `Mat34` runs at 0.5–0.75× the hand-written quaternion path's
cost (i.e. faster); motor composition runs at 1.9–3.7× the quat+vector
compose cost (64 fused monomorphized terms vs ~46 flops — the price of
carrying screw semantics and drift correctability); the raw per-point
sandwich is ~14× quat and is NOT the bulk path — it exists for
correctness cross-checks and one-off transforms.

## No-claim boundaries

- **No SIMD capsules yet.** The kernels are scalar (compiler-unrollable
  fixed loops); fs-simd-tiered versions are follow-up work once fs-render
  or fs-time demand them.
- **Raw-sandwich parity is not claimed.** Parity holds for the lowered
  bulk path and near-parity (≲2–4×) for composition; the dense 16- and
  32-component products are correctness machinery, not hot paths.
- **CGA covers construction, not kinematics.** Conformal versors
  (rotors/translators/dilators in Cl(4,1)) are out of scope for this
  bead; rigid motion lives in PGA. Apollonius-type solvers are expressed
  through meets/tangency residuals, not packaged as a dedicated solver.
- **Principal branch only.** `motor_log` folds the ±M double cover and
  returns screws with θ ∈ [0, π]; callers interpolating past π must
  compose increments.
- Timing numbers above are from a shared, sometimes-loaded worker;
  ratios (not absolute ns) are the claim.
- The generated dense G0 battery covers coefficients in `[−1, 1)` and is a
  replayable defect detector, not a formal proof over arbitrary-magnitude f64
  multivectors. Exhaustive blade-table tests remain the exact algebraic
  structure check.
