# CONTRACT: fs-scene

## Purpose and layer

Layer L3. Static robot environments share yawed boxes and solid half-spaces,
so a floor or table participates in geometric checks instead of existing only
in a renderer. Dependencies are `fs-ga`, `fs-geom`, and `fs-query`.

## Public types and semantics

- `SceneBody::Box` stores a center, positive half-extents and rotation about
  world +Z. Coordinates, extents and overlap depths are metres; yaw is radians.
- `SceneBody::HalfSpace` denotes `normal dot position <= offset`. The normal
  points into free space. Admission normalizes its direction; the offset is
  already a distance and is **not** divided by the original normal length.
  `ground_plane(h)` is solid below world height `h`.
- `BodyRole::{KeepOut, Support}` labels a body's purpose. Both use the same
  caller-declared nonnegative `skin_m`; role does not change the overlap law.
- `StaticScene::push` validates one body before appending it and returns its
  declaration-order index. `entries`, `len`, and `is_empty` inspect the scene.
- `deepest_sphere_penetration` scans spherical colliders against all bodies.
  It returns the greatest excess over skin, with body and collider indices,
  raw depth and role. Equal excess retains the first body/collider pair.
  Empty scenes or collider lists return `None`.
- `SceneBody::sphere_overlap_depth` and `sphere_box_overlap_depth` evaluate
  closed-form nominal overlap. A sphere outside a box uses nearest-point
  distance; a center inside uses radius plus distance to the nearest face.
- `bounding_radius_m` and `center_m` return `None` for an unbounded half-space.
  `convex_support_map` constructs `fs-query::ConvexOrientedBox` for a usable
  bounded box; `None` also covers a refused support-map construction.
- `prepare` allocates a vector of `PreparedBody` values with cached yaw and
  bounding-sphere data. `may_breach` provides a preliminary nominal sphere
  reject; half-spaces always return true. Prepared overlap uses cached yaw.

## Invariants

Scene admission requires finite centers, extents, yaw, plane components,
offset and skin; positive box half-extents; nonnegative skin; and a normal
whose computed Euclidean length exceeds `1e-12`. A refused push leaves the
existing entries unchanged. Ordinary admitted half-space normals are unit
directions, subject to the numerical limits below.

Queries do not mutate the scene. Penetration requires depth strictly greater
than skin. Skin applies to KeepOut as well as Support, so KeepOut with a
positive skin intentionally tolerates that much overlap.

Direct body/helper queries accept publicly constructible values without
calling scene admission. Callers must supply finite coordinates, nonnegative
finite radii, positive finite extents and consistent yaw sine/cosine values.
Do not mutate a prepared body's public geometry after preparation: its cached
private values would no longer describe that geometry.

## Error model

`SceneError::NonFinite` and `NonPositive` name the rejected field;
`DegenerateNormal` rejects an unusable normal. Zero skin is allowed despite
the shared NonPositive diagnostic wording. Query methods return nominal
numbers or `Option`, not typed arithmetic-error or certificate objects.

Arithmetic is ordinary `f64`. Closed form does not mean outward-rounded or
certified. Extreme finite components can overflow norm, distance or bounding
radius calculations; normal admission does not currently reject an infinite
computed norm. Unchecked malformed query values can yield misleading results.
Neither prepared nor direct broad-phase rejection has a certified rounding
margin. Claims of separation require a separately admitted `fs-query` path.

## Determinism class

Fixed declaration order, collider order and strict-greater tie breaking.
There is no RNG or scheduling. Repeated queries on an unchanged scene and
the same floating-point environment follow the same operations. No cross-ISA
bit identity is claimed for trigonometry, square roots or compiler arithmetic.

## Cancellation behavior

Synchronous queries have no `Cx`, cancellation poll, deadline or explicit
work budget. A scene scan performs at most bodies times colliders pair tests.
Queries allocate no per-pair buffers; `push` and `prepare` allocate storage.
The caller owns batching and cancellation around scans.

## Unsafe boundary

`#![forbid(unsafe_code)]`; no unsafe implementation or FFI.

## Feature flags

None. These primitive scene operations are available in the default build.

## Conformance tests

Existing tests in `src/lib.rs` cover clear/touching/interior box queries,
yaw, ground planes, support indentation, deepest-hit identity, refused body
admission, normal direction normalization and the convex support bridge.
Sampled tests compare broad-phase and exhaustive scans, prepared and direct
depths, and the unbounded-body reject. They are nominal finite-fixture
regressions, not exhaustive interval or adversarial-magnitude proofs.

## No-claim boundaries

No contact-force integration, friction, dynamics, continuous collision
detection, swept volumes, deformable bodies, arbitrary mesh collision,
distance enclosure, authenticated geometry, physical validation or general
collision-free trajectory certificate. BodyRole is metadata plus a skin
policy; it does not supply a support-force or contact-stability model.
