# Authored constructive domains for native 3-D studies

The native `sdf3-study` path accepts `constructive-implicit` domains. The same
CutFEM operator builder consumes them for independent-load compliance/SIMP,
goal-guided octree refinement, and fixed-background stress-constrained
minimum-volume studies. The shape is fixed while material densities are
optimized. This is not moving-boundary or CAD parameter optimization.

```sh
cargo run --release -p fs-cli --features sdf3-study --bin frankensim -- \
  --json study examples/marquee/bracket-3d-bored.fsim bored-study.db
```

The supplied one-stage example is a numerically scaled plate with a transverse
bore. It retains SI declarations and explicit work limits; its modulus and loads
are demonstration values, not a validated material/design specification. Add
SIMP stages to its existing schedule to request goal-guided background refinement.
Each completed stage uses the ordinary ledger/checkpoint path. Compliance resume
replays the accepted prefix under the same executable and charges replay to the
original allowances; stress-study cross-process resume remains unsupported.

Replace only the domain section of an existing version-1 3-D study:

```lisp
(domain
  :type constructive-implicit
  :bounds ((0.0 0.0 0.0) (1.0 1.0 1.0))
  :shape (difference
    :blend-m 0.0
    :left (half-space :normal (0.0 0.0 1.0) :offset-m 0.63)
    :right (cylinder :axis y :center-m (0.55 0.0 0.28) :radius-m 0.12)))
```

## Shape forms

Fields are required in the shown order; unknown or repeated fields are refused.
All coordinates, semi-axes, radii and offsets use the declared SI frame.

```lisp
(half-space :normal (0.0 0.0 1.0) :offset-m 0.63)
(sphere :center-m (0.5 0.5 0.5) :radius-m 0.3)
(ellipsoid :center-m (0.5 0.5 0.5) :semi-axes-m (0.4 0.3 0.2))
(cylinder :axis z :center-m (0.5 0.5 0.0) :radius-m 0.12)
(box :center-m (0.5 0.5 0.5) :half-extents-m (0.3 0.2 0.1))
(union :blend-m 0.0 :left SHAPE :right SHAPE)
(intersection :blend-m 0.0 :left SHAPE :right SHAPE)
(difference :blend-m 0.0 :left SHAPE :right SHAPE)
```

Cylinders are infinite along x, y or z and clipped by the background bounds or
another shape. A half-space normal is used literally, not silently normalized.
Zero blend width selects exact hard set operations. Positive width deliberately
changes the boundary through a polynomial implicit blend; it is a width in the
length-valued **field**, not a guaranteed Euclidean fillet radius. Composite and
quadric magnitudes are not exact signed distances. Sharp child creases remain
sharp even when a parent operation is blended.

Recipes are capped at 128 nodes and 24 nesting levels. The usual physical box,
material, load, field-memory and quadrature/solver limits still apply. Interval
ranges and partial-slope ranges guide the existing quadrature; unresolved cuts,
unsupported/disconnected mechanics, and work exhaustion refuse rather than
publishing a fabricated successful design. Primitive support is not a promise
that every Boolean arrangement will resolve at every requested grid/budget.

Reference pressure/traction patches act on **all retained implicit surfaces in
the declared x-fraction patch**, including hole/cavity walls, not clipping-box
faces. There is no surface-name selector or follower-load semantics in this path.
Body loads remain independent cases. Geometry expressions are retained in the
canonical source and affect study identity and replay automatically.

The old `curved-height-sdf` and `physical-curved-height-sdf` declarations retain
their previous arithmetic and envelopes. No new continuum safety, exact-distance,
physical-validation, or globally optimal-design claim is introduced.
