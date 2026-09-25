# Physical shaft bending during a hi-hat performance

`--flexible-stick INPUT.fst` and `--second-flexible-stick INPUT.fst` replace the
selected rigid effective stick with its geometry-derived pin-supported shaft.
Both flags work with `hihat`, `hihat-wav` and `hihat-mic`. The second requires
`--second-stick-position-m X Y`; either stick can remain on the original rigid
path. These flags are not yet admitted by the single-drum/splash commands.

Tip contact now acts on the shaft's rigid rotation **and every retained bending
mode**, with the same reciprocal reaction on the upper cymbal. A strike can
store energy in shaft bending and return it through later contact. Nothing
adds a recorded attack, authored bending frequencies, duplicate time stepper,
output filter, or independently synthesized stick sound.

The source is the existing `fs-plate` round Hermite/Rayleigh beam pencil and
`fs-couple::render::plate::impact::striker::flexible` owner. The complete shaft,
including its butt behind the pin, contributes consistent translational and
section-rotation inertia. Shaft modes are appended after the original rigid
stick/shell/carriage prefix. Cymbal acoustic source addresses stay unchanged;
only the actual resulting shell motion reaches the pressure observer.

## Explicit shaft input

Every singleton record is required once; comments and blank lines are allowed.
The strict UTF-8 format has an 8 KiB ceiling. The included file is a **synthetic
uniform test shaft**, not manufacturer geometry or measured wood properties:

```text
frankensim-flexible-stick-v1
material,12000000000,800,0.001
support,0.1,0.39,0.16
basis,8,3000,17
station,0,0.005
station,0.4,0.005
```

`material` specifies Young's modulus [Pa], density [kg/m3], and a modal damping
ratio in [0,1). The true rigid rotation has zero stiffness and zero damping.
`support` specifies pin, tip-contact, and hand-force axial stations [m] on the
same shaft, with positive tip and hand levers. The pin fixes transverse
translation, not rotation; it is not an inferred grip impedance or constraint
on the tip. `station` supplies 2..33 ordered axial positions and radii [m],
with linearly tapered radius between them. No adjacent pair can both be zero.

`basis` supplies subdivisions per segment, maximum retained frequency [Hz],
and total mode ceiling including rigid rotation. The existing bounds are
1..16 subdivisions, 66 source nodes, and 2..17 modes. Every elastic mode within
the requested frequency window must fit; overflow refuses instead of dropping
modes. The unchanged mechanical Nyquist guard applies. A complete discrete
slice is not a spatial/bandwidth convergence certificate.

## Hand force, not tip force

With a flexible stick selected, its existing `--stick-force-file` or
`--second-stick-force-file` acts at the supplied **hand station**. One force
program is integrated once over each mechanical tick, then projected through
the full signed hand row. Modal signs and lever ratios are retained. Work is
`F_hand * v_hand`, not `F_hand * v_tip`. Without a flexible selection, the old
scalar force-port behavior remains unchanged. Pedal and both hands share the
same accepted clock, and a refused/cancelled step consumes none of them.

Initial launch still uses the explicit strike speed and 0.2 mm tip clearance.
The entire shaft starts in rigid motion; no initial bending energy or arbitrary
mode phase is injected. Both sticks can be launched and driven independently.
CSV adds tip/hand displacement, tip/hand velocity and bending energy only for
selected flexible sticks; that energy is already included in `total_energy_j`.

```sh
cargo run --release -p fs-couple --example percussion -- \
  hihat crates/fs-couple/examples/percussion/estimated-hihat.fshh 12000 \
  --flexible-stick crates/fs-couple/examples/percussion/estimated-flexible-stick.fst \
  --second-stick-position-m -0.06 0.01 --second-stick-speed-m-s 0.6 \
  --second-flexible-stick crates/fs-couple/examples/percussion/estimated-flexible-stick.fst \
  --analytic-newton --impact-substeps 8 511
```

This is a supported invocation, not a claimed completed render. Existing
squeeze-film, radiation, receiver and force-file options remain composable.
The tip's existing estimated Hertz law is separate and unchanged: importing a
shaft does **not** calibrate tip curvature, transverse contact modulus or wood
anisotropy. Planar small-deflection bending does not model three-dimensional
stick orientation, torsion, shear, changing grips, or hand biomechanics. Large
bending and impact bandwidth still need physical validity/refinement studies.
Native regressions cover signed hand work and retry, two-shaft mode/source
addresses, rigid-only launch, and actual contact excitation with shared energy.
Native execution, measured force/response comparison and throughput remain
required before instrument-fidelity or real-time claims.
