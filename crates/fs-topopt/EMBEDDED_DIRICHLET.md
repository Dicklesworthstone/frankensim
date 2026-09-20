# Supports and prescribed motion on a 3-D implicit surface

The 3-D Cartesian and adaptive elasticity builders now accept
`build_with_embedded_dirichlet`. A pure `(position, outward_normal) -> bool`
callback selects the full-vector displacement patch. The box-clamp callback
may select no nodes: supports can lie strictly inside the background box.
The old builders still reject missing box support.

The symmetric Nitsche cell contribution is
`-v dot sigma(u)n - sigma(v)n dot u + gamma v dot u`, with
`gamma = beta*(lambda+2*mu)/min(background cell spans)`. It is added to the
existing cell matrices, not a second assembled solver. Density multiplies the
bulk AND boundary cell terms; ghost terms retain their original density law.
Hanging-node elimination, Jacobi/Galerkin preparation and density contractions
therefore see the same modified operator. The default beta is 32, not a theorem
that every cut, material distribution or disconnected component is stable.

## Runnable consumers

```sh
# Homogeneous embedded support; two independent right-face loads; two grids.
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example embedded_supported_sdf3 -- 3 250000 2
# Nonzero imposed translation plus external traction; prints physical fields.
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example embedded_supported_sdf3 -- --motion
```

Use the repository's DSR/RCH execution lane when available. The optimization
example uses the dimensionless slab `0.17 < x < 0.83`, supported on its LEFT
implicit face. Neither support nor loaded face coincides with a grid boundary.
The compression/shear callbacks explicitly vanish on the supported patch;
Neumann loads are not silently clipped to avoid contradictory user inputs.
Optional refinement uses the existing surface-load DWR estimator, transfers
raw densities, restores volume feasibility and solves a fresh baseline.
It does not assume objective descent across different grids.

## Nonzero displacement is not a fixed force

`prescribed_displacement_load(g, checkpoint)` computes the CURRENT-density
lifting `b_g = integral(-sigma(v)n dot g + gamma v dot g)` on the retained patch.
Add it to the external load before solving. Recompute after EVERY density
change. The callback is a displacement, not a traction; it must agree with any
additional strongly imposed zero box clamps.

`prescribed_displacement_scale_work(g, v, checkpoint)` supplies
`v^T db_g/dscale_c`. For the augmented functional `(f+b_g)^T u`, the scale
derivative is `2*load_work(g,u)-scale_quadratic_forms(u)`. That functional is
NOT physical external compliance `f^T u`, a reaction force or actuator work.
An external-only observation requires its own adjoint and both stiffness and
lifting derivatives. The existing fixed-load OC driver and reference-load DWR
entry points are used here only with homogeneous `g=0`; this example does not
silently treat a nonzero lifting as density independent.

Support selection and normals use retained numerical surface quadrature, not
a certified patch measure or a geometric boundary partition. Resolve narrow
patches and junctions and support every connected component. Inter-grid
admission detects changed boundary method/penalty; the same actual patch and
implicit domain remain caller obligations. DWR cell-operator terms include the
retained Nitsche contribution. No coercivity theorem, nonzero-lifting DWR,
follower pressure, shape derivative or continuum-error bound is claimed.
