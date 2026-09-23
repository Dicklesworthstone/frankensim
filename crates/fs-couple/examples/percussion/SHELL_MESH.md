# Explicit three-dimensional cymbal geometry

`--shell-mesh instrument.fss` replaces the revolved profile with an explicit
triangular reference midsurface, nodal physical thickness and per-facet material
assignments. It is intended for supplied surveys or prepared CAD meshes whose
asymmetry, holes and local geometry must not be replaced by a radial average.
It works with `splash`, `splash-wav` and `splash-mic` through the same shell FEM,
nonlinear reduction, contact and two-sided radiation owners as the profile path.

## Start from the existing estimated instrument

The cold export command writes the declared profile as explicit samples, without
running an eigensolve or playback. It refuses to overwrite an existing file.

```sh
cargo run --release -p fs-couple --example percussion -- \
  export-shell-mesh cymbal.fss

# The optional third argument replaces the built-in estimated profile.
cargo run --release -p fs-couple --example percussion -- \
  export-shell-mesh supplied-cymbal.fss measured-meridian.profile

# Edit or replace the mesh with actual supplied geometry and physical fields.
cargo run --release -p fs-couple --example percussion -- \
  splash 4096 --shell-mesh cymbal.fss --analytic-newton --impact-substeps 8 511 \
  --strike-position-m 0.06 0.01 --second-stick-position-m -0.05 0.02 \
  --second-stick-speed-m-s 0.8 > mesh-cymbal.csv

# The same geometry also reaches both actual fixed-receiver BEM observations.
cargo run --release -p fs-couple --example percussion -- \
  splash-mic 4800 20 --shell-mesh cymbal.fss --analytic-newton \
  --impact-substeps 8 511 --microphone-right -0.08,0.05,0.35 > mesh-cymbal.wav
```

The default export is the **existing estimated splash**, not a newly acquired
scan, measured digital twin or manufacturer's geometry. A supplied file keeps
its input status; parsing never confers measurement or calibration evidence.
Usage examples do not establish a converged or auditioned audio render.

## Physical input format

The bounded UTF-8 format begins with `frankensim-shell-mesh-v1`. Blank lines and
`#` comments are allowed. All physical quantities use SI units. Records are:

```text
band_hz,lower_Hz,upper_Hz
strike,default_x_m,default_y_m
material,material_ID,Young_Pa,Poisson_ratio,density_kg_m3
node,node_ID,x_m,y_m,z_m,thickness_m
triangle,triangle_ID,node_a_ID,node_b_ID,node_c_ID,material_ID
```

`band_hz` and `strike` occur exactly once. The explicit default strike is needed
because a general mesh has no unique radial meridian from which to infer one.
`--strike-position-m X Y` still overrides it for a particular performance.
Node, triangle and material IDs are distinct nonnegative integers in their own
namespaces; they may be sparse and records may be reordered. Dense indexing is
deterministic by source ID, and each triangle keeps its corner order. Every
material and vertex must be referenced. Unknown records, duplicate IDs, missing
assignments, invalid physical values, and inputs over 8 MiB refuse.

Positions are **midsurface** coordinates in the instrument's reference frame,
with positive z upward. Thickness is a positive physical field, not a visual
normal map. Each facet uses its three thickness samples' arithmetic mean in the
existing isotropic section law; mass uses its actual three-dimensional area.
This retains the profile path's section approximation, rather than claiming
exact integration of bending stiffness through an arbitrary thickness field.
Young's modulus, Poisson's ratio and density belong to the assigned material.
No alloy-name lookup or unprovided residual stress/hardness is inferred.

One connected, consistently oriented manifold shell is admitted. Boundaries and
multiple holes are retained. Exact duplicate positions/faces, disconnected
pieces, unused vertices, nonmanifold edges and pinched vertex fans refuse. The
adapter does not weld, smooth, resample, fill holes or repair winding. Near-
coincident geometry and global self-intersections are **not certified**.

This percussion host still uses vertical XY contact stations. Facets must face
upward and have nondegenerate XY projection; vertical walls and local folds do
not fit that chart. A station on a shared topological edge or vertex is legal;
multiple unrelated projected facets at a queried station refuse instead of
selecting one by input order. This local query check is not a global overlap
certificate. Supply a physical midsurface rather than a closed solid CAD skin.

## One geometry throughout the instrument

The exact admitted vertices and connectivity reach shell assembly and nonlinear
strain reduction. The nodal thickness also defines the finite-thickness
radiation lift, including both faces and boundary walls. Per-facet sections
control actual stiffness and mass. Nothing is substituted into an output EQ.

Both sticks, fixed mufflers and compliant mutes sample the same supplied shell.
Independent force programs, analytic Newton, internal recovery, felt history and
mono/stereo observers remain available. `--shell-profile` and `--shell-mesh` are
mutually exclusive and reject before either input is opened; drum/snare commands
reject both shell selectors rather than consuming an unrelated geometry.

The estimated stand remains unchanged: its pads at 12 mm radius must lie on the
supplied surface. A differently mounted instrument needs different hardware;
this importer does not invent it. The retained vertical rigid translation and
elastic basis are not a full six-degree-of-freedom suspension model.

All existing limits remain: 10,000 nodes, 20,000 triangles, 32 shell coordinates,
physical and timestep guards, reduction budgets, and the 2,048-panel radiation
budget. The fixed-receiver bake still uses 40–1640 Hz. Wider input windows may
exceed these limits and refuse; they do not silently truncate modes or establish
full-band cymbal fidelity. Spatial/modal/temporal convergence, manufacturing
stress, measured damping, two-way radiation loading and real-time qualification
remain separate work. There is no new time integrator or runtime dependency.

Focused native checks:

```sh
cargo test --release -p fs-plate --test percussion_shell survey
cargo test --release -p fs-couple --example percussion mesh_input -- --test-threads=1
cargo test --release -p fs-couple --example percussion playing::tests
```
