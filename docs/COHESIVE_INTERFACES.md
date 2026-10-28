# Cohesive interface failure and structural load sharing

`fs-material::cohesive` and `fs-solid::cohesive` implement a bounded part of
issue #1 R4: explicitly bonded surfaces can soften, unload, irreversibly fail,
and redistribute forces through a supplied elastic structure. They do not
infer adhesion parameters or whether an interface should initially be bonded.

## Run the physical example

```sh
cargo run -p fs-solid --example cohesive_support_removal
```

The example puts a support-bond cohesive plane and a declared notched-feature
fracture plane in series with a compliant displacement-controlled actuator.
One run assigns the weaker law to the support interface; the second assigns
it to the feature plane. It prints actual reaction force, both openings and
both damage variables, and checks equal series forces and the expected
failure location. Parameters are illustrative, not measured print properties.
A predeclared cohesive feature plane is not automatic bulk crack prediction.

## Constitutive and structural API

Create `BilinearCohesive` with stiffness (Pa/m), peak traction (Pa), fracture
energy (J/m2), compression penalty (Pa/m) and a separation limit (m). The
onset and complete-failure openings are `strength/stiffness` and
`2*fracture_energy/strength`. Normal tension and both local shear modes share
one effective-separation history and the same parameters. Compression alone
does not damage the bond. After complete failure, the normal compression
penalty remains; frictionless penalty reclosure is not exact nonpenetration.

`bonded_point()` starts a zero-jump, stress-free bond. `point.trial(jump, gate)`
returns a candidate without mutating the accepted history. Unloading is
linear at the retained damaged stiffness, with no healing. The returned
traction, stored-energy change and fracture dissipation are per original area.
Work is analytically integrated along the declared straight separation path;
unreported intermediate excursions cannot be reconstructed from endpoints.

For structural coupling, bind a constant `fs_solver::LinearOp` to
`CohesiveSystem`. Each `CohesivePatch` has a stable ID, original area, law and
sparse `JumpTerm` columns. Columns map reduced displacements to local jumps.
Paired nodes use equal/opposite columns; supplied trace interpolation may use
more. Resolve fixed local frames, quadrature and eliminated constraints before
binding. The operator and maps must stay unchanged through a state's lifetime.

`state.solve(force_n, prescribed_jumps_m, config, cx)` solves

```text
K u + sum_i B_i^T A_i traction_i(B_i u + prescribed_jump_i; old_history_i) = force
```

The existing Newton/FGMRES driver handles the generally indefinite softening
tangent. Every nonlinear probe uses the original accepted damage, not the last
Newton trial. A successful result publishes displacement and all interface
histories together only after convergence, a fresh force-residual check and
local/aggregate cohesive energy checks. Bounds cover dimensions, sparse map
size, Krylov basis entries, iterations and operator evaluations. Cancellation
and any refused trial leave the accepted state unchanged. An invalid trial
refuses the increment rather than being used to extrapolate the material law.

Inspect `result.state.points()` in the same order as `system.patches()` for
jump, traction, damage and failure. These observables can inform an explicit
thermal/contact update. Thermal conductance and cohesive strength remain
independent; this module does not convert one into the other.

The separate `fs-conduction` `FractureCooling::synchronize_with_heat` adapter
can deposit an explicitly supplied fraction of NEW cohesive dissipation into
the original thermal cells. It retains nonthermal fracture energy separately,
uses the original enthalpy law for latent heat, and commits contacts, damage and
heat together. Prior damage, unloading and numerical defects do not supply heat;
the adapter neither infers conversion fractions nor changes this constitutive law.

## Applicability and verification

The model is quasistatic, fixed-topology and fixed-frame. It does not establish
unique/stable equilibrium, resolve dynamic snap-back, create new bonds, update
geometry, supply healing/friction/temperature-dependent adhesion, or couple
fracture heat into a thermal solver. Subdivide loading paths to resolve events
and reversals. Reported energy is cohesive-only, not a complete actuator,
structural or thermodynamic balance. No whole-print qualification is claimed.

Six material and seven structural Rust regressions cover fracture energy,
normal/shear unloading and reclosure, tangents, nonproportional paths,
work-conjugate force assembly, prescribed motion, area partition, load sharing,
refusal/cancellation atomicity and replay. Native Rust execution is still
unverified in the authoring environment. Independent Python/SciPy checks
exercise the equations, not the Rust compiler or the native solver runtime.
