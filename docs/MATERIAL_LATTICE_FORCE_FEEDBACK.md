# Material viscosity and buoyancy on the shared 3-D lattice

`fs_lbm::lattice3::sources::material::MaterialForceGrid3` combines the existing
D3Q19/D3Q27 collision and streaming owner with the existing `MaterialScaling`
and `BoussinesqMaterial` owners. Local material changes now drive spatial body
forces and viscosity on an evolving grid; no second fluid solver is introduced.

## Physical inputs and ownership

Construct an ordinary, unstepped `Grid3` with zero uniform force and an explicit
collision model, topology, populations and cell/Mach budgets. Give the material
wrapper that grid, its original `MaterialScaling`, a `BoussinesqConfig<3>` and
one `Option<F>` per cell. `F: FluidTransport` can be an existing gas, mixture,
liquid or other concrete source response. Every fluid needs a source; every
solid needs `None`. Sources and their query receipts remain in `samples()`.

The wrapper owns the scaling for its lifetime. Reference density, dx and dt
cannot be replaced during stepping. Its velocity gate also enforces the
scaling's declared scenario speed, retaining any stricter original grid gate.

Each original source is sampled exactly once per supplied endpoint. The shared
Boussinesq projection checks density/Cp departure bounds and retains rho0/Cp0.
The same response supplies local viscosity and buoyancy:

```
nu = mu_source / rho0
F_SI = (rho_source - rho0) * gravity
F_lattice = (F_SI / rho0) * dt^2 / dx
```

The last conversion uses `MaterialScaling::force_density_lattice`, not an
acceleration shortcut. No second beta*g force or hydrostatic reference load is
added. The chosen collision family is unchanged: D3Q27 preserves independent
bulk and higher-order rates; D3Q19 retains its existing tied second-order rates
and unforced-only restriction for its moment collision families. Force-bearing
endpoints are refused for those unforced families before publication, rather
than accepting a half-kick and failing only on the following collision.

## Endpoint sequencing

Initialize with sources at step zero. Then call:

```rust,ignore
run.advance_serial(&cx, next_step, sources_at_next_endpoint)?;
// Or use the SAME source timing on the fixed-tile executor:
run.advance_pooled(&pool, &gate, next_step, sources_at_next_endpoint)?;
```

The step must be exactly `run.flow().completed_steps() + 1`. A repeated, skipped
or unrepresentable endpoint refuses before a source update. The source field
may be evaluated by an external thermal/composition model; this API does not
invent temperatures or assert that those external values are physically valid.

Collision consumes the current viscosity and current force F_n. The streamed
state is interpreted with endpoint force F_(n+1), using the original half-force
velocity convention:

```
rho * u_n = sum_q(f_q * e_q) + F_n / 2
```

For a periodic grid the total momentum increment is therefore the trapezoidal
source impulse `(sum(F_n) + sum(F_(n+1))) / 2`. New viscosity becomes current
only after acceptance and drives the next collision. This is an explicit,
lagged coefficient update, not a stiff implicit material-feedback solve.

`step_serial` and `step_pooled` hold the current source snapshot without
re-querying it. `sampled_at_step()` does not advance during that hold: a frozen
source is not falsely presented as a newly evaluated material epoch.

Initial populations must follow the documented force-shifted convention.
Starting with all material states equal to the reference gives zero initial
buoyancy and permits ordinary equilibrium-at-rest populations. For a nonzero
initial force, supply populations with the intended force-corrected velocity;
attaching a field does not silently recenter them.

## Publication and limits

Material samples, collision coefficients, force history, populations and step
count publish together only after source admission, both native passes and
final cancellation checks. A source-domain, ownership, conversion, allocation,
Mach, epoch or cancellation refusal leaves all accepted state intact. Producer
callbacks must not mutate shared external simulation history; their side
effects are not rolled back by this owner.

This route advances weakly compressible Boussinesq momentum with prescribed
material samples. It does not advance temperature, enthalpy, composition,
phase fraction, a physical acoustic equation of state or moving geometry. It
neither adds viscous/pressure-work heat nor closes a total-energy balance.
Density/Cp bounds are model-admission limits, not experimental qualification.

## Focused checks

The eight force-owner and seven material-driver Rust tests cover original-path
parity, independent push streaming, source impulse/reversal, worker-count
agreement, atomic cancellation/retry, source sampling, actual gas-temperature
substitution, epoch/domain/speed refusals and manufactured buoyancy-driven flow.
The pre-existing ten collision-field regressions remain unchanged.

The Stokes fixture retains its original error bounds. The high-viscosity
16-cell control misses them (12.244% relative velocity error in the independent
NumPy reference), so that arm uses 32/64 cells rather than a relaxed gate.
The reference errors are 1.635% -> 0.411% for mu=0.01 Pa*s at 16/32 cells and
3.059% -> 0.765% for mu=0.04 Pa*s at 32/64 cells. The expected second-order
interval is 1.8..=2.2. These are independent numerical-reference results, not
native Rust execution or experimental validation.

Native verification remains required; Cargo/rustc were unavailable during
implementation. On the project's configured build lane, the focused commands
are `cargo test -p fs-lbm --lib lattice3::sources -- --nocapture` and
`cargo test -p fs-lbm --lib lattice3::collision_field_tests`, followed by the
existing D3Q19/D3Q27 integration batteries. No test is claimed executed merely
because its source has been committed.
