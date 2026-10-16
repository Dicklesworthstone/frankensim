# Supplied hereditary plate bending (material memory)

The retained plate can now carry regional Generalized Maxwell bending laws in
its nonlinear valve/contact/air solve. This adds **stress relaxation and a
history-dependent restoring force**, not an output damping envelope. Two valves
with identical opening and velocity but different viscous history can respond
differently to the same pressure.

Run the explicitly illustrative supplied spectrum with spatial lay closure:

```bash
cargo run --release -p fs-couple --example plate_aperture -- \
  --profiled \
  --relaxation crates/fs-couple/examples/plate_aperture/illustrative-relaxation.txt \
  4000000000 900 0.0003 > /tmp/plate-memory.csv
```

The spectrum file supplies equilibrium Young modulus [Pa], constant Poisson
ratio, a declared material-use band [Hz], temporal resolution allowances,
initial history, and `(additional modulus [Pa], relaxation time [s])` pairs.
Records are complete and ordered as in the example; unknown/trailing records
refuse. The file is bounded at 64 KiB and 64 supplied terms. The supplied
illustrative law has E_inf=4 GPa and additional arms (2 GPa, 1 ms) and
(1 GPa, 10 ms). It is **not measured wet-cane or polymer data**. Modulus,
thickness and density arguments still define the actual equilibrium specimen;
the spectrum is refused if it does not match that specimen.

Selecting `--relaxation` replaces the example's illustrative 0.02 viscous damping
ratio with the supplied hereditary law, and prints that choice. The library
requires zero separate modal damping rather than silently counting both losses.
Unselected runs retain the old physical solver and parameters. Both profiled
and uniform slit geometry are supported. Pressure remains the same 5 Pa program;
there is no change-dependent normalization, new sound source or gain.

For regional models, call
`DynamicAperture::with_plate_relaxation(PlateRelaxationSpec, InitialApertureMemory)`
after either plate constructor and **before** any accepted step or tube/network
construction. The regions must partition every source triangle exactly once.
Each region supplies its own isotropic `GeneralizedMaxwell`, Poisson ratio,
material-use band and source/condition label. `E_inf` and Poisson ratio must
reproduce each covered triangle's actual equilibrium DKT bending matrix;
orthotropic/nonproportional viscoelastic tensors are not inferred from a scalar
E(omega). Original plate geometry, material receipts and lay are retained.

For the retained opening-coordinate shape psi, each region's existing element
matrices produce `K_region = sum psi_e^T K_e psi_e`. An arm's stiffness is
`K_region E_j/E_inf`. This excludes membrane prestress and does not replace or
double-count the equilibrium stiffness. Regional thickness/material variation
therefore changes the actual restoring force. One positive arm per regional
term suffices for this one structural coordinate and a uniform equivalent
initial viscous displacement within that region; arbitrary spatial initial
memory and additional structural modes are not represented.

Runtime arms use the existing fs-phs internal-memory storage and resistance:
`H_arm=(sqrt(k) x-z)^2/2`, `z_dot=(sqrt(k) x-z)/tau`, with x measured from the
original plate rest geometry. The linear z equation is eliminated at the same
midpoint as the nonlinear pressure/contact/momentum solve. There is no lagged
stress correction or replacement integrator. All candidate z values remain
private until the complete valve/tube/network step is accepted. A failed
contact, slope, penetration, propagation or finite-value check consumes no
material history; cancellation and budget extension retain it.

`InitialApertureMemory` requires relaxed, unrelaxed, or explicit per-arm viscous
displacement [m]. Relaxed means zero initial arm stress at the current bent
shape, not zero plate equilibrium energy. Unrelaxed means zero prior viscous
displacement, not zero initial arm energy. Initial memory is not a way to reset
an already vibrating instrument; repeated or late attachment refuses.

Admission requires caller-selected dt/tau and dt*omega_instantaneous bounds.
The former is at most 2, excluding midpoint's sign-alternating held-strain
memory regime; choose a smaller allowance for accuracy. The latter is at most
1 radian. Retained equilibrium and instantaneous frequencies must fit every
supplied material-use band. These are resolution/domain checks, **not** proof
that an arbitrary nonlinear transient has no out-of-band content. Midpoint
memory is not exact exponential relaxation, and the supplied spectrum is not
automatically validated against experiment.

`relaxation()` exposes immutable source laws, projected regional stiffness,
accepted energy-normalized history and material energy/loss. The CSV adds
`material_energy_j` and `material_step_loss_j`; both are already included in
its total energy/loss columns and must not be added again. Contact and jet loss
remain separate physical mechanisms. Memory state/energy validation allocates
bounded candidate storage per step; no allocation-free or real-time claim is
made. The one-mode, linear-geometry, fixed-lay and internal-pressure limitations
in `PLATE_APERTURES.md` still apply.
