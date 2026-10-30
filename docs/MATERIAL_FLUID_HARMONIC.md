# Forced vibration of source-defined frames in regional fluids

`MaterialFluidFrame::prepare_harmonic` prepares a steady forced-response solve
at one positive angular frequency. It borrows the exact original frame/fluid
owner, including all solid material ports, loaded-region identities, spatial
traces, added inertia and spectral memory. It does not create new material
coefficients, retune modes, or alter an accepted time-domain state.

Peak complex amplitudes use `exp(+i omega t)`. With the original structural
velocity v and a station trace b, each fluid branch obeys

```
v_j = lambda_j/(lambda_j + i*omega) * b^T v
slip_j = i*omega/(lambda_j + i*omega) * b^T v
```

The solve eliminates these independent equations exactly, factors only the
structural dynamic stiffness, and then reconstructs every branch momentum.
The full loaded mass is included once. Only the low-tail and memory resistance
are added afterward; potential-flow inertia is not counted again. Small slip
is calculated directly, not by subtracting nearly equal fluid/body velocities.

```rust,ignore
use fs_couple::air_path::cylinder::frame::harmonic::FrameHarmonicBudget;

let response = model.prepare_harmonic(
    omega_rad_s,
    FrameHarmonicBudget {
        max_modal_entries: 4096,
        max_setup_terms: 2_000_000,
        max_response_values: 100_000,
        equation_relative_tolerance: 1e-10,
        power_absolute_tolerance_w: 1e-10,
        power_relative_tolerance: 1e-8,
    },
    &gate,
)?;
let first = response.solve_modal(&first_generalized_force_phasor, &gate)?;
let second = response.solve_modal(&second_generalized_force_phasor, &gate)?;
// Both loads reuse the SAME factorization at this frequency.
```

Each result includes original-basis displacement/velocity, the complete state
phasor, cycle-mean stored energy, separate solid-port and regional-fluid loss
powers, input power and equation/power residuals. Forces may have independent
phases. Different regions of one element remain distinct accounts. Heating is
not automatically deposited into caloric cells by a steady harmonic query.

Preparation work is bounded by `m^3 + m^2*(1+ports+stations) + memories`, with
separate matrix and complete-result limits before allocation. The existing
cancellable `fs-la` complex LU is reused; no new linear factorization or time
integrator is introduced. Source viscosity frequency support is checked at the
actual angular frequency. DC, singular factors, incomplete inputs, nonfinite
arithmetic, exceeded budgets, cancellation and failed power/residual gates
refuse without changing either the source model or previous results.

The original model still owns its dense preparation/state budget: this path
reduces the harmonic factorization, not that earlier allocation. Native elapsed
time, total byte use, condition numbers and real-time speed are not certified.
A small residual does not bound forward error near resonance. Modal, spatial,
spectral-cutoff and physical-model refinement remain independent obligations.

This is frozen small-amplitude Newtonian local-cylinder forced vibration,
not transient initialization, live thermal constitutive replacement, background
flow, free-surface interaction, acoustic radiation or measured-material validation.

## Qualification

Six native regression groups compare the reduced response with the original
full `i*omega*I - (J-R)*Q` pencil, check independent stored energy/loss, analytic
dry response, complex amplitude scaling, uniform-region partition invariance,
data-only source swaps, source bands and cancellation/retry.

```
rch exec -- cargo test -p fs-couple --lib air_path::cylinder::frame::harmonic
```

Rust/Cargo/DSR/RCH are unavailable in the editing environment. Native compilation,
tests and formatting are UNEXECUTED. Independent NumPy/SciPy manufactured
Euler-Bernoulli/Stokes equations passed 30 full-versus-reduced cases: maximum
normalized state difference 1.91e-15, scaled full-equation residual 5.41e-17,
power discrepancy 1.14e-17 W and mean-energy difference 2.85e-16 J. A four-mode,
two-region reference factors a 4-by-4 rather than 40-by-40 matrix. This is an
algebraic dimension comparison, not a native performance benchmark.

## Physical force, moment and member-point ports

`harmonic_ports::FrameHarmonicPorts` borrows the original model and projects
physical nodal force/moment columns once for an entire frequency sweep.
`harmonic_transverse_ports` derives columns directly from the original member's
cubic transverse trace at `(element, xi)`. The corresponding nodal moments are
retained; a point inside a member is not replaced by endpoint-only forcing.

```rust,ignore
let ports = model.harmonic_transverse_ports(
    &[(0, 0.3), (1, 0.6)], port_budget, &gate,
)?;
for omega in frequencies_rad_s {
    let prepared = model.prepare_harmonic(omega, harmonic_budget, &gate)?;
    let mobility = ports.mobility(&prepared, &gate)?;
    let driven = ports.solve(&prepared, &complex_drive_amplitudes, &gate)?;
    // driven retains actual nodal displacement/velocity, source losses and
    // signed work at each drive. No physical clock or accepted history moves.
}
```

`FrameHarmonicPortBudget` bounds port count, `ports*modes*nodal_dofs` projection
work, the complete `ports*ports` matrix, and independent nodal/port power error.
General columns use original node-major x/y/rotation order. Caller-declared
units must make each column times its effort a nodal force or moment; the dual
motion observation preserves that work. No column is normalized or made positive.
Zero, repeated and dependent columns remain valid; singular port mobility is not
repaired. Source identity is checked: another same-sized model cannot reuse a
port map accidentally.

The returned matrix is row-major `velocity_i = sum_j Y_ij effort_j`, with one
shared factorization for all columns at that frequency. Reciprocity is not
forced by averaging entries. No time-domain realization, measured accuracy or
unconditional forward-error bound is implied by a matrix or residual report.

Coherent drive phases are combined BEFORE solving motion and computing heat.
Adding separate unit-drive heat accounts would incorrectly omit interference.
Individual drive work can be negative while the complete passive system absorbs
positive net power. The API retains those signs and independently checks nodal,
port and modal work. Returned reports and triangular solves allocate; no
allocation-free solve or native speedup is claimed.

Five additional Rust regressions cover complex mobility/reciprocity, coherent
cancellation/reinforcement, original point geometry and rotation work, signed
coordinate scaling, source mismatch and cancellation/retry. They are authored
but UNEXECUTED under the toolchain limitation above. Independent NumPy/SciPy
physical force/moment tests at 1, 30 and 100 rad/s passed: reciprocity discrepancy
below 7.11e-15, nodal work mismatch below 4.34e-19 W, exact cancellation for
opposing duplicate drives, and fourfold power for reinforcing duplicate drives.
These are synthetic equation checks, not native Rust or experimental evidence.
