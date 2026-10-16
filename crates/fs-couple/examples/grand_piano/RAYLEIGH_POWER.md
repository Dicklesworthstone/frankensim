# Match soundboard radiation loss to outgoing pressure power

The piano's explicit `--rayleigh-loading` preparation uses
`fs_bem::helmholtz::rayleigh_modal::PreparedSurface::power_consistent_impedance`.
Both the full-load and `--radiation-ports` paths sample every admitted loaded
surface field through this operator before the existing passive fitter. The
old `PreparedSurface::impedance` remains an unchanged comparison image; the
scalar baffled-piston API also remains unchanged.

No new CLI control or output EQ is needed. Existing renders that omit
`--rayleigh-loading` are unaffected. One/two half-space selection, physical
medium, source geometry, string-loaded board basis, fitting thresholds,
mechanical clock, damping, and microphone pressure scale are unchanged.
The report identifies the resistance quadrature and retained self reactance.

## Why the self term matters

For area-point source strengths `b_i = area_i * normal_velocity_i`, the real
part of the Rayleigh loading kernel is smooth, even when the two points meet:

```text
K_R(i,j) = rho * omega^2 / (2*pi*c) * sinc(k * distance(i,j))
sinc(0)  = 1
R_modal  = (A Phi)^T K_R (A Phi)
P_mean   = 1/2 * real(v* R_modal v)        (peak complex velocity amplitudes)
```

The same quadrature points generate the outgoing pressure field. Integrating
its squared magnitude over a far hemisphere gives precisely this mean power:
the angular integral of `exp(i*k*direction.(x_i-x_j))` is
`2*pi*sinc(k*distance(i,j))` for coplanar points. Thus the resistance matrix is
a Gram matrix and is positive semidefinite in exact arithmetic for arbitrary
signed fields and complex modal phases. Signed off-diagonal entries can still
be negative; zero fields still have exactly zero loss. No diagonal floor,
eigenvalue clipping, fitted damping addition or receiver-gain change is used.

The legacy image integrates only the self-cell real part over an equivalent
disc while evaluating the other real entries at points. Its diagonal is
smaller than the regular-kernel limit. Subtracting that diagonal correction
from the Gram matrix can turn a weakly radiating cancelling field into an
acoustic energy source. Wavelength resolution alone does not prevent this.

An explicit three-cell counterexample uses centres at -20, 0 and +20 mm,
100 mm^2 per cell, and relative velocities 1, -2, 1. Equivalent discs do not
overlap. At 400 Hz in rho=1.2 kg/m^3, c=343 m/s, the legacy rule gives
-1.34080123e-8 W; consistent quadrature gives +1.61800384e-9 W. An independent
positive angular integral of `16*A^2*sin(k*d*direction_x/2)^4` gives the latter
value. This is a manufactured physics regression, not a measured piano result.

## What is and is not changed

The singular imaginary self term still uses the existing equivalent-disc
integral. All other imaginary entries, work budgets, cell admission and full
reciprocal projection are unchanged. The resistance consistency result does
not certify that this mixed spatial approximation is an exact finite-cell
Galerkin discretization, converged on a particular mesh, or a causal
positive-real impedance over the entire frequency axis. The existing passive
realization and complex/resistance fit gates remain mandatory and can still
refuse a requested load. No test failure or fit threshold is suppressed.

Very weak radiation directions can still have ordinary floating-point
nullspace error. The implementation does not relabel that roundoff as extra
material loss. This correction does not resolve missing soundboard modes,
transverse shear, measured wood properties, rim/lid/room effects or perceptual
calibration. The high-band piano fidelity task remains open.

## Verification

```sh
cargo test --release -p fs-bem --lib rayleigh_modal
cargo test --release -p fs-couple --example grand_piano rayleigh_loading
```

Five new kernel regressions cover the negative-power counterexample, an
independent outgoing-power integral with 3 and 128 fields, reciprocal signed
projection, unchanged imaginary entries, physical scaling and admission.
Two piano regressions exercise the actual pre-fit sweep, both half-spaces,
loaded shapes, the unchanged passive fitter and explicit acoustic projection.
Existing piston, radiation, and combined string/board/air energy tests remain.

Independent Python/NumPy/SciPy references were executed for three cancelling
cases, six complete-field angular-power comparisons and 96 randomized surfaces.
The maximum randomized relative power discrepancy was 4.56e-15; the smallest
relative resistance eigenvalue was -1.54e-16 (roundoff). The legacy matrices
had materially negative eigenvalues in 49 of those randomized cases. These
are separate numerical references, not execution of the Rust tests. Native
compilation, fitted-load rendering and listening remain unverified here.
