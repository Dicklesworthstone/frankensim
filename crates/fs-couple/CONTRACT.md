

## Harmonic Rayleigh pressure and vector receivers

`pcm_wav::baffled::harmonic::Prepared` is the continuous-frequency image of the
existing `BaffledPressure` surface basis. It consumes the same positive physical
areas and complete signed `SurfaceSample::mode_shape` coefficients in the z=0
infinite rigid baffle. A 1/area basis still denotes volume flow. Prepared transfer
columns expose pressure, Cartesian pressure gradient and fluid particle velocity
for every input and receiver; repeated coherent source evaluations allocate no
memory. The shared `fs_bem::near_field::FirstOrder` produces ideal omni, cardioid
and figure-eight Pa-equivalent observations with the same resolved medium.

All phasors use `exp(-i omega t)`. For one volume-velocity element Q,
`p = -i rho omega Q exp(ikR)/(2 pi R)` and
`u = Q exp(ikR) (1/R^2 - i k/R) r_hat/(2 pi)`.
The gradient is `i rho omega u`. Retaining the velocity's reactive `1/R^2` term
gives `u_r/(p/(rho c)) = 1 + i/(kR)`, so proximity response and phase follow from
the field itself. Signed complex surface/input contributions are summed before
microphone projection; a scalar pressure null may retain transverse particle
velocity and directional output.

Preparation admits at most the original 128 input coordinates, the shared BEM
64 receivers and 120,000 surface points. Complete sample/input/receiver work is
checked against the declared budget (default 2^21, hard cap 2^24 contributions),
with cancellation before receiver rows and every 64 surface samples. The default
receiver clearance is the existing 5 cm; an explicitly smaller positive clearance
requires finer caller-supplied quadrature. Invalid/nonfinite inputs, unrepresentable
fields, clearance violations and incomplete arrays refuse. Source evaluation and
microphone projection leave output buffers unchanged on refusal.

Focused native tests compare the analytic baffled monopole and ideal cardioid
for kR=0.05..50, the far-field velocity limit, complex spatial finite differences,
a pressure-cancelled but velocity-live source pair, the original circular surface
quadrature against the analytic piston, and the original differentiating causal
pressure owner through 21.6 kHz at 48 kHz (relative RMS difference <5e-4 after its
reported numerical latency). Full 64-receiver/128-input order, work and cancellation
checks use the same production kernels and original `FirstOrder` implementation.

This provides exact Green kernels at the supplied quadrature points, not an error
bound for the continuous surface. Structural and spatial refinement remain the
caller's responsibility. The harmonic matrix has no sampling clock or temporal
interpolation filter. It adds no causal vector receiver, measured capsule response,
molecular absorption, finite-baffle/room scattering or reciprocal radiation load,
and does not close the wider M2.6-W1 work item.
