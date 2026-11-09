

## Causal pressure-driven molecular compressibility

`RelaxationState::compressibility_relaxation` retains the original evaluated
species/reservoir state and lowers its complete `Cp(s)` to a positive finite
pressure-memory law `beta(s) = beta_inf + sum delta_beta/(1+s*tau)`.
The bounded original symmetric Jacobi owner diagonalizes the rank-one inverse
of `Cp_inf + sum C_j/(1+s*tau_j)`. These pressure-driven times differ from the
caloric and compression-driven times. Equal source times combine exactly;
zero-strength reservoirs remain in the retained source but allocate no arm.
No frequency fit, attenuation conversion or extra bulk viscosity is introduced.

Each arm obeys `tau*q_dot + q = p`, stores `delta_beta*q^2/2` per volume and
dissipates `delta_beta*(p-q)^2/tau`. A spatial pressure consumer must use the
frozen compressibility for its instantaneous mass; the memories supply the
relaxed increment once. Original gas coherence, finite positivity, eigenpairs
and relaxed/frozen limits are checked before accepting a realization.
Harmonic response uses the existing peak `exp(-i*omega*t)` convention.

Four focused regressions compare the independent one-reservoir pressure time
and strength, multi-reservoir response against the reciprocal original caloric
modulus over frequencies from zero to `f64::MAX`, equal/zero/empty reservoirs,
physical pressure scaling and original-state admission. They ran in a bounded
native build of unchanged original gas/caloric/Jacobi numerical owners; the
largest sampled relative compressibility discrepancy was `2.144e-16`.
This is a frozen-state, linear molecular constitutive realization; it supplies
no shear viscosity, Fourier heat conduction, material calibration or full
workspace-build claim.
