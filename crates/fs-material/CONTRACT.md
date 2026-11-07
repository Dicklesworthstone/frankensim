

## Analytic frequency dispersion of original Maxwell spectra

`GeneralizedMaxwell::modulus_frequency_derivative` returns the original
storage and loss slopes with respect to angular frequency [Pa s]. Scalar
frequency retains signed parity: storage is even, loss is odd, storage slope
is odd and loss slope is even. `TensorProny::complex_modulus_frequency_derivative`
returns those full Mandel tensor slopes plus `d(C''/omega)/domega` [Pa s²],
including signed anisotropic coupling entries. The original nonnegative
applicability band is enforced, including its endpoints; the algebraic law
supplies the admitted one-sided endpoint slopes. At admitted DC, storage and
effective-viscosity slopes are zero and loss slope is `C_kv + sum(C_j tau_j)`.
Kelvin viscosity contributes its original constant loss slope and exactly
zero effective-viscosity slope.

The derivatives use the authored branch stiffnesses and relaxation times,
with binary-scaled products and compensated accumulation in branch order.
No spectral refit, perturbation solve or frequency differencing is performed.
Reciprocal evaluation and delayed physical scaling retain representable weak
storage, loss and derivative entries when `omega*tau` or its square cannot
be represented. Nonfinite inputs or required unrepresentable derivatives
return `ViscoError::Parameters`; tensor band violations remain
`ViscoError::OutOfBand`. A finite primal query does not require a finite
tangent, and the scalar modulus slope does not require the additional
effective-viscosity slope to fit. The small fixed tensor queries allocate no
workspace and introduce no worker-dependent reduction order.

`tests/visco_frequency.rs` contains ten synthetic G0/G1/G3 regressions for
closed-form branch slopes, full anisotropy and Kelvin terms, independent
native-primal finite differences, signed scalar parity, tensor band admission,
DC and representability boundaries, weak mutual entries over extreme scales,
subnormal relaxation times and cancellation between relaxation slopes.
These integration tests are authored but have not run in the implementation
environment: the real dependency graph is blocked while compiling
`asupersync`. Six inline tests of the actual shared spectrum kernel passed
with the pinned native compiler. Neither result establishes measured-material
qualification or replaces the original hereditary states during time
integration; `C''/omega` here is a harmonic constitutive quantity.
