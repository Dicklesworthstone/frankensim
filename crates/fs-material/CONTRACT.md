

## Complex plane-stress harmonic frequency derivatives

`TensorProny::plane_stress_complex_frequency_derivative` differentiates the
original complex transverse-stress constraint in all four ordered material
planes. It returns signed storage/loss slopes in the same engineering-shear
coordinates as `plane_stress_complex`. The primal material evaluation,
transverse solve and reduced passivity admission are shared once; the public
primal query retains its original pair result and admission semantics.

For `X=C_oo^-1 C_oi` and `T=[I;-X]`, the exact complex-symmetric Schur derivative
is `Q_omega=T^T C_omega T`, with ordinary transpose and no conjugation. This
retains the changing transverse-strain coupling between storage and loss
without a second factorization or a derivative solve. The shared material
aggregation requests only modulus slopes: an unrepresentable derivative of
effective viscosity does not invalidate a finite plane-stress slope.

Each reduced entry has its own factored, compensated complex contraction.
A large eliminated normal derivative therefore cannot globally rescale away
an unrelated weak in-plane derivative. Physical units are restored after
signed accumulation; finite-precision cancellation is not a relative-error
certificate. The result preserves exact symmetry and the original frequency
band. Derivatives need not themselves be passive or positive matrices.

Six integration regressions cover independent full stress constraints with
noncommuting anisotropic branches in every plane, native frequency differences,
complex isotropic formulas, independent viscosity-slope admission, an in-plane
`1e-100` slope beside an eliminated `1e300` slope, and typed source refusals.
They remain unexecuted with the full dependency graph. Four native tests of the
actual production projection module passed against the real `fs-math` library,
including complex product signs, overflowed products that cancel, weak retained
responses and minimum-subnormal coefficient scaling. This is a material query
for plate/shell consumers; it does not add an assembled plate/shell frequency
solver or replace the original hereditary material history.
