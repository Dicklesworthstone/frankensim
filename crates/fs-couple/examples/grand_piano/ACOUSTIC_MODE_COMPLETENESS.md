# Complete soundboard modes in a wide acoustic band

The shared `fs-modal` shift-invert Lanczos kernel no longer treats `beta <=
1e-14` as an unconditional invariant-subspace breakdown. The shifted inverse
has units of inverse eigenvalue. In a 20 kHz soundboard window its eigenvalues
can cluster around 1e-10, with useful expansion directions below that absolute
cutoff. Restarting cannot reliably recover those directions if the same cutoff
rejects them again. Merely raising a mode cap does not repair this defect.

The stopping test now compares beta with 16 machine epsilons times the largest
observed tridiagonal coefficient in that sweep. There is no floor of one on the
scale, no increased iteration budget and no relaxation of Ritz convergence or
inertia-count admission. Exact breakdown still permits deflation restarts;
insufficient budgets still refuse rather than returning incomplete modes.
Both assembled soundboard preparation and the matrix-free modal API use this
same kernel. Nothing changes the physical stiffness, mass, material damping,
bridge projections, or radiation mapping.

The existing piano regression
`board_geometry::tests::bridge_ritz_preparation_preserves_sampled_transfer_and_one_physical_motion_basis`
previously refused with 15 of 36 inertia-counted modes in its 20 kHz source
window. An independent diagonal numerical reference with 36 known resonances
at 50..1500 Hz inside a 20 kHz window similarly recovered only 10 with the old
cutoff and all 36 with the new cutoff, within the same 36 iterations. This is
a solver reference, not a measured piano or a full soundboard test result.

Focused native checks:

```sh
cargo test -p fs-modal
cargo test -p fs-couple --example grand_piano bridge_ritz_preparation
```

New tests cover wide acoustic windows with nonuniform mass, matrix-free scaling
across 26 orders of magnitude, actual degeneracy, mass orthogonality, explicit
residuals and budget refusal. The existing native piano workflow also exercises
the real board, its reduced bridge transfer and both acoustic renderers.

Recovering all modes in the selected discrete model is necessary but does not
certify continuum mesh convergence, high-frequency bridge impedance, a complete
acoustic band at the microphone, or perceptual similarity to a real piano.
