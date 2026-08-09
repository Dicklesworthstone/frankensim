# fs-nlmodal — Contract

## Purpose and layer

Layer L3. Geometric nonlinearity for thin structures in modal
coordinates (bead frankensim-fsim-vonkarman-u5nme): von Karman plates
and Kirchhoff-Carrier strings as port-Hamiltonian systems with a
quartic SUM-OF-SQUARES potential
`H = 1/2 sum (p^2 + w^2 q^2) + 1/4 sum_j c_j (q^T E_j q)^2` — the
Airy-stress-eliminated modal form. Pitch glide, mode coupling, and
the crash cascade with structurally exact energy accounting.

## Public types and semantics

- `SosModalStorage` (implements fs-phs `Storage`): frequencies +
  `StressChannel`s (`c_j >= 0`, symmetric `E_j`). Gradient is the
  EXACT analytic gradient of the coded `H`; state layout matches
  fs-phs `modal_bank` (`[q, p]` interleaved).
- `assemble(storage, zetas, strike_weights)` — the pHS: symplectic
  pair blocks, per-mode viscous damping `R = diag(0, 2 zeta w)`
  (zetas come from the caller; the visco-damping facility's per-mode
  output slots in here — no second damping representation invented),
  one strike port. REFUSES asymmetric couplings, negative
  coefficients/damping, non-positive/non-finite frequencies, and
  non-finite couplings by name (NaN-proof negated comparisons);
  duplicate modes in either constructor list are refused.
- `von_karman_ss_plate` — simply-supported rectangle with ANALYTIC
  sine modes for displacement (mass-normalized) and Airy stress
  (unit-normalized; both are biharmonic eigenfunctions on the SS
  rectangle, in-plane movable edges). Channel coefficients
  `E h / (2 xi_j^4)`. Coupling integrals by Gauss-Legendre quadrature
  whose order SCALES with the highest half-wave sum (a fixed order
  left ~5 points per wave — executed), certified by a second
  independent order judged against max(channel scale, 1e-12 x global
  scale) — entrywise relative comparison falsely refuses
  analytically-zero entries, and an ALL-zero channel's own scale is
  pure roundoff (both executed); the residual is returned, not
  hidden. The two orders share the quadrature engine, so the
  certificate is a CONVERGENCE check, not an independent-route one. Stress-mode count is a SEPARATE truncation from the
  displacement count (both explicit inputs).
- `kirchhoff_carrier_string` — one diagonal channel
  `E[k,k] = (k pi/L)^2 * 2/(mu L)`, coefficient `E A L / 8`: exactly
  the averaged-tension Kirchhoff-Carrier form in mass-normalized
  coordinates (hand-derivation pinned in the battery).
- `duffing_backbone` / `single_mode_beta` — the analytic
  perturbation pins.
- Time stepping, striking, damping, and the energy ledger are
  fs-phs's (Gonzalez discrete gradients — implicit, energy-exact,
  stable at crash amplitudes). No explicit integrator exists in this
  crate to misuse; the "explicit scheme refusal" is by construction.

## Invariants

1. Quadrature certificate: two independent Gauss orders agree to
   1e-8 of the channel scale, or construction refuses
   (`QuadratureMismatch`).
2. Coupling symmetry is an admission requirement
   (`AsymmetricCoupling`), and the tensor is non-vacuous for the SS
   plate.
3. Duffing backbone: the measured amplitude-dependent frequency of a
   single-mode system tracks `w0 (1 + 3 beta a^2 / (8 w0^2))` within
   10% of the SHIFT at small amplitude.
4. Kirchhoff-Carrier exactness: `single_mode_beta` equals the
   hand-derived `E A (k pi/L)^4 / (2 mu^2 L)` to 1e-12, and the
   backbone reproduces the classical glide
   `dw/w = (3/32)(EA/T0)(k pi A/L)^2` in physical amplitude.
5. Energy conservation: undamped struck-level free run holds `H` to
   1e-9 relative over 20k steps.
6. ARCHITECTURAL FINDING (executed, recorded): under Gonzalez
   discrete-gradient stepping the correction term forces
   `dg.(x1-x0) = H(x1)-H(x0)` for ANY gradient function, so
   force-vs-energy divergence (the unsymmetrized-tensor bug class)
   CANNOT break energy conservation — unlike force-side integrators.
   The guards here are admission (2) and the TRAJECTORY: a divergent
   mutant conserves its coded `H` to ~1e-7 while departing the true
   dynamics by > 5% — both pinned.
7. Amplitude-scaling covariance: the nonlinear frequency shift
   scales as amplitude^2 (measured ratio 4 +/- 10%).
8. Mode coupling cross-check: the two-mode tension-modulated string
   transfers energy into a seeded second mode (weak, second-order
   parametric — correct string physics) and the discrete-gradient
   trajectory agrees with an independent in-test RK4 integration of
   the same vector field to 2%.
9. Cascade onset: a fundamental-only initial condition leaks energy
   into other modes only through the nonlinearity; the leaked
   fraction grows strongly with `w/h` (0.1 -> 0.6 -> 2.0), ledger
   closed throughout.

## Error model

Typed `NlModalError` (parameter, asymmetric coupling, quadrature
mismatch, pHS passthrough). No silent degradation.

## Determinism class

Deterministic: Newton-on-recurrence Gauss nodes, fixed orders, fs-phs
stepping; no RNG or time dependence.

## Cancellation behavior

Construction is milliseconds at analytic-mode scale; stepping is the
long pole (implicit Newton per step). No `Cx` integration (workspace
`frankensim-ccmn`).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/nlmodal.rs` (10): construction certificates (quadrature
residual, mode counts, assembly); non-vacuous coupling + asymmetry
refusal; Duffing backbone vs perturbation formula; Kirchhoff-Carrier
beta + glide hand formulas; energy conservation; the
unsymmetrized-tensor architectural falsifier (conservation survives,
trajectory diverges); amplitude-scaling covariance; parametric
coupling + RK4 cross-check; plate cascade casebook (JSON lines with
modal energies, leaked fractions, tensor size); typed refusals.

## No-claim boundaries

- Moderate-rotation von Karman regime only: no damage, plasticity,
  wrinkling, or full shell curvature (quadratic coupling from
  curvature is the shells follow-up).
- Simply-supported in-plane-movable boundary in v1; FE-mode input is
  the follow-up (trigger: consumers with non-analytic geometry —
  which is also when coupling-tensor CACHING starts to matter; at
  analytic-mode scale the tensor build is milliseconds and a cache
  manifest would be ceremony, deliberately not built).
- The Airy potential form `V = (Eh/8) sum S_j^2/xi_j^4` is taken
  from the standard modal literature (Ducceschi-Touze class); its
  validation here is through the Duffing/KC analytic pins and the
  structural energy identities, not a from-scratch continuum
  re-derivation.
- Literature-value comparisons for internal-resonance energy
  exchange are deferred (measured tables need verified
  transcription); the executed cross-check is an independent
  INTEGRATOR over the same coded vector field — it validates the
  stepper, not the tensor physics (review-corrected wording).
- Damping is per-mode viscous `zeta_k` supplied by the caller;
  frequency-dependent loss models remain the visco-damping
  facility's.
