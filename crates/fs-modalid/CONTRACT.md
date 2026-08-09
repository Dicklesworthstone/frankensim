# fs-modalid — Contract

## Purpose and layer

Layer L3. Experimental modal identification (bead
frankensim-fsim-modal-id-1qxdh): measured FRFs -> modal parameters
`(f_k, eta_k, phi_k)` with machine-readable quality gates — the
measurement side of the simulate-vs-measure calibration loop.

## Public types and semantics

- `FrfData` — multi-channel FRF on a shared grid; `new` CONJUGATES
  engineering (`e^{-i*omega*t}`) data onto the Laplace axis (the
  load-bearing convention bridge from the vector-fitting bead);
  `new_laplace` skips it; `parse_csv` for the v1 header-free CSV
  table (freq, re, im[, coherence] per channel; comment lines `#`).
  No UFF/proprietary formats in v1 (stated).
- `estimate_snr` — coherence-driven `gamma^2/(1-gamma^2)` when
  present; otherwise a SECOND-DIFFERENCE roughness estimator (white
  noise amplifies 6x in `d2` power; a smooth FRF contributes only
  curvature). The naive lowest-decile magnitude floor is WRONG for
  FRFs (antiresonance valleys are signal — executed failure pinned in
  the doc).
- `identify` — the pipeline: SNR gate (typed `SnrTooLow` refusal) ->
  vector-fit stabilization ladder (per-order pole tables) ->
  `stabilization` verdicts (frequency/damping match tolerances,
  required consecutive runs; every final-order pole gets an
  accept/reject + reason) -> near-duplicate merge at 1/5 the match
  tolerance (the two tolerance roles differ — executed: merging at
  the match tolerance swallowed a genuine 0.5%-separated pair) ->
  per-channel residues at the SHARED accepted poles (fs-vfit
  `residue_fit_at_poles`, conjugate-EXPANDED pole list) -> residue
  significance gate (a stable spare pole modeling the noise floor
  carries residues orders below physical modes) -> split-sample
  (even/odd half-grid) confidence intervals -> optional
  exponential-window damping correction with raw values logged.
- `rfp_fit` — classical rational-fraction-polynomial identifier on a
  Forsythe-orthogonal basis (three-term recurrence, orthonormal under
  the fit weights), denominator roots via the companion matrix,
  residues through the SAME fs-vfit residue pass. The third
  independent identifier (vector fitting and Loewner re-exported from
  fs-vfit are the other two); disagreement is a diagnostic.
- `mac` / `mac_matrix` / `mac_pairing` — modal assurance criterion
  and greedy floor-gated pairing: the calibration diff primitive.
- `correct_exponential_window` — `zeta_true = zeta_meas -
  1/(tau*omega_n)`, clamped at zero, delta reported.

## Invariants

1. Known-answer: 10 synthetic modes recovered across noise levels
   1e-6..1e-3 within the reported split-sample intervals plus an
   authored noise-scaled floor; shapes MAC > 0.99 against truth.
2. Close-mode resolution: a 0.5%-separated pair resolves into exactly
   two modes where an order-2 single fit cannot hold both; the RFP
   cross-check finds the same pair.
3. Stabilization verdicts are machine-readable (frequency, damping,
   run length, accepted, reason) and the diagram (per-order pole
   tables) ships in the result.
4. Window-correction mutation is visible: disabling the correction on
   exponentially windowed data biases damping by 10x the corrected
   error against the known answer.
5. SNR refusal fires by name below the floor; no fabricated modes.
   Zero surviving poles is the typed `NoStablePoles`, not an empty
   table.
6. RFP conditioning: the monomial-basis Gram matrix at degree 20 is
   numerically singular (condition > 1e8 asserted) while the Forsythe
   basis is orthonormal by construction — the reason the orthogonal
   basis is load-bearing, demonstrated not asserted.
7. Published benchmark: identification recovers Carcagno et al. JASA
   144(6):3533 Table I modal parameters (Brazilian-rosewood guitar,
   F/Q verified verbatim from the saved CC-BY PDF) from an FRF
   SYNTHESIZED from those published values under measurement-grade
   noise — parameters are published, the trace is not (mobility
   magnitudes are figures-only across the CC-BY literature; honest
   label).

## Error model

Typed `ModalIdError` (shape, SNR refusal, fit passthrough, CSV line,
no stable poles). No silent degradation.

## Determinism class

Deterministic: fs-vfit/fs-la kernels, no RNG (test noise is
golden-angle-stride pseudo-noise), fixed ladder orders.

## Cancellation behavior

Synchronous; the ladder is the long pole (seconds-class at order 28
on 1600 samples). No `Cx` integration (workspace `frankensim-ccmn`).

## Unsafe boundary

`#![deny(unsafe_code)]` via workspace lints; no capsules.

## Feature flags

None.

## Conformance tests

`tests/modalid.rs` (11): known-answer noise sweep inside intervals;
MAC shape recovery; close-mode pair + RFP cross-check + low-order
failure demonstration; MAC identity/orthogonality/scale-phase
invariance; perturbed-plate calibration diff (intact modes pair
> 0.99, perturbed flagged); window-correction mutation visibility +
formula unit test; SNR refusal + coherence-driven estimates; RFP
Gram-conditioning demonstration; CSV ingest round-trip with typed
line refusal; Carcagno published-parameter benchmark.

## No-claim boundaries

- MIMO identification (matrix-valued Loewner over simultaneous
  references) is a follow-up; v1 is common-pole multi-channel from a
  reference channel.
- No measured-trace benchmark: no license-compatible published FRF
  TRACE was available (figures-only across the CC-BY literature);
  the published-parameter round-trip is the honest v1 benchmark, and
  ingesting a real measured dataset is the stated trigger for the
  follow-up.
- Uncertainty is split-sample spread, not a statistical confidence
  interval with coverage guarantees (recorded; a bootstrap/posterior
  treatment is a follow-up).
- Complex-mode-indicator function (CMIF) for repeated roots at ONE
  frequency is a follow-up; the close-mode battery covers separated
  close pairs.
- UFF/universal-file ingestion is out of scope (stated in the bead's
  polish round).
