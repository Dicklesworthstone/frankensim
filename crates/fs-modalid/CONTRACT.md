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
  present; otherwise a FOURTH-DIFFERENCE roughness estimator (white
  noise passes with power gain 70; the median sits at ln2 of the
  mean — both divided out). The naive lowest-decile magnitude floor
  is WRONG for FRFs (antiresonance valleys are signal — executed
  failure). Channels under 5 samples vote ZERO (refusal), never
  infinity.
- `identify` — the pipeline: SNR gate (typed `SnrTooLow` refusal) ->
  vector-fit stabilization ladder (per-order pole tables) ->
  `stabilization` verdicts (frequency/damping match tolerances,
  required consecutive runs; every final-order pole gets an
  accept/reject + reason) -> near-duplicate merge at 1/5 the match
  tolerance (the two tolerance roles differ — executed: merging at
  the match tolerance swallowed a genuine 0.5%-separated pair;
  within a cluster the LONGEST-stable-run pole wins — keeping the
  first kept a spare 0.1 Hz off the true mode, executed) ->
  per-channel residues at the SHARED accepted poles (fs-vfit
  `residue_fit_at_poles`, conjugate-EXPANDED pole list) -> residue
  significance gate (dropped frequencies RECORDED in
  `insignificant_freqs_hz`, never silent) -> split-sample (even/odd
  half-grid) confidence intervals with a capped nearest-match (an
  unmatched half-grid pole yields a NaN interval, never zero) ->
  optional exponential-window damping correction with raw values
  logged. Ladder options are validated (`order_step >= 1` etc.).
- `rfp_fit` — classical rational-fraction-polynomial identifier on a
  Forsythe-recurrence basis. HONEST MECHANISM (review-corrected): on
  a positive-frequencies-only grid the basis is NOT orthonormal in
  the complex inner product; its i^d PARITY structure makes the real
  Re/Im-stacked LS Gram near-identity (condition < 100 asserted,
  measured 32). Denominator roots via the companion matrix; residues
  through the SAME fs-vfit residue pass. The third
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
6. RFP conditioning: the monomial-basis Gram at degree 20 is
   numerically singular (> 1e8) while the real-stacked Forsythe Gram
   conditions below 100 — both DEMONSTRATED in the battery, with the
   honest mechanism statement above.
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
- The d4 SNR estimator assumes the grid resolves the peaks (several
  samples per half-width); a legitimate coarse survey grid can still
  refuse — conservative, but a stated limitation (coherence data
  sidesteps it).
- Stabilization verdicts cover FINAL-order poles only: a mode stable
  through intermediate orders but relocated at the top order gets no
  verdict (stated; longest-run-over-all-orders diagrams are a
  follow-up).
