# Certificate-regime doctrine

Status: ratified schema v1 for
`frankensim-extreal-program-f85xj.9.1`.

This doctrine routes a requested claim to the kind of evidence object that can
honestly support it. It does not rank every rigorous object above every
statistical object. The proposition, horizon, domain, population, and intended
decision determine the route.

The table is descriptive governance data. Looking up a row does not mint
evidence, authenticate an artifact, prove a scientific statement, or admit a
runtime claim. The executable router is separately owned by
`frankensim-extreal-program-f85xj.9.3`.

## Closed claim-to-evidence table

The section between the code-derived markers is rendered from
`fs_govern::CERTIFICATE_REGIMES` and checked byte-for-byte by the crate's
integration tests. Edit the Rust table and this section together; neither is a
fallback authority for the other.

<!-- BEGIN CODE-DERIVED CERTIFICATE REGIME TABLE -->
| ID | Claim class | Required evidence object | Current capability map | Applicability and no-claim boundary |
| --- | --- | --- | --- | --- |
| `CR-01` | Root or event time | `interval-root-or-taylor-enclosure` — Interval root isolation or Taylor enclosure | `fs-ivl` / `interval-root-isolation` (available)<br>`fs-ivl` / `univariate-taylor-enclosure` (available) | Scope: one declared finite parameter/time domain with explicit function and derivative enclosure assumptions<br>No claim: does not establish behavior after the isolated event or outside the admitted box |
| `CR-02` | Short-horizon collision or reachability | `validated-reachability-tube` — Validated finite-horizon tube | `fs-ivl` / `validated-reachability-tube` (staged) | Scope: one declared finite horizon, dynamics version, initial set, disturbance set, and event geometry<br>No claim: validated reachability tubes are staged; ordinary interval propagation is not a tube certificate |
| `CR-03` | Conserved quantity | `discrete-balance-certificate` — Discrete conservation or balance certificate | `fs-evidence` / `discrete-balance-defect-receipt` (available) | Scope: one discrete operator, mesh/complex, boundary partition, source convention, units, and time slab<br>No claim: a balance receipt does not validate the constitutive model or measurement chain |
| `CR-04` | Local stability | `spectral-or-lyapunov-certificate` — Spectral or Lyapunov certificate in a stated domain | `fs-spectral` / `residual-enclosed-spectral-service` (available)<br>`fs-sos` / `quadratic-lyapunov-check` (available) | Scope: one stated equilibrium/orbit, model and linearization version, parameter domain, norm, and spectral or Lyapunov assumptions<br>No claim: local stability evidence does not imply global attraction, nonlinear robustness, or long-time predictive accuracy |
| `CR-05` | Long-horizon mean load | `statistical-observable-with-model-evidence` — Statistical observable with sampling and model evidence | `fs-eproc` / `time-uniform-mean-evidence` (available)<br>`fs-uq` / `model-uncertainty-context` (thin) | Scope: one declared observable, population/regime, sampling design, dependence model, stopping rule, and model-form evidence<br>No claim: sampling evidence does not prove the simulator model, and a trajectory enclosure is not required or implied |
| `CR-06` | Turbulent or broadband spectrum | `distributional-spectral-validation` — Distributional and spectral validation | `fs-spectral` / `spectral-computation-evidence` (thin)<br>`fs-uq` / `distributional-spectrum-validation` (staged) | Scope: one declared spectrum/statistic, frequency band, windowing convention, operating regime, observation process, and comparison metric<br>No claim: computed eigenvalue or FFT evidence alone does not validate a turbulent distribution or broadband field |
| `CR-07` | Reliability over duty cycles | `sequential-rare-event-statistics` — Sequential or rare-event statistical evidence | `fs-eproc` / `anytime-sequential-evidence` (available)<br>`fs-uq` / `rare-event-duty-cycle-model` (staged) | Scope: one declared duty-cycle population, failure event, censoring/dependence model, stopping policy, and model/field evidence<br>No claim: anytime validity under the stated sampling law does not validate the failure model or transfer to an undeclared population |
| `CR-08` | Exact long chaotic trajectory | `no-useful-bound` — NoUsefulBound | `fs-govern` / `certificate-regime-no-useful-bound-policy` (available) | Scope: a request for one exact trajectory materially beyond the admitted predictability horizon of the stated chaotic model<br>No claim: NoUsefulBound is an honest routing result, not proof that all statistical, local, event, or finite-horizon claims are impossible |
<!-- END CODE-DERIVED CERTIFICATE REGIME TABLE -->

## Intervals still have strong roles in chaotic systems

`NoUsefulBound` for an exact long chaotic trajectory is not a retreat from
certified arithmetic. Interval and Taylor evidence remain the appropriate
objects for bounded propositions:

- short-time event and root isolation, inside the declared finite search domain;
- local exclusions and collision-free boxes, inside the admitted geometry and
  finite horizon;
- parameter and constitutive-law bounds, for the stated law version and domain;
- discrete conservation and balance audits, for the exact operator, boundary,
  source convention, units, and time slab;
- ingredients of finite-horizon tubes, while the validated reachability-tube
  construction itself remains explicitly staged.

The fifth role is deliberately narrow. Ordinary repeated interval propagation
does not become a validated reachability tube merely because every arithmetic
operation was outward rounded.

## Report-language boundary

A mathematically valid enclosure remains a valid result even when it is too
wide to support the engineering decision. It must be reported as inconclusive
or `NoUsefulBound` for that claim, never promoted to an engineering certificate
by changing the adjective.

This is a routing result, not silence. A report must retain the requested claim,
domain/horizon, evidence route, observed bound or refusal, missing capability,
and the narrower claims that remain supportable. The doctrine's
`NoUsefulBound` row does not prove that statistical observables, local
stability, short-time events, or finite-horizon exclusions are impossible.

## Worked thermal routes

### Steady maximum temperature

Claim: maximum steady-state component temperature over declared material, load,
boundary, and discretization uncertainty.

Route: finite interval/Taylor or residual evidence for that stated steady
problem, with every model and discretization assumption retained.

Refusal boundary: the resulting bound does not establish reliability over a
population of duty cycles.

### Duty-cycle reliability

Claim: probability of thermal-limit exceedance over a declared duty-cycle
population.

Route: sequential or rare-event statistics plus sampling, dependence,
model-form, calibration, and field-transfer evidence.

Refusal boundary: a steady temperature interval or one long simulated
trajectory cannot mint this reliability claim.

## Capability maturity and future routing

`available` means a narrow serving implementation exists at the source locator;
it does not mean the complete scientific claim is established. `thin` means
useful machinery exists but the full route still needs additional evidence.
`staged` means the named capability has no implementation locator and cannot be
inferred from adjacent primitives.

The v1 table is closed and ordered. New claim families or evidence objects
require a schema-version change, explicit migration semantics, updated
no-claim boundaries, and tests that reject the old schema where the new
semantics matter. The future router must consume this table as data and remain
unable to widen a staged, thin, inconclusive, or `NoUsefulBound` route into
stronger authority.
