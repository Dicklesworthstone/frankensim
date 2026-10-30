# Three-day integration audit: flute and related instrument physics

## Scope and depth

The audit covers **1,304 individual commits** on the reviewed main history, using the approximately 72-hour window 2026-10-27 17:05:00 UTC through 2026-10-30 17:05:33 UTC and fixed head `90245fbea53699d5de59ff4ead78a016ec08e5d3`. [The per-commit inventory](2026-10-30-instrument-integration.tsv) contains one row per commit, its title, instrument relevance, integration assessment, changed-path evidence and review depth. Coverage was checked against the original commit manifest: no missing SHAs, duplicate SHAs or mismatched titles.

In the initial inventory pass, every commit received metadata and changed-path assessment. **73 selected rows additionally received relevant patch, interface, current-callsite or current-family trace review; the other 1,231 are scope assessments, not claims of full implementation validation.** The TSV preserves that pass's review depth. Its appended `followup_*` columns retain the separate review's per-commit relevance and integration assessments; they do not upgrade unrelated-domain screening into numerical qualification. The comprehensive simulation plan, instrument realism plan, open realism beads and current instrument status informed prioritization.

## Principal finding

The published flute performances still instantiate the characteristic-network `frankensim-flue-performance-v1` route with a reduced jet card. They do not select the newer resolved NSF source-to-DG path. Many new gas, wall, radiation and acoustic capabilities are connected to real native owners, but adding a capability elsewhere does not make these flute cards use it.

The repository already records a concrete source failure. At the current authored gain of 90, the D4/G4 jet displacement reaches approximately 117/204 labium half-widths RMS. The cut is effectively switching between two saturated sides. The recorded even-harmonic deficit is about −51/−58 dB for H2/H1; changing the labium offset did not repair it. A5 remains an octave low. These are **existing repository observations**, not fresh measurements made in this audit. The optional forced-Bickley route has already failed the flute regime campaign, and the Nolle/Price transfer failed its whistle control; neither is an evidence-backed default switch. The source needs admitted viscous amplification/profile evolution and physically resolved mouth geometry, followed by register and spectral validation. An output gain or timbre adjustment would not establish that physics.

Evidence: `data/instruments/descriptions/flute/status.json`, `flute-d4.performance`, `flute-g4.performance`, `flute-a5.performance`; diagnosis commit `46666de1d8`; `crates/fs-aeroac/src/forced_jet.rs`; `crates/fs-couple/src/flue_loop/forced_jet.rs`.

## Most consequential integration gaps at the audit baseline

| Area | Existing integrated capability | Remaining flute gap |
| --- | --- | --- |
| Resolved source and acoustic feedback | Material-bound `NsfGrid3` trials exchange live characteristics with DG; native parsing, renderer, independent audio clock, observer and Instrument Lab consume this path. Commits `e6c89d399f`, `0abbb9e3d9`, `c739a07afe`, `ca66b3eba0`, `a6e5cdc753`. | Published flute cards do not instantiate a resolved jet. The existing source boundary topology offers one complete characteristic axis, insufficient by itself for separate mouth atmosphere, breath/flue inflow and bore feedback. A physical lip/flue fixture and a validated flow campaign are still absent. |
| Turbulence and heat | NSF SGS combines subgrid shear, molecular bulk viscosity, stress work and heat diffusion (`eeeacdf9e9` and follow-ups). | SGS explicitly rejects characteristic and prescribed-reservoir boundaries. It therefore cannot simply be enabled on the native resolved source-to-DG path. A resolved flow run must respect resolution/transport limits; new apertures alone do not qualify playing-Re turbulence. |
| Jet models and windway | Shared `FlueFlowInput`/`FlueJetInput` already connects network and DG paths to slit, developing, unsteady and contraction flows, Nolle and forced-Bickley transfers (`27cb4106aa`). Variable-speed causal convection and physical speed semantics are also wired. | The free flute lip jet has `windway none`; developing channel laws do not address its spreading. The Rayleigh solver is a temporal eigenmode oracle; the forced jet is a parallel Bickley perturbation solver. Neither is a drop-in spatially spreading viscous source. This is a model/qualification gap, not an undiscovered parser switch. |
| Full flute bore and side openings | DG supports finite chimneys, sourced material walls, molecular gas, side channels and **closed terminal channels**. | The current source grammar is limited to seven side openings plus the axial outlet, with each closed chimney still consuming a slot. Sixteen Boehm holes plus the cork cavity cannot be carried intact by that bounded route. A wholesale card conversion would lose geometry unless this is addressed. |
| Regional and changing breath gas | Native DG live rebind preserves pressure/flow, viscothermal and wall histories. Finite gas transport and regional/molecular acoustic paths have genuine consumers. | Published flute graph playback uses one fixed moist-air snapshot. It does not select transported warm breath or continuously changing regional gas. Live gas rebinding alone does not supply the missing mass/enthalpy transport or jet constitutive update. |
| Mutual radiation and receiver sound | Full matrix BEM/network exterior owners, modal outlets and bounded passive rational fits exist. Shared physical exterior propagation is wired to flute/reed observers. | Flute cards retain compact independent terminal loads and an approximately 2.8 kHz receiver band. Complete mutual load topology and practical high-clock fitting remain blockers. The earlier clipped automatic mutual-radiation heuristic was explicitly removed by `3e14984693`; it must not be counted as current full mutual coupling. |
| Window, tone holes and walls | Spatial thermoviscous/hybrid FE interiors, molecular calorics and finite material walls now compose through harmonic multiports and duct scattering. Shared wind readers also carry directional wall response into real consumers. | The active flute does not instantiate the new spatial interior owners; finite junction volumes, undercut/pad details and a same-specimen measured geometry remain incomplete. Harmonic capabilities still need an admitted causal coupling and actual geometry qualification before performance playback can use them. |

The relevant native owners are `crates/fs-couple/src/nsf_dg.rs`, `crates/fs-lbm/src/nsf3.rs` and its boundary modules, `crates/fs-couple/src/dg_port/input.rs`, `crates/fs-couple/src/flue_loop/network/performance/file/excitation.rs`, and the native resolved-performance input/render paths. The per-commit TSV gives the changed evidence paths for each cited commit.

The tone-hole limit is deeper than the parser. Although the underlying DG owner
supports 32 side ports, native preparation builds radiation for every opening
before attaching terminal closures. Closed chimneys retain exterior matrix rows
and histories; finite branches require both original flow ports; the observer
requires matching terminal/aperture dimensions. Admitting a full flute requires
separating sealed internal branches from radiating exterior-port indices while
preserving chimney storage and bore coupling. Raising the input count or
filtering closed holes alone would break those existing ownership invariants.

## Session implementation selected from the audit

The highest-value bounded extension is a **localized prescribed-reservoir aperture on a stationary transverse NSF wall**, coexisting with the original characteristic source/bore axis. It reuses the original source material, physical ghost state, population reconstruction, thermal-pressure force and conservative total-energy flux. Geometry must avoid characteristic collars, domain corners, moving walls and solids. Native file input and NSF-to-DG construction need to expose that same owner, rather than introducing a second jet equation.

This enables a separately ventilated resolved source domain. It does not establish an acoustically transparent mouth, a finite atmospheric reservoir, a calibrated flue, a realistic flute recording, or SGS compatibility. Reservoir exchange must remain distinct from characteristic/DG exchange in the ledgers. Tests must demonstrate simultaneous exchange and conservative ownership, not only a parser accepting a new token.

**Concurrent implementation after the audit window:** commit
[`2fdf53b2b`](https://github.com/Dicklesworthstone/frankensim/commit/2fdf53b2b0ee8a77dc23d86d05c701496eaa490b),
with its boundary-owner prerequisite, now exposes selected characteristic bore
patches and transverse reservoir apertures through the original NSF/DG owners,
native inputs, two example performances and Instrument Lab. It retains separate
complete-face, selected-patch and reservoir accounts. This advances the topology
gap identified above. Its commit records passing native-client checks and
explicitly blocked executable numerical qualification. The manufactured examples
remain distinct from a resolved, validated flute mouth; escaping mouth flow
still has no separate exterior radiation observer. The native results in this
supplement predate those late boundary changes and do not qualify them.

## Prioritized continuation

1. Complete and qualify an actual resolved mouth/flue geometry with separate atmosphere and breath boundaries; couple it through the existing characteristic/DG owner. Check mass, momentum, energy, onset, flow structure and interface sensitivity before judging sound.
2. Resolve either source-resolution limits or compatible SGS/open-boundary treatment for the admitted playing regime. Validate the reduced viscous jet alternative independently if it remains the practical playback path; do not fit a gain to mask the relay source.
3. Preserve the complete Boehm bore, all closed/open chimneys and cork when moving to DG, then supply one physically consistent exterior and receiver model at a useful band.
4. Integrate transported warm breath and the relevant material state after the source topology and full acoustic graph are supported. Complete same-specimen passive impedance, register, harmonic and listening tests through the existing realism beads.

These priorities map to the instrument realism plan’s M3.2–M3.4 resolved source/coupling ladder and flute section 6.8, and the open `.133`–`.141`, `.19`, `.207`, `.319`–`.323` realism beads. Generic high-temperature chemistry, phase conversion, optical radiation and contact mechanics were assessed for scope but are not appropriate additions merely because they are recent. Bowed/plucked strings, reeds, percussion and transducers have real separate consumers documented in the per-commit inventory; this flute-focused audit does not claim those families are comprehensively qualified.

## Supplement: source-to-consumer repairs from the parallel review

A separate review used the exact UTC interval **2026-10-27T17:04:57Z through 2026-10-30T17:04:57Z**, at the same fixed baseline `90245fbea53699d5de59ff4ead78a016ec08e5d3`. It resolves to the **same 1,304-commit set**. Its assessments supplement the existing inventory rather than creating additional reviewed commits. Current-callsite tracing identified three bounded defects in existing wind consumers, alongside the larger source-topology work selected above.

### Canonical flute observation now consumes the exterior gas law

The baseline has shared physical gas propagation in DG wind and other pressure observers, but the characteristic-network `FlueFilePerformance` route lacks the complete attachment. The patch connects `GasPropagation` through `DuctPressure`, `FlueOutletPerformance`, native source parsing, and `music_render` receipts.

Each terminal opening and the mouth retains its own distance, phase, and absorption. Common numerical padding aligns the separately prepared mouth and terminal pressure histories before summation. The explicit gas law therefore reaches the actual listener-pressure calculation rather than remaining a parsed record or applying one filter after unrelated paths have been mixed.

The optional record leaves source pressure, jet dynamics, and bore state under their existing owners. Descriptions that omit it retain their original observation behavior. Relevant changes and regressions are in `baffled_pressure_duct.rs`, `flue_loop/network/performance.rs`, its `file.rs` and tests, and `tests/duct_pressure_packets.rs` under `crates/fs-couple`.

Concurrent commit
[`da45542bd`](https://github.com/Dicklesworthstone/frankensim/commit/da45542bd3b3bba88acd0d21ca6036db435765d1)
then selected `observer_gas stokes-kirchhoff 4096 1000000 50000000` in the
published D4, G4 and A5 descriptions. That activation and its extra direct-owner
tests are preserved. The integration is therefore selected by those presets;
the estimated source, authored gain and outstanding register/realism problems
remain separate.

### Resolved observation now parses molecular propagation against its gas source

The resolved source reader already has an observation attachment, but called `GasPropagation::read_record` without the exterior gas context. `nsf_dg/performance/file.rs` now calls `read_record_for_gas`. Molecular propagation can consequently use the original exterior gas record and its source-defined relaxation law, with the same admission rules as other native consumers.

Concurrent commit `1cede3a1b` delivered the same parser repair during this review; that implementation and its molecular-pressure tests are retained. The follow-up adds a distinct-upstream-exterior regression in `nsf_dg/performance/file/tests/pressure.rs`, so selecting the last parsed gas record cannot accidentally bind the listener to the upstream atmosphere. This closes a source-record-to-observer gap. It does not change the resolved flow’s composition transport, add a lip boundary, or qualify an atmospheric interface.

### Jet-card centre speed now participates in pressure admission

Commit `84bcc88c0` connects exit centre speed to card convection, growth, Reynolds number, and onset seeding, while bulk speed continues to carry volume flux. The associated pressure controls and schedule admission still checked only bulk speed. The fully developed slit’s Poiseuille centre was also checked against its inviscid bound only at initial attachment’s supplied reference pressure.

`card_speed_at` now receives the proposed pressure explicitly, and shared `validated_jet_speed` checks the corresponding centre speed before accepting a control. The declared Mach limit and `sqrt(2P/rho)` bound apply consistently. Cold attachment validates both its reference pressure and any already accepted input; bulk flow retains its original physical law.

New `flue_windway` regressions exercise a legal bulk speed with an illegal centre speed, pressures on both sides of a limit derived from the slit law, and a developing core proposed from a zero-pressure state. A refusal test compares subsequent output with an untouched equivalent voice, covering retained pressure and source history as well as immediate state. This matters directly to shared flue/recorder source semantics; the shipped free flute jet has `windway none`, so this correction alone should not be credited with improving its timbre.

### Recent interface changes also left callers behind

The native build exposed additional integration breakages. The staged double-acting acoustic `StagedPort` had not followed `FiniteGasTransportNode`’s new molecular-heat preparation, advance argument, and fourth accounting term. It now delegates preparation to the complete-load owner, forwards the accounting shape, and explicitly refuses a nonzero molecular heat integral that this passive-load wrapper does not own. A focused regression checks refusal before caloric trials; silently discarding the new heat term would conceal a physical accounting error.

Other repairs restore the actual ownership/visibility of the articulated-contact `Term` import, circuit unknown admission, and thermoelastic conductor-geometry refresh/restore helpers. New wet-frame tests now call the existing gate-aware `step_under_gate` interface, and reed contact tests import their solver types. These changes repair concrete trait/caller seams encountered while building the wind changes; they add no new flute excitation model.

### Validation and remaining physical boundary

Validation used the repository-pinned `nightly-2026-07-06` toolchain and locked,
offline Cargo dependencies. DSR and RCH were unavailable in this environment.

- `cargo test --locked --offline -p fs-couple --test flue_windway --test duct_pressure_packets`
  passed **18 tests**: all five duct pressure tests and all thirteen windway
  tests, including the new source-speed and material-path regressions. This
  compiled the production coupling crate and native renderer through
  `ddf9f2100`. Later concurrent topology/ensemble changes are not covered by
  that completed build.
- Eight direct `music_render wind` cases passed their expected outcomes:
  legacy observation, classical propagation at callback sizes 37 and 257,
  explicit molecular exterior gas, three invalid-input refusals, and a resolved
  source with distinct upstream/downstream atmospheres. The two callback sizes
  produced identical PCM; classical propagation produced PCM different from
  the legacy case. Receipts retained 32 terminal surface paths, eight mouth
  paths, one exterior gas, and common additional numerical latency.
- Missing molecular species, a zero preparation budget, and unexpected trailing
  fields refused before WAV or sidecar creation. The resolved case used the
  explicitly requested 96-to-48 kHz decimation and accepted the downstream
  N2/CO2 relaxation law despite the separate upstream N2/O2 exterior.
- The full `fs-couple --lib` test-binary build was stopped while contending at
  the shared 8 GiB memory ceiling. Its in-module tests are **not reported as
  passing**. The narrower native suites and renderer cases above completed.
  The shortened observer fixture was corrected to use an in-window held
  pressure after the native renderer correctly refused its old, longer gesture.
- Additional short-onset checks using the full published D4/G4/A5 geometries
  were started after their observer activation. The shared filesystem reached
  zero available bytes during the first case; that extra render was stopped.
  No full-preset pass, full-note render or audible-realism result is claimed.
- Changed-range formatting and `git diff --check` passed. The constellation
  checkout preflight was not fully clean (NumPy legacy-file raw bytes and
  NetworkX/Pandas checkout-index materialization); this is not a claim that the
  repository-wide CI/DSR gates passed.

Published implementation checkpoints are
[`4aa8d7522`](https://github.com/Dicklesworthstone/frankensim/commit/4aa8d7522291688edc394e6a7d8e091788e1d4c5)
(source-speed admission),
[`68cf23ca9`](https://github.com/Dicklesworthstone/frankensim/commit/68cf23ca94bc32741b00ac529ec40476f74b1e5f)
(complete canonical-flue gas observation), and
[`ddf9f2100`](https://github.com/Dicklesworthstone/frankensim/commit/ddf9f2100539cc60d7597ba5639b8c07a7c2a08c)
(cross-module heat/force API repairs). The concurrent resolved-observer repair
`1cede3a1b` is retained and covered by the distinct-atmosphere check. This work
contains no fresh graded listening result or recording comparison.

The late concurrent boundary commits implement the separate aperture topology described earlier; these three consumer repairs and their completed native checks do not establish its numerical qualification. The dominant source limitations remain those documented above and in realism beads `.135`, `.138`, `.188`, `.207`, and `.320`–`.323`: playing-Re excitation, complete embouchure/bore geometry, correct closed-loop regimes, and measured listener output. Plan M3.2–M3.4 and §6.8 require that physical composition. Existing forced-Bickley register failures and the Nolle/Price whistle failures recorded in `.320` comments 7303/7325 still preclude treating a default-model switch as the completed repair.

## Reservoir owner checkpoint and complementary composition coverage

Commit [`9296d4fe2`](https://github.com/Dicklesworthstone/frankensim/commit/9296d4fe28e814b0175103e2cee421a632972153) added the localized transverse `NsfReservoirAperture` owner selected by this audit. It retains the original gas/material law, population streaming, pressure force and finite-volume energy flux. Sparse sorted tile ownership keeps aperture-edge streaming and independently counted mass consistent; cold admission, reversible trials and cancellation retain the same physical owner. Mixed-boundary internal/kinetic exchange is reported as unresolved instead of being relabeled viscous heating.

The subsequent selected-patch/native implementation cited above also connected this reservoir to the native reader and both NSF/DG frame constructors. The final reconciliation preserves that implementation and its source grammar: the reservoir declaration precedes `initial`, the count follows the SI state, and tiles use `reservoir_cell`. It preserves both native examples, selected-patch feedback, separate reservoir mass, the existing dry-air/conservation tests and the parallel review supplement.

Complementary regressions in `crates/fs-couple/src/nsf_dg/performance/file/tests/aperture.rs` add three distinct composition checks:

- Exact replay through the original fluid trial and original DG load, including reservoir mass and energy, with a warmer reservoir and cancellation/retry.
- Retained thermochemical species, NASA9 caloric and transport sources after the borrowed inputs are dropped; changing only calorics at fixed geometry, mechanical clock and reservoir state must change solved pressure.
- A finite upstream duct together with the side reservoir; upstream termination must affect the main bore, all three physical openings retain their accounts, and cancellation or an unbracketed pressure rise must preserve the fluid and both ducts before retry.

Native reservoir provenance now explicitly names the mass/energy account fields, SI units and opposite positive directions. The native CLI regression also consumes the existing side-aperture example through the same original pressure/decimation comparison and checks those receipt fields. The current native grammar, characteristic patch and independently maintained observer integrations are preserved.

### Verification boundary for this checkpoint

The pinned compiler and genuine production dependencies passed **23 focused native runtime tests** on the reservoir-owner checkpoint: six new aperture tests, fourteen existing open-boundary tests and three existing characteristic-trial tests. These were rebuilt and passed after rebuilding the current material dependency. The aperture tests compile the exact original NSF module tree; they do not substitute numerical models or stubs.

The patched full `fs-couple` production library and actual `music_render` binary source also passed metadata compilation on the pre-reconciliation verification tree. Those checks establish type integration for that tree; they do not cover the later selected-characteristic-patch implementation now preserved on main. The broad local unit-test check did not complete and supplies no passing result. Its old local import/API diagnostics were already corrected on concurrent main and are not reported as current-main defects.

The three complementary coupled cases and the extended native CLI case are **added but not runtime-verified**. The linked build was constrained by exhausted shared storage. The native-client suite passed 21 tests, including exact forwarding of resolved sources. These checks are not a clean Cargo/lockfile workspace rebuild, full-preset render, physical-mouth convergence study or new listening result. The published flute's source/regime problems remain open even though its physical gas observation and the resolved-source boundary topology have advanced.
