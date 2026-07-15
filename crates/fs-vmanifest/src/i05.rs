//! The I05 HIL/WCET/fixed-point deployment-twin compiler
//! VerificationManifest draft (bead
//! frankensim-leapfrog-2026-program-i94v.3.1.7.1).
//!
//! The baseline lattice freezes bounded DeploymentTwinIR admission,
//! range- and error-accounted fixed-point synthesis, reproducible safe
//! `no_std` binaries, non-confusable timing-evidence classes, explicit
//! multi-clock HIL semantics, fail-closed fault/safe-state behavior, and
//! bounded target-refinement evidence. The stronger lattice attempts
//! compositional deadline closure, static/measurement reconciliation,
//! machine-checked floating-to-target refinement, and a complete WCET
//! argument for a deliberately finite target model. A separate refutation
//! lane actively hunts false equivalence certificates. None of these
//! preregistrations is itself target, timing, safety, or theorem authority.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

const CAMPAIGN_POLICY_FIXTURE: &str = "i05-campaign-policy-v1";

/// Build the I05 draft. Consumers freeze it themselves; focused
/// conformance tests prove that the authored seed satisfies schema v2.
#[must_use]
pub fn i05_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "I05",
        title: "Deployment-twin compiler gate: bounded fixed-point code, target-specific timing, \
                explicit HIL clocks and faults, and theorem-qualified floating-to-target refinement",
        version: 1,
        explicits: FiveExplicits {
            units: "SI base units at model and HIL ports; target integer encodings carry an exact \
                    rational scale and offset; time is seconds plus integer target cycles; frequency \
                    is hertz; memory is bytes; probability and normalized error are unit '1'; \
                    exact/refusal/coverage verdicts use unit 'bit'",
            seeds: "Philox 4x32-10 streams keyed by BLAKE3 domain \
                    org.frankensim.i05.fixture-stream.v1 over the exact fixture id, target-profile \
                    digest, and case index; development indices 0..=16383, core held-out indices \
                    65536..=81919, and maximal held-out indices 131072..=147455 are disjoint; \
                    hardware noise/fault schedules retain raw authenticated draws",
            budgets: "smoke <= 90 s and 2 GiB; core <= 45 min and 16 GiB; maximal <= 12 h and \
                      32 GiB per pinned target; code-size <= 2 MiB per image; bounded-horizon \
                      refinement <= 10^7 target ticks; timing exploration <= 10^9 analyzed states \
                      per finite target profile; every cap exhaustion is BudgetExhausted/Unknown, \
                      never evidence of a deadline or equivalence",
            versions: "fs-vmanifest schema v2; DeploymentTwinIR v1; TargetProfile v1; \
                       FixedPointSemantics v1; TimingEvidenceReceipt v1; HilTranscript v1; \
                       rust-toolchain.toml and constellation.lock digests pinned by campaign receipt",
            capabilities: "safe Rust only in generated production controller code; no allocator, network, \
                           unwinding, FFI, inline assembly, hidden floating point, or wall-clock reads \
                           in that admitted no_std subset; any boot/vector/HAL/interrupt/probe runtime capsule \
                           is non-generated, separately content-addressed, capability-allowlisted and audited \
                           behind safe facades; deterministic mode mandatory; target, \
                           linker, compiler, binary, probe, clock, and HIL-rig identities mandatory; \
                           maximal theorem and exhaustive-WCET lanes remain feature-gated",
        },
        claims: i05_claims(),
        fixtures: i05_fixtures(),
        obligations: i05_obligations(),
        waivers: i05_waivers(),
        amendment_rules: "Any semantic change creates the exact next manifest version through \
                          FrozenManifest::amend. Candidate/model/toolchain/target/profile, claim \
                          text, arithmetic, timing model, clocks, partitions, bands, seeds, budgets, \
                          checker identities, and policy freeze before evidence. The amendment \
                          receipt invalidates affected predecessor claim and leaf authority; \
                          observed results never edit this version in place or tune a successor \
                          without a new campaign lineage.",
    }
}

#[allow(clippy::too_many_lines)]
fn i05_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "i05-deployment-ir-admission",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "DeploymentTwinIR canonical encoding round-trips byte-identically and \
                        admission rejects every enumerated malformed, dimensionally inconsistent, \
                        unbounded, nondeterministic, or target-incompatible construct with a stable \
                        typed diagnostic before code generation",
            hypotheses: &[
                "the IR belongs to the frozen bounded grammar and every state, port, clock, task, \
                 fault, safe state, numeric format, loop bound, and target resource has a stable id",
                "admission receives the exact TargetProfile, compiler-policy, and capability digests",
            ],
            qoi: "canonical_roundtrip_and_admission_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "independent schema-v1 decoder plus a separately authored malformed-class \
                           recognizer over canonical bytes",
                independent: true,
                tcb_overlap: "shares only primitive UTF-8 and integer decoding, not the production \
                              admission or encoder implementation",
            },
            activation: "the I05 DeploymentTwinIR implementation leaf opens",
            kill: "one accepted prohibited construct or one noncanonical round trip blocks every \
                   downstream I05 lane",
            fallback: "refuse code generation and emit the bounded explicit-interpreter target plan",
            no_claim: "schema admission does not prove controller stability, timing, physical fidelity, \
                       safe behavior, or target equivalence",
        },
        ClaimSpec {
            id: "i05-fixed-point-range-soundness",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "For admitted bounded expressions, synthesized signedness, word lengths, \
                        binary/rational scales, rounding, saturation/trap policy, and accumulator \
                        widths enclose every reachable mathematical value and the reported local \
                        quantization-error bound encloses the exact rational fixed-point error",
            hypotheses: &[
                "all input/state/parameter ranges and loop horizons are frozen finite sets or certified \
                 intervals; NaN, infinity, subnormal, and implicit cast behavior is absent",
                "every operation has an explicit infinite-precision semantics and one explicit target \
                 rounding and overflow semantics; wraparound is admitted only when the source law \
                 explicitly specifies modular arithmetic",
                "range propagation either returns a checker-verifiable enclosure or Unknown; sampling \
                 alone can never close this claim",
            ],
            qoi: "range_and_local_error_containment_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G1,
            oracle: OracleRoute {
                identity: "exact-rational exhaustive evaluator for enumerable microprograms plus \
                           independently implemented interval/range certificate checker",
                independent: true,
                tcb_overlap: "the checker shares the frozen operation semantics but no synthesis or \
                              range-propagation code",
            },
            activation: "DeploymentTwinIR admission is green on the development corpus",
            kill: "one reachable out-of-range value, unreported saturation/trap, or underestimated \
                   exact-rational error refutes the synthesized numeric format",
            fallback: "increase precision or use a slower checked integer kernel; if no admitted format \
                       closes, retain the floating model and refuse deployment",
            no_claim: "local arithmetic containment is not closed-loop refinement, stability, WCET, or \
                       freedom from analog ADC/DAC and sensor error",
        },
        ClaimSpec {
            id: "i05-reproducible-no-std-binary",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "The admitted compiler emits safe no_std Rust whose source, object graph, linked \
                        image, memory map, ABI description, and symbol/section inventory reproduce \
                        byte-identically from the frozen IR and toolchain on the same target/toolchain \
                        fingerprint, with no undeclared runtime capability",
            hypotheses: &[
                "the exact rustc/LLVM target specification, linker, linker script, flags, environment, \
                 source-date epoch, and build inputs are content-addressed",
                "the generated controller subset excludes allocator, panic unwinding, recursion, dynamic \
                 dispatch, unbounded loops, unsafe code, FFI, assembly, and hidden target-feature dispatch; \
                 separately pinned boot/vector/HAL/interrupt capsules are explicit TCB inputs",
            ],
            qoi: "binary_and_capability_replay_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G5,
            oracle: OracleRoute {
                identity: "binary-format parser and capability allowlist checker independent of the code \
                           generator, plus a clean-room rebuild comparator",
                independent: true,
                tcb_overlap: "shares the pinned compiler/linker binaries but not generator logic or the \
                              first build directory",
            },
            activation: "fixed-point synthesis returns a checker-accepted format assignment",
            kill: "a byte mismatch under the same fingerprint, undeclared import/section/instruction, or \
                   memory-map overflow blocks target loading",
            fallback: "emit source and proof obligations only; no flashable image is promoted",
            no_claim: "same-input binary reproducibility does not establish semantic equivalence, WCET, \
                       security hardening, correctness of the separately audited target-runtime capsule, \
                       or certification of the toolchain",
        },
        ClaimSpec {
            id: "i05-timing-evidence-separation",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Every timing result is typed as StaticUpperBound, CompositionalUpperBound, \
                        MeasuredSampleMaximum, or Unavailable; it binds binary, target, clock, cache/bus/ \
                        interrupt assumptions, task context, path coverage, and units, and no measured \
                        maximum is relabeled as WCET",
            hypotheses: &[
                "cycle counters and probes have calibrated overhead and monotonicity receipts",
                "static bounds name loop/call bounds and a finite microarchitectural model; compositional \
                 bounds additionally name interference and scheduling contracts",
            ],
            qoi: "timing_evidence_kind_and_binding_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independent receipt decoder/binary-hash checker plus adversarial timing-kind \
                           confusion corpus",
                independent: true,
                tcb_overlap: "none beyond the public TimingEvidenceReceipt schema",
            },
            activation: "a reproducible target image exists",
            kill: "a stale binary/target binding, unit confusion, hidden assumption, or measurement-as- \
                   WCET promotion is a fail-closed integrity error",
            fallback: "report the measured distribution and Unavailable WCET separately; do not admit a \
                       hard real-time deployment",
            no_claim: "finite testing cannot prove an unobserved worst case; even a valid static bound \
                       may be too loose for deadline admission",
        },
        ClaimSpec {
            id: "i05-hil-clock-causality",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "HIL transcripts preserve distinct simulation, host monotonic, DAQ, bus, and \
                        device-clock domains with calibrated affine/piecewise clock maps and uncertainty; \
                        every sample, actuation, deadline, dropout, duplicate, reordering, and reset has \
                        a causality interval, and lag/catch-up cannot masquerade as on-time execution",
            hypotheses: &[
                "each clock has a stable identity, tick rate, wrap rule, calibration interval, and \
                 synchronization uncertainty bound",
                "device and host I/O paths expose bounded timestamp placement or report Unknown; \
                 interpolation never invents physical samples",
            ],
            qoi: "clock_mapping_and_causality_containment_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent logic-analyzer/PTP trace aligner with injected wrap, drift, \
                           jitter, dropout, duplicate, and reorder schedules",
                independent: true,
                tcb_overlap: "shares the physical timestamp sources but not the HIL transcript builder \
                              or clock-map fitter",
            },
            activation: "the target image and rig profile pass identity checks",
            kill: "one event outside its claimed causality interval or one hidden missed deadline blocks \
                   HIL evidence",
            fallback: "widen intervals or mark the run TimingUnknown; open-loop offline replay remains \
                       available but carries no real-time claim",
            no_claim: "clock containment does not prove plant fidelity, deterministic Ethernet, or \
                       deadline satisfaction outside the pinned rig and load",
        },
        ClaimSpec {
            id: "i05-safe-state-fault-semantics",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Every injected sensor, actuator, compute, bus, clock, power, and task fault has \
                        a typed activation/detection/recovery interval and expected state transition; \
                        the runtime reaches the fixture-declared safe output or an explicit containment \
                        failure without publishing a normal-success receipt",
            hypotheses: &[
                "safe states, interlocks, watchdog windows, fault observability, simultaneous-fault \
                 policy, and de-energization semantics are fixture inputs approved before execution",
                "injection fidelity is separately classified as software, electrical, protocol, or \
                 physical; one class is never silently generalized to another",
            ],
            qoi: "fault_transition_and_safe_output_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G4,
            oracle: OracleRoute {
                identity: "independent rig interlock/logger and physical output probe checked against \
                           the fixture fault-state automaton",
                independent: true,
                tcb_overlap: "shares the rig wiring and target power supply but not target control code \
                              or in-band telemetry",
            },
            activation: "clock-causality calibration is green and the external emergency stop is armed",
            kill: "unsafe output, missing containment transition, silent fault, or success publication \
                   after containment failure halts the rig and refutes this lane",
            fallback: "de-energize through the independent interlock and retain a FailureBundle; no \
                       automatic retry of the same physical fault",
            no_claim: "this is a scoped engineering test, not IEC/ISO regulatory certification, a \
                       complete hazard analysis, or proof against unmodeled common-cause faults",
        },
        ClaimSpec {
            id: "i05-bounded-target-refinement",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "On the frozen bounded-horizon controller/plant interface, target traces refine \
                        the floating reference under the declared input/noise/fault relation: discrete \
                        mode labels agree and every named state/output/safety-margin QoI stays inside its \
                        directed interval after clock and quantization uncertainty are accounted",
            hypotheses: &[
                "the reference uses the same sampled-data semantics, saturation, delay, reset, and \
                 finite input relation as the target rather than an idealized continuous controller",
                "the plant/interface model and HIL I/O error envelopes are frozen; the claim horizon and \
                 initial set are finite",
                "comparison aligns traces through the clock-causality intervals and never picks the \
                 most favorable timestamp inside an uncertainty window",
            ],
            qoi: "worst_normalized_trace_refinement_error",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "separately implemented exact-rational target emulator plus floating-reference \
                           replay and interval trace aligner",
                independent: true,
                tcb_overlap: "shares fixture inputs and frozen arithmetic semantics, not emitted code, \
                              production emulator, or target telemetry decoder",
            },
            activation: "numeric, binary, clock, and fault baseline leaves are green at core tier",
            kill: "one held-out trace exceeds normalized error 1, changes an undeclared mode, or violates \
                   a safety margin; report the minimized trace and stop promotion",
            fallback: "increase precision/rate, narrow ContextOfUse, or keep the floating controller in \
                       simulation-only status",
            no_claim: "bounded empirical refinement is not an unbounded theorem, plant-model validation, \
                       stability proof, WCET proof, or authority outside the frozen target and horizon",
        },
        ClaimSpec {
            id: "i05-compositional-deadline-contracts",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Assume-guarantee timing contracts compose task-local execution bounds with \
                        scheduler, preemption, interrupt, DMA, cache, bus, and shared-resource interference \
                        into end-to-end response-time and age-of-data upper bounds for the frozen task graph",
            hypotheses: &[
                "every local bound is an actual upper bound for the same binary/target/context and every \
                 shared resource has one non-overlapping ownership/interference model",
                "release jitter, priority, blocking, preemption cost, interrupt masks, DMA masters, and \
                 clock conversion bounds are complete for the admitted profile",
                "every normalized response/age denominator is a finite strictly positive frozen deadline \
                 or age limit; a missing or zero limit is Unsupported rather than favorable arithmetic",
            ],
            qoi: "maximum_normalized_response_or_age_ratio",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "independent event-driven schedule enumerator on bounded adversarial task graphs \
                           plus trace-based lower-bound falsifier",
                independent: true,
                tcb_overlap: "shares frozen task/resource semantics only; response-time equations and \
                              enumeration implementation are disjoint",
            },
            activation: "timing evidence is correctly typed and all required upper-bound inputs exist",
            kill: "an enumerated or measured response exceeds the composed upper bound, or an interference \
                   owner is absent/double-counted",
            fallback: "serialize tasks, reserve resources, lower rates, or report DeadlineUnknown; measured \
                       slack alone never repairs an invalid upper bound",
            no_claim: "no portability across binaries, silicon revisions, clock/power states, schedulers, \
                       or interference profiles not named in the receipt",
        },
        ClaimSpec {
            id: "i05-static-measurement-reconciliation",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Static/compositional timing bounds and authenticated measurements are reconciled \
                        without collapsing their evidence kinds: every observation lies below its applicable \
                        upper bound, coverage gaps remain explicit, and statistically surprising slack or \
                        exceedance triggers model/profile investigation rather than post-hoc rebasing",
            hypotheses: &[
                "binary, target state, workload, probe, clock, and task-context identities match exactly",
                "measurement partitions and sequential diagnostics are frozen before the held-out campaign; \
                 optional stopping uses anytime-valid e-process thresholds",
            ],
            qoi: "bound_exceedance_count_and_coverage_gap_verdict",
            unit: "1",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent receipt joiner and e-process implementation over retained raw timing \
                           traces, separate from static analyzer and target harness",
                independent: true,
                tcb_overlap: "shares authenticated trace bytes and public receipt schema only",
            },
            activation: "compositional timing contracts and held-out rig campaign are admitted",
            kill: "one unexplained upper-bound exceedance refutes the applicable bound; identity mismatch \
                   or coverage gap yields IntegrityFailed/Unknown, not a favorable comparison",
            fallback: "retain both evidence classes and localize the mismatch; hard real-time promotion stays \
                       blocked until a valid upper bound is re-earned",
            no_claim: "absence of observed exceedance does not tighten an upper bound or estimate an extreme \
                       quantile beyond the preregistered statistical claim",
        },
        ClaimSpec {
            id: "i05-machine-checked-target-refinement",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "A machine-checked theorem proves bounded-horizon trace refinement from the frozen \
                        sampled-data real-arithmetic controller to the emitted fixed-point binary semantics, \
                        including rounding, saturation/traps, scheduling delay, clock uncertainty, mode \
                        transitions, and the declared fault relation",
            hypotheses: &[
                "a successor manifest freezes a canonical proposition AST, all definitions, target decoder \
                 semantics, proof-assistant translation, axiom policy, and total runtime-premise map before proof",
                "all generated range/error/timing premises are checker-accepted and bound to the exact binary, \
                 target profile, initial/input set, plant abstraction, and finite horizon",
                "the proof's observable relation is no weaker than the public QoI/safety contract and carries \
                 explicit error accumulation rather than equating real and finite arithmetic",
            ],
            qoi: "kernel_checked_refinement_theorem_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G5,
            oracle: OracleRoute {
                identity: "pinned Lean kernel replay plus independent canonical-AST-to-Lean round-trip and \
                           runtime-premise checker",
                independent: true,
                tcb_overlap: "the theorem consumes exported binary semantics and certificates but shares no \
                              production compiler proof-search code; kernel/axioms are declared TCB",
            },
            activation: "feature i05-target-refinement-theorem; a pre-proof successor has frozen machine \
                         proposition/definition bytes and every baseline/frontier premise is green",
            kill: "translation mismatch, unresolved premise, disallowed axiom/postulate, kernel rejection, \
                   or one admitted counterexample leaves the theorem Refuted/Unknown as applicable",
            fallback: "retain bounded empirical target-refinement evidence with its existing no-claim boundary",
            no_claim: "version-1 prose mints no theorem authority; no unbounded time, liveness, physical-plant, \
                       compiler-correctness, or target behavior outside the exact decoded binary semantics",
        },
        ClaimSpec {
            id: "i05-complete-finite-target-wcet",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "For a deliberately finite target semantics, a checker-verified exhaustive/symbolic \
                        exploration proves the exact maximum cycle count of every admitted task entry over all \
                        declared inputs, initial microstates, cache/bus/interrupt states, and bounded loop paths",
            hypotheses: &[
                "a pre-campaign successor freezes complete executable instruction, pipeline, memory, cache, \
                 bus, interrupt, DMA, clock-state, reset, and arbitration semantics for one exact silicon profile",
                "the finite state grammar, valid-state predicate, initial set, transition relation, cost model, \
                 symmetry quotient, rank/unrank or symbolic coverage method, exclusion order, and completeness \
                 certificate are machine-readable and independently decoded",
                "the analyzed image has no self-modification, unbounded input, asynchronous source outside the \
                 profile, or undocumented microarchitectural state",
            ],
            qoi: "exact_worst_case_cycle_count_and_argmax",
            unit: "cycle",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G5,
            oracle: OracleRoute {
                identity: "independent target-semantics decoder, coverage/Merkle checker, maximizing-trace replay, \
                           and hardware lower-bound falsifier",
                independent: true,
                tcb_overlap: "shares frozen target semantics and binary bytes only; search, quotient, and coverage \
                              checkers are separately implemented",
            },
            activation: "feature i05-finite-target-wcet; machine-readable finite semantics and a whole-campaign \
                         cost preflight are frozen before candidate exploration",
            kill: "unmodeled state, incomplete coverage, quotient unsoundness, cost mismatch, hardware trace above \
                   the result, or budget exhaustion prevents an exact-WCET claim",
            fallback: "use a conservative static/compositional upper bound or report WCET Unavailable; preserve \
                       explored traces as non-exhaustive lower-bound evidence",
            no_claim: "version-1 target prose is not an executable state grammar and mints no exhaustive authority; \
                       a model-exact maximum is silicon WCET only after independent model-fidelity qualification",
        },
        ClaimSpec {
            id: "i05-false-equivalence-certificate-falsifier",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Refutation,
            statement: "Adversarial numeric, scheduling, clock, mode-boundary, compiler, and target-model mutations \
                        try to obtain a green end-to-end refinement or WCET certificate for a known-nonrefining or \
                        underbounded program; any accepted mutant refutes the corresponding certifier",
            hypotheses: &[
                "mutants preserve syntax and enough superficial checks to exercise semantic binders: signedness, \
                 scale, rounding, saturation, stale binary hash, omitted interrupt, clock wrap, timestamp choice, \
                 dropped mode, weakened relation, unsound symmetry quotient, and target-semantic drift",
                "each mutant has an independently generated witness trace or proof of omitted state and is hidden \
                 from certifier development until adjudication",
            ],
            qoi: "false_certificate_count",
            unit: "1",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "independent mutation generator plus exact-rational trace witnesses and finite-state \
                           countermodel checker",
                independent: true,
                tcb_overlap: "none with production certificate generation; public schemas only",
            },
            activation: "a maximal certifier candidate is frozen and its checker interfaces are sealed",
            kill: "the intended lane succeeds when any false certificate is accepted: retain the minimized mutant \
                   and mark the certifier Refuted; zero false certificates is merely falsification evidence",
            fallback: "disable the affected certifier and fall back to the strongest independently surviving lower \
                       rung; never relabel a failed maximal certificate as baseline proof",
            no_claim: "a finite mutation campaign cannot prove certifier soundness or completeness and never upgrades \
                       Unknown to Proved",
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i05_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: CAMPAIGN_POLICY_FIXTURE,
            source: FixtureSource::AuthoredSpec {
                spec: "I05_CAMPAIGN_POLICY_V1\n\
TARGET_IDENTITY=TargetProfile digest + silicon/board revision + ISA/ABI/endianness + clock/power state + memory/cache/TCM/bus/interrupt/DMA model + debug/probe state + compiler/linker/linker-script/flags + emitted source/object/link-map/binary digests; a mismatch is IntegrityFailed, never tolerated drift\n\
ARITHMETIC=each value binds signedness,width,exact rational scale,offset,rounding,overflow and unit; operations have infinite-precision reference semantics and explicit target semantics; saturation/trap/modular wrap are distinct; range/error certificates return Accepted,Refuted,Unknown,IntegrityFailed; testing alone cannot accept a universal range\n\
TIMING_KINDS={StaticUpperBound,CompositionalUpperBound,MeasuredSampleMaximum,Unavailable}; only the first two may carry upper-bound authority under their exact assumptions; a measurement is always a lower bound on the true maximum observed search space and never WCET; cycles and seconds remain distinct with a clock-bound conversion receipt\n\
CLOCKS={simulation,host_monotonic,DAQ,bus,device,target_cycle,external_probe}; clock maps bind calibration interval,wrap,reset,rate/offset/drift enclosure and timestamp-placement uncertainty; missed,late,dropped,duplicated,reordered,interpolated and caught-up events remain distinct\n\
FAULTS=fault record binds source,class,activation interval,detection interval,recovery/containment interval,simultaneity,observability and physical-injection fidelity; normal completion, safe-state completion, containment failure, infrastructure failure and emergency-stop termination are non-confusable\n\
REFINEMENT=reference and target share sampled-data,delay,saturation,reset and fault semantics; align through whole clock intervals without favorable point selection; compare mode trace plus directed per-QoI intervals over a frozen finite horizon; empirical, kernel-checked and physical-model-qualified authority are separate\n\
THEOREM_AUTHORITY=version 1 has prose cards only and mints no proof/exhaustive-search authority; a pre-proof successor must freeze canonical proposition/definition/target-semantics ASTs, total runtime-premise schema, deterministic AST-to-Lean translation, exact axiom allowlist {propext,Quot.sound,Classical.choice} and transitive closure; sorryAx, custom postulates, native/unsafe oracle shortcuts and unbound generated premises are IntegrityFailed; a pre-search successor must freeze the executable finite-state grammar, validity and initial predicates, transition/cost semantics, enumeration or symbolic coverage algorithm, quotient proof obligations, independent decoder, whole-campaign preflight and completeness root\n\
EVIDENCE_STATES=Execution{Succeeded,Failed,Cancelled,TimedOut,BudgetExhausted,InfrastructureFailed,IntegrityFailed}; Predicate{Accepted,Refuted,Unknown,NotEvaluated}; Claim{NoClaim,EvidencePartial,Reproduced,Refuted,Proved}; Promotion{Blocked,CoreEligible,MaxEligible}; one axis never substitutes for another and no green badge collapses them\n\
HOLDOUT=the exact candidate, model, target/toolchain, thresholds, seeds, policy and checker versions freeze before held-out access; core and maximal partitions have one named consumer stage; premature/cross-stage read, replacement, retry-after-result, label leak or post-result amendment is IntegrityFailed\n\
LIFECYCLE=request->drain->finalize; cancellation polls at IR node, basic block, task release, HIL sample and proof-search tile boundaries; losing/speculative work drains; checkpoints bind manifest/candidate/target and completed logical ids; resume/fork never republishes accepted ids\n\
LOGGING=bounded schema-versioned fs-obs JSONL with stable run/case/claim/fixture/leaf/target/binary/task/clock/fault/attempt ids, units, seeds, budgets, versions, capabilities and redaction reason; large/raw streams content-addressed; stdout is not evidence\n\
RETENTION=retain manifest, exact command, code/contract/schema/toolchain/target/binary hashes, canonical inputs, raw authorized HIL/probe traces, clock calibrations, timing assumptions, range/proof/coverage certificates, all terminal states, minimized counterexamples and replay-verifier result; promotion evidence/refutations are durable; inaccessible/redacted/expired evidence constrains replay and never becomes silently verified\n\
FAILURE_BUNDLE=first semantic divergence plus bounded context, exact inputs/state/task/clock/fault, binary and target ids, expected/actual intervals, causal predecessors, minimized reproducer, checker disagreement, terminal state and replay command; partial success cannot publish normal authority\n\
PROMOTION=baseline claims consumed only by I05.G4 after independent G2 reproduction and G3 falsification; maximal claims consumed only by I05.G7 after independent G5 reproduction and G6 red-team; missing/stale/waived/integrity-failed evidence cannot promote\n\
LEAF_REQUIREMENT=every obligation references this policy fixture, declares all nine unit-case classes, smoke/core/max tier, DSR lane, events, replay command, G4 drain/checkpoint schedule, G5 matrix, independent adjudication receipt, performance envelope, accessibility/agent-parity check and consuming gate",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i05-ir-admission-corpus",
            source: FixtureSource::AuthoredSpec {
                spec: "I05_FIXTURE_V1 DeploymentTwinIR corpus. VALID: empty task graph with explicit \
                       no-output semantics; single/multirate motor-control graphs; bounded state \
                       machines; typed ADC/PWM/CAN ports; periodic/sporadic tasks; declared faults and \
                       safe states. INVALID CLASSES: dangling ids, duplicate authority ids, unit/clock \
                       mismatch, algebraic task cycle without solver, unbounded recursion/loop/allocation, \
                       missing range/rounding/overflow, implicit float/cast, invalid schedule, overlapping \
                       memory, inaccessible safe state, nondeterministic tie, undeclared capability, and \
                       target-resource overflow. Independent generator emits validity witnesses and one \
                       minimal counterexample class tag per invalid case.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i05-fixed-point-exact-corpus",
            source: FixtureSource::AuthoredSpec {
                spec: "I05_FIXTURE_V1 exact-rational fixed-point microprograms. WIDTHS={2,3,4,7,8,15,16, \
                       24,31,32,63}; signed/unsigned; binary and rational scales; ties under nearest-even, \
                       toward-zero, floor and ceiling; saturating and trapping add/sub/mul/FMA/div; modular \
                       wrap only in explicitly modular programs; accumulator/cast edges; affine recurrences \
                       and lookup interpolation. Enumerable cases exhaust every input; larger cases carry \
                       rational interval witnesses. Includes negative zero source normalization, min/-1 \
                       division, shift boundaries, intermediate-width overflow and double-rounding traps.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i05-synthetic-target-profiles",
            source: FixtureSource::AuthoredSpec {
                spec: "I05_TARGET_PROFILE_V1 fully public synthetic targets, not vendor-silicon claims: \
                       T0=rv32imc-like in-order 3-stage, little-endian, no cache, 64 KiB SRAM, fixed \
                       instruction/memory costs, priority interrupts; T1=cortex-m7-like dual-issue model, \
                       little-endian, I/D cache sets/ways/replacement explicitly enumerated, ITCM/DTCM, \
                       flash wait states, branch predictor, DMA and bus arbitration finite automata. Every \
                       opcode/pipeline/memory/interrupt/reset transition either has executable semantics or \
                       is Unsupported. These profiles validate tooling only; hardware authority requires the \
                       waived independently characterized pack.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i05-motor-controller-suite",
            source: FixtureSource::AuthoredSpec {
                spec: "I05_FIXTURE_V1 sampled-data motor controls: PI current loop with anti-windup and \
                       Clarke/Park/SVPWM; six-step BLDC commutation; speed loop cascade; safe torque off; \
                       rates {10,20,40} kHz inner and {0.5,1,2} kHz outer; ADC/PWM delays, dead time, \
                       saturation and encoder wrap explicit. Plant family is bounded PMSM/BLDC dq plus \
                       inverter delay; initial/input/noise sets and 2 s horizon frozen. Analytic equilibria \
                       and high-precision floating traces provide development checks, not physical validation.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i05-reproducible-build-matrix",
            source: FixtureSource::AuthoredSpec {
                spec: "I05_FIXTURE_V1 hermetic build matrix over both synthetic profiles and controller \
                       suite: clean directories, path/order/locale/timezone perturbations, worker counts \
                       {1,2,7}, identical pinned toolchain/container bytes and SOURCE_DATE_EPOCH=0. Compare \
                       emitted source, metadata-normalized and raw objects, link map, final image, symbols, \
                       sections, relocations and capability inventory. A declared nondeterministic tool \
                       field must be removed at source; golden normalization cannot hide semantic bytes.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i05-timing-task-graphs",
            source: FixtureSource::AuthoredSpec {
                spec: "I05_FIXTURE_V1 enumerable timing graphs: periodic/sporadic fixed-priority and cyclic \
                       executive tasks; nested critical sections; ISR release; DMA/cache/bus interference; \
                       end-to-end sensor-to-actuator chains. Task counts 1..8, integer cycle costs, periods \
                       50..10000 cycles, deadlines constrained/arbitrary, jitter/blocking/preemption costs \
                       explicit. Small graphs exhaust all release/interference schedules; larger graphs have \
                       independently derived response-time upper bounds and trace lower bounds.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i05-timing-adversaries-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I05_FIXTURE_V1 CORE HOLDOUT timing adversaries generated after candidate freeze from \
                       indices 65536..=69631: cache conflict aliases, branch-history worst paths, flash/TCM \
                       boundary placement, interrupt phasing, DMA bursts, bus contention, priority inversion, \
                       wraparound counters and stale-binary/profile substitutions. One I05.G3 consumer; raw \
                       schedules and labels withheld until final one-shot core adjudication.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i05-compositional-timing-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I05_FIXTURE_V1 MAX HOLDOUT compositional timing graphs, indices \
                       131072..=135167: task/resource ownership joins with correlated releases, nested \
                       critical sections, ISR storms, DMA/cache/bus interference, mixed clock domains, \
                       age-of-data chains and authenticated trace partitions. Small projections have exact \
                       schedule-enumeration maxima; larger cases retain independently generated violating \
                       traces or conservative bounds. Includes stale identity, missing/double-counted owner, \
                       invalid optional-stopping and measurement-as-WCET adversaries. One I05.G6 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i05-hil-clock-faults-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I05_FIXTURE_V1 CORE HOLDOUT multi-clock HIL campaigns: independent host/DAQ/bus/device \
                       oscillators with offset, affine drift, piecewise thermal drift, tick wrap/reset, timestamp \
                       before/after transport, jitter, delay, drop, duplicate, reorder and burst catch-up. Logic- \
                       analyzer edges provide external intervals. Includes adversarial trace whose pointwise-best \
                       alignment passes but whole-interval directed refinement fails. Indices 69632..=73727.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i05-safe-state-faults-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I05_FIXTURE_V1 CORE HOLDOUT fault automata: sensor stuck/bias/noise/saturation/disconnect; \
                       actuator stuck/open/short/derate; CAN loss/corruption/reorder; task overrun/deadlock; clock \
                       jump/freeze; brownout/reset; watchdog and dual faults. Each case freezes injection fidelity, \
                       activation/detection/containment windows, permitted transient output, terminal safe state and \
                       external interlock action. Physical and software injection classes are never conflated. \
                       Indices 73728..=77823; one I05.G3 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i05-target-refinement-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I05_FIXTURE_V1 CORE HOLDOUT target-refinement traces over the motor-controller suite: \
                       boundary initial states, quantization ties, saturation entry/exit, encoder wrap, mode changes, \
                       delay/jitter extremes and declared single faults. Reference and target share exact sampled- \
                       data/fault semantics. QoIs: current/speed/torque/energy/state/safety margins normalized by \
                       frozen directed bands; trace aligner evaluates the entire clock interval. Indices \
                       77824..=81919; one I05.G3 consumer.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i05-refinement-theorem-card",
            source: FixtureSource::AuthoredSpec {
                spec: "I05_THEOREM_CARD_V1 TARGET. Intended successor proposition quantifies over a finite horizon \
                       N, initial/input/noise/fault sets, sampled-data reference transition R, decoded binary target \
                       transition T, clock-delay relation C and observation relation O_epsilon, and proves every T \
                       trace simulates an R trace under O_epsilon with invariant/safety preservation. REQUIRED \
                       MACHINE ARTIFACTS before proof: proposition and definition AST bytes/digests, total decoder \
                       semantics, premise schema, deterministic Lean translation and structural round-trip, exact \
                       axiom allowlist {propext,Quot.sound,Classical.choice}, transitive axiom closure, explicit rejection \
                       of sorryAx/custom postulates/native or unsafe oracle shortcuts, and retained kernel replay. This prose card grants \
                       no theorem authority.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i05-finite-target-wcet-card",
            source: FixtureSource::AuthoredSpec {
                spec: "I05_EXHAUSTIVENESS_CARD_V1 TARGET. Intended successor covers exact binary x finite initial \
                       architectural/microarchitectural/input state x bounded interrupt/DMA schedules, applies the \
                       executable transition/cycle-cost relation, and returns maximum cost plus argmax trace. REQUIRED \
                       MACHINE ARTIFACTS before search: complete grammar and canonical encoding, validity/initial/ \
                       transition/cost predicates, unsupported-opcode closure, total enumeration or sound symbolic \
                       coverage and exclusion order, symmetry quotient with proof obligations, rank/unrank/sharding, \
                       independent decoder, source digests, preflight and Merkle completeness root. This prose card \
                       grants no exhaustive or silicon-WCET authority.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i05-certificate-mutants-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I05_FIXTURE_V1 MAX HOLDOUT certifier mutants, indices 135168..=147455: signedness/scale/rounding/ \
                       overflow changes; stale source/object/binary/profile binders; omitted task/ISR/DMA/cache/bus \
                       state; clock wrap and favorable timestamp selection; dropped hybrid mode/fault; weakened \
                       observation relation; false loop bound; cost off-by-one; unsound state merge/symmetry quotient; \
                       incomplete shard; target-semantic drift; disallowed theorem axiom. Each mutant has a separately \
                       retained exact-rational trace, schedule, hardware lower-bound, or finite countermodel witness. \
                       One I05.G6 consumer; identities/labels withheld until the maximal candidate is sealed.",
            },
            partition: Partition::HeldOut,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i05_obligations() -> Vec<ObligationRow> {
    const UNIT_CASES: &[&str] = &[
        "happy",
        "empty",
        "boundary",
        "max",
        "error",
        "unit-dimension",
        "tie-break",
        "cancellation",
        "migration",
    ];
    vec![
        ObligationRow {
            leaf: "i05-ir-numeric-admission",
            claims_covered: &[
                "i05-deployment-ir-admission",
                "i05-fixed-point-range-soundness",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: valid/malformed bounded DeploymentTwinIR and exact-rational fixed-point \
                 microprograms; predicates: independent grammar validity, target-capability closure, and \
                 rational range/error containment; laws: canonical roundtrip, admission/refusal totality, \
                 synthesis idempotence, range monotonicity, exact scale/unit covariance, no unreported \
                 overflow, and checker acceptance iff every certificate premise binds; shrink by IR node/ \
                 basic-block removal while preserving the witness; replay from frozen Philox identity",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i05-ir-admission-corpus",
                "i05-fixed-point-exact-corpus",
                "i05-synthetic-target-profiles",
            ],
            g3_relations: &[
                "canonical id relabeling and insertion-order invariance",
                "exact unit/scale covariance with correspondingly transformed range and error bounds",
                "widening an input range cannot shrink an accepted reachable-range enclosure",
                "independent exact-rational evaluation either stays contained or refutes the certificate",
            ],
            g4_schedule: "inject cancel/panic/allocation failure after decode, validation, range propagation, \
                          format search, checker pass and canonical assembly; poll at IR/basic-block tiles; \
                          request->drain->finalize publishes no partial admitted IR or format; checkpoint after \
                          canonical node ids and resume/fork must equal one-shot identity; retain typed refusal \
                          and bounded FailureBundle",
            g5_matrix: "logical shards {1,2,7} x input orders {forward,reverse,permuted} x deterministic mode x \
                        same toolchain/ISA; require byte-identical canonical IR, format assignment, certificates, \
                        event JSON and digest; cross-ISA runs compare semantic exact-rational records, not floats",
            entry_point: "scripts/e2e/leapfrog/i05_ir_numeric.sh",
            tier: CampaignTier::Smoke,
            dsr_lane: "dsr quality --tool frankensim (i05-ir-numeric slice)",
            obs_events: &[
                "request.received",
                "ir.admitted",
                "ir.refused",
                "numeric.synthesized",
                "numeric.refuted",
                "numeric.unknown",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i05_ir_numeric.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i05-reproducible-codegen",
            claims_covered: &["i05-reproducible-no-std-binary"],
            unit_cases: UNIT_CASES,
            g0: "generators: admitted controller/task graphs x synthetic profiles; predicates: generated source \
                 stays inside the safe bounded no_std subset, link/resource maps fit, capability inventory is an \
                 allowlisted subset; laws: codegen purity, clean-build byte identity, source/object/binary/ABI \
                 binder consistency, zero undeclared symbols/sections/instructions, deterministic diagnostic order; \
                 shrink by task/function/basic-block removal; replay uses frozen build-input Merkle root",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i05-motor-controller-suite",
                "i05-synthetic-target-profiles",
                "i05-reproducible-build-matrix",
            ],
            g3_relations: &[
                "build path, locale, timezone and directory enumeration changes preserve raw image bytes",
                "worker-count changes preserve source, object graph, link map and final image",
                "one-byte IR or target-profile semantic mutation moves every downstream binder",
                "undeclared import/section/instruction injection is refused by the independent inventory checker",
            ],
            g4_schedule: "cancel or fault between source emission, per-object compile, link, binary parse, inventory \
                          check and artifact publication; drain all compiler/linker children; atomic finalize publishes \
                          the complete bound artifact set or none; checkpoint only content-addressed completed objects; \
                          resume cannot reuse an object whose full build-input key differs; request->drain->finalize is \
                          the only publication lifecycle",
            g5_matrix: "clean directories {A,B} x workers {1,2,7} x input orders {forward,reverse} x deterministic \
                        mode for each exact target/toolchain fingerprint; raw source/object/link-map/image bytes, \
                        inventories, logs and artifact graph must match bitwise",
            entry_point: "scripts/e2e/leapfrog/i05_codegen.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr build frankensim --target darwin/arm64 (i05 hermetic cross-target fixture lane)",
            obs_events: &[
                "request.received",
                "codegen.source_emitted",
                "codegen.object_completed",
                "codegen.link_completed",
                "binary.inventory_checked",
                "artifact.published",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i05_codegen.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i05-timing-evidence",
            claims_covered: &["i05-timing-evidence-separation"],
            unit_cases: UNIT_CASES,
            g0: "generators: timing receipt kinds, binaries/profiles/task contexts, valid/stale/malformed assumption \
                 graphs and task micrographs; predicates: exact identity binding, dimension/cycle-clock conversion, \
                 upper-bound premise completeness; laws: kind disjointness, receipt canonicality, no measurement-to- \
                 upper-bound coercion, bound >= independently replayed trace cost when applicable, stale binder refusal; \
                 shrink to the minimal omitted assumption or distinguishing trace",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i05-synthetic-target-profiles",
                "i05-timing-task-graphs",
                "i05-compositional-timing-max-holdout",
            ],
            g3_relations: &[
                "cycle-to-second rescaling follows only the frozen clock interval and preserves kind",
                "binary/profile/task-context substitution is IntegrityFailed before numeric comparison",
                "duplicating or permuting measurements preserves sample maximum and never creates WCET authority",
                "adding a possible interference source cannot improve a sound upper bound without a discharged exclusion",
            ],
            g4_schedule: "cancel/fault during binary decode, CFG/loop analysis, target-state exploration, probe capture, \
                          receipt join and retention; completed static states and raw samples checkpoint separately by \
                          logical id; request->drain->finalize prevents a partial bound/sample set from becoming complete; \
                          timeout/budget/infrastructure/integrity states remain distinct",
            g5_matrix: "analysis shards {1,2,7} x sample ingestion orders {capture,reverse,permuted} x deterministic mode \
                        on each identical target/toolchain fingerprint; kind, assumptions, cycles, clock interval, sample \
                        multiset root, counterexample and canonical receipt must match bitwise",
            entry_point: "scripts/e2e/leapfrog/i05_timing.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i05-timing slice)",
            obs_events: &[
                "request.received",
                "timing.static_bound",
                "timing.compositional_bound",
                "timing.sample",
                "timing.unavailable",
                "timing.identity_failed",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i05_timing.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i05-hil-fault-safety",
            claims_covered: &["i05-hil-clock-causality", "i05-safe-state-fault-semantics"],
            unit_cases: UNIT_CASES,
            g0: "generators: multi-clock maps plus typed sensor/actuator/compute/bus/clock/power/task faults and \
                 fixture fault automata; predicates: calibration containment, causal interval ordering, injection \
                 identity/fidelity, safe-output transition and external interlock agreement; laws: clock-wrap/reset \
                 totality, interval containment, no favorable alignment, fault terminal-state exclusivity, no normal \
                 publication after containment failure; shrink time/fault graph while retaining the violation",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i05-motor-controller-suite",
                "i05-synthetic-target-profiles",
                "i05-hil-clock-faults-core-holdout",
                "i05-safe-state-faults-core-holdout",
            ],
            g3_relations: &[
                "common clock-origin translation preserves causal intervals and deadline verdicts",
                "clock-unit rescaling with exact rational tick conversion preserves physical intervals",
                "external-id relabeling and record-order permutation preserve the canonical fault transcript",
                "widening timestamp uncertainty cannot turn TimingUnknown or late into proven on-time",
                "software injection cannot acquire electrical/physical injection-fidelity authority",
            ],
            g4_schedule: "arm independent emergency stop, then inject cancel/panic/timeout and every frozen fault at \
                          pre-activation, active, detection, containment and recovery boundaries; poll per HIL sample/task \
                          release; request->drain commands declared safe output, drains I/O, finalizes terminal receipt and \
                          never auto-retries physical faults; request->drain->finalize is the only publication lifecycle; \
                          checkpoint/resume uses a fresh attempt id and cannot repeat \
                          an actuation without explicit idempotency authorization",
            g5_matrix: "logical ingestion shards {1,2,7} x record orders {capture,reverse,permuted} x deterministic \
                        mode offline replay for identical rig/clock/probe fingerprints; raw physical timing may vary and is \
                        retained, while canonical intervals, classifications and replay verdict must match",
            entry_point: "scripts/e2e/leapfrog/i05_hil_faults.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i05-hil-fault slice on admitted rig)",
            obs_events: &[
                "request.received",
                "clock.calibrated",
                "sample.causality_bounded",
                "fault.activated",
                "fault.detected",
                "safe_state.entered",
                "containment.failed",
                "emergency_stop.activated",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i05_hil_faults.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i05-target-refinement",
            claims_covered: &["i05-bounded-target-refinement"],
            unit_cases: UNIT_CASES,
            g0: "generators: bounded controller/plant/input/noise/fault traces with exact reference and target arithmetic; \
                 predicates: shared sampled-data semantics, clock-interval alignment, mode equality and directed QoI/safety \
                 interval acceptance; laws: reference self-refines at zero error, normalized bands are unit-rescaling \
                 invariant, target emulator/physical transcript agreement stays separately typed, worst error is monotone \
                 under fixture union; shrink to earliest divergent tick and smallest active state/input/fault slice",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i05-fixed-point-exact-corpus",
                "i05-motor-controller-suite",
                "i05-target-refinement-core-holdout",
            ],
            g3_relations: &[
                "consistent physical-unit rescaling preserves normalized directed-error verdict",
                "refining clock uncertainty can preserve or improve but widening cannot improve acceptance",
                "held-out trace reorder/id relabel leaves per-trace result and worst normalized error unchanged",
                "one-bit arithmetic or one-tick schedule mutation must move identity and be exposed if it changes a trace",
            ],
            g4_schedule: "cancel/panic/timeout at reference tick, target tick, mode transition, clock-alignment and QoI \
                          reduction boundaries; request->drain->finalize publishes no accepted aggregate without every \
                          scheduled trace terminal state; checkpoints bind completed trace/tick prefixes and arithmetic/ \
                          target digests; resume/fork equals one-shot or is IntegrityFailed",
            g5_matrix: "trace shards {1,2,7} x worker counts {1,2,7} x input orders {forward,reverse,permuted} x \
                        deterministic mode under the same ISA/target-emulator fingerprint; exact rational target state, \
                        mode trace, interval alignment, QoIs, worst-case tie break, logs and receipt match bitwise",
            entry_point: "scripts/e2e/leapfrog/i05_refinement.sh",
            tier: CampaignTier::Core,
            dsr_lane: "dsr quality --tool frankensim (i05-refinement slice)",
            obs_events: &[
                "request.received",
                "reference.tick",
                "target.tick",
                "trace.aligned",
                "refinement.accepted",
                "refinement.refuted",
                "refinement.unknown",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i05_refinement.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i05-compositional-timing",
            claims_covered: &[
                "i05-compositional-deadline-contracts",
                "i05-static-measurement-reconciliation",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: bounded task/resource/interference graphs, upper-bound receipts and authenticated measured \
                 traces; predicates: ownership completeness/non-overlap, assumption compatibility, schedule enumeration, \
                 bound applicability and frozen anytime-valid diagnostic; laws: composition contains exhaustive small-graph \
                 maximum, observed trace <= applicable bound, kind separation, interference monotonicity, deterministic \
                 owner attribution; shrink to minimal task/resource/trace causing bound or identity failure",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i05-synthetic-target-profiles",
                "i05-timing-task-graphs",
                "i05-timing-adversaries-core-holdout",
            ],
            g3_relations: &[
                "independent schedule enumeration falsifies but never manufactures an upper bound",
                "task-id permutation conjugates the schedule and preserves chain response/age bounds",
                "adding release/interference possibilities cannot improve a valid worst-case bound",
                "sample duplication/order changes do not alter maxima or anytime-valid evidence under logical ids",
                "profile/binary/workload mismatch fails integrity before reconciliation",
            ],
            g4_schedule: "cancel/fault during local analysis, contract join, schedule search, trace ingestion and \
                          sequential diagnostic; drain child analyses and persist completed logical ids/e-process state; \
                          timeout/budget leaves DeadlineUnknown; an exceedance finalizes Refuted before any retry; resume \
                          cannot double-count a sample or change the frozen stopping rule; request->drain->finalize and the \
                          checkpoint state are receipt-bound",
            g5_matrix: "analysis/trace shards {1,2,7} x worker counts {1,2,7} x ingestion orders {forward,reverse} x \
                        deterministic mode on identical fingerprints; owner map, applicable assumptions, bounds, trace roots, \
                        e-process path, tie breaks, terminal states and receipt bytes match",
            entry_point: "scripts/e2e/leapfrog/i05_compositional_timing.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i05-compositional-timing slice)",
            obs_events: &[
                "request.received",
                "timing.contract_local",
                "timing.contract_composed",
                "timing.bound_refuted",
                "timing.coverage_gap",
                "timing.eprocess_updated",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i05_compositional_timing.sh --manifest <manifest-id> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i05-maximal-certifiers",
            claims_covered: &[
                "i05-machine-checked-target-refinement",
                "i05-complete-finite-target-wcet",
                "i05-false-equivalence-certificate-falsifier",
            ],
            unit_cases: UNIT_CASES,
            g0: "generators: proposition/definition/target-semantics ASTs, runtime premise records, finite target states/ \
                 transitions/costs and certifier mutants; predicates: canonical decode/roundtrip, total premise binding, \
                 axiom closure, finite grammar validity, coverage/root verification, quotient obligations, maximizing trace \
                 replay and countermodel witness; laws: kernel acceptance only under allowlist, rank/unrank bijection where \
                 used, shard disjoint union completeness, maximum >= every replayed cost, every known mutant refused; shrink \
                 proof term/state grammar/program/mutant while preserving the failure",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                "i05-synthetic-target-profiles",
                "i05-refinement-theorem-card",
                "i05-finite-target-wcet-card",
                "i05-certificate-mutants-max-holdout",
                "i05-industrial-target-pack",
            ],
            g3_relations: &[
                "AST alpha-renaming with canonical rebinding preserves proposition semantics and proof verdict",
                "state-id permutation conjugates transitions/costs and preserves exact maximum/argmax equivalence class",
                "shard permutation and worker count preserve completeness root and maximizing trace tie break",
                "adding valid initial states or transitions cannot reduce an exact worst-case maximum",
                "disallowed axiom, omitted premise/state/shard, stale binary/profile or unsound quotient mutation fails closed",
                "every accepted false-certificate mutant is retained as a certifier refutation",
            ],
            g4_schedule: "whole-campaign preflight refuses impossible budgets before held-out reveal; then inject cancel/ \
                          panic/timeout/corrupt-checkpoint at AST translation, kernel replay, state expansion, quotient, shard, \
                          Merkle join, maximizing replay and mutant adjudication boundaries; drain all proof/search children; \
                          checkpoints bind immutable frontier/visited/completeness roots; BudgetExhausted stays Unknown and \
                          cannot publish theorem/exhaustive authority; request->drain->finalize governs every proof/search \
                          terminal state",
            g5_matrix: "proof/search shards {1,2,7,31} x workers {1,2,7} x frontier orders {lex,reverse,permuted} x \
                        deterministic mode on identical toolchain/target-semantics fingerprint; AST bytes, premise/axiom \
                        closure, kernel result, explored/valid/excluded counts, roots, exact maximum/argmax, mutant verdicts, \
                        terminal states and receipts match bitwise",
            entry_point: "scripts/e2e/leapfrog/i05_maximal_certifiers.sh",
            tier: CampaignTier::Max,
            dsr_lane: "dsr quality --tool frankensim (i05-maximal-certifiers isolated lane)",
            obs_events: &[
                "request.received",
                "theorem.ast_checked",
                "theorem.premise_bound",
                "theorem.kernel_checked",
                "wcet.state_expanded",
                "wcet.shard_completed",
                "wcet.completeness_checked",
                "wcet.maximum_replayed",
                "certifier.mutant_refused",
                "certifier.refuted",
                "cancel.requested",
                "drain.completed",
                "finalize.completed",
                "failure_bundle.retained",
                "adjudication.receipt",
            ],
            replay_command: "scripts/e2e/leapfrog/i05_maximal_certifiers.sh --manifest <manifest-id> --replay <artifact-id>",
        },
    ]
}

fn i05_waivers() -> Vec<Waiver> {
    vec![Waiver {
        subject: "i05-industrial-target-pack",
        reason: "the maximal finite-target/WCET lane currently has only fully public synthetic target \
                 semantics; cycle-accurate silicon characterization, errata, proprietary bus/cache details, \
                 calibrated probes and independent lab traces for a production MCU are not yet licensed and \
                 content-addressed for third-party replay",
        owner: "I05 target-characterization and V&V registry owners",
        predicate: "an independently governed production-target pack is admitted with exact silicon/board/ \
                    toolchain/probe identities, executable supported-state semantics, errata coverage, model- \
                    discrepancy qualification, raw authorized traces and a replayable retention receipt",
        expiry: "before the first I05 maximal campaign submission; review at every manifest amendment and \
                 immediately on silicon, board, compiler, linker, clock/power or errata revision",
        promotion_effect: "machine-checked model refinement may be attempted on synthetic semantics, but no \
                           silicon WCET, industrial HIL equivalence, safety certification or maximal I05 \
                           promotion may be claimed while this waiver is live; baseline synthetic/core \
                           engineering evidence remains separately eligible",
    }]
}
