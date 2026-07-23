//! fs-obs — structured observability: the ONE event schema for kernels,
//! solvers, test suites, and (once it lands) the ledger `events` table.
//!
//! Roughly forty beads specify per-suite logging ("structured JSON records
//! sufficient to reproduce any failure from logs alone"). Without a single
//! owned schema, forty suites invent forty dialects and that promise dies.
//! This crate makes "diagnosable from logs alone" a CHECKABLE property:
//! every emitter produces [`Event`]s, every event serializes to one JSON-line
//! dialect, and the [`validate_line`] / [`lint_failure_record`] gates run in
//! CI.
//!
//! Determinism split (Decalogue P2): an event's CONTENT (kind, payload,
//! scope, seq) is deterministic in deterministic mode and is what
//! [`Event::content_hash`] covers; wall-clock time lives in the envelope
//! only, excluded from the hash, so logs from two runs of the same seed diff
//! cleanly.
//!
//! Serialization is in-house (Decalogue P1: std + constellation only — serde
//! is not on that list). The wire format is JSON-lines with CANONICAL field
//! order; the strict validator treats deviation as corruption, not dialect.
//!
//! The [`ident`] module owns the CANONICAL REPLAY IDENTITY encoding
//! (bead gp3.14): versioned, typed, length-prefixed — the shared
//! replacement for ad hoc delimiter-concatenation identities.

pub mod containment;
pub mod ident;
pub mod privacy;
pub mod process;

use core::fmt;
use core::num::NonZeroU64;
use std::fmt::Write as _;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Schema version stamped into every line; bump on any non-additive change.
pub const SCHEMA_VERSION: u32 = 1;

/// Exact typed event-content identity semantics. Version 2 replaced the
/// legacy display-JSON hash, which collapsed distinct NaN payload bits;
/// version 10 adds exact bounded-process frame and gap projections.
pub const EVENT_CONTENT_IDENTITY_VERSION: u32 = 10;

/// Domain-separated artifact kind framed into the typed event identity.
pub const EVENT_CONTENT_IDENTITY_DOMAIN: &str = "org.frankensim.fs-obs.event-content.v10";

/// Owner-local declaration consumed by `xtask check-identities`.
pub const EVENT_CONTENT_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-obs:event-content",
    "version_const=EVENT_CONTENT_IDENTITY_VERSION",
    "version=10",
    "domain=org.frankensim.fs-obs.event-content.v10",
    "domain_const=EVENT_CONTENT_IDENTITY_DOMAIN",
    "encoder=Event::content_identity",
    "encoder_helpers=Event::content_identity_with_versions,Event::content_identity_with_schema,bind_outcome_identity,bind_u128_identity,bind_process_artifact_identity,Severity::name,EventKind::kind_name,ReceiptScope::name,ExecutionDisposition::name,PredicateOutcome::name,EpistemicGrade::name,DomainApplicability::name,OperationalSupport::name,EvidenceCompleteness::name,EvidenceIntegrity::name,PromotionEffect::name,SubmissionLane::name,SubmissionOutcome::name,AttachAction::name,CapabilityDecision::name,LifecycleTransitionKind::name,ArtifactAction::name",
    "schema_constants=EVENT_CONTENT_IDENTITY_VERSION,EVENT_CONTENT_IDENTITY_DOMAIN,SCHEMA_VERSION",
    "schema_functions=check_event_content_identity_version,Event::content_identity_receipt,Event::admit_content_identity,fnv1a64,crates/fs-obs/src/process.rs#ProcessStream::name,crates/fs-obs/src/process.rs#ProcessEventClass::name,crates/fs-obs/src/process.rs#LossReason::name",
    "schema_dependencies=fs-obs:replay-identity-frame",
    "digest=fnv1a64",
    "encoding=typed-binary",
    "sources=Event,EventIdentityReceipt",
    "source_fields=Event.session:semantic,Event.scope:semantic,Event.seq:semantic,Event.severity:semantic,Event.kind:semantic,Event.wall_ns:nonsemantic:wall-clock-envelope-only,EventIdentityReceipt.declared_identity_version:semantic,EventIdentityReceipt.canonical_bytes:semantic,EventIdentityReceipt.root:derived:validated-fnv-root-of-retained-canonical-bytes",
    "source_bindings=Event.session>session,Event.scope>scope,Event.seq>seq,Event.severity>severity,Event.kind>kind+solver-residual-solver+solver-residual-iter+solver-residual-residual+tile-complete-tile+tile-complete-kernel+cancellation-reason+budget-delta-resource+budget-delta-spent+budget-delta-remaining+gradient-check-op+gradient-check-max-rel-err+gradient-check-pass+conformance-case-suite+conformance-case-case+conformance-case-pass+conformance-case-detail+conformance-case-seed+benchmark-result-kernel+benchmark-result-metric+benchmark-result-value+benchmark-result-machine+storm-assertion-name+storm-assertion-pass+storm-assertion-seed+manifest-selection-manifest+manifest-selection-source-snapshot+stratum-expansion-stratum+stratum-expansion-profile+stratum-expansion-cases+dsr-run-run+dsr-run-outcome-scope+dsr-run-outcome-receipt+dsr-run-outcome-disposition+dsr-run-outcome-predicate+dsr-run-outcome-evidence-methods+dsr-run-outcome-grade+dsr-run-outcome-applicability+dsr-run-outcome-support+dsr-run-outcome-completeness+dsr-run-outcome-integrity+dsr-run-outcome-promotion+dsr-run-outcome-detail+campaign-run-campaign+campaign-run-outcome-scope+campaign-run-outcome-receipt+campaign-run-outcome-disposition+campaign-run-outcome-predicate+campaign-run-outcome-evidence-methods+campaign-run-outcome-grade+campaign-run-outcome-applicability+campaign-run-outcome-support+campaign-run-outcome-completeness+campaign-run-outcome-integrity+campaign-run-outcome-promotion+campaign-run-outcome-detail+submission-decision-job+submission-decision-lane+submission-decision-outcome+submission-decision-detail+work-identity-binding-job+work-identity-binding-attempt+work-identity-binding-operation+lease-cursor-queue+lease-cursor-lease+lease-cursor-holder+lease-cursor-cursor+attach-detach-observer+attach-detach-target+attach-detach-action+journey-phase-journey+journey-phase-phase+journey-phase-ordinal+scope-tile-progress-completed+scope-tile-progress-total+heartbeat-worker+heartbeat-beat+observation-gap-expected-seq+observation-gap-resumed-seq+observation-gap-reason+process-frame-process+process-frame-stream+process-frame-ordinal+process-frame-class+process-frame-payload+process-frame-original-bytes+process-frame-artifact-pointer+process-gap-process+process-gap-stream+process-gap-class+process-gap-first-ordinal+process-gap-last-ordinal+process-gap-original-count+process-gap-emitted-count+process-gap-dropped-count+process-gap-omitted-inline-bytes+process-gap-reason+process-gap-policy-version+process-gap-artifact-pointer+oracle-comparison-oracle+oracle-comparison-subject+oracle-comparison-metric+oracle-comparison-observed+oracle-comparison-tolerance+oracle-comparison-pass+oracle-comparison-detail+tolerance-derivation-quantity+tolerance-derivation-tolerance+tolerance-derivation-basis+claim-adjudication-claim+claim-adjudication-outcome-scope+claim-adjudication-outcome-receipt+claim-adjudication-outcome-disposition+claim-adjudication-outcome-predicate+claim-adjudication-outcome-evidence-methods+claim-adjudication-outcome-grade+claim-adjudication-outcome-applicability+claim-adjudication-outcome-support+claim-adjudication-outcome-completeness+claim-adjudication-outcome-integrity+claim-adjudication-outcome-promotion+claim-adjudication-outcome-detail+capability-domain-decision-capability+capability-domain-decision-domain+capability-domain-decision-decision+capability-domain-decision-detail+lifecycle-transition-entity+lifecycle-transition-transition+lifecycle-transition-detail+artifact-lifecycle-artifact+artifact-lifecycle-action+artifact-lifecycle-actor+artifact-lifecycle-detail+visualization-transform-view+visualization-transform-transform+visualization-transform-source+diagnostic-repair-subject+diagnostic-repair-action+diagnostic-repair-pass+diagnostic-repair-detail+containment-node-attempt+containment-node-role+containment-node-node+containment-node-parent-role+containment-node-parent+containment-node-seq+containment-node-dsr-run+containment-node-campaign-run+containment-node-shard+containment-node-journey+containment-node-case+containment-gap-attempt+containment-gap-node-role+containment-gap-node+containment-gap-missing-parent-role+containment-gap-missing-parent+race-record-resource+race-record-schedule+race-record-pass+race-record-seed+degradation-event-resource+degradation-event-limit+degradation-event-observed+degradation-event-action+import-receipt-format+import-receipt-artifact+import-receipt-accepted+import-receipt-detail+certificate-verdict-certificate+certificate-verdict-pass+certificate-verdict-bound+certificate-verdict-detail+custom-name+custom-json-exact-opaque-utf8,EventIdentityReceipt.declared_identity_version>retained-producer-version,EventIdentityReceipt.canonical_bytes>retained-canonical-bytes",
    "external_semantic_fields=artifact-domain,identity-version,wire-schema",
    "semantic_fields=artifact-domain,identity-version,wire-schema,session,scope,seq,severity,kind,solver-residual-solver,solver-residual-iter,solver-residual-residual,tile-complete-tile,tile-complete-kernel,cancellation-reason,budget-delta-resource,budget-delta-spent,budget-delta-remaining,gradient-check-op,gradient-check-max-rel-err,gradient-check-pass,conformance-case-suite,conformance-case-case,conformance-case-pass,conformance-case-detail,conformance-case-seed,benchmark-result-kernel,benchmark-result-metric,benchmark-result-value,benchmark-result-machine,storm-assertion-name,storm-assertion-pass,storm-assertion-seed,manifest-selection-manifest,manifest-selection-source-snapshot,stratum-expansion-stratum,stratum-expansion-profile,stratum-expansion-cases,dsr-run-run,dsr-run-outcome-scope,dsr-run-outcome-receipt,dsr-run-outcome-disposition,dsr-run-outcome-predicate,dsr-run-outcome-evidence-methods,dsr-run-outcome-grade,dsr-run-outcome-applicability,dsr-run-outcome-support,dsr-run-outcome-completeness,dsr-run-outcome-integrity,dsr-run-outcome-promotion,dsr-run-outcome-detail,campaign-run-campaign,campaign-run-outcome-scope,campaign-run-outcome-receipt,campaign-run-outcome-disposition,campaign-run-outcome-predicate,campaign-run-outcome-evidence-methods,campaign-run-outcome-grade,campaign-run-outcome-applicability,campaign-run-outcome-support,campaign-run-outcome-completeness,campaign-run-outcome-integrity,campaign-run-outcome-promotion,campaign-run-outcome-detail,submission-decision-job,submission-decision-lane,submission-decision-outcome,submission-decision-detail,work-identity-binding-job,work-identity-binding-attempt,work-identity-binding-operation,lease-cursor-queue,lease-cursor-lease,lease-cursor-holder,lease-cursor-cursor,attach-detach-observer,attach-detach-target,attach-detach-action,journey-phase-journey,journey-phase-phase,journey-phase-ordinal,scope-tile-progress-completed,scope-tile-progress-total,heartbeat-worker,heartbeat-beat,observation-gap-expected-seq,observation-gap-resumed-seq,observation-gap-reason,process-frame-process,process-frame-stream,process-frame-ordinal,process-frame-class,process-frame-payload,process-frame-original-bytes,process-frame-artifact-pointer,process-gap-process,process-gap-stream,process-gap-class,process-gap-first-ordinal,process-gap-last-ordinal,process-gap-original-count,process-gap-emitted-count,process-gap-dropped-count,process-gap-omitted-inline-bytes,process-gap-reason,process-gap-policy-version,process-gap-artifact-pointer,oracle-comparison-oracle,oracle-comparison-subject,oracle-comparison-metric,oracle-comparison-observed,oracle-comparison-tolerance,oracle-comparison-pass,oracle-comparison-detail,tolerance-derivation-quantity,tolerance-derivation-tolerance,tolerance-derivation-basis,claim-adjudication-claim,claim-adjudication-outcome-scope,claim-adjudication-outcome-receipt,claim-adjudication-outcome-disposition,claim-adjudication-outcome-predicate,claim-adjudication-outcome-evidence-methods,claim-adjudication-outcome-grade,claim-adjudication-outcome-applicability,claim-adjudication-outcome-support,claim-adjudication-outcome-completeness,claim-adjudication-outcome-integrity,claim-adjudication-outcome-promotion,claim-adjudication-outcome-detail,capability-domain-decision-capability,capability-domain-decision-domain,capability-domain-decision-decision,capability-domain-decision-detail,lifecycle-transition-entity,lifecycle-transition-transition,lifecycle-transition-detail,artifact-lifecycle-artifact,artifact-lifecycle-action,artifact-lifecycle-actor,artifact-lifecycle-detail,visualization-transform-view,visualization-transform-transform,visualization-transform-source,diagnostic-repair-subject,diagnostic-repair-action,diagnostic-repair-pass,diagnostic-repair-detail,containment-node-attempt,containment-node-role,containment-node-node,containment-node-parent-role,containment-node-parent,containment-node-seq,containment-node-dsr-run,containment-node-campaign-run,containment-node-shard,containment-node-journey,containment-node-case,containment-gap-attempt,containment-gap-node-role,containment-gap-node,containment-gap-missing-parent-role,containment-gap-missing-parent,race-record-resource,race-record-schedule,race-record-pass,race-record-seed,degradation-event-resource,degradation-event-limit,degradation-event-observed,degradation-event-action,import-receipt-format,import-receipt-artifact,import-receipt-accepted,import-receipt-detail,certificate-verdict-certificate,certificate-verdict-pass,certificate-verdict-bound,certificate-verdict-detail,custom-name,custom-json-exact-opaque-utf8,retained-producer-version,retained-canonical-bytes",
    "excluded_fields=to-jsonl:display-transport-only",
    "consumers=Event::content_hash,EventIdentityReceipt,Event::admit_content_identity,ledger-event-sinks,replay-comparison",
    "mutations=artifact-domain:crates/fs-obs/src/lib.rs#event_content_identity_domain_version_and_wire_schema_bytes_are_independent,identity-version:crates/fs-obs/src/lib.rs#event_content_identity_domain_version_and_wire_schema_bytes_are_independent,wire-schema:crates/fs-obs/src/lib.rs#event_content_identity_domain_version_and_wire_schema_bytes_are_independent,session:crates/fs-obs/src/lib.rs#event_content_identity_mutation_battery,scope:crates/fs-obs/src/lib.rs#event_content_identity_mutation_battery,seq:crates/fs-obs/src/lib.rs#event_content_identity_mutation_battery,severity:crates/fs-obs/src/lib.rs#event_content_identity_mutation_battery,kind:crates/fs-obs/src/lib.rs#event_content_identity_mutation_battery,solver-residual-solver:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,solver-residual-iter:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,solver-residual-residual:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,tile-complete-tile:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,tile-complete-kernel:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,cancellation-reason:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,budget-delta-resource:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,budget-delta-spent:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,budget-delta-remaining:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,gradient-check-op:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,gradient-check-max-rel-err:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,gradient-check-pass:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,conformance-case-suite:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,conformance-case-case:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,conformance-case-pass:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,conformance-case-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,conformance-case-seed:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,benchmark-result-kernel:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,benchmark-result-metric:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,benchmark-result-value:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,benchmark-result-machine:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,storm-assertion-name:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,storm-assertion-pass:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,storm-assertion-seed:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,manifest-selection-manifest:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,manifest-selection-source-snapshot:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,stratum-expansion-stratum:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,stratum-expansion-profile:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,stratum-expansion-cases:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-run:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-scope:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-receipt:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-disposition:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-predicate:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-evidence-methods:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-grade:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-applicability:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-support:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-completeness:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-integrity:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-promotion:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,dsr-run-outcome-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,campaign-run-campaign:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,campaign-run-outcome-scope:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,campaign-run-outcome-receipt:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,campaign-run-outcome-disposition:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,campaign-run-outcome-predicate:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,campaign-run-outcome-evidence-methods:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,campaign-run-outcome-grade:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,campaign-run-outcome-applicability:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,campaign-run-outcome-support:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,campaign-run-outcome-completeness:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,campaign-run-outcome-integrity:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,campaign-run-outcome-promotion:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,campaign-run-outcome-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,submission-decision-job:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,submission-decision-lane:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,submission-decision-outcome:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,submission-decision-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,work-identity-binding-job:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,work-identity-binding-attempt:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,work-identity-binding-operation:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,lease-cursor-queue:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,lease-cursor-lease:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,lease-cursor-holder:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,lease-cursor-cursor:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,attach-detach-observer:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,attach-detach-target:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,attach-detach-action:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,journey-phase-journey:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,journey-phase-phase:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,journey-phase-ordinal:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,scope-tile-progress-completed:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,scope-tile-progress-total:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,heartbeat-worker:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,heartbeat-beat:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,observation-gap-expected-seq:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,observation-gap-resumed-seq:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,observation-gap-reason:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-frame-process:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-frame-stream:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-frame-ordinal:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-frame-class:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-frame-payload:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-frame-original-bytes:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-frame-artifact-pointer:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-gap-process:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-gap-stream:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-gap-class:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-gap-first-ordinal:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-gap-last-ordinal:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-gap-original-count:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-gap-emitted-count:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-gap-dropped-count:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-gap-omitted-inline-bytes:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-gap-reason:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-gap-policy-version:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,process-gap-artifact-pointer:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,oracle-comparison-oracle:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,oracle-comparison-subject:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,oracle-comparison-metric:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,oracle-comparison-observed:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,oracle-comparison-tolerance:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,oracle-comparison-pass:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,oracle-comparison-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,tolerance-derivation-quantity:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,tolerance-derivation-tolerance:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,tolerance-derivation-basis:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,claim-adjudication-claim:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,claim-adjudication-outcome-scope:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,claim-adjudication-outcome-receipt:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,claim-adjudication-outcome-disposition:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,claim-adjudication-outcome-predicate:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,claim-adjudication-outcome-evidence-methods:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,claim-adjudication-outcome-grade:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,claim-adjudication-outcome-applicability:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,claim-adjudication-outcome-support:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,claim-adjudication-outcome-completeness:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,claim-adjudication-outcome-integrity:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,claim-adjudication-outcome-promotion:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,claim-adjudication-outcome-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,capability-domain-decision-capability:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,capability-domain-decision-domain:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,capability-domain-decision-decision:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,capability-domain-decision-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,lifecycle-transition-entity:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,lifecycle-transition-transition:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,lifecycle-transition-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,artifact-lifecycle-artifact:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,artifact-lifecycle-action:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,artifact-lifecycle-actor:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,artifact-lifecycle-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,visualization-transform-view:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,visualization-transform-transform:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,visualization-transform-source:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,diagnostic-repair-subject:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,diagnostic-repair-action:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,diagnostic-repair-pass:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,diagnostic-repair-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-node-attempt:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-node-role:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-node-node:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-node-parent-role:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-node-parent:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-node-seq:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-node-dsr-run:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-node-campaign-run:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-node-shard:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-node-journey:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-node-case:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-gap-attempt:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-gap-node-role:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-gap-node:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-gap-missing-parent-role:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,containment-gap-missing-parent:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,race-record-resource:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,race-record-schedule:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,race-record-pass:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,race-record-seed:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,degradation-event-resource:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,degradation-event-limit:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,degradation-event-observed:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,degradation-event-action:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,import-receipt-format:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,import-receipt-artifact:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,import-receipt-accepted:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,import-receipt-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,certificate-verdict-certificate:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,certificate-verdict-pass:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,certificate-verdict-bound:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,certificate-verdict-detail:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,custom-name:crates/fs-obs/src/lib.rs#every_event_kind_payload_field_moves_identity,custom-json-exact-opaque-utf8:crates/fs-obs/src/lib.rs#custom_payload_identity_is_exact_opaque_utf8,retained-producer-version:crates/fs-obs/src/lib.rs#retained_event_identity_receipts_admit_exactly_or_fail_closed,retained-canonical-bytes:crates/fs-obs/src/lib.rs#retained_event_identity_receipts_admit_exactly_or_fail_closed",
    "nonsemantic_mutations=Event.wall_ns:crates/fs-obs/src/lib.rs#wall_clock_is_envelope_only,to-jsonl:crates/fs-obs/src/lib.rs#content_identity_preserves_bits_that_display_json_collapses",
    "field_guard=classify_event_identity_fields",
    "transport_guard=Event::admit_content_identity",
    "version_guard=crates/fs-obs/src/lib.rs#event_content_identity_versions_fail_closed",
    "coupling_surface=fs-obs:event-content-identity",
];

/// Structured refusal for event identities produced under unknown semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EventIdentityVersionError {
    /// Version declared by retained evidence.
    pub declared: u32,
    /// Exact version supported by this build.
    pub supported: u32,
}

impl fmt::Display for EventIdentityVersionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "event content identity v{} is unsupported; this build accepts exactly v{}",
            self.declared, self.supported
        )
    }
}

impl core::error::Error for EventIdentityVersionError {}

/// Retained proof of one event's exact content identity.
///
/// A receipt is deliberately more than a naked root: it carries the declared
/// event-identity semantics and the complete canonical preimage. Consumers
/// must call [`Event::admit_content_identity`] before trusting retained data.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EventIdentityReceipt {
    declared_identity_version: u32,
    canonical_bytes: Vec<u8>,
    root: u64,
}

impl EventIdentityReceipt {
    /// Reconstruct a receipt loaded from retained evidence.
    ///
    /// This constructor intentionally does not bless the parts. Call
    /// [`Event::admit_content_identity`] to verify the declared version, root,
    /// exact bytes, and event content together.
    #[must_use]
    pub fn from_retained_parts(
        declared_identity_version: u32,
        canonical_bytes: Vec<u8>,
        root: u64,
    ) -> Self {
        Self {
            declared_identity_version,
            canonical_bytes,
            root,
        }
    }

    /// Event-content identity semantics declared by the retained producer.
    #[must_use]
    pub fn declared_identity_version(&self) -> u32 {
        self.declared_identity_version
    }

    /// Complete canonical typed preimage retained by the producer.
    #[must_use]
    pub fn canonical_bytes(&self) -> &[u8] {
        &self.canonical_bytes
    }

    /// Retained FNV-1a root over [`Self::canonical_bytes`].
    #[must_use]
    pub fn root(&self) -> u64 {
        self.root
    }
}

/// Fail-closed refusal when retained event-identity evidence is not exactly
/// the identity of the event being admitted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EventIdentityAdmissionError {
    /// The producer declared identity semantics unknown to this build.
    UnsupportedVersion(EventIdentityVersionError),
    /// The retained root is not derived from the retained canonical bytes.
    RootMismatch {
        /// Root recorded in the receipt.
        declared: u64,
        /// Root recomputed from the receipt's canonical bytes.
        computed: u64,
    },
    /// The self-consistent receipt names different exact event content.
    CanonicalBytesMismatch {
        /// Root recorded in the receipt.
        declared_root: u64,
        /// Root computed from the event supplied to admission.
        expected_root: u64,
    },
}

impl fmt::Display for EventIdentityAdmissionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedVersion(error) => fmt::Display::fmt(error, f),
            Self::RootMismatch { declared, computed } => write!(
                f,
                "retained event identity root {declared:016x} does not match canonical bytes root {computed:016x}"
            ),
            Self::CanonicalBytesMismatch {
                declared_root,
                expected_root,
            } => write!(
                f,
                "retained event identity {declared_root:016x} does not bind the admitted event identity {expected_root:016x}"
            ),
        }
    }
}

impl core::error::Error for EventIdentityAdmissionError {}

/// Refuse retained event identities whose exact typed semantics are unknown.
///
/// # Errors
/// [`EventIdentityVersionError`] for any version other than the current one.
pub fn check_event_content_identity_version(
    declared: u32,
) -> Result<(), EventIdentityVersionError> {
    if declared == EVENT_CONTENT_IDENTITY_VERSION {
        Ok(())
    } else {
        Err(EventIdentityVersionError {
            declared,
            supported: EVENT_CONTENT_IDENTITY_VERSION,
        })
    }
}

/// Receipt scope for one verification outcome (i94v.7.3.1): a successful
/// status query, a failed worker attempt, and a refuted scientific job can
/// coexist, so every outcome names the exact scope its receipt covers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReceiptScope {
    /// One scoped operation.
    Operation,
    /// One worker attempt.
    Attempt,
    /// One job across attempts.
    Job,
    /// One whole campaign.
    Campaign,
}

impl ReceiptScope {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Operation => "operation",
            Self::Attempt => "attempt",
            Self::Job => "job",
            Self::Campaign => "campaign",
        }
    }
}

/// How the scoped execution ended, independent of what it proved.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExecutionDisposition {
    /// Ran to completion.
    Completed,
    /// Ran and failed.
    Failed,
    /// Cancelled before completion.
    Cancelled,
    /// Exceeded its budget/deadline.
    TimedOut,
    /// Refused before execution.
    Refused,
    /// Ended in a state the harness cannot classify; closure incomplete.
    Indeterminate,
}

impl ExecutionDisposition {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Completed => "completed",
            Self::Failed => "failed",
            Self::Cancelled => "cancelled",
            Self::TimedOut => "timed_out",
            Self::Refused => "refused",
            Self::Indeterminate => "indeterminate",
        }
    }
}

/// What the requested predicate concluded, independent of how it ran.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PredicateOutcome {
    /// The predicate held.
    Satisfied,
    /// The predicate was refuted.
    Refuted,
    /// Ran but could not decide.
    Indeterminate,
    /// This outcome carries no predicate.
    NotApplicable,
}

impl PredicateOutcome {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Satisfied => "satisfied",
            Self::Refuted => "refuted",
            Self::Indeterminate => "indeterminate",
            Self::NotApplicable => "not_applicable",
        }
    }
}

/// Epistemic grade of the produced evidence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EpistemicGrade {
    /// Machine-checked enclosure/certificate.
    Verified,
    /// Validated against admitted external data.
    Validated,
    /// Estimated by an uncertified method.
    Estimated,
    /// Reported without evidence machinery.
    Reported,
    /// No evidence claim in this outcome.
    NotApplicable,
}

impl EpistemicGrade {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Verified => "verified",
            Self::Validated => "validated",
            Self::Estimated => "estimated",
            Self::Reported => "reported",
            Self::NotApplicable => "not_applicable",
        }
    }
}

/// Whether the run stayed inside its declared applicability domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DomainApplicability {
    /// In-domain.
    InDomain,
    /// Out of the declared domain.
    OutOfDomain,
    /// Applicability was not evaluated.
    Unevaluated,
    /// No domain claim applies.
    NotApplicable,
}

impl DomainApplicability {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::InDomain => "in_domain",
            Self::OutOfDomain => "out_of_domain",
            Self::Unevaluated => "unevaluated",
            Self::NotApplicable => "not_applicable",
        }
    }
}

/// Operational support level under which the outcome was produced.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OperationalSupport {
    /// Fully supported configuration.
    Supported,
    /// Degraded configuration; outcome stands with caveats.
    Degraded,
    /// Unsupported configuration.
    Unsupported,
    /// Support classification does not apply.
    NotApplicable,
}

impl OperationalSupport {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Supported => "supported",
            Self::Degraded => "degraded",
            Self::Unsupported => "unsupported",
            Self::NotApplicable => "not_applicable",
        }
    }
}

/// Whether the evidence set is complete (SEPARATE from integrity).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceCompleteness {
    /// Every required artifact/receipt is present.
    Complete,
    /// Required evidence is missing; closure incomplete.
    Incomplete,
    /// Completeness does not apply.
    NotApplicable,
}

impl EvidenceCompleteness {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Complete => "complete",
            Self::Incomplete => "incomplete",
            Self::NotApplicable => "not_applicable",
        }
    }
}

/// Whether the evidence set is intact (SEPARATE from completeness).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceIntegrity {
    /// All present evidence authenticated.
    Intact,
    /// Tamper/corruption detected.
    Compromised,
    /// Integrity not yet checked.
    Unchecked,
    /// Integrity does not apply.
    NotApplicable,
}

impl EvidenceIntegrity {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Intact => "intact",
            Self::Compromised => "compromised",
            Self::Unchecked => "unchecked",
            Self::NotApplicable => "not_applicable",
        }
    }
}

/// Effect of this outcome on promotion state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PromotionEffect {
    /// Promoted evidence/authority.
    Promoted,
    /// Demoted evidence/authority.
    Demoted,
    /// No promotion movement.
    Unchanged,
    /// Promotion does not apply.
    NotApplicable,
}

impl PromotionEffect {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Promoted => "promoted",
            Self::Demoted => "demoted",
            Self::Unchanged => "unchanged",
            Self::NotApplicable => "not_applicable",
        }
    }
}

/// Submission lane under the two-tier durable-submission doctrine
/// (i94v.7.3.1 F2): durable submissions are flushable and survive recovery;
/// plain submissions are in-memory only.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubmissionLane {
    /// Flushable, recovery-ledger-backed submission.
    Durable,
    /// In-memory submission with no recovery claim.
    Plain,
}

impl SubmissionLane {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Durable => "durable",
            Self::Plain => "plain",
        }
    }
}

/// Queue admission decision for one submission (i94v.7.3.1 F2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubmissionOutcome {
    /// Admitted into the queue.
    Admitted,
    /// Refused at admission.
    Refused,
    /// Deferred pending capacity or dependency.
    Deferred,
}

impl SubmissionOutcome {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Admitted => "admitted",
            Self::Refused => "refused",
            Self::Deferred => "deferred",
        }
    }
}

/// Observer attach/detach action on a running scope (i94v.7.3.1 F2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AttachAction {
    /// Observer attached.
    Attach,
    /// Observer detached.
    Detach,
}

impl AttachAction {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Attach => "attach",
            Self::Detach => "detach",
        }
    }
}

/// Capability admission decision for one declared domain of use
/// (i94v.7.3.1 F4).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CapabilityDecision {
    /// Capability admitted in the domain.
    Admitted,
    /// Capability refused in the domain.
    Refused,
    /// Capability admitted with restrictions.
    Restricted,
}

impl CapabilityDecision {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Admitted => "admitted",
            Self::Refused => "refused",
            Self::Restricted => "restricted",
        }
    }
}

/// Durable-lifecycle transition kind (i94v.7.3.1 F5).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LifecycleTransitionKind {
    /// A fault was contained.
    Fault,
    /// A checkpoint was taken.
    Checkpoint,
    /// Execution paused.
    Pause,
    /// Execution resumed.
    Resume,
    /// Work migrated to another host/worker.
    Migrate,
    /// The scope is draining (no new work admitted).
    Drain,
    /// The scope finalized.
    Finalize,
}

impl LifecycleTransitionKind {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Fault => "fault",
            Self::Checkpoint => "checkpoint",
            Self::Pause => "pause",
            Self::Resume => "resume",
            Self::Migrate => "migrate",
            Self::Drain => "drain",
            Self::Finalize => "finalize",
        }
    }
}

/// Artifact lifecycle action (i94v.7.3.1 F5).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArtifactAction {
    /// Artifact committed to retained storage.
    Commit,
    /// Retention decision applied.
    Retain,
    /// Content redacted.
    Redact,
    /// Artifact shared outward.
    Share,
    /// Artifact integrity verified.
    Verify,
}

impl ArtifactAction {
    /// Stable wire name.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Self::Commit => "commit",
            Self::Retain => "retain",
            Self::Redact => "redact",
            Self::Share => "share",
            Self::Verify => "verify",
        }
    }
}

/// One scoped verification outcome (i94v.7.3.1): the shared core every
/// outcome-bearing event kind embeds. Not-applicable states are explicit
/// variants, never absent-by-convention fields.
#[derive(Debug, Clone, PartialEq)]
pub struct ScopedReceiptOutcome {
    /// Scope the receipt covers.
    pub scope: ReceiptScope,
    /// Receipt identity (hex) at that scope.
    pub receipt: String,
    /// How execution ended.
    pub disposition: ExecutionDisposition,
    /// What the requested predicate concluded.
    pub predicate: PredicateOutcome,
    /// Comma-joined registered evidence-method names ("" = none).
    pub evidence_methods: String,
    /// Epistemic grade of the produced evidence.
    pub grade: EpistemicGrade,
    /// Domain-applicability classification.
    pub applicability: DomainApplicability,
    /// Operational-support classification.
    pub support: OperationalSupport,
    /// Evidence completeness (separate from integrity).
    pub completeness: EvidenceCompleteness,
    /// Evidence integrity (separate from completeness).
    pub integrity: EvidenceIntegrity,
    /// Promotion effect.
    pub promotion: PromotionEffect,
    /// Diagnostic detail; REQUIRED non-empty for failed/refused/
    /// indeterminate dispositions (lint-enforced).
    pub detail: String,
}

/// Severity ladder. `Error` events MUST satisfy [`lint_failure_record`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Severity {
    /// High-volume diagnostics (tile completions).
    Trace,
    /// Normal progress (solver iterations, case verdicts).
    Info,
    /// Something degraded but the run continues.
    Warn,
    /// A failure record — must be reproducible from its own payload.
    Error,
}

impl Severity {
    fn name(self) -> &'static str {
        match self {
            Severity::Trace => "trace",
            Severity::Info => "info",
            Severity::Warn => "warn",
            Severity::Error => "error",
        }
    }
}

/// Deterministic sampling policy for one severity's logical opportunities.
///
/// The representation is intentionally opaque: callers choose either
/// [`Self::never`] or [`Self::every`] and cannot construct a zero cadence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SamplingCadence(Option<NonZeroU64>);

impl SamplingCadence {
    /// Reject every logical opportunity.
    #[must_use]
    pub const fn never() -> Self {
        Self(None)
    }

    /// Admit logical ordinal zero and every `period`-th opportunity after it.
    #[must_use]
    pub const fn every(period: NonZeroU64) -> Self {
        Self(Some(period))
    }

    fn admits(self, logical_ordinal: u64) -> bool {
        self.0
            .is_some_and(|period| logical_ordinal.is_multiple_of(period.get()))
    }
}

/// Independent deterministic sampling policy for trace and info events.
///
/// Warn and error events bypass both cadences and are always admitted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EmissionGate {
    trace: SamplingCadence,
    info: SamplingCadence,
}

impl EmissionGate {
    /// Build a gate with explicit trace and info cadences.
    #[must_use]
    pub const fn new(trace: SamplingCadence, info: SamplingCadence) -> Self {
        Self { trace, info }
    }

    fn admits(self, severity: Severity, logical_ordinal: u64) -> bool {
        match severity {
            Severity::Trace => self.trace.admits(logical_ordinal),
            Severity::Info => self.info.admits(logical_ordinal),
            Severity::Warn | Severity::Error => true,
        }
    }
}

/// Typed payload registry, v1. Additive evolution only: new kinds may be
/// added; existing fields may never change meaning (validator-enforced by
/// the golden lines in the conformance suite).
#[derive(Debug, Clone, PartialEq)]
pub enum EventKind {
    /// One iteration of an iterative solver.
    SolverResidual {
        /// Solver identity (e.g. "cg", "p-mg").
        solver: String,
        /// Iteration index.
        iter: u64,
        /// Residual norm.
        residual: f64,
    },
    /// A tile finished executing.
    TileComplete {
        /// Logical tile identity (the SAME identity that keys RNG streams
        /// and deterministic reductions — plan §5.2).
        tile: u64,
        /// Kernel name.
        kernel: String,
    },
    /// A scope was cancelled.
    Cancellation {
        /// Why (kill-handle, budget, panic containment, ...).
        reason: String,
    },
    /// Budget accounting delta (P4: budgets first).
    BudgetDelta {
        /// Resource name ("wall_s", "mem_bytes", "core_s", "energy_j").
        resource: String,
        /// Amount spent in this step.
        spent: f64,
        /// Remaining grant.
        remaining: f64,
    },
    /// A gradient verification outcome (the merge-gate evidence, §8.7).
    GradientCheck {
        /// Operator under test.
        op: String,
        /// Max relative error across probed directions.
        max_rel_err: f64,
        /// Verdict.
        pass: bool,
    },
    /// One conformance-suite case verdict (plan §13.3).
    ConformanceCase {
        /// Suite id (e.g. "fs-qty/conformance").
        suite: String,
        /// Case id (e.g. "qty-001/0.12Pa*s").
        case: String,
        /// Verdict.
        pass: bool,
        /// Human/agent-readable detail; REQUIRED non-empty when `pass=false`
        /// (lint-enforced) so failures reproduce from the log alone.
        detail: String,
        /// Replay seed when the case is randomized (lint: required on fail
        /// for randomized cases; 0 = not randomized).
        seed: u64,
    },
    /// A performance measurement (roofline harness rows).
    BenchmarkResult {
        /// Kernel name.
        kernel: String,
        /// Metric name ("glups", "gflops", "bandwidth_gbs", "mrays").
        metric: String,
        /// Measured value.
        value: f64,
        /// Machine fingerprint hash (fs-substrate probe).
        machine: u64,
    },
    /// A chaos/storm assertion outcome (G4).
    StormAssertion {
        /// Assertion name ("no-arena-leak", "cancel-latency-p99").
        name: String,
        /// Verdict.
        pass: bool,
        /// Storm seed for replay.
        seed: u64,
    },
    /// Manifest + source-snapshot selection for one verification run
    /// (i94v.7.3.1 F1).
    ManifestSelection {
        /// Selected manifest identity (hex).
        manifest: String,
        /// Source-snapshot identity (hex) the manifest was resolved against.
        source_snapshot: String,
    },
    /// Core/Max stratum and campaign-profile expansion (i94v.7.3.1 F1).
    StratumExpansion {
        /// Stratum name ("core", "max").
        stratum: String,
        /// Campaign profile the stratum expanded under.
        profile: String,
        /// Number of expanded cases.
        cases: u64,
    },
    /// One DSR run outcome carrying the full scoped receipt core
    /// (i94v.7.3.1 F1).
    DsrRun {
        /// DSR run identity.
        run: String,
        /// The scoped verification outcome.
        outcome: ScopedReceiptOutcome,
    },
    /// One campaign-level outcome carrying the full scoped receipt core
    /// (i94v.7.3.1 F2).
    CampaignRun {
        /// Campaign identity.
        campaign: String,
        /// The scoped verification outcome.
        outcome: ScopedReceiptOutcome,
    },
    /// Queue admission decision for one submission under the two-tier
    /// durable-submission doctrine (i94v.7.3.1 F2).
    SubmissionDecision {
        /// Submitted job identity (hex).
        job: String,
        /// Durable (flushable) or plain (in-memory) lane.
        lane: SubmissionLane,
        /// Admission decision.
        outcome: SubmissionOutcome,
        /// Decision detail; REQUIRED non-empty when refused or deferred
        /// (lint-enforced) so queue refusals reproduce from the log alone.
        detail: String,
    },
    /// Identity-chain binding correlating one operation to its attempt and
    /// job (i94v.7.3.1 F2) so scoped receipts at different [`ReceiptScope`]s
    /// name the same work.
    WorkIdentityBinding {
        /// Job identity (hex).
        job: String,
        /// Attempt ordinal within the job (1-based).
        attempt: u64,
        /// Operation identity (hex) within the attempt.
        operation: String,
    },
    /// Queue lease-cursor progress for one worker lease (i94v.7.3.1 F2).
    LeaseCursor {
        /// Queue identity.
        queue: String,
        /// Lease identity (hex).
        lease: String,
        /// Lease holder (worker identity).
        holder: String,
        /// Cursor position within the leased range.
        cursor: u64,
    },
    /// Observer attach/detach on a running scope (i94v.7.3.1 F2).
    AttachDetach {
        /// Observer identity.
        observer: String,
        /// Scope or campaign the observer targets.
        target: String,
        /// Attach or detach.
        action: AttachAction,
    },
    /// One CORE user journey entering a named phase (i94v.7.3.1 F3). The
    /// envelope scope carries the journey's scope-tree position; the payload
    /// names the journey and phase so replays reconstruct the storyboard.
    JourneyPhase {
        /// Journey identity (e.g. "thermal-fatigue-gearbox").
        journey: String,
        /// Phase name within the journey (e.g. "mesh", "solve", "adjudicate").
        phase: String,
        /// Phase ordinal within the journey (1-based, monotone).
        ordinal: u64,
    },
    /// Tile-level progress within the emitting scope (i94v.7.3.1 F3).
    ScopeTileProgress {
        /// Tiles completed so far in this scope.
        completed: u64,
        /// Total tiles planned for this scope.
        total: u64,
    },
    /// Worker liveness heartbeat (i94v.7.3.1 F3).
    Heartbeat {
        /// Worker identity emitting the beat.
        worker: String,
        /// Monotone beat counter per worker.
        beat: u64,
    },
    /// An honest declaration that observation was lost for a span of
    /// sequence numbers (i94v.7.3.1 F3): readers must treat the gap as
    /// unknown execution, never as absence of events.
    ObservationGap {
        /// First sequence number that was expected but not observed.
        expected_seq: u64,
        /// Sequence number at which observation resumed.
        resumed_seq: u64,
        /// Why the gap happened; REQUIRED non-empty (lint-enforced).
        reason: String,
    },
    /// One bounded inline frame captured from a child process stream
    /// (i94v.7.3.3). Opaque bytes stay exact; presentation renders them as
    /// lowercase hexadecimal rather than guessing UTF-8.
    ProcessFrame {
        /// Validated child-process identity.
        process: process::ProcessId,
        /// Source stream.
        stream: process::ProcessStream,
        /// One-based ordinal within that stream.
        ordinal: u64,
        /// Loss/backpressure class.
        class: process::ProcessEventClass,
        /// Exact inline bytes retained by the bounded capture policy.
        payload: Vec<u8>,
        /// Byte count before any bounded inline truncation.
        original_bytes: u64,
        /// Already-committed full-detail artifact, when inline bytes are a
        /// bounded projection.
        artifact: Option<process::DurableArtifactPointer>,
    },
    /// One quantified omission range from bounded process capture
    /// (i94v.7.3.3). Counts and the exact policy version make every sampling,
    /// truncation, spill, or sink-loss decision replay-visible.
    ProcessGap {
        /// Validated child-process identity.
        process: process::ProcessId,
        /// Source stream.
        stream: process::ProcessStream,
        /// Loss/backpressure class.
        class: process::ProcessEventClass,
        /// First affected source ordinal.
        first_ordinal: u64,
        /// Last affected source ordinal.
        last_ordinal: u64,
        /// Frames presented in the range.
        original_count: u128,
        /// Frames still represented inline.
        emitted_count: u128,
        /// Whole frames omitted.
        dropped_count: u128,
        /// Bytes omitted from inline payloads.
        omitted_inline_bytes: u128,
        /// Stable omission reason.
        reason: process::LossReason,
        /// Exact capture-policy version.
        policy_version: u32,
        /// Already-committed retained detail, when present.
        artifact: Option<process::DurableArtifactPointer>,
    },
    /// One oracle comparison verdict (i94v.7.3.1 F4): a run output measured
    /// against an independent oracle under a named metric and tolerance.
    OracleComparison {
        /// Oracle identity (e.g. "mms-poisson", "nafems-le1").
        oracle: String,
        /// Subject artifact/receipt identity (hex) being compared.
        subject: String,
        /// Comparison metric name (e.g. "l2-rel", "max-abs").
        metric: String,
        /// Observed metric value.
        observed: f64,
        /// Tolerance the comparison ran under.
        tolerance: f64,
        /// Verdict.
        pass: bool,
        /// Detail; REQUIRED non-empty when `pass=false` (lint-enforced).
        detail: String,
    },
    /// How a tolerance was derived (i94v.7.3.1 F4): every tolerance names
    /// its basis so no comparison rests on an unexplained constant.
    ToleranceDerivation {
        /// Quantity the tolerance governs.
        quantity: String,
        /// Derived tolerance value.
        tolerance: f64,
        /// Derivation basis; REQUIRED non-empty (lint-enforced) — an
        /// unexplained tolerance is the silent-constant bug class.
        basis: String,
    },
    /// One independent claim adjudication carrying the full scoped receipt
    /// core (i94v.7.3.1 F4).
    ClaimAdjudication {
        /// Adjudicated claim identity (hex).
        claim: String,
        /// The scoped adjudication outcome.
        outcome: ScopedReceiptOutcome,
    },
    /// Capability admission decision in a declared domain of use
    /// (i94v.7.3.1 F4).
    CapabilityDomainDecision {
        /// Capability identity (e.g. "time-varying-flux").
        capability: String,
        /// Declared domain of use the decision covers.
        domain: String,
        /// Admission decision.
        decision: CapabilityDecision,
        /// Decision detail; REQUIRED non-empty when refused or restricted
        /// (lint-enforced).
        detail: String,
    },
    /// One durable-lifecycle transition (i94v.7.3.1 F5): faults, checkpoints,
    /// pause/resume, migration, drain, and finalization all record the same
    /// typed shape so recovery replays the exact lifecycle path.
    LifecycleTransition {
        /// Entity making the transition (job/campaign/worker identity).
        entity: String,
        /// Typed transition.
        transition: LifecycleTransitionKind,
        /// Transition detail; REQUIRED non-empty for faults (lint-enforced).
        detail: String,
    },
    /// One artifact lifecycle action (i94v.7.3.1 F5).
    ArtifactLifecycle {
        /// Artifact content identity (hex).
        artifact: String,
        /// Typed action.
        action: ArtifactAction,
        /// Acting principal (user/agent/policy identity).
        actor: String,
        /// Action detail; REQUIRED non-empty for redact/share
        /// (lint-enforced): redactions need a rationale, shares a target.
        detail: String,
    },
    /// A visualization derived from retained data (i94v.7.3.1 F5): every
    /// rendered view names its source and transform so plots trace to data.
    VisualizationTransform {
        /// View identity (e.g. "stress-contour-p1").
        view: String,
        /// Named transform applied (e.g. "warp-by-displacement x50").
        transform: String,
        /// Source artifact/receipt identity (hex).
        source: String,
    },
    /// A diagnostic-driven repair attempt (i94v.7.3.1 F5).
    DiagnosticRepair {
        /// Subject that was diagnosed (artifact/config/receipt identity).
        subject: String,
        /// Repair action applied.
        action: String,
        /// Whether the repair verified as effective.
        pass: bool,
        /// Detail; REQUIRED non-empty when `pass=false` (lint-enforced).
        detail: String,
    },
    /// One node of a sealed local execution-containment tree
    /// (i94v.7.3.2; see [`containment`]). Role strings come from
    /// [`containment::LocalNodeId::role`]; empty context strings mean the
    /// edge is absent.
    ContainmentNode {
        /// Propagated Attempt root the tree is sealed under.
        attempt: String,
        /// Node role ("op", "scope", "tile").
        role: String,
        /// Node raw identity.
        node: String,
        /// Primary parent role ("attempt", "op", "scope", "tile").
        parent_role: String,
        /// Primary parent raw identity.
        parent: String,
        /// Deterministic sibling sequence.
        seq: u64,
        /// Hosting DSR invocation ("" = none).
        dsr_run: String,
        /// Selecting campaign execution ("" = none).
        campaign_run: String,
        /// Deterministic shard membership ("" = none).
        shard: String,
        /// Logical journey definition ("" = none).
        journey: String,
        /// Logical case definition ("" = none).
        case: String,
    },
    /// One explicit lineage gap in a sealed containment tree
    /// (i94v.7.3.2): the parent never arrived, so closure is incomplete.
    ContainmentGap {
        /// Propagated Attempt root the tree is sealed under.
        attempt: String,
        /// Role of the node whose lineage is incomplete.
        node_role: String,
        /// Raw identity of that node.
        node: String,
        /// Role of the parent that never arrived.
        missing_parent_role: String,
        /// Raw identity of the parent that never arrived.
        missing_parent: String,
    },
    /// A concurrency-race observation from a G4 race/storm harness.
    RaceRecord {
        /// Contended resource or invariant ("arena-slot", "tune-row", ...).
        resource: String,
        /// Interleaving/schedule descriptor sufficient to replay the race.
        schedule: String,
        /// Verdict.
        pass: bool,
        /// Replay seed for the schedule (lint: required non-zero on fail).
        seed: u64,
    },
    /// A graceful-degradation decision at a bounded resource.
    DegradationEvent {
        /// Resource that hit its bound ("retained_bytes_per_scope", ...).
        resource: String,
        /// Bound in force.
        limit: u64,
        /// Observed demand when the bound fired.
        observed: u64,
        /// Action taken ("rotate-lane", "flush", "refuse", ...).
        action: String,
    },
    /// An external-artifact import decision (STEP/mesh/dataset quarantine).
    ImportReceipt {
        /// Import format ("step", "gltf", "csv-dataset", ...).
        format: String,
        /// Content hash of the imported artifact (hex).
        artifact: String,
        /// Whether the importer admitted the artifact.
        accepted: bool,
        /// Refusal/quarantine detail; REQUIRED non-empty when
        /// `accepted=false` (lint-enforced).
        detail: String,
    },
    /// A certificate check outcome (enclosure/trim/conversion verdicts).
    CertificateVerdict {
        /// Certificate family ("ivl-enclosure", "trim-cert", ...).
        certificate: String,
        /// Verdict.
        pass: bool,
        /// Certified numeric bound (0.0 when the certificate is non-numeric).
        bound: f64,
        /// Detail; REQUIRED non-empty when `pass=false` (lint-enforced).
        detail: String,
    },
    /// Escape hatch for kinds not yet in the registry. `json` is exact opaque
    /// UTF-8: whitespace and object-member order are semantic identity bytes;
    /// fs-obs never claims to canonicalize unchecked JSON. The caller must
    /// still supply one valid pre-serialized JSON object for `to_jsonl`.
    Custom {
        /// Kind name (kebab-case).
        name: String,
        /// Exact opaque UTF-8 bytes of the pre-serialized JSON object.
        json: String,
    },
}

impl EventKind {
    /// Stable wire/ledger kind name. This is the documented export surface
    /// for ledger `events`-table sinks (huq.16): sinks store
    /// (`Event::to_jsonl` bytes, `kind_name`) without re-encoding.
    #[must_use]
    pub fn kind_name(&self) -> &'static str {
        match self {
            EventKind::SolverResidual { .. } => "solver_residual",
            EventKind::TileComplete { .. } => "tile_complete",
            EventKind::Cancellation { .. } => "cancellation",
            EventKind::BudgetDelta { .. } => "budget_delta",
            EventKind::GradientCheck { .. } => "gradient_check",
            EventKind::ConformanceCase { .. } => "conformance_case",
            EventKind::BenchmarkResult { .. } => "benchmark_result",
            EventKind::StormAssertion { .. } => "storm_assertion",
            EventKind::ManifestSelection { .. } => "manifest_selection",
            EventKind::StratumExpansion { .. } => "stratum_expansion",
            EventKind::DsrRun { .. } => "dsr_run",
            EventKind::CampaignRun { .. } => "campaign_run",
            EventKind::SubmissionDecision { .. } => "submission_decision",
            EventKind::WorkIdentityBinding { .. } => "work_identity_binding",
            EventKind::LeaseCursor { .. } => "lease_cursor",
            EventKind::AttachDetach { .. } => "attach_detach",
            EventKind::JourneyPhase { .. } => "journey_phase",
            EventKind::ScopeTileProgress { .. } => "scope_tile_progress",
            EventKind::Heartbeat { .. } => "heartbeat",
            EventKind::ObservationGap { .. } => "observation_gap",
            EventKind::ProcessFrame { .. } => "process_frame",
            EventKind::ProcessGap { .. } => "process_gap",
            EventKind::OracleComparison { .. } => "oracle_comparison",
            EventKind::ToleranceDerivation { .. } => "tolerance_derivation",
            EventKind::ClaimAdjudication { .. } => "claim_adjudication",
            EventKind::CapabilityDomainDecision { .. } => "capability_domain_decision",
            EventKind::LifecycleTransition { .. } => "lifecycle_transition",
            EventKind::ArtifactLifecycle { .. } => "artifact_lifecycle",
            EventKind::VisualizationTransform { .. } => "visualization_transform",
            EventKind::DiagnosticRepair { .. } => "diagnostic_repair",
            EventKind::ContainmentNode { .. } => "containment_node",
            EventKind::ContainmentGap { .. } => "containment_gap",
            EventKind::RaceRecord { .. } => "race_record",
            EventKind::DegradationEvent { .. } => "degradation_event",
            EventKind::ImportReceipt { .. } => "import_receipt",
            EventKind::CertificateVerdict { .. } => "certificate_verdict",
            EventKind::Custom { .. } => "custom",
        }
    }
}

/// One observability event: envelope + typed payload.
#[derive(Debug, Clone, PartialEq)]
pub struct Event {
    /// Session identity (study/run scope; content-hashed).
    pub session: String,
    /// Slash-separated scope path mirroring the asupersync scope tree
    /// (e.g. "study-x/op-3/kernel-lbm/tile-42"); content-hashed.
    pub scope: String,
    /// Per-emitter monotone sequence number; content-hashed (gives logs a
    /// deterministic total order per scope without wall-clock).
    pub seq: u64,
    /// Severity.
    pub severity: Severity,
    /// Typed payload.
    pub kind: EventKind,
    /// Wall-clock nanoseconds since the unix epoch. ENVELOPE ONLY: excluded
    /// from `content_hash`, always serialized LAST, `None` in deterministic
    /// replay comparisons.
    pub wall_ns: Option<u64>,
}

#[allow(dead_code)]
fn classify_event_identity_fields(event: &Event, receipt: &EventIdentityReceipt) {
    let Event {
        session: _,
        scope: _,
        seq: _,
        severity: _,
        kind: _,
        wall_ns: _,
    } = event;
    let EventIdentityReceipt {
        declared_identity_version: _,
        canonical_bytes: _,
        root: _,
    } = receipt;
}

fn esc(out: &mut String, s: &str) {
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            '\r' => out.push_str("\\r"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
}

fn push_str_field(out: &mut String, key: &str, val: &str) {
    let _ = write!(out, "\"{key}\":\"");
    esc(out, val);
    out.push('"');
}

fn push_hex_field(out: &mut String, key: &str, bytes: &[u8]) {
    let _ = write!(out, "\"{key}\":\"");
    for byte in bytes {
        let _ = write!(out, "{byte:02x}");
    }
    out.push('"');
}

fn push_process_artifact_fields(
    out: &mut String,
    artifact: Option<&process::DurableArtifactPointer>,
) {
    match artifact {
        Some(artifact) => {
            push_str_field(out, "artifact", artifact.artifact());
            out.push(',');
            push_str_field(out, "artifact_content_hash", artifact.content_hash());
        }
        None => out.push_str("\"artifact\":null,\"artifact_content_hash\":null"),
    }
}

/// Serialize a float for the wire: finite → shortest-round-trip; non-finite
/// → tagged string (JSON has no NaN/Inf; readers must handle both shapes).
fn push_f64(out: &mut String, key: &str, v: f64) {
    if v.is_finite() {
        let _ = write!(out, "\"{key}\":{v}");
    } else {
        let _ = write!(out, "\"{key}\":\"non-finite:{v}\"");
    }
}

/// Serialize one scoped receipt outcome as the `"outcome":{...}` member
/// shared by every outcome-bearing kind (i94v.7.3.1).
fn push_outcome_json(out: &mut String, outcome: &ScopedReceiptOutcome) {
    out.push_str(",\"outcome\":{");
    push_str_field(out, "scope", outcome.scope.name());
    out.push(',');
    push_str_field(out, "receipt", &outcome.receipt);
    out.push(',');
    push_str_field(out, "disposition", outcome.disposition.name());
    out.push(',');
    push_str_field(out, "predicate", outcome.predicate.name());
    out.push(',');
    push_str_field(out, "evidence_methods", &outcome.evidence_methods);
    out.push(',');
    push_str_field(out, "grade", outcome.grade.name());
    out.push(',');
    push_str_field(out, "applicability", outcome.applicability.name());
    out.push(',');
    push_str_field(out, "support", outcome.support.name());
    out.push(',');
    push_str_field(out, "completeness", outcome.completeness.name());
    out.push(',');
    push_str_field(out, "integrity", outcome.integrity.name());
    out.push(',');
    push_str_field(out, "promotion", outcome.promotion.name());
    out.push(',');
    push_str_field(out, "detail", &outcome.detail);
    out.push('}');
}

/// Bind one scoped receipt outcome into the typed identity frame, the exact
/// field-for-field mirror of [`push_outcome_json`].
fn bind_outcome_identity(
    builder: ident::IdentityBuilder,
    outcome: &ScopedReceiptOutcome,
) -> ident::IdentityBuilder {
    builder
        .str("outcome_scope", outcome.scope.name())
        .str("outcome_receipt", &outcome.receipt)
        .str("outcome_disposition", outcome.disposition.name())
        .str("outcome_predicate", outcome.predicate.name())
        .str("outcome_evidence_methods", &outcome.evidence_methods)
        .str("outcome_grade", outcome.grade.name())
        .str("outcome_applicability", outcome.applicability.name())
        .str("outcome_support", outcome.support.name())
        .str("outcome_completeness", outcome.completeness.name())
        .str("outcome_integrity", outcome.integrity.name())
        .str("outcome_promotion", outcome.promotion.name())
        .str("outcome_detail", &outcome.detail)
}

fn bind_u128_identity(
    builder: ident::IdentityBuilder,
    key: &str,
    value: u128,
) -> ident::IdentityBuilder {
    builder.bytes(key, &value.to_le_bytes())
}

fn bind_process_artifact_identity(
    builder: ident::IdentityBuilder,
    artifact: Option<&process::DurableArtifactPointer>,
) -> ident::IdentityBuilder {
    match artifact {
        Some(artifact) => builder
            .flag("artifact_present", true)
            .str("artifact", artifact.artifact())
            .str("artifact_content_hash", artifact.content_hash()),
        None => builder
            .flag("artifact_present", false)
            .str("artifact", "")
            .str("artifact_content_hash", ""),
    }
}

impl Event {
    /// Serialize the CONTENT portion (everything except `wall_ns`) in
    /// canonical field order.
    #[allow(clippy::too_many_lines)] // One exhaustive table preserves canonical JSON field order.
    fn content_json(&self) -> String {
        let mut s = String::with_capacity(160);
        let _ = write!(s, "{{\"v\":{SCHEMA_VERSION},");
        push_str_field(&mut s, "session", &self.session);
        s.push(',');
        push_str_field(&mut s, "scope", &self.scope);
        let _ = write!(s, ",\"seq\":{},", self.seq);
        push_str_field(&mut s, "severity", self.severity.name());
        s.push(',');
        push_str_field(&mut s, "kind", self.kind.kind_name());
        s.push_str(",\"payload\":{");
        match &self.kind {
            EventKind::SolverResidual {
                solver,
                iter,
                residual,
            } => {
                push_str_field(&mut s, "solver", solver);
                let _ = write!(s, ",\"iter\":{iter},");
                push_f64(&mut s, "residual", *residual);
            }
            EventKind::TileComplete { tile, kernel } => {
                let _ = write!(s, "\"tile\":{tile},");
                push_str_field(&mut s, "kernel", kernel);
            }
            EventKind::Cancellation { reason } => {
                push_str_field(&mut s, "reason", reason);
            }
            EventKind::BudgetDelta {
                resource,
                spent,
                remaining,
            } => {
                push_str_field(&mut s, "resource", resource);
                s.push(',');
                push_f64(&mut s, "spent", *spent);
                s.push(',');
                push_f64(&mut s, "remaining", *remaining);
            }
            EventKind::GradientCheck {
                op,
                max_rel_err,
                pass,
            } => {
                push_str_field(&mut s, "op", op);
                s.push(',');
                push_f64(&mut s, "max_rel_err", *max_rel_err);
                let _ = write!(s, ",\"pass\":{pass}");
            }
            EventKind::ConformanceCase {
                suite,
                case,
                pass,
                detail,
                seed,
            } => {
                push_str_field(&mut s, "suite", suite);
                s.push(',');
                push_str_field(&mut s, "case", case);
                let _ = write!(s, ",\"pass\":{pass},");
                push_str_field(&mut s, "detail", detail);
                let _ = write!(s, ",\"seed\":{seed}");
            }
            EventKind::BenchmarkResult {
                kernel,
                metric,
                value,
                machine,
            } => {
                push_str_field(&mut s, "kernel", kernel);
                s.push(',');
                push_str_field(&mut s, "metric", metric);
                s.push(',');
                push_f64(&mut s, "value", *value);
                let _ = write!(s, ",\"machine\":{machine}");
            }
            EventKind::StormAssertion { name, pass, seed } => {
                push_str_field(&mut s, "name", name);
                let _ = write!(s, ",\"pass\":{pass},\"seed\":{seed}");
            }
            EventKind::ManifestSelection {
                manifest,
                source_snapshot,
            } => {
                push_str_field(&mut s, "manifest", manifest);
                s.push(',');
                push_str_field(&mut s, "source_snapshot", source_snapshot);
            }
            EventKind::StratumExpansion {
                stratum,
                profile,
                cases,
            } => {
                push_str_field(&mut s, "stratum", stratum);
                s.push(',');
                push_str_field(&mut s, "profile", profile);
                let _ = write!(s, ",\"cases\":{cases}");
            }
            EventKind::DsrRun { run, outcome } => {
                push_str_field(&mut s, "run", run);
                push_outcome_json(&mut s, outcome);
            }
            EventKind::CampaignRun { campaign, outcome } => {
                push_str_field(&mut s, "campaign", campaign);
                push_outcome_json(&mut s, outcome);
            }
            EventKind::SubmissionDecision {
                job,
                lane,
                outcome,
                detail,
            } => {
                push_str_field(&mut s, "job", job);
                s.push(',');
                push_str_field(&mut s, "lane", lane.name());
                s.push(',');
                push_str_field(&mut s, "outcome", outcome.name());
                s.push(',');
                push_str_field(&mut s, "detail", detail);
            }
            EventKind::WorkIdentityBinding {
                job,
                attempt,
                operation,
            } => {
                push_str_field(&mut s, "job", job);
                let _ = write!(s, ",\"attempt\":{attempt},");
                push_str_field(&mut s, "operation", operation);
            }
            EventKind::LeaseCursor {
                queue,
                lease,
                holder,
                cursor,
            } => {
                push_str_field(&mut s, "queue", queue);
                s.push(',');
                push_str_field(&mut s, "lease", lease);
                s.push(',');
                push_str_field(&mut s, "holder", holder);
                let _ = write!(s, ",\"cursor\":{cursor}");
            }
            EventKind::AttachDetach {
                observer,
                target,
                action,
            } => {
                push_str_field(&mut s, "observer", observer);
                s.push(',');
                push_str_field(&mut s, "target", target);
                s.push(',');
                push_str_field(&mut s, "action", action.name());
            }
            EventKind::JourneyPhase {
                journey,
                phase,
                ordinal,
            } => {
                push_str_field(&mut s, "journey", journey);
                s.push(',');
                push_str_field(&mut s, "phase", phase);
                let _ = write!(s, ",\"ordinal\":{ordinal}");
            }
            EventKind::ScopeTileProgress { completed, total } => {
                let _ = write!(s, "\"completed\":{completed},\"total\":{total}");
            }
            EventKind::Heartbeat { worker, beat } => {
                push_str_field(&mut s, "worker", worker);
                let _ = write!(s, ",\"beat\":{beat}");
            }
            EventKind::ObservationGap {
                expected_seq,
                resumed_seq,
                reason,
            } => {
                let _ = write!(
                    s,
                    "\"expected_seq\":{expected_seq},\"resumed_seq\":{resumed_seq},"
                );
                push_str_field(&mut s, "reason", reason);
            }
            EventKind::ProcessFrame {
                process,
                stream,
                ordinal,
                class,
                payload,
                original_bytes,
                artifact,
            } => {
                push_str_field(&mut s, "process", process.as_str());
                s.push(',');
                push_str_field(&mut s, "stream", stream.name());
                let _ = write!(s, ",\"ordinal\":{ordinal},");
                push_str_field(&mut s, "class", class.name());
                s.push(',');
                push_hex_field(&mut s, "payload_hex", payload);
                let _ = write!(s, ",\"original_bytes\":{original_bytes},");
                push_process_artifact_fields(&mut s, artifact.as_ref());
            }
            EventKind::ProcessGap {
                process,
                stream,
                class,
                first_ordinal,
                last_ordinal,
                original_count,
                emitted_count,
                dropped_count,
                omitted_inline_bytes,
                reason,
                policy_version,
                artifact,
            } => {
                push_str_field(&mut s, "process", process.as_str());
                s.push(',');
                push_str_field(&mut s, "stream", stream.name());
                s.push(',');
                push_str_field(&mut s, "class", class.name());
                let _ = write!(
                    s,
                    ",\"first_ordinal\":{first_ordinal},\"last_ordinal\":{last_ordinal},\
                     \"original_count\":{original_count},\"emitted_count\":{emitted_count},\
                     \"dropped_count\":{dropped_count},\
                     \"omitted_inline_bytes\":{omitted_inline_bytes},"
                );
                push_str_field(&mut s, "reason", reason.name());
                let _ = write!(s, ",\"policy_version\":{policy_version},");
                push_process_artifact_fields(&mut s, artifact.as_ref());
            }
            EventKind::OracleComparison {
                oracle,
                subject,
                metric,
                observed,
                tolerance,
                pass,
                detail,
            } => {
                push_str_field(&mut s, "oracle", oracle);
                s.push(',');
                push_str_field(&mut s, "subject", subject);
                s.push(',');
                push_str_field(&mut s, "metric", metric);
                s.push(',');
                push_f64(&mut s, "observed", *observed);
                s.push(',');
                push_f64(&mut s, "tolerance", *tolerance);
                let _ = write!(s, ",\"pass\":{pass},");
                push_str_field(&mut s, "detail", detail);
            }
            EventKind::ToleranceDerivation {
                quantity,
                tolerance,
                basis,
            } => {
                push_str_field(&mut s, "quantity", quantity);
                s.push(',');
                push_f64(&mut s, "tolerance", *tolerance);
                s.push(',');
                push_str_field(&mut s, "basis", basis);
            }
            EventKind::ClaimAdjudication { claim, outcome } => {
                push_str_field(&mut s, "claim", claim);
                push_outcome_json(&mut s, outcome);
            }
            EventKind::CapabilityDomainDecision {
                capability,
                domain,
                decision,
                detail,
            } => {
                push_str_field(&mut s, "capability", capability);
                s.push(',');
                push_str_field(&mut s, "domain", domain);
                s.push(',');
                push_str_field(&mut s, "decision", decision.name());
                s.push(',');
                push_str_field(&mut s, "detail", detail);
            }
            EventKind::LifecycleTransition {
                entity,
                transition,
                detail,
            } => {
                push_str_field(&mut s, "entity", entity);
                s.push(',');
                push_str_field(&mut s, "transition", transition.name());
                s.push(',');
                push_str_field(&mut s, "detail", detail);
            }
            EventKind::ArtifactLifecycle {
                artifact,
                action,
                actor,
                detail,
            } => {
                push_str_field(&mut s, "artifact", artifact);
                s.push(',');
                push_str_field(&mut s, "action", action.name());
                s.push(',');
                push_str_field(&mut s, "actor", actor);
                s.push(',');
                push_str_field(&mut s, "detail", detail);
            }
            EventKind::VisualizationTransform {
                view,
                transform,
                source,
            } => {
                push_str_field(&mut s, "view", view);
                s.push(',');
                push_str_field(&mut s, "transform", transform);
                s.push(',');
                push_str_field(&mut s, "source", source);
            }
            EventKind::DiagnosticRepair {
                subject,
                action,
                pass,
                detail,
            } => {
                push_str_field(&mut s, "subject", subject);
                s.push(',');
                push_str_field(&mut s, "action", action);
                let _ = write!(s, ",\"pass\":{pass},");
                push_str_field(&mut s, "detail", detail);
            }
            EventKind::ContainmentNode {
                attempt,
                role,
                node,
                parent_role,
                parent,
                seq,
                dsr_run,
                campaign_run,
                shard,
                journey,
                case,
            } => {
                push_str_field(&mut s, "attempt", attempt);
                s.push(',');
                push_str_field(&mut s, "role", role);
                s.push(',');
                push_str_field(&mut s, "node", node);
                s.push(',');
                push_str_field(&mut s, "parent_role", parent_role);
                s.push(',');
                push_str_field(&mut s, "parent", parent);
                let _ = write!(s, ",\"seq\":{seq},");
                push_str_field(&mut s, "dsr_run", dsr_run);
                s.push(',');
                push_str_field(&mut s, "campaign_run", campaign_run);
                s.push(',');
                push_str_field(&mut s, "shard", shard);
                s.push(',');
                push_str_field(&mut s, "journey", journey);
                s.push(',');
                push_str_field(&mut s, "case", case);
            }
            EventKind::ContainmentGap {
                attempt,
                node_role,
                node,
                missing_parent_role,
                missing_parent,
            } => {
                push_str_field(&mut s, "attempt", attempt);
                s.push(',');
                push_str_field(&mut s, "node_role", node_role);
                s.push(',');
                push_str_field(&mut s, "node", node);
                s.push(',');
                push_str_field(&mut s, "missing_parent_role", missing_parent_role);
                s.push(',');
                push_str_field(&mut s, "missing_parent", missing_parent);
            }
            EventKind::RaceRecord {
                resource,
                schedule,
                pass,
                seed,
            } => {
                push_str_field(&mut s, "resource", resource);
                s.push(',');
                push_str_field(&mut s, "schedule", schedule);
                let _ = write!(s, ",\"pass\":{pass},\"seed\":{seed}");
            }
            EventKind::DegradationEvent {
                resource,
                limit,
                observed,
                action,
            } => {
                push_str_field(&mut s, "resource", resource);
                let _ = write!(s, ",\"limit\":{limit},\"observed\":{observed},");
                push_str_field(&mut s, "action", action);
            }
            EventKind::ImportReceipt {
                format,
                artifact,
                accepted,
                detail,
            } => {
                push_str_field(&mut s, "format", format);
                s.push(',');
                push_str_field(&mut s, "artifact", artifact);
                let _ = write!(s, ",\"accepted\":{accepted},");
                push_str_field(&mut s, "detail", detail);
            }
            EventKind::CertificateVerdict {
                certificate,
                pass,
                bound,
                detail,
            } => {
                push_str_field(&mut s, "certificate", certificate);
                let _ = write!(s, ",\"pass\":{pass},");
                push_f64(&mut s, "bound", *bound);
                s.push(',');
                push_str_field(&mut s, "detail", detail);
            }
            EventKind::Custom { name, json } => {
                push_str_field(&mut s, "name", name);
                let _ = write!(s, ",\"data\":{json}");
            }
        }
        s.push('}');
        s.push('}');
        s
    }

    /// Serialize one JSON-line: content with `wall_ns` spliced in as the
    /// LAST field before the closing brace (envelope position).
    #[must_use]
    pub fn to_jsonl(&self) -> String {
        let mut content = self.content_json();
        if let Some(ns) = self.wall_ns {
            content.pop(); // strip trailing '}'
            let _ = write!(content, ",\"wall_ns\":{ns}}}");
        }
        content
    }

    fn content_identity_with_versions(
        &self,
        event_identity_version: u32,
        event_wire_schema_version: u32,
    ) -> ident::ReplayIdentity {
        self.content_identity_with_schema(
            EVENT_CONTENT_IDENTITY_DOMAIN,
            event_identity_version,
            event_wire_schema_version,
        )
    }

    #[allow(clippy::too_many_lines)] // One exhaustive mirror keeps identity field order auditable against JSON.
    fn content_identity_with_schema(
        &self,
        artifact_domain: &str,
        event_identity_version: u32,
        event_wire_schema_version: u32,
    ) -> ident::ReplayIdentity {
        let Event {
            session,
            scope,
            seq,
            severity,
            kind,
            wall_ns: _,
        } = self;
        let builder = ident::IdentityBuilder::new(artifact_domain)
            .u64(
                "event_content_identity_version",
                u64::from(event_identity_version),
            )
            .u64(
                "event_wire_schema_version",
                u64::from(event_wire_schema_version),
            )
            .str("session", session)
            .str("scope", scope)
            .u64("seq", *seq)
            .str("severity", severity.name())
            .str("kind", kind.kind_name());
        let builder = match kind {
            EventKind::SolverResidual {
                solver,
                iter,
                residual,
            } => builder
                .str("solver", solver)
                .u64("iter", *iter)
                .f64_bits("residual", *residual),
            EventKind::TileComplete { tile, kernel } => {
                builder.u64("tile", *tile).str("kernel", kernel)
            }
            EventKind::Cancellation { reason } => builder.str("reason", reason),
            EventKind::BudgetDelta {
                resource,
                spent,
                remaining,
            } => builder
                .str("resource", resource)
                .f64_bits("spent", *spent)
                .f64_bits("remaining", *remaining),
            EventKind::GradientCheck {
                op,
                max_rel_err,
                pass,
            } => builder
                .str("op", op)
                .f64_bits("max_rel_err", *max_rel_err)
                .flag("pass", *pass),
            EventKind::ConformanceCase {
                suite,
                case,
                pass,
                detail,
                seed,
            } => builder
                .str("suite", suite)
                .str("case", case)
                .flag("pass", *pass)
                .str("detail", detail)
                .u64("seed", *seed),
            EventKind::BenchmarkResult {
                kernel,
                metric,
                value,
                machine,
            } => builder
                .str("kernel", kernel)
                .str("metric", metric)
                .f64_bits("value", *value)
                .u64("machine", *machine),
            EventKind::StormAssertion { name, pass, seed } => builder
                .str("name", name)
                .flag("pass", *pass)
                .u64("seed", *seed),
            EventKind::ManifestSelection {
                manifest,
                source_snapshot,
            } => builder
                .str("manifest", manifest)
                .str("source_snapshot", source_snapshot),
            EventKind::StratumExpansion {
                stratum,
                profile,
                cases,
            } => builder
                .str("stratum", stratum)
                .str("profile", profile)
                .u64("cases", *cases),
            EventKind::DsrRun { run, outcome } => {
                bind_outcome_identity(builder.str("run", run), outcome)
            }
            EventKind::CampaignRun { campaign, outcome } => {
                bind_outcome_identity(builder.str("campaign", campaign), outcome)
            }
            EventKind::SubmissionDecision {
                job,
                lane,
                outcome,
                detail,
            } => builder
                .str("job", job)
                .str("lane", lane.name())
                .str("outcome", outcome.name())
                .str("detail", detail),
            EventKind::WorkIdentityBinding {
                job,
                attempt,
                operation,
            } => builder
                .str("job", job)
                .u64("attempt", *attempt)
                .str("operation", operation),
            EventKind::LeaseCursor {
                queue,
                lease,
                holder,
                cursor,
            } => builder
                .str("queue", queue)
                .str("lease", lease)
                .str("holder", holder)
                .u64("cursor", *cursor),
            EventKind::AttachDetach {
                observer,
                target,
                action,
            } => builder
                .str("observer", observer)
                .str("target", target)
                .str("action", action.name()),
            EventKind::JourneyPhase {
                journey,
                phase,
                ordinal,
            } => builder
                .str("journey", journey)
                .str("phase", phase)
                .u64("ordinal", *ordinal),
            EventKind::ScopeTileProgress { completed, total } => {
                builder.u64("completed", *completed).u64("total", *total)
            }
            EventKind::Heartbeat { worker, beat } => {
                builder.str("worker", worker).u64("beat", *beat)
            }
            EventKind::ObservationGap {
                expected_seq,
                resumed_seq,
                reason,
            } => builder
                .u64("expected_seq", *expected_seq)
                .u64("resumed_seq", *resumed_seq)
                .str("reason", reason),
            EventKind::ProcessFrame {
                process,
                stream,
                ordinal,
                class,
                payload,
                original_bytes,
                artifact,
            } => bind_process_artifact_identity(
                builder
                    .str("process", process.as_str())
                    .str("stream", stream.name())
                    .u64("ordinal", *ordinal)
                    .str("class", class.name())
                    .bytes("payload", payload)
                    .u64("original_bytes", *original_bytes),
                artifact.as_ref(),
            ),
            EventKind::ProcessGap {
                process,
                stream,
                class,
                first_ordinal,
                last_ordinal,
                original_count,
                emitted_count,
                dropped_count,
                omitted_inline_bytes,
                reason,
                policy_version,
                artifact,
            } => {
                let builder = builder
                    .str("process", process.as_str())
                    .str("stream", stream.name())
                    .str("class", class.name())
                    .u64("first_ordinal", *first_ordinal)
                    .u64("last_ordinal", *last_ordinal);
                let builder = bind_u128_identity(builder, "original_count", *original_count);
                let builder = bind_u128_identity(builder, "emitted_count", *emitted_count);
                let builder = bind_u128_identity(builder, "dropped_count", *dropped_count);
                let builder =
                    bind_u128_identity(builder, "omitted_inline_bytes", *omitted_inline_bytes)
                        .str("reason", reason.name())
                        .u64("policy_version", u64::from(*policy_version));
                bind_process_artifact_identity(builder, artifact.as_ref())
            }
            EventKind::OracleComparison {
                oracle,
                subject,
                metric,
                observed,
                tolerance,
                pass,
                detail,
            } => builder
                .str("oracle", oracle)
                .str("subject", subject)
                .str("metric", metric)
                .f64_bits("observed", *observed)
                .f64_bits("tolerance", *tolerance)
                .flag("pass", *pass)
                .str("detail", detail),
            EventKind::ToleranceDerivation {
                quantity,
                tolerance,
                basis,
            } => builder
                .str("quantity", quantity)
                .f64_bits("tolerance", *tolerance)
                .str("basis", basis),
            EventKind::ClaimAdjudication { claim, outcome } => {
                bind_outcome_identity(builder.str("claim", claim), outcome)
            }
            EventKind::CapabilityDomainDecision {
                capability,
                domain,
                decision,
                detail,
            } => builder
                .str("capability", capability)
                .str("domain", domain)
                .str("decision", decision.name())
                .str("detail", detail),
            EventKind::LifecycleTransition {
                entity,
                transition,
                detail,
            } => builder
                .str("entity", entity)
                .str("transition", transition.name())
                .str("detail", detail),
            EventKind::ArtifactLifecycle {
                artifact,
                action,
                actor,
                detail,
            } => builder
                .str("artifact", artifact)
                .str("action", action.name())
                .str("actor", actor)
                .str("detail", detail),
            EventKind::VisualizationTransform {
                view,
                transform,
                source,
            } => builder
                .str("view", view)
                .str("transform", transform)
                .str("source", source),
            EventKind::DiagnosticRepair {
                subject,
                action,
                pass,
                detail,
            } => builder
                .str("subject", subject)
                .str("action", action)
                .flag("pass", *pass)
                .str("detail", detail),
            EventKind::ContainmentNode {
                attempt,
                role,
                node,
                parent_role,
                parent,
                seq,
                dsr_run,
                campaign_run,
                shard,
                journey,
                case,
            } => builder
                .str("attempt", attempt)
                .str("role", role)
                .str("node", node)
                .str("parent_role", parent_role)
                .str("parent", parent)
                .u64("seq", *seq)
                .str("dsr_run", dsr_run)
                .str("campaign_run", campaign_run)
                .str("shard", shard)
                .str("journey", journey)
                .str("case", case),
            EventKind::ContainmentGap {
                attempt,
                node_role,
                node,
                missing_parent_role,
                missing_parent,
            } => builder
                .str("attempt", attempt)
                .str("node_role", node_role)
                .str("node", node)
                .str("missing_parent_role", missing_parent_role)
                .str("missing_parent", missing_parent),
            EventKind::RaceRecord {
                resource,
                schedule,
                pass,
                seed,
            } => builder
                .str("resource", resource)
                .str("schedule", schedule)
                .flag("pass", *pass)
                .u64("seed", *seed),
            EventKind::DegradationEvent {
                resource,
                limit,
                observed,
                action,
            } => builder
                .str("resource", resource)
                .u64("limit", *limit)
                .u64("observed", *observed)
                .str("action", action),
            EventKind::ImportReceipt {
                format,
                artifact,
                accepted,
                detail,
            } => builder
                .str("format", format)
                .str("artifact", artifact)
                .flag("accepted", *accepted)
                .str("detail", detail),
            EventKind::CertificateVerdict {
                certificate,
                pass,
                bound,
                detail,
            } => builder
                .str("certificate", certificate)
                .flag("pass", *pass)
                .f64_bits("bound", *bound)
                .str("detail", detail),
            EventKind::Custom { name, json } => builder
                .str("name", name)
                .bytes("custom_json_opaque_utf8", json.as_bytes()),
        };
        builder
            .exclude(
                "wall_ns",
                "wall-clock is observability envelope, not replay identity",
            )
            .finish()
    }

    /// Canonical typed identity for the deterministic event content.
    ///
    /// This is deliberately independent of the display-oriented JSON line:
    /// floats bind by their exact bit patterns, every payload variant has a
    /// closed typed encoder, and wall-clock remains an explicit exclusion.
    #[must_use]
    pub fn content_identity(&self) -> ident::ReplayIdentity {
        self.content_identity_with_versions(EVENT_CONTENT_IDENTITY_VERSION, SCHEMA_VERSION)
    }

    /// Capture a retained receipt containing the declared event-identity
    /// version, exact canonical bytes, and their root.
    #[must_use]
    pub fn content_identity_receipt(&self) -> EventIdentityReceipt {
        let identity = self.content_identity();
        EventIdentityReceipt::from_retained_parts(
            EVENT_CONTENT_IDENTITY_VERSION,
            identity.canonical_bytes().to_vec(),
            identity.root(),
        )
    }

    /// Admit a retained identity receipt only when all three proof surfaces
    /// agree: declared semantics, root-over-retained-bytes, and the exact
    /// canonical identity of this event.
    ///
    /// # Errors
    /// [`EventIdentityAdmissionError`] if any proof surface differs.
    pub fn admit_content_identity(
        &self,
        receipt: &EventIdentityReceipt,
    ) -> Result<(), EventIdentityAdmissionError> {
        check_event_content_identity_version(receipt.declared_identity_version)
            .map_err(EventIdentityAdmissionError::UnsupportedVersion)?;

        let computed = fnv1a64(&receipt.canonical_bytes);
        if computed != receipt.root {
            return Err(EventIdentityAdmissionError::RootMismatch {
                declared: receipt.root,
                computed,
            });
        }

        let expected = self.content_identity();
        if receipt.canonical_bytes.as_slice() != expected.canonical_bytes()
            || receipt.root != expected.root()
        {
            return Err(EventIdentityAdmissionError::CanonicalBytesMismatch {
                declared_root: receipt.root,
                expected_root: expected.root(),
            });
        }
        Ok(())
    }

    /// Deterministic FNV-1a root over [`Event::content_identity`]'s exact
    /// typed bytes. Not cryptographic; ledger-grade content addressing uses
    /// the same canonical identity bytes under a stronger digest.
    #[must_use]
    pub fn content_hash(&self) -> u64 {
        self.content_identity().root()
    }
}

/// FNV-1a 64-bit (in-house, deterministic across platforms).
#[must_use]
pub fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// A per-scope emitter handing out monotone sequence numbers.
#[derive(Debug)]
pub struct Emitter {
    session: String,
    scope: String,
    /// Byte length of the immutable base scope; `exit_scope` never
    /// truncates past it.
    base_scope_len: usize,
    seq: u64,
}

impl Emitter {
    /// Create an emitter for one (session, scope) pair.
    #[must_use]
    pub fn new(session: impl Into<String>, scope: impl Into<String>) -> Self {
        let scope = scope.into();
        let base_scope_len = scope.len();
        Emitter {
            session: session.into(),
            scope,
            base_scope_len,
            seq: 0,
        }
    }

    /// Enter a child scope segment: the asupersync scope tree IS the trace
    /// tree (huq.16), so the emitter mirrors it with explicit enter/exit
    /// rather than thread-local magic. Entering emits nothing by itself;
    /// `seq` stays one monotone stream per emitter so interleaved child
    /// scopes replay in exact emission order.
    ///
    /// # Errors
    /// Refuses empty segments and segments containing `/` (the path
    /// separator) or control characters — a forged segment must not be able
    /// to impersonate a different tree position.
    pub fn enter_scope(&mut self, segment: &str) -> Result<(), SchemaError> {
        if segment.is_empty() {
            return Err(SchemaError {
                at: 0,
                message: "scope segment must be non-empty".to_string(),
            });
        }
        if segment.contains('/') || segment.chars().any(char::is_control) {
            return Err(SchemaError {
                at: 0,
                message: format!(
                    "scope segment {segment:?} must not contain '/' or control characters \
                     (path forgery guard)"
                ),
            });
        }
        self.scope.push('/');
        self.scope.push_str(segment);
        Ok(())
    }

    /// Exit the innermost child scope.
    ///
    /// # Errors
    /// Refuses at the base scope: exits must pair with entries, and a
    /// mismatched exit is a harness bug worth surfacing, not masking.
    pub fn exit_scope(&mut self) -> Result<(), SchemaError> {
        let Some(cut) = self.scope[self.base_scope_len..].rfind('/') else {
            return Err(SchemaError {
                at: 0,
                message: "exit_scope at the base scope: unbalanced enter/exit".to_string(),
            });
        };
        self.scope.truncate(self.base_scope_len + cut);
        Ok(())
    }

    /// The live slash-separated scope path events are currently stamped with.
    #[must_use]
    pub fn current_scope(&self) -> &str {
        &self.scope
    }

    /// Depth of entered child segments below the base scope.
    #[must_use]
    pub fn scope_depth(&self) -> usize {
        self.scope[self.base_scope_len..]
            .chars()
            .filter(|&c| c == '/')
            .count()
    }

    /// Build the next event (seq auto-increments). `wall_ns` is supplied by
    /// the caller because THIS crate must stay deterministic and I/O-free;
    /// runtime layers pass real clocks, tests pass `None`.
    pub fn emit(&mut self, severity: Severity, kind: EventKind, wall_ns: Option<u64>) -> Event {
        let e = Event {
            session: self.session.clone(),
            scope: self.scope.clone(),
            seq: self.seq,
            severity,
            kind,
            wall_ns,
        };
        self.seq += 1;
        e
    }

    /// Build and emit an event only when `gate` admits its logical opportunity.
    ///
    /// Trace and info decisions use the caller-supplied, zero-based
    /// `logical_ordinal`, not the emitter sequence or invocation order. Warn
    /// and error events are always admitted. The builder supplies both the
    /// typed payload and optional caller clock so a rejected opportunity does
    /// not construct either. `Some` events consume the same sequence stream as
    /// [`Self::emit`]; `None` does not advance it.
    pub fn emit_gated(
        &mut self,
        gate: EmissionGate,
        severity: Severity,
        logical_ordinal: u64,
        build: impl FnOnce() -> (EventKind, Option<u64>),
    ) -> Option<Event> {
        if !gate.admits(severity, logical_ordinal) {
            return None;
        }

        let (kind, wall_ns) = build();
        Some(self.emit(severity, kind, wall_ns))
    }
}

/// Validation failure for a JSON-line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SchemaError {
    /// Byte position (approximate for structural errors).
    pub at: usize,
    /// What is wrong and how to fix it.
    pub message: String,
}

impl fmt::Display for SchemaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "event schema violation at byte {}: {}",
            self.at, self.message
        )
    }
}

impl core::error::Error for SchemaError {}

/// Registry of known kinds (validator uses it; keep in sync with EventKind).
pub const KNOWN_KINDS: &[&str] = &[
    "solver_residual",
    "tile_complete",
    "cancellation",
    "budget_delta",
    "gradient_check",
    "conformance_case",
    "benchmark_result",
    "storm_assertion",
    "manifest_selection",
    "stratum_expansion",
    "dsr_run",
    "campaign_run",
    "submission_decision",
    "work_identity_binding",
    "lease_cursor",
    "attach_detach",
    "journey_phase",
    "scope_tile_progress",
    "heartbeat",
    "observation_gap",
    "process_frame",
    "process_gap",
    "oracle_comparison",
    "tolerance_derivation",
    "claim_adjudication",
    "capability_domain_decision",
    "lifecycle_transition",
    "artifact_lifecycle",
    "visualization_transform",
    "diagnostic_repair",
    "containment_node",
    "containment_gap",
    "race_record",
    "degradation_event",
    "import_receipt",
    "certificate_verdict",
    "custom",
];

/// Strict structural validation of one JSON-line: required envelope keys in
/// canonical order, known kind, balanced payload object. (The writer is
/// ours; deviation means corruption, so the checks are cheap and strict
/// rather than a full JSON parser.)
///
/// # Errors
/// Returns [`SchemaError`] naming the first violated requirement.
pub fn validate_line(line: &str) -> Result<(), SchemaError> {
    let need = |cond: bool, at: usize, msg: &str| -> Result<(), SchemaError> {
        if cond {
            Ok(())
        } else {
            Err(SchemaError {
                at,
                message: msg.to_string(),
            })
        }
    };
    need(
        line.starts_with('{') && line.ends_with('}'),
        0,
        "line must be one JSON object",
    )?;
    let ver = format!("{{\"v\":{SCHEMA_VERSION},");
    need(
        line.starts_with(&ver),
        0,
        "first field must be the schema version \"v\"",
    )?;
    for key in [
        "\"session\":",
        "\"scope\":",
        "\"seq\":",
        "\"severity\":",
        "\"kind\":",
        "\"payload\":{",
    ] {
        need(
            line.contains(key),
            0,
            &format!("missing required field {key}"),
        )?;
    }
    let kind_pos = line.find("\"kind\":\"").ok_or(SchemaError {
        at: 0,
        message: "kind must be a string".to_string(),
    })?;
    let kind_rest = &line[kind_pos + 8..];
    let kind_end = kind_rest.find('"').unwrap_or(0);
    let kind = &kind_rest[..kind_end];
    need(
        KNOWN_KINDS.contains(&kind),
        kind_pos,
        &format!("unknown kind {kind:?}; register new kinds in the payload registry"),
    )?;
    // Braces balance (escaped quotes handled by the writer's escaping).
    let mut depth = 0i64;
    let mut in_str = false;
    let mut prev_backslash = false;
    for (i, c) in line.char_indices() {
        match c {
            '"' if !prev_backslash => in_str = !in_str,
            '{' if !in_str => depth += 1,
            '}' if !in_str => depth -= 1,
            _ => {}
        }
        prev_backslash = c == '\\' && !prev_backslash;
        if depth < 0 {
            return Err(SchemaError {
                at: i,
                message: "unbalanced braces".to_string(),
            });
        }
    }
    need(
        depth == 0 && !in_str,
        line.len(),
        "unbalanced braces or unterminated string",
    )?;
    Ok(())
}

/// The failure-record completeness lint, v1 (the "log-replay lint" of the
/// observability bead): failure verdicts must carry enough to reproduce.
///
/// # Errors
/// Returns [`SchemaError`] describing the missing reproduction ingredient.
#[allow(clippy::too_many_lines)] // One fail-closed policy table keeps every refusal and precedence visible.
pub fn lint_failure_record(event: &Event) -> Result<(), SchemaError> {
    match &event.kind {
        EventKind::ConformanceCase {
            pass: false,
            detail,
            ..
        } if detail.is_empty() => Err(SchemaError {
            at: 0,
            message: "failing conformance case must carry a non-empty detail \
                          (reproduce-from-log-alone doctrine)"
                .to_string(),
        }),
        EventKind::StormAssertion {
            pass: false, seed, ..
        } if *seed == 0 => Err(SchemaError {
            at: 0,
            message: "failing storm assertion must carry its replay seed".to_string(),
        }),
        EventKind::GradientCheck {
            pass: false, op, ..
        } if op.is_empty() => Err(SchemaError {
            at: 0,
            message: "failing gradient check must name its operator".to_string(),
        }),
        EventKind::DsrRun { outcome, .. }
        | EventKind::CampaignRun { outcome, .. }
        | EventKind::ClaimAdjudication { outcome, .. }
            if matches!(
                outcome.disposition,
                ExecutionDisposition::Failed
                    | ExecutionDisposition::Refused
                    | ExecutionDisposition::Indeterminate
            ) && outcome.detail.is_empty() =>
        {
            Err(SchemaError {
                at: 0,
                message: "failed/refused/indeterminate scoped outcome must carry a \
                          non-empty diagnostic detail (reproduce-from-log-alone doctrine)"
                    .to_string(),
            })
        }
        EventKind::SubmissionDecision {
            outcome: SubmissionOutcome::Refused | SubmissionOutcome::Deferred,
            detail,
            ..
        } if detail.is_empty() => Err(SchemaError {
            at: 0,
            message: "refused/deferred submission must carry a non-empty decision \
                          detail (reproduce-from-log-alone doctrine)"
                .to_string(),
        }),
        EventKind::OracleComparison {
            pass: false,
            detail,
            ..
        } if detail.is_empty() => Err(SchemaError {
            at: 0,
            message: "failing oracle comparison must carry a non-empty detail \
                          (reproduce-from-log-alone doctrine)"
                .to_string(),
        }),
        EventKind::ToleranceDerivation { basis, .. } if basis.is_empty() => Err(SchemaError {
            at: 0,
            message: "a tolerance derivation must name its basis: an unexplained \
                          tolerance is the silent-constant bug class"
                .to_string(),
        }),
        EventKind::CapabilityDomainDecision {
            decision: CapabilityDecision::Refused | CapabilityDecision::Restricted,
            detail,
            ..
        } if detail.is_empty() => Err(SchemaError {
            at: 0,
            message: "refused/restricted capability decision must carry a non-empty \
                          detail (reproduce-from-log-alone doctrine)"
                .to_string(),
        }),
        EventKind::LifecycleTransition {
            transition: LifecycleTransitionKind::Fault,
            detail,
            ..
        } if detail.is_empty() => Err(SchemaError {
            at: 0,
            message: "a contained fault must carry a non-empty detail \
                          (reproduce-from-log-alone doctrine)"
                .to_string(),
        }),
        EventKind::ArtifactLifecycle {
            action: ArtifactAction::Redact | ArtifactAction::Share,
            detail,
            ..
        } if detail.is_empty() => Err(SchemaError {
            at: 0,
            message: "redact/share artifact actions must carry a non-empty detail: \
                          redactions need a rationale, shares a target"
                .to_string(),
        }),
        EventKind::DiagnosticRepair {
            pass: false,
            detail,
            ..
        } if detail.is_empty() => Err(SchemaError {
            at: 0,
            message: "an ineffective repair must carry a non-empty detail \
                          (reproduce-from-log-alone doctrine)"
                .to_string(),
        }),
        EventKind::ObservationGap { reason, .. } if reason.is_empty() => Err(SchemaError {
            at: 0,
            message: "an observation gap must carry a non-empty reason: readers treat \
                          the gap as unknown execution, never as absence of events"
                .to_string(),
        }),
        EventKind::ProcessFrame { ordinal: 0, .. } => Err(SchemaError {
            at: 0,
            message: "a process frame ordinal must be one-based".to_string(),
        }),
        EventKind::ProcessFrame {
            payload,
            original_bytes,
            ..
        } if u64::try_from(payload.len()).is_err()
            || *original_bytes < u64::try_from(payload.len()).unwrap_or(u64::MAX) =>
        {
            Err(SchemaError {
                at: 0,
                message: "a process frame's original byte count cannot be smaller than its \
                          retained inline payload"
                    .to_string(),
            })
        }
        EventKind::ProcessGap {
            first_ordinal,
            last_ordinal,
            policy_version,
            ..
        } if *first_ordinal == 0 || last_ordinal < first_ordinal || *policy_version == 0 => {
            Err(SchemaError {
                at: 0,
                message:
                    "a process gap needs a one-based ordered range and non-zero policy version"
                        .to_string(),
            })
        }
        EventKind::ProcessGap {
            first_ordinal,
            last_ordinal,
            original_count,
            emitted_count,
            dropped_count,
            ..
        } if emitted_count.checked_add(*dropped_count) != Some(*original_count)
            || *original_count
                != u128::from(*last_ordinal)
                    .checked_sub(u128::from(*first_ordinal))
                    .and_then(|span| span.checked_add(1))
                    .unwrap_or(0) =>
        {
            Err(SchemaError {
                at: 0,
                message:
                    "a process gap must balance emitted+dropped counts over its exact ordinal range"
                        .to_string(),
            })
        }
        EventKind::ProcessGap {
            reason: process::LossReason::DurableSpill,
            artifact: None,
            ..
        } => Err(SchemaError {
            at: 0,
            message: "a durable-spill process gap must carry its committed artifact pointer"
                .to_string(),
        }),
        EventKind::RaceRecord {
            pass: false, seed, ..
        } if *seed == 0 => Err(SchemaError {
            at: 0,
            message: "failing race record must carry its schedule replay seed".to_string(),
        }),
        EventKind::ImportReceipt {
            accepted: false,
            detail,
            ..
        } if detail.is_empty() => Err(SchemaError {
            at: 0,
            message: "refused import must carry a non-empty quarantine detail \
                          (reproduce-from-log-alone doctrine)"
                .to_string(),
        }),
        EventKind::CertificateVerdict {
            pass: false,
            detail,
            ..
        } if detail.is_empty() => Err(SchemaError {
            at: 0,
            message: "failing certificate verdict must carry a non-empty detail \
                          (reproduce-from-log-alone doctrine)"
                .to_string(),
        }),
        _ => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_events() -> Vec<Event> {
        let mut em = Emitter::new("study-x", "op-1/kernel-cg");
        vec![
            em.emit(
                Severity::Info,
                EventKind::SolverResidual {
                    solver: "cg".into(),
                    iter: 3,
                    residual: 1.5e-7,
                },
                Some(1_000),
            ),
            em.emit(
                Severity::Trace,
                EventKind::TileComplete {
                    tile: 42,
                    kernel: "lbm_d3q19".into(),
                },
                None,
            ),
            em.emit(
                Severity::Warn,
                EventKind::Cancellation {
                    reason: "budget".into(),
                },
                Some(1_500),
            ),
            em.emit(
                Severity::Warn,
                EventKind::BudgetDelta {
                    resource: "wall_s".into(),
                    spent: 12.5,
                    remaining: 7187.5,
                },
                Some(2_000),
            ),
            em.emit(
                Severity::Info,
                EventKind::GradientCheck {
                    op: "poisson".into(),
                    max_rel_err: 2.5e-8,
                    pass: true,
                },
                None,
            ),
            em.emit(
                Severity::Error,
                EventKind::ConformanceCase {
                    suite: "fs-qty/conformance".into(),
                    case: "qty-001".into(),
                    pass: false,
                    detail: "value mismatch: got 0.13, want 0.12".into(),
                    seed: 7,
                },
                Some(3_000),
            ),
            em.emit(
                Severity::Info,
                EventKind::BenchmarkResult {
                    kernel: "gemm".into(),
                    metric: "gflops".into(),
                    value: 123.5,
                    machine: 0x1234,
                },
                None,
            ),
            em.emit(
                Severity::Info,
                EventKind::StormAssertion {
                    name: "no-arena-leak".into(),
                    pass: true,
                    seed: 99,
                },
                None,
            ),
            em.emit(
                Severity::Info,
                EventKind::Custom {
                    name: "regime-report".into(),
                    json: r#"{"re":100.5,"we":0.3}"#.into(),
                },
                None,
            ),
        ]
    }

    fn read_identity_len(bytes: &[u8], cursor: &mut usize) -> usize {
        let end = cursor
            .checked_add(core::mem::size_of::<u64>())
            .expect("test identity length offset fits usize");
        let encoded: [u8; 8] = bytes[*cursor..end]
            .try_into()
            .expect("identity length has eight bytes");
        *cursor = end;
        usize::try_from(u64::from_le_bytes(encoded)).expect("test identity length fits usize")
    }

    fn identity_field<'a>(canonical: &'a [u8], wanted: &str) -> (u8, &'a [u8]) {
        assert!(
            canonical.starts_with(ident::REPLAY_IDENTITY_DOMAIN.as_bytes()),
            "identity frame must start with its declared replay domain"
        );
        let mut cursor = ident::REPLAY_IDENTITY_DOMAIN.len() + core::mem::size_of::<u32>();
        let kind_len = read_identity_len(canonical, &mut cursor);
        cursor += kind_len;
        while cursor < canonical.len() {
            let tag = canonical[cursor];
            cursor += 1;
            let key_len = read_identity_len(canonical, &mut cursor);
            let key_end = cursor + key_len;
            let key = core::str::from_utf8(&canonical[cursor..key_end])
                .expect("identity field keys are UTF-8");
            cursor = key_end;
            let value_len = read_identity_len(canonical, &mut cursor);
            let value_end = cursor + value_len;
            if key == wanted {
                return (tag, &canonical[cursor..value_end]);
            }
            cursor = value_end;
        }
        panic!("identity field {wanted:?} was not encoded");
    }

    fn event_with_kind(kind: EventKind) -> Event {
        Event {
            session: "identity-session".into(),
            scope: "identity-scope".into(),
            seq: 17,
            severity: Severity::Info,
            kind,
            wall_ns: None,
        }
    }

    fn assert_payload_mutations(
        base: EventKind,
        mutations: Vec<(&'static str, EventKind)>,
        observed: &mut Vec<&'static str>,
    ) {
        let base_root = event_with_kind(base).content_hash();
        for (field, mutation) in mutations {
            assert_ne!(
                event_with_kind(mutation).content_hash(),
                base_root,
                "mutating {field} must move the exact event identity"
            );
            observed.push(field);
        }
    }

    #[test]
    fn every_kind_serializes_and_validates() {
        for e in sample_events() {
            let line = e.to_jsonl();
            validate_line(&line).unwrap_or_else(|err| panic!("{line}: {err}"));
        }
    }

    #[test]
    fn wall_clock_is_envelope_only() {
        let mut a = sample_events().remove(0);
        let h1 = a.content_hash();
        a.wall_ns = Some(999_999_999);
        let h2 = a.content_hash();
        assert_eq!(h1, h2, "content hash must exclude wall clock");
        // ...but the serialized line DOES carry it, as the last field.
        assert!(a.to_jsonl().ends_with(",\"wall_ns\":999999999}"));
    }

    #[test]
    fn content_hash_is_sensitive_to_content() {
        let events = sample_events();
        let mut hashes: Vec<u64> = events.iter().map(Event::content_hash).collect();
        hashes.sort_unstable();
        hashes.dedup();
        assert_eq!(
            hashes.len(),
            events.len(),
            "distinct events must hash distinctly"
        );
    }

    #[test]
    fn content_identity_preserves_bits_that_display_json_collapses() {
        let event = |residual| Event {
            session: "same-session".into(),
            scope: "same-scope".into(),
            seq: 7,
            severity: Severity::Info,
            kind: EventKind::SolverResidual {
                solver: "cg".into(),
                iter: 3,
                residual,
            },
            wall_ns: None,
        };
        let first = event(f64::from_bits(0x7ff8_0000_0000_0001));
        let second = event(f64::from_bits(0x7ff8_0000_0000_0002));
        assert_eq!(
            first.to_jsonl(),
            second.to_jsonl(),
            "tagged display JSON intentionally does not expose NaN payload bits"
        );
        assert_ne!(first.content_hash(), second.content_hash());
        assert_ne!(
            first.content_identity().canonical_bytes(),
            second.content_identity().canonical_bytes()
        );
        assert_eq!(
            first.content_hash(),
            fnv1a64(first.content_identity().canonical_bytes()),
            "the stored root is derived from the exact canonical bytes"
        );
    }

    #[test]
    // This exhaustive mutation table is deliberately kept in one test so a
    // newly added event field cannot escape the single identity-coverage audit.
    #[allow(clippy::too_many_lines)]
    fn every_event_kind_payload_field_moves_identity() {
        let mut observed = Vec::new();

        assert_payload_mutations(
            EventKind::SolverResidual {
                solver: "cg".into(),
                iter: 3,
                residual: 0.25,
            },
            vec![
                (
                    "solver_residual.solver",
                    EventKind::SolverResidual {
                        solver: "gmres".into(),
                        iter: 3,
                        residual: 0.25,
                    },
                ),
                (
                    "solver_residual.iter",
                    EventKind::SolverResidual {
                        solver: "cg".into(),
                        iter: 4,
                        residual: 0.25,
                    },
                ),
                (
                    "solver_residual.residual",
                    EventKind::SolverResidual {
                        solver: "cg".into(),
                        iter: 3,
                        residual: 0.5,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::TileComplete {
                tile: 7,
                kernel: "lbm".into(),
            },
            vec![
                (
                    "tile_complete.tile",
                    EventKind::TileComplete {
                        tile: 8,
                        kernel: "lbm".into(),
                    },
                ),
                (
                    "tile_complete.kernel",
                    EventKind::TileComplete {
                        tile: 7,
                        kernel: "gemm".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::Cancellation {
                reason: "budget".into(),
            },
            vec![(
                "cancellation.reason",
                EventKind::Cancellation {
                    reason: "panic".into(),
                },
            )],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::BudgetDelta {
                resource: "wall_s".into(),
                spent: 1.0,
                remaining: 9.0,
            },
            vec![
                (
                    "budget_delta.resource",
                    EventKind::BudgetDelta {
                        resource: "mem_bytes".into(),
                        spent: 1.0,
                        remaining: 9.0,
                    },
                ),
                (
                    "budget_delta.spent",
                    EventKind::BudgetDelta {
                        resource: "wall_s".into(),
                        spent: 2.0,
                        remaining: 9.0,
                    },
                ),
                (
                    "budget_delta.remaining",
                    EventKind::BudgetDelta {
                        resource: "wall_s".into(),
                        spent: 1.0,
                        remaining: 8.0,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::GradientCheck {
                op: "poisson".into(),
                max_rel_err: 0.125,
                pass: true,
            },
            vec![
                (
                    "gradient_check.op",
                    EventKind::GradientCheck {
                        op: "elasticity".into(),
                        max_rel_err: 0.125,
                        pass: true,
                    },
                ),
                (
                    "gradient_check.max_rel_err",
                    EventKind::GradientCheck {
                        op: "poisson".into(),
                        max_rel_err: 0.25,
                        pass: true,
                    },
                ),
                (
                    "gradient_check.pass",
                    EventKind::GradientCheck {
                        op: "poisson".into(),
                        max_rel_err: 0.125,
                        pass: false,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::ConformanceCase {
                suite: "suite-a".into(),
                case: "case-a".into(),
                pass: true,
                detail: "ok".into(),
                seed: 11,
            },
            vec![
                (
                    "conformance_case.suite",
                    EventKind::ConformanceCase {
                        suite: "suite-b".into(),
                        case: "case-a".into(),
                        pass: true,
                        detail: "ok".into(),
                        seed: 11,
                    },
                ),
                (
                    "conformance_case.case",
                    EventKind::ConformanceCase {
                        suite: "suite-a".into(),
                        case: "case-b".into(),
                        pass: true,
                        detail: "ok".into(),
                        seed: 11,
                    },
                ),
                (
                    "conformance_case.pass",
                    EventKind::ConformanceCase {
                        suite: "suite-a".into(),
                        case: "case-a".into(),
                        pass: false,
                        detail: "ok".into(),
                        seed: 11,
                    },
                ),
                (
                    "conformance_case.detail",
                    EventKind::ConformanceCase {
                        suite: "suite-a".into(),
                        case: "case-a".into(),
                        pass: true,
                        detail: "different".into(),
                        seed: 11,
                    },
                ),
                (
                    "conformance_case.seed",
                    EventKind::ConformanceCase {
                        suite: "suite-a".into(),
                        case: "case-a".into(),
                        pass: true,
                        detail: "ok".into(),
                        seed: 12,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::BenchmarkResult {
                kernel: "gemm".into(),
                metric: "gflops".into(),
                value: 100.0,
                machine: 7,
            },
            vec![
                (
                    "benchmark_result.kernel",
                    EventKind::BenchmarkResult {
                        kernel: "spmv".into(),
                        metric: "gflops".into(),
                        value: 100.0,
                        machine: 7,
                    },
                ),
                (
                    "benchmark_result.metric",
                    EventKind::BenchmarkResult {
                        kernel: "gemm".into(),
                        metric: "bandwidth_gbs".into(),
                        value: 100.0,
                        machine: 7,
                    },
                ),
                (
                    "benchmark_result.value",
                    EventKind::BenchmarkResult {
                        kernel: "gemm".into(),
                        metric: "gflops".into(),
                        value: 101.0,
                        machine: 7,
                    },
                ),
                (
                    "benchmark_result.machine",
                    EventKind::BenchmarkResult {
                        kernel: "gemm".into(),
                        metric: "gflops".into(),
                        value: 100.0,
                        machine: 8,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::StormAssertion {
                name: "no-leak".into(),
                pass: true,
                seed: 23,
            },
            vec![
                (
                    "storm_assertion.name",
                    EventKind::StormAssertion {
                        name: "cancel-latency".into(),
                        pass: true,
                        seed: 23,
                    },
                ),
                (
                    "storm_assertion.pass",
                    EventKind::StormAssertion {
                        name: "no-leak".into(),
                        pass: false,
                        seed: 23,
                    },
                ),
                (
                    "storm_assertion.seed",
                    EventKind::StormAssertion {
                        name: "no-leak".into(),
                        pass: true,
                        seed: 24,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::ManifestSelection {
                manifest: "ab12".into(),
                source_snapshot: "cd34".into(),
            },
            vec![
                (
                    "manifest_selection.manifest",
                    EventKind::ManifestSelection {
                        manifest: "ef56".into(),
                        source_snapshot: "cd34".into(),
                    },
                ),
                (
                    "manifest_selection.source_snapshot",
                    EventKind::ManifestSelection {
                        manifest: "ab12".into(),
                        source_snapshot: "0789".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::StratumExpansion {
                stratum: "core".into(),
                profile: "nightly".into(),
                cases: 40,
            },
            vec![
                (
                    "stratum_expansion.stratum",
                    EventKind::StratumExpansion {
                        stratum: "max".into(),
                        profile: "nightly".into(),
                        cases: 40,
                    },
                ),
                (
                    "stratum_expansion.profile",
                    EventKind::StratumExpansion {
                        stratum: "core".into(),
                        profile: "release".into(),
                        cases: 40,
                    },
                ),
                (
                    "stratum_expansion.cases",
                    EventKind::StratumExpansion {
                        stratum: "core".into(),
                        profile: "nightly".into(),
                        cases: 41,
                    },
                ),
            ],
            &mut observed,
        );
        {
            let base_outcome = || ScopedReceiptOutcome {
                scope: ReceiptScope::Job,
                receipt: "ab12".into(),
                disposition: ExecutionDisposition::Completed,
                predicate: PredicateOutcome::Satisfied,
                evidence_methods: "oracle,enclosure".into(),
                grade: EpistemicGrade::Verified,
                applicability: DomainApplicability::InDomain,
                support: OperationalSupport::Supported,
                completeness: EvidenceCompleteness::Complete,
                integrity: EvidenceIntegrity::Intact,
                promotion: PromotionEffect::Unchanged,
                detail: "ok".into(),
            };
            let with = |mutate: &dyn Fn(&mut ScopedReceiptOutcome)| {
                let mut o = base_outcome();
                mutate(&mut o);
                EventKind::DsrRun {
                    run: "run-1".into(),
                    outcome: o,
                }
            };
            assert_payload_mutations(
                EventKind::DsrRun {
                    run: "run-1".into(),
                    outcome: base_outcome(),
                },
                vec![
                    (
                        "dsr_run.run",
                        EventKind::DsrRun {
                            run: "run-2".into(),
                            outcome: base_outcome(),
                        },
                    ),
                    (
                        "dsr_run.outcome.scope",
                        with(&|o| o.scope = ReceiptScope::Attempt),
                    ),
                    (
                        "dsr_run.outcome.receipt",
                        with(&|o| o.receipt = "cd34".into()),
                    ),
                    (
                        "dsr_run.outcome.disposition",
                        with(&|o| o.disposition = ExecutionDisposition::Failed),
                    ),
                    (
                        "dsr_run.outcome.predicate",
                        with(&|o| o.predicate = PredicateOutcome::Refuted),
                    ),
                    (
                        "dsr_run.outcome.evidence_methods",
                        with(&|o| o.evidence_methods = "oracle".into()),
                    ),
                    (
                        "dsr_run.outcome.grade",
                        with(&|o| o.grade = EpistemicGrade::Estimated),
                    ),
                    (
                        "dsr_run.outcome.applicability",
                        with(&|o| o.applicability = DomainApplicability::OutOfDomain),
                    ),
                    (
                        "dsr_run.outcome.support",
                        with(&|o| o.support = OperationalSupport::Degraded),
                    ),
                    (
                        "dsr_run.outcome.completeness",
                        with(&|o| o.completeness = EvidenceCompleteness::Incomplete),
                    ),
                    (
                        "dsr_run.outcome.integrity",
                        with(&|o| o.integrity = EvidenceIntegrity::Unchecked),
                    ),
                    (
                        "dsr_run.outcome.promotion",
                        with(&|o| o.promotion = PromotionEffect::Promoted),
                    ),
                    (
                        "dsr_run.outcome.detail",
                        with(&|o| o.detail = "different".into()),
                    ),
                ],
                &mut observed,
            );
        }
        {
            let base_outcome = || ScopedReceiptOutcome {
                scope: ReceiptScope::Campaign,
                receipt: "ef56".into(),
                disposition: ExecutionDisposition::Completed,
                predicate: PredicateOutcome::Satisfied,
                evidence_methods: "oracle,enclosure".into(),
                grade: EpistemicGrade::Verified,
                applicability: DomainApplicability::InDomain,
                support: OperationalSupport::Supported,
                completeness: EvidenceCompleteness::Complete,
                integrity: EvidenceIntegrity::Intact,
                promotion: PromotionEffect::Unchanged,
                detail: "ok".into(),
            };
            let with = |mutate: &dyn Fn(&mut ScopedReceiptOutcome)| {
                let mut o = base_outcome();
                mutate(&mut o);
                EventKind::CampaignRun {
                    campaign: "camp-1".into(),
                    outcome: o,
                }
            };
            assert_payload_mutations(
                EventKind::CampaignRun {
                    campaign: "camp-1".into(),
                    outcome: base_outcome(),
                },
                vec![
                    (
                        "campaign_run.campaign",
                        EventKind::CampaignRun {
                            campaign: "camp-2".into(),
                            outcome: base_outcome(),
                        },
                    ),
                    (
                        "campaign_run.outcome.scope",
                        with(&|o| o.scope = ReceiptScope::Job),
                    ),
                    (
                        "campaign_run.outcome.receipt",
                        with(&|o| o.receipt = "0011".into()),
                    ),
                    (
                        "campaign_run.outcome.disposition",
                        with(&|o| o.disposition = ExecutionDisposition::Cancelled),
                    ),
                    (
                        "campaign_run.outcome.predicate",
                        with(&|o| o.predicate = PredicateOutcome::Indeterminate),
                    ),
                    (
                        "campaign_run.outcome.evidence_methods",
                        with(&|o| o.evidence_methods = "enclosure".into()),
                    ),
                    (
                        "campaign_run.outcome.grade",
                        with(&|o| o.grade = EpistemicGrade::Validated),
                    ),
                    (
                        "campaign_run.outcome.applicability",
                        with(&|o| o.applicability = DomainApplicability::Unevaluated),
                    ),
                    (
                        "campaign_run.outcome.support",
                        with(&|o| o.support = OperationalSupport::Unsupported),
                    ),
                    (
                        "campaign_run.outcome.completeness",
                        with(&|o| o.completeness = EvidenceCompleteness::Incomplete),
                    ),
                    (
                        "campaign_run.outcome.integrity",
                        with(&|o| o.integrity = EvidenceIntegrity::Compromised),
                    ),
                    (
                        "campaign_run.outcome.promotion",
                        with(&|o| o.promotion = PromotionEffect::Demoted),
                    ),
                    (
                        "campaign_run.outcome.detail",
                        with(&|o| o.detail = "different".into()),
                    ),
                ],
                &mut observed,
            );
        }
        assert_payload_mutations(
            EventKind::SubmissionDecision {
                job: "job-1".into(),
                lane: SubmissionLane::Durable,
                outcome: SubmissionOutcome::Admitted,
                detail: "ok".into(),
            },
            vec![
                (
                    "submission_decision.job",
                    EventKind::SubmissionDecision {
                        job: "job-2".into(),
                        lane: SubmissionLane::Durable,
                        outcome: SubmissionOutcome::Admitted,
                        detail: "ok".into(),
                    },
                ),
                (
                    "submission_decision.lane",
                    EventKind::SubmissionDecision {
                        job: "job-1".into(),
                        lane: SubmissionLane::Plain,
                        outcome: SubmissionOutcome::Admitted,
                        detail: "ok".into(),
                    },
                ),
                (
                    "submission_decision.outcome",
                    EventKind::SubmissionDecision {
                        job: "job-1".into(),
                        lane: SubmissionLane::Durable,
                        outcome: SubmissionOutcome::Deferred,
                        detail: "ok".into(),
                    },
                ),
                (
                    "submission_decision.detail",
                    EventKind::SubmissionDecision {
                        job: "job-1".into(),
                        lane: SubmissionLane::Durable,
                        outcome: SubmissionOutcome::Admitted,
                        detail: "different".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::WorkIdentityBinding {
                job: "job-1".into(),
                attempt: 1,
                operation: "op-a".into(),
            },
            vec![
                (
                    "work_identity_binding.job",
                    EventKind::WorkIdentityBinding {
                        job: "job-2".into(),
                        attempt: 1,
                        operation: "op-a".into(),
                    },
                ),
                (
                    "work_identity_binding.attempt",
                    EventKind::WorkIdentityBinding {
                        job: "job-1".into(),
                        attempt: 2,
                        operation: "op-a".into(),
                    },
                ),
                (
                    "work_identity_binding.operation",
                    EventKind::WorkIdentityBinding {
                        job: "job-1".into(),
                        attempt: 1,
                        operation: "op-b".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::LeaseCursor {
                queue: "dsr".into(),
                lease: "lease-1".into(),
                holder: "worker-1".into(),
                cursor: 10,
            },
            vec![
                (
                    "lease_cursor.queue",
                    EventKind::LeaseCursor {
                        queue: "storm".into(),
                        lease: "lease-1".into(),
                        holder: "worker-1".into(),
                        cursor: 10,
                    },
                ),
                (
                    "lease_cursor.lease",
                    EventKind::LeaseCursor {
                        queue: "dsr".into(),
                        lease: "lease-2".into(),
                        holder: "worker-1".into(),
                        cursor: 10,
                    },
                ),
                (
                    "lease_cursor.holder",
                    EventKind::LeaseCursor {
                        queue: "dsr".into(),
                        lease: "lease-1".into(),
                        holder: "worker-2".into(),
                        cursor: 10,
                    },
                ),
                (
                    "lease_cursor.cursor",
                    EventKind::LeaseCursor {
                        queue: "dsr".into(),
                        lease: "lease-1".into(),
                        holder: "worker-1".into(),
                        cursor: 11,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::AttachDetach {
                observer: "workbench-1".into(),
                target: "camp-1".into(),
                action: AttachAction::Attach,
            },
            vec![
                (
                    "attach_detach.observer",
                    EventKind::AttachDetach {
                        observer: "workbench-2".into(),
                        target: "camp-1".into(),
                        action: AttachAction::Attach,
                    },
                ),
                (
                    "attach_detach.target",
                    EventKind::AttachDetach {
                        observer: "workbench-1".into(),
                        target: "camp-2".into(),
                        action: AttachAction::Attach,
                    },
                ),
                (
                    "attach_detach.action",
                    EventKind::AttachDetach {
                        observer: "workbench-1".into(),
                        target: "camp-1".into(),
                        action: AttachAction::Detach,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::JourneyPhase {
                journey: "thermal-fatigue-gearbox".into(),
                phase: "solve".into(),
                ordinal: 3,
            },
            vec![
                (
                    "journey_phase.journey",
                    EventKind::JourneyPhase {
                        journey: "bracket-topology".into(),
                        phase: "solve".into(),
                        ordinal: 3,
                    },
                ),
                (
                    "journey_phase.phase",
                    EventKind::JourneyPhase {
                        journey: "thermal-fatigue-gearbox".into(),
                        phase: "adjudicate".into(),
                        ordinal: 3,
                    },
                ),
                (
                    "journey_phase.ordinal",
                    EventKind::JourneyPhase {
                        journey: "thermal-fatigue-gearbox".into(),
                        phase: "solve".into(),
                        ordinal: 4,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::ScopeTileProgress {
                completed: 10,
                total: 64,
            },
            vec![
                (
                    "scope_tile_progress.completed",
                    EventKind::ScopeTileProgress {
                        completed: 11,
                        total: 64,
                    },
                ),
                (
                    "scope_tile_progress.total",
                    EventKind::ScopeTileProgress {
                        completed: 10,
                        total: 65,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::Heartbeat {
                worker: "worker-1".into(),
                beat: 7,
            },
            vec![
                (
                    "heartbeat.worker",
                    EventKind::Heartbeat {
                        worker: "worker-2".into(),
                        beat: 7,
                    },
                ),
                (
                    "heartbeat.beat",
                    EventKind::Heartbeat {
                        worker: "worker-1".into(),
                        beat: 8,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::ObservationGap {
                expected_seq: 40,
                resumed_seq: 52,
                reason: "sink rotated under memory pressure".into(),
            },
            vec![
                (
                    "observation_gap.expected_seq",
                    EventKind::ObservationGap {
                        expected_seq: 41,
                        resumed_seq: 52,
                        reason: "sink rotated under memory pressure".into(),
                    },
                ),
                (
                    "observation_gap.resumed_seq",
                    EventKind::ObservationGap {
                        expected_seq: 40,
                        resumed_seq: 53,
                        reason: "sink rotated under memory pressure".into(),
                    },
                ),
                (
                    "observation_gap.reason",
                    EventKind::ObservationGap {
                        expected_seq: 40,
                        resumed_seq: 52,
                        reason: "different".into(),
                    },
                ),
            ],
            &mut observed,
        );
        {
            let process_id =
                |raw: &str| process::ProcessId::new(raw).expect("identity fixture process id");
            let artifact = || {
                process::DurableArtifactPointer::committed(
                    "artifact://process/worker-1/stdout/7",
                    "0123456789abcdef0123456789abcdef",
                )
                .expect("identity fixture artifact")
            };
            let base = || EventKind::ProcessFrame {
                process: process_id("worker-1"),
                stream: process::ProcessStream::Stdout,
                ordinal: 7,
                class: process::ProcessEventClass::Diagnostic,
                payload: vec![0, 0xff, b'\n'],
                original_bytes: 8,
                artifact: Some(artifact()),
            };
            let field = |mutate: &dyn Fn(&mut EventKind)| {
                let mut kind = base();
                mutate(&mut kind);
                kind
            };
            assert_payload_mutations(
                base(),
                vec![
                    (
                        "process_frame.process",
                        field(&|kind| {
                            if let EventKind::ProcessFrame { process, .. } = kind {
                                *process = process_id("worker-2");
                            }
                        }),
                    ),
                    (
                        "process_frame.stream",
                        field(&|kind| {
                            if let EventKind::ProcessFrame { stream, .. } = kind {
                                *stream = process::ProcessStream::Stderr;
                            }
                        }),
                    ),
                    (
                        "process_frame.ordinal",
                        field(&|kind| {
                            if let EventKind::ProcessFrame { ordinal, .. } = kind {
                                *ordinal = 8;
                            }
                        }),
                    ),
                    (
                        "process_frame.class",
                        field(&|kind| {
                            if let EventKind::ProcessFrame { class, .. } = kind {
                                *class = process::ProcessEventClass::Telemetry;
                            }
                        }),
                    ),
                    (
                        "process_frame.payload",
                        field(&|kind| {
                            if let EventKind::ProcessFrame { payload, .. } = kind {
                                payload.push(1);
                            }
                        }),
                    ),
                    (
                        "process_frame.original_bytes",
                        field(&|kind| {
                            if let EventKind::ProcessFrame { original_bytes, .. } = kind {
                                *original_bytes = 9;
                            }
                        }),
                    ),
                    (
                        "process_frame.artifact_pointer",
                        field(&|kind| {
                            if let EventKind::ProcessFrame { artifact, .. } = kind {
                                *artifact = None;
                            }
                        }),
                    ),
                ],
                &mut observed,
            );
        }
        {
            let process_id =
                |raw: &str| process::ProcessId::new(raw).expect("identity fixture process id");
            let artifact = || {
                process::DurableArtifactPointer::committed(
                    "artifact://process/worker-1/stderr/9-11",
                    "fedcba9876543210fedcba9876543210",
                )
                .expect("identity fixture artifact")
            };
            let base = || EventKind::ProcessGap {
                process: process_id("worker-1"),
                stream: process::ProcessStream::Stderr,
                class: process::ProcessEventClass::Telemetry,
                first_ordinal: 9,
                last_ordinal: 11,
                original_count: 3,
                emitted_count: 1,
                dropped_count: 2,
                omitted_inline_bytes: 144,
                reason: process::LossReason::PolicySampling,
                policy_version: 7,
                artifact: Some(artifact()),
            };
            let field = |mutate: &dyn Fn(&mut EventKind)| {
                let mut kind = base();
                mutate(&mut kind);
                kind
            };
            assert_payload_mutations(
                base(),
                vec![
                    (
                        "process_gap.process",
                        field(&|kind| {
                            if let EventKind::ProcessGap { process, .. } = kind {
                                *process = process_id("worker-2");
                            }
                        }),
                    ),
                    (
                        "process_gap.stream",
                        field(&|kind| {
                            if let EventKind::ProcessGap { stream, .. } = kind {
                                *stream = process::ProcessStream::Stdout;
                            }
                        }),
                    ),
                    (
                        "process_gap.class",
                        field(&|kind| {
                            if let EventKind::ProcessGap { class, .. } = kind {
                                *class = process::ProcessEventClass::Diagnostic;
                            }
                        }),
                    ),
                    (
                        "process_gap.first_ordinal",
                        field(&|kind| {
                            if let EventKind::ProcessGap { first_ordinal, .. } = kind {
                                *first_ordinal = 8;
                            }
                        }),
                    ),
                    (
                        "process_gap.last_ordinal",
                        field(&|kind| {
                            if let EventKind::ProcessGap { last_ordinal, .. } = kind {
                                *last_ordinal = 12;
                            }
                        }),
                    ),
                    (
                        "process_gap.original_count",
                        field(&|kind| {
                            if let EventKind::ProcessGap { original_count, .. } = kind {
                                *original_count = 4;
                            }
                        }),
                    ),
                    (
                        "process_gap.emitted_count",
                        field(&|kind| {
                            if let EventKind::ProcessGap { emitted_count, .. } = kind {
                                *emitted_count = 2;
                            }
                        }),
                    ),
                    (
                        "process_gap.dropped_count",
                        field(&|kind| {
                            if let EventKind::ProcessGap { dropped_count, .. } = kind {
                                *dropped_count = 3;
                            }
                        }),
                    ),
                    (
                        "process_gap.omitted_inline_bytes",
                        field(&|kind| {
                            if let EventKind::ProcessGap {
                                omitted_inline_bytes,
                                ..
                            } = kind
                            {
                                *omitted_inline_bytes = 145;
                            }
                        }),
                    ),
                    (
                        "process_gap.reason",
                        field(&|kind| {
                            if let EventKind::ProcessGap { reason, .. } = kind {
                                *reason = process::LossReason::QueueCapacity;
                            }
                        }),
                    ),
                    (
                        "process_gap.policy_version",
                        field(&|kind| {
                            if let EventKind::ProcessGap { policy_version, .. } = kind {
                                *policy_version = 8;
                            }
                        }),
                    ),
                    (
                        "process_gap.artifact_pointer",
                        field(&|kind| {
                            if let EventKind::ProcessGap { artifact, .. } = kind {
                                *artifact = None;
                            }
                        }),
                    ),
                ],
                &mut observed,
            );
        }
        {
            let base = || EventKind::OracleComparison {
                oracle: "mms-poisson".into(),
                subject: "ab12".into(),
                metric: "l2-rel".into(),
                observed: 1.5e-9,
                tolerance: 1e-8,
                pass: true,
                detail: "ok".into(),
            };
            let field = |mutate: &dyn Fn(&mut EventKind)| {
                let mut k = base();
                mutate(&mut k);
                k
            };
            assert_payload_mutations(
                base(),
                vec![
                    (
                        "oracle_comparison.oracle",
                        field(&|k| {
                            if let EventKind::OracleComparison { oracle, .. } = k {
                                *oracle = "nafems-le1".into();
                            }
                        }),
                    ),
                    (
                        "oracle_comparison.subject",
                        field(&|k| {
                            if let EventKind::OracleComparison { subject, .. } = k {
                                *subject = "cd34".into();
                            }
                        }),
                    ),
                    (
                        "oracle_comparison.metric",
                        field(&|k| {
                            if let EventKind::OracleComparison { metric, .. } = k {
                                *metric = "max-abs".into();
                            }
                        }),
                    ),
                    (
                        "oracle_comparison.observed",
                        field(&|k| {
                            if let EventKind::OracleComparison { observed, .. } = k {
                                *observed = 2.5e-9;
                            }
                        }),
                    ),
                    (
                        "oracle_comparison.tolerance",
                        field(&|k| {
                            if let EventKind::OracleComparison { tolerance, .. } = k {
                                *tolerance = 1e-7;
                            }
                        }),
                    ),
                    (
                        "oracle_comparison.pass",
                        field(&|k| {
                            if let EventKind::OracleComparison { pass, .. } = k {
                                *pass = false;
                            }
                        }),
                    ),
                    (
                        "oracle_comparison.detail",
                        field(&|k| {
                            if let EventKind::OracleComparison { detail, .. } = k {
                                *detail = "different".into();
                            }
                        }),
                    ),
                ],
                &mut observed,
            );
        }
        assert_payload_mutations(
            EventKind::ToleranceDerivation {
                quantity: "displacement-l2".into(),
                tolerance: 1e-8,
                basis: "mms convergence order".into(),
            },
            vec![
                (
                    "tolerance_derivation.quantity",
                    EventKind::ToleranceDerivation {
                        quantity: "stress-max".into(),
                        tolerance: 1e-8,
                        basis: "mms convergence order".into(),
                    },
                ),
                (
                    "tolerance_derivation.tolerance",
                    EventKind::ToleranceDerivation {
                        quantity: "displacement-l2".into(),
                        tolerance: 1e-7,
                        basis: "mms convergence order".into(),
                    },
                ),
                (
                    "tolerance_derivation.basis",
                    EventKind::ToleranceDerivation {
                        quantity: "displacement-l2".into(),
                        tolerance: 1e-8,
                        basis: "different".into(),
                    },
                ),
            ],
            &mut observed,
        );
        {
            let base_outcome = || ScopedReceiptOutcome {
                scope: ReceiptScope::Job,
                receipt: "9a0b".into(),
                disposition: ExecutionDisposition::Completed,
                predicate: PredicateOutcome::Satisfied,
                evidence_methods: "oracle".into(),
                grade: EpistemicGrade::Verified,
                applicability: DomainApplicability::InDomain,
                support: OperationalSupport::Supported,
                completeness: EvidenceCompleteness::Complete,
                integrity: EvidenceIntegrity::Intact,
                promotion: PromotionEffect::Promoted,
                detail: "ok".into(),
            };
            let with = |mutate: &dyn Fn(&mut ScopedReceiptOutcome)| {
                let mut o = base_outcome();
                mutate(&mut o);
                EventKind::ClaimAdjudication {
                    claim: "claim-1".into(),
                    outcome: o,
                }
            };
            assert_payload_mutations(
                EventKind::ClaimAdjudication {
                    claim: "claim-1".into(),
                    outcome: base_outcome(),
                },
                vec![
                    (
                        "claim_adjudication.claim",
                        EventKind::ClaimAdjudication {
                            claim: "claim-2".into(),
                            outcome: base_outcome(),
                        },
                    ),
                    (
                        "claim_adjudication.outcome.scope",
                        with(&|o| o.scope = ReceiptScope::Operation),
                    ),
                    (
                        "claim_adjudication.outcome.receipt",
                        with(&|o| o.receipt = "1c2d".into()),
                    ),
                    (
                        "claim_adjudication.outcome.disposition",
                        with(&|o| o.disposition = ExecutionDisposition::TimedOut),
                    ),
                    (
                        "claim_adjudication.outcome.predicate",
                        with(&|o| o.predicate = PredicateOutcome::NotApplicable),
                    ),
                    (
                        "claim_adjudication.outcome.evidence_methods",
                        with(&|o| o.evidence_methods = "enclosure".into()),
                    ),
                    (
                        "claim_adjudication.outcome.grade",
                        with(&|o| o.grade = EpistemicGrade::Reported),
                    ),
                    (
                        "claim_adjudication.outcome.applicability",
                        with(&|o| o.applicability = DomainApplicability::NotApplicable),
                    ),
                    (
                        "claim_adjudication.outcome.support",
                        with(&|o| o.support = OperationalSupport::NotApplicable),
                    ),
                    (
                        "claim_adjudication.outcome.completeness",
                        with(&|o| o.completeness = EvidenceCompleteness::NotApplicable),
                    ),
                    (
                        "claim_adjudication.outcome.integrity",
                        with(&|o| o.integrity = EvidenceIntegrity::NotApplicable),
                    ),
                    (
                        "claim_adjudication.outcome.promotion",
                        with(&|o| o.promotion = PromotionEffect::NotApplicable),
                    ),
                    (
                        "claim_adjudication.outcome.detail",
                        with(&|o| o.detail = "different".into()),
                    ),
                ],
                &mut observed,
            );
        }
        assert_payload_mutations(
            EventKind::CapabilityDomainDecision {
                capability: "time-varying-flux".into(),
                domain: "thermal-transient".into(),
                decision: CapabilityDecision::Admitted,
                detail: "ok".into(),
            },
            vec![
                (
                    "capability_domain_decision.capability",
                    EventKind::CapabilityDomainDecision {
                        capability: "adaptive-remesh".into(),
                        domain: "thermal-transient".into(),
                        decision: CapabilityDecision::Admitted,
                        detail: "ok".into(),
                    },
                ),
                (
                    "capability_domain_decision.domain",
                    EventKind::CapabilityDomainDecision {
                        capability: "time-varying-flux".into(),
                        domain: "structural-static".into(),
                        decision: CapabilityDecision::Admitted,
                        detail: "ok".into(),
                    },
                ),
                (
                    "capability_domain_decision.decision",
                    EventKind::CapabilityDomainDecision {
                        capability: "time-varying-flux".into(),
                        domain: "thermal-transient".into(),
                        decision: CapabilityDecision::Restricted,
                        detail: "ok".into(),
                    },
                ),
                (
                    "capability_domain_decision.detail",
                    EventKind::CapabilityDomainDecision {
                        capability: "time-varying-flux".into(),
                        domain: "thermal-transient".into(),
                        decision: CapabilityDecision::Admitted,
                        detail: "different".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::LifecycleTransition {
                entity: "job-1".into(),
                transition: LifecycleTransitionKind::Checkpoint,
                detail: "ok".into(),
            },
            vec![
                (
                    "lifecycle_transition.entity",
                    EventKind::LifecycleTransition {
                        entity: "job-2".into(),
                        transition: LifecycleTransitionKind::Checkpoint,
                        detail: "ok".into(),
                    },
                ),
                (
                    "lifecycle_transition.transition",
                    EventKind::LifecycleTransition {
                        entity: "job-1".into(),
                        transition: LifecycleTransitionKind::Drain,
                        detail: "ok".into(),
                    },
                ),
                (
                    "lifecycle_transition.detail",
                    EventKind::LifecycleTransition {
                        entity: "job-1".into(),
                        transition: LifecycleTransitionKind::Checkpoint,
                        detail: "different".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::ArtifactLifecycle {
                artifact: "ab12".into(),
                action: ArtifactAction::Commit,
                actor: "policy-retention".into(),
                detail: "ok".into(),
            },
            vec![
                (
                    "artifact_lifecycle.artifact",
                    EventKind::ArtifactLifecycle {
                        artifact: "cd34".into(),
                        action: ArtifactAction::Commit,
                        actor: "policy-retention".into(),
                        detail: "ok".into(),
                    },
                ),
                (
                    "artifact_lifecycle.action",
                    EventKind::ArtifactLifecycle {
                        artifact: "ab12".into(),
                        action: ArtifactAction::Verify,
                        actor: "policy-retention".into(),
                        detail: "ok".into(),
                    },
                ),
                (
                    "artifact_lifecycle.actor",
                    EventKind::ArtifactLifecycle {
                        artifact: "ab12".into(),
                        action: ArtifactAction::Commit,
                        actor: "user-7".into(),
                        detail: "ok".into(),
                    },
                ),
                (
                    "artifact_lifecycle.detail",
                    EventKind::ArtifactLifecycle {
                        artifact: "ab12".into(),
                        action: ArtifactAction::Commit,
                        actor: "policy-retention".into(),
                        detail: "different".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::VisualizationTransform {
                view: "stress-contour-p1".into(),
                transform: "warp-by-displacement x50".into(),
                source: "ab12".into(),
            },
            vec![
                (
                    "visualization_transform.view",
                    EventKind::VisualizationTransform {
                        view: "temp-slice-z0".into(),
                        transform: "warp-by-displacement x50".into(),
                        source: "ab12".into(),
                    },
                ),
                (
                    "visualization_transform.transform",
                    EventKind::VisualizationTransform {
                        view: "stress-contour-p1".into(),
                        transform: "isosurface 0.5".into(),
                        source: "ab12".into(),
                    },
                ),
                (
                    "visualization_transform.source",
                    EventKind::VisualizationTransform {
                        view: "stress-contour-p1".into(),
                        transform: "warp-by-displacement x50".into(),
                        source: "cd34".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::DiagnosticRepair {
                subject: "mesh-9".into(),
                action: "collapse slivers".into(),
                pass: true,
                detail: "ok".into(),
            },
            vec![
                (
                    "diagnostic_repair.subject",
                    EventKind::DiagnosticRepair {
                        subject: "mesh-10".into(),
                        action: "collapse slivers".into(),
                        pass: true,
                        detail: "ok".into(),
                    },
                ),
                (
                    "diagnostic_repair.action",
                    EventKind::DiagnosticRepair {
                        subject: "mesh-9".into(),
                        action: "split edges".into(),
                        pass: true,
                        detail: "ok".into(),
                    },
                ),
                (
                    "diagnostic_repair.pass",
                    EventKind::DiagnosticRepair {
                        subject: "mesh-9".into(),
                        action: "collapse slivers".into(),
                        pass: false,
                        detail: "ok".into(),
                    },
                ),
                (
                    "diagnostic_repair.detail",
                    EventKind::DiagnosticRepair {
                        subject: "mesh-9".into(),
                        action: "collapse slivers".into(),
                        pass: true,
                        detail: "different".into(),
                    },
                ),
            ],
            &mut observed,
        );
        {
            let base = || EventKind::ContainmentNode {
                attempt: "attempt-9".into(),
                role: "tile".into(),
                node: "tile-42".into(),
                parent_role: "scope".into(),
                parent: "kernel-lbm".into(),
                seq: 3,
                dsr_run: "dsr-7".into(),
                campaign_run: "camp-1".into(),
                shard: "shard-3".into(),
                journey: "thermal-fatigue".into(),
                case: "case-12".into(),
            };
            let field = |mutate: &dyn Fn(&mut EventKind)| {
                let mut k = base();
                mutate(&mut k);
                k
            };
            let set = |name: &'static str, value: &'static str| {
                field(&move |k| {
                    if let EventKind::ContainmentNode {
                        attempt,
                        role,
                        node,
                        parent_role,
                        parent,
                        dsr_run,
                        campaign_run,
                        shard,
                        journey,
                        case,
                        ..
                    } = k
                    {
                        *match name {
                            "attempt" => attempt,
                            "role" => role,
                            "node" => node,
                            "parent_role" => parent_role,
                            "parent" => parent,
                            "dsr_run" => dsr_run,
                            "campaign_run" => campaign_run,
                            "shard" => shard,
                            "journey" => journey,
                            "case" => case,
                            _ => unreachable!("unknown containment field {name}"),
                        } = value.into();
                    }
                })
            };
            assert_payload_mutations(
                base(),
                vec![
                    ("containment_node.attempt", set("attempt", "attempt-10")),
                    ("containment_node.role", set("role", "op")),
                    ("containment_node.node", set("node", "tile-43")),
                    ("containment_node.parent_role", set("parent_role", "op")),
                    ("containment_node.parent", set("parent", "solve")),
                    (
                        "containment_node.seq",
                        field(&|k| {
                            if let EventKind::ContainmentNode { seq, .. } = k {
                                *seq = 4;
                            }
                        }),
                    ),
                    ("containment_node.dsr_run", set("dsr_run", "dsr-8")),
                    (
                        "containment_node.campaign_run",
                        set("campaign_run", "camp-2"),
                    ),
                    ("containment_node.shard", set("shard", "shard-4")),
                    ("containment_node.journey", set("journey", "bracket")),
                    ("containment_node.case", set("case", "case-13")),
                ],
                &mut observed,
            );
        }
        assert_payload_mutations(
            EventKind::ContainmentGap {
                attempt: "attempt-9".into(),
                node_role: "tile".into(),
                node: "tile-1".into(),
                missing_parent_role: "scope".into(),
                missing_parent: "gone".into(),
            },
            vec![
                (
                    "containment_gap.attempt",
                    EventKind::ContainmentGap {
                        attempt: "attempt-10".into(),
                        node_role: "tile".into(),
                        node: "tile-1".into(),
                        missing_parent_role: "scope".into(),
                        missing_parent: "gone".into(),
                    },
                ),
                (
                    "containment_gap.node_role",
                    EventKind::ContainmentGap {
                        attempt: "attempt-9".into(),
                        node_role: "op".into(),
                        node: "tile-1".into(),
                        missing_parent_role: "scope".into(),
                        missing_parent: "gone".into(),
                    },
                ),
                (
                    "containment_gap.node",
                    EventKind::ContainmentGap {
                        attempt: "attempt-9".into(),
                        node_role: "tile".into(),
                        node: "tile-2".into(),
                        missing_parent_role: "scope".into(),
                        missing_parent: "gone".into(),
                    },
                ),
                (
                    "containment_gap.missing_parent_role",
                    EventKind::ContainmentGap {
                        attempt: "attempt-9".into(),
                        node_role: "tile".into(),
                        node: "tile-1".into(),
                        missing_parent_role: "op".into(),
                        missing_parent: "gone".into(),
                    },
                ),
                (
                    "containment_gap.missing_parent",
                    EventKind::ContainmentGap {
                        attempt: "attempt-9".into(),
                        node_role: "tile".into(),
                        node: "tile-1".into(),
                        missing_parent_role: "scope".into(),
                        missing_parent: "elsewhere".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::RaceRecord {
                resource: "arena-slot".into(),
                schedule: "a-b-a".into(),
                pass: true,
                seed: 31,
            },
            vec![
                (
                    "race_record.resource",
                    EventKind::RaceRecord {
                        resource: "tune-row".into(),
                        schedule: "a-b-a".into(),
                        pass: true,
                        seed: 31,
                    },
                ),
                (
                    "race_record.schedule",
                    EventKind::RaceRecord {
                        resource: "arena-slot".into(),
                        schedule: "b-a-b".into(),
                        pass: true,
                        seed: 31,
                    },
                ),
                (
                    "race_record.pass",
                    EventKind::RaceRecord {
                        resource: "arena-slot".into(),
                        schedule: "a-b-a".into(),
                        pass: false,
                        seed: 31,
                    },
                ),
                (
                    "race_record.seed",
                    EventKind::RaceRecord {
                        resource: "arena-slot".into(),
                        schedule: "a-b-a".into(),
                        pass: true,
                        seed: 32,
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::DegradationEvent {
                resource: "retained_bytes_per_scope".into(),
                limit: 4096,
                observed: 5000,
                action: "rotate-lane".into(),
            },
            vec![
                (
                    "degradation_event.resource",
                    EventKind::DegradationEvent {
                        resource: "events_per_scope".into(),
                        limit: 4096,
                        observed: 5000,
                        action: "rotate-lane".into(),
                    },
                ),
                (
                    "degradation_event.limit",
                    EventKind::DegradationEvent {
                        resource: "retained_bytes_per_scope".into(),
                        limit: 8192,
                        observed: 5000,
                        action: "rotate-lane".into(),
                    },
                ),
                (
                    "degradation_event.observed",
                    EventKind::DegradationEvent {
                        resource: "retained_bytes_per_scope".into(),
                        limit: 4096,
                        observed: 6000,
                        action: "rotate-lane".into(),
                    },
                ),
                (
                    "degradation_event.action",
                    EventKind::DegradationEvent {
                        resource: "retained_bytes_per_scope".into(),
                        limit: 4096,
                        observed: 5000,
                        action: "refuse".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::ImportReceipt {
                format: "step".into(),
                artifact: "ab12".into(),
                accepted: true,
                detail: "clean".into(),
            },
            vec![
                (
                    "import_receipt.format",
                    EventKind::ImportReceipt {
                        format: "gltf".into(),
                        artifact: "ab12".into(),
                        accepted: true,
                        detail: "clean".into(),
                    },
                ),
                (
                    "import_receipt.artifact",
                    EventKind::ImportReceipt {
                        format: "step".into(),
                        artifact: "cd34".into(),
                        accepted: true,
                        detail: "clean".into(),
                    },
                ),
                (
                    "import_receipt.accepted",
                    EventKind::ImportReceipt {
                        format: "step".into(),
                        artifact: "ab12".into(),
                        accepted: false,
                        detail: "clean".into(),
                    },
                ),
                (
                    "import_receipt.detail",
                    EventKind::ImportReceipt {
                        format: "step".into(),
                        artifact: "ab12".into(),
                        accepted: true,
                        detail: "quarantined".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::CertificateVerdict {
                certificate: "ivl-enclosure".into(),
                pass: true,
                bound: 0.125,
                detail: "ok".into(),
            },
            vec![
                (
                    "certificate_verdict.certificate",
                    EventKind::CertificateVerdict {
                        certificate: "trim-cert".into(),
                        pass: true,
                        bound: 0.125,
                        detail: "ok".into(),
                    },
                ),
                (
                    "certificate_verdict.pass",
                    EventKind::CertificateVerdict {
                        certificate: "ivl-enclosure".into(),
                        pass: false,
                        bound: 0.125,
                        detail: "ok".into(),
                    },
                ),
                (
                    "certificate_verdict.bound",
                    EventKind::CertificateVerdict {
                        certificate: "ivl-enclosure".into(),
                        pass: true,
                        bound: 0.25,
                        detail: "ok".into(),
                    },
                ),
                (
                    "certificate_verdict.detail",
                    EventKind::CertificateVerdict {
                        certificate: "ivl-enclosure".into(),
                        pass: true,
                        bound: 0.125,
                        detail: "different".into(),
                    },
                ),
            ],
            &mut observed,
        );
        assert_payload_mutations(
            EventKind::Custom {
                name: "opaque".into(),
                json: r#"{"a":1,"b":2}"#.into(),
            },
            vec![
                (
                    "custom.name",
                    EventKind::Custom {
                        name: "other".into(),
                        json: r#"{"a":1,"b":2}"#.into(),
                    },
                ),
                (
                    "custom.json",
                    EventKind::Custom {
                        name: "opaque".into(),
                        json: r#"{"b":2,"a":1}"#.into(),
                    },
                ),
            ],
            &mut observed,
        );

        let expected = [
            "solver_residual.solver",
            "solver_residual.iter",
            "solver_residual.residual",
            "tile_complete.tile",
            "tile_complete.kernel",
            "cancellation.reason",
            "budget_delta.resource",
            "budget_delta.spent",
            "budget_delta.remaining",
            "gradient_check.op",
            "gradient_check.max_rel_err",
            "gradient_check.pass",
            "conformance_case.suite",
            "conformance_case.case",
            "conformance_case.pass",
            "conformance_case.detail",
            "conformance_case.seed",
            "benchmark_result.kernel",
            "benchmark_result.metric",
            "benchmark_result.value",
            "benchmark_result.machine",
            "storm_assertion.name",
            "storm_assertion.pass",
            "storm_assertion.seed",
            "manifest_selection.manifest",
            "manifest_selection.source_snapshot",
            "stratum_expansion.stratum",
            "stratum_expansion.profile",
            "stratum_expansion.cases",
            "dsr_run.run",
            "dsr_run.outcome.scope",
            "dsr_run.outcome.receipt",
            "dsr_run.outcome.disposition",
            "dsr_run.outcome.predicate",
            "dsr_run.outcome.evidence_methods",
            "dsr_run.outcome.grade",
            "dsr_run.outcome.applicability",
            "dsr_run.outcome.support",
            "dsr_run.outcome.completeness",
            "dsr_run.outcome.integrity",
            "dsr_run.outcome.promotion",
            "dsr_run.outcome.detail",
            "campaign_run.campaign",
            "campaign_run.outcome.scope",
            "campaign_run.outcome.receipt",
            "campaign_run.outcome.disposition",
            "campaign_run.outcome.predicate",
            "campaign_run.outcome.evidence_methods",
            "campaign_run.outcome.grade",
            "campaign_run.outcome.applicability",
            "campaign_run.outcome.support",
            "campaign_run.outcome.completeness",
            "campaign_run.outcome.integrity",
            "campaign_run.outcome.promotion",
            "campaign_run.outcome.detail",
            "submission_decision.job",
            "submission_decision.lane",
            "submission_decision.outcome",
            "submission_decision.detail",
            "work_identity_binding.job",
            "work_identity_binding.attempt",
            "work_identity_binding.operation",
            "lease_cursor.queue",
            "lease_cursor.lease",
            "lease_cursor.holder",
            "lease_cursor.cursor",
            "attach_detach.observer",
            "attach_detach.target",
            "attach_detach.action",
            "journey_phase.journey",
            "journey_phase.phase",
            "journey_phase.ordinal",
            "scope_tile_progress.completed",
            "scope_tile_progress.total",
            "heartbeat.worker",
            "heartbeat.beat",
            "observation_gap.expected_seq",
            "observation_gap.resumed_seq",
            "observation_gap.reason",
            "process_frame.process",
            "process_frame.stream",
            "process_frame.ordinal",
            "process_frame.class",
            "process_frame.payload",
            "process_frame.original_bytes",
            "process_frame.artifact_pointer",
            "process_gap.process",
            "process_gap.stream",
            "process_gap.class",
            "process_gap.first_ordinal",
            "process_gap.last_ordinal",
            "process_gap.original_count",
            "process_gap.emitted_count",
            "process_gap.dropped_count",
            "process_gap.omitted_inline_bytes",
            "process_gap.reason",
            "process_gap.policy_version",
            "process_gap.artifact_pointer",
            "oracle_comparison.oracle",
            "oracle_comparison.subject",
            "oracle_comparison.metric",
            "oracle_comparison.observed",
            "oracle_comparison.tolerance",
            "oracle_comparison.pass",
            "oracle_comparison.detail",
            "tolerance_derivation.quantity",
            "tolerance_derivation.tolerance",
            "tolerance_derivation.basis",
            "claim_adjudication.claim",
            "claim_adjudication.outcome.scope",
            "claim_adjudication.outcome.receipt",
            "claim_adjudication.outcome.disposition",
            "claim_adjudication.outcome.predicate",
            "claim_adjudication.outcome.evidence_methods",
            "claim_adjudication.outcome.grade",
            "claim_adjudication.outcome.applicability",
            "claim_adjudication.outcome.support",
            "claim_adjudication.outcome.completeness",
            "claim_adjudication.outcome.integrity",
            "claim_adjudication.outcome.promotion",
            "claim_adjudication.outcome.detail",
            "capability_domain_decision.capability",
            "capability_domain_decision.domain",
            "capability_domain_decision.decision",
            "capability_domain_decision.detail",
            "lifecycle_transition.entity",
            "lifecycle_transition.transition",
            "lifecycle_transition.detail",
            "artifact_lifecycle.artifact",
            "artifact_lifecycle.action",
            "artifact_lifecycle.actor",
            "artifact_lifecycle.detail",
            "visualization_transform.view",
            "visualization_transform.transform",
            "visualization_transform.source",
            "diagnostic_repair.subject",
            "diagnostic_repair.action",
            "diagnostic_repair.pass",
            "diagnostic_repair.detail",
            "containment_node.attempt",
            "containment_node.role",
            "containment_node.node",
            "containment_node.parent_role",
            "containment_node.parent",
            "containment_node.seq",
            "containment_node.dsr_run",
            "containment_node.campaign_run",
            "containment_node.shard",
            "containment_node.journey",
            "containment_node.case",
            "containment_gap.attempt",
            "containment_gap.node_role",
            "containment_gap.node",
            "containment_gap.missing_parent_role",
            "containment_gap.missing_parent",
            "race_record.resource",
            "race_record.schedule",
            "race_record.pass",
            "race_record.seed",
            "degradation_event.resource",
            "degradation_event.limit",
            "degradation_event.observed",
            "degradation_event.action",
            "import_receipt.format",
            "import_receipt.artifact",
            "import_receipt.accepted",
            "import_receipt.detail",
            "certificate_verdict.certificate",
            "certificate_verdict.pass",
            "certificate_verdict.bound",
            "certificate_verdict.detail",
            "custom.name",
            "custom.json",
        ];
        assert_eq!(
            observed.as_slice(),
            expected.as_slice(),
            "all 173 payload fields stay enumerated"
        );
    }

    #[test]
    fn scope_tree_mirrors_enter_exit_and_refuses_forgery_and_imbalance() {
        let mut em = Emitter::new("study-x", "op-1");
        assert_eq!(em.current_scope(), "op-1");
        assert_eq!(em.scope_depth(), 0);

        em.enter_scope("kernel-cg").expect("child scope");
        em.enter_scope("tile-7").expect("grandchild scope");
        assert_eq!(em.current_scope(), "op-1/kernel-cg/tile-7");
        assert_eq!(em.scope_depth(), 2);
        let deep = em.emit(
            Severity::Info,
            EventKind::TileComplete {
                tile: 7,
                kernel: "cg".into(),
            },
            None,
        );
        assert_eq!(deep.scope, "op-1/kernel-cg/tile-7");

        em.exit_scope().expect("exit grandchild");
        assert_eq!(em.current_scope(), "op-1/kernel-cg");
        let mid = em.emit(
            Severity::Info,
            EventKind::Cancellation {
                reason: "budget".into(),
            },
            None,
        );
        assert_eq!(mid.scope, "op-1/kernel-cg");
        assert_eq!(
            mid.seq,
            deep.seq + 1,
            "seq is one monotone stream across scope moves"
        );

        em.exit_scope().expect("exit child");
        assert_eq!(em.current_scope(), "op-1");
        assert!(
            em.exit_scope().is_err(),
            "exit at the base scope must refuse (unbalanced enter/exit)"
        );

        assert!(em.enter_scope("").is_err(), "empty segment must refuse");
        assert!(
            em.enter_scope("a/b").is_err(),
            "path separator inside a segment is scope forgery"
        );
        assert!(
            em.enter_scope("a\u{7}b").is_err(),
            "control characters inside a segment must refuse"
        );
        assert_eq!(
            em.current_scope(),
            "op-1",
            "refused entries must not move the scope"
        );
    }

    #[test]
    fn scope_tree_replay_is_deterministic_and_wire_valid() {
        let run = || {
            let mut em = Emitter::new("study-x", "op-1");
            let mut lines = Vec::new();
            em.enter_scope("solve").expect("scope");
            lines.push(
                em.emit(
                    Severity::Info,
                    EventKind::SolverResidual {
                        solver: "cg".into(),
                        iter: 0,
                        residual: 0.5,
                    },
                    None,
                )
                .to_jsonl(),
            );
            em.enter_scope("tile-3").expect("scope");
            lines.push(
                em.emit(
                    Severity::Trace,
                    EventKind::TileComplete {
                        tile: 3,
                        kernel: "cg".into(),
                    },
                    None,
                )
                .to_jsonl(),
            );
            em.exit_scope().expect("exit");
            em.exit_scope().expect("exit");
            lines
        };
        let first = run();
        let second = run();
        assert_eq!(
            first, second,
            "identical scope walks replay bit-identically"
        );
        for line in &first {
            validate_line(line).expect("scope-tree events stay wire-valid");
        }
    }

    #[test]
    fn verification_outcome_kinds_are_wire_valid_and_lint_enforced() {
        let mut em = Emitter::new("study-x", "dsr-1");
        let selection = em.emit(
            Severity::Info,
            EventKind::ManifestSelection {
                manifest: "ab12".into(),
                source_snapshot: "cd34".into(),
            },
            None,
        );
        let expansion = em.emit(
            Severity::Info,
            EventKind::StratumExpansion {
                stratum: "core".into(),
                profile: "nightly".into(),
                cases: 40,
            },
            None,
        );
        let mut outcome = ScopedReceiptOutcome {
            scope: ReceiptScope::Job,
            receipt: "ab12".into(),
            disposition: ExecutionDisposition::Completed,
            predicate: PredicateOutcome::Satisfied,
            evidence_methods: "oracle,enclosure".into(),
            grade: EpistemicGrade::Verified,
            applicability: DomainApplicability::InDomain,
            support: OperationalSupport::Supported,
            completeness: EvidenceCompleteness::Complete,
            integrity: EvidenceIntegrity::Intact,
            promotion: PromotionEffect::Unchanged,
            detail: "ok".into(),
        };
        let run = em.emit(
            Severity::Info,
            EventKind::DsrRun {
                run: "run-1".into(),
                outcome: outcome.clone(),
            },
            Some(9),
        );
        for event in [&selection, &expansion, &run] {
            validate_line(&event.to_jsonl()).expect("verification kinds stay wire-valid");
            lint_failure_record(event).expect("conclusive outcomes pass the lint");
        }

        // A failed outcome without diagnostic detail refuses.
        outcome.disposition = ExecutionDisposition::Failed;
        outcome.detail = String::new();
        let bad = em.emit(
            Severity::Error,
            EventKind::DsrRun {
                run: "run-1".into(),
                outcome: outcome.clone(),
            },
            None,
        );
        assert!(
            lint_failure_record(&bad).is_err(),
            "failed outcome without detail must refuse"
        );
        outcome.detail = "worker lost lease".into();
        let repaired = em.emit(
            Severity::Error,
            EventKind::DsrRun {
                run: "run-1".into(),
                outcome,
            },
            None,
        );
        lint_failure_record(&repaired).expect("detailed failure passes");
        validate_line(&repaired.to_jsonl()).expect("failure line stays wire-valid");
    }

    #[test]
    // The end-to-end wire/lint cases share one emitter and outcome vocabulary;
    // splitting them would obscure which queue transitions are covered.
    #[allow(clippy::too_many_lines)]
    fn run_queue_identity_kinds_are_wire_valid_and_lint_enforced() {
        let mut em = Emitter::new("study-x", "campaign-1");
        let outcome = ScopedReceiptOutcome {
            scope: ReceiptScope::Campaign,
            receipt: "ef56".into(),
            disposition: ExecutionDisposition::Completed,
            predicate: PredicateOutcome::Satisfied,
            evidence_methods: "oracle".into(),
            grade: EpistemicGrade::Verified,
            applicability: DomainApplicability::InDomain,
            support: OperationalSupport::Supported,
            completeness: EvidenceCompleteness::Complete,
            integrity: EvidenceIntegrity::Intact,
            promotion: PromotionEffect::Promoted,
            detail: "ok".into(),
        };
        let conclusive = [
            em.emit(
                Severity::Info,
                EventKind::CampaignRun {
                    campaign: "camp-1".into(),
                    outcome: outcome.clone(),
                },
                Some(11),
            ),
            em.emit(
                Severity::Info,
                EventKind::SubmissionDecision {
                    job: "job-1".into(),
                    lane: SubmissionLane::Durable,
                    outcome: SubmissionOutcome::Admitted,
                    detail: String::new(),
                },
                None,
            ),
            em.emit(
                Severity::Trace,
                EventKind::WorkIdentityBinding {
                    job: "job-1".into(),
                    attempt: 1,
                    operation: "op-a".into(),
                },
                None,
            ),
            em.emit(
                Severity::Trace,
                EventKind::LeaseCursor {
                    queue: "dsr".into(),
                    lease: "lease-1".into(),
                    holder: "worker-1".into(),
                    cursor: 10,
                },
                None,
            ),
            em.emit(
                Severity::Info,
                EventKind::AttachDetach {
                    observer: "workbench-1".into(),
                    target: "camp-1".into(),
                    action: AttachAction::Attach,
                },
                None,
            ),
        ];
        for event in &conclusive {
            validate_line(&event.to_jsonl()).expect("run/queue/identity kinds stay wire-valid");
            lint_failure_record(event).expect("conclusive records pass the lint");
        }

        // An indeterminate campaign outcome without detail refuses.
        let mut indeterminate = outcome;
        indeterminate.disposition = ExecutionDisposition::Indeterminate;
        indeterminate.detail = String::new();
        let bad_campaign = em.emit(
            Severity::Error,
            EventKind::CampaignRun {
                campaign: "camp-1".into(),
                outcome: indeterminate,
            },
            None,
        );
        assert!(
            lint_failure_record(&bad_campaign).is_err(),
            "indeterminate campaign outcome without detail must refuse"
        );

        // Refused and deferred submissions without detail refuse; admitted
        // never demands one.
        for refused_like in [SubmissionOutcome::Refused, SubmissionOutcome::Deferred] {
            let bad = em.emit(
                Severity::Error,
                EventKind::SubmissionDecision {
                    job: "job-1".into(),
                    lane: SubmissionLane::Plain,
                    outcome: refused_like,
                    detail: String::new(),
                },
                None,
            );
            assert!(
                lint_failure_record(&bad).is_err(),
                "{} submission without detail must refuse",
                refused_like.name()
            );
        }
        let repaired = em.emit(
            Severity::Error,
            EventKind::SubmissionDecision {
                job: "job-1".into(),
                lane: SubmissionLane::Plain,
                outcome: SubmissionOutcome::Refused,
                detail: "queue at durable-flush capacity".into(),
            },
            None,
        );
        lint_failure_record(&repaired).expect("detailed refusal passes");
        validate_line(&repaired.to_jsonl()).expect("refusal line stays wire-valid");
    }

    #[test]
    fn execution_progress_kinds_ride_the_scope_tree_and_gaps_demand_reasons() {
        // One journey walk: phases and progress are emitted INSIDE scope-tree
        // spans, so the envelope scope carries the storyboard position.
        let walk = || -> Vec<String> {
            let mut em = Emitter::new("study-x", "journey-thermal-fatigue");
            let mut lines = Vec::new();
            lines.push(
                em.emit(
                    Severity::Info,
                    EventKind::JourneyPhase {
                        journey: "thermal-fatigue-gearbox".into(),
                        phase: "solve".into(),
                        ordinal: 2,
                    },
                    None,
                )
                .to_jsonl(),
            );
            em.enter_scope("kernel-lbm").expect("child scope");
            lines.push(
                em.emit(
                    Severity::Trace,
                    EventKind::ScopeTileProgress {
                        completed: 10,
                        total: 64,
                    },
                    None,
                )
                .to_jsonl(),
            );
            lines.push(
                em.emit(
                    Severity::Trace,
                    EventKind::Heartbeat {
                        worker: "worker-1".into(),
                        beat: 7,
                    },
                    None,
                )
                .to_jsonl(),
            );
            em.exit_scope().expect("balanced exit");
            lines
        };
        let first = walk();
        let second = walk();
        assert_eq!(first, second, "progress walks replay bit-identically");
        assert!(
            first[1].contains("journey-thermal-fatigue/kernel-lbm"),
            "tile progress carries its scope-tree position in the envelope"
        );
        for line in &first {
            validate_line(line).expect("progress kinds stay wire-valid");
        }

        // A gap with no reason refuses; a reasoned gap passes.
        let mut em = Emitter::new("study-x", "journey-thermal-fatigue");
        let bad = em.emit(
            Severity::Warn,
            EventKind::ObservationGap {
                expected_seq: 40,
                resumed_seq: 52,
                reason: String::new(),
            },
            None,
        );
        assert!(
            lint_failure_record(&bad).is_err(),
            "an unexplained observation gap must refuse"
        );
        let good = em.emit(
            Severity::Warn,
            EventKind::ObservationGap {
                expected_seq: 40,
                resumed_seq: 52,
                reason: "sink rotated under memory pressure".into(),
            },
            None,
        );
        lint_failure_record(&good).expect("reasoned gap passes");
        validate_line(&good.to_jsonl()).expect("gap line stays wire-valid");
    }

    #[test]
    // Keep the adjudication variants together as one exhaustive wire contract.
    #[allow(clippy::too_many_lines)]
    fn adjudication_kinds_are_wire_valid_and_lint_enforced() {
        let mut em = Emitter::new("study-x", "adjudicate-1");
        let outcome = ScopedReceiptOutcome {
            scope: ReceiptScope::Job,
            receipt: "9a0b".into(),
            disposition: ExecutionDisposition::Completed,
            predicate: PredicateOutcome::Refuted,
            evidence_methods: "oracle".into(),
            grade: EpistemicGrade::Verified,
            applicability: DomainApplicability::InDomain,
            support: OperationalSupport::Supported,
            completeness: EvidenceCompleteness::Complete,
            integrity: EvidenceIntegrity::Intact,
            promotion: PromotionEffect::Demoted,
            detail: "refuted by oracle within tolerance basis".into(),
        };
        let conclusive = [
            em.emit(
                Severity::Info,
                EventKind::OracleComparison {
                    oracle: "mms-poisson".into(),
                    subject: "ab12".into(),
                    metric: "l2-rel".into(),
                    observed: 1.5e-9,
                    tolerance: 1e-8,
                    pass: true,
                    detail: String::new(),
                },
                None,
            ),
            em.emit(
                Severity::Info,
                EventKind::ToleranceDerivation {
                    quantity: "displacement-l2".into(),
                    tolerance: 1e-8,
                    basis: "mms convergence order".into(),
                },
                None,
            ),
            em.emit(
                Severity::Info,
                EventKind::ClaimAdjudication {
                    claim: "claim-1".into(),
                    outcome,
                },
                Some(13),
            ),
            em.emit(
                Severity::Info,
                EventKind::CapabilityDomainDecision {
                    capability: "time-varying-flux".into(),
                    domain: "thermal-transient".into(),
                    decision: CapabilityDecision::Admitted,
                    detail: String::new(),
                },
                None,
            ),
        ];
        for event in &conclusive {
            validate_line(&event.to_jsonl()).expect("adjudication kinds stay wire-valid");
            lint_failure_record(event).expect("conclusive adjudications pass the lint");
        }

        // Every refusal-shaped adjudication record demands its detail.
        let bad_oracle = em.emit(
            Severity::Error,
            EventKind::OracleComparison {
                oracle: "mms-poisson".into(),
                subject: "ab12".into(),
                metric: "l2-rel".into(),
                observed: 3.0e-8,
                tolerance: 1e-8,
                pass: false,
                detail: String::new(),
            },
            None,
        );
        assert!(
            lint_failure_record(&bad_oracle).is_err(),
            "failing oracle comparison without detail must refuse"
        );
        let bad_basis = em.emit(
            Severity::Info,
            EventKind::ToleranceDerivation {
                quantity: "displacement-l2".into(),
                tolerance: 1e-8,
                basis: String::new(),
            },
            None,
        );
        assert!(
            lint_failure_record(&bad_basis).is_err(),
            "a tolerance without a named basis must refuse"
        );
        for decision in [CapabilityDecision::Refused, CapabilityDecision::Restricted] {
            let bad = em.emit(
                Severity::Error,
                EventKind::CapabilityDomainDecision {
                    capability: "time-varying-flux".into(),
                    domain: "thermal-transient".into(),
                    decision,
                    detail: String::new(),
                },
                None,
            );
            assert!(
                lint_failure_record(&bad).is_err(),
                "{} capability decision without detail must refuse",
                decision.name()
            );
        }
    }

    #[test]
    fn lifecycle_and_artifact_kinds_are_wire_valid_and_lint_enforced() {
        let mut em = Emitter::new("study-x", "campaign-1");
        let conclusive = [
            em.emit(
                Severity::Info,
                EventKind::LifecycleTransition {
                    entity: "job-1".into(),
                    transition: LifecycleTransitionKind::Checkpoint,
                    detail: String::new(),
                },
                None,
            ),
            em.emit(
                Severity::Info,
                EventKind::ArtifactLifecycle {
                    artifact: "ab12".into(),
                    action: ArtifactAction::Commit,
                    actor: "policy-retention".into(),
                    detail: String::new(),
                },
                None,
            ),
            em.emit(
                Severity::Info,
                EventKind::VisualizationTransform {
                    view: "stress-contour-p1".into(),
                    transform: "warp-by-displacement x50".into(),
                    source: "ab12".into(),
                },
                None,
            ),
            em.emit(
                Severity::Info,
                EventKind::DiagnosticRepair {
                    subject: "mesh-9".into(),
                    action: "collapse slivers".into(),
                    pass: true,
                    detail: String::new(),
                },
                None,
            ),
        ];
        for event in &conclusive {
            validate_line(&event.to_jsonl()).expect("lifecycle kinds stay wire-valid");
            lint_failure_record(event).expect("conclusive lifecycle records pass the lint");
        }

        // Faults, redactions/shares, and failed repairs demand their detail.
        let bad_fault = em.emit(
            Severity::Error,
            EventKind::LifecycleTransition {
                entity: "job-1".into(),
                transition: LifecycleTransitionKind::Fault,
                detail: String::new(),
            },
            None,
        );
        assert!(
            lint_failure_record(&bad_fault).is_err(),
            "a contained fault without detail must refuse"
        );
        for action in [ArtifactAction::Redact, ArtifactAction::Share] {
            let bad = em.emit(
                Severity::Warn,
                EventKind::ArtifactLifecycle {
                    artifact: "ab12".into(),
                    action,
                    actor: "user-7".into(),
                    detail: String::new(),
                },
                None,
            );
            assert!(
                lint_failure_record(&bad).is_err(),
                "{} without detail must refuse",
                action.name()
            );
        }
        let bad_repair = em.emit(
            Severity::Error,
            EventKind::DiagnosticRepair {
                subject: "mesh-9".into(),
                action: "collapse slivers".into(),
                pass: false,
                detail: String::new(),
            },
            None,
        );
        assert!(
            lint_failure_record(&bad_repair).is_err(),
            "an ineffective repair without detail must refuse"
        );
    }

    #[test]
    fn custom_payload_identity_is_exact_opaque_utf8() {
        let opaque = r#"{ "b":2, "a":1 }"#;
        let event = event_with_kind(EventKind::Custom {
            name: "opaque-json".into(),
            json: opaque.into(),
        });
        let identity = event.content_identity();
        let (tag, retained) = identity_field(identity.canonical_bytes(), "custom_json_opaque_utf8");
        assert_eq!(
            tag, 0x05,
            "custom JSON is bound as exact bytes, not text claiming canonical JSON"
        );
        assert_eq!(
            retained,
            opaque.as_bytes(),
            "opaque UTF-8 round-trips byte-for-byte through the identity frame"
        );
        assert!(event.to_jsonl().contains(&format!("\"data\":{opaque}")));
        validate_line(&event.to_jsonl())
            .expect("the valid opaque object remains a valid event line");

        let reordered = event_with_kind(EventKind::Custom {
            name: "opaque-json".into(),
            json: r#"{"a":1,"b":2}"#.into(),
        });
        assert_ne!(
            event.content_hash(),
            reordered.content_hash(),
            "whitespace and member order are honestly semantic under opaque-byte identity"
        );
    }

    #[test]
    fn event_content_identity_mutation_battery() {
        let base = Event {
            session: "session-a".into(),
            scope: "scope-a".into(),
            seq: 11,
            severity: Severity::Info,
            kind: EventKind::GradientCheck {
                op: "poisson".into(),
                max_rel_err: 0.25,
                pass: true,
            },
            wall_ns: None,
        };
        let base_hash = base.content_hash();
        let mutations = [
            Event {
                session: "session-b".into(),
                ..base.clone()
            },
            Event {
                scope: "scope-b".into(),
                ..base.clone()
            },
            Event {
                seq: 12,
                ..base.clone()
            },
            Event {
                severity: Severity::Warn,
                ..base.clone()
            },
            Event {
                kind: EventKind::Cancellation {
                    reason: "budget".into(),
                },
                ..base.clone()
            },
            Event {
                kind: EventKind::GradientCheck {
                    op: "poisson".into(),
                    max_rel_err: 0.25_f64.next_up(),
                    pass: true,
                },
                ..base.clone()
            },
        ];
        assert!(
            mutations
                .iter()
                .all(|mutation| mutation.content_hash() != base_hash),
            "every mutable semantic event field must move the typed identity"
        );
        let mut envelope = base;
        envelope.wall_ns = Some(u64::MAX);
        assert_eq!(
            envelope.content_hash(),
            base_hash,
            "the declared wall-clock exclusion must not move identity"
        );
    }

    #[test]
    // The assertions independently perturb every framing/version component and
    // therefore intentionally remain a single identity-separation proof.
    #[allow(clippy::too_many_lines)]
    fn event_content_identity_domain_version_and_wire_schema_bytes_are_independent() {
        const FRAME_HEADER_LEN: usize = ident::REPLAY_IDENTITY_DOMAIN.len()
            + core::mem::size_of::<u32>()
            + core::mem::size_of::<u64>();

        let event = event_with_kind(EventKind::Cancellation {
            reason: "budget".into(),
        });
        let current =
            event.content_identity_with_versions(EVENT_CONTENT_IDENTITY_VERSION, SCHEMA_VERSION);
        let domain_mutation = event.content_identity_with_schema(
            "org.frankensim.fs-obs.event-content.v2.alternate",
            EVENT_CONTENT_IDENTITY_VERSION,
            SCHEMA_VERSION,
        );
        let identity_version_mutation = event
            .content_identity_with_versions(EVENT_CONTENT_IDENTITY_VERSION + 1, SCHEMA_VERSION);
        let wire_schema_mutation = event
            .content_identity_with_versions(EVENT_CONTENT_IDENTITY_VERSION, SCHEMA_VERSION + 1);

        let (identity_tag, current_identity_version) =
            identity_field(current.canonical_bytes(), "event_content_identity_version");
        let (schema_tag, current_wire_schema) =
            identity_field(current.canonical_bytes(), "event_wire_schema_version");
        assert_eq!(identity_tag, 0x02);
        assert_eq!(schema_tag, 0x02);
        // Canonical bytes carry the shared replay-identity framing —
        // "fsid" + u32 ident-schema version + u64 kind length — and bind the
        // event artifact domain as the length-framed identity kind, so the
        // domain begins after the fixed 16-byte frame header rather than at
        // byte zero (length framing is what makes prefix collisions
        // unrepresentable; see ident_003).
        assert!(
            current
                .canonical_bytes()
                .starts_with(ident::REPLAY_IDENTITY_DOMAIN.as_bytes())
        );
        assert!(
            current.canonical_bytes()[FRAME_HEADER_LEN..]
                .starts_with(EVENT_CONTENT_IDENTITY_DOMAIN.as_bytes())
        );
        assert!(
            domain_mutation.canonical_bytes()[FRAME_HEADER_LEN..]
                .starts_with(b"org.frankensim.fs-obs.event-content.v2.alternate")
        );
        assert_eq!(
            current_identity_version,
            u64::from(EVENT_CONTENT_IDENTITY_VERSION)
                .to_le_bytes()
                .as_slice()
        );
        assert_eq!(
            current_wire_schema,
            u64::from(SCHEMA_VERSION).to_le_bytes().as_slice()
        );
        assert_eq!(
            identity_field(
                domain_mutation.canonical_bytes(),
                "event_content_identity_version"
            )
            .1,
            current_identity_version,
            "artifact-domain mutation must leave identity-version bytes unchanged"
        );
        assert_eq!(
            identity_field(
                domain_mutation.canonical_bytes(),
                "event_wire_schema_version"
            )
            .1,
            current_wire_schema,
            "artifact-domain mutation must leave wire-schema bytes unchanged"
        );

        assert_eq!(
            identity_field(
                identity_version_mutation.canonical_bytes(),
                "event_wire_schema_version"
            )
            .1,
            current_wire_schema,
            "identity-version mutation must leave the wire-schema bytes unchanged"
        );
        assert_eq!(
            identity_field(
                wire_schema_mutation.canonical_bytes(),
                "event_content_identity_version"
            )
            .1,
            current_identity_version,
            "wire-schema mutation must leave the identity-version bytes unchanged"
        );
        assert_ne!(
            identity_field(
                identity_version_mutation.canonical_bytes(),
                "event_content_identity_version"
            )
            .1,
            current_identity_version
        );
        assert_ne!(
            identity_field(
                wire_schema_mutation.canonical_bytes(),
                "event_wire_schema_version"
            )
            .1,
            current_wire_schema
        );
        assert_ne!(
            current.canonical_bytes(),
            identity_version_mutation.canonical_bytes()
        );
        assert_ne!(current.canonical_bytes(), domain_mutation.canonical_bytes());
        assert_ne!(current.root(), domain_mutation.root());
        assert_ne!(current.root(), identity_version_mutation.root());
        assert_ne!(
            current.canonical_bytes(),
            wire_schema_mutation.canonical_bytes()
        );
        assert_ne!(current.root(), wire_schema_mutation.root());
    }

    #[test]
    fn retained_event_identity_receipts_admit_exactly_or_fail_closed() {
        let event = event_with_kind(EventKind::StormAssertion {
            name: "no-leak".into(),
            pass: true,
            seed: 23,
        });
        let captured = event.content_identity_receipt();
        assert_eq!(
            captured.declared_identity_version(),
            EVENT_CONTENT_IDENTITY_VERSION
        );
        assert_eq!(captured.root(), fnv1a64(captured.canonical_bytes()));

        let retained = EventIdentityReceipt::from_retained_parts(
            captured.declared_identity_version(),
            captured.canonical_bytes().to_vec(),
            captured.root(),
        );
        assert_eq!(event.admit_content_identity(&retained), Ok(()));

        let stale = EventIdentityReceipt::from_retained_parts(
            0,
            captured.canonical_bytes().to_vec(),
            captured.root(),
        );
        assert_eq!(
            event.admit_content_identity(&stale),
            Err(EventIdentityAdmissionError::UnsupportedVersion(
                EventIdentityVersionError {
                    declared: 0,
                    supported: EVENT_CONTENT_IDENTITY_VERSION,
                }
            ))
        );

        let wrong_root = EventIdentityReceipt::from_retained_parts(
            EVENT_CONTENT_IDENTITY_VERSION,
            captured.canonical_bytes().to_vec(),
            captured.root() ^ 1,
        );
        assert!(matches!(
            event.admit_content_identity(&wrong_root),
            Err(EventIdentityAdmissionError::RootMismatch { .. })
        ));

        let mut foreign_bytes = captured.canonical_bytes().to_vec();
        let last = foreign_bytes
            .last_mut()
            .expect("event identity canonical bytes are non-empty");
        *last ^= 1;
        let foreign_root = fnv1a64(&foreign_bytes);
        let self_consistent_but_foreign = EventIdentityReceipt::from_retained_parts(
            EVENT_CONTENT_IDENTITY_VERSION,
            foreign_bytes,
            foreign_root,
        );
        assert_eq!(
            event.admit_content_identity(&self_consistent_but_foreign),
            Err(EventIdentityAdmissionError::CanonicalBytesMismatch {
                declared_root: foreign_root,
                expected_root: captured.root(),
            })
        );
    }

    #[test]
    fn event_content_identity_versions_fail_closed() {
        assert!(check_event_content_identity_version(EVENT_CONTENT_IDENTITY_VERSION).is_ok());
        for declared in [0, EVENT_CONTENT_IDENTITY_VERSION + 1] {
            assert_eq!(
                check_event_content_identity_version(declared),
                Err(EventIdentityVersionError {
                    declared,
                    supported: EVENT_CONTENT_IDENTITY_VERSION,
                })
            );
        }
    }

    #[test]
    fn sequence_numbers_are_monotone_per_emitter() {
        let events = sample_events();
        for (i, e) in events.iter().enumerate() {
            assert_eq!(e.seq, i as u64);
        }
    }

    #[test]
    fn sampling_cadence_is_zero_based_and_total_at_u64_max() {
        let every_one = SamplingCadence::every(NonZeroU64::MIN);
        let every_two =
            SamplingCadence::every(NonZeroU64::new(2).expect("two is a non-zero test cadence"));

        assert!(!SamplingCadence::never().admits(0));
        assert!(!SamplingCadence::never().admits(u64::MAX));
        assert!(every_one.admits(0));
        assert!(every_one.admits(u64::MAX));
        assert!(every_two.admits(0));
        assert!(!every_two.admits(1));
        assert!(every_two.admits(2));
        assert!(!every_two.admits(u64::MAX));
    }

    #[test]
    fn emission_gate_samples_levels_independently() {
        let gate = EmissionGate::new(
            SamplingCadence::every(NonZeroU64::new(2).expect("two is a non-zero trace cadence")),
            SamplingCadence::every(NonZeroU64::new(3).expect("three is a non-zero info cadence")),
        );

        assert!(gate.admits(Severity::Trace, 0));
        assert!(!gate.admits(Severity::Trace, 1));
        assert!(gate.admits(Severity::Trace, 2));
        assert!(!gate.admits(Severity::Info, 2));
        assert!(gate.admits(Severity::Info, 3));

        let quiet = EmissionGate::new(SamplingCadence::never(), SamplingCadence::never());
        assert!(!quiet.admits(Severity::Trace, u64::MAX));
        assert!(!quiet.admits(Severity::Info, u64::MAX));
        assert!(quiet.admits(Severity::Warn, u64::MAX));
        assert!(quiet.admits(Severity::Error, u64::MAX));
    }

    #[test]
    fn emit_gated_is_lazy_and_sequences_count_only_emissions() {
        use core::cell::Cell;

        let quiet = EmissionGate::new(SamplingCadence::never(), SamplingCadence::never());
        let trace_even = EmissionGate::new(
            SamplingCadence::every(NonZeroU64::new(2).expect("two is a non-zero trace cadence")),
            SamplingCadence::never(),
        );
        let builds = Cell::new(0_u64);
        let mut emitter = Emitter::new("sampling", "gate");

        let rejected = emitter.emit_gated(quiet, Severity::Trace, u64::MAX, || {
            panic!("a rejected opportunity must not build its payload or caller clock")
        });
        assert!(rejected.is_none());
        assert_eq!(builds.get(), 0);

        let trace = emitter
            .emit_gated(trace_even, Severity::Trace, 2, || {
                builds.set(builds.get() + 1);
                (
                    EventKind::TileComplete {
                        tile: 2,
                        kernel: "sampled".into(),
                    },
                    Some(20),
                )
            })
            .expect("even trace opportunity is admitted");
        let info = emitter.emit_gated(trace_even, Severity::Info, 0, || {
            panic!("the disabled info cadence must stay lazy")
        });
        let warn = emitter
            .emit_gated(quiet, Severity::Warn, u64::MAX, || {
                builds.set(builds.get() + 1);
                (
                    EventKind::Cancellation {
                        reason: "budget".into(),
                    },
                    Some(30),
                )
            })
            .expect("warn events bypass sampling");
        let direct = emitter.emit(
            Severity::Info,
            EventKind::Custom {
                name: "direct".into(),
                json: "{}".into(),
            },
            None,
        );

        assert!(info.is_none());
        assert_eq!(builds.get(), 2);
        assert_eq!((trace.seq, warn.seq, direct.seq), (0, 1, 2));
        assert_eq!((trace.wall_ns, warn.wall_ns), (Some(20), Some(30)));
    }

    #[test]
    fn golden_line_shape_is_stable() {
        // Schema evolution is additive-only; this golden line is the contract.
        // Changing it requires a SCHEMA_VERSION bump and a semantic justification
        // (golden-evidence policy, AGENTS.md).
        let e = Event {
            session: "s".into(),
            scope: "a/b".into(),
            seq: 5,
            severity: Severity::Info,
            kind: EventKind::GradientCheck {
                op: "poisson".into(),
                max_rel_err: 1e-9,
                pass: true,
            },
            wall_ns: None,
        };
        // Note: Rust's shortest-round-trip float Display never uses scientific
        // notation, so 1e-9 serializes as 0.000000001 — that IS the contract.
        assert_eq!(
            e.to_jsonl(),
            "{\"v\":1,\"session\":\"s\",\"scope\":\"a/b\",\"seq\":5,\"severity\":\"info\",\
             \"kind\":\"gradient_check\",\"payload\":{\"op\":\"poisson\",\
             \"max_rel_err\":0.000000001,\"pass\":true}}"
        );
    }

    #[test]
    fn validator_rejects_corruption() {
        let good = sample_events()[0].to_jsonl();
        for bad in [
            String::new(),
            "not json".to_string(),
            good.replace("\"v\":1", "\"v\":99"),
            good.replace("solver_residual", "mystery_kind"),
            good.replace("\"session\"", "\"sesion\""),
            good[..good.len() - 1].to_string(),
        ] {
            assert!(validate_line(&bad).is_err(), "should reject: {bad}");
        }
    }

    #[test]
    fn failure_lint_demands_reproduction_ingredients() {
        let mut em = Emitter::new("s", "x");
        let bad = em.emit(
            Severity::Error,
            EventKind::ConformanceCase {
                suite: "s".into(),
                case: "c".into(),
                pass: false,
                detail: String::new(),
                seed: 0,
            },
            None,
        );
        assert!(lint_failure_record(&bad).is_err());
        let good = em.emit(
            Severity::Error,
            EventKind::StormAssertion {
                name: "no-leak".into(),
                pass: false,
                seed: 42,
            },
            None,
        );
        assert!(lint_failure_record(&good).is_ok());
        let bad_storm = em.emit(
            Severity::Error,
            EventKind::StormAssertion {
                name: "no-leak".into(),
                pass: false,
                seed: 0,
            },
            None,
        );
        assert!(lint_failure_record(&bad_storm).is_err());
    }

    #[test]
    fn escaping_handles_hostile_strings() {
        let mut em = Emitter::new("s\"es\\sion\n", "sc\tope");
        let e = em.emit(
            Severity::Error,
            EventKind::Cancellation {
                reason: "quote\" backslash\\ newline\n tab\t".into(),
            },
            None,
        );
        let line = e.to_jsonl();
        validate_line(&line).unwrap_or_else(|err| panic!("{line}: {err}"));
        assert!(!line.contains('\n'), "JSONL lines must be single-line");
    }

    #[test]
    fn non_finite_floats_are_tagged_not_invalid() {
        let mut em = Emitter::new("s", "x");
        let e = em.emit(
            Severity::Info,
            EventKind::SolverResidual {
                solver: "cg".into(),
                iter: 1,
                residual: f64::NAN,
            },
            None,
        );
        let line = e.to_jsonl();
        validate_line(&line).expect("tagged non-finite must stay valid");
        assert!(line.contains("non-finite:NaN"));
    }

    #[test]
    fn fnv_matches_known_answers() {
        // Published FNV-1a 64 test vectors.
        assert_eq!(fnv1a64(b""), 0xcbf2_9ce4_8422_2325);
        assert_eq!(fnv1a64(b"a"), 0xaf63_dc4c_8601_ec8c);
        assert_eq!(fnv1a64(b"foobar"), 0x85944171f73967e8);
    }
}
