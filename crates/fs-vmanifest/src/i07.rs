//! The I07 (AP242 semantic manufacturing round trip) VerificationManifest
//! draft (bead frankensim-leapfrog-2026-program-i94v.2.2.8.1).
//!
//! The baseline lattice ([S]) freezes a bounded ISO 10303-21/AP242 profile,
//! product/configuration/occurrence lineage, representation-qualified geometry,
//! semantic PMI/GD&T, material/process/lot references, harness/EWIS graphs,
//! canonical deterministic export, and an exhaustive semantic-difference
//! receipt. The maximal lattice ([F]/[M]) adds safe opaque-extension capsules,
//! kinematic-pair semantics, partial bidirectional edition migrations,
//! governed real-corpus/standards/interoperability validation with a narrowly
//! scoped external-artifact waiver, proof-carrying semantic equivalence and
//! sheaf/stack descent, plus an independent counterexample lane. A weaker
//! receipt closes only its own lattice element.
//!
//! This is preregistration, not standards, manufacturing, drawing-approval, or
//! regulatory authority. In particular, visual similarity, parser acceptance,
//! and byte-stable output are never substitutes for semantic preservation.

use crate::{
    Ambition, CampaignTier, ClaimPolarity, ClaimSpec, FiveExplicits, FixturePin, FixtureSource,
    GauntletTier, ManifestDraft, ObligationRow, OracleRoute, Partition, ToleranceSemantics, Waiver,
};

const CAMPAIGN_POLICY_FIXTURE: &str = "i07-campaign-policy-v1";
const PROFILE_REGISTRY_FIXTURE: &str = "i07-profile-registry-v1";
const REAL_CORPUS_SLOT: &str = "i07-governed-real-ap242-corpus";

/// Build the I07 draft. Consumers freeze it themselves; campaign-time proof
/// and authority remain outside this constructor.
#[must_use]
pub fn i07_draft() -> ManifestDraft {
    ManifestDraft {
        initiative: "I07",
        title: "AP242 semantic manufacturing round-trip gate: bounded Part-21 profiles, \
                configuration lineage, geometry, PMI, materials, EWIS, kinematics, \
                semantic-loss receipts, governed real-corpus validation, edition migration, \
                and proof-carrying sheaf/stack descent",
        version: 1,
        explicits: FiveExplicits {
            units: "six-base quantity system L M T I Theta N; every AP242 numeric value \
                    retains source unit, semantic quantity kind, uncertainty/qualifier, \
                    coordinate frame, and exact conversion path; angles declare radians \
                    or degrees as semantic angle kinds; entity identities, exact \
                    receipt/reconciliation outcomes, and admission verdicts use unit 'bit'; only \
                    explicitly normalized numeric geometry/error defect scores use unit '1'",
            seeds: "Philox 4x32-10 aliases are exact, case-sensitive UTF-8 'i07/<stream>'; \
                    the ten-round Random123 word convention uses M0=0xD2511F53, \
                    M1=0xCD9E8D57, W0=0x9E3779B9, and W1=0xBB67AE85; \
                    d=BLAKE3::derive_key('org.frankensim.i07.fixture-stream.v1', alias_utf8); \
                    k0=LE32(d[0..4]); k1=LE32(d[4..8]); c0=low32(case_index); \
                    c1=high32(case_index); c2=low32(output_block_ordinal); \
                    c3=high32(output_block_ordinal). Lane is the Philox output-word index \
                    0..3 and is never folded into the counter; each output word is serialized \
                    LE32, so byte offset n selects block floor(n/16), lane floor((n mod 16)/4), \
                    and byte n mod 4. Native-endian casts and Unicode normalization are forbidden. \
                    development indices 0..=4095, core stage-heldout 65536..=69631, \
                    maximal stage-heldout 131072..=135167. These public deterministic \
                    ranges provide replay and stage separation only, never IID, secrecy, \
                    untouched-population, or statistical validation authority. KAT_GATE: before \
                    evidence generation the concrete successor freezes independently reproduced \
                    development and heldout-format alias/digest/key/counter/four-output-word \
                    vectors from two implementations. The set covers every counter word and output \
                    lane, exact LE serialized bytes, nonzero high32 case and block words, block \
                    0-to-1 and byte-offset 15-to-16 boundaries, and case-sensitive/UTF-8 alias \
                    twins. Cross-alias derived-key/counter collisions fail admission; version-1 \
                    prose does not invent those computed words",
            budgets: "smoke <= 120 s on one host and <= 4 GiB; core <= 60 min and \
                      <= 24 GiB; max <= 16 h on a quiet campaign host and <= 64 GiB; \
                      parser/entity/reference caps, text/blob caps, graph depth, SCC, \
                      geometry subdivision, proof, and migration-step budgets are fixed \
                      by fixture. Core polls occur at entity/reference/semantic-atom batches \
                      <=4096 items or geometry tiles <=16384 primitives and demonstrate \
                      request observation <=250 ms on the admitted host. Max polls also bound \
                      opaque scans and migration steps <=4096 items, proof declaration/checker \
                      steps <=256 declarations, falsifier batches <=1024 candidates, and governed-\
                      corpus work <=256 items or <=4096 semantic atoms per cancellable batch with \
                      request observation <=1 s. Every external reader/vendor/checker process must \
                      heartbeat and support supervised terminate-drain-finalize within the same bound \
                      or preflight refuses; no opaque tool call is an unbounded tile. Accuracy is fixed by claim \
                      tolerances",
            versions: "fs-vmanifest schema v2; I07 profile-registry/profile-definition/\
                       required-assignment-universe schemas v1; governed receipt-schema-set and \
                       authority-head transition schemas v1; \
                       semantic-atom/loss-receipt schema v1; Machine-IR identity/lineage \
                       schema v1; material/process/lot/reference envelope schema v1 with \
                       every referenced revision/content root pinned per case; AP242 \
                       edition/profile identities are opaque pinned \
                       artifact digests rather than guessed prose; toolchain pinned by \
                       rust-toolchain.toml and sibling revisions by constellation.lock",
            capabilities: "baseline default-off frontier-ap242-semantic-roundtrip: bounded clear-text \
                           ISO-10303-21 parsing, profile admission, product/assembly/\
                           configuration lineage, exact-and-tessellated geometry, \
                           semantic-PMI, material-process-lot, harness-EWIS, canonical \
                           export, semantic-loss receipt; maximal default-off \
                           moonshot-ap242-proof-carrying-semantics: safe opaque capsules, \
                           kinematic pairs, edition lenses, proof-carrying semantic \
                           equivalence and descent; no network, FFI, foreign code, or \
                           executable AP242 payload in the production import/export path. \
                           Development-only proof infrastructure may execute only a pinned, \
                           sandboxed, capability-limited formal checker or ingest its \
                           independently authenticated receipt; no licensed-text embedding, \
                           drawing approval, or regulatory authority; deterministic mode \
                           mandatory for G5",
        },
        claims: i07_claims(),
        fixtures: i07_fixtures(),
        obligations: i07_obligations(),
        waivers: i07_waivers(),
        amendment_rules: "After campaign admission every semantic change is a successor \
                          through FrozenManifest::amend. Changes to edition/profile/\
                          conformance class, schema digest, supported entity or atom set, \
                          unit/frame rule, identity/lineage rule, fixture partition, \
                          geometry norm/band, loss taxonomy, oracle/threat model, \
                          material/lot/evidence pin, capability, or terminal-state policy \
                          invalidate exactly their reverse-dependency descendants. No \
                          result, vendor exception, waiver, or legal interpretation may \
                          edit this version in place",
    }
}

#[allow(clippy::too_many_lines)]
fn i07_claims() -> Vec<ClaimSpec> {
    vec![
        ClaimSpec {
            id: "i07-part21-profile-bounded-parse",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "For an exact admitted ProfileKeyDigest/ProfileEntryRoot pair, the bounded clear-text \
                        ISO 10303-21 reader reconstructs the complete token, section, \
                        instance, complex-entity, attribute, aggregate, select, and \
                        reference graph prescribed by the pinned supported subset, \
                        retaining source spans and refusing malformed, ambiguous, \
                        over-budget, wrong-edition, wrong-schema, or unsupported input \
                        before any engineering semantic authority is emitted",
            hypotheses: &[
                "ProfileKeyDigest binds the Part-21 encoding edition, AP242 schema artifact digest, application-protocol edition, conformance class/profile, implementation method, character encoding, signature policy, external-reference policy, validation-rule inventory, semantic-atom inventory, unit/frame policy, canonical-writer policy, no-claim policy, coverage floors, and resource limits; ProfileEntryRoot additionally binds ProfileDefinitionSchemaRoot and the canonical definition bytes and is the complete expanded authority identity",
                "the supported entity/type/attribute/WHERE/UNIQUE rule inventory and all implementation-defined limits are content-pinned before parsing",
                "forward references and legal strongly connected reference components are retained as graph structure; dangling, duplicate, type-invalid, or profile-forbidden references refuse, while cycles are never rejected merely because they are cycles",
                "all strings, binary/hex encodings, comments, omitted/derived markers, enumeration/select tags, numeric lexical forms, and source spans stay bounded and loss-accounted",
            ],
            qoi: "exact_profile_token_entity_reference_and_refusal_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G0,
            oracle: OracleRoute {
                identity: "i07.oracle.part21.v1 at fs-vmanifest-oracles/i07/part21.rs::parse_generated_ast_and_refusals",
                independent: true,
                tcb_overlap: "shares only fixture bytes and profile-registry bytes; the generator AST, independent reader, and resource-bound auditor do not use production lexer/parser code",
            },
            activation: "a pre-candidate successor has discharged PROFILE_INSTANCE_GATE and KAT_GATE and binds the exact required ProfileDefinitionSchemaRoot, ProfileKeyDigest/ProfileEntryRoot, registry/required-assignment-universe/assignment roots, and claim/leaf set; the bounded reader and typed refusal surface then exist for that registry",
            kill: "one accepted wrong profile, token/attribute/reference mismatch, silent numeric/string normalization, legal-cycle rejection, dangling-reference admission, unbounded allocation, or partial semantic publication after refusal kills the parser claim",
            fallback: "quarantine the entire exchange as Unsupported with exact profile/rule/span/budget diagnostics and a content digest",
            no_claim: "parser agreement proves neither AP242 semantic conformance nor engineering correctness; no arbitrary STEP application protocol, vendor dialect, XML encoding, signature, or external fetch is implied",
        },
        ClaimSpec {
            id: "i07-product-assembly-configuration-lineage",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Within the frozen product-structure subset, import-export-\
                        reimport preserves products, versions, definitions/views, \
                        occurrences, assembly usage, configuration/effectivity, external \
                        identities, provenance, sharing, and instance transforms as a \
                        typed attributed graph, independent of exchange-file #numbering",
            hypotheses: &[
                "source-local Part-21 instance numbers are provenance locators, never global semantic identities",
                "stable BodyId and occurrence/configuration identities derive from pinned semantic keys plus explicit lineage and collision/refusal rules, not traversal order or geometry proximity",
                "assembly sharing, multiplicity, parent-child direction, transform composition order, frame handedness, units, configuration alternatives, and effectivity domains are explicit",
                "split, merge, suppression, substitution, and ambiguous rebinding are represented by typed lineage morphisms or refusal; no heuristic topological naming is silently promoted",
            ],
            qoi: "exact_product_occurrence_configuration_lineage_graph_isomorphism_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i07.oracle.product_graph.v1 at fs-vmanifest-oracles/i07/product_graph.rs::check_typed_graph_isomorphism",
                independent: true,
                tcb_overlap: "shares canonical semantic-atom vocabulary and input fixture bytes; independent graph construction and isomorphism do not consume production Machine IR ids",
            },
            activation: "bounded Part-21/profile parsing is green and Machine-IR lineage schema v1 is admitted",
            kill: "one merged distinct occurrence, duplicated shared definition, lost configuration/effectivity edge, frame/unit inversion, lineage collision, traversal-order identity, or ambiguous rebind admitted as preserved kills promotion",
            fallback: "preserve the parse graph and emit Ambiguous/Unsupported semantic-loss atoms; refuse stable native identity for the affected component",
            no_claim: "no arbitrary vendor feature-history recovery, CAD topological naming through undeclared edits, bill-of-material approval, or configuration-management authority",
        },
        ClaimSpec {
            id: "i07-representation-qualified-geometry-roundtrip",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Exact B-rep/NURBS and tessellated AP242 representations remain \
                        separately typed. Every admitted conversion preserves units, \
                        frames, orientation, representation role, validation properties, \
                        and topology predicates it actually certifies. Acceptance is the \
                        conjunction of independently enclosed bidirectional geometric/integral \
                        defects inside preregistered normalized bands and every applicable \
                        discrete orientation, closure, manifoldness, self-intersection, \
                        validation-property, and topology certificate being Satisfied; Unknown \
                        or unsupported certificate state blocks that property's promotion",
            hypotheses: &[
                "source and target representation classes, trims, seams, shell/solid role, tolerances, validation properties, coordinate systems, and transforms are frozen",
                "exact-to-exact paths retain analytic/control data unless a named certified refit is explicitly selected; tessellated input is never relabeled exact and an exact-to-tessellated path is a BoundedApproximation loss disposition",
                "directed Hausdorff, boundary/volume measure, moment, normal-orientation, closure/manifold, self-intersection, and topology checks use property-specific complete oracles; sampled distance or visual overlays are falsifiers only",
                "a topology-preservation claim additionally binds reach/feature-separation and coverage premises sufficient for the exact property; absent premises yield Unknown, not inference from a small distance",
            ],
            qoi: "maximum_preregistered_normalized_bidirectional_geometry_validation_and_topology_defect",
            unit: "1",
            tolerance: ToleranceSemantics::Interval { lo: 0.0, hi: 1.0 },
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "i07.oracle.geometry.v1 composite fs-vmanifest-oracles/i07/geometry.rs::enclose_roundtrip and independent interval subdivision/oriented-intersection/winding checkers",
                independent: true,
                tcb_overlap: "shares exact source geometry bytes and unit/frame schema only; independent analytic, interval, and exact-predicate paths do not use production conversion meshes or tolerances",
            },
            activation: "product/occurrence frames are green and at least one certified exact or tessellated chart route exists",
            kill: "representation laundering, one out-of-band directed defect, orientation/frame/unit error, validation-property fabrication, missed seam/self-intersection/topology counterexample, or sampled-only positive certificate kills the affected route",
            fallback: "retain source representation, emit BoundedApproximation or Unknown with every measured/certified defect, and refuse downstream geometry authority above that color",
            no_claim: "no universal B-rep healing, continuum equivalence from tessellation, arbitrary trimmed-surface watertightness, feature-history equivalence, or manufacturing fitness",
        },
        ClaimSpec {
            id: "i07-semantic-pmi-gdt-datum-texture",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "For the pinned semantic-PMI subset, dimensions, tolerances, \
                        datum features/systems/targets and precedence, geometric \
                        characteristics/modifiers/zones, surface texture, annotations, \
                        requirements, and their associations reimport to the same typed \
                        semantic graph and target identities with exact unit/frame and \
                        value semantics",
            hypotheses: &[
                "semantic PMI and presentation-only graphical PMI are disjoint types; presentation geometry may be retained but never supplies missing semantic authority",
                "every numeric value retains nominal/limit/plus-minus semantics, semantic quantity kind, source unit, conversion, decimal/rational intent where available, uncertainty/qualifier, and applicable frame",
                "datum precedence, target/feature association, tolerance-zone geometry, material-condition/modifier semantics, and requirement scope are explicit in the supported profile",
                "an unresolved target, unsupported modifier, conflicting duplicate, or presentation-only annotation receives a named loss/refusal disposition",
            ],
            qoi: "exact_semantic_pmi_datum_texture_and_association_graph_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i07.oracle.pmi.v1 at fs-vmanifest-oracles/i07/pmi.rs::check_semantic_atoms_and_associations",
                independent: true,
                tcb_overlap: "shares the supported semantic-atom vocabulary; generator truth graph and dimensional checker are independent of production PMI mapping",
            },
            activation: "geometry/feature target identities and quantity-kind schema are admitted",
            kill: "one value/modifier/datum-order/unit/frame/target change, semantic/presentation laundering, silently dropped unsupported atom, or ambiguous target reported Preserved kills the PMI route",
            fallback: "retain graphical/source payload and emit itemized Unsupported or Ambiguous atoms; require human engineering review",
            no_claim: "no drawing approval, inspection-plan adequacy, legal interpretation of GD&T, tolerance-stack validation, manufacturability, or meaning reconstructed from pixels/text alone",
        },
        ClaimSpec {
            id: "i07-material-process-lot-requirement-evidence",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Named material/specification, process/heat-treatment/coating, \
                        lot/batch, requirement, document, source/license, evidence, and \
                        supersession references preserve exact identity, scope, \
                        applicability, lineage, and association to product occurrences \
                        without inventing inaccessible property values or validation",
            hypotheses: &[
                "every material/process/lot/evidence reference has a typed identity namespace, edition/revision when applicable, content digest or explicit unpinned state, license/access class, and occurrence/feature scope",
                "external documents and property cards are data references only; the no-network reader never dereferences a URI and unavailable content remains unavailable",
                "material name, grade, temper, process state, lot, specimen/population, validity domain, covariance, and evidence color remain distinct fields",
                "conflicting references coexist with a conflict receipt; selection/fusion requires a separately versioned policy and cannot arise from map overwrite",
            ],
            qoi: "exact_material_process_lot_requirement_evidence_lineage_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i07.oracle.material_lineage.v1 at fs-vmanifest-oracles/i07/material_lineage.rs::check_reference_graph",
                independent: true,
                tcb_overlap: "shares identity namespace and canonical fixture bytes only; independent reference-graph checker does not query production material selection",
            },
            activation: "product occurrence identities and fs-matdb reference envelope are available",
            kill: "one lost lot/process/revision/license/evidence link, silent URI fetch, grade-only property inference, conflict overwrite, or scope rebind kills the affected material lineage claim",
            fallback: "preserve opaque/reference bytes and mark the property/material decision Unsupported or Unknown pending an authenticated card",
            no_claim: "reference preservation does not validate material truth, supplier conformity, lot acceptance, process completion, property applicability, service life, or regulatory compliance",
        },
        ClaimSpec {
            id: "i07-harness-ewis-connectivity-semantics",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "Within the pinned harness/EWIS subset, connectors, cavities/\
                        pins, contacts, wires/cores, splices, shields, bundles, branches, \
                        route segments, terminals, grounds/bonds, labels, gauges, \
                        materials, occurrence placement, and requirement/evidence links \
                        round-trip to an isomorphic typed HarnessGraph with no inferred \
                        connectivity",
            hypotheses: &[
                "electrical connectivity, physical containment/bundling, geometric routing, shielding, grounding/bonding, and logical signal naming are separate edge families",
                "pin/cavity/contact direction, conductor/core identity, splice multiplicity, shield termination, wire gauge/unit/system, occurrence context, and route frame are explicit",
                "open circuits, intentionally unconnected pins, multi-drop nets, duplicate labels, connector reuse, and shared shields are represented without name-based merging",
                "unsupported EWIS constructs, missing endpoints, or ambiguous ground/bond ownership refuse or receive semantic-loss atoms",
            ],
            qoi: "exact_harness_connectivity_containment_route_and_identity_graph_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i07.oracle.harness.v1 at fs-vmanifest-oracles/i07/harness.rs::check_multigraph_isomorphism",
                independent: true,
                tcb_overlap: "shares typed HarnessGraph vocabulary only; generator netlist/containment/route graphs and independent isomorphism checker do not consume production mappings",
            },
            activation: "product occurrence identity, frames, and harness schema v1 are admitted",
            kill: "one false/missing connection, label-based merge, lost shield/ground/bond termination, route-frame/unit error, containment-connectivity conflation, or ambiguous endpoint called Preserved kills the harness route",
            fallback: "retain source entities plus an itemized Unknown/Unsupported HarnessGraph fragment and block electrical decisions that depend on it",
            no_claim: "no ampacity, voltage-drop, EMC, creepage/clearance, installation, maintainability, airworthiness, or EWIS regulatory compliance claim",
        },
        ClaimSpec {
            id: "i07-exhaustive-semantic-loss-receipt",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "For every admitted input, the semantic-difference receipt \
                        reconciles the complete bounded source and target atom/edge \
                        inventories: each source atom has exactly one DispositionRecord, each \
                        target atom has exactly one OriginRecord, each record carries a typed \
                        zero/one/many counterpart set plus its derivation, and every PreservedExact, \
                        PreservedNormalized, SemanticallyMapped, BoundedApproximation, \
                        OpaquePreservedNoSemanticAuthority, Ambiguous, UnsupportedWithRefusal, \
                        or DroppedWithRefusal item names subject, \
                        rule, reason, bound, owner, authority impact, and provenance",
            hypotheses: &[
                "the ProfileKey freezes a total semantic-atom and relationship inventory for every supported entity plus a catch-all token/reference inventory for unsupported entities",
                "receipt records are disjoint and exhaustive over stable atom ids; a source DispositionRecord may name zero, one, or many target atoms and a target OriginRecord may name zero, one, or many source atoms, so split, merge, normalization, synthesis, and deletion remain representable without duplicate authority",
                "forward and inverse source-target incidence relations are exact transposes; synthesized targets with zero source atoms name a derivation rule, dropped sources with zero targets name a refusal, and every nonempty many-to-many relation names its semantic mapping rule",
                "reconciliation includes source/target counts and roots, record roots, orphan/dangling checks, approximation/error ledgers, opaque capsules, warnings, refusals, and all semantic domains in this manifest",
                "Failed execution, Refuted claim, Unknown meaning, Unsupported construct, PartialEvidence, IntegrityFailed, and InfrastructureFailed are orthogonal terminal states and never collapse into an absent or falsely clean empty loss receipt; legitimately empty admitted inventories retain both roots and zero reconciliation counts",
            ],
            qoi: "exact_source_target_atom_reconciliation_and_loss_classification_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i07.oracle.semantic_loss.v1 at fs-vmanifest-oracles/i07/semantic_loss.rs::reconcile_independent_inventory",
                independent: true,
                tcb_overlap: "shares profile-registry and semantic-atom schema bytes; independent source/target walkers, coverage auditor, and receipt reconciler do not consume production loss rows",
            },
            activation: "all baseline semantic atomizers exist; this receipt is mandatory even when every other baseline claim refuses",
            kill: "one unclassified source/target atom, duplicate or missing DispositionRecord/OriginRecord, non-transpose incidence relation, unowned synthesis/split/merge/drop, missing loss owner/bound/authority effect, empty-success laundering, silent unsupported drop, or receipt produced from the same production bookkeeping without independent reconciliation kills baseline promotion",
            fallback: "refuse the round trip and retain parse graph, partial artifacts, and a FailureBundle that explicitly says receipt completeness is Unknown",
            no_claim: "completeness is only over the exact pinned atom inventory/profile and bounded parse; a clean receipt is not arbitrary vendor-dialect, feature-history, legal, drawing, or manufacturing equivalence",
        },
        ClaimSpec {
            id: "i07-canonical-export-determinism-and-containment",
            ambition: Ambition::Solid,
            polarity: ClaimPolarity::Affirmative,
            statement: "The canonical writer emits one bounded profile-admitted \
                        clear-text exchange whose ordering, numbering, whitespace, numeric \
                        rendering, string encoding, and content root are deterministic \
                        across schedules on the same admitted ISA/toolchain profile; it \
                        performs no network access or foreign execution, publishes only \
                        after complete validation, and reimports to the frozen semantic \
                        receipt",
            hypotheses: &[
                "canonical ordering and tie-breaks operate on semantic identity plus complete content digests, never hash-map or traversal arrival order",
                "numeric rendering is round-trip safe under the pinned lexical policy and preserves semantic quantity/unit records even when presentation normalizes",
                "external identifiers, URIs, signatures, and documents are inert data; no resolver, plugin, macro, script, or foreign code executes",
                "request-cancel-drain-finalize uses a transactional temporary artifact and atomically publishes only a fully reimported/adjudicable file and receipt",
            ],
            qoi: "exact_canonical_bytes_reimport_and_no_partial_publication_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G5,
            oracle: OracleRoute {
                identity: "i07.oracle.canonical_writer.v1 at fs-vmanifest-oracles/i07/canonical_writer.rs::independent_reparse_and_compare",
                independent: true,
                tcb_overlap: "shares deterministic scalar formatting specification and BLAKE3; independent parser and semantic-root builder do not share production writer traversal",
            },
            activation: "baseline semantic mappings and transactional artifact sink exist",
            kill: "same-profile schedule drift, non-round-tripping numeric/string output, unresolved reference, side effect/network access, partial publication after fault/cancel, or canonical bytes that reimport to a different semantic receipt kills canonical export",
            fallback: "retain the native content-addressed package and semantic-loss receipt; emit no AP242 file",
            no_claim: "canonical bytes are a FrankenSim profile convention, not vendor byte identity, digital-signature continuity, compression equivalence, or semantic proof by themselves",
        },
        ClaimSpec {
            id: "i07-safe-opaque-extension-capsules",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "A syntactically valid but semantically unsupported entity or \
                        attribute can cross the round trip only inside a SafeOpaqueCapsule \
                        that retains exact raw octets, token AST, source span/context, \
                        entity-reference closure and remap proof, profile/edition, digest, \
                        hazard classification, and explicit NoSemanticAuthority. Injection \
                        into output occurs only when structural safety and reference-closure \
                        premises reproduce; otherwise the capsule remains sidecar evidence \
                        and the exchange refuses lossless status",
            hypotheses: &[
                "the bounded Part-21 lexer recognizes the entire unsupported payload without interpreting vendor semantics and preserves exact source octets separately from normalized tokens",
                "all #references, nested aggregate/select/function tokens, external-reference strings, and context dependencies are structurally enumerated; reference remapping is total, injective on live identities, and independently replayed",
                "capsules are inert data: network access, URI resolution, signature trust, plugin/macro/script execution, executable foreign payloads, and unbounded decompression are forbidden",
                "byte-exact reinsertion and token-equivalent canonical reinsertion are distinct dispositions; any changed context, unresolved dependency, profile mismatch, token failing the frozen reinsertion-safety policy, or signature-covered region refuses reinsertion",
            ],
            qoi: "exact_opaque_octet_token_reference_closure_and_safety_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i07.oracle.opaque_capsule.v1 at fs-vmanifest-oracles/i07/opaque_capsule.rs::check_bytes_tokens_refs_and_hazards",
                independent: true,
                tcb_overlap: "shares the bounded lexical token schema and fixture bytes; independent byte ledger, reference walker, remapper, and deny-all side-effect harness do not use production capsule code",
            },
            activation: "all baseline claims are adjudicated and the safe-opaque maximal feature is enabled",
            kill: "one byte/token mismatch, missing dependency, noninjective remap, active side effect, reinsertion without a valid structural-safety receipt, signature laundering, or OpaquePreserved item granted semantic authority kills the capsule route",
            fallback: "retain the exact capsule only in the native evidence package, mark the construct Unsupported, and emit no AP242 artifact requiring it",
            no_claim: "opaque preservation proves syntax/bytes/reference custody only; it never proves vendor meaning, safety, conformance, semantic equivalence, or digital-signature validity",
        },
        ClaimSpec {
            id: "i07-kinematic-pair-semantic-roundtrip",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "For the pinned AP242 kinematic subset, mechanisms, links, \
                        pair/joint types, ordered attachment frames, admissible \
                        configuration relations, freedoms, coordinates, limits, sense, \
                        coupling laws, initial branch, and occurrence associations \
                        round-trip to an equivalent native constraint graph and back",
            hypotheses: &[
                "fixed, revolute, prismatic, cylindrical, spherical, planar, screw, gear-ratio, rack-pinion, and declared compound fixtures use exact profile-defined semantics; unsupported pairs refuse",
                "attachment frames, motor composition order, handedness, axis direction, coordinate zero, pitch/ratio sign, angular/linear units, limit openness/closure, and configuration branch are explicit",
                "equivalence compares the admitted configuration relation and tangent/Pfaffian relation on the pinned regular domain, not labels, sampled motion, or equal nominal degree counts",
                "holonomic position constraints, nonholonomic velocity distributions, unilateral limits/contact, friction, clearance, compliance, actuation, and dynamics remain separately typed",
            ],
            qoi: "exact_kinematic_relation_frame_freedom_limit_and_coupling_equivalence_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i07.oracle.kinematics.v1 at fs-vmanifest-oracles/i07/kinematics.rs::check_relation_and_tangent_graph",
                independent: true,
                tcb_overlap: "shares PGA motor arithmetic specification and fixture semantic graph; independent symbolic relation and interval tangent checkers do not use production AP242 mapping",
            },
            activation: "baseline product occurrence/frame semantics are green and the kinematic maximal feature is enabled",
            kill: "one joint/pair/frame/sense/unit/limit/coupling/branch change, sampled-only equivalence, holonomic/nonholonomic conflation, silent clearance/contact inference, or ambiguous occurrence attachment called Preserved kills the route",
            fallback: "preserve kinematic entities as opaque/loss atoms and export the assembly without native motion authority",
            no_claim: "no mobility completeness, finite motion, collision freedom, force/reaction, friction, clearance, compliance, dynamics, controller, or mechanism safety claim",
        },
        ClaimSpec {
            id: "i07-partial-bidirectional-edition-migration",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "Each automated migration between two exact pinned AP242 \
                        edition/profile identities is a partial bidirectional lens over \
                        explicitly declared common, refined, deprecated, split, merged, \
                        and unrepresentable semantic atoms. On its admitted domain the \
                        lens laws and identity/lineage preservation reproduce; outside it \
                        migration returns itemized information loss, ambiguity, or refusal \
                        rather than inventing an inverse",
            hypotheses: &[
                "both schema/profile artifact roots, conformance classes, mapping-table root, semantic-atom crosswalk, defaults, deprecations, and migration policy are frozen",
                "GetPut, PutGet, and PutPut laws are stated only on the law's exact admitted domain/quotient; noninjective split/merge/defaulting and edition-specific semantics carry residual information or refuse",
                "migration composes product/configuration lineage, geometry bounds, PMI, material/evidence, harness, kinematic, and opaque/loss receipts rather than checking schema validity alone",
                "a migration path records every intermediate edition/profile and composed loss; direct and multi-hop routes compare under a pinned coherence defect",
            ],
            qoi: "exact_admitted_lens_laws_lineage_and_migration_loss_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i07.oracle.edition_migration.v1 at fs-vmanifest-oracles/i07/edition_migration.rs::replay_lens_and_crosswalk",
                independent: true,
                tcb_overlap: "shares frozen mapping-table and atom-schema bytes; independent lens-law enumerator and graph reconciler do not use production migration code",
            },
            activation: "baseline round-trip and semantic-loss claims are green for both endpoint profiles",
            kill: "one unpinned schema/default, false inverse, law failure inside the admitted domain, lost lineage, hidden intermediate loss, path-coherence defect outside band, or Unsupported case reported migrated kills that mapping revision",
            fallback: "single-profile operation plus an explicit human-reviewed semantic-difference package; no automatic migration",
            no_claim: "no universal cross-edition compatibility, legal equivalence, semantic recovery of removed constructs, migration of proprietary history, or guarantee for an unpinned schema/profile",
        },
        ClaimSpec {
            id: "i07-governed-real-ap242-corpus-validation",
            ambition: Ambition::Frontier,
            polarity: ClaimPolarity::Affirmative,
            statement: "For a successor-bound ProfileRegistryRoot and ProfileAssignmentRoot, a \
                        governed content-pinned real AP242 corpus campaign validates every declared \
                        item and stratum against its assigned exact ProfileKeyDigest/ProfileEntryRoot, \
                        independently owned expected semantic inventory, named reader/writer/vendor/\
                        version matrix, complete semantic-loss receipt, and scoped conformance and \
                        interoperability predicates. StandardAxis, InteroperabilityAxis, and \
                        IndustrialCorpusAxis remain separately serialized and no axis is inferred \
                        from parser acceptance, synthetic evidence, or another axis",
            hypotheses: &[
                "fs-vvreg discharges the external corpus waiver only through a typed I07GovernedCorpusDischargeEnvelopeV1 and verified DischargeReceipt binding the exact waiver predicate, license/access/custody, standard/application-protocol/Part-21 edition and schema/profile roots, immutable corpus item roots, redaction policy, independent review/revocation receipts, predecessor manifest and transaction intent; the verified AmendmentRecord separately binds the final successor and the atomic transaction advances the authority head",
                "ProfileAssignmentRoot transitively binds ProfileDefinitionSchemaRoot and RequiredAssignmentUniverseRoot; its decoded logical-key domain exactly equals every required leaf, claim/domain, fixture/case/corpus-item/stratum and semantic profile role, each with exactly one immutable ProfileKeyDigest/ProfileEntryRoot assignment; every logical row runs and the caller cannot select a favorable profile or corpus subset",
                "expected semantic inventories, acceptable loss dispositions, QoIs, bands, standard clauses or scoped predicates, vendor/tool versions, independence threat graph, and adjudication code are frozen before candidate access",
                "StandardAxis=ScopedClauseConformance requires its exact licensed/public clause-pack and conformance-class root; InteroperabilityAxis=ScopedRequiredMatrixValidated requires every preregistered tool/version/direction cell; IndustrialCorpusAxis=GovernedBlindScopedValidation additionally requires an untouched-population custody and joint candidate/checker/expectation commitment, while a previously public corpus remains PublicCorpusOnly",
                "each component result binds the exact manifest, claim, checker, profile registry, profile assignment, corpus, toolchain, and output roots; missing, inaccessible, malformed, selectively omitted, or stale components are Unknown/IntegrityFailed rather than success",
            ],
            qoi: "exact_governed_corpus_coverage_semantic_receipt_and_scoped_validation_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G2,
            oracle: OracleRoute {
                identity: "i07.oracle.governed_real_corpus.v1 at fs-vmanifest-oracles/i07/governed_real_corpus.rs::check_custody_profiles_semantics_and_scoped_axes",
                independent: true,
                tcb_overlap: "shares only approved corpus/profile/expected-inventory digests and public predicate definitions; independent custody, reference readers, semantic adjudicators, reconciliation checker, and review signers are disjoint from the production importer/exporter",
            },
            activation: "one atomic fs-vvreg authority transaction installs the same-ID typed external discharge-envelope root, removes the external-corpus waiver, freezes realized governed corpus roots plus exact profile assignments, verifies the AmendmentRecord, advances the authority head, and only then grants one-shot capabilities to the named validation leaf",
            kill: "one custody/profile/semantic-role/assignment mismatch, repeated entry disguised by another ordinal, omitted required item or matrix cell, semantic inventory or loss-receipt defect, reference-reader disagreement outside its frozen band, public-as-blind relabeling, axis collapse, stale expectation, or favorable result from missing/inaccessible evidence kills the scoped validation claim",
            fallback: "retain the governed evidence and component disagreements with StandardAxis/InteroperabilityAxis/IndustrialCorpusAxis at their honest nonpromoting states; synthetic lanes remain independently adjudicable",
            no_claim: "one governed campaign proves only its exact corpus, profiles, predicates, tools, versions, strata, and bands; it grants no universal AP242 conformance, vendor universality, arbitrary industrial-data support, manufacturing fitness, drawing approval, supplier acceptance, safety, airworthiness, or regulatory authority",
        },
        ClaimSpec {
            id: "i07-proof-carrying-bidirectional-semantic-equivalence",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "A machine-checked theorem card proves that the frozen AP242 \
                        import functor and export functor form an equivalence on the exact-\
                        preservation supported semantic subcategory, modulo only the exact \
                        declared normalization congruence. Separately, approximate paths \
                        lift to a typed directed-error/loss quantaloid whose comparison \
                        morphisms compose certified geometry and semantic-loss budgets \
                        without treating a finite tolerance as equality. Unit/frame, \
                        identity/lineage, product/configuration, PMI, material/evidence, \
                        harness, kinematic, geometry-budget, and loss naturality squares \
                        commute; exact inverse laws and enriched comparison laws construct \
                        the complete semantic-difference receipt rather than assuming it",
            hypotheses: &[
                "the exact source AP242 and target Machine-IR semantic categories are finite typed attributed relational structures with explicit objects, morphisms, equations, partiality/refusal, and a proved normalization congruence",
                "the approximation layer is a typed directed-error/loss quantaloid whose homs are complete lattices in a frozen authority order and whose unital associative composition preserves arbitrary joins in each variable; identity/zero budgets, explicit infinity/refusal, property-specific geometry components, semantic-loss colors, and independently checked composition bounds are machine data; tolerance closeness is never silently made symmetric, transitive, or equal",
                "import/export object and morphism maps, unit/frame conversions, identity/lineage maps, and every domain-specific naturality square are complete machine ASTs with total serialization and runtime-premise bindings",
                "ordinary equivalence is restricted to the exact-preservation supported subcategory and exact quotient; opaque capsules, ambiguous atoms, unsupported constructs, and migration residuals are outside or carried as explicit loss objects, while certified approximations inhabit only the enriched comparison layer",
                "a canonical binding receipt proves byte/semantic equality among manifest claim, theorem-card AST, definitions, generated Lean proposition, elaborated declaration, proof term, environment, and complete transitive axiom closure under exact policy i07.lean-axioms.v1={propext, Quot.sound, Classical.choice}",
            ],
            qoi: "independent_bidirectional_semantic_equivalence_theorem_checker_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i07.oracle.semantic_equivalence_lean.v1 at proofs/i07/SemanticEquivalence.lean::supportedRoundTripEquivalence checked by pinned Lean4 kernel receipt",
                independent: true,
                tcb_overlap: "shares bound theorem-card AST, semantic signatures, and rational geometry-relation data only; proof kernel is disjoint from production import/export",
            },
            activation: "the target card is frozen and a pre-candidate manifest successor freezes the complete proposition/definition/runtime-premise/translation AST required by the formal-projection fixture before any theorem attempt",
            kill: "binding or kernel rejection, sorryAx, custom postulate/theorem-equivalent axiom, native-oracle authority outside the admitted kernel/checker TCB, transitive axiom outside exactly {propext, Quot.sound, Classical.choice}, failed naturality square, nonassociative/nonmonotone or under-enclosing budget composition, tolerance-as-equality laundering, or one independently verified premise-satisfying counterexample refutes exactly the bound theorem revision",
            fallback: "retain independently checked per-domain receipts with no global semantic-equivalence theorem",
            no_claim: "manifest version 1 prose mints no theorem color; checker success covers only the frozen exact category/congruence and separately typed budget-enriched approximation layer, never continuum equality from tolerance, and does not prove runtime premises, arbitrary vendor AP242, legal equivalence, or manufacturing/regulatory authority",
        },
        ClaimSpec {
            id: "i07-semantic-sheaf-descent-obstruction-theorem",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Affirmative,
            statement: "A machine-checked finite semantic sheaf-and-stack theorem organizes local \
                        AP242 modules, product/configuration views, geometry/PMI scopes, \
                        material/evidence contexts, harness subgraphs, and kinematic \
                        neighborhoods over the complete nerve of a pinned cover. In the proved \
                        effective normalization quotient, compatible strict local sections glue \
                        uniquely. Gauge-bearing objects instead satisfy separately typed groupoid/\
                        stack descent, including morphism descent and uniqueness only up to the \
                        exact retained equivalence/automorphism structure, while every \
                        failed restriction square localizes a typed obstruction; an \
                        additive pinned-cover Cech H1 obstruction is named only from a retained closed \
                        non-exact cochain witness in an admitted abelian coefficient \
                        complex, while nonabelian/categorical descent records a separately \
                        typed Cech cocycle, torsor, or groupoid obstruction",
            hypotheses: &[
                "cover, complete nerve with every nonempty finite intersection, stalk signatures, restriction maps, restriction-composition identities, pair/triple/higher-overlap cocycle orientations, coefficient category, gauges, normalization quotient, and domain-specific semantic equations are complete finite machine data",
                "restriction maps preserve typing, units/frames, identities/lineage, associations, and loss ownership; local compatibility is not inferred from equal labels or geometry proximity",
                "the strict H0 gluing theorem proves existence and literal uniqueness only after effectiveness of the normalization quotient; the groupoid/stack theorem separately proves object and morphism descent up to its precise equivalence while retaining stabilizers/automorphisms; failure produces localized restriction defects rather than a decorative cohomology label",
                "any H1 statement supplies exact closure and independent non-exactness against the full image of the previous coboundary, with basis/change-of-cover maps and coefficient semantics pinned",
                "identification of pinned-cover Cech H1 with derived sheaf H1 requires a separately kernel-checked Leray/acyclic-cover comparison, cofinal direct-limit theorem, or another exact comparison theorem; without it the result remains cover-relative",
                "an additive H1 label requires an admitted abelian-group/module coefficient complex and independently checked d^2=0; nonabelian or category-valued data use typed Cech cocycle/torsor/groupoid descent and are never laundered into vector-space H1",
                "a canonical binding receipt links manifest bytes, theorem-card AST, definitions, generated/elaborated Lean declaration, proof term, environment, and complete transitive axiom closure under i07.lean-axioms.v1={propext, Quot.sound, Classical.choice}",
            ],
            qoi: "independent_semantic_descent_gluing_and_obstruction_checker_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i07.oracle.semantic_descent_lean.v1 at proofs/i07/SemanticDescent.lean::finiteDescentAndObstruction checked by pinned Lean4 kernel receipt",
                independent: true,
                tcb_overlap: "shares theorem-card signature/restriction bytes only; independent finite matrix/cochain checker and Lean kernel are disjoint from production merge/mapping code",
            },
            activation: "baseline semantic-loss localization is green and a pre-candidate successor freezes the complete sheaf/cochain proposition and total AST-to-Lean translation",
            kill: "binding/kernel rejection, a non-functorial restriction, omitted higher intersection, pair/triple/higher-overlap coherence failure, false strict gluing before an effective quotient, erased stack automorphism, missed restriction defect, unproved non-exactness, d^2 != 0, abelian-H1 laundering of nonabelian data, cover-relative Cech H1 relabeled as derived sheaf H1 without a comparison theorem, basis/change-of-cover failure, forbidden axiom, or an independently verified in-domain counterexample refutes the bound theorem",
            fallback: "retain the ordinary typed graph diff and localized restriction failures with no H0/H1 authority",
            no_claim: "local agreement does not imply global consistency without full-nerve descent premises; uniqueness up to gauge is not ordinary sheaf uniqueness unless an effective quotient is proved, a mismatch is not automatically H1, pinned-cover Cech H1 is not derived sheaf H1 without a proved comparison, nonabelian descent is not additive cohomology, and sheaf/stack/cohomology authority grants no legal, manufacturing, safety, or arbitrary vendor semantic equivalence",
        },
        ClaimSpec {
            id: "i07-semantic-equivalence-counterexample-search",
            ambition: Ambition::Moonshot,
            polarity: ClaimPolarity::Refutation,
            statement: "Falsifier lane: exhaust a cardinality-audited canonical \
                        microgrammar and search a separately seeded larger supergrammar \
                        for supported-profile entity graphs, reference SCCs, assemblies, \
                        geometry relations, PMI, material/evidence, harness, kinematic, \
                        opaque, migration, and descent cases that violate a bound maximal \
                        theorem, lens law, loss reconciliation, or nonvacuity obligation. \
                        Every valid candidate and minimized counterexample is retained",
            hypotheses: &[
                "the adversary target freezes a raw 16x16x16x16x16x16 decorated-record space (N=16,777,216) plus a larger seeded supergrammar, but manifest-version-1 prose has no exhaustive authority",
                "before candidate generation a successor manifest freezes the full grammar/encoding/validity/stratum/tag/parameter AST, canonical quotient, total enumeration order, rank/unrank/sharding algorithms, source digests, independent decoder, bijection/cardinality proof, cost preflight, and Merkle completeness root",
                "independent premise and domain checkers bind each candidate to exact manifest/profile/theorem/migration/axiom-policy roots; out-of-domain candidates never count as refutation",
                "nonvacuity floors require witnesses for every baseline semantic domain, every maximal domain, nontrivial reference SCCs, lossy migration residuals, nonzero restriction defects, nontrivial directed-error/loss joins and compositions, genuine closed exact and closed non-exact abelian cochains, nonabelian cocycle/torsor/groupoid descent cases, and both theorem-satisfying and theorem-violating raw candidates",
            ],
            qoi: "exact_nonvacuity_completeness_and_zero_verified_in_domain_counterexample_verdict",
            unit: "bit",
            tolerance: ToleranceSemantics::Exact,
            evidence_tier: GauntletTier::G3,
            oracle: OracleRoute {
                identity: "i07.oracle.semantic_falsifier.v1 at fs-vmanifest-oracles/i07/semantic_falsifier.rs::verify_membership_cover_and_minimize",
                independent: true,
                tcb_overlap: "shares canonical candidate and declaration digests; independent enumerator, decoder, premise checker, proof kernels, high-precision geometry checker, coverage adjudicator, and minimizer are separately pinned",
            },
            activation: "all maximal declarations and profile roots are frozen and a pre-candidate FrozenManifest::amend successor has discharged the complete microgrammar formalization gate",
            kill: "first independently reproduced in-domain counterexample refutes exactly its bound immutable revision; empty domain, missed floor, grammar escape, rank/unrank/canonicalization defect, incomplete shard/root, or declaration mismatch fails the campaign rather than passing",
            fallback: "restrict or replace the theorem/lens/profile through authenticated amendment while preserving candidates, counterexamples, and tombstone/defect classifications",
            no_claim: "even a complete, nonvacuous, empty bounded search is not a proof; kernel acceptance of the exactly bound theorem remains separate, and a candidate caused by specification/checker/TCB defect invalidates that artifact rather than being mislabeled mathematical refutation",
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i07_fixtures() -> Vec<FixturePin> {
    vec![
        FixturePin {
            id: CAMPAIGN_POLICY_FIXTURE,
            source: FixtureSource::AuthoredSpec {
                spec: "I07_CAMPAIGN_POLICY_V1. AUTHORITY: freeze is preregistration only; \
no evidence color, AP242 conformance, drawing approval, supplier acceptance, manufacturing \
fitness, airworthiness, regulatory approval, or legal conclusion is minted. PROFILE_BINDING: \
manifest version 1 freezes the registry schema and synthetic target only, not a selectable \
execution profile. Before any candidate or evidence case, a FrozenManifest::amend successor \
must freeze a complete machine ProfileRegistryInstance, ProfileDefinitionSchemaRoot, \
ProfileRegistryRoot, RequiredAssignmentUniverseRoot, ProfileAssignmentRoot, schema/entity/type/\
attribute/semantic-atom inventories, numeric resource limits, canonical writer \
policy, positive per-domain nonvacuity floors, and required claim/leaf set. PROFILE_HASH_ALGEBRA: \
H(context,payload)=BLAKE3::derive_key(context,payload) with the context's exact ASCII bytes and raw \
32-byte output. U64LE and U32LE are fixed-width unsigned little-endian. FRAME_BYTES(x)=U64LE(len(x)) \
|| x; FRAME_UTF8(s)=FRAME_BYTES(exact_case_sensitive_UTF8_without_normalization(s)); \
FRAME_LIST(xs)=U64LE(count(xs)) || concat(FRAME_BYTES(x) for x in xs). No native-endian integer, \
hex-text digest, Unicode normalization, case folding, omitted empty field, or ambiguous concatenation \
is accepted. TAGGED_FIELD(tag,payload)=U32LE(tag)||FRAME_BYTES(payload), where tag is the one-based \
field ordinal fixed by the named schema. ProfileKeyBytes is FRAME_UTF8('i07-profile-key-v1') \
followed by exactly these nineteen TAGGED_FIELD entries, once each in increasing ordinal order: \
1 part21_encoding_edition_id Utf8; 2 part21_encoding_artifact_digest Digest; 3 \
ap242_application_protocol_id Utf8; 4 ap242_edition_id Utf8; 5 ap242_schema_artifact_digest \
Digest; 6 conformance_class_id Utf8; 7 implementation_method_id Utf8; 8 character_encoding_id \
Utf8; 9 signature_policy_id Utf8; 10 external_reference_policy_id Utf8; 11 \
supported_entity_inventory_digest Digest; 12 supported_semantic_atom_inventory_digest Digest; 13 \
validation_rule_inventory_digest Digest; 14 allowed_reference_scc_digest Digest; 15 \
unit_frame_policy_digest Digest; 16 canonical_writer_policy_digest Digest; 17 \
no_claim_policy_digest Digest; 18 coverage_floor_digest Digest; 19 limit_policy_digest Digest. \
Utf8 payloads use FRAME_UTF8, Digest payloads are exactly 32 nonzero raw bytes, and no field may be \
missing, empty, duplicated, unknown, reordered, normalized, case-folded or trailed. \
PROFILE_BOOTSTRAP_CAPS: every untrusted length/count is checked against remaining input, its local \
cap and the cumulative projection cap with checked arithmetic before allocation, recursion or hash \
finalization. MAX_PROFILE_UTF8_BYTES=65536; MAX_PROFILE_NESTED_COLLECTION_ITEMS=65536; \
MAX_PROFILE_NESTED_PAYLOAD_BYTES=16777216; MAX_PROFILE_DEFINITION_SCHEMA_BYTES=4194304; \
MAX_CANONICAL_PROFILE_DEFINITION_BYTES=16777216; MAX_PROFILE_DEFINITION_FIELDS=4096; \
MAX_PROFILE_REGISTRY_ENTRIES=4096; MAX_PROFILE_REGISTRY_PROJECTION_BYTES=67108864; \
MAX_REQUIRED_ASSIGNMENT_KEYS=1048576; MAX_REQUIRED_ASSIGNMENT_UNIVERSE_BYTES=268435456; \
MAX_PROFILE_ASSIGNMENT_BYTES=536870912. These immutable bootstrap limits are enforced before the \
decoder can trust or consult limit_policy_digest. MAX_PROFILE_NESTED_PAYLOAD_BYTES bounds byte_len(x) \
inside FRAME_BYTES(x), excluding that frame's outer eight-byte U64LE length; the complete frame size is \
checked_add(8,byte_len(x)), and every enclosing cumulative projection cap applies to that complete \
framed size. Top-level ProfileKeyBytes, ProfileDefinitionSchemaBytes, \
CanonicalProfileDefinitionBytes, ProfileRegistry, RequiredAssignmentUniverse and ProfileAssignment \
projections each have depth 0; entering any framed nested record or collection increments depth by \
exactly one, depth 64 is admitted, and an attempt to enter depth 65 is refused before recursion or \
allocation. Registry and assignment \
decoders stream into their domain-separated hashes instead of requiring whole-projection allocation, \
and decoded row counts must exactly match their framed counts. A concrete successor must ship two \
independent decoders plus exact-payload-cap, payload-cap-plus-one, exact-depth-64, depth-65, U64::MAX \
length/count, checked-add/multiply overflow, truncated nested payload, duplicate/order and trailing-byte \
mutation tests before any profile authority is used. \
ProfileDefinitionSchemaBytes is the exact content-addressed canonical schema artifact whose bytes \
start with FRAME_UTF8('i07-profile-definition-schema-v1') and freeze the complete definition-field \
ordinal/name/type/requiredness/default/constraint-AST table plus the canonical encodings of every \
nested type. ProfileDefinitionSchemaRoot=H('org.frankensim.i07.profile-definition-schema.v1',\
ProfileDefinitionSchemaBytes). CanonicalProfileDefinitionBytes is \
FRAME_UTF8('i07-profile-definition-v1')||raw(ProfileDefinitionSchemaRoot)||U32LE(field_count) \
followed by exactly field_count definition-AST fields, each as one TAGGED_FIELD in exact strictly \
increasing schema ordinal order; field_count must equal the exact schema table cardinality. Every \
empty field is present, nested sequences use FRAME_LIST, nested records are themselves framed, maps \
are lexicographically raw-key sorted and unique, booleans are one byte 0 or 1, and every integer width \
is fixed by the bound schema and serialized unsigned little-endian. The schema bytes/root and complete \
definition-field table are mandatory content-addressed inputs; a missing, mismatched or unavailable \
schema root is IntegrityFailed, and an encoder may not substitute host serialization, an omitted \
default or a friendly-name projection. ProfileKeyDigest=H(\
'org.frankensim.i07.profile-key.v1',ProfileKeyBytes). ProfileEntryRoot=H(\
'org.frankensim.i07.profile-entry.v1',raw(ProfileKeyDigest)||raw(ProfileDefinitionSchemaRoot)||\
FRAME_BYTES(CanonicalProfileDefinitionBytes)). ProfileRegistryRoot=H(\
'org.frankensim.i07.profile-registry.v1',FRAME_UTF8('i07-profile-registry-v1')||\
raw(ProfileDefinitionSchemaRoot)||FRAME_LIST(\
lexicographically_raw_sorted_unique(raw(ProfileKeyDigest)||raw(ProfileEntryRoot)) pairs)). A duplicate \
ProfileKeyDigest is IntegrityFailed even if paired with another entry root. AssignmentTargetKind=\
{Fixture=1,Case=2,CorpusItem=3,Stratum=4}; ProfileRole={Input=1,Output=2,SourceEndpoint=3,\
TargetEndpoint=4,ReferenceOracle=5,MigrationIntermediate=6,TheoremInstantiation=7}. \
AssignmentLogicalKey=(leaf_id,claim_or_domain_id,target_kind,target_id,profile_role), encoded as \
exactly 1 leaf_id Utf8; 2 claim_or_domain_id Utf8; 3 target_kind U32LE; 4 target_id Utf8; 5 \
profile_role U32LE, each as one TAGGED_FIELD in increasing ordinal order; target kind and profile \
role payloads are their U32LE discriminants. \
AssignmentRow=FRAME_BYTES(encoded_logical_key)||raw(ProfileKeyDigest)||raw(ProfileEntryRoot). \
RequiredAssignmentUniverseRoot=H('org.frankensim.i07.required-assignment-universe.v1',\
raw(ProfileDefinitionSchemaRoot)||FRAME_LIST(lexicographically_raw_sorted_unique(\
encoded_logical_key))). ProfileAssignmentRoot=H('org.frankensim.i07.profile-assignment.v1',\
raw(ProfileRegistryRoot)||raw(RequiredAssignmentUniverseRoot)||FRAME_LIST(\
lexicographically_raw_sorted_unique(AssignmentRow))). The decoded AssignmentRow logical-key set must \
equal, not merely contain, the exact RequiredAssignmentUniverseRoot set, and every key has exactly one \
row. Duplicate, missing or unexpected logical keys, unknown roles/kinds, a repeated entry disguised by \
another ordinal, mismatched embedded keys, or an entry absent from the exact registry are \
IntegrityFailed. Ordinals have no profile-authority semantics. --profile <profile-root> means exactly \
the assigned ProfileEntryRoot, not the registry root. Every request/replay receipt binds manifest \
root, profile-definition-schema root, registry root, required-assignment-universe root, assignment \
root, leaf, claim/domain, target kind/id, semantic profile role, ProfileKeyDigest, ProfileEntryRoot, \
and schema/membership/assignment proofs. Every required assignment row \
runs; no caller-selected favorable member, role, item, stratum, or subset exists. FILE_SCHEMA must \
match the already bound profile; profile/edition inference from FILE_NAME, host filename, vendor, \
FILE_SCHEMA guessing, or entity coincidence is forbidden. REQUEST_SET: the required claim and leaf \
set is frozen before execution; exit 0 means all required authorities passed, never a favorable \
caller-selected subset. ORTHOGONAL_STATES: \
ExecutionState={Succeeded,Failed,Cancelled,TimedOut,BudgetExhausted,InfrastructureFailed,\
IntegrityFailed}; ClaimState={Pending,Passed,Failed,Refuted,Unknown,Unsupported}; \
EvidenceCompleteness={CompleteEvidence,PartialEvidence,NoEvidence}; \
LossDisposition={PreservedExact,PreservedNormalized,SemanticallyMapped,BoundedApproximation,\
OpaquePreservedNoSemanticAuthority,Ambiguous,UnsupportedWithRefusal,DroppedWithRefusal}; \
StandardAxis={NotAssessed,SyntheticProfileOnly,ScopedClauseConformance}; \
InteroperabilityAxis={NotAssessed,ScopedRequiredMatrixValidated}; \
IndustrialCorpusAxis={NotAssessed,PublicCorpusOnly,GovernedBlindScopedValidation}; \
RegulatoryAxis={NotGranted}. AXIS_AGGREGATION: every non-NotAssessed state carries its exact scope, \
predicate/clause/matrix/corpus, profile assignment, checker, result and adjudication roots. Scoped \
clause conformance requires its clause pack; interoperability requires every preregistered matrix \
cell; only jointly committed untouched-population custody may mint GovernedBlindScopedValidation. \
Public/previously accessed corpora remain PublicCorpusOnly. Missing, inaccessible, failed, stale, or \
selectively omitted components cannot mint a favorable axis. Failed means the runner did not establish its acceptance \
predicate; Refuted requires an independently reproduced in-domain counterexample to the exact \
immutable claim; Unknown means evidence cannot decide; Unsupported is outside admitted \
capability. CROSS_AXIS_VALIDITY: case.finalized forbids ClaimState=Pending. ClaimState=Passed and every \
favorable non-NotAssessed validation axis require ExecutionState=Succeeded, \
EvidenceCompleteness=CompleteEvidence, valid integrity, admitted applicability/support and every \
frozen acceptance predicate Passed. ClaimState=Refuted likewise requires CompleteEvidence plus the \
independently reproduced exact immutable counterexample and ExecutionState=Succeeded. A complete \
supported acceptance failure is ClaimState=Failed; ClaimState=Unsupported requires CompleteEvidence \
of exact non-applicability. PartialEvidence or NoEvidence forces ClaimState=Unknown. \
ExecutionState=Failed, Cancelled, TimedOut, BudgetExhausted, InfrastructureFailed or IntegrityFailed \
cannot carry Passed, Refuted or a favorable validation axis. PartialEvidence and NoEvidence have zero \
promotion authority. INTEGRITY and \
INFRASTRUCTURE never become scientific failure/refutation. EXIT_CODES: \
0=all_requested_claims_passed; 10=failed; 11=refuted; \
12=unknown_or_partial; 13=unsupported; 20=cancelled; 21=timed_out; 22=budget_exhausted; \
30=integrity_failed; 31=infrastructure_failed. Mixed cases return the highest-severity \
fail-closed code in order 30,31,21,22,11,10,20,13,12,0 while retaining every per-claim state. \
LIFECYCLE: request -> admit/refuse -> start -> bounded work -> request_cancel? -> drain -> finalize -> adjudicate -> \
atomic_publish. A successful publication requires one complete atomic AP242/loss-receipt pair; a \
non-success publishes no AP242 artifact and records the explicit absent-artifact disposition in at \
most one atomic terminal/FailureBundle evidence transaction. No path publishes one member of a \
success pair, a partial member, or an unrecorded absence. \
LIFECYCLE_EVENT_CONTRACT: every submitted attempt emits exactly one campaign.requested. A refused \
admission emits no campaign.admitted, case.started or work event and retains one typed refusal/\
finalization receipt. Every admitted attempt emits exactly one campaign.admitted and then exactly one \
case.started before work. Every observer or descendant registration atomically commits its catalog \
membership receipt and acquires its worst-case terminal row-and-encoded-byte reservation before it is \
admitted; a registration with no capacity cannot become live. The first accepted cancellation for one \
open scope/epoch is the primary request: exactly one case.cancel.requested binds request id, scope \
root, request-event identity, primary_request_logical_sequence, supervisor-clock request interval, \
calibration-artifact root and conservative timestamp-uncertainty enclosure, then atomically seals the \
request-time observer-set root/count, descendant-frontier root/count, registration epoch and exact \
terminal reservation. Child admission closes in the same transition: a racing registration either \
commits before the seal and is included with deterministic catalog membership and reservation \
receipts, or is refused and cannot escape either reconciliation. Exact-key byte-identical replay \
returns the primary receipt without another event; a conflicting replay is IntegrityFailed. A later \
distinct request for the sealed, draining or finalized epoch receives deterministic AlreadySealed, \
AlreadyDraining or AlreadyFinalized response bound to the primary event and cannot reopen drain or \
reset latency. On a cancellation path, exactly one case.drain.started follows the primary seal and \
initiates propagation. On a non-cancellation path, exactly one case.drain.started follows the \
earliest applicable work-completion, domain-failure, budget, campaign-timeout or infrastructure \
drain cause without inventing a cancellation seal or observation. \
Every live cancellable logical work unit in the sealed observer set emits exactly one \
case.cancel.observed binding the primary request id/event identity/primary_request_logical_sequence, \
sealed observer root/count, descendant-frontier root/count, registration epoch, admitted unit/tile id, \
deterministic membership proof, optional last-completed boundary, a distinct globally unique \
observation_event_logical_sequence and supervisor-domain calibrated observation interval with an \
outward uncertainty enclosure. Observation forbids every new domain, candidate, protected-access, \
adjudication or publication effect in that unit; only already reserved bounded cleanup, descendant \
join, exit-receipt, drain and finalization work may continue, and none may mint new scientific \
authority. \
On a cancellation path, case.drain.completed carries the sealed roots/counts and a complete \
exit-reconciliation root and may occur only after every sealed unit has observed cancellation and \
every admitted descendant has one matching membership-plus-exit receipt. On a non-cancellation path, \
it instead carries the exact admitted-frontier root/count and complete exit reconciliation and may \
occur only after every admitted descendant exits; it neither requires nor fabricates a cancellation \
observation. A legitimately empty sealed observer set requires a present, \
independently reconciled root with cardinality zero and exactly zero observations; its \
request-to-observation maximum is exactly zero at seal, but descendant drain is separately mandatory. \
DRAIN_TRIGGER: DrainTrigger={CancellationObserved,EmptyObserverSeal,ObservationTimeoutDrain,\
InfrastructureFailure,NonCancellationDrain}. Its calibrated trigger interval is respectively the \
earliest effective on-time observation interval, the primary-seal interval for an independently \
verified empty observer set, the singleton first nanosecond after the missed inclusive observation \
deadline, the independently receipt-verified infrastructure-failure-onset interval, or the explicit \
drain-start interval. Trigger choice orders effective interval lower endpoint first, causal rank \
InfrastructureFailure=0,CancellationObserved=1,EmptyObserverSeal=2,ObservationTimeoutDrain=3,\
NonCancellationDrain=4 \
second, causal logical sequence third and stable identity fourth. The inclusive observation deadline \
is checked_add(request_interval_lower,250000000 ns) for Core or \
checked_add(request_interval_lower,1000000000 ns) for Max, and its timeout onset is \
checked_add(deadline,1 ns); arithmetic overflow is IntegrityFailed. case.drain.completed carries a \
calibrated drained interval in the same supervisor domain, and case.finalized carries its calibrated \
finalized interval. Conservative trigger-to-drained latency is \
checked_sub(drained_interval_upper,trigger_interval_lower); conservative drained-to-finalized latency \
is checked_sub(finalized_interval_upper,drained_interval_lower). Core caps are inclusively 2000000000 \
ns and 2000000000 ns respectively; Max caps are inclusively 8000000000 ns and 8000000000 ns. The \
failure onset for either cap is exactly the first nanosecond outside its inclusive bound. An \
observation whose interval upper equals the inclusive observation deadline, or a drain/finalization \
conservative latency equal to its inclusive cap, is on time. At the derived timeout onset \
checked_add(deadline,1 ns), an observation completion loses to and latches TimedOut; likewise a \
conservative latency of checked_add(cap,1 ns) latches TimedOut. Overall admitted wall-budget expiry is \
TimedOut. \
TERMINAL_ARBITRATION: any malformed canonical prefix, invalid integrity binding or impossible interval \
selects IntegrityFailed regardless of an operational candidate. Otherwise the earliest \
terminal-eligible logical boundary wins; at the same boundary the exact cause order is \
InfrastructureFailed > TimedOut > Failed > Cancelled > BudgetExhausted > Succeeded. A cancellation \
request makes Succeeded and BudgetExhausted ineligible until the request is observed and the sealed \
frontier drains and finalizes. Cancelled is eligible only when every required observation, descendant \
exit, drain and finalization receipt is complete inside every inclusive cap and no higher cause is \
eligible. A domain/acceptance effect whose own durable commit point precedes the primary cancellation \
may supply Failed; no post-cancellation effect may do so. A receipt-verified supervisor, external \
process, sink, authentication, drain-protocol or publication-protocol failure supplies \
InfrastructureFailed. A missed observation, trigger-to-drained, drained-to-finalized or campaign \
deadline supplies TimedOut. An unverifiable failure without the required receipt sets \
ExecutionState=IntegrityFailed and EvidenceCompleteness=NoEvidence rather than a guessed scientific \
result. case.cancelled is the exactly-once \
terminal-selection event when Cancelled wins; case.failed is exactly once for selected Failed, \
TimedOut, BudgetExhausted, InfrastructureFailed or IntegrityFailed and binds that exact disposition. \
Every admitted attempt emits exactly one case.finalized binding final execution/input/integrity/\
domain/support state and immutable evidence payload, then exactly one claim.adjudicated per requested \
claim binds its claim/evidence axes. Successful adjudication requires exactly one \
evidence.atomic_publish_committed followed by exactly one artifact.published. A non-success may \
publish at most one paired atomic terminal/failure evidence transaction with zero promotion \
authority; a partial/unpaired publication is IntegrityFailed. Duplicate, missing, out-of-order or \
inapplicable lifecycle events, or work/cancellation/drain after case.finalized, never promote. \
For each exact sealed observer member i, conservative latency_i is \
checked_sub(observation_interval_upper_i,request_interval_lower) in nanoseconds in one supervisor \
monotonic domain after applying the receipt-bound calibration and outward drift enclosure; the reported \
latency is max_i(latency_i), inclusively capped at 250000000 ns for Core and 1000000000 ns for Max. \
Overflow, an invalid/expired calibration, an inward enclosure, incomparable domains, an observation \
interval wholly before its request interval, or unknown/duplicate/extra/root-, count- or \
epoch-mismatched observations is IntegrityFailed. A missing or late observation crossing its inclusive \
deadline is TimedOut; preflight inability to establish a common clock or bounded observer is an \
Unsupported admission refusal; a receipt-verified observer/supervisor loss is InfrastructureFailed. \
None can succeed, and incomplete persistence separately lowers EvidenceCompleteness. Raw calibrated \
timestamps, uncertainty and watchdog arrival remain telemetry; request \
event identities, logical sequences, unit/boundary ids, sealed root/epoch and selected disposition \
remain canonical. A legitimately empty admitted \
source and target inventory still publishes a present independently reconciled receipt with both \
roots and zero counts; failure or incomplete inventory can never publish a clean empty item set. \
CHECKPOINT: canonical \
version/type/profile/manifest/case/root/checksum envelope; resume and fork are idempotent; a stale \
or corrupt envelope is IntegrityFailed. LOSS_RECEIPT: complete source and target atom roots, exactly \
one DispositionRecord per source atom and exactly one OriginRecord per target atom; each record has a \
typed zero/one/many counterpart set and derivation/refusal rule, and the forward/inverse incidence \
relations must be exact transposes. Split, merge, synthesis and deletion are never forced into false \
one-to-one lineage. Approximation bounds, opaque capsules, warnings/refusals, owners, authority effects, and \
reconciliation counts. OBSERVABILITY: typed canonical fs-obs JSONL is a non-authoritative mirror. \
One event kind is at most 256 UTF-8 bytes; one complete record including envelope and LF is at most \
1048576 bytes; one durable flush is at most 1024 rows and 4194304 bytes; one run/scope is at most \
65536 rows and 67108864 bytes; the governor retains at most 268435456 bytes across concurrent runs. \
Attempt admission uses checked arithmetic to reserve base ordinary and terminal rows plus their \
worst-case encoded bytes including envelope, redaction and framing overhead. Each observer/descendant \
registration then atomically acquires its incremental worst-case cancellation/drain/finalize rows and \
bytes from the same global ledger before catalog admission; the cancellation seal freezes both catalogs \
and the exact non-borrowable reservation. Transition, observation, drain, adjudication, publication and \
                      failure rows are reserved separately. Every reservation follows exactly \
                      Provisional->Released for refused admission, or \
                      Provisional->Reserved->Consumed|Released for admitted capacity: a refused admission releases every \
                      provisional unit atomically; a durable row consumes its exact row/byte unit; finalization emits one \
reservation-settlement root reconciling every unit as consumed or released with a reason; and crash \
recovery replays the durable ledger to the same settlement without double consumption or leak. A unit \
cannot be released before its possible terminal row is durably made unnecessary, borrowed by another \
run, or retained after settlement. Priority evidence uses an independently serviced bounded writer \
and bounded priority-flush segments; ordinary traffic cannot consume that lane. Admission performs a \
checked global earliest-deadline-first demand-bound test: for every priority deadline d, the sum of \
remaining worst-case durable service time of all already admitted and proposed priority segments with \
deadline <=d must not exceed the remaining writer service capacity through d computed from the \
writer's frozen service rate, elapsed service and nonpreemptive blocking bound. Per-run \
segment_count*worst_case_durable_service_time is necessary but never sufficient. The \
writer's frozen nonpreemptive segment maximum also fits the smallest tier request-observation SLO, so \
priority traffic from concurrent runs cannot queue cancellation observations behind ordinary or other \
priority work beyond its admitted deadline. The global reservation ledger refuses per-run, governor, \
byte, row or service-time overcommit before allocation. High-volume progress is \
deterministically aggregated before event creation by frozen logical tile/shard/window with detail \
in content-addressed Merkle artifacts. No created event is silently dropped, overwritten, truncated, \
arrival-time sampled or spilled without bound. Sink refusal sets \
ExecutionState=InfrastructureFailed and EvidenceCompleteness=PartialEvidence when a durable partial \
FailureBundle exists; if no terminal evidence can persist it sets ExecutionState=InfrastructureFailed \
and EvidenceCompleteness=NoEvidence with zero authority. A slow sink consumes the wall budget; \
admission requires a durable-service bound compatible with the exact observation, drain and \
finalization SLOs. Redact licensed/proprietary payload, personal data, secrets, signatures and absolute host \
paths; log stable ids, domains, digests, counts, bounds, decisions and ranked fixes. RETENTION: canonical input/output native/AP242 \
artifacts, manifest/profile/oracle/TCB roots, semantic loss, proof/falsifier/adjudication receipts, \
checkpoints, logs, and FailureBundle are content-addressed; success and all terminal failures retain \
the minimum replay pack. INDEPENDENCE: threat graph declares shared parser, lexer, schema, semantic \
vocabulary, geometry kernel, units, graph algorithms, generator, data, compiler/toolchain, and human \
review; a second function name is not independence. HOLDOUT: public deterministic ranges are stage \
separation/replay only and have no IID/secrecy authority; access before the named consumer begins is \
IntegrityFailed by campaign policy but the manifest itself cannot enforce custody. \
HOLDOUT_REALIZATION: every version-1 HeldOut AuthoredSpec is only a generator/stratum schema. Before \
access, a FrozenManifest::amend successor replaces it with a content-addressed realized artifact root \
and binds generator/toolchain, canonical case count/order, Merkle root, custodian signatures, \
candidate/checker/ProfileRegistryRoot/ProfileAssignmentRoot, scales, tolerances, and decision code; \
until then no holdout capability or promotion authority exists. JOINT_REVEAL: a leaf consuming \
multiple protected strata gets one atomic signed commitment over every realized root plus candidate, \
checker, profile registry/assignment, scales, tolerances, and decision code before any root, byte, \
summary, or derived label is revealed. Sequential reveal, tuning between strata, selective opening, \
or retry is IntegrityFailed and splitting strata never creates another attempt. GOVERNED_EXTERNAL: \
public standard/conformance corpora may support only their exact StandardAxis or PublicCorpusOnly \
scope; industrial blind authority additionally requires untouched-population custody and the same \
joint commitment discipline over candidate, checker, expectations, profiles, corpus roots and decision \
code. The unresolved i07-governed-real-ap242-corpus deck and its Waiver remain \
NoPromotionAuthority until fs-vvreg verifies an I07GovernedCorpusDischargeReceiptV1 precommit \
authorization binding the exact \
waiver subject/predicate, typed artifact kind, license/access/custody, standard and profile identities, \
realized item/corpus/expectation roots, candidate/checker/toolchain/decision roots, review and \
revocation receipts, predecessor manifest digest, successor version and a domain-separated \
transaction-intent digest. One atomic FrozenManifest::amend authority transaction installs at that \
same deck id an I07GovernedCorpusDischargeEnvelopeV1 root containing the realized-root references and \
precommit authorization receipt, removes the Waiver row, verifies the AmendmentRecord and advances \
the authority head before protected access or adjudication. TRANSACTION_INTENT_V1 defines \
TransactionIntentDigest=BLAKE3::derive_key(\
'org.frankensim.i07.governed-transaction-intent.v1',P), where P is the canonical successor intent \
for exactly one fenced authority transaction. It binds the predecessor manifest and governance-stage \
receipt; expected authority-head digest/generation; immutable CandidateFrozen candidate/checker/\
toolchain/decision/ProfileRegistryRoot/ProfileAssignmentRoot commitment; lexically ordered retired \
waiver subjects; role-addressed (slot_id,artifact_role,artifact_schema_digest,protected_root) records; \
role-addressed future envelope slots that bind each output artifact_schema_digest before Pending; \
coupled transaction group; governance stage/scope; idempotency key/attempt/capability epoch; future \
envelope/final-successor/AmendmentRecord slots; an exact receipt-schema set root; and the exact closed \
GovernedOutputSchemaSetRoot. Separately sorted slot-id and root \
lists are forbidden because any role swap must change P. Cross-field validity requires \
initiative_id='I07', schema_identity='i07-governed-transaction-intent-v1', checked \
successor_version=predecessor_version+1, governance_stage=RealizationCommitted, the exact governed \
authority scope, an expected-head compare-and-swap with generation incremented exactly once and an \
idempotency/capability epoch that prevents ABA replay, and equality of tags 0x0012/0x0013 with the \
opened CandidateFrozen profile-root commitment; because ProfileRegistryRoot and \
ProfileAssignmentRoot transitively bind ProfileDefinitionSchemaRoot and \
RequiredAssignmentUniverseRoot, opening either transitive binding to a different root is \
IntegrityFailed. GovernedOutputSchemaSetRoot is frozen by GovernanceCommitted before candidate \
execution, copied byte-identically into CandidateFrozen, and must byte-equal tag 0x0014 in P. The \
retired-waiver list must equal the exact predecessor-minus-successor waiver-set \
difference. Protected bindings, retired subjects and envelope future artifacts form an exact \
role-addressed bijection. Each future output artifact schema is a mandatory content-addressed canonical \
schema whose bytes are available and independently rehashed before realization, and the realized \
envelope must decode under the byte-identical precommitted schema digest; a caller-selected or merely \
nonzero schema digest has no authority. Every envelope-slot id is unique and the envelope-slot-id set, \
final-successor slot id and AmendmentRecord slot id are pairwise disjoint. \
GOVERNED_OUTPUT_SCHEMA_AUTHORITY_V1: the protocol authority is closed to exactly one output and has no \
extension row. GovernedOutputSchemaMatrixBytesV1 is byte-for-byte \
FRAME_UTF8('i07-governed-output-schema-matrix-v1')||U16LE(1)||U16LE(1)||U16LE(3)||\
FRAME_BYTES(Utf8('i07.governed-real-ap242-corpus-discharge'))||\
FRAME_BYTES(Utf8('i07-governed-real-ap242-corpus'))||\
FRAME_BYTES(Utf8('GovernedRealAp242CorpusDischarge/'))||U64LE(1)||U64LE(1)||U16LE(1)||\
U16LE(1): row_count=1, protocol_kind=GovernedRealAp242CorpusDischarge=1, \
governance_stage=RealizationCommitted=3, \
authority_scope='i07.governed-real-ap242-corpus-discharge', \
target_slot_id='i07-governed-real-ap242-corpus', \
role_prefix='GovernedRealAp242CorpusDischarge/', min_count=max_count=1, \
related_waiver_rule=ExactRetiredSubject=1, and protected_binding_rule=Required=1. The exact matrix \
length is (8+36)+3*2+(8+8+40)+(8+8+30)+(8+8+33)+2*8+2*2=221 bytes. The sole \
retired waiver subject, protected slot_id, target_slot_id and FutureArtifact related_waiver_subject must all \
byte-equal i07-governed-real-ap242-corpus. Its only artifact_role is derived, never supplied by a \
caller: artifact_role='GovernedRealAp242CorpusDischarge/'||related_waiver_subject, exactly \
'GovernedRealAp242CorpusDischarge/i07-governed-real-ap242-corpus'. \
GovernedOutputSchemaMatrixDigest=BLAKE3::derive_key(\
'org.frankensim.i07.governed-output-schema-matrix.v1',GovernedOutputSchemaMatrixBytesV1). \
I07GovernedCorpusDischargeEnvelopeSchemaBytesV1 starts with \
FRAME_UTF8('i07-governed-corpus-discharge-envelope-schema-v1')||U16LE(34) and then exactly thirty-four \
SchemaField records in increasing ordinal order. SchemaField=U16LE(ordinal)||\
FRAME_BYTES(Utf8(field_name))||U16LE(field_type)||U16LE(1)||FRAME_BYTES(ConstraintAst), where the \
fourth component is required=1; field_type is Utf8=1,U16LE=2,U64LE=3,Digest=4; and ConstraintAst is \
exactly EqUtf8=U16LE(1)||FRAME_BYTES(Utf8(literal)), EqU16=U16LE(2)||U16LE(value), \
NonzeroDigest=U16LE(3), EqIntentTag=U16LE(4)||U16LE(tag), \
RangeU64=U16LE(5)||U64LE(min)||U64LE(max), or \
EqCandidateFrozenComponent=U16LE(6)||U16LE(component), with \
Candidate=1,Checker=2,Toolchain=3,Decision=4; EqTransactionIntentDigest=U16LE(7); \
EqPrecommitAuthorizationDigest=U16LE(8); \
EqAuthorizationBinding=U16LE(9)||FRAME_BYTES(Utf8(binding_name)); or \
EqAuthorizationBoundU64Range=U16LE(10)||FRAME_BYTES(Utf8(binding_name))||U64LE(min)||U64LE(max). \
The authorization-binding namespace is closed to exactly the seventeen case-sensitive binding_name \
values used by fields 14..23 and 28..34 below; each resolves one unique mandatory field of the \
verified precommit authorization, and a missing, extra, duplicate, aliased or differently typed \
binding is IntegrityFailed. \
The exact field table is: \
1 schema_identity Utf8 EqUtf8('i07-governed-corpus-discharge-envelope-v1'); \
2 initiative_id Utf8 EqUtf8('I07'); 3 protocol_kind U16LE EqU16(1); \
4 governance_stage U16LE EqU16(3); \
5 authority_scope Utf8 EqUtf8('i07.governed-real-ap242-corpus-discharge'); \
6 slot_id Utf8 EqUtf8('i07-governed-real-ap242-corpus'); \
7 related_waiver_subject Utf8 EqUtf8('i07-governed-real-ap242-corpus'); \
8 artifact_role Utf8 \
EqUtf8('GovernedRealAp242CorpusDischarge/i07-governed-real-ap242-corpus'); \
9 transaction_intent_digest Digest EqTransactionIntentDigest; \
10 precommit_authorization_digest Digest EqPrecommitAuthorizationDigest; \
11 governed_output_schema_set_root Digest EqIntentTag(0x0014); \
12 profile_registry_root Digest EqIntentTag(0x0012); \
13 profile_assignment_root Digest EqIntentTag(0x0013); \
14 standard_application_protocol_part21_identity_root Digest \
EqAuthorizationBinding('standard_application_protocol_part21_identity_root'); \
15 conformance_class_root Digest EqAuthorizationBinding('conformance_class_root'); \
16 license_access_custody_root Digest EqAuthorizationBinding('license_access_custody_root'); \
17 corpus_item_root Digest EqAuthorizationBinding('corpus_item_root'); \
18 corpus_item_count U64LE \
EqAuthorizationBoundU64Range('corpus_item_count',1,1048576); \
19 expectation_root Digest EqAuthorizationBinding('expectation_root'); \
20 semantic_inventory_root Digest EqAuthorizationBinding('semantic_inventory_root'); \
21 qoi_band_root Digest EqAuthorizationBinding('qoi_band_root'); \
22 clause_interoperability_matrix_root Digest \
EqAuthorizationBinding('clause_interoperability_matrix_root'); \
23 public_governed_blind_authority_root Digest \
EqAuthorizationBinding('public_governed_blind_authority_root'); \
24 candidate_root Digest EqCandidateFrozenComponent(1); \
25 checker_root Digest EqCandidateFrozenComponent(2); \
26 toolchain_root Digest EqCandidateFrozenComponent(3); \
27 decision_root Digest EqCandidateFrozenComponent(4); \
28 joint_candidate_checker_expectation_root Digest \
EqAuthorizationBinding('joint_candidate_checker_expectation_root'); \
29 oracle_tcb_threat_graph_root Digest EqAuthorizationBinding('oracle_tcb_threat_graph_root'); \
30 independent_review_receipt_root Digest \
EqAuthorizationBinding('independent_review_receipt_root'); \
31 revocation_receipt_root Digest EqAuthorizationBinding('revocation_receipt_root'); \
32 semantic_loss_and_component_result_root Digest \
EqAuthorizationBinding('semantic_loss_and_component_result_root'); \
33 redaction_policy_root Digest EqAuthorizationBinding('redaction_policy_root'); \
34 validation_axis_root Digest EqAuthorizationBinding('validation_axis_root'). \
I07GovernedCorpusDischargeEnvelopeBytesV1 is exactly \
FRAME_UTF8('i07-governed-corpus-discharge-envelope-v1')||U16LE(34)||\
concat(U16LE(ordinal)||U64LE(payload_byte_len)||payload in increasing ordinal order); Utf8 payload is \
U64LE(byte_len)||the exact nonempty bytes, U16LE and U64LE payloads are exactly 2 and 8 bytes, and a \
Digest payload is exactly 32 nonzero raw bytes. Every one of the 34 schema constraints is checked \
before accepting the record. The exact envelope length is \
(8+41)+2+34*(2+8)+(6*8+41+3+40+30+30+63)+2*2+8+25*32=1458 bytes; every \
other length is refused. Unknown, optional, duplicate, omitted, reordered or trailing fields and \
unknown field types or constraint kinds are refused. \
I07GovernedCorpusDischargeEnvelopeSchemaDigest=BLAKE3::derive_key(\
'org.frankensim.i07.governed-corpus-discharge-envelope-schema.v1',\
I07GovernedCorpusDischargeEnvelopeSchemaBytesV1), and \
I07GovernedCorpusDischargeEnvelopeDigest=BLAKE3::derive_key(\
'org.frankensim.i07.governed-corpus-discharge-envelope.v1',\
I07GovernedCorpusDischargeEnvelopeBytesV1). The realized FutureArtifact digest and successor-installed \
same-ID External root must both byte-equal I07GovernedCorpusDischargeEnvelopeDigest. \
GovernedOutputSchemaRecordV1=U16LE(1)||U16LE(3)||\
FRAME_BYTES(Utf8('i07.governed-real-ap242-corpus-discharge'))||\
FRAME_BYTES(Utf8('GovernedRealAp242CorpusDischarge/i07-governed-real-ap242-corpus'))||\
raw32(I07GovernedCorpusDischargeEnvelopeSchemaDigest). \
GovernedOutputSchemaRecordV1 is exactly 2*2+(8+8+40)+(8+8+63)+32=171 bytes. \
GovernedOutputSchemaSetRoot=BLAKE3::derive_key(\
'org.frankensim.i07.governed-output-schema-set.v1',raw32(GovernedOutputSchemaMatrixDigest)||\
U64LE(1)||FRAME_BYTES(GovernedOutputSchemaRecordV1)). GovernedOutputSchemaMembershipProofV1 is exactly \
raw32(GovernedOutputSchemaMatrixDigest)||U64LE(0)||U64LE(1)||\
FRAME_BYTES(GovernedOutputSchemaRecordV1); verification recomputes the singleton set root, so a \
proof is exactly 32+2*8+8+171=227 bytes and a different record cannot masquerade as membership. \
GovernanceCommitted freezes the exact matrix bytes, \
matrix digest, schema bytes, schema digest and set root before candidate execution; CandidateFrozen \
copies the same set root; P binds it at tag 0x0014; and the precommit authorization carries the exact \
membership proof for the sole FutureArtifact. Authorization verifies protocol/stage/scope, the derived \
role, exactly-one cardinality, exact-waiver rule, protected/future bijection, all three root equalities, \
schema availability and independent rehash before granting realization. The realized envelope must \
decode under the exact thirty-four-field schema and every EqIntentTag and EqCandidateFrozenComponent \
constraint, both exact digest constraints, and every authorization binding must byte-equal its already \
committed source. A permissive substitute schema, matrix \
extension, role alias or swap, missing or extra output, cross-protocol, cross-stage or cross-scope \
schema, absent or forged membership proof, GovernanceCommitted/CandidateFrozen/P root mismatch, or \
schema/output mutation is IntegrityFailed. Untrusted matrix bytes are preflight-capped at 4096 before \
the exact 221-byte comparison, schema bytes at 4194304, membership-proof bytes at 8192 before the \
exact 227-byte comparison, realized-envelope bytes are exactly 1458, and every Utf8 payload is capped \
at 65536 before its exact literal or schema constraint is checked; \
schema field count is exactly 34 and schema-set count exactly 1. Top-level matrix, schema, set, proof \
and envelope projections have depth 0; each framed child increments depth once, depth 64 is admitted \
and an attempt to enter depth 65 is refused before recursion or allocation. Counts, lengths and \
cumulative sizes use checked arithmetic against remaining input and caps before allocation or hash \
finalization. Before authority use, two independently implemented matrix/schema/set/proof/envelope \
encoders and decoders must match published exact-byte and exact-root KATs and pass schema-byte, \
field-table, matrix-row, role-derivation, membership, permissive-schema, missing/extra-output, \
cross-stage/scope, root-mismatch, exact-cap/cap-plus-one, exact-envelope-length/length-plus-one, \
depth-64/depth-65, overflow, truncation, \
reorder and trailing-byte mutation twins. \
PrecommitAuthorizationSchemaBytes and CommitReceiptSchemaBytes are mandatory content-addressed \
canonical schema artifacts starting respectively with \
FRAME_UTF8('i07-governed-corpus-discharge-precommit-authorization-schema-v1') and \
FRAME_UTF8('i07-governed-corpus-discharge-commit-receipt-schema-v1'); each freezes its complete \
ordinal/name/type/requiredness/constraint-AST table and every nested canonical encoding. \
PrecommitAuthorizationSchemaDigest=H(\
'org.frankensim.i07.precommit-authorization-schema.v1',PrecommitAuthorizationSchemaBytes), and \
CommitReceiptSchemaDigest=H('org.frankensim.i07.commit-receipt-schema.v1',\
CommitReceiptSchemaBytes). ReceiptSchemaSetRoot=H('org.frankensim.i07.receipt-schema-set.v1',\
FRAME_LIST(lexicographically_raw_sorted_unique(\
FRAME_UTF8('PrecommitAuthorization')||raw(PrecommitAuthorizationSchemaDigest),\
FRAME_UTF8('CommitReceipt')||raw(CommitReceiptSchemaDigest)))) and contains exactly those two \
case-sensitive role/schema pairs; arbitrary nonzero schema digests and role aliases are forbidden. A \
decoder recomputes both schema digests and the set root from the exact available schema bytes before \
accepting either receipt; absence or mismatch is IntegrityFailed. Every \
authorization/commit receipt embeds its role, schema digest and ReceiptSchemaSetRoot membership proof; \
a missing, extra, swapped or mismatched role/schema is IntegrityFailed. \
I07_ENVELOPE_DIGEST_PENDING_V1 and \
I07_AMENDMENT_RECORD_PENDING_V1 name semantic future slots only; their encoded union is Pending and \
cannot alias a digest. P may bind the already-existing predecessor-stage receipt digest but never \
serializes the newly created authorization/commit receipt digest, realized envelope digest, final \
successor digest or realized AmendmentRecord digest. The embedded \
I07GovernedCorpusDischargeReceiptV1 precommit authorization binds the predecessor/stage receipt, \
expected old head, ReceiptSchemaSetRoot, GovernedOutputSchemaSetRoot and P, but never the realized \
envelope, final successor or new head. The envelope \
embeds only that authorization; the successor embeds the envelope digest; the AmendmentRecord binds \
predecessor and final successor. A separate non-embedded \
I07GovernedCorpusDischargeCommitReceiptV1 atomically binds P, authorization, realized envelope, \
predecessor/final successor, verified AmendmentRecord, ReceiptSchemaSetRoot, \
GovernedOutputSchemaSetRoot and observed old/new \
authority heads. NewAuthorityHeadDigest=H('org.frankensim.i07.authority-head.v1',\
raw(old_head_digest)||U64LE(old_generation)||raw(TransactionIntentDigest)||\
raw(authorization_digest)||raw(realized_envelope_digest)||raw(final_successor_digest)||\
raw(amendment_record_digest)); NewAuthorityHead is \
NewAuthorityHeadDigest||U64LE(checked_add(old_generation,1)). The fenced transaction succeeds only if \
the current head byte-equals expected_authority_head, and it durably writes exactly the derived new \
head plus the commit receipt in the same atomic commit; CAS failure, generation overflow or a \
different proposed new head performs no waiver retirement, successor publication, capability grant \
or other effect. This strict P -> authorization -> envelope -> successor -> AmendmentRecord -> \
NewAuthorityHead -> commit-receipt DAG prevents a \
content-hash cycle. \
TRANSACTION_INTENT_ENCODING_V1: P starts with the 25 exact ASCII bytes \
I07_TRANSACTION_INTENT_V1 followed by one trailing byte 0x00, then U16LE(20), then exactly twenty \
U16LE(tag)||U64LE(payload_byte_len)||payload fields in strictly increasing order: 0x0001 initiative_id \
Utf8; 0x0002 schema_identity Utf8; 0x0003 successor_version U64LE; 0x0004 predecessor_manifest_digest \
Digest; 0x0005 expected_authority_head AuthorityHead; 0x0006 predecessor_stage_receipt_digest Digest; \
0x0007 candidate_freeze_commitment_digest Digest; 0x0008 retired_waiver_subjects Utf8List; 0x0009 \
protected_bindings ProtectedBindingList; 0x000a coupled_transaction_group Utf8; 0x000b governance_stage \
U16LE where GovernanceCommitted=1,CandidateFrozen=2,RealizationCommitted=3,\
RevealedForAdjudication=4,Closed=5; 0x000c authority_scope Utf8; 0x000d mutation_fence MutationFence; \
0x000e envelope_slots FutureArtifactList; 0x000f final_successor_slot FutureDigest; 0x0010 \
amendment_record_slot FutureDigest; 0x0011 receipt_schema_set_root Digest; 0x0012 profile_registry_root \
Digest; 0x0013 profile_assignment_root Digest; 0x0014 governed_output_schema_set_root Digest. Digest is \
exactly 32 nonzero raw bytes. Utf8 is \
U64LE(byte_len)||exact nonempty case-sensitive bytes without normalization. AuthorityHead is \
Digest||U64LE(generation), exactly 40 bytes. MutationFence is nonzero Digest(idempotency_key)||\
U64LE(attempt_id)||U64LE(capability_epoch), exactly 48 bytes. Utf8List=U64LE(count)||\
concat(FRAME_BYTES(Utf8)) with 1..4096 lexically encoded-element-ordered unique entries here. \
ProtectedBinding=FRAME_BYTES(Utf8(slot_id))||FRAME_BYTES(Utf8(artifact_role))||\
Digest(artifact_schema)||Digest(protected_root). ProtectedBindingList=U64LE(count)||\
concat(FRAME_BYTES(ProtectedBinding)); it has 1..4096 records ordered by the full encoded record, with \
unique slot_id, artifact_role and (slot_id,artifact_role) identities. \
FutureArtifact=FRAME_BYTES(Utf8(slot_id))||FRAME_BYTES(Utf8(related_waiver_subject))||\
FRAME_BYTES(Utf8(artifact_role))||Digest(artifact_schema)||digest_union. FutureArtifactList=U64LE(count)||\
concat(FRAME_BYTES(FutureArtifact)); it has 1..4096 full-record-ordered entries with unique slot, \
waiver-subject and role identities. FutureDigest=FRAME_BYTES(Utf8(slot_id))||digest_union. \
digest_union is exactly byte 0=Pending with no following bytes or byte 1=Digest followed by 32 raw \
bytes; tags 0x000e..0x0010 require Pending. P is at most 16777216 bytes, every Utf8 payload is at most \
65536 bytes, and every U64 count/length is range-checked with checked arithmetic against remaining \
input and these caps before allocation. Every scalar occurs once with exact length; missing, extra, \
duplicate, unknown, out-of-order, empty-required, non-minimal, overflowed or trailing bytes are \
refused. A raw or all-zero digest, role or artifact-schema swap, unavailable/mismatched future artifact \
schema, permissive or unregistered output schema, arbitrary same-ID External fixture, receipt alone, \
waiver removal alone, split successor, post-reveal amendment or structurally frozen successor without \
the verified fs-vvreg transaction grants no discharge or promotion authority. Before discharge, two \
independent encoders must match a published 20-field transaction-intent KAT and \
header/NUL/count/every-tag/type/length/cardinality/order/role-swap/bijection/union/cap/overflow/\
trailing-byte mutation suite. SECURITY: the \
production import/export path performs no network or URI resolution and executes no plugin, macro, \
script, foreign code, or embedded payload; inputs cannot trigger proof execution or unbounded \
decompression. Development-only theorem lanes may invoke only the manifest-pinned sandboxed formal \
checker under a capability-limited supervisor, or consume an independently authenticated checker \
receipt; that checker never receives production network or file authority. \
ACCESSIBILITY_AGENT_PARITY: every verdict and loss is available in deterministic structured form; \
color, plots, proprietary UI, and visual inspection are never the sole carrier. \
DETERMINISM_PROJECTION: the bitwise-comparable logical-event/adjudication projection contains only \
logical ids, deterministic sequence keys, roots, states, bounds and decisions. Raw wall/host \
timestamps, scheduler traces, process ids, durations, percentiles and machine telemetry are retained \
under a separate content root and compared only under preregistered bands; their bytes/root never \
participate in G5 equality. PERFORMANCE: p50/\
p95/p99 wall, peak resident bytes, bytes/entities/edges per second, and cancellation latency bind \
exact model/profile/corpus/machine fingerprints; initially Unmeasured and no speed claim. PROMOTION: \
generated no-orphan/waiver/drift lint, independent manifest adjudication, exact terminal states, \
semantic-loss reconciliation, required falsifier receipts, and DSR proof are all load-bearing; any \
missing/stale field fails closed.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: PROFILE_REGISTRY_FIXTURE,
            source: FixtureSource::AuthoredSpec {
                spec: "I07_PROFILE_REGISTRY_V1. schema_id is the exact 23-byte ASCII string i07-profile-registry-v1. The canonical decoder refuses a missing, extra, duplicate, unknown, out-of-order, case-folded, normalized, or trailing schema/profile field. ProfileKey={part21_encoding_edition_id,\
part21_encoding_artifact_digest,ap242_application_protocol_id,ap242_edition_id,\
ap242_schema_artifact_digest,conformance_class_id,implementation_method_id,\
character_encoding_id,signature_policy_id,external_reference_policy_id,\
supported_entity_inventory_digest,supported_semantic_atom_inventory_digest,\
validation_rule_inventory_digest,allowed_reference_scc_digest,unit_frame_policy_digest,\
canonical_writer_policy_digest,no_claim_policy_digest,coverage_floor_digest,limit_policy_digest}. \
All ids are opaque exact strings plus content roots; no code infers edition/profile from a friendly \
name. The initial synthetic profile admits only the fixture-declared clear-text sections/entities/\
rules and inert external-reference strings. Anchors, digital signatures, XML encodings, foreign \
application protocols, unregistered FILE_SCHEMA values, executable extensions, and network \
resolution are Unsupported unless a successor registry explicitly admits them. Each profile lists \
entity/type/attribute/derived/omitted/select/aggregate/reference validation rules; semantic-atom \
coverage; allowed reference SCC classes; unit/frame policies; resource limits; canonical writer \
policy; exact no-claims; and positive required coverage for every claim/leaf semantic domain. \
ProfileDefinitionSchemaRoot, ProfileKeyDigest, ProfileEntryRoot, ProfileRegistryRoot, \
RequiredAssignmentUniverseRoot and ProfileAssignmentRoot use exactly the framed derive-key algebra \
frozen by I07_CAMPAIGN_POLICY_V1. ProfileEntryRoot, not a friendly key, is the complete expanded \
profile authority identity; every entry envelope embeds its recomputed ProfileDefinitionSchemaRoot, \
ProfileKeyDigest and canonical definition bytes with exact U32LE field count. The assignment is a \
total immutable mapping whose decoded logical-key domain exactly equals RequiredAssignmentUniverseRoot, from \
every required (leaf_id,claim_or_domain_id,AssignmentTargetKind,target_id,ProfileRole) logical key to \
exactly one registered (ProfileKeyDigest,ProfileEntryRoot). Target kinds include Fixture, Case, \
CorpusItem and Stratum; roles include Input, Output, SourceEndpoint, TargetEndpoint, ReferenceOracle, \
MigrationIntermediate and TheoremInstantiation. Multiple required logical keys all execute, and \
duplicate/missing/unexpected logical keys, unknown kinds/roles, repeated entries hidden by ordinal aliases, \
mismatched embedded keys, or a favorable replay-time choice are IntegrityFailed. \
PROFILE_INSTANCE_GATE: this text is a schema target, not a concrete execution registry. Before any \
candidate, a successor freezes the full machine registry AST, exact ProfileDefinitionSchemaBytes/\
ProfileDefinitionSchemaRoot and ProfileKeyDigest/ProfileEntryRoot pairs, ProfileRegistryRoot, \
RequiredAssignmentUniverseRoot, ProfileAssignmentRoot, every inventory/rule/limit/canonicalization byte, \
required claim/leaf set, target/role assignment matrix, nonvacuity matrix, independent decoder, \
canonical registry digest, and independently reproduced KAT_GATE vector. \
--profile <profile-root> is exactly the assigned ProfileEntryRoot, and replay refuses any request \
without the bound schema/registry/universe/assignment membership proofs. Real ISO artifacts \
are not embedded here and enter only through the governed external corpus slot with exact edition, \
license/access, digest, and review receipt.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i07-part21-synthetic-corpus",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 part21-synthetic-corpus. GENERATOR emits bounded clear-text \
exchanges from a grammar with HEADER/DATA/ENDSEC/END-ISO structure, simple and complex entity \
instances, forward/back references, legal SCCs, aggregates, selects, enumerations, omitted '$', \
derived '*', escaped strings, Unicode/code-page declarations allowed by the synthetic profile, \
binary/hex literals, comments, signed exponents, negative zero, and numeric lexical twins. NEGATIVE: \
wrong FILE_SCHEMA/profile, duplicate/dangling/type-invalid refs, illegal complex-entity components, \
unterminated strings/comments, nonfinite/overflowing numeric values, invalid encodings, over-depth/\
over-count/over-byte graphs, forbidden external/anchor/signature sections, and trailing payload. \
GROUND_TRUTH: complete token+span AST, typed instance/reference graph, SCC inventory, refusal rule/span, \
and source root generated without production parser. SEEDS: alias 'i07/part21', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i07-part21-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 part21 CORE STAGE-HOLDOUT. Same exact generator, grammar, \
negative taxonomy, and ground-truth schema as i07-part21-synthetic-corpus with unseen token layouts, \
reference SCCs, complex entities, encodings, and limit boundaries. SEEDS: alias 'i07/part21', \
k=65536..=69631. Public deterministic replay; no IID/secrecy authority. Sole consumer \
i07-parser-writer-core; access before that leaf begins is campaign IntegrityFailed.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i07-product-assembly-configurations",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 product-assembly-configurations. GENERATOR creates 1..64 \
products, versions/definitions/views, 1..512 occurrences, shared definitions, nested assemblies, \
configuration alternatives, serial/date/lot effectivities, external identities, provenance, \
suppression/substitution, and exact rational rigid transforms in mixed SI/source units. ADVERSARIES: \
#renumber/permutation, equal names for distinct occurrences, reused definition in many parents, \
diamond sharing, traversal-order changes, transform-order/handedness/unit twins, collision attempts, \
ambiguous split/merge/rebind, forbidden assembly cycle, and separately legal non-assembly reference \
cycles. GROUND_TRUTH is a typed attributed multigraph plus lineage morphisms and exact refusal set. \
SEEDS: alias 'i07/assembly', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i07-product-assembly-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 product/assembly CORE STAGE-HOLDOUT. Same generator and \
ground-truth graph schema as i07-product-assembly-configurations with unseen sharing, effectivity, \
identity-collision, ambiguity, and transform cases. SEEDS: alias 'i07/assembly', \
k=65536..=69631. Public deterministic replay only. Sole consumer i07-structure-geometry-core.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i07-mixed-exact-tessellated-geometry",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 mixed-exact-tessellated-geometry. CASES combine analytic \
boxes/cylinders/cones/spheres/tori, rational NURBS patches, trimmed multi-shell B-reps, open shells, \
voids, seams, near-tangent faces, slivers, mixed exact+tessellated representations, independent \
tessellations, and representation/validation properties under nested units/frames/occurrences. \
NEGATIVE: flipped shells/normals, broken trim loop, seam gap, self-intersection, duplicate shell, \
lost void, exact/tessellated laundering, falsified validation property, unit/frame transform error, \
and sub-feature below separation premise. ORACLES: exact closed-form values where available and \
independently checked outward-rounded interval/ball enclosures otherwise; coverage-complete interval \
subdivision for bidirectional distance; exact oriented intersections, \
winding/degree and property-specific topology checks. SCORES: every case preregisters L_ref>0 m, \
A_ref=max(A_source,A_target), V_ref=max(V_source,V_target), and \
J_ref=max(||J_source||_F,||J_target||_F), where J is the centroidal density-free volume \
second-moment tensor expressed in the preregistered comparison frame with unit m^5. \
The two directed d terms are independently certified upper endpoints, never sample maxima. \
s_H=max(d_source_to_target,d_target_to_source)/(1e-12 m+1e-8 L_ref); \
s_A=|A_source-A_target|/(1e-18 m^2+1e-8 A_ref); \
s_V=|V_source-V_target|/(1e-24 m^3+1e-8 V_ref); \
s_J=||J_source-J_target||_F/(1e-30 m^5+1e-8 J_ref). Every A,V,J, reference scale, \
absolute difference, maximum, norm, sum and quotient is exact or is evaluated as one dependency-aware \
outward-rounded expression; s_H,s_A,s_V,s_J denote certified finite upper endpoints, never \
point estimates or independently rounded numerator/denominator samples. The pre-candidate applicability mask \
selects only defined properties (for example, no V/J score for an open shell); \
	s_geometry=max of every applicable s_H,s_A,s_V,s_J. GEOMETRY_ACCEPTANCE is the conjunction \
	of finite s_geometry<=1 and Satisfied for every pre-applicable discrete orientation, closure, \
	manifoldness, self-intersection-freedom, validation-property, and topology certificate. Unknown, \
	Indeterminate, Unsupported, missing, or malformed certificate state blocks promotion for that \
	property and can never be averaged into s_geometry. Every scale, \
support, norm, applicability predicate, quadrature/subdivision, reach/separation premise, formula/unit, \
and oracle root freezes before candidate execution. SEEDS: alias 'i07/geometry', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i07-mixed-geometry-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 geometry CORE STAGE-HOLDOUT. Same exact generator, oracle, \
score arithmetic, and topology/refusal policy as i07-mixed-exact-tessellated-geometry with unseen \
trims, seams, mixed representations, unit/frame nests, and negative twins. SEEDS: alias \
'i07/geometry', k=65536..=69631. Public deterministic replay only. Sole consumer \
i07-structure-geometry-core.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i07-semantic-pmi-gdt-texture",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 semantic-pmi-gdt-texture. GENERATOR builds typed semantic \
feature/annotation graphs spanning linear/angular/radial dimensions; plus-minus, limit, and \
compound values; form/orientation/location/runout/profile characteristics; datum features, targets, \
systems and precedence; projected/unequally-disposed/free-state/tangent-plane/material-condition \
modifiers; tolerance zones; surface texture direction/process/limits; requirement/document links; \
and separate presentation-only curves/text. Numeric values use exact rationals/decimals with source \
units and conversion roots. ADVERSARIES: modifier deletion/substitution, datum reorder, target rebind, \
wrong zone frame, degree/radian and inch/mm twins, conflicting duplicates, orphan target, graphical-only \
lookalike, and unsupported enumeration. GROUND_TRUTH: complete typed atom/association graph and exact \
loss dispositions. SEEDS: alias 'i07/pmi', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i07-semantic-pmi-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 semantic-PMI CORE STAGE-HOLDOUT. Same generator, typed \
ground-truth graph, numeric semantics, and adversary taxonomy as i07-semantic-pmi-gdt-texture with \
unseen datum systems, modifiers, zones, target associations, textures, and unit/frame twins. SEEDS: \
alias 'i07/pmi', k=65536..=69631. Public deterministic replay only. Sole consumer \
i07-pmi-material-core.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i07-material-process-lot-references",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 material-process-lot-references. GENERATOR associates \
occurrences/features with named material identities, specifications+editions, temper/heat treatment, \
casting/forging/machining/additive/coating processes, lot/heat/batch, supplier namespace, requirement \
and evidence/document refs, license/access classes, supersession/conflict edges, validity-domain \
stubs, and content digests or explicit Unpinned states. No external content is fetched or invented. \
ADVERSARIES: same marketing grade/different process, same lot string/different namespace, conflicting \
spec editions, missing digest, inaccessible URI, superseded reference, duplicate-key overwrite, \
occurrence/feature scope drift, and fake property inferred from name. GROUND_TRUTH: typed reference/\
lineage graph and loss/refusal inventory. SEEDS: alias 'i07/material', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i07-material-process-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 material/process/lot CORE STAGE-HOLDOUT. Same generator, \
ground-truth graph, no-fetch policy, and adversaries as i07-material-process-lot-references with unseen \
namespace, conflict, supersession, license/access, and scope cases. SEEDS: alias 'i07/material', \
k=65536..=69631. Public deterministic replay only. Sole consumer i07-pmi-material-core.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i07-harness-ewis-graphs",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 harness-ewis-graphs. GENERATOR creates connector occurrences, \
cavities/pins/contacts, wires and multicore cables, splices, branches, shields and shield terminations, \
bundles, route segments, labels, gauges, materials, grounds/bonds, intentionally open pins, multi-drop \
nets, requirement/evidence links, and nested assembly placement. Separate edge colors encode electrical \
connectivity, containment/bundling, route adjacency, shield coverage, ground/bond, and logical signal. \
ADVERSARIES: duplicate labels for distinct nets, connector reuse, pin/cavity swap, lost core, false \
name-based merge, missing splice leg, shield grounded at wrong end, ground-vs-bond conflation, gauge unit \
swap, route-frame inversion, and ambiguous occurrence endpoint. GROUND_TRUTH: typed attributed multigraph \
per edge family and exact loss/refusal set. SEEDS: alias 'i07/harness', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i07-harness-ewis-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 harness/EWIS CORE STAGE-HOLDOUT. Same generator, colored \
ground-truth multigraph, and adversary taxonomy as i07-harness-ewis-graphs with unseen reuse, splice, \
shield, multi-drop, open-pin, ground/bond, and occurrence ambiguity. SEEDS: alias 'i07/harness', \
k=65536..=69631. Public deterministic replay only. Sole consumer i07-harness-loss-core.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i07-semantic-loss-unsupported-corpus",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 semantic-loss-unsupported-corpus. CASES cross every supported \
semantic atom with PreservedExact/PreservedNormalized/SemanticallyMapped/BoundedApproximation \
dispositions and inject unsupported entity, \
attribute, select alternative, enumeration, representation, PMI modifier, material/evidence link, \
harness edge, kinematic construct, external reference, signature-covered region, and vendor extension. \
NEGATIVE receipts omit one source atom, duplicate a DispositionRecord, fabricate an OriginRecord, \
break forward/inverse incidence transposition, collapse a split/merge into false one-to-one lineage, hide an \
approximation bound, relabel presentation PMI semantic, report unsupported as success, or publish empty \
loss after failure. GROUND_TRUTH: independent complete source/target atom inventories, forward+inverse \
incidence maps, counts/roots, exact records/counterpart sets/derivations/reasons/owners/authority \
effects, and expected terminal state. \
SEEDS: alias 'i07/loss', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i07-semantic-loss-core-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 semantic-loss CORE STAGE-HOLDOUT. Same semantic-atom cross, \
DispositionRecord/OriginRecord reconciliation, unsupported-domain taxonomy, independent inventory oracle, and \
negative receipt grammar as i07-semantic-loss-unsupported-corpus, with unseen omission, \
duplicate-record, fabricated-origin, non-transpose-incidence, split/merge-lineage, hidden-bound, \
presentation-as-semantic, failure-with-clean-empty, and \
legitimate-empty-inventory cases. A legitimate empty case must publish both independently reconciled \
roots and zero counts; absence or partial publication never passes. SEEDS: alias 'i07/loss', \
k=65536..=69631. Public deterministic replay only. Sole consumer i07-harness-loss-core.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i07-adversarial-reference-graphs",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 adversarial-reference-graphs. Finite typed graphs contain long \
forward chains, wide fanout, diamond sharing, self/cross references, legal non-assembly SCCs, forbidden \
assembly/configuration cycles, dangling/type-wrong/duplicate refs, external-document identifiers, \
opaque refs crossing supported/unsupported boundaries, repeated labels, hash-prefix collisions, depth/\
count/byte cap boundaries, and canonical-order ties. Each case records exact SCCs, topological quotient, \
allowed/refused rule, stable semantic identities, resource ceiling, and minimized fault. SEEDS: alias \
'i07/reference-graph', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i07-safe-opaque-extensions",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 safe-opaque-extensions. GENERATOR embeds bounded unknown \
entities/attributes with exact raw octets, normalized token AST, nested aggregates/select/function \
tokens, strings/escapes, references within the unknown closure and across supported entities, and \
source context/profile. CASES cover identity #maps, canonical #renumbering with total ref remap, \
sidecar-only preservation, and safe refusal. NEGATIVE: missing ref, ref collision, context/profile drift, \
signature-covered edit, forbidden external/anchor section, URI resolver attempt, executable/plugin \
payload declaration, over-budget blob/depth, invalid tokenization, changed raw byte, and semantic \
authority laundering. GROUND_TRUTH: octet/token/ref/hazard ledgers and expected reinsertion \
disposition. SEEDS: alias 'i07/opaque', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i07-safe-opaque-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 safe-opaque MAX STAGE-HOLDOUT. Same generator, byte/token/ref/\
hazard ground truth, and negative taxonomy as i07-safe-opaque-extensions with unseen nested tokens, \
cross-boundary closures, remaps, signature/context hazards, and cap boundaries. SEEDS: alias \
'i07/opaque', k=131072..=135167. Public deterministic replay only. Sole consumer \
i07-opaque-max.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i07-kinematic-pair-assemblies",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 kinematic-pair-assemblies. GENERATOR builds mechanisms over \
revolute, prismatic, cylindrical, spherical, planar, screw, fixed, gear-ratio, rack-pinion, and declared \
compound pairs within the synthetic supported profile. Every pair records ordered attachment \
occurrences/frames, axis/sense, coordinate zero, angular/linear units, admissible configuration \
relation, tangent/Pfaffian relation, freedoms, open/closed limits, pitch/ratio/coupling signs, branch, \
and occurrence/configuration scope. ADVERSARIES: axis flip, frame order, motor composition order, \
degree/radian and mm/m twins, pitch/ratio inverse/sign, swapped limits, branch alias, sampled-motion \
lookalike with different relation, holonomic/nonholonomic conflation, inferred contact/clearance, and \
ambiguous attachment. GROUND_TRUTH: symbolic relation graph plus interval domain checks. SEEDS: alias \
'i07/kinematics', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i07-kinematic-pair-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 kinematic-pair MAX STAGE-HOLDOUT. Same generator, symbolic \
relation/tangent ground truth, and adversary taxonomy as i07-kinematic-pair-assemblies with unseen pair \
compositions, frame/sign/unit/limit/branch and sampled-lookalike cases. SEEDS: alias \
'i07/kinematics', k=131072..=135167. Public deterministic replay only. Sole consumer \
i07-kinematics-max.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i07-edition-profile-migration-pairs",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 edition-profile-migration-pairs. SYNTHETIC schema/profile pairs \
freeze endpoint roots, common/refined/deprecated/split/merged/unrepresentable atom crosswalks, defaults, \
residual-information carriers, direct and multi-hop mappings, and exact admitted domains/quotients for \
GetPut, PutGet, PutPut and path-coherence laws. Cases span all baseline domains plus opaque/kinematic \
atoms. NEGATIVE: implicit default, removed residual, false inverse after merge, cross-edition leakage, \
lost stable lineage, hidden intermediate loss, direct-vs-multihop incoherence, unpinned endpoint, and \
Unsupported reported migrated. GROUND_TRUTH: independent finite lens tables and exact loss graphs. No \
licensed schema text. SEEDS: alias 'i07/migration', k=0..=4095.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i07-edition-migration-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 edition/profile migration MAX STAGE-HOLDOUT. Same synthetic \
endpoint/crosswalk/lens generator and ground truth as i07-edition-profile-migration-pairs with unseen \
split/merge/default/residual/path-coherence cases across all semantic domains. SEEDS: alias \
'i07/migration', k=131072..=135167. Public deterministic replay only. Sole consumer \
i07-migration-max.",
            },
            partition: Partition::HeldOut,
        },
        FixturePin {
            id: "i07-semantic-equivalence-theorem-card",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_THEOREM_TARGET_V1 supported-roundtrip-equivalence. TARGET only: finite \
typed attributed relational exact source/target semantic categories; proved normalization congruence; \
separate typed directed-error/loss quantaloid whose homs are complete lattices in a frozen authority \
order and whose identity/zero, infinity/refusal, unital associative composition, and arbitrary-join \
preservation in each variable are machine laws, with property-specific geometry components and \
semantic-loss colors; partial/refusal and loss \
objects; import/export exact functors plus enriched comparison maps; object/morphism maps; unit/frame, \
identity/lineage, product/configuration, geometry-budget, PMI, material/evidence, harness, kinematic and \
loss naturality squares; unit/counit or explicit natural isomorphisms on the exact subcategory; enriched \
comparison/composition laws; no finite-tolerance-as-equality rule; runtime premise binding. AXIOM_POLICY \
target i07.lean-axioms.v1 admits exactly propext, Quot.sound, \
Classical.choice and rejects sorryAx, custom postulates, theorem-equivalent axioms, native-oracle \
declarations outside the admitted kernel/checker TCB. FORMAL_PROJECTION_GATE: manifest version 1 is prose and mints no theorem authority. Before \
candidate/proof execution a FrozenManifest::amend successor must freeze canonical proposition AST, \
all symbol/definition bytes+digests, total runtime-premise schema, total deterministic AST-to-Lean \
translation, structural round-trip checks, fully qualified declaration/environment, complete \
transitive axiom-closure algorithm, proof/TCB roots, and independent binding receipt.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i07-semantic-descent-theorem-card",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_THEOREM_TARGET_V1 finite-semantic-sheaf-and-stack-descent. TARGET only: \
finite cover of module/configuration/geometry-PMI/material/harness/kinematic contexts; complete nerve \
with every nonempty finite intersection; typed stalk signatures; functorial restriction maps and \
composition identities; pair/triple/higher-overlap orientations and coherence; abelian coefficient \
complexes and separately typed nonabelian Cech cocycle/torsor/groupoid data; a proved effective \
normalization quotient carrying strict sheaf H0 existence+literal uniqueness; separately typed \
groupoid/stack object-and-morphism descent up to the exact retained equivalence and automorphism \
structure; localized restriction defects; coboundaries with independently checked d^2=0; \
exact closedness and non-exactness witness for any additive pinned-cover Cech H1 obstruction; \
nonabelian descent laws; basis/change-of-cover maps. SHEAF_H1_RATCHET: identification with derived \
sheaf H1 requires a separately bound and kernel-checked Leray/acyclic-cover comparison, cofinal \
direct-limit theorem, or another exact comparison; absent that receipt the conclusion remains \
cover-relative. AXIOM_POLICY \
and FORMAL_PROJECTION_GATE are exactly i07-semantic-equivalence-theorem-card. Version-1 prose mints \
neither H0 nor H1 authority; a pre-candidate amendment must freeze the complete machine proposition, \
definitions, runtime premise schema, translation, declaration/environment, axiom closure, and \
independent binding receipt.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i07-semantic-adversary-microgrammar",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_M0_TARGET_V1 semantic-adversary-microgrammar. TARGET raw grammar has six \
16-way axes for entity/edge structure, reference SCC, semantic domain decoration, unit/frame/identity \
mutation, loss/migration disposition, and descent/theorem stratum: 16^6=N=16,777,216 raw decorated \
records before validity/canonical quotient. It spans every baseline/max domain, positive/negative \
premises, nontrivial SCCs, split/merge residuals, restriction defects, closed exact and closed \
non-exact cochains, geometry boundary cases, and counterexample-minimization moves. SUPERGRAMMAR uses \
alias 'i07/adversary', k=0..=4095 with max 4096 trajectories x 4096 evaluations=16,777,216 \
evaluations and whole-campaign preflight. M0_FORMALIZATION_GATE: version-1 prose has no cardinality, \
enumeration, quotient, or exhaustive authority. Before any candidate a successor must freeze full \
decorated grammar/encoding/domain/validity/stratum/tag/parameter/event ASTs, total enumeration/exclusion \
order, rank/unrank/sharding algorithms, source digests, independent decoder, bijection/cardinality/\
canonical-quotient proofs, nonvacuity floors, cost preflight, and Merkle completeness root.",
            },
            partition: Partition::Development,
        },
        FixturePin {
            id: "i07-semantic-adversary-max-holdout",
            source: FixtureSource::AuthoredSpec {
                spec: "I07_FIXTURE_V1 semantic adversary MAX STAGE-HOLDOUT. Same future formalized \
supergrammar, predicates, nonvacuity floors, and independent membership/minimization protocol as \
i07-semantic-adversary-microgrammar, but alias 'i07/adversary', k=131072..=135167. Public \
deterministic replay only. Sole consumer i07-falsifier-max. No candidate generation is authorized \
until M0_FORMALIZATION_GATE is discharged by amendment.",
            },
            partition: Partition::HeldOut,
        },
    ]
}

#[allow(clippy::too_many_lines)]
fn i07_obligations() -> Vec<ObligationRow> {
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
    const COMMON_EVENTS: &[&str] = &[
        "artifact.published",
        "campaign.admitted",
        "campaign.requested",
        "case.cancel.requested",
        "case.cancel.observed",
        "case.cancelled",
        "case.drain.started",
        "case.drain.completed",
        "case.failed",
        "case.finalized",
        "case.started",
        "checkpoint.committed",
        "claim.adjudicated",
        "failure_bundle.retained",
        "evidence.atomic_publish_committed",
        "integrity.failed",
        "loss.itemized",
        "resume.completed",
    ];
    vec![
        ObligationRow {
            leaf: "i07-parser-writer-core",
            claims_covered: &[
                "i07-part21-profile-bounded-parse",
                "i07-canonical-export-determinism-and-containment",
            ],
            unit_cases: UNIT_CASES,
            g0: "GENERATORS: bounded Part-21 grammar, profile registry, reference graphs, \
and negative/refusal twins. VALIDITY: exact ProfileKey, token/AST/type/ref graph and cap predicates. \
LAWS: full parse-ground-truth equality; span/ref/SCC preservation; wrong-profile/dangling/type/cap \
refusal; canonical parse(write(parse(x))) semantic-root equality; canonical write idempotence; \
transactional no-partial-publication. SHRINKERS: token/section/entity/reference/SCC delta-debugging that \
preserves the failing rule and ProfileKey. REPLAY: exact fixture aliases/ranges in explicits; every \
case binds manifest/profile/generator/oracle/TCB/input roots and terminal exit code",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                PROFILE_REGISTRY_FIXTURE,
                "i07-part21-synthetic-corpus",
                "i07-part21-core-holdout",
                "i07-adversarial-reference-graphs",
                "i07-semantic-loss-unsupported-corpus",
            ],
            g3_relations: &[
                "entity #renumbering and record permutation preserve semantic graph/root while canonical output bytes converge",
                "whitespace/comment/layout changes preserve semantics but remain visible in source provenance",
                "resource-cap tightening can turn success into named Unsupported/Budget refusal but never partial success",
                "numeric lexical twins preserve exact declared value semantics and canonicalize deterministically",
                "equivalent reference-SCC presentation preserves SCC quotient and stable semantic identities",
            ],
            g4_schedule: "For empty/boundary/max graphs inject cancellation and faults after every \
lexer block, entity batch, SCC pass, canonical-sort shard, temporary write, reparse, and pre-publish \
barrier. Every path is request->drain->finalize; no AP242/receipt is public before validation. Commit a \
canonical checkpoint at 1/3 and 2/3, interrupt, corrupt/stale negative twins, resume/fork idempotently, \
and compare terminal state plus retained FailureBundle. TimedOut, BudgetExhausted, InfrastructureFailed, \
IntegrityFailed, Failed, Unknown and Unsupported remain distinct per campaign policy",
            g5_matrix: "ISA families {aarch64-apple, x86_64}; threads {1,2,7,16}; shards \
{1,3,8}; schedules {forward,reverse,work-steal}; mode deterministic; profiles {debug,release}. Bitwise \
canonical bytes/root/receipt/canonical-logical-event-order equality is required only for identical admitted \
ISA+toolchain+profile fingerprints. Cross-ISA compares parsed semantic roots, exact graph/status fields, \
and declared numeric bands; no cross-ISA byte/digest authority until retained two-host evidence lands",
            entry_point: "scripts/e2e/leapfrog/i07_parser_writer.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i07-parser-writer-core dsr quality --tool frankensim",
            obs_events: COMMON_EVENTS,
            replay_command: "scripts/e2e/leapfrog/i07_parser_writer.sh --manifest <manifest-root> --profile <profile-root> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i07-structure-geometry-core",
            claims_covered: &[
                "i07-product-assembly-configuration-lineage",
                "i07-representation-qualified-geometry-roundtrip",
            ],
            unit_cases: UNIT_CASES,
            g0: "GENERATORS: product/configuration/occurrence multigraphs and mixed exact/tessellated \
geometry with positive/negative twins. VALIDITY: typed graph, lineage, frame/unit, representation, \
geometry-oracle and topology-premise predicates. LAWS: #number/permutation invariance; sharing/\
multiplicity/effectivity/lineage isomorphism; rigid-transform covariance; exact/tessellated type \
separation; bidirectional certified distance/integral bands; property-specific topology refusals. \
SHRINKERS: graph edge/node pruning and geometry face/trim/subdivision reduction that retain identity or \
geometry failure. REPLAY binds every score scale/norm/support/oracle root",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                PROFILE_REGISTRY_FIXTURE,
                "i07-product-assembly-configurations",
                "i07-product-assembly-core-holdout",
                "i07-mixed-exact-tessellated-geometry",
                "i07-mixed-geometry-core-holdout",
                "i07-adversarial-reference-graphs",
            ],
            g3_relations: &[
                "source-local #renumbering and traversal reorder leave semantic identities and lineage fixed",
                "global rigid frame change preserves graph semantics and covaries geometry/frame records",
                "unit rescaling preserves physical geometry and quantity kinds while changing source-unit provenance",
                "shared-definition fanout changes only occurrence edges, never duplicates definition identity",
                "refinement/tessellation cannot strengthen exactness or topology evidence and must compose directed bounds",
            ],
            g4_schedule: "Inject faults/cancel after occurrence identity, lineage collision scan, \
frame composition, each interval geometry tile, topology oracle, loss reconciliation, checkpoint, and \
pre-publish. Drain all geometry/oracle tasks, finalize one terminal state, retain unfinished bounds as \
PartialEvidence only, and publish no preserved/topology claim from partial traversal. Resume/fork twice \
from checksummed canonical graph+subdivision checkpoints and require identical outcomes",
            g5_matrix: "ISA {aarch64-apple,x86_64} x threads {1,2,7,16} x graph shards \
{1,4,9} x geometry tiles {1,8,32} x deterministic debug/release. Same fingerprint requires bitwise \
stable semantic ids, graph roots, exact predicates, interval endpoints, receipt ordering and terminal \
state; cross-ISA is exact graph equality plus certified geometric bands, not unearned digest equality",
            entry_point: "scripts/e2e/leapfrog/i07_structure_geometry.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i07-structure-geometry-core dsr quality --tool frankensim",
            obs_events: &[
                "artifact.published",
                "campaign.admitted",
                "campaign.requested",
                "case.cancel.requested",
                "case.cancel.observed",
                "case.cancelled",
                "case.drain.started",
                "case.drain.completed",
                "case.failed",
                "case.finalized",
                "case.started",
                "checkpoint.committed",
                "claim.adjudicated",
                "failure_bundle.retained",
                "evidence.atomic_publish_committed",
                "geometry.bound.enclosed",
                "identity.bound",
                "integrity.failed",
                "lineage.refused",
                "loss.itemized",
                "resume.completed",
            ],
            replay_command: "scripts/e2e/leapfrog/i07_structure_geometry.sh --manifest <manifest-root> --profile <profile-root> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i07-pmi-material-core",
            claims_covered: &[
                "i07-semantic-pmi-gdt-datum-texture",
                "i07-material-process-lot-requirement-evidence",
            ],
            unit_cases: UNIT_CASES,
            g0: "GENERATORS: semantic PMI/datum/texture graphs and material/process/lot/evidence \
reference graphs. VALIDITY: quantity-kind/unit/frame/target association, identity namespace, \
revision/digest/license/access/scope and no-fetch predicates. LAWS: exact typed graph isomorphism; \
datum precedence and target stability; semantic/presentation separation; dimensional covariance; \
reference/lineage/conflict preservation; no name-based property inference. SHRINKERS: prune one \
modifier/datum/target/reference/conflict edge while preserving the mismatch. REPLAY uses declared \
Philox ranges and exact generator roots",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                PROFILE_REGISTRY_FIXTURE,
                "i07-semantic-pmi-gdt-texture",
                "i07-semantic-pmi-core-holdout",
                "i07-material-process-lot-references",
                "i07-material-process-core-holdout",
                "i07-semantic-loss-unsupported-corpus",
            ],
            g3_relations: &[
                "unit rescaling and frame change preserve physical PMI while provenance/canonical representation changes",
                "presentation-only annotation deletion cannot change semantic-PMI verdict and cannot fill a semantic gap",
                "datum/target id relabeling with a full graph isomorphism preserves semantics; reorder without precedence mapping does not",
                "material/reference input order never chooses a winner among conflicting claims",
                "unavailable external content stays unavailable under URI/string normalization",
            ],
            g4_schedule: "Cancel/fault after each PMI atom family, association pass, dimension check, \
material namespace/reference edge, conflict audit, loss row, checkpoint and pre-publish. External refs \
remain inert under all failures. Drain, finalize, retain redacted FailureBundle; no orphan target or \
partial material graph may publish as Preserved. Resume/fork canonical atom/reference checkpoints and \
require identical graphs and terminal states",
            g5_matrix: "ISA {aarch64-apple,x86_64} x threads {1,2,7,16} x atom shards \
{1,3,11} x deterministic debug/release. Same fingerprint requires bitwise semantic/reference roots, \
loss rows, refusals, canonical logical-event/adjudication projections and exit states; raw timing/\
performance telemetry is separately rooted and banded. Cross-ISA requires exact rational/unit/id graph equality and \
declared numeric conversion bands only",
            entry_point: "scripts/e2e/leapfrog/i07_pmi_material.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i07-pmi-material-core dsr quality --tool frankensim",
            obs_events: &[
                "artifact.published",
                "campaign.admitted",
                "campaign.requested",
                "case.cancel.requested",
                "case.cancel.observed",
                "case.cancelled",
                "case.drain.started",
                "case.drain.completed",
                "case.failed",
                "case.finalized",
                "case.started",
                "checkpoint.committed",
                "claim.adjudicated",
                "failure_bundle.retained",
                "evidence.atomic_publish_committed",
                "integrity.failed",
                "loss.itemized",
                "material.reference.bound",
                "pmi.association.checked",
                "resume.completed",
            ],
            replay_command: "scripts/e2e/leapfrog/i07_pmi_material.sh --manifest <manifest-root> --profile <profile-root> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i07-harness-loss-core",
            claims_covered: &[
                "i07-harness-ewis-connectivity-semantics",
                "i07-exhaustive-semantic-loss-receipt",
            ],
            unit_cases: UNIT_CASES,
            g0: "GENERATORS: colored HarnessGraph fixtures, semantic-loss/unsupported corpus, and \
reference adversaries. VALIDITY: edge-family typing, endpoint/occurrence/unit/frame, source/target atom \
inventory and one-DispositionRecord/one-OriginRecord reconciliation. LAWS: multigraph isomorphism without label \
merging; connectivity/containment/route/shield/ground separation; exact-transpose forward/inverse \
incidence with typed zero/one/many split/merge/synthesis/drop lineage; \
failure never yields empty loss; unsupported remains Unsupported. SHRINKERS: preserve the false/missing \
harness edge or unreconciled atom. REPLAY binds profile/atom inventory/generator/oracle roots",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                PROFILE_REGISTRY_FIXTURE,
                "i07-harness-ewis-graphs",
                "i07-harness-ewis-core-holdout",
                "i07-semantic-loss-unsupported-corpus",
                "i07-semantic-loss-core-holdout",
                "i07-adversarial-reference-graphs",
            ],
            g3_relations: &[
                "connector/wire/net relabeling under typed graph isomorphism preserves semantics",
                "permuting source/target atom traversal preserves the canonical loss receipt",
                "adding presentation-only data cannot discharge a semantic loss",
                "splitting a bundle changes containment only and cannot invent electrical connectivity",
                "unknown construct insertion adds exactly one unsupported/opaque/refused inventory branch and never disappears",
            ],
            g4_schedule: "Fault/cancel after every harness edge-color pass, endpoint resolver, \
source/target walker, reconciliation shard, inverse-map audit, checkpoint and pre-publish. Drain/finalize \
before atomic receipt publication; on any interrupted/incomplete inventory publish only FailureBundle \
with ClaimState=Unknown and EvidenceCompleteness=PartialEvidence, never a clean loss list. Resume/fork independently checks \
canonical inventories and reconciliation roots",
            g5_matrix: "ISA {aarch64-apple,x86_64} x threads {1,2,7,16} x graph/inventory \
shards {1,5,13} x deterministic debug/release. Same fingerprint requires bitwise graph/atom/loss roots, \
row order, terminal states and redacted canonical logical-event fields; raw timing/performance \
telemetry is separately rooted and banded. Cross-ISA requires exact discrete graphs/states and \
declared route-coordinate bands only",
            entry_point: "scripts/e2e/leapfrog/i07_harness_loss.sh",
            tier: CampaignTier::Core,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i07-harness-loss-core dsr quality --tool frankensim",
            obs_events: &[
                "artifact.published",
                "campaign.admitted",
                "campaign.requested",
                "case.cancel.requested",
                "case.cancel.observed",
                "case.cancelled",
                "case.drain.started",
                "case.drain.completed",
                "case.failed",
                "case.finalized",
                "case.started",
                "checkpoint.committed",
                "claim.adjudicated",
                "failure_bundle.retained",
                "evidence.atomic_publish_committed",
                "harness.edge.checked",
                "integrity.failed",
                "loss.inventory.reconciled",
                "loss.itemized",
                "resume.completed",
            ],
            replay_command: "scripts/e2e/leapfrog/i07_harness_loss.sh --manifest <manifest-root> --profile <profile-root> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i07-opaque-max",
            claims_covered: &["i07-safe-opaque-extension-capsules"],
            unit_cases: UNIT_CASES,
            g0: "GENERATORS: bounded reinsertion-admissible/reinsertion-refused opaque token/reference/context capsules and cap \
twins. VALIDITY: exact raw octets, token AST, source context/profile, complete ref closure, injective \
remap, inert-data and signature policies. LAWS: round-trip octet custody; token/ref equivalence; total \
remap or refusal; sidecar-vs-reinsert distinction; no semantic authority; no side effects. SHRINKERS: \
minimize token/ref/context/hazard while retaining the fault. REPLAY binds byte/token/ref roots and \
deny-all side-effect trace",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                PROFILE_REGISTRY_FIXTURE,
                "i07-safe-opaque-extensions",
                "i07-safe-opaque-max-holdout",
                "i07-adversarial-reference-graphs",
                "i07-semantic-loss-unsupported-corpus",
            ],
            g3_relations: &[
                "entity #renumbering preserves opaque token semantics only through the independently checked total ref map",
                "whitespace normalization never changes retained raw-octet custody and changes reinsertion disposition when byte exactness is requested",
                "adding an unresolved dependency monotonically weakens reinsertion to refusal/sidecar",
                "context/profile/signature change cannot preserve an earlier safety admission",
                "opaque payload permutation outside reference closure leaves other capsule identities fixed",
            ],
            g4_schedule: "Cancel/fault at byte capture, tokenization, ref walk, hazard scan, remap, \
reinsertion, independent reparse and pre-publish; inject network/URI/plugin attempts and require deny \
events. Drain/finalize and retain exact redacted capsule root plus FailureBundle; never publish partial \
capsule as safe. Checkpoint/resume/fork at token and remap boundaries with corruption/stale-profile \
negative twins",
            g5_matrix: "ISA {aarch64-apple,x86_64} x threads {1,2,7,16} x capsule shards \
{1,4,16} x deterministic debug/release. Same fingerprint requires bitwise raw/token/ref/hazard/remap/\
receipt roots and terminal states; cross-ISA requires exact byte/token/ref equality, not a blanket \
writer digest claim",
            entry_point: "scripts/e2e/leapfrog/i07_opaque.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i07-opaque-max dsr quality --tool frankensim",
            obs_events: &[
                "artifact.published",
                "campaign.admitted",
                "campaign.requested",
                "case.cancel.requested",
                "case.cancel.observed",
                "case.cancelled",
                "case.drain.started",
                "case.drain.completed",
                "case.failed",
                "case.finalized",
                "case.started",
                "checkpoint.committed",
                "claim.adjudicated",
                "failure_bundle.retained",
                "evidence.atomic_publish_committed",
                "integrity.failed",
                "loss.itemized",
                "opaque.capsule.retained",
                "opaque.reinsertion.refused",
                "resume.completed",
                "side_effect.denied",
            ],
            replay_command: "scripts/e2e/leapfrog/i07_opaque.sh --manifest <manifest-root> --profile <profile-root> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i07-kinematics-max",
            claims_covered: &["i07-kinematic-pair-semantic-roundtrip"],
            unit_cases: UNIT_CASES,
            g0: "GENERATORS: supported pair/joint and compound mechanisms with relation/tangent \
ground truth and adversarial lookalikes. VALIDITY: occurrence/frame/axis/unit/sense/limit/branch plus \
holonomic-vs-Pfaffian relation predicates. LAWS: exact relation-graph equivalence; motor/frame \
covariance; tangent agreement on regular domains; sign/unit/limit preservation; unsupported contact/\
clearance/dynamics refusal. SHRINKERS: remove links/pairs/constraints/parameters while retaining the \
semantic mismatch. REPLAY binds symbolic relation and interval-domain roots",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                PROFILE_REGISTRY_FIXTURE,
                "i07-product-assembly-configurations",
                "i07-kinematic-pair-assemblies",
                "i07-kinematic-pair-max-holdout",
                "i07-semantic-loss-unsupported-corpus",
            ],
            g3_relations: &[
                "global rigid frame conjugation preserves the kinematic relation and covaries attachments",
                "axis reversal paired with coordinate/sense transformation preserves semantics only under the exact declared map",
                "unit rescaling preserves admissible configurations and tangent distribution",
                "joint-name and #number permutations cannot change relation identity",
                "sample-equivalent trajectories do not imply relation equivalence and must be separated by adversarial probes",
            ],
            g4_schedule: "Cancel/fault after occurrence binding, motor/frame normalization, symbolic \
relation construction, tangent/Pfaffian check, interval domain, loss row and pre-publish. Drain every \
relation/oracle task, then finalize; partial sampling can only be PartialEvidence. Resume/fork from canonical graph+\
domain checkpoints and reproduce relation/refusal/terminal roots",
            g5_matrix: "ISA {aarch64-apple,x86_64} x threads {1,2,7,16} x mechanism \
shards {1,3,9} x deterministic debug/release. Same fingerprint requires bitwise pair graph, symbolic \
relation, exact rank/incidence, interval endpoints, loss rows and states; cross-ISA remains within \
declared transcendental/frame bands",
            entry_point: "scripts/e2e/leapfrog/i07_kinematics.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i07-kinematics-max dsr quality --tool frankensim",
            obs_events: &[
                "artifact.published",
                "campaign.admitted",
                "campaign.requested",
                "case.cancel.requested",
                "case.cancel.observed",
                "case.cancelled",
                "case.drain.started",
                "case.drain.completed",
                "case.failed",
                "case.finalized",
                "case.started",
                "checkpoint.committed",
                "claim.adjudicated",
                "failure_bundle.retained",
                "evidence.atomic_publish_committed",
                "integrity.failed",
                "kinematic.relation.checked",
                "loss.itemized",
                "resume.completed",
            ],
            replay_command: "scripts/e2e/leapfrog/i07_kinematics.sh --manifest <manifest-root> --profile <profile-root> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i07-migration-max",
            claims_covered: &["i07-partial-bidirectional-edition-migration"],
            unit_cases: UNIT_CASES,
            g0: "GENERATORS: synthetic endpoint schemas/profiles, finite partial lenses, direct/\
multi-hop paths and split/merge/default/residual negatives. VALIDITY: exact endpoint roots, atom \
crosswalk, admitted law domains/quotients, residual information and lineage predicates. LAWS: GetPut, \
PutGet, PutPut only on their bound domains; source/target/loss reconciliation; identity/lineage \
preservation; direct/path coherence; Unsupported/refusal outside domain. SHRINKERS: minimize mapping \
table/atom/path while retaining a law/coherence fault. REPLAY binds all endpoint/mapping roots",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                PROFILE_REGISTRY_FIXTURE,
                "i07-edition-profile-migration-pairs",
                "i07-edition-migration-max-holdout",
                "i07-semantic-loss-unsupported-corpus",
            ],
            g3_relations: &[
                "identity migration is neutral only for byte-identical profile/schema roots and policy",
                "adding an unrepresentable atom monotonically adds loss/refusal and cannot improve lens authority",
                "direct and multi-hop routes agree only modulo the frozen coherence relation and composed loss",
                "renaming edition labels without root equality changes endpoint identity and refuses",
                "mapping-table input order and shard order preserve canonical lens/receipt roots",
            ],
            g4_schedule: "Cancel/fault after endpoint admission, each crosswalk/lens stage, residual \
capture, domain proof, path composition, loss reconciliation, checkpoint and pre-publish. Drain and \
finalize with no half-migrated target. Resume/fork each hop; corrupt mapping/root/checkpoint twins are \
IntegrityFailed and cannot fall back to a guessed edition",
            g5_matrix: "ISA {aarch64-apple,x86_64} x threads {1,2,7,16} x mapping \
shards {1,4,12} x path orders {direct,forward-composed,balanced-composed} x deterministic debug/release. \
Same fingerprint requires bitwise lens/loss/lineage roots and states; cross-ISA exact discrete equality \
plus declared geometry bands",
            entry_point: "scripts/e2e/leapfrog/i07_migration.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i07-migration-max dsr quality --tool frankensim",
            obs_events: &[
                "artifact.published",
                "campaign.admitted",
                "campaign.requested",
                "case.cancel.requested",
                "case.cancel.observed",
                "case.cancelled",
                "case.drain.started",
                "case.drain.completed",
                "case.failed",
                "case.finalized",
                "case.started",
                "checkpoint.committed",
                "claim.adjudicated",
                "failure_bundle.retained",
                "evidence.atomic_publish_committed",
                "integrity.failed",
                "loss.itemized",
                "migration.lens.checked",
                "migration.path.reconciled",
                "resume.completed",
            ],
            replay_command: "scripts/e2e/leapfrog/i07_migration.sh --manifest <manifest-root> --profile <profile-root> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i07-governed-real-corpus-validation-max",
            claims_covered: &["i07-governed-real-ap242-corpus-validation"],
            unit_cases: UNIT_CASES,
            g0: "GENERATORS: governed corpus package/semantic-role assignment/custody/expected-inventory metadata plus arbitrary-digest, same-ID-fixture-only, receipt-only, waiver-only, split-transaction and negative twins that mutate one root, edition, item, stratum, role, vendor/version, expected atom, loss disposition, band, signature or axis. VALIDITY: exact license/access/custody, corpus completeness, item-addressable ProfileRegistryRoot/ProfileAssignmentRoot membership, immutable pre-access expectations, independent-reader/tool identity, public-vs-governed-blind authority, typed I07GovernedCorpusDischargeEnvelopeV1 and receipt, atomic amendment/authority-head transition, redaction and no-favorable-missingness predicates. LAWS: corpus/item order and worker sharding preserve roots; every required logical assignment and matrix cell executes exactly once; duplicate entries under ordinal aliases refuse; ScopedClauseConformance requires its clause pack, ScopedRequiredMatrixValidated requires every declared cell, and GovernedBlindScopedValidation requires untouched-population joint commitment; separately serialized StandardAxis, InteroperabilityAxis and IndustrialCorpusAxis never infer one another; generic manifest syntax or partial discharge grants no authority; missing/inaccessible evidence cannot pass. SHRINKERS: preserve the first custody, assignment, semantic, reconciliation, omission, discharge, independence or axis-collapse defect. REPLAY binds the predecessor/successor manifests, transaction intent, authority head, discharge envelope/receipt, claim/checker, registry/assignment, corpus, expectation, toolchain, output and adjudication roots without exposing restricted bytes",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                PROFILE_REGISTRY_FIXTURE,
                REAL_CORPUS_SLOT,
            ],
            g3_relations: &[
                "permuting corpus items, declarations, independently owned expected inventories or worker shards preserves the canonical coverage and component-verdict roots",
                "adding a required corpus item, profile assignment, vendor/version cell or predicate cannot improve the all-required verdict until that component passes",
                "renaming friendly editions/vendors/tools without identical bound roots changes identity and refuses rather than preserving a result",
                "redaction that preserves approved content roots and adjudication fields preserves verdicts while unauthorized byte disclosure fails custody",
                "a synthetic success or one scoped interoperability success cannot strengthen any separately serialized governed-corpus axis",
                "replacing the unresolved deck with an arbitrary same-ID External digest or removing its waiver without the typed envelope, verified receipt and one atomic authority-head transaction cannot preserve discharge authority",
            ],
            g4_schedule: "Preflight access and whole-campaign bounds, then execute at most 256 corpus items or 4096 semantic atoms per cancellable batch and inject fault/cancel at discharge-envelope/receipt verification, atomic amendment and authority-head commit, custody capability issue, semantic-role assignment proof, every independent reader/import/export heartbeat, semantic inventory, loss reconciliation, component signature, axis adjudication, checkpoint and pre-publish. Every external process must support supervised terminate-drain-finalize with request observation <=1 s or admission refuses. Drain every granted capability and worker. After the primary cancellation no new capability issue/use, protected read, reveal, adjudication, publication or authority-head advance may begin; only an operation whose own durable commit point preceded that cancellation remains committed, and every in-flight precommit reconciles conservatively without performing a new effect. Finalize exactly one terminal state, retain the redacted authenticated transaction/FailureBundle plus component receipts, and publish no partial favorable axis. Resume/fork only from canonical signed custody-free checkpoints; stale roots, missing items, partial discharge and duplicate attempts are IntegrityFailed",
            g5_matrix: "Corpus order {canonical,reverse,seeded-permutation} x vendors/tools {declared matrix} x ISA {aarch64-apple,x86_64} x workers {1,2,7,16} x shards {1,4,17} x deterministic debug/release. Same fingerprint requires bitwise coverage, profile assignment, semantic/loss/component/axis, canonical logical-event/adjudication projection and terminal roots. Raw timestamps, durations, scheduler/process telemetry and performance roots are retained separately and compared under bands, never for bitwise equality. Cross-implementation and Cross-ISA comparisons use exact discrete receipts plus preregistered certified numeric bands and never assume byte identity",
            entry_point: "scripts/e2e/leapfrog/i07_governed_real_corpus.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i07-governed-real-corpus-validation-max dsr quality --tool frankensim",
            obs_events: &[
                "artifact.published",
                "campaign.admitted",
                "campaign.requested",
                "case.cancel.requested",
                "case.cancel.observed",
                "case.cancelled",
                "case.drain.started",
                "case.drain.completed",
                "case.failed",
                "case.finalized",
                "case.started",
                "checkpoint.committed",
                "claim.adjudicated",
                "corpus.custody.checked",
                "failure_bundle.retained",
                "evidence.atomic_publish_committed",
                "industrial_corpus.axis.adjudicated",
                "integrity.failed",
                "interoperability.axis.adjudicated",
                "loss.itemized",
                "profile.assignment.checked",
                "resume.completed",
                "standard.axis.adjudicated",
            ],
            replay_command: "scripts/e2e/leapfrog/i07_governed_real_corpus.sh --manifest <successor-manifest-root> --profile <profile-root> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i07-theorem-descent-max",
            claims_covered: &[
                "i07-proof-carrying-bidirectional-semantic-equivalence",
                "i07-semantic-sheaf-descent-obstruction-theorem",
            ],
            unit_cases: UNIT_CASES,
            g0: "GENERATORS: theorem-card target AST fixtures after formalization, finite semantic \
categories, functors/naturality squares, covers/stalks/restrictions, exact cochains and positive/negative \
premise bindings. VALIDITY: complete machine proposition/definition/translation/declaration/axiom \
closure plus runtime-premise maps. LAWS: AST serialization and Lean translation round trips; exact \
axiom policy; exact equivalence unit/counit/naturality; directed-error/loss quantaloid complete-lattice, \
identity, monotonicity, associativity, arbitrary-join preservation, enclosure and enriched-comparison \
                laws; complete-nerve restriction composition and higher-overlap coherence; strict \
                H0 equalizer/descent only after an effective quotient; groupoid/stack object-and-\
                morphism descent with retained automorphisms; localized restriction defects; \
	abelian d^2=0 plus pinned-cover Cech H1 closure/non-exactness; separately bound \
	Leray/acyclic or direct-limit comparison before any derived-sheaf-H1 authority; nonabelian Cech \
	cocycle/torsor/groupoid descent laws. \
SHRINKERS: preserve proof/binding/naturality/cochain fault. REPLAY: bind \
the successor manifest, proposition/definition/translation/declaration/environment/axiom-policy, proof, \
runtime-premise, fixture, oracle, and TCB roots. VERSION-1 TARGET ONLY: no candidate may run before \
amendment discharges formal-projection gates",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                PROFILE_REGISTRY_FIXTURE,
                "i07-semantic-equivalence-theorem-card",
                "i07-semantic-descent-theorem-card",
                "i07-product-assembly-configurations",
                "i07-semantic-loss-unsupported-corpus",
            ],
            g3_relations: &[
                "category/object/record presentation reorder preserves the elaborated bound proposition",
                "semantic-id relabeling under a certified isomorphism transports naturality and descent witnesses",
                "cover refinement preserves the theorem only through an explicit full-nerve comparison/refinement map",
                "strict sheaf uniqueness and groupoid/stack uniqueness up to equivalence remain distinct conclusions and transport retained automorphisms",
                "basis change transports cochains contragrediently and preserves exact closure/non-exactness verdict",
                "removing one required semantic naturality square or runtime premise blocks theorem consumption",
                "adding an opaque/unsupported/loss atom outside the supported subcategory cannot be erased by quotient normalization",
            ],
            g4_schedule: "After formalization only, cancel/fault at AST decode, definition/root bind, \
translation, elaboration, transitive axiom walk, proof-kernel check, runtime-premise check, finite cochain \
oracle, and pre-publish. Drain, then finalize, proof/checker subprocess-equivalent tasks without foreign production \
dependency, retain proof/TCB/FailureBundle roots, and never publish theorem authority from partial \
elaboration. Checkpoint/resume/fork exact immutable theorem environments; stale/mismatched environments \
are IntegrityFailed",
            g5_matrix: "Pinned theorem environment/toolchain/axiom policy x host ISA \
{aarch64-apple,x86_64} x checker workers {1,2,7} x deterministic orderings {source,reverse,hash}. The \
fully qualified declaration, elaborated type, proof term, transitive axiom set, finite matrix/cochain \
results, binding roots and terminal state must be identical under the same environment. Cross-host \
authority requires independently retained proof-kernel and binding receipts, not assumed digest parity",
            entry_point: "scripts/e2e/leapfrog/i07_theorems.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i07-theorem-descent-max dsr quality --tool frankensim",
            obs_events: &[
                "artifact.published",
                "campaign.admitted",
                "campaign.requested",
                "case.cancel.requested",
                "case.cancel.observed",
                "case.cancelled",
                "case.drain.started",
                "case.drain.completed",
                "case.failed",
                "case.finalized",
                "case.started",
                "checkpoint.committed",
                "claim.adjudicated",
                "failure_bundle.retained",
                "evidence.atomic_publish_committed",
                "integrity.failed",
                "loss.itemized",
                "proof.axioms.checked",
                "proof.binding.checked",
                "proof.kernel.checked",
                "resume.completed",
                "sheaf.obstruction.checked",
            ],
            replay_command: "scripts/e2e/leapfrog/i07_theorems.sh --manifest <successor-manifest-root> --profile <profile-root> --replay <artifact-id>",
        },
        ObligationRow {
            leaf: "i07-falsifier-max",
            claims_covered: &["i07-semantic-equivalence-counterexample-search"],
            unit_cases: UNIT_CASES,
            g0: "GENERATORS: formalized M0 rank/unrank exhaustive microgrammar and independent \
stage-heldout supergrammar after amendment. VALIDITY: canonical decode/encode/rank/unrank bijection; \
profile/theorem/lens domain membership; nonvacuity floors; shard/root completeness; candidate \
classification. LAWS: every raw rank visited exactly once; canonical quotient coverage; all floors met; \
out-of-domain never refutes; independently reproduced exact countermodel binds immutable roots; \
specification/checker/TCB defects invalidate their own artifacts. SHRINKERS: preserve membership and \
counterexample while minimizing entity/edge/domain/parameter structure. REPLAY: bind successor manifest, \
formal grammar/decoder/bijection/cardinality/canonical-quotient/nonvacuity/preflight/Merkle, profile, \
theorem, migration, oracle, and TCB roots plus exact completed shard set. VERSION-1 TARGET ONLY",
            decks: &[
                CAMPAIGN_POLICY_FIXTURE,
                PROFILE_REGISTRY_FIXTURE,
                "i07-semantic-adversary-microgrammar",
                "i07-semantic-adversary-max-holdout",
                "i07-semantic-equivalence-theorem-card",
                "i07-semantic-descent-theorem-card",
                "i07-edition-profile-migration-pairs",
            ],
            g3_relations: &[
                "rank/unrank and shard permutation preserve the exact candidate multiset and Merkle root",
                "external candidate-id relabeling cannot affect validity, membership, minimization or classification",
                "canonical semantic isomorphs collapse only under the frozen decorated-object quotient",
                "adding a premise violation changes classification to out-of-domain and cannot manufacture refutation",
                "minimization is monotone in the frozen cost order and retains exact immutable claim/profile roots",
                "microgrammar and supergrammar outcomes remain separate coverage axes",
            ],
            g4_schedule: "After M0 formalization only, preflight the entire 16,777,216-record \
microgrammar and every 4096x4096 supergrammar campaign budget before start. Cancel/fault at rank \
boundaries, shard commits, decode/validity/premise/geometry/proof checks, candidate retention, minimizer, \
Merkle merge, checkpoint and pre-publish barrier. Drain all workers, then finalize: retain completed \
shard roots and every \
candidate, but incomplete campaign state has EvidenceCompleteness=PartialEvidence and \
ClaimState=Unknown and never zero-counterexample \
success. Resume/fork exact disjoint shards; duplicate/missing/corrupt shards are IntegrityFailed",
            g5_matrix: "ISA {aarch64-apple,x86_64} x workers {1,2,7,16} x shards \
{1,16,257,4096} x traversal {ascending,descending,bit-reversed} x deterministic debug/release. Same \
fingerprint requires identical raw/canonical counts, floor counts, candidate bytes/classes/minima, shard/\
Merkle roots and terminal states. Cross-host authority requires independent decoder/enumerator/root \
receipts and exact discrete equality; floating geometry checks use certified bands",
            entry_point: "scripts/e2e/leapfrog/i07_falsifier.sh",
            tier: CampaignTier::Max,
            dsr_lane: "env FRANKENSIM_VMANIFEST_LEAF=i07-falsifier-max dsr quality --tool frankensim",
            obs_events: &[
                "artifact.published",
                "campaign.admitted",
                "campaign.requested",
                "candidate.classified",
                "candidate.retained",
                "case.cancel.requested",
                "case.cancel.observed",
                "case.cancelled",
                "case.drain.started",
                "case.drain.completed",
                "case.failed",
                "case.finalized",
                "case.started",
                "checkpoint.committed",
                "claim.adjudicated",
                "counterexample.minimized",
                "failure_bundle.retained",
                "evidence.atomic_publish_committed",
                "integrity.failed",
                "loss.itemized",
                "microgrammar.shard.committed",
                "nonvacuity.checked",
                "resume.completed",
            ],
            replay_command: "scripts/e2e/leapfrog/i07_falsifier.sh --manifest <successor-manifest-root> --profile <profile-root> --replay <artifact-id>",
        },
    ]
}

fn i07_waivers() -> Vec<Waiver> {
    vec![Waiver {
        subject: REAL_CORPUS_SLOT,
        reason: "real ISO 10303-21/AP242 schema/profile artifacts and public/licensed/vendor \
                 manufacturing corpora cannot be invented or embedded here; all authored fixtures \
                 are synthetic AP242-shaped semantic decks with no licensed standard text",
        owner: "I07 manifest/corpus owner frankensim-leapfrog-2026-program-i94v.2.2.8.1 and the named G2/G5 independent-reproduction adjudicators",
        predicate: "fs-vvreg issues a typed I07GovernedCorpusDischargeReceiptV1 for the exact waiver \
                    subject and predicate, and one atomic FrozenManifest::amend authority transaction \
                    installs its same-ID I07GovernedCorpusDischargeEnvelopeV1 External root, removes \
                    this Waiver row, verifies the AmendmentRecord and advances the authority head. \
                    The envelope/receipt pin every exact standard/application-protocol/Part-21 \
                    edition, schema/profile/conformance-class artifact digest, license/access/custody, \
                    item-addressable semantic-role profile assignment, corpus item and expectation \
                    root, QoI/band, standard clause pack or interoperability matrix when claimed, \
                    public-versus-governed-blind authority, joint candidate/checker/expectation \
                    commitment for untouched-population authority, oracle/TCB threat graph, \
                    independent review and revocation receipts, predecessor manifest and exact \
                    transaction intent without placing forbidden text in the production dependency \
                    graph. The verified AmendmentRecord separately binds the final successor and the \
                    atomic transaction advances the authority head. A raw/all-zero digest, fixture-only replacement, \
                    receipt-only change, waiver-only removal or split transaction is not discharge",
        expiry: "before any protected corpus byte, aggregate, derived label or expectation is accessed \
                 and before activation of i07-governed-real-ap242-corpus-validation or its sole \
                 i07-governed-real-corpus-validation-max leaf; review again for every corpus item or \
                 expectation root, standard/application-protocol/Part-21 edition, schema/profile/\
                 conformance class, license/access/custody, semantic-role assignment, QoI/band, \
                 clause pack, interoperability matrix, candidate/checker/toolchain/decision, \
                 oracle/TCB, reviewer/revocation, envelope/receipt, transaction-intent or authority-\
                 head change",
        promotion_effect: "while live, only i07-governed-real-ap242-corpus-validation and \
                           i07-governed-real-corpus-validation-max have NoPromotionAuthority. \
                           Synthetic decks retain their separately scoped code/schema/profile \
                           mechanics authority, but no synthetic result may be relabeled as \
                           real-corpus, scoped-standard, vendor-interoperability, physical-\
                           manufacturing, supplier, drawing, safety, or regulatory validation",
    }]
}
