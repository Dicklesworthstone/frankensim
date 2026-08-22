//! Integration battery for the mechanism-registry contract binding and the
//! retained registry log (bead
//! frankensim-euler-disc-emergent-flagship-t6314.1.2, remaining slices).
//!
//! - Contract-side binding: a declared `HypothesisSource` (id + locator)
//!   joins its version-pinned registry row, failing closed on unregistered
//!   citation names and on locator drift (a locator edit IS a new source
//!   declaration).
//! - Retained log: schema-versioned, deterministic, bounded,
//!   content-addressed projection of the full registry state; any semantic
//!   row or card edit moves the identity.

use fs_euler_disc_e2e::mechanism_registry::{
    ExponentFamily, MechanismCard, MechanismClass, MechanismRegistry, RegisteredSource, SourceType,
};

fn campaign() -> MechanismRegistry {
    fs_euler_disc_e2e::mechanism_registry::campaign_registry().expect("shipped rows register")
}

const MOFFATT_ID: &str = "moffatt-2000-nature";

fn pinned_row(id: &str) -> RegisteredSource {
    RegisteredSource {
        id: id.into(),
        source_type: SourceType::ExperimentalStudy,
        version_or_date: "2020".into(),
        locator: "Journal of Disc Studies 1, 1-2 (2020)".into(),
        license: "cited".into(),
        model_equations: Vec::new(),
        assumptions: vec!["assumption".into()],
        measured_configurations: Vec::new(),
        qois: vec!["qoi".into()],
        reported_uncertainty: "unreported".into(),
        transfer_limitations: vec!["apparatus-specific".into()],
    }
}

fn citing_card(id: &str, source_id: &str) -> MechanismCard {
    MechanismCard {
        id: id.into(),
        class: MechanismClass::WedgeThinGapFlow,
        observables: vec!["terminal exponent".into()],
        distinguishing_interventions: vec!["vacuum run".into()],
        exponent_family: Some(ExponentFamily { lo: 0.3, hi: 0.5 }),
        source_ids: vec![source_id.into()],
    }
}

#[test]
fn binding_joins_a_declared_source_to_its_pinned_row() {
    let registry = campaign();
    let locator = registry
        .source(MOFFATT_ID)
        .expect("campaign row")
        .locator
        .clone();
    let binding = registry
        .bind_contract_source(MOFFATT_ID, &locator)
        .expect("exact declaration binds");
    assert_eq!(binding.source_id(), MOFFATT_ID);
    assert_eq!(binding.locator(), locator);
    assert_eq!(
        binding.row_identity(),
        registry.source(MOFFATT_ID).expect("row").identity(),
        "the binding carries the rich row identity, not just the name"
    );
    assert!(
        binding.citing_cards() >= 1,
        "the Moffatt gap-flow card cites this source"
    );
}

#[test]
fn an_unregistered_citation_name_refuses_to_bind() {
    let registry = campaign();
    let error = registry
        .bind_contract_source("not-a-registered-row", "anything")
        .expect_err("names alone never authorize");
    assert_eq!(error.rule, "euler-registry-unresolved-source");
}

#[test]
fn locator_drift_refuses_as_a_new_declaration() {
    let registry = campaign();
    let pinned = registry.source(MOFFATT_ID).expect("row").locator.clone();
    let drifted = format!("{pinned} [retrieved 2026-08-22]");
    let error = registry
        .bind_contract_source(MOFFATT_ID, &drifted)
        .expect_err("a locator edit is a new declaration");
    assert_eq!(error.rule, "euler-registry-locator-drift");
}

#[test]
fn duplicate_declarations_refuse_and_order_is_preserved() {
    let registry = campaign();
    let moffatt = registry.source(MOFFATT_ID).expect("row").locator.clone();
    let engh_id = "van-den-engh-2000-nature";
    let engh = registry.source(engh_id).expect("row").locator.clone();

    let duplicated = registry.bind_contract_sources(&[
        (MOFFATT_ID, moffatt.as_str()),
        (MOFFATT_ID, moffatt.as_str()),
    ]);
    assert_eq!(
        duplicated.expect_err("duplicates refuse").rule,
        "euler-registry-duplicate-binding"
    );

    let bound = registry
        .bind_contract_sources(&[(engh_id, engh.as_str()), (MOFFATT_ID, moffatt.as_str())])
        .expect("distinct declarations bind");
    let ids: Vec<_> = bound
        .iter()
        .map(fs_euler_disc_e2e::mechanism_registry::SourceBinding::source_id)
        .collect();
    assert_eq!(ids, vec![engh_id, MOFFATT_ID], "declaration order is kept");
}

#[test]
fn an_empty_declaration_list_binds_nothing() {
    let registry = campaign();
    assert!(
        registry
            .bind_contract_sources(&[])
            .expect("vacuous")
            .is_empty()
    );
}

#[test]
fn registry_log_is_deterministic_bounded_and_content_addressed() {
    let first = campaign().retained_log();
    let replay = campaign().retained_log();
    assert_eq!(
        first.canonical_bytes(),
        replay.canonical_bytes(),
        "rebuilds must be byte-identical"
    );
    assert_eq!(first.identity(), replay.identity());
    assert_eq!(first.source_count(), 7);
    assert_eq!(first.card_count(), 8);
    assert!(
        first.canonical_bytes().len() < 65_536,
        "the retained log stays bounded"
    );
    let reproduction = first.reproduction_command();
    assert!(reproduction.starts_with("cargo test -p fs-euler-disc-e2e"));
    assert!(!reproduction.contains("/Users") && !reproduction.contains("/home"));
}

#[test]
fn a_semantic_row_edit_moves_the_log_identity() {
    let base_registry = {
        let mut registry = MechanismRegistry::new();
        registry.register_source(pinned_row("src-a")).expect("row");
        registry
            .register_card(citing_card("card-a", "src-a"))
            .expect("card");
        registry
    };
    let edited_registry = {
        let mut registry = MechanismRegistry::new();
        let mut row = pinned_row("src-a");
        row.assumptions = vec!["assumption".into(), "refined assumption".into()];
        registry.register_source(row).expect("row");
        registry
            .register_card(citing_card("card-a", "src-a"))
            .expect("card");
        registry
    };

    let base = base_registry.retained_log();
    let edited = edited_registry.retained_log();
    assert_ne!(
        base.identity(),
        edited.identity(),
        "an assumption edit is a new declaration and must move the log"
    );
    // The edit is attributable to the source row: the binding identity moves.
    let locator = "Journal of Disc Studies 1, 1-2 (2020)";
    assert_ne!(
        base_registry
            .bind_contract_source("src-a", locator)
            .expect("binds")
            .row_identity(),
        edited_registry
            .bind_contract_source("src-a", locator)
            .expect("binds")
            .row_identity()
    );
}
