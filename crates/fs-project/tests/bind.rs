//! Battery for validate-time material/interface binding resolution
//! (bead f85xj.6.4): the reference project binds fully with retained,
//! replayable receipts; envelope-outside-domain refuses BEFORE any
//! solve; Unstated uncertainty warns up front without refusing; and a
//! card with coexisting conflicting claims resolves only through the
//! explicit claim pin recorded in the project file — auto-pick is
//! impossible by construction.
#![allow(clippy::float_cmp)] // bit-exact propagation of stored card values is the contract under test

use std::collections::BTreeMap;

use fs_blake3::hash_bytes;
use fs_evidence::{Color, ValidityDomain};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    PINNED_CLAIM_POLICY_TAG, PropertyClaim, PropertyKey, PropertyValue, Provenance, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_project::{
    BindingRequirements, CONTACT_RESISTANCE_DIMS, CONTACT_RESISTANCE_PROPERTY, CardLibrary,
    EntityDecl, Envelope, InterfaceCardBinding, MaterialBinding, MaterialResolution, ProjectSpec,
    TEMPERATURE_AXIS, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY, ThermalLimit,
    resolve_bindings,
};
use fs_qty::QtyAny;
use fs_regime::{
    AxisViolationKind, EnvelopeCoverage, OperatingPoint, QoiClaim, audit_product_output_with_cards,
};
use fs_scenario::Violation;

const KELVIN: fs_qty::Dims = fs_qty::Dims([0, 0, 0, 1, 0, 0]);
const PASCAL: fs_qty::Dims = fs_qty::Dims([-1, 1, -2, 0, 0, 0]);

fn kelvin(value: f64) -> QtyAny {
    QtyAny {
        value,
        dims: KELVIN,
    }
}

fn provenance(source: &str) -> Provenance {
    Provenance {
        source: source.to_string(),
        license: "internal-use".to_string(),
        artifact: None,
    }
}

fn conductivity_claim(
    value: f64,
    lo: f64,
    hi: f64,
    source: &str,
    uncertainty: UncertaintyModel,
) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS),
        value: PropertyValue::Scalar {
            value,
            dims: THERMAL_CONDUCTIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained().with(TEMPERATURE_AXIS, lo, hi),
        uncertainty,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: provenance(source),
    }
}

fn resistance_claim(value: f64, lo: f64, hi: f64, source: &str) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(CONTACT_RESISTANCE_PROPERTY, CONTACT_RESISTANCE_DIMS),
        value: PropertyValue::Scalar {
            value,
            dims: CONTACT_RESISTANCE_DIMS,
        },
        validity: ValidityDomain::unconstrained().with(TEMPERATURE_AXIS, lo, hi),
        uncertainty: UncertaintyModel::RelativeHalfWidth {
            fraction: 0.2,
            confidence: 0.9,
        },
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: provenance(source),
    }
}

fn stated() -> UncertaintyModel {
    UncertaintyModel::HalfWidth {
        half_width: 0.02,
        confidence: 0.95,
    }
}

fn fr4_state() -> MaterialStateId {
    MaterialStateId {
        chemistry: "FR4".to_string(),
        phase: "laminate".to_string(),
        process: "cured".to_string(),
        revision: 0,
    }
}

fn copper_state() -> MaterialStateId {
    MaterialStateId {
        chemistry: "Cu-OFE".to_string(),
        phase: "wrought".to_string(),
        process: "annealed".to_string(),
        revision: 0,
    }
}

fn material_card(state: MaterialStateId, claims: Vec<PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).expect("fixture claim inserts");
    }
    MaterialCard::assemble(state, set, Vec::new()).expect("fixture card assembles")
}

fn interface_card(claims: Vec<PropertyClaim>) -> InterfaceSystemCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).expect("fixture claim inserts");
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: copper_state(),
            texture_frame: "lapped-frame-1".to_string(),
        },
        SurfaceSpec {
            material: fr4_state(),
            texture_frame: "as-cured-frame-2".to_string(),
        },
        SystemContext {
            medium: "tim-paste-x".to_string(),
            third_body: None,
            environment: "air".to_string(),
            history: "cured-once".to_string(),
        },
        set,
        Vec::new(),
    )
    .expect("fixture interface card assembles")
}

/// The reference library: FR4 board, copper spreader, one TIM system,
/// every claim valid over [200, 450] K.
fn reference_library() -> (CardLibrary, String, String, String) {
    let mut library = CardLibrary::new();
    let board = library.insert_material(material_card(
        fr4_state(),
        vec![conductivity_claim(
            0.3,
            200.0,
            450.0,
            "laminate handbook",
            stated(),
        )],
    ));
    let spreader = library.insert_material(material_card(
        copper_state(),
        vec![conductivity_claim(
            390.0,
            200.0,
            450.0,
            "copper handbook",
            stated(),
        )],
    ));
    let tim = library.insert_interface(interface_card(vec![resistance_claim(
        2.0e-5,
        200.0,
        450.0,
        "tim datasheet",
    )]));
    (library, board, spreader, tim)
}

/// The reference project skeleton for resolution: two regions, one
/// interface, envelope 263.15..318.15 K, one 371.15 K junction limit on
/// the spreader region.
fn reference_spec(board_card: &str, spreader_card: &str, tim_card: &str) -> ProjectSpec {
    ProjectSpec {
        assembly: Some(vec![
            EntityDecl::Assembly {
                name: "rig".to_string(),
                display: "Reference rig".to_string(),
                expect_id: None,
            },
            EntityDecl::Part {
                parent: "rig".to_string(),
                name: "stack".to_string(),
                display: "Board stack".to_string(),
                expect_id: None,
            },
            EntityDecl::Region {
                parent: "stack".to_string(),
                name: "board".to_string(),
                display: "PCB".to_string(),
                expect_id: None,
            },
            EntityDecl::Region {
                parent: "stack".to_string(),
                name: "spreader".to_string(),
                display: "Heat spreader".to_string(),
                expect_id: None,
            },
            EntityDecl::Interface {
                parent: "rig".to_string(),
                name: "tim".to_string(),
                display: "Spreader TIM".to_string(),
                from: "board".to_string(),
                to: "spreader".to_string(),
                expect_id: None,
            },
        ]),
        envelope: Some(Envelope {
            ambient_lo: kelvin(263.15),
            ambient_hi: kelvin(318.15),
            pressure: QtyAny {
                value: 101_325.0,
                dims: PASCAL,
            },
        }),
        requirements: Some(vec![ThermalLimit {
            class: "junction".to_string(),
            region: "spreader".to_string(),
            limit: kelvin(371.15),
            margin: kelvin(10.0),
        }]),
        materials: Some(vec![
            MaterialBinding {
                region: "board".to_string(),
                card: board_card.to_string(),
                claim: None,
                state: "FR4/laminate/cured rev 0".to_string(),
                temp_lo: kelvin(233.15),
                temp_hi: kelvin(398.15),
                source: "seed-v1".to_string(),
            },
            MaterialBinding {
                region: "spreader".to_string(),
                card: spreader_card.to_string(),
                claim: None,
                state: "Cu-OFE/wrought/annealed rev 0".to_string(),
                temp_lo: kelvin(233.15),
                temp_hi: kelvin(398.15),
                source: "seed-v1".to_string(),
            },
        ]),
        interface_cards: Some(vec![InterfaceCardBinding {
            interface: "tim".to_string(),
            card: tim_card.to_string(),
            claim: None,
            source: "seed-v1".to_string(),
        }]),
        ..ProjectSpec::default()
    }
}

fn assert_code<'v>(
    resolution: &'v MaterialResolution,
    code: &str,
    what_contains: &str,
    fix_contains: &str,
) -> &'v Violation {
    let found = resolution
        .violations
        .iter()
        .find(|violation| violation.code == code)
        .unwrap_or_else(|| {
            panic!(
                "expected violation `{code}`; got: {:?}",
                resolution
                    .violations
                    .iter()
                    .map(|violation| violation.code)
                    .collect::<Vec<_>>()
            )
        });
    assert!(
        found.what.contains(what_contains),
        "violation `{code}` what {:?} must name {what_contains:?}",
        found.what
    );
    assert!(
        found.fix.contains(fix_contains),
        "violation `{code}` fix {:?} must offer {fix_contains:?}",
        found.fix
    );
    found
}

#[test]
fn the_reference_project_binds_fully_with_replayable_receipts() {
    let (library, board, spreader, tim) = reference_library();
    let spec = reference_spec(&board, &spreader, &tim);
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());

    assert!(
        resolution.admissible(),
        "reference project must bind: {:?}",
        resolution.violations
    );
    assert!(resolution.advisories.is_empty());
    assert_eq!(resolution.bindings.len(), 3);
    let receipts: Vec<_> = resolution.receipts().collect();
    assert_eq!(receipts.len(), 6, "3 bindings x 1 property x 2 endpoints");

    // Every retained receipt replays against the claim set that made it,
    // and the retained bytes decode back to the retained receipt.
    for retained in &receipts {
        let claims = if retained.context.starts_with("interface") {
            library.interface(&tim).expect("tim card").claims()
        } else if retained.context.contains("`board`") {
            library.material(&board).expect("board card").claims()
        } else {
            library.material(&spreader).expect("spreader card").claims()
        };
        claims
            .verify_receipt(&retained.receipt)
            .expect("retained receipt replays");
        let hash = fs_matdb::PropertyUsageReceipt::from_bytes_verified(&retained.bytes, {
            retained.receipt.content_hash()
        })
        .expect("retained bytes decode");
        assert_eq!(&hash, &retained.receipt);
    }

    // The logged table carries the complete chain per row.
    let table = resolution.render_table();
    for needle in [
        "region `board`",
        "region `spreader`",
        "interface `tim`",
        THERMAL_CONDUCTIVITY_PROPERTY,
        CONTACT_RESISTANCE_PROPERTY,
        "FR4/laminate/cured rev 0",
        "Cu-OFE/wrought/annealed rev 0",
        "tim-paste-x",
        "laminate handbook",
        "seed-v1",
        "±0.02 @ confidence 0.95",
    ] {
        assert!(
            table.contains(needle),
            "table must log {needle:?}:\n{table}"
        );
    }
}

#[test]
fn actual_matdb_claim_demotes_after_reference_project_runtime_drift() {
    let (library, board, spreader, tim) = reference_library();
    let spec = reference_spec(&board, &spreader, &tim);
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());
    assert!(
        resolution.admissible(),
        "reference project must bind: {:?}",
        resolution.violations
    );

    let board_binding = resolution
        .bindings
        .iter()
        .find(|binding| binding.target == fs_project::BindingTarget::Region("board".to_string()))
        .expect("reference board binding");
    let conductivity = board_binding
        .properties
        .iter()
        .find(|property| property.property == THERMAL_CONDUCTIVITY_PROPERTY)
        .expect("board conductivity");
    let registry = resolution.regime_audit_cards();
    assert_eq!(registry.len(), 3, "board, spreader, and TIM claims");
    assert_eq!(
        conductivity.regime_card.name,
        format!(
            "fs-matdb:claim:{}:{}",
            board_binding.card, conductivity.selected_claim
        )
    );
    assert_eq!(
        conductivity.regime_card.version,
        format!(
            "property-usage-receipt-v{}",
            conductivity.receipt_lo.receipt.schema_version
        )
    );
    assert_eq!(
        conductivity.regime_card.validity.bound(TEMPERATURE_AXIS),
        Some((200.0, 450.0))
    );

    let audit = audit_product_output_with_cards(
        &registry,
        &[
            OperatingPoint {
                id: "reference".to_string(),
                groups: BTreeMap::from([(TEMPERATURE_AXIS.to_string(), 371.15)]),
            },
            OperatingPoint {
                id: "runtime-drift".to_string(),
                groups: BTreeMap::from([(TEMPERATURE_AXIS.to_string(), 451.0)]),
            },
        ],
        &[QoiClaim {
            qoi: "board-temperature".to_string(),
            color: Color::Estimated {
                estimator: "thermal-reference-project".to_string(),
                dispersion: 0.02,
            },
            model_cards: vec![conductivity.regime_card.name.clone()],
            override_acknowledgement: None,
        }],
    )
    .expect("actual matdb card audits");
    let receipt = &audit.receipts[0];
    assert_eq!(receipt.coverage, EnvelopeCoverage::Partial);
    assert_eq!(receipt.in_domain_points, ["reference"]);
    assert_eq!(receipt.out_of_domain_points, ["runtime-drift"]);
    assert_eq!(receipt.model_cards.len(), 1);
    assert_eq!(receipt.model_cards[0].name, conductivity.regime_card.name);
    assert_eq!(
        receipt.model_cards[0].version,
        conductivity.regime_card.version
    );
    assert_eq!(receipt.violations.len(), 1);
    let violation = &receipt.violations[0];
    assert_eq!(violation.point, "runtime-drift");
    assert_eq!(violation.card, conductivity.regime_card.name);
    assert_eq!(violation.axis, TEMPERATURE_AXIS);
    assert_eq!(violation.observed, Some(451.0));
    assert_eq!(violation.lo, 200.0);
    assert_eq!(violation.hi, 450.0);
    assert_eq!(violation.kind, AxisViolationKind::Above);
    assert!(violation.distance > 0.0);
    assert!(matches!(
        receipt.effective_color,
        Color::Estimated { dispersion, .. } if dispersion.is_infinite()
    ));
    let canonical = receipt.to_canonical_json();
    assert!(canonical.contains(&board_binding.card));
    assert!(canonical.contains(&conductivity.selected_claim));
}

#[test]
fn envelope_outside_the_card_domain_refuses_at_validate() {
    let (mut library, board, _, tim) = reference_library();
    // A spreader card whose data ends at 350 K, below the 371.15 K limit.
    let narrow = library.insert_material(material_card(
        copper_state(),
        vec![conductivity_claim(
            390.0,
            250.0,
            350.0,
            "narrow source",
            stated(),
        )],
    ));
    let spec = reference_spec(&board, &narrow, &tim);
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());

    assert!(!resolution.admissible());
    assert_code(
        &resolution,
        "project-binding-domain-uncovered",
        "region `spreader`",
        "never extrapolates",
    );
}

#[test]
fn an_admitted_range_narrower_than_the_envelope_refuses() {
    let (library, board, spreader, tim) = reference_library();
    let mut spec = reference_spec(&board, &spreader, &tim);
    // The user admits the board material only above the envelope floor.
    spec.materials.as_mut().expect("materials")[0].temp_lo = kelvin(283.15);
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());

    let found = assert_code(
        &resolution,
        "project-material-envelope-uncovered",
        "region `board`",
        "qualified for this envelope",
    );
    // The message states both ranges so the fix is decidable from it.
    assert!(found.what.contains("263.15") && found.what.contains("283.15"));
}

#[test]
fn unstated_uncertainty_warns_up_front_without_refusing() {
    let (mut library, _, spreader, tim) = reference_library();
    let unstated = library.insert_material(material_card(
        fr4_state(),
        vec![conductivity_claim(
            0.3,
            200.0,
            450.0,
            "uncertainty-free datasheet",
            UncertaintyModel::Unstated,
        )],
    ));
    let spec = reference_spec(&unstated, &spreader, &tim);
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());

    assert!(
        resolution.admissible(),
        "Unstated uncertainty warns, never refuses: {:?}",
        resolution.violations
    );
    let advisory = resolution
        .advisories
        .iter()
        .find(|advisory| advisory.code == "binding-uncertainty-unstated")
        .expect("unstated advisory");
    assert!(advisory.what.contains("region `board`"));
    assert!(advisory.note.contains("Estimated"));
    let board_row = &resolution.bindings[0];
    assert!(board_row.properties[0].unstated_uncertainty);
    assert!(resolution.render_table().contains("UNSTATED"));
}

#[test]
fn conflicting_claims_refuse_without_a_pin_and_resolve_only_through_one() {
    let (mut library, _, spreader, tim) = reference_library();
    let handbook = conductivity_claim(0.29, 200.0, 450.0, "handbook", stated());
    let vendor = conductivity_claim(0.35, 200.0, 450.0, "vendor datasheet", stated());
    let vendor_pin = vendor.content_hash().to_hex();
    let conflicted = library.insert_material(material_card(fr4_state(), vec![handbook, vendor]));

    // Without a pin: a typed refusal naming every candidate, and NO
    // resolved row for the region — auto-pick is impossible.
    let spec = reference_spec(&conflicted, &spreader, &tim);
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());
    assert!(!resolution.admissible());
    let conflict = assert_code(
        &resolution,
        "project-binding-claims-conflict",
        "region `board`",
        ":claim <hex>",
    );
    assert!(
        conflict.fix.contains(&vendor_pin),
        "fix lists the candidates"
    );
    assert!(
        !resolution
            .bindings
            .iter()
            .any(|binding| binding.target == fs_project::BindingTarget::Region("board".into())),
        "no value may be produced for a conflicted, unpinned binding"
    );

    // With the pin recorded in the project file: resolves to exactly the
    // pinned claim, through the receipted pinned-claim policy.
    let mut pinned_spec = reference_spec(&conflicted, &spreader, &tim);
    pinned_spec.materials.as_mut().expect("materials")[0].claim = Some(vendor_pin.clone());
    let pinned = resolve_bindings(
        &pinned_spec,
        &library,
        &BindingRequirements::thermal_steady_v1(),
    );
    assert!(pinned.admissible(), "{:?}", pinned.violations);
    let board_row = &pinned.bindings[0];
    assert_eq!(board_row.pinned_claim.as_deref(), Some(vendor_pin.as_str()));
    let property = &board_row.properties[0];
    assert_eq!(property.selected_claim, vendor_pin);
    assert_eq!(property.value_lo, 0.35);
    assert_eq!(property.receipt_lo.receipt.policy, PINNED_CLAIM_POLICY_TAG);
    assert_eq!(property.provenance_source, "vendor datasheet");
}

#[test]
fn pin_refusals_are_typed_and_never_bypass_the_domain() {
    let (mut library, _, spreader, tim) = reference_library();
    let wide = conductivity_claim(0.29, 200.0, 450.0, "handbook", stated());
    let narrow = conductivity_claim(0.35, 250.0, 300.0, "narrow vendor", stated());
    let narrow_pin = narrow.content_hash().to_hex();
    let card = library.insert_material(material_card(fr4_state(), vec![wide, narrow]));

    // A pin that names no claim on the card.
    let mut spec = reference_spec(&card, &spreader, &tim);
    spec.materials.as_mut().expect("materials")[0].claim =
        Some(hash_bytes(b"not a claim").to_hex());
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());
    assert_code(
        &resolution,
        "project-binding-pin-unknown",
        "region `board`",
        "stale",
    );

    // A pinned claim that does not cover the admitted range: the wide
    // claim covers it, but a pin never silently substitutes or extends.
    let mut spec = reference_spec(&card, &spreader, &tim);
    spec.materials.as_mut().expect("materials")[0].claim = Some(narrow_pin);
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());
    assert_code(
        &resolution,
        "project-binding-pin-domain",
        "region `board`",
        "never bypasses the extrapolation refusal",
    );

    // A malformed pin.
    let mut spec = reference_spec(&card, &spreader, &tim);
    spec.materials.as_mut().expect("materials")[0].claim = Some("zz".repeat(32));
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());
    assert_code(
        &resolution,
        "project-binding-pin-malformed",
        "region `board`",
        "full content hash",
    );
}

#[test]
fn a_range_stitched_from_two_claims_is_not_a_resolution() {
    let (mut library, _, spreader, tim) = reference_library();
    let cold = conductivity_claim(0.29, 200.0, 300.0, "cold source", stated());
    let hot = conductivity_claim(0.31, 300.0, 450.0, "hot source", stated());
    let stitched = library.insert_material(material_card(fr4_state(), vec![cold, hot]));
    let spec = reference_spec(&stitched, &spreader, &tim);
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());

    assert_code(
        &resolution,
        "project-binding-domain-split",
        "region `board`",
        "whole admitted range",
    );
}

#[test]
fn structural_refusals_cover_cards_states_targets_and_coverage() {
    let (mut library, board, spreader, tim) = reference_library();

    // Unknown card hash.
    let mut spec = reference_spec(&board, &spreader, &tim);
    spec.materials.as_mut().expect("materials")[0].card = "ab".repeat(32);
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());
    assert_code(
        &resolution,
        "project-material-card-unknown",
        "region `board`",
        "correct the card hash",
    );

    // Manufactured-state mismatch.
    let mut spec = reference_spec(&board, &spreader, &tim);
    spec.materials.as_mut().expect("materials")[0].state = "FR4".to_string();
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());
    assert_code(
        &resolution,
        "project-material-state-mismatch",
        "`FR4`",
        "manufactured state, never a name",
    );

    // Binding to a non-region entity.
    let mut spec = reference_spec(&board, &spreader, &tim);
    spec.materials.as_mut().expect("materials")[0].region = "stack".to_string();
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());
    assert_code(
        &resolution,
        "project-material-target-kind",
        "`stack`",
        "region entities",
    );
    // ... which also leaves `board` unbound.
    assert_code(
        &resolution,
        "project-material-unbound-region",
        "`board`",
        "exactly one material card",
    );

    // Duplicate binding.
    let mut spec = reference_spec(&board, &spreader, &tim);
    let duplicate = spec.materials.as_ref().expect("materials")[0].clone();
    spec.materials.as_mut().expect("materials").push(duplicate);
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());
    assert_code(
        &resolution,
        "project-material-binding-duplicate",
        "`board`",
        "exactly once",
    );

    // Unbound interface.
    let mut spec = reference_spec(&board, &spreader, &tim);
    spec.interface_cards = Some(Vec::new());
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());
    assert_code(
        &resolution,
        "project-interface-unbound",
        "`tim`",
        "exactly one TIM/contact system card",
    );

    // A card claim under the right name with the wrong dimensions.
    let wrong_dims = library.insert_material(material_card(
        fr4_state(),
        vec![PropertyClaim {
            key: PropertyKey::new(THERMAL_CONDUCTIVITY_PROPERTY, KELVIN),
            value: PropertyValue::Scalar {
                value: 0.3,
                dims: KELVIN,
            },
            validity: ValidityDomain::unconstrained().with(TEMPERATURE_AXIS, 200.0, 450.0),
            uncertainty: stated(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: provenance("dimensionally confused source"),
        }],
    ));
    let spec = reference_spec(&wrong_dims, &spreader, &tim);
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());
    assert_code(
        &resolution,
        "project-binding-property-dims",
        "region `board`",
        "fix the card data",
    );

    // A card without the required property at all.
    let empty = library.insert_material(material_card(fr4_state(), Vec::new()));
    let spec = reference_spec(&empty, &spreader, &tim);
    let resolution = resolve_bindings(&spec, &library, &BindingRequirements::thermal_steady_v1());
    assert_code(
        &resolution,
        "project-binding-property-missing",
        "region `board`",
        "bind a card that states this property",
    );
}

#[test]
fn missing_sections_name_the_precondition_instead_of_panicking() {
    let resolution = resolve_bindings(
        &ProjectSpec::default(),
        &CardLibrary::new(),
        &BindingRequirements::thermal_steady_v1(),
    );
    assert_code(
        &resolution,
        "project-binding-preconditions",
        "material resolution needs",
        "structural validation",
    );
}

/// The product-level constants must track `fs-conduction`'s consumption:
/// a drift here would validate against one property spelling and solve
/// against another.
#[test]
fn resolution_constants_track_fs_conduction() {
    assert_eq!(
        THERMAL_CONDUCTIVITY_DIMS,
        fs_conduction::material::CONDUCTIVITY_DIMS
    );
    assert_eq!(TEMPERATURE_AXIS, fs_conduction::material::TEMPERATURE_AXIS);
    assert_eq!(
        CONTACT_RESISTANCE_PROPERTY,
        fs_conduction::interface::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY
    );
    assert_eq!(
        CONTACT_RESISTANCE_DIMS,
        fs_conduction::interface::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS
    );
}
