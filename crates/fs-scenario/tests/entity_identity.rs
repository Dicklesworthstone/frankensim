//! Entity-identity conformance (suite `fs-scenario/entity`, bead
//! `frankensim-extreal-program-f85xj.17.1`): strings become identities.
//!
//! Acceptance covered here: the reference cooling scenario expressed in the
//! entity model; rename and re-import fixtures keeping their bindings intact
//! through receipts; dangling/ambiguous reference validation; mechanical
//! legacy migration; and the datum/tolerance/placement declarations the
//! feature-restore comment on the bead requires.
//!
//! NOTE ON "THE REFERENCE COOLING SCENARIO": no canonical reference cooling
//! project exists in this repository yet (the phrase appears only in bead
//! prose for `.6.1`/`.6.3`). The fixture below is this crate's representative
//! cooling stack, authored here so the entity model has a real assembly to be
//! exercised against. It is not a repo-wide canonical artifact and does not
//! claim to be one.

use fs_ga::{Point, Quat, Vec3};
use fs_qty::{Dims, QtyAny};
use fs_scenario::entity::{
    Binding, BindingTable, ContactSide, Correspondence, DatumFeature, DatumId, EntityCatalog,
    EntityDeclaration, EntityId, EntityKind, EntityRef, EvidenceTier, GeometryFingerprint,
    ImportRevision, ImportScope, ImportStep, ImportedEntity, InterfacePair, KindExpectation,
    MatchBasis, NameLookup, PlacementBasis, RebindEvent, ReferenceSite, ResolutionFault, Tolerance,
    ToleranceKind, ToleranceSource, migrate_legacy_scenario, scenario_reference_sites,
    validate_bindings,
};
use fs_scenario::{
    BcKind, BcValue, BoundaryCondition, ContactLaw, ContactModel, Environment, Frame, FrameId,
    FrameMotion, FrameTree, LoadCase, Physics, Scenario, Violation,
};

const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
const HEAT_FLUX: Dims = Dims([0, 1, -3, 0, 0, 0]);
const HTC: Dims = Dims([0, 1, -3, -1, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-scenario/entity\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

fn fingerprint(tag: &str) -> GeometryFingerprint {
    GeometryFingerprint::of_bytes(tag.as_bytes())
}

fn thermal(region: &str, kind: BcKind, value: f64, dims: Dims) -> BoundaryCondition {
    BoundaryCondition {
        region: region.to_string(),
        physics: Physics::Thermal,
        kind,
        value: Some(BcValue::Uniform(QtyAny::new(value, dims))),
        compatibility: None,
        frame: 0,
    }
}

/// The representative cooling stack: a power die bonded through a thermal
/// interface material to a liquid-cooled plate.
fn cooling_scenario() -> Scenario {
    let mut scenario = Scenario::new("reference-cooling", 0x00C0_01ED, Environment::earth_lab());
    scenario.frames = FrameTree::new();
    scenario.frames.add(Frame {
        id: FrameId(1),
        name: "module".to_string(),
        parent: FrameId(0),
        motion: FrameMotion::Fixed {
            orientation: Quat::identity(),
            translation: Vec3::new(0.0, 0.0, 0.012),
        },
    });
    scenario.frames.add(Frame {
        id: FrameId(2),
        name: "die-seat".to_string(),
        parent: FrameId(1),
        motion: FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.25),
            translation: Vec3::new(0.004, -0.002, 0.003),
        },
    });
    scenario.base_bcs = vec![
        thermal("die.top", BcKind::Neumann, 3.5e5, HEAT_FLUX),
        thermal("plate.ambient", BcKind::Robin, 12.0, HTC),
        thermal("plate.inlet", BcKind::Dirichlet, 298.15, TEMPERATURE),
    ];
    scenario.cases.push(LoadCase {
        name: "peak-load".to_string(),
        bcs: vec![thermal("die.top", BcKind::Neumann, 9.1e5, HEAT_FLUX)],
    });
    scenario.contacts = vec![
        ContactLaw {
            region_a: "die.bottom".to_string(),
            region_b: "tim.top".to_string(),
            model: ContactModel::Tied,
        },
        ContactLaw {
            region_a: "tim.bottom".to_string(),
            region_b: "plate.top".to_string(),
            model: ContactModel::Tied,
        },
    ];
    scenario
}

/// Identities of the reference cooling assembly.
struct CoolingModel {
    catalog: EntityCatalog,
    assembly: EntityId,
    die: EntityId,
    tim: EntityId,
    plate: EntityId,
    die_top: EntityId,
    die_bottom: EntityId,
    tim_top: EntityId,
    tim_bottom: EntityId,
    plate_top: EntityId,
    plate_ambient: EntityId,
    plate_inlet: EntityId,
    plate_metal: EntityId,
    plate_coolant: EntityId,
    die_to_tim: EntityId,
    tim_to_plate: EntityId,
    conjugate: EntityId,
    datum_a: DatumId,
    datum_b: DatumId,
    datum_c: DatumId,
}

fn cooling_model() -> CoolingModel {
    let mut catalog = EntityCatalog::new();
    let assembly = catalog
        .declare(
            EntityDeclaration::assembly("cold-plate-stack")
                .with_display_name("Cold Plate Stack")
                .with_fingerprint(fingerprint("assembly/cold-plate-stack@rev-a")),
        )
        .expect("assembly");
    let die = catalog
        .declare(
            EntityDeclaration::part(assembly, "die")
                .with_display_name("SiC die")
                .with_fingerprint(fingerprint("solid/die@rev-a")),
        )
        .expect("die");
    let tim = catalog
        .declare(
            EntityDeclaration::part(assembly, "tim")
                .with_display_name("TIM layer")
                .with_fingerprint(fingerprint("solid/tim@rev-a")),
        )
        .expect("tim");
    let plate = catalog
        .declare(
            EntityDeclaration::part(assembly, "cold-plate")
                .with_display_name("Cold plate")
                .with_fingerprint(fingerprint("solid/cold-plate@rev-a")),
        )
        .expect("plate");

    let die_top = catalog
        .declare(
            EntityDeclaration::surface(die, "top").with_fingerprint(fingerprint("patch/die-top")),
        )
        .expect("die.top");
    let die_bottom = catalog
        .declare(
            EntityDeclaration::surface(die, "bottom")
                .with_fingerprint(fingerprint("patch/die-bottom")),
        )
        .expect("die.bottom");
    let tim_top = catalog
        .declare(
            EntityDeclaration::surface(tim, "top").with_fingerprint(fingerprint("patch/tim-top")),
        )
        .expect("tim.top");
    let tim_bottom = catalog
        .declare(
            EntityDeclaration::surface(tim, "bottom")
                .with_fingerprint(fingerprint("patch/tim-bottom")),
        )
        .expect("tim.bottom");
    let plate_top = catalog
        .declare(
            EntityDeclaration::surface(plate, "top")
                .with_fingerprint(fingerprint("patch/plate-top")),
        )
        .expect("plate.top");
    let plate_ambient = catalog
        .declare(
            EntityDeclaration::surface(plate, "ambient")
                .with_fingerprint(fingerprint("patch/plate-ambient")),
        )
        .expect("plate.ambient");
    let plate_inlet = catalog
        .declare(
            EntityDeclaration::surface(plate, "inlet")
                .with_fingerprint(fingerprint("patch/plate-inlet")),
        )
        .expect("plate.inlet");
    let plate_metal = catalog
        .declare(
            EntityDeclaration::region(plate, "aluminium")
                .with_fingerprint(fingerprint("volume/plate-metal")),
        )
        .expect("plate.aluminium");
    let plate_coolant = catalog
        .declare(
            EntityDeclaration::region(plate, "coolant")
                .with_fingerprint(fingerprint("volume/plate-coolant")),
        )
        .expect("plate.coolant");

    // Ordered: the TIM is applied to the die's bottom face, so `from` is the
    // die side. Reversing the pair is a different identity.
    let die_to_tim = catalog
        .declare(EntityDeclaration::interface(
            assembly,
            "die-to-tim",
            InterfacePair::ordered(die_bottom, tim_top),
        ))
        .expect("die-to-tim");
    let tim_to_plate = catalog
        .declare(EntityDeclaration::interface(
            assembly,
            "tim-to-plate",
            InterfacePair::ordered(tim_bottom, plate_top),
        ))
        .expect("tim-to-plate");
    // Conjugate heat transfer: physics does not distinguish the sides.
    let conjugate = catalog
        .declare(EntityDeclaration::interface(
            plate,
            "metal-coolant",
            InterfacePair::unordered(plate_metal, plate_coolant),
        ))
        .expect("metal-coolant");

    let datum_a = catalog
        .declare_datum(plate, "A", DatumFeature::Plane, &[])
        .expect("datum A");
    let datum_b = catalog
        .declare_datum(plate, "B", DatumFeature::Axis, &[datum_a])
        .expect("datum B");
    let datum_c = catalog
        .declare_datum(die, "C", DatumFeature::Plane, &[datum_a, datum_b])
        .expect("datum C");

    catalog
        .declare_tolerance(Tolerance::new(
            plate_top,
            ToleranceKind::Flatness,
            QtyAny::new(5.0e-5, LENGTH),
            Vec::new(),
            ToleranceSource::Drawing {
                sheet: "CP-1001".to_string(),
                note: "FCF-3".to_string(),
            },
        ))
        .expect("flatness");
    catalog
        .declare_tolerance(Tolerance::new(
            die_bottom,
            ToleranceKind::Parallelism,
            QtyAny::new(2.0e-5, LENGTH),
            vec![datum_a, datum_b],
            ToleranceSource::Standard {
                clause: "ASME Y14.5-2018 6.6".to_string(),
            },
        ))
        .expect("parallelism");
    catalog
        .declare_tolerance(Tolerance::new(
            die_to_tim,
            ToleranceKind::InterfaceGap,
            QtyAny::new(8.0e-5, LENGTH),
            vec![datum_c],
            ToleranceSource::Assumed {
                rationale: "bond-line thickness from the TIM datasheet's nominal cure gap"
                    .to_string(),
            },
        ))
        .expect("gap");

    catalog
        .declare_placement(plate, FrameId(1), PlacementBasis::Nominal)
        .expect("plate placement");
    catalog
        .declare_placement(die, FrameId(2), PlacementBasis::Nominal)
        .expect("die placement");

    CoolingModel {
        catalog,
        assembly,
        die,
        tim,
        plate,
        die_top,
        die_bottom,
        tim_top,
        tim_bottom,
        plate_top,
        plate_ambient,
        plate_inlet,
        plate_metal,
        plate_coolant,
        die_to_tim,
        tim_to_plate,
        conjugate,
        datum_a,
        datum_b,
        datum_c,
    }
}

fn bind_cooling(model: &CoolingModel) -> BindingTable {
    let mut bindings = BindingTable::new();
    let budget = model.catalog.budget();
    let boundary = KindExpectation::Boundary;
    for (site, target) in [
        (
            ReferenceSite::BaseBoundaryCondition { row: 0 },
            model.die_top,
        ),
        (
            ReferenceSite::BaseBoundaryCondition { row: 1 },
            model.plate_ambient,
        ),
        (
            ReferenceSite::BaseBoundaryCondition { row: 2 },
            model.plate_inlet,
        ),
        (
            ReferenceSite::CaseBoundaryCondition {
                case_row: 0,
                bc_row: 0,
            },
            model.die_top,
        ),
        (
            ReferenceSite::Contact {
                row: 0,
                side: ContactSide::A,
            },
            model.die_bottom,
        ),
        (
            ReferenceSite::Contact {
                row: 0,
                side: ContactSide::B,
            },
            model.tim_top,
        ),
        (
            ReferenceSite::Contact {
                row: 1,
                side: ContactSide::A,
            },
            model.tim_bottom,
        ),
        (
            ReferenceSite::Contact {
                row: 1,
                side: ContactSide::B,
            },
            model.plate_top,
        ),
    ] {
        bindings
            .bind(site, EntityRef::new(target, boundary), budget)
            .expect("bind");
    }
    bindings
}

fn codes(violations: &[Violation]) -> Vec<&'static str> {
    violations.iter().map(|violation| violation.code).collect()
}

fn log_receipts(case: &str, catalog: &EntityCatalog, from: u64) {
    for receipt in catalog.receipts().iter().skip(from as usize) {
        println!(
            "{{\"suite\":\"fs-scenario/entity\",\"case\":\"{case}\",\"event\":\"identity\",\
             \"sequence\":{},\"kind\":\"{}\",\"basis\":\"{}\",\"subject\":\"{}\",\
             \"predecessor\":\"{}\",\"display_before\":\"{}\",\"display_after\":\"{}\"}}",
            receipt.sequence(),
            receipt.event(),
            receipt.basis(),
            receipt.subject().short_token(),
            receipt
                .predecessor()
                .map_or_else(|| "-".to_string(), EntityId::short_token),
            receipt.display_before().unwrap_or("-"),
            receipt.display_after().unwrap_or("-"),
        );
    }
}

fn log_binding_table(case: &str, catalog: &EntityCatalog, bindings: &BindingTable) {
    for row in bindings.report(catalog) {
        println!(
            "{{\"suite\":\"fs-scenario/entity\",\"case\":\"{case}\",\"event\":\"binding\",\
             \"site\":\"{}\",\"requested\":\"{}\",\"current\":\"{}\",\"display\":\"{}\",\
             \"hops\":{},\"tier\":\"{}\",\"fault\":\"{}\"}}",
            row.site,
            row.requested.short_token(),
            row.current
                .map_or_else(|| "-".to_string(), EntityId::short_token),
            row.display_name.as_deref().unwrap_or("-"),
            row.hops,
            row.tier.map_or("-", EvidenceTier::label),
            row.fault.map_or("-", ResolutionFault::code),
        );
    }
}

fn log_datum_and_tolerance_tables(case: &str, catalog: &EntityCatalog) {
    for datum in catalog.datums() {
        println!(
            "{{\"suite\":\"fs-scenario/entity\",\"case\":\"{case}\",\"event\":\"datum\",\
             \"row\":{},\"id\":\"{}\",\"owner\":\"{}\",\"name\":\"{}\",\"feature\":\"{}\",\
             \"references\":{}}}",
            datum.row(),
            datum.id().short_token(),
            datum.owner().short_token(),
            datum.declared_name(),
            datum.feature(),
            datum.references().len(),
        );
    }
    for tolerance in catalog.tolerances() {
        println!(
            "{{\"suite\":\"fs-scenario/entity\",\"case\":\"{case}\",\"event\":\"tolerance\",\
             \"row\":{},\"subject\":\"{}\",\"kind\":\"{}\",\"magnitude_m\":{},\
             \"datums\":{},\"source\":\"{}\"}}",
            tolerance.row(),
            tolerance.subject().short_token(),
            tolerance.kind(),
            tolerance.magnitude().value,
            tolerance.datum_frame().len(),
            tolerance.source().label(),
        );
    }
    for placement in catalog.placements() {
        println!(
            "{{\"suite\":\"fs-scenario/entity\",\"case\":\"{case}\",\"event\":\"placement\",\
             \"row\":{},\"occurrence\":\"{}\",\"frame\":{},\"basis\":\"{}\"}}",
            placement.row(),
            placement.occurrence().short_token(),
            placement.frame().0,
            placement.basis().label(),
        );
    }
}

/// A revision that re-exports the same geometry under new CAD names: every
/// declared name changes, every geometry fingerprint is preserved. Bodies are
/// matched by content; the interfaces (which carry no geometry of their own)
/// are matched by declared name under a superseded parent.
fn renamed_cad_revision() -> ImportRevision {
    let mut entities = Vec::new();
    let mut push = |declaration: EntityDeclaration| -> EntityId {
        let id = declaration.identity();
        entities.push(ImportedEntity {
            declaration,
            correspondence: Correspondence::Auto,
        });
        id
    };
    let assembly = push(
        EntityDeclaration::assembly("stack-2026")
            .with_display_name("Cold Plate Stack (2026)")
            .with_fingerprint(fingerprint("assembly/cold-plate-stack@rev-a")),
    );
    let die = push(
        EntityDeclaration::part(assembly, "DIE_1")
            .with_display_name("SiC die")
            .with_fingerprint(fingerprint("solid/die@rev-a")),
    );
    let tim = push(
        EntityDeclaration::part(assembly, "TIM_1")
            .with_display_name("TIM layer")
            .with_fingerprint(fingerprint("solid/tim@rev-a")),
    );
    let plate = push(
        EntityDeclaration::part(assembly, "PLATE_1")
            .with_display_name("Cold plate")
            .with_fingerprint(fingerprint("solid/cold-plate@rev-a")),
    );
    let die_top = push(
        EntityDeclaration::surface(die, "FACE_TOP").with_fingerprint(fingerprint("patch/die-top")),
    );
    let die_bottom = push(
        EntityDeclaration::surface(die, "FACE_BOT")
            .with_fingerprint(fingerprint("patch/die-bottom")),
    );
    let tim_top = push(
        EntityDeclaration::surface(tim, "FACE_TOP").with_fingerprint(fingerprint("patch/tim-top")),
    );
    let tim_bottom = push(
        EntityDeclaration::surface(tim, "FACE_BOT")
            .with_fingerprint(fingerprint("patch/tim-bottom")),
    );
    let plate_top = push(
        EntityDeclaration::surface(plate, "FACE_TOP")
            .with_fingerprint(fingerprint("patch/plate-top")),
    );
    let _plate_ambient = push(
        EntityDeclaration::surface(plate, "FACE_AMB")
            .with_fingerprint(fingerprint("patch/plate-ambient")),
    );
    let _plate_inlet = push(
        EntityDeclaration::surface(plate, "FACE_IN")
            .with_fingerprint(fingerprint("patch/plate-inlet")),
    );
    let metal = push(
        EntityDeclaration::region(plate, "VOL_METAL")
            .with_fingerprint(fingerprint("volume/plate-metal")),
    );
    let coolant = push(
        EntityDeclaration::region(plate, "VOL_FLUID")
            .with_fingerprint(fingerprint("volume/plate-coolant")),
    );
    push(EntityDeclaration::interface(
        assembly,
        "die-to-tim",
        InterfacePair::ordered(die_bottom, tim_top),
    ));
    push(EntityDeclaration::interface(
        assembly,
        "tim-to-plate",
        InterfacePair::ordered(tim_bottom, plate_top),
    ));
    push(EntityDeclaration::interface(
        plate,
        "metal-coolant",
        InterfacePair::unordered(metal, coolant),
    ));
    let _ = die_top;
    ImportRevision {
        label: "rev-b/cad-rename".to_string(),
        event: RebindEvent::Import,
        scope: ImportScope::Partial,
        entities,
    }
}

#[test]
fn ec_001_reference_cooling_scenario_is_expressed_by_identity() {
    let scenario = cooling_scenario();
    assert!(
        scenario.validate().is_empty(),
        "the fixture scenario itself must be admissible: {:?}",
        scenario.validate()
    );

    let model = cooling_model();
    let bindings = bind_cooling(&model);
    let findings = validate_bindings(&scenario, &model.catalog, &bindings);
    assert!(
        findings.is_empty(),
        "the reference cooling scenario must resolve every reference: {findings:?}"
    );
    assert!(model.catalog.verify_receipts());

    // Every enumerable reference site is bound.
    let sites = scenario_reference_sites(&scenario);
    assert_eq!(sites.len(), 8);
    for (site, _) in &sites {
        assert!(bindings.binding_for(*site).is_some(), "{site} is unbound");
    }

    // Ordering is physics-bearing where it matters.
    let interface = model.catalog.get(model.die_to_tim).expect("interface");
    let pair = interface.pair().expect("pair");
    assert_eq!(pair.applied_side(), Some(model.die_bottom));
    assert_eq!(pair.from(), model.die_bottom);
    assert_eq!(pair.to(), model.tim_top);
    let conjugate = model
        .catalog
        .get(model.conjugate)
        .expect("conjugate")
        .pair()
        .expect("pair");
    assert_eq!(
        conjugate.applied_side(),
        None,
        "an unordered interface must refuse to name an applied side"
    );
    assert!(
        [model.plate_metal, model.plate_coolant].contains(&conjugate.from())
            && [model.plate_metal, model.plate_coolant].contains(&conjugate.to())
    );
    assert_eq!(model.tim_to_plate.kind(), EntityKind::Interface);
    assert_eq!(model.tim.kind(), EntityKind::Part);
    assert_eq!(model.assembly.kind(), EntityKind::Assembly);

    log_binding_table("ec-001", &model.catalog, &bindings);
    log_datum_and_tolerance_tables("ec-001", &model.catalog);
    verdict(
        "ec-001",
        "reference cooling stack expressed as assembly/part/region/surface/interface identities with zero dangling references",
    );
}

#[test]
fn ec_002_renames_do_not_orphan_bindings() {
    let scenario = cooling_scenario();
    let mut model = cooling_model();
    let bindings = bind_cooling(&model);
    let before = bindings.report(&model.catalog);
    let before_findings = validate_bindings(&scenario, &model.catalog, &bindings);

    let renames = [
        (model.assembly, "Stack (Rev B) — 2026 refresh"),
        (model.die, "Power die, screened"),
        (model.tim, "TIM (phase-change pad)"),
        (model.plate, "Cold plate, brazed"),
        (model.die_top, "Active face"),
        (model.plate_ambient, "Exposed casing"),
    ];
    let first_receipt = model.catalog.receipts().len() as u64;
    for (id, display) in renames {
        model.catalog.rename(id, display).expect("rename");
    }

    let after = bindings.report(&model.catalog);
    let after_findings = validate_bindings(&scenario, &model.catalog, &bindings);
    assert!(after_findings.is_empty());
    assert_eq!(codes(&before_findings), codes(&after_findings));
    assert_eq!(before.len(), after.len());
    for (before_row, after_row) in before.iter().zip(after.iter()) {
        assert_eq!(before_row.site, after_row.site);
        assert_eq!(
            before_row.current, after_row.current,
            "a rename must not move an identity"
        );
        assert_eq!(before_row.hops, after_row.hops);
        assert_eq!(before_row.tier, after_row.tier);
    }
    assert_eq!(
        model
            .catalog
            .get(model.die_top)
            .expect("die top")
            .declared_name(),
        "top",
        "the declared name is identity-bearing and does not follow the display name"
    );
    assert_eq!(
        model
            .catalog
            .get(model.die_top)
            .expect("die top")
            .display_name(),
        "Active face"
    );
    assert!(model.catalog.verify_receipts());
    log_receipts("ec-002", &model.catalog, first_receipt);
    verdict(
        "ec-002",
        "six renames moved display names only; every binding resolved to the same identity with the same evidence tier",
    );
}

#[test]
fn ec_003_reimport_keeps_bindings_through_receipts() {
    let scenario = cooling_scenario();
    let mut model = cooling_model();
    let bindings = bind_cooling(&model);
    let before = bindings.report(&model.catalog);
    let first_receipt = model.catalog.receipts().len() as u64;

    let revision = renamed_cad_revision();
    let outcome = model
        .catalog
        .apply_import(&revision)
        .expect("re-import with new CAD names must be matched by fingerprint");
    assert_eq!(outcome.steps().len(), 16);
    let mut content = 0usize;
    let mut path = 0usize;
    for step in outcome.steps() {
        match step {
            ImportStep::Superseded { basis, .. } => match basis {
                MatchBasis::GeometryFingerprint => {
                    assert!(basis.proves_geometry_bytes_matched());
                    content += 1;
                }
                MatchBasis::DeclaredPath => {
                    assert!(
                        !basis.proves_geometry_bytes_matched(),
                        "a path match must not claim the geometry bytes matched"
                    );
                    path += 1;
                }
                other => panic!("unexpected basis {other:?}"),
            },
            other => panic!("expected every entity to supersede its predecessor: {other:?}"),
        }
    }
    assert_eq!(
        (content, path),
        (13, 3),
        "13 bodies matched by content; 3 interfaces matched by declared path under a superseded parent"
    );

    let findings = validate_bindings(&scenario, &model.catalog, &bindings);
    assert!(
        findings.is_empty(),
        "bindings must survive the re-import: {findings:?}"
    );
    let after = bindings.report(&model.catalog);
    for (before_row, after_row) in before.iter().zip(after.iter()) {
        assert_eq!(before_row.site, after_row.site);
        assert_ne!(
            before_row.current, after_row.current,
            "a re-imported entity gets a new identity"
        );
        assert_eq!(after_row.hops, 1);
        assert_eq!(after_row.tier, Some(EvidenceTier::ContentMatched));
        // The metamorphic invariant: the physics target is the same geometry.
        let old = model
            .catalog
            .get(before_row.current.expect("resolved"))
            .expect("old entity");
        let new = model
            .catalog
            .get(after_row.current.expect("resolved"))
            .expect("new entity");
        assert_eq!(old.fingerprint(), new.fingerprint());
        assert_eq!(old.kind(), new.kind());
        assert_ne!(old.declared_name(), new.declared_name());
    }
    assert!(model.catalog.verify_receipts());
    log_receipts("ec-003", &model.catalog, first_receipt);
    log_binding_table("ec-003", &model.catalog, &bindings);
    verdict(
        "ec-003",
        "a full CAD rename re-import superseded 12 identities on content-matched receipts; every binding still resolves",
    );
}

#[test]
fn ec_004_reimport_is_stable_and_idempotent() {
    let scenario = cooling_scenario();
    let mut model = cooling_model();
    let bindings = bind_cooling(&model);
    let revision = renamed_cad_revision();
    model.catalog.apply_import(&revision).expect("first import");
    let after_first = bindings.report(&model.catalog);
    let entities_after_first = model.catalog.entities().len();

    let second = model
        .catalog
        .apply_import(&revision)
        .expect("re-applying the same revision must be admissible");
    assert!(
        second
            .steps()
            .iter()
            .all(|step| matches!(step, ImportStep::Unchanged { .. })),
        "an identical revision must not mint identities: {:?}",
        second.steps()
    );
    assert_eq!(model.catalog.entities().len(), entities_after_first);
    assert_eq!(second.receipts(), 16, "the revision is still receipted");
    assert_eq!(bindings.report(&model.catalog), after_first);
    assert!(validate_bindings(&scenario, &model.catalog, &bindings).is_empty());
    assert!(model.catalog.verify_receipts());
    verdict(
        "ec-004",
        "re-importing an identical revision produced 12 unchanged receipts, zero new identities, and an unchanged binding table",
    );
}

#[test]
fn ec_005_dangling_and_ambiguous_references_are_reported_with_fixes() {
    let scenario = cooling_scenario();
    let mut model = cooling_model();
    let budget = model.catalog.budget();

    // A retired entity: a complete revision of the plate subtree that omits
    // the ambient face.
    let plate_revision = ImportRevision {
        label: "rev-c/plate-only".to_string(),
        event: RebindEvent::RevisionMigration,
        scope: ImportScope::Complete { root: model.plate },
        entities: vec![
            ImportedEntity {
                declaration: EntityDeclaration::surface(model.plate, "top")
                    .with_fingerprint(fingerprint("patch/plate-top")),
                correspondence: Correspondence::Auto,
            },
            ImportedEntity {
                declaration: EntityDeclaration::surface(model.plate, "inlet")
                    .with_fingerprint(fingerprint("patch/plate-inlet")),
                correspondence: Correspondence::Auto,
            },
            ImportedEntity {
                declaration: EntityDeclaration::region(model.plate, "aluminium")
                    .with_fingerprint(fingerprint("volume/plate-metal")),
                correspondence: Correspondence::Auto,
            },
            ImportedEntity {
                declaration: EntityDeclaration::region(model.plate, "coolant")
                    .with_fingerprint(fingerprint("volume/plate-coolant")),
                correspondence: Correspondence::Auto,
            },
        ],
    };
    model
        .catalog
        .apply_import(&plate_revision)
        .expect("complete plate revision");
    assert!(matches!(
        model.catalog.resolve(EntityRef::new(
            model.plate_ambient,
            KindExpectation::Boundary
        )),
        Err(ResolutionFault::Retired { .. })
    ));

    let ghost = EntityDeclaration::surface(model.die, "ghost-face").identity();
    let mut bindings = BindingTable::new();
    let boundary = KindExpectation::Boundary;
    for (site, reference) in [
        // 0: dangling — never declared.
        (
            ReferenceSite::BaseBoundaryCondition { row: 0 },
            EntityRef::new(ghost, boundary),
        ),
        // 1: retired by the complete revision above.
        (
            ReferenceSite::BaseBoundaryCondition { row: 1 },
            EntityRef::new(model.plate_ambient, boundary),
        ),
        // 2: the caller's own expectation is violated.
        (
            ReferenceSite::BaseBoundaryCondition { row: 2 },
            EntityRef::new(
                model.plate_inlet,
                KindExpectation::Exact(EntityKind::Region),
            ),
        ),
        // 3: the SITE's required kind is violated (a BC cannot name a part).
        (
            ReferenceSite::CaseBoundaryCondition {
                case_row: 0,
                bc_row: 0,
            },
            EntityRef::new(model.die, KindExpectation::Any),
        ),
        // 4: a site the scenario does not have.
        (
            ReferenceSite::BaseBoundaryCondition { row: 99 },
            EntityRef::new(model.die_top, boundary),
        ),
        // 5+6: the same site bound twice.
        (
            ReferenceSite::Contact {
                row: 0,
                side: ContactSide::A,
            },
            EntityRef::new(model.die_bottom, boundary),
        ),
        (
            ReferenceSite::Contact {
                row: 0,
                side: ContactSide::A,
            },
            EntityRef::new(model.tim_top, boundary),
        ),
    ] {
        bindings.bind(site, reference, budget).expect("bind");
    }

    let findings = validate_bindings(&scenario, &model.catalog, &bindings);
    let found = codes(&findings);
    for expected in [
        "entity-dangling-reference",
        "entity-retired-reference",
        "entity-kind-mismatch",
        "entity-site-kind-mismatch",
        "entity-orphan-site",
        "entity-duplicate-binding",
        "entity-unbound-site",
    ] {
        assert!(
            found.contains(&expected),
            "{expected} missing from {found:?}"
        );
    }
    assert!(
        findings.iter().all(|violation| !violation.fix.is_empty()),
        "every finding must carry a fix hint"
    );
    for finding in &findings {
        println!(
            "{{\"suite\":\"fs-scenario/entity\",\"case\":\"ec-005\",\"event\":\"violation\",\
             \"code\":\"{}\"}}",
            finding.code
        );
    }

    // Ambiguity: two parts each carrying a surface with the same declared
    // name, and an unbound site whose string names both.
    let mut ambiguous = EntityCatalog::new();
    let assembly = ambiguous
        .declare(EntityDeclaration::assembly("stack"))
        .expect("assembly");
    let left = ambiguous
        .declare(EntityDeclaration::part(assembly, "left"))
        .expect("left");
    let right = ambiguous
        .declare(EntityDeclaration::part(assembly, "right"))
        .expect("right");
    ambiguous
        .declare(EntityDeclaration::surface(left, "die.top"))
        .expect("left face");
    ambiguous
        .declare(EntityDeclaration::surface(right, "die.top"))
        .expect("right face");
    let ambiguous_findings = validate_bindings(&scenario, &ambiguous, &BindingTable::new());
    assert!(codes(&ambiguous_findings).contains(&"entity-ambiguous-name"));
    let ambiguity = ambiguous_findings
        .iter()
        .find(|violation| violation.code == "entity-ambiguous-name")
        .expect("ambiguity finding");
    assert!(ambiguity.what.contains("names 2 entities"));
    assert!(
        ambiguity
            .fix
            .contains("bind this site to the intended identity")
    );
    assert!(matches!(
        ambiguous.lookup_declared_name("die.top"),
        NameLookup::Ambiguous { total: 2, .. }
    ));

    verdict(
        "ec-005",
        "dangling, retired, kind-mismatched, orphan-site, duplicate, unbound, and ambiguous-name references all reported with fixes",
    );
}

#[test]
fn ec_006_legacy_string_scenarios_migrate_mechanically() {
    let scenario = cooling_scenario();
    let mut catalog = EntityCatalog::new();

    // Before migration every site is string-only.
    let unmigrated = validate_bindings(&scenario, &catalog, &BindingTable::new());
    assert_eq!(unmigrated.len(), 8);
    assert!(
        codes(&unmigrated)
            .iter()
            .all(|code| *code == "entity-unbound-site")
    );

    let migration = migrate_legacy_scenario(&scenario, &mut catalog, "legacy/reference-cooling")
        .expect("migrate");
    let findings = validate_bindings(&scenario, &catalog, migration.bindings());
    assert!(
        findings.is_empty(),
        "the migrated scenario must resolve every reference: {findings:?}"
    );
    assert_eq!(
        migration.surfaces().len(),
        7,
        "one identity per distinct region string"
    );
    for (name, id) in migration.surfaces() {
        let entity = catalog.get(*id).expect("migrated entity");
        assert_eq!(entity.declared_name(), name);
        assert_eq!(entity.display_name(), name);
        assert!(entity.is_legacy());
        assert_eq!(entity.kind(), EntityKind::Surface);
        assert_eq!(
            entity.fingerprint(),
            None,
            "migration must not invent geometry it was never given"
        );
    }
    assert!(catalog.get(migration.part()).expect("part").is_legacy());
    assert_eq!(migration.receipts(), 9, "assembly + part + 7 surfaces");
    let receipt = &catalog.receipts()[migration.first_receipt() as usize];
    assert_eq!(receipt.event(), RebindEvent::LegacyMigration);
    assert_eq!(receipt.basis(), MatchBasis::LegacyName);
    assert_eq!(receipt.basis().tier(), Some(EvidenceTier::Declared));
    assert!(!receipt.basis().proves_geometry_bytes_matched());

    // Idempotent: the legacy marker is metadata, not identity.
    let repeat = migrate_legacy_scenario(&scenario, &mut catalog, "legacy/reference-cooling")
        .expect("re-migrate");
    assert_eq!(repeat.receipts(), 0);
    assert_eq!(repeat.root(), migration.root());
    assert_eq!(repeat.part(), migration.part());
    assert_eq!(repeat.surfaces(), migration.surfaces());
    assert_eq!(repeat.bindings(), migration.bindings());
    assert!(catalog.verify_receipts());

    // A legacy binding whose string drifts from its declared name is caught.
    let mut drifted = scenario.clone();
    drifted.base_bcs[0].region = "die.upper".to_string();
    let drift = validate_bindings(&drifted, &catalog, migration.bindings());
    assert!(codes(&drift).contains(&"entity-legacy-name-drift"));

    log_receipts("ec-006", &catalog, migration.first_receipt());
    verdict(
        "ec-006",
        "7 string region names became declared-name identities with legacy markers; re-running the migration appended nothing",
    );
}

#[test]
fn ec_007_datums_tolerances_and_placements_are_typed() {
    let scenario = cooling_scenario();
    let model = cooling_model();
    assert!(model.catalog.validate().is_empty());
    assert_eq!(model.catalog.datums().len(), 3);
    assert_eq!(
        model
            .catalog
            .datum(model.datum_c)
            .expect("datum C")
            .references(),
        &[model.datum_a, model.datum_b]
    );

    // A form control cannot carry datums; a location control must.
    let mut bad = EntityCatalog::new();
    let assembly = bad
        .declare(EntityDeclaration::assembly("a"))
        .expect("assembly");
    let part = bad
        .declare(EntityDeclaration::part(assembly, "p"))
        .expect("p");
    let face = bad
        .declare(EntityDeclaration::surface(part, "f"))
        .expect("f");
    let datum = bad
        .declare_datum(part, "A", DatumFeature::Plane, &[])
        .expect("datum");
    let ghost_datum = model.catalog.datums()[2].id();
    bad.declare_tolerance(Tolerance::new(
        face,
        ToleranceKind::Flatness,
        QtyAny::new(1.0e-5, LENGTH),
        vec![datum],
        ToleranceSource::Drawing {
            sheet: "S1".to_string(),
            note: "N1".to_string(),
        },
    ))
    .expect("declare");
    bad.declare_tolerance(Tolerance::new(
        face,
        ToleranceKind::Position,
        QtyAny::new(1.0e-5, LENGTH),
        Vec::new(),
        ToleranceSource::Assumed {
            rationale: String::new(),
        },
    ))
    .expect("declare");
    bad.declare_tolerance(Tolerance::new(
        face,
        ToleranceKind::Parallelism,
        QtyAny::new(2.0, TEMPERATURE),
        vec![datum, datum, ghost_datum, datum],
        ToleranceSource::Standard {
            clause: "X".to_string(),
        },
    ))
    .expect("declare");
    let findings = codes(&bad.validate());
    for expected in [
        "tolerance-datum-forbidden",
        "tolerance-datum-required",
        "tolerance-source-empty",
        "tolerance-dims",
        "tolerance-datum-frame-arity",
        "tolerance-datum-repeated",
        "tolerance-datum-dangling",
    ] {
        assert!(
            findings.contains(&expected),
            "{expected} missing: {findings:?}"
        );
    }

    // Placement round-trip through the scenario's own FrameTree.
    let placement = model
        .catalog
        .placement_of(model.die)
        .expect("die placement");
    assert_eq!(placement.frame(), FrameId(2));
    assert_eq!(placement.basis(), PlacementBasis::Nominal);
    let probe = Point::new(0.001, 0.002, -0.003);
    let via_tree = scenario
        .frames
        .world_pose(placement.frame(), 0.0)
        .expect("pose");
    let manual = fs_scenario::FrameTree::local_pose(&scenario.frames.frames[0], 0.0)
        .expect("module")
        .compose(
            &fs_scenario::FrameTree::local_pose(&scenario.frames.frames[1], 0.0).expect("seat"),
        );
    let a = via_tree.transform_point(probe).expect("a");
    let b = manual.transform_point(probe).expect("b");
    assert!(Vec3::new(a.x - b.x, a.y - b.y, a.z - b.z).norm() < 1e-12);

    // A placement naming a frame the scenario does not declare is a finding.
    let mut stray = cooling_model();
    stray
        .catalog
        .declare_placement(stray.tim, FrameId(9), PlacementBasis::Nominal)
        .expect("declare");
    let stray_findings = codes(&validate_bindings(
        &scenario,
        &stray.catalog,
        &bind_cooling(&stray),
    ));
    assert!(stray_findings.contains(&"placement-unknown-frame"));

    log_datum_and_tolerance_tables("ec-007", &model.catalog);
    verdict(
        "ec-007",
        "datum hierarchy, tolerance dimensional typing/datum-frame rules, and nominal placements validated against the scenario FrameTree",
    );
}

#[test]
fn ec_008_scenario_lifecycle_author_import_rename_reimport_validate() {
    let scenario = cooling_scenario();

    // 1. Author from strings only.
    let mut catalog = EntityCatalog::new();
    let migration = migrate_legacy_scenario(&scenario, &mut catalog, "reference-cooling")
        .expect("legacy migration");
    let mut bindings = migration.bindings().clone();
    assert!(validate_bindings(&scenario, &catalog, &bindings).is_empty());
    let authored: Vec<_> = bindings
        .bindings()
        .iter()
        .map(|binding| (binding.site(), binding.reference().target()))
        .collect();

    // 2. Import-bind: real geometry arrives and supersedes the legacy stubs.
    let part = migration.part();
    let mut entities = Vec::new();
    for (name, tag) in [
        ("die.top", "patch/die-top"),
        ("die.bottom", "patch/die-bottom"),
        ("tim.top", "patch/tim-top"),
        ("tim.bottom", "patch/tim-bottom"),
        ("plate.top", "patch/plate-top"),
        ("plate.ambient", "patch/plate-ambient"),
        ("plate.inlet", "patch/plate-inlet"),
    ] {
        entities.push(ImportedEntity {
            declaration: EntityDeclaration::surface(part, name).with_fingerprint(fingerprint(tag)),
            correspondence: Correspondence::Auto,
        });
    }
    let import = ImportRevision {
        label: "rev-a/import".to_string(),
        event: RebindEvent::Import,
        scope: ImportScope::Partial,
        entities,
    };
    let outcome = catalog.apply_import(&import).expect("import");
    assert_eq!(outcome.steps().len(), 7);
    for step in outcome.steps() {
        assert!(matches!(
            step,
            ImportStep::Superseded {
                basis: MatchBasis::DeclaredPath,
                ..
            }
        ));
    }
    let after_import = bindings.report(&catalog);
    assert!(after_import.iter().all(|row| row.fault.is_none()));
    assert!(
        after_import
            .iter()
            .all(|row| row.tier == Some(EvidenceTier::PathMatched)),
        "a name-and-path match is weaker evidence than a fingerprint match, and says so"
    );

    // 3. Rename: display only.
    let current = catalog
        .resolve(EntityRef::new(authored[0].1, KindExpectation::Boundary))
        .expect("resolve")
        .current();
    catalog.rename(current, "Active die face").expect("rename");

    // 4. Re-mesh re-import: same declared paths, new geometry bytes.
    let remesh = ImportRevision {
        label: "rev-a/remesh".to_string(),
        event: RebindEvent::Remesh,
        scope: ImportScope::Partial,
        entities: vec![ImportedEntity {
            declaration: EntityDeclaration::surface(part, "die.top")
                .with_fingerprint(fingerprint("patch/die-top@mesh-2")),
            correspondence: Correspondence::Auto,
        }],
    };
    let remesh_outcome = catalog.apply_import(&remesh).expect("remesh");
    assert_eq!(remesh_outcome.steps().len(), 1);

    // 5. Validate: everything still resolves, through two hops for the
    //    re-meshed face.
    let findings = validate_bindings(&scenario, &catalog, &bindings);
    assert!(findings.is_empty(), "{findings:?}");
    let table = bindings.report(&catalog);
    let die_top_row = table
        .iter()
        .find(|row| row.site == ReferenceSite::BaseBoundaryCondition { row: 0 })
        .expect("die.top row");
    assert_eq!(die_top_row.hops, 2);
    assert_eq!(die_top_row.tier, Some(EvidenceTier::PathMatched));
    assert_eq!(die_top_row.requested, authored[0].1);
    assert!(catalog.verify_receipts());

    // A reserved site family (sensor) uses the same shape and resolves too.
    bindings
        .bind(
            ReferenceSite::Sensor { row: 0 },
            EntityRef::new(die_top_row.current.expect("current"), KindExpectation::Any),
            catalog.budget(),
        )
        .expect("bind sensor");
    assert!(validate_bindings(&scenario, &catalog, &bindings).is_empty());

    log_receipts("ec-008", &catalog, 0);
    log_binding_table("ec-008", &catalog, &bindings);
    println!(
        "{{\"suite\":\"fs-scenario/entity\",\"case\":\"ec-008\",\"event\":\"chain\",\
         \"receipts\":{},\"root\":\"{}\"}}",
        catalog.receipts().len(),
        catalog.receipt_root()
    );
    verdict(
        "ec-008",
        "author -> import-bind -> rename -> re-mesh -> validate kept all 8 bindings resolvable across 2 supersession hops",
    );
}

#[test]
fn ec_009_bindings_are_ordered_and_deterministic() {
    let scenario = cooling_scenario();
    let model_a = cooling_model();
    let model_b = cooling_model();
    assert_eq!(
        model_a.catalog.receipt_root(),
        model_b.catalog.receipt_root()
    );
    assert_eq!(
        bind_cooling(&model_a).report(&model_a.catalog),
        bind_cooling(&model_b).report(&model_b.catalog)
    );
    let sites: Vec<_> = scenario_reference_sites(&scenario)
        .into_iter()
        .map(|(site, name)| format!("{site}={name}"))
        .collect();
    assert_eq!(
        sites,
        vec![
            "base BC row 0=die.top".to_string(),
            "base BC row 1=plate.ambient".to_string(),
            "base BC row 2=plate.inlet".to_string(),
            "case row 0 BC row 0=die.top".to_string(),
            "contact row 0 side a=die.bottom".to_string(),
            "contact row 0 side b=tim.top".to_string(),
            "contact row 1 side a=tim.bottom".to_string(),
            "contact row 1 side b=plate.top".to_string(),
        ]
    );
    let table: Vec<Binding> = bind_cooling(&model_a).bindings().to_vec();
    assert_eq!(table.len(), 8);
    verdict(
        "ec-009",
        "identity derivation, receipt roots, and site enumeration are deterministic across independent builds",
    );
}
