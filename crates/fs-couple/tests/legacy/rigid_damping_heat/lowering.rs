//! The pairing declared as a draft and lowered by the owner constructors
//! (`fs_system::lowering::owners`), against the declared graph above.
//!
//! Four 1 m boxes of density 2 (2 kg each) at x = 0, 1, 2, 3 m, moving along
//! x at 1, -1, 0.5 and -0.5 m/s, joined in two pairs by bushings of 0.7 and
//! 1.1 N s/m; two cells of 2 and 3 kg on the fixture laws, each on a bed
//! with the fixture's reservoir (to the face of a stationary `lab` frame);
//! loss routes of the fixture's fractions, each from the face of a bushing's
//! first body; and a constant source on cell 0. Every interface joins faces. The oracle is `graph()` on the lowered mass
//! properties, read back from the lowering's receipts. Only the source is
//! declared: the fixture's contact and its 4 W at 0.35 s are thermal events
//! that a draft does not declare. Fixture values; not a seed-dataset
//! authority.

use super::*;
use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterfaceSystemCard, InterpolationPolicy,
    LawId, LawParameter, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::draft::{
    DraftAssembly, DraftAssemblyInput, DraftEnvironment, DraftGeometrySource, DraftInitialState,
    DraftInterfaceDecl, DraftInterfaceKind, DraftPartDecl, DraftPartKind, FidelityErrorBudget,
    FidelityPolicy, FidelityRung, PhysicsDomain, PrimitiveShape, SourceDecl, SourceTarget,
    StepPolicy, TimeHistory, TimePlan,
};
use fs_scenario::assembly::{MaterialRef, ParameterKey, ReferenceState};
use fs_scenario::{EntityCatalog, EntityDeclaration, EntityId, FrameId, InterfacePair, Scenario};
use fs_system::lowering::HEAT_ROUTE_FAMILY;
use fs_system::lowering::owners::thermal::{
    self, CaloricCell, CaloricLaws, CaloricRequest, HEAT_SOURCE_PORT, SENSIBLE_HEAT_LAW,
};
use fs_system::lowering::owners::{
    LOSS_ROUTES, THERMAL_BOUNDARY_FAMILY, VISCOUS_BUSHING_FAMILY, mbd, register_owner_defaults,
    register_owners,
};
use fs_system::{
    CardStore, ConstitutiveModelRef, LoweredSystem, LoweringError, LoweringRegistry,
    ParameterSource, ProvenanceRecorder,
};
use std::collections::BTreeSet;
use std::rc::Rc;

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const SPECIFIC_HEAT: Dims = Dims([2, 0, -2, -1, 0, 0]);
const DAMPING: Dims = Dims([0, 1, -1, 0, 0, 0]);
const TORSIONAL_DAMPING: Dims = Dims([2, 1, -1, 0, 0, 0]);
const HEAT_TRANSFER: Dims = Dims([0, 1, -3, -1, 0, 0]);
const AREA: Dims = Dims([2, 0, 0, 0, 0, 0]);
const CONDUCTANCE: Dims = Dims([2, 1, -3, -1, 0, 0]);
const POWER: Dims = Dims([2, 1, -3, 0, 0, 0]);

fn provenance(what: &str) -> Provenance {
    let source = format!("fixture value; not a seed-dataset authority ({what})");
    let artifact = Some(hash_bytes(source.as_bytes()));
    Provenance {
        source,
        license: "CC-BY-4.0".to_owned(),
        artifact,
    }
}

fn claims(values: &[(&str, f64, Dims)]) -> ClaimSet {
    let mut set = ClaimSet::new();
    for &(name, value, dims) in values {
        let claim = PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::HalfWidth {
                half_width: 1e-6,
                confidence: 0.95,
            },
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: provenance(name),
        };
        set.insert_claim(claim).unwrap();
    }
    set
}

fn state(name: &str) -> MaterialStateId {
    let fixture = || "fixture".to_owned();
    MaterialStateId {
        chemistry: name.to_owned(),
        phase: fixture(),
        process: fixture(),
        revision: 0,
    }
}

fn model(law: &str, parameters: &[(&str, f64, Dims)]) -> ConstitutiveModelCard {
    let parameters = parameters.iter();
    let parameters =
        parameters.map(|&(name, value, dims)| (name.to_owned(), LawParameter { value, dims }));
    ConstitutiveModelCard {
        law: LawId(law.to_owned()),
        law_version: 1,
        parameters: parameters.collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained(),
        sources: Vec::new(),
        provenance: provenance(law),
    }
}

fn material(values: &[(&str, f64, Dims)], models: Vec<ConstitutiveModelCard>) -> MaterialCard {
    MaterialCard::assemble(state("fixture solid"), claims(values), models).unwrap()
}

fn interface_card(family: &str, values: &[(&str, f64, Dims)]) -> InterfaceSystemCard {
    let surface = |side: &str| SurfaceSpec {
        material: state(side),
        texture_frame: side.to_owned(),
    };
    let fixture = || "fixture".to_owned();
    let context = SystemContext {
        medium: family.to_owned(),
        third_body: None,
        environment: fixture(),
        history: fixture(),
    };
    InterfaceSystemCard::assemble(
        surface("a"),
        surface("b"),
        context,
        claims(values),
        Vec::new(),
    )
    .unwrap()
}

/// The fixture's two laws on cards: law ids by storage basis, `a` and `b`
/// as model parameters.
fn law_id(basis: ThermalStorageBasis) -> String {
    format!("fixture-quadratic-{basis:?}")
}

fn fixture_cards(laws: &[Caloric; 2]) -> [MaterialCard; 2] {
    core::array::from_fn(|c| {
        let law = &laws[c];
        let ab = [
            ("a", law.a, SPECIFIC_HEAT),
            ("b", law.b, Dims([2, 0, -2, -2, 0, 0])),
        ];
        material(
            &[("density", 2.0 + c as f64, DENSITY)],
            vec![model(&law_id(law.basis), &ab)],
        )
    })
}

/// The fixture laws, registered from the test: `a` and `b` from the model,
/// and the energy of the initial temperature from the law's own 300 K datum.
fn fixture_laws() -> CaloricLaws {
    let mut laws = CaloricLaws::new();
    for basis in [
        ThermalStorageBasis::Enthalpy,
        ThermalStorageBasis::InternalEnergy,
    ] {
        laws.register(
            law_id(basis),
            move |request: &CaloricRequest<'_>, r: &mut ProvenanceRecorder| {
                let model = request.model;
                let mut parameter = |name: &str| {
                    let LawParameter { value, dims } = model.parameters[name];
                    let reference = ConstitutiveModelRef {
                        card_hash: request.material.card,
                        model_hash: model.content_hash(),
                        law: &model.law.0,
                        law_version: model.law_version,
                        parameter: name,
                    };
                    r.record_model_parameter(name, value, dims, reference);
                    value
                };
                let law = Caloric {
                    a: parameter("a"),
                    b: parameter("b"),
                    maximum: 1e6,
                    basis,
                };
                let t = request.temperature_k - 300.0;
                let specific_energy_j_kg = law.a * t + 0.5 * law.b * t * t;
                Ok(CaloricCell {
                    law: Rc::new(law),
                    specific_energy_j_kg,
                })
            },
        );
    }
    laws
}

/// The pairing as a draft, its cells on `cell_cards`.
#[allow(clippy::too_many_lines)] // One draft literal: entities, parts, interfaces, cards.
fn pairing(cell_cards: [MaterialCard; 2], source_w: f64) -> (DraftAssembly, CardStore) {
    let mut cards = CardStore::new();
    let body_card = cards.insert_material(material(&[("density", 2.0, DENSITY)], Vec::new()));
    let cell_cards = cell_cards.map(|card| cards.insert_material(card));
    let mut catalog = EntityCatalog::new();
    let root = catalog
        .declare(EntityDeclaration::assembly("pairing"))
        .unwrap();
    // Each part and its face: only surfaces and regions join interfaces.
    let mut part = |name: String| {
        let part = catalog.declare(EntityDeclaration::part(root, &name));
        let part = part.unwrap();
        let face = catalog.declare(EntityDeclaration::surface(part, "face"));
        (part, face.unwrap())
    };
    let bodies: [(EntityId, EntityId); 4] = core::array::from_fn(|b| part(format!("body{b}")));
    let cells: [(EntityId, EntityId); 2] = core::array::from_fn(|c| part(format!("cell{c}")));
    let lab = part("lab".to_owned());
    let mut iface = |name: &str, [a, b]: [EntityId; 2]| {
        let pair = InterfacePair::ordered(a, b);
        let entity = catalog
            .declare(EntityDeclaration::interface(root, name, pair))
            .unwrap();
        (name.to_owned(), entity, a, b)
    };
    let bushings = [
        iface("bushing0", [bodies[0].1, bodies[1].1]),
        iface("bushing1", [bodies[2].1, bodies[3].1]),
    ];
    let beds = [
        iface("bed0", [cells[0].1, lab.1]),
        iface("bed1", [cells[1].1, lab.1]),
    ];
    let route_ends = [(0, 0), (0, 1), (1, 0), (1, 1)];
    let route = |(b, c): (usize, usize)| [bodies[2 * b].1, cells[c].1];
    let losses = route_ends.map(|ends| iface(&format!("route{}{}", ends.0, ends.1), route(ends)));

    let m = |value| QtyAny::new(value, LENGTH);
    let geometry = DraftGeometrySource::Primitive(PrimitiveShape::Box {
        lx: m(1.0),
        ly: m(1.0),
        lz: m(1.0),
    });
    let fidelity = FidelityPolicy {
        rungs: vec![FidelityRung::Lumped],
        error_budget: FidelityErrorBudget {
            relative_tol: 1e-5,
            absolute_tol: QtyAny::dimensionless(1e-10),
            max_work_units: Some(10_000),
        },
    };
    let initial = |x: f64, v: f64| DraftInitialState {
        position: [m(x), m(0.0), m(0.0)],
        orientation: [1.0, 0.0, 0.0, 0.0],
        linear_velocity: [v, 0.0, 0.0].map(|v| QtyAny::new(v, SPEED)),
        angular_velocity: [QtyAny::new(0.0, Dims([0, 0, -1, 0, 0, 0])); 3],
        temperature: None,
    };
    let declare =
        |name: String, entity, kind, card, physics, start: DraftInitialState| DraftPartDecl {
            name,
            entity,
            kind,
            geometry: geometry.clone(),
            materials: vec![MaterialRef {
                card,
                query_point: vec![(
                    ParameterKey::new("T", TEMPERATURE),
                    QtyAny::new(300.0, TEMPERATURE),
                )],
                frame: FrameId(0),
                layer: None,
                unknown: Vec::new(),
            }],
            physics: BTreeSet::from([physics]),
            fidelity: fidelity.clone(),
            uncertainty: Vec::new(),
            initial: Some(start),
        };
    let mut parts = Vec::new();
    for (b, v) in [1.0, -1.0, 0.5, -0.5].into_iter().enumerate() {
        let start = initial(b as f64, v);
        parts.push(declare(
            format!("body{b}"),
            bodies[b].0,
            DraftPartKind::RigidBody,
            body_card,
            PhysicsDomain::Mechanics,
            start,
        ));
    }
    for c in 0..2 {
        let start = initial(c as f64 + 0.5, 0.0);
        parts.push(declare(
            format!("cell{c}"),
            cells[c].0,
            DraftPartKind::Solid,
            cell_cards[c],
            PhysicsDomain::Thermal,
            start,
        ));
    }
    // The stationary frame the beds join, after the cells (receipts 4 and 5).
    let (kind, physics) = (DraftPartKind::DrivenBody, PhysicsDomain::Mechanics);
    let frame = declare(
        "lab".to_owned(),
        lab.0,
        kind,
        body_card,
        physics,
        initial(0.0, 0.0),
    );
    parts.push(DraftPartDecl {
        initial: None,
        ..frame
    });

    let contact = |(name, entity, feature_a, feature_b), kind| DraftInterfaceDecl {
        name,
        entity,
        feature_a,
        feature_b,
        kind,
    };
    let mut interfaces = Vec::new();
    for (bushing, c) in bushings.into_iter().zip([0.7, 1.1]) {
        let mut damping = vec![
            ("damping_x", c, DAMPING),
            ("damping_y", 0.0, DAMPING),
            ("damping_z", 0.0, DAMPING),
        ];
        damping.extend(
            ["damping_rx", "damping_ry", "damping_rz"].map(|name| (name, 0.0, TORSIONAL_DAMPING)),
        );
        let law_card = cards
            .insert_interface(
                VISCOUS_BUSHING_FAMILY,
                interface_card(VISCOUS_BUSHING_FAMILY, &damping),
            )
            .unwrap();
        let kind = DraftInterfaceKind::MechanicalContact {
            law_card,
            third_body: None,
            state: Vec::new(),
        };
        interfaces.push(contact(bushing, kind));
    }
    let exchange = [
        ("convection_coefficient", 0.0, HEAT_TRANSFER),
        ("reservoir_area", 0.1, AREA),
        ("reservoir_conductance", 0.0, CONDUCTANCE),
    ];
    let bed_card = cards
        .insert_interface(
            THERMAL_BOUNDARY_FAMILY,
            interface_card(THERMAL_BOUNDARY_FAMILY, &exchange),
        )
        .unwrap();
    for bed in beds {
        let kind = DraftInterfaceKind::ThermalContact {
            law_card: bed_card,
            third_body: None,
            state: Vec::new(),
        };
        interfaces.push(contact(bed, kind));
    }
    let fractions = routes();
    for (route, (b, c)) in losses.into_iter().zip(route_ends) {
        let fraction = QtyAny::dimensionless(fractions[b][&CellId(c as u64)]);
        let state = vec![(ParameterKey::new("fraction", Dims::NONE), fraction)];
        interfaces.push(contact(route, DraftInterfaceKind::HeatRoute { state }));
    }

    let mut environment = DraftEnvironment::earth_lab_air();
    environment.base.gravity = [QtyAny::new(0.0, Dims([1, 0, -2, 0, 0, 0])); 3];
    environment.base.ambient_temperature = QtyAny::new(300.0, TEMPERATURE);
    let heater = (source_w > 0.0).then(|| SourceDecl {
        name: "heater".to_owned(),
        target: SourceTarget::Port {
            entity: cells[0].0,
            port: HEAT_SOURCE_PORT.to_owned(),
        },
        history: TimeHistory::Constant(QtyAny::new(source_w, POWER)),
    });
    let scenario = Scenario::new("pairing", 112, environment.base.clone());
    let input = DraftAssemblyInput {
        name: "pairing".to_owned(),
        root,
        catalog,
        reference_state: Some(ReferenceState::Unloaded),
        environment,
        time_plan: TimePlan {
            horizon: QtyAny::new(0.5, TIME),
            step_policy: StepPolicy::Fixed {
                dt: QtyAny::new(0.125, TIME),
            },
            default_sample_rate: QtyAny::new(8.0, Dims([0, 0, -1, 0, 0, 0])),
        },
        parts,
        joints: Vec::new(),
        gear_meshes: Vec::new(),
        interfaces,
        sources: heater.into_iter().collect(),
        controllers: Vec::new(),
        probes: Vec::new(),
    };
    (DraftAssembly::admit(input, &scenario).unwrap(), cards)
}

fn lower(draft: &DraftAssembly, cards: &CardStore) -> LoweredSystem {
    let mut registry = LoweringRegistry::with_kernel_defaults();
    register_owners(&mut registry, fixture_laws());
    registry.lower(draft, cards).unwrap()
}

/// The bodies' mass properties as the lowering recorded them.
fn lowered_mass(lowered: &LoweredSystem) -> [MassProperties; 4] {
    core::array::from_fn(|b| {
        let sources = &lowered.receipt.components[b].parameter_sources;
        let value = |name: &str| sources.iter().find(|p| p.parameter == name).unwrap().value;
        let inertia = Vec3::new(value("inertia_x"), value("inertia_y"), value("inertia_z"));
        MassProperties::new(value("mass"), Vec3::ZERO, inertia).unwrap()
    })
}

fn bodies_of(system: &System, node: NodeId) -> &[RigidBodyState] {
    system.state::<Vec<RigidBodyState>>(node).unwrap()
}

fn temperature_bits(system: &System, node: NodeId) -> [u64; 2] {
    let timeline = system.state::<TimelineState<'static>>(node).unwrap();
    [0, 1].map(|c| {
        timeline
            .simulation()
            .temperature_k(CellId(c))
            .unwrap()
            .to_bits()
    })
}

#[test]
fn the_lowered_draft_steps_with_the_declared_graph_bit_for_bit() {
    let laws = static_laws();
    for source_w in [0.0, 2.0] {
        let (draft, cards) = pairing(fixture_cards(laws), source_w);
        let mut bodies = Bodies::new([0.7, 1.1]);
        bodies.properties = lowered_mass(&lower(&draft, &cards));
        // The oracle takes these from the receipts, so hold them to the box:
        // 2 kg, and 2 kg (1/12 + 1/12) m^2 about each axis.
        let i = 2.0 * (1.0 / 12.0 + 1.0 / 12.0);
        let declared = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(i, i, i)).unwrap();
        assert_eq!(bodies.properties, [declared; 4]);
        let schedule = vec![birth(0, 0.0, source_w), birth(1, 0.0, 0.0)];
        let heat = DepositionSimulation::new(
            schedule,
            &[(MaterialId(0), &laws[0]), (MaterialId(1), &laws[1])],
            budget(),
        );
        let heat = heat
            .unwrap()
            .advance_to(0.0, config().thermal, || Ok(()))
            .unwrap()
            .0;
        for target_s in [0.5, 0.45] {
            let case = format!("{source_w} W to {target_s} s");
            let (mut oracle, network, cells) = graph(&bodies, &routes(), &heat);
            let mut lowered = lower(&draft, &cards);
            let nodes = [
                lowered.node_map["viscous-network"],
                lowered.node_map["thermal-timeline"],
            ];
            let system = &mut lowered.system;
            let gate = CancelGate::new_clock_free();
            with_cx(&gate, |cx| {
                while oracle.time_s() < target_s {
                    let theirs = oracle.step(target_s - oracle.time_s(), cx).unwrap().dt_s;
                    let step = system.step(target_s - system.time_s(), cx).unwrap().dt_s;
                    assert_eq!(step.to_bits(), theirs.to_bits(), "{case}: the step");
                    assert_eq!(
                        bodies_of(system, nodes[0]),
                        bodies_of(&oracle, network),
                        "{case}"
                    );
                    let [a, b] = [
                        temperature_bits(system, nodes[1]),
                        temperature_bits(&oracle, cells),
                    ];
                    assert_eq!(a, b, "{case}: the cells");
                }
            });
            assert_eq!(
                system.time_s().to_bits(),
                oracle.time_s().to_bits(),
                "{case}"
            );
            // The same accounts, and the lowered ledger closes.
            let (got, want) = (system.ledger(), oracle.ledger());
            let dissipated = node(want, "network").dissipated_j.value();
            assert!(
                dissipated > 0.5,
                "{case}: enough heat to compare, {dissipated} J"
            );
            // The step histories are equal to the bit, so each account is too.
            let same = |a: f64, b: f64, what: &str| {
                assert_eq!(a.to_bits(), b.to_bits(), "{case}: {what}: {a:e} vs {b:e}");
            };
            same(
                node(got, "viscous-network").dissipated_j.value(),
                dissipated,
                "dissipated",
            );
            let stored = |ledger: &Ledger, name: &str| node(ledger, name).stored_change_j.value();
            same(
                stored(got, "thermal-timeline"),
                stored(want, "cells"),
                "stored heat",
            );
            for name in ["unallocated", "environment"] {
                same(sink(got, name), sink(want, name), name);
            }
            // The residual is a sum of accounts that cancel, so it closes to
            // rounding: 1e-12 of the dissipated energy.
            let residual = got.energy_residual_j();
            assert!(
                residual.abs() <= 1e-12 * dissipated,
                "{case}: the ledger closes: {residual:e} J"
            );
        }
    }
}

#[test]
fn cells_on_the_product_law_lower_with_their_card_receipts() {
    let card = |c: f64| {
        let values = [
            ("density", 2.0 + c, DENSITY),
            ("specific_heat_capacity", 500.0, SPECIFIC_HEAT),
        ];
        material(
            &values,
            vec![model(
                SENSIBLE_HEAT_LAW,
                &[("datum_temperature", 0.0, TEMPERATURE)],
            )],
        )
    };
    let (draft, cards) = pairing([card(0.0), card(1.0)], 2.0);
    let mut registry = LoweringRegistry::with_kernel_defaults();
    register_owner_defaults(&mut registry);
    let mut lowered = registry.lower(&draft, &cards).unwrap();
    for cell in &lowered.receipt.components[4..6] {
        assert!(
            cell.property_receipts
                .iter()
                .any(|r| r.property == "specific_heat_capacity")
        );
        let datum = cell
            .parameter_sources
            .iter()
            .find(|p| p.parameter == "datum_temperature");
        assert!(matches!(
            datum.unwrap().source,
            ParameterSource::ConstitutiveModelParam { .. }
        ));
    }
    let (intervals, _) = run_counting(&mut lowered.system, 0.5);
    assert!(
        intervals.iter().all(|&dt| dt <= 0.125),
        "the time plan bounds the step: {intervals:?}"
    );
    let ledger = lowered.system.ledger();
    assert!(node(ledger, "thermal-timeline").stored_change_j.value() > 0.0);
    assert!(
        ledger.energy_residual_j().abs() <= 1e-11,
        "{}",
        ledger.energy_residual_j()
    );
}

#[test]
fn the_product_laws_band_refuses_a_cell_outside_it_and_never_clamps() {
    // The cells start at the 300 K ambient.
    let lower = |lo: f64, hi: f64| {
        let datum = [("datum_temperature", 0.0, TEMPERATURE)];
        let mut law = model(SENSIBLE_HEAT_LAW, &datum);
        law.validity = ValidityDomain::unconstrained().with("T", lo, hi);
        let c = ("specific_heat_capacity", 500.0, SPECIFIC_HEAT);
        let card = || material(&[("density", 2.0, DENSITY), c], vec![law.clone()]);
        let (draft, cards) = pairing([card(), card()], 2.0);
        let mut registry = LoweringRegistry::with_kernel_defaults();
        register_owner_defaults(&mut registry);
        registry.lower(&draft, &cards)
    };
    let refusal = lower(301.0, 310.0).err().unwrap().to_string();
    assert!(refusal.contains("outside the law's band"), "{refusal}");
    let mut lowered = lower(290.0, 300.000_1).unwrap();
    let sources = lowered.receipt.components[4].parameter_sources.iter();
    let band = sources.filter(|p| p.parameter.starts_with("validity.T."));
    let band: Vec<f64> = band.map(|p| p.value).collect();
    assert_eq!(band, [290.0, 300.000_1]);
    // 1 J into 1000 J/K by 0.5 s takes cell 0 past the band's top: the law
    // refuses the step, and no temperature is clamped to the band.
    let gate = CancelGate::new_clock_free();
    let refused = with_cx(&gate, |cx| {
        (0..1000).find_map(|_| lowered.system.step(0.5 - lowered.system.time_s(), cx).err())
    });
    assert!(refused.is_some() && lowered.system.time_s() < 0.5);
}

#[test]
fn each_owner_group_deleted_refuses_by_what_it_lowered() {
    let (draft, cards) = pairing(fixture_cards(static_laws()), 2.0);
    let deleted = [
        (mbd::GROUP, "must be a lowered viscous bushing"),
        (
            thermal::GROUP,
            "feature_b must be a cell of the thermal timeline",
        ),
        (LOSS_ROUTES, "heat-route"),
    ];
    for (group, needle) in deleted {
        let mut registry = LoweringRegistry::with_kernel_defaults();
        register_owners(&mut registry, fixture_laws());
        registry.register_group(group, |_, _| Ok(Vec::new()), |_, _| unreachable!());
        let refusal = registry.lower(&draft, &cards).err().unwrap();
        assert!(refusal.to_string().contains(needle), "{group}: {refusal}");
        if group == LOSS_ROUTES {
            let LoweringError::MissingInterfaceConstructor {
                interface,
                card_family,
                ..
            } = refusal
            else {
                panic!("{refusal:?}");
            };
            assert_eq!(
                (interface.as_str(), card_family.as_str()),
                ("route00", HEAT_ROUTE_FAMILY)
            );
        }
    }
}
