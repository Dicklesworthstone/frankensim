//! Integration tests for bead `frankensim-rc-1109-root-3gc12.11.2`:
//! owner-extensible lowering registry, Slice-A graph equivalence, unseen-device
//! lowering with zero Rust changes, named refusals, receipt completeness, and
//! `AcousticAssembly` family registry entries.

use std::collections::{BTreeMap, BTreeSet};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_bytes};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_ports::PortKind;
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::draft::{
    ControllerDecl, ControllerLaw, DraftAssembly, DraftAssemblyInput, DraftEnvironment,
    DraftGeometrySource, DraftInterfaceDecl, DraftInterfaceKind, DraftPartDecl, DraftPartKind,
    FidelityErrorBudget, FidelityPolicy, FidelityRung, GearGenerator, GearMeshDecl, JointDecl,
    JointKind, PhysicsDomain, PrimitiveShape, SourceDecl, SourceTarget, StepPolicy, ThirdBodyDecl,
    TimeHistory, TimePlan,
};
use fs_scenario::assembly::{MaterialRef, ParameterKey, ReferenceState, UnknownProperty};
use fs_scenario::{
    EntityCatalog, EntityDeclaration, Environment, FrameId, InterfacePair, Scenario,
};
use fs_system::elements::{
    Admittance, Damper, Force, HeatCapacity, Mass, Motion, SourceAccount, Spring, Storage, Waveform,
};
use fs_system::{
    AcousticAssemblySpec, AcousticMemberSpec, CardStore, FidelityRungKind, FidelitySelectionPolicy,
    InterfaceKindTag, InterfaceLoweringKey, LoweredAcousticMemberOutput, LoweredInterfaceOutput,
    LoweredPartOutput, LoweringError, LoweringRegistry, ParameterProvenance, ParameterSource,
    PartLoweringKey, ProvenanceRecorder, RouteTemperature, SystemBuilder, heat_sink,
};

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const FREQUENCY: Dims = Dims([0, 0, -1, 0, 0, 0]);
const FORCE: Dims = Dims([1, 1, -2, 0, 0, 0]);
const POWER: Dims = Dims([2, 1, -3, 0, 0, 0]);
const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const SPECIFIC_HEAT: Dims = Dims([2, 0, -2, -1, 0, 0]);
const MODULUS: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const THERMAL_COND: Dims = Dims([1, 1, -3, -1, 0, 0]);
const CONDUCTANCE: Dims = Dims([2, 1, -3, -1, 0, 0]);
const STIFFNESS: Dims = Dims([0, 1, -2, 0, 0, 0]);
const TORSIONAL_STIFFNESS: Dims = Dims([2, 1, -2, 0, 0, 0]);
const DAMPING: Dims = Dims([0, 1, -1, 0, 0, 0]);
const TORSIONAL_DAMPING: Dims = Dims([2, 1, -1, 0, 0, 0]);
const AREA: Dims = Dims([2, 0, 0, 0, 0, 0]);
const INDUCTANCE: Dims = Dims([2, 1, -2, 0, -2, 0]);
const CAPACITANCE: Dims = Dims([-2, -1, 4, 0, 2, 0]);
const ELEC_CONDUCTANCE: Dims = Dims([-2, -1, 3, 0, 2, 0]);

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 112,
            kernel_id: 1,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

fn prov(source: &str) -> Provenance {
    Provenance {
        source: source.to_owned(),
        license: "CC-BY-4.0".to_owned(),
        artifact: Some(hash_bytes(source.as_bytes())),
    }
}

fn scalar_claim(name: &str, value: f64, dims: Dims, source: &str) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with("T", 200.0, 600.0),
        uncertainty: UncertaintyModel::HalfWidth {
            half_width: (value.abs() * 0.01).max(1.0e-6),
            confidence: 0.95,
        },
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: prov(source),
    }
}

fn steel_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(scalar_claim("density", 8000.0, DENSITY, "ASM-304-density"))
        .unwrap();
    claims
        .insert_claim(scalar_claim(
            "specific_heat_capacity",
            500.0,
            SPECIFIC_HEAT,
            "ASM-304-cp",
        ))
        .unwrap();
    claims
        .insert_claim(scalar_claim(
            "youngs_modulus",
            193.0e9,
            MODULUS,
            "ASM-304-E",
        ))
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "AISI-304".to_owned(),
            phase: "austenitic".to_owned(),
            process: "annealed".to_owned(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn aluminum_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(scalar_claim("density", 2700.0, DENSITY, "ASM-6061-density"))
        .unwrap();
    claims
        .insert_claim(scalar_claim(
            "specific_heat_capacity",
            896.0,
            SPECIFIC_HEAT,
            "ASM-6061-cp",
        ))
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "AA6061".to_owned(),
            phase: "wrought".to_owned(),
            process: "T6".to_owned(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn tim_paste_card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(scalar_claim(
            "thermal_conductivity",
            4.0,
            THERMAL_COND,
            "TIM-silicone-k",
        ))
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "silicone-alumina".to_owned(),
            phase: "paste".to_owned(),
            process: "dispensed".to_owned(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn interface_card(claims: ClaimSet, medium: &str) -> InterfaceSystemCard {
    let surface_a = SurfaceSpec {
        material: MaterialStateId {
            chemistry: "AISI-304".to_owned(),
            phase: "austenitic".to_owned(),
            process: "annealed".to_owned(),
            revision: 0,
        },
        texture_frame: "ground-Ra0.4".to_owned(),
    };
    let surface_b = SurfaceSpec {
        material: MaterialStateId {
            chemistry: "AA6061".to_owned(),
            phase: "wrought".to_owned(),
            process: "T6".to_owned(),
            revision: 0,
        },
        texture_frame: "milled-Ra0.8".to_owned(),
    };
    InterfaceSystemCard::assemble(
        surface_a,
        surface_b,
        SystemContext {
            medium: medium.to_owned(),
            third_body: None,
            environment: "air".to_owned(),
            history: "run-in".to_owned(),
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn mat_ref(card: ContentHash) -> MaterialRef {
    MaterialRef {
        card,
        query_point: vec![(
            ParameterKey::new("T", TEMPERATURE),
            QtyAny::new(293.15, TEMPERATURE),
        )],
        frame: FrameId(0),
        layer: None,
        unknown: Vec::new(),
    }
}

fn policy_rungs(rungs: Vec<FidelityRung>) -> FidelityPolicy {
    FidelityPolicy {
        rungs,
        error_budget: FidelityErrorBudget {
            relative_tol: 1.0e-5,
            absolute_tol: QtyAny::dimensionless(1.0e-12),
            max_work_units: Some(1_000_000),
        },
    }
}

fn fixed_time_plan(horizon_s: f64, dt_s: f64) -> TimePlan {
    TimePlan {
        horizon: QtyAny::new(horizon_s, TIME),
        step_policy: StepPolicy::Fixed {
            dt: QtyAny::new(dt_s, TIME),
        },
        default_sample_rate: QtyAny::new(1.0 / dt_s, FREQUENCY),
    }
}

fn base_scenario(name: &str) -> Scenario {
    Scenario::new(name, 112, Environment::earth_lab())
}

#[test]
#[allow(clippy::too_many_lines)]
fn slice_a_draft_lowers_to_exact_hand_built_graph_match() {
    let mut cards = CardStore::new();
    let steel_hash = cards.insert_material(steel_card());

    let mut catalog = EntityCatalog::new();
    let root = catalog
        .declare(EntityDeclaration::assembly("slice_a_rig"))
        .unwrap();
    let slider = catalog
        .declare(EntityDeclaration::part(root, "slider"))
        .unwrap();

    // Box 0.1 x 0.05 x 0.05 m => V = 2.5e-4 m^3.
    // With rho = 8000 kg/m^3 => m = 2.0 kg; cp = 500 J/(kg K) => C = 1000 J/K.
    let draft_input = DraftAssemblyInput {
        name: "slice_a_rigid_damper_heat".to_owned(),
        root,
        catalog,
        reference_state: Some(ReferenceState::Unloaded),
        environment: DraftEnvironment::earth_lab_air(),
        time_plan: fixed_time_plan(0.1, 1.0e-3),
        parts: vec![DraftPartDecl {
            name: "slider".to_owned(),
            entity: slider,
            kind: DraftPartKind::RigidBody,
            geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                lx: QtyAny::new(0.1, LENGTH),
                ly: QtyAny::new(0.05, LENGTH),
                lz: QtyAny::new(0.05, LENGTH),
            }),
            materials: vec![mat_ref(steel_hash)],
            physics: BTreeSet::from([PhysicsDomain::Mechanics, PhysicsDomain::Thermal]),
            fidelity: policy_rungs(vec![FidelityRung::Lumped]),
            uncertainty: Vec::new(),
        }],
        joints: vec![JointDecl {
            name: "suspension".to_owned(),
            parent: root,
            child: slider,
            frame: FrameId(0),
            kind: JointKind::Prismatic,
            axis: [1.0, 0.0, 0.0],
            pitch: None,
            limits: None,
            damping: Some(QtyAny::new(12.0, DAMPING)),
            stiffness: Some(QtyAny::new(400.0, STIFFNESS)),
        }],
        gear_meshes: Vec::new(),
        interfaces: Vec::new(),
        sources: vec![SourceDecl {
            name: "drive_force".to_owned(),
            target: SourceTarget::Port {
                entity: slider,
                port: "mechanics".to_owned(),
            },
            history: TimeHistory::Harmonic {
                amplitude: QtyAny::new(15.0, FORCE),
                frequency: QtyAny::new(3.0, FREQUENCY),
                phase: QtyAny::dimensionless(0.25),
                offset: QtyAny::new(5.0, FORCE),
            },
        }],
        controllers: Vec::new(),
        probes: Vec::new(),
    };

    let scenario = base_scenario("slice_a_rig");
    let draft = DraftAssembly::admit(draft_input, &scenario).expect("slice A draft admits");

    let registry = LoweringRegistry::with_kernel_defaults();
    let mut lowered = registry
        .lower(&draft, &cards)
        .expect("slice A draft lowers cleanly");

    // Verify planner ran and receipt is complete and deterministic.
    assert!(lowered.plan.block_count >= 1);
    assert_eq!(lowered.receipt.components.len(), 1);
    assert_eq!(lowered.receipt.joints.len(), 1);
    assert_eq!(lowered.receipt.sources.len(), 1);
    assert_eq!(lowered.receipt.total_property_receipts(), 2); // density + cp
    lowered
        .receipt
        .verify_completeness()
        .expect("receipt complete");
    let lowered_again = registry.lower(&draft, &cards).unwrap();
    assert_eq!(
        lowered.receipt.receipt_hash,
        lowered_again.receipt.receipt_hash
    );

    // Hand-build the exact equivalent kernel graph for step-for-step comparison.
    let mut hb = SystemBuilder::new();
    let t0 = draft.input().environment.base.ambient_temperature.value;
    let m_node = hb
        .add_with_tangent(
            "slider.mechanics",
            Mass::new(2.0).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let c_node = hb
        .add_with_tangent("slider.thermal", HeatCapacity::new(1000.0).unwrap(), t0)
        .unwrap();
    let s_node = hb
        .add_with_tangent(
            "joint.suspension.spring",
            Spring::new(400.0, 0.0).unwrap(),
            0.0,
        )
        .unwrap();
    let d_node = hb
        .add_with_tangent("joint.suspension.damper", Damper::new(12.0).unwrap(), ())
        .unwrap();
    let f_node = hb
        .add_with_tangent(
            "source.drive_force",
            Force::new(Waveform::sine(15.0, 3.0, 0.25, 5.0).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();

    let m_p = hb.port(m_node, "p").unwrap();
    let c_heat = hb.port(c_node, "heat").unwrap();
    let (sa, sb) = (hb.port(s_node, "a").unwrap(), hb.port(s_node, "b").unwrap());
    let (da, db, dh) = (
        hb.port(d_node, "a").unwrap(),
        hb.port(d_node, "b").unwrap(),
        hb.port(d_node, "heat").unwrap(),
    );
    let f_p = hb.port(f_node, "p").unwrap();
    hb.common_flow(m_p, &[sa, da, f_p]).unwrap();
    hb.common_effort(c_heat, &[dh]).unwrap();
    hb.boundary(sb, Box::new(|_, _, out| out[0] = 0.0)).unwrap();
    hb.boundary(db, Box::new(|_, _, out| out[0] = 0.0)).unwrap();
    let mut hand_built = hb.seal().expect("hand-built system seals");
    assert!(hand_built.guess(dh, &[t0]));

    let gate = CancelGate::new_clock_free();
    let low_m = lowered.node_map["slider.mechanics"];
    let low_c = lowered.node_map["slider.thermal"];
    with_cx(&gate, |cx| {
        for step_idx in 0..100 {
            let _ = lowered.system.step(1.0e-3, cx).expect("lowered step");
            let _ = hand_built.step(1.0e-3, cx).expect("hand-built step");
            let m_low: &Motion = lowered.system.state(low_m).unwrap();
            let m_hb: &Motion = hand_built.state(m_node).unwrap();
            assert!(
                (m_low.position_m - m_hb.position_m).abs() <= 1.0e-14,
                "step {step_idx} position mismatch"
            );
            assert!(
                (m_low.velocity_m_s - m_hb.velocity_m_s).abs() <= 1.0e-14,
                "step {step_idx} velocity mismatch"
            );
            let t_low: &f64 = lowered.system.state(low_c).unwrap();
            let t_hb: &f64 = hand_built.state(c_node).unwrap();
            assert!(
                (t_low - t_hb).abs() <= 1.0e-13,
                "step {step_idx} temperature mismatch"
            );
        }
        assert!(
            (lowered.system.ledger().total_dissipated_j()
                - hand_built.ledger().total_dissipated_j())
            .abs()
                <= 1.0e-13
        );
        assert!(
            (lowered.system.ledger().entropy_generated_j_per_k()
                - hand_built.ledger().entropy_generated_j_per_k())
            .abs()
                <= 1.0e-15
        );
        assert!(lowered.system.ledger().energy_residual_j().abs() <= 1.0e-11);
    });
    let t_final: &f64 = lowered.system.state(low_c).unwrap();
    assert!(
        *t_final > t0,
        "dissipated mechanical energy heated the slider"
    );

    // Every history must retain its units when it becomes a force or boundary.
    let kelvin = QtyAny::new(300.0, TEMPERATURE);
    for (port, history) in [
        ("mechanics", TimeHistory::Constant(kelvin)),
        (
            "mechanics",
            TimeHistory::Harmonic {
                amplitude: kelvin,
                frequency: QtyAny::new(3.0, FREQUENCY),
                phase: QtyAny::dimensionless(0.0),
                offset: kelvin,
            },
        ),
        (
            "mechanics",
            TimeHistory::PiecewiseLinear {
                knots: vec![
                    (QtyAny::new(0.0, TIME), kelvin),
                    (QtyAny::new(0.1, TIME), kelvin),
                ],
            },
        ),
        (
            "mechanics",
            TimeHistory::Step {
                t_switch: QtyAny::new(0.05, TIME),
                before: kelvin,
                after: kelvin,
            },
        ),
        ("thermal", TimeHistory::Constant(kelvin)),
    ] {
        let mut wrong_units = draft.input().clone();
        wrong_units.sources[0].target = SourceTarget::Port {
            entity: slider,
            port: port.to_owned(),
        };
        wrong_units.sources[0].history = history;
        let wrong_units = DraftAssembly::admit(wrong_units, &scenario).unwrap();
        let error = registry
            .lower(&wrong_units, &cards)
            .err()
            .expect("wrong source units refuse");
        assert!(matches!(&error, LoweringError::Invalid { item, detail }
            if item == "source:drive_force" && detail.contains("dimensions") && detail.contains(port)));
    }

    let mut controlled = draft.input().clone();
    controlled.controllers.push(ControllerDecl {
        name: "velocity_feedback".to_owned(),
        sensor_ports: vec![(slider, "velocity".to_owned())],
        actuator_port: (slider, "command".to_owned()),
        sample_rate: QtyAny::new(1000.0, FREQUENCY),
        law: ControllerLaw::Pid {
            kp: QtyAny::dimensionless(1.0),
            ki: QtyAny::new(0.0, FREQUENCY),
            kd: QtyAny::new(0.0, TIME),
            setpoint: QtyAny::dimensionless(1.0),
        },
    });
    let controlled = DraftAssembly::admit(controlled, &scenario).unwrap();
    let error = registry
        .lower(&controlled, &cards)
        .err()
        .expect("feedback cannot silently disappear");
    assert!(matches!(error, LoweringError::Invalid { ref item, .. }
        if item == "controller:velocity_feedback"));
}

#[test]
#[allow(clippy::too_many_lines)]
fn unseen_multi_domain_device_lowers_with_zero_rust_changes() {
    let mut cards = CardStore::new();
    let steel_hash = cards.insert_material(steel_card());
    let alu_hash = cards.insert_material(aluminum_card());
    let tim_hash = cards.insert_material(tim_paste_card());

    let mut kv_claims = ClaimSet::new();
    kv_claims
        .insert_claim(scalar_claim(
            "contact_stiffness",
            250.0,
            STIFFNESS,
            "KV-contact-k",
        ))
        .unwrap();
    kv_claims
        .insert_claim(scalar_claim(
            "damping_coefficient",
            8.0,
            DAMPING,
            "KV-contact-c",
        ))
        .unwrap();
    let kv_hash = cards
        .insert_interface("kelvin-voigt", interface_card(kv_claims, "elastomer-pad"))
        .unwrap();

    let mut therm_claims = ClaimSet::new();
    therm_claims
        .insert_claim(scalar_claim(
            "thermal_conductance",
            15.0,
            CONDUCTANCE,
            "TIM-contact-G",
        ))
        .unwrap();
    let therm_hash = cards
        .insert_interface("thermal-conductance", interface_card(therm_claims, "tim"))
        .unwrap();

    let mut mesh_claims = ClaimSet::new();
    mesh_claims
        .insert_claim(scalar_claim(
            "mesh_stiffness",
            500.0,
            TORSIONAL_STIFFNESS,
            "AGMA-gear-mesh-k",
        ))
        .unwrap();
    let mesh_hash = cards
        .insert_interface(
            "involute-gear-mesh",
            interface_card(mesh_claims, "oil-film"),
        )
        .unwrap();

    let mut catalog = EntityCatalog::new();
    let root = catalog
        .declare(EntityDeclaration::assembly("gear_and_thermal_bench"))
        .unwrap();
    let pinion = catalog
        .declare(EntityDeclaration::part(root, "pinion"))
        .unwrap();
    let wheel = catalog
        .declare(EntityDeclaration::part(root, "wheel"))
        .unwrap();
    let carriage = catalog
        .declare(EntityDeclaration::part(root, "carriage"))
        .unwrap();
    let heatsink = catalog
        .declare(EntityDeclaration::part(root, "heatsink"))
        .unwrap();
    let frame = catalog
        .declare(EntityDeclaration::part(root, "bench_frame"))
        .unwrap();
    let ground_surf = catalog
        .declare(EntityDeclaration::surface(frame, "ground_pad"))
        .unwrap();
    let carriage_stop_surf = catalog
        .declare(EntityDeclaration::surface(carriage, "carriage_stop_face"))
        .unwrap();
    let carriage_tim_surf = catalog
        .declare(EntityDeclaration::surface(carriage, "carriage_tim_face"))
        .unwrap();
    let heatsink_tim_surf = catalog
        .declare(EntityDeclaration::surface(heatsink, "heatsink_tim_face"))
        .unwrap();
    let kv_iface = catalog
        .declare(EntityDeclaration::interface(
            root,
            "carriage_stop",
            InterfacePair::ordered(carriage_stop_surf, ground_surf),
        ))
        .unwrap();
    let tim_iface = catalog
        .declare(EntityDeclaration::interface(
            root,
            "carriage_to_sink",
            InterfacePair::ordered(carriage_tim_surf, heatsink_tim_surf),
        ))
        .unwrap();

    // Unseen composite device combining:
    // - 2 meshed gears (pinion 20T -> wheel 40T) with revolute bearings;
    // - thermo-mechanical carriage with multi-rung ladder [FiniteElement, Lumped];
    // - Kelvin-Voigt contact + TIM thermal conductance (with third-body paste!)
    //   into an aluminum heatsink;
    // - piecewise-linear force + constant Joule heater.
    let draft_input = DraftAssemblyInput {
        name: "unseen_gear_thermal_rig".to_owned(),
        root,
        catalog,
        reference_state: Some(ReferenceState::Unloaded),
        environment: DraftEnvironment::earth_lab_air(),
        time_plan: fixed_time_plan(0.05, 1.0e-3),
        parts: vec![
            DraftPartDecl {
                name: "pinion".to_owned(),
                entity: pinion,
                kind: DraftPartKind::RigidBody,
                geometry: DraftGeometrySource::Gear(GearGenerator {
                    teeth: 20,
                    module: QtyAny::new(0.002, LENGTH),
                    face_width: QtyAny::new(0.015, LENGTH),
                    pressure_angle: QtyAny::dimensionless(20.0_f64.to_radians()),
                    profile_shift: 0.0,
                    bore_radius: QtyAny::new(0.004, LENGTH),
                }),
                materials: vec![mat_ref(steel_hash)],
                physics: BTreeSet::from([PhysicsDomain::Mechanics]),
                fidelity: policy_rungs(vec![FidelityRung::Lumped]),
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "wheel".to_owned(),
                entity: wheel,
                kind: DraftPartKind::RigidBody,
                geometry: DraftGeometrySource::Gear(GearGenerator {
                    teeth: 40,
                    module: QtyAny::new(0.002, LENGTH),
                    face_width: QtyAny::new(0.015, LENGTH),
                    pressure_angle: QtyAny::dimensionless(20.0_f64.to_radians()),
                    profile_shift: 0.0,
                    bore_radius: QtyAny::new(0.006, LENGTH),
                }),
                materials: vec![mat_ref(steel_hash)],
                physics: BTreeSet::from([PhysicsDomain::Mechanics]),
                fidelity: policy_rungs(vec![FidelityRung::Lumped]),
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "carriage".to_owned(),
                entity: carriage,
                kind: DraftPartKind::RigidBody,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Cylinder {
                    radius: QtyAny::new(0.03, LENGTH),
                    height: QtyAny::new(0.08, LENGTH),
                }),
                materials: vec![mat_ref(steel_hash)],
                physics: BTreeSet::from([PhysicsDomain::Mechanics, PhysicsDomain::Thermal]),
                fidelity: policy_rungs(vec![
                    FidelityRung::FiniteElement {
                        order: 2,
                        h_target: QtyAny::new(0.005, LENGTH),
                    },
                    FidelityRung::Lumped,
                ]),
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "heatsink".to_owned(),
                entity: heatsink,
                kind: DraftPartKind::Solid,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                    lx: QtyAny::new(0.12, LENGTH),
                    ly: QtyAny::new(0.08, LENGTH),
                    lz: QtyAny::new(0.02, LENGTH),
                }),
                materials: vec![mat_ref(alu_hash)],
                physics: BTreeSet::from([PhysicsDomain::Thermal]),
                fidelity: policy_rungs(vec![FidelityRung::Lumped]),
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "bench_frame".to_owned(),
                entity: frame,
                kind: DraftPartKind::DrivenBody,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                    lx: QtyAny::new(0.5, LENGTH),
                    ly: QtyAny::new(0.5, LENGTH),
                    lz: QtyAny::new(0.05, LENGTH),
                }),
                materials: vec![mat_ref(steel_hash)],
                physics: BTreeSet::from([PhysicsDomain::Mechanics]),
                fidelity: policy_rungs(vec![FidelityRung::Lumped]),
                uncertainty: Vec::new(),
            },
        ],
        joints: vec![
            JointDecl {
                name: "pinion_bearing".to_owned(),
                parent: root,
                child: pinion,
                frame: FrameId(0),
                kind: JointKind::Revolute,
                axis: [0.0, 0.0, 1.0],
                pitch: None,
                limits: None,
                damping: Some(QtyAny::new(0.02, TORSIONAL_DAMPING)),
                stiffness: None,
            },
            JointDecl {
                name: "wheel_bearing".to_owned(),
                parent: root,
                child: wheel,
                frame: FrameId(0),
                kind: JointKind::Revolute,
                axis: [0.0, 0.0, 1.0],
                pitch: None,
                limits: None,
                damping: Some(QtyAny::new(0.05, TORSIONAL_DAMPING)),
                stiffness: None,
            },
        ],
        gear_meshes: vec![GearMeshDecl {
            name: "stage_1_mesh".to_owned(),
            driver: pinion,
            driven: wheel,
            teeth_driver: 20,
            teeth_driven: 40,
            pressure_angle: QtyAny::dimensionless(20.0_f64.to_radians()),
            center_distance: QtyAny::new(0.06, LENGTH),
            backlash: None,
            mesh_card: mesh_hash,
        }],
        interfaces: vec![
            DraftInterfaceDecl {
                name: "carriage_stop".to_owned(),
                entity: kv_iface,
                feature_a: carriage_stop_surf,
                feature_b: ground_surf,
                kind: DraftInterfaceKind::MechanicalContact {
                    law_card: kv_hash,
                    third_body: None,
                    state: Vec::new(),
                },
            },
            DraftInterfaceDecl {
                name: "carriage_to_sink".to_owned(),
                entity: tim_iface,
                feature_a: carriage_tim_surf,
                feature_b: heatsink_tim_surf,
                kind: DraftInterfaceKind::ThermalContact {
                    law_card: therm_hash,
                    third_body: Some(ThirdBodyDecl {
                        material_card: tim_hash,
                        film_thickness: QtyAny::new(0.0005, LENGTH),
                        query_point: vec![(
                            ParameterKey::new("T", TEMPERATURE),
                            QtyAny::new(293.15, TEMPERATURE),
                        )],
                    }),
                    state: vec![(
                        ParameterKey::new("contact_area", AREA),
                        QtyAny::new(0.002, AREA),
                    )],
                },
            },
        ],
        sources: vec![
            SourceDecl {
                name: "carriage_pulse".to_owned(),
                target: SourceTarget::Port {
                    entity: carriage,
                    port: "mechanics".to_owned(),
                },
                history: TimeHistory::PiecewiseLinear {
                    knots: vec![
                        (QtyAny::new(0.0, TIME), QtyAny::new(0.0, FORCE)),
                        (QtyAny::new(0.01, TIME), QtyAny::new(20.0, FORCE)),
                        (QtyAny::new(0.05, TIME), QtyAny::new(5.0, FORCE)),
                    ],
                },
            },
            SourceDecl {
                name: "sink_heater".to_owned(),
                target: SourceTarget::Port {
                    entity: heatsink,
                    port: "thermal".to_owned(),
                },
                history: TimeHistory::Constant(QtyAny::new(25.0, POWER)),
            },
        ],
        controllers: Vec::new(),
        probes: Vec::new(),
    };

    let scenario = base_scenario("unseen_gear_thermal_rig");
    let draft = DraftAssembly::admit(draft_input, &scenario).expect("unseen draft admits");
    let registry = LoweringRegistry::with_kernel_defaults()
        .with_fidelity_policy(FidelitySelectionPolicy::FirstRegistered);
    let mut lowered = registry
        .lower(&draft, &cards)
        .expect("unseen multi-domain device lowers with zero Rust changes");

    // Verify carriage skipped unregistered FiniteElement and selected Lumped with explicit rationale.
    let carriage_rec = lowered
        .receipt
        .components
        .iter()
        .find(|c| c.part_name == "carriage")
        .unwrap();
    assert_eq!(carriage_rec.chosen_fidelity, FidelityRung::Lumped);
    assert!(
        carriage_rec.fidelity_rationale.contains("ladder[1]"),
        "expected ladder[1] rationale: {}",
        carriage_rec.fidelity_rationale
    );
    assert!(carriage_rec.fidelity_rationale.contains("finite_element"));

    // Verify TIM third-body receipt was recorded on carriage_to_sink interface.
    let tim_rec = lowered
        .receipt
        .interfaces
        .iter()
        .find(|i| i.interface_name == "carriage_to_sink")
        .unwrap();
    assert_eq!(tim_rec.third_body_card, Some(tim_hash));
    assert_eq!(tim_rec.property_receipts.len(), 2); // thermal_conductance + third-body thermal_conductivity

    // Verify gear ratio = -20/40 = -0.5.
    assert_eq!(lowered.receipt.gear_meshes.len(), 1);
    assert!((lowered.receipt.gear_meshes[0].ratio - (-0.5)).abs() <= 1.0e-15);

    // Step the lowered multi-domain system and confirm energy + entropy closure.
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for _ in 0..50 {
            let _ = lowered.system.step(1.0e-3, cx).expect("step succeeds");
        }
    });
    assert!(lowered.system.ledger().entropy_generated_j_per_k() >= 0.0);
    assert!(lowered.system.ledger().energy_residual_j().abs() <= 1.0e-10);
    let json = lowered.receipt.canonical_json();
    assert!(json.contains("\"schema\":\"frankensim.lowering_receipt.v1\""));
    assert!(json.contains("\"third_body_thermal_conductivity\""));
}

#[test]
#[allow(clippy::too_many_lines)]
fn unknown_combinations_refuse_by_name_stating_missing_constructor() {
    let mut cards = CardStore::new();
    let steel_hash = cards.insert_material(steel_card());
    let registry = LoweringRegistry::with_kernel_defaults();

    // 1. Missing part constructor.
    let mut catalog = EntityCatalog::new();
    let root = catalog
        .declare(EntityDeclaration::assembly("bad_part_rig"))
        .unwrap();
    let flex = catalog
        .declare(EntityDeclaration::part(root, "exotic_coil_flex"))
        .unwrap();
    let scenario = base_scenario("bad_part_rig");
    let draft = DraftAssembly::admit(
        DraftAssemblyInput {
            name: "bad_part_rig".to_owned(),
            root,
            catalog,
            reference_state: Some(ReferenceState::Unloaded),
            environment: DraftEnvironment::earth_lab_air(),
            time_plan: fixed_time_plan(0.01, 1.0e-3),
            parts: vec![DraftPartDecl {
                name: "exotic_coil_flex".to_owned(),
                entity: flex,
                kind: DraftPartKind::FlexibleBody,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Sphere {
                    radius: QtyAny::new(0.02, LENGTH),
                }),
                materials: vec![mat_ref(steel_hash)],
                physics: BTreeSet::from([
                    PhysicsDomain::Mechanics,
                    PhysicsDomain::Electromagnetics,
                ]),
                fidelity: policy_rungs(vec![
                    FidelityRung::FiniteElement {
                        order: 2,
                        h_target: QtyAny::new(0.002, LENGTH),
                    },
                    FidelityRung::Modal { mode_count: 8 },
                ]),
                uncertainty: Vec::new(),
            }],
            joints: Vec::new(),
            gear_meshes: Vec::new(),
            interfaces: Vec::new(),
            sources: Vec::new(),
            controllers: Vec::new(),
            probes: Vec::new(),
        },
        &scenario,
    )
    .unwrap();
    let err = registry.lower(&draft, &cards).err().expect("must refuse");
    let msg = err.to_string();
    assert!(matches!(err, LoweringError::MissingPartConstructor { .. }));
    assert!(msg.contains("exotic_coil_flex"), "{msg}");
    assert!(msg.contains("flexible_body"), "{msg}");
    assert!(msg.contains("electromagnetics"), "{msg}");
    assert!(msg.contains("finite_element"), "{msg}");
    assert!(msg.contains("modal"), "{msg}");

    // 2. Missing joint constructor (Helical).
    let mut catalog = EntityCatalog::new();
    let root = catalog
        .declare(EntityDeclaration::assembly("helical_rig"))
        .unwrap();
    let nut = catalog
        .declare(EntityDeclaration::part(root, "lead_nut"))
        .unwrap();
    let draft_joint = DraftAssembly::admit(
        DraftAssemblyInput {
            name: "helical_rig".to_owned(),
            root,
            catalog,
            reference_state: Some(ReferenceState::Unloaded),
            environment: DraftEnvironment::earth_lab_air(),
            time_plan: fixed_time_plan(0.01, 1.0e-3),
            parts: vec![DraftPartDecl {
                name: "lead_nut".to_owned(),
                entity: nut,
                kind: DraftPartKind::RigidBody,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                    lx: QtyAny::new(0.05, LENGTH),
                    ly: QtyAny::new(0.05, LENGTH),
                    lz: QtyAny::new(0.05, LENGTH),
                }),
                materials: vec![mat_ref(steel_hash)],
                physics: BTreeSet::from([PhysicsDomain::Mechanics]),
                fidelity: policy_rungs(vec![FidelityRung::Lumped]),
                uncertainty: Vec::new(),
            }],
            joints: vec![JointDecl {
                name: "lead_screw".to_owned(),
                parent: root,
                child: nut,
                frame: FrameId(0),
                kind: JointKind::Helical,
                axis: [0.0, 0.0, 1.0],
                pitch: Some(QtyAny::new(0.002, LENGTH)),
                limits: None,
                damping: None,
                stiffness: None,
            }],
            gear_meshes: Vec::new(),
            interfaces: Vec::new(),
            sources: Vec::new(),
            controllers: Vec::new(),
            probes: Vec::new(),
        },
        &scenario,
    )
    .unwrap();
    let err = registry
        .lower(&draft_joint, &cards)
        .err()
        .expect("helical joint must refuse");
    assert!(matches!(err, LoweringError::MissingJointConstructor { .. }));
    assert!(err.to_string().contains("lead_screw"));
    assert!(err.to_string().contains("joint(helical)"));

    // 3. Missing interface constructor (unknown card family).
    let mut exotic_claims = ClaimSet::new();
    exotic_claims
        .insert_claim(scalar_claim(
            "conductance",
            1.0,
            ELEC_CONDUCTANCE,
            "exotic-contact",
        ))
        .unwrap();
    let exotic_hash = cards
        .insert_interface(
            "magneto-strictive-gap",
            interface_card(exotic_claims, "vacuum"),
        )
        .unwrap();
    let mut catalog = EntityCatalog::new();
    let root = catalog
        .declare(EntityDeclaration::assembly("iface_rig"))
        .unwrap();
    let body_a = catalog
        .declare(EntityDeclaration::part(root, "body_a"))
        .unwrap();
    let body_b = catalog
        .declare(EntityDeclaration::part(root, "body_b"))
        .unwrap();
    let surf_a = catalog
        .declare(EntityDeclaration::surface(body_a, "surf_a"))
        .unwrap();
    let surf_b = catalog
        .declare(EntityDeclaration::surface(body_b, "surf_b"))
        .unwrap();
    let iface_id = catalog
        .declare(EntityDeclaration::interface(
            root,
            "mystery_gap",
            InterfacePair::ordered(surf_a, surf_b),
        ))
        .unwrap();
    let draft_iface = DraftAssembly::admit(
        DraftAssemblyInput {
            name: "iface_rig".to_owned(),
            root,
            catalog,
            reference_state: Some(ReferenceState::Unloaded),
            environment: DraftEnvironment::earth_lab_air(),
            time_plan: fixed_time_plan(0.01, 1.0e-3),
            parts: vec![
                DraftPartDecl {
                    name: "body_a".to_owned(),
                    entity: body_a,
                    kind: DraftPartKind::RigidBody,
                    geometry: DraftGeometrySource::Primitive(PrimitiveShape::Sphere {
                        radius: QtyAny::new(0.02, LENGTH),
                    }),
                    materials: vec![mat_ref(steel_hash)],
                    physics: BTreeSet::from([PhysicsDomain::Mechanics]),
                    fidelity: policy_rungs(vec![FidelityRung::Lumped]),
                    uncertainty: Vec::new(),
                },
                DraftPartDecl {
                    name: "body_b".to_owned(),
                    entity: body_b,
                    kind: DraftPartKind::RigidBody,
                    geometry: DraftGeometrySource::Primitive(PrimitiveShape::Sphere {
                        radius: QtyAny::new(0.02, LENGTH),
                    }),
                    materials: vec![mat_ref(steel_hash)],
                    physics: BTreeSet::from([PhysicsDomain::Mechanics]),
                    fidelity: policy_rungs(vec![FidelityRung::Lumped]),
                    uncertainty: Vec::new(),
                },
            ],
            joints: Vec::new(),
            gear_meshes: Vec::new(),
            interfaces: vec![DraftInterfaceDecl {
                name: "mystery_gap".to_owned(),
                entity: iface_id,
                feature_a: surf_a,
                feature_b: surf_b,
                kind: DraftInterfaceKind::ElectricalContact {
                    law_card: exotic_hash,
                    third_body: None,
                    state: Vec::new(),
                },
            }],
            sources: Vec::new(),
            controllers: Vec::new(),
            probes: Vec::new(),
        },
        &scenario,
    )
    .unwrap();
    let err = registry
        .lower(&draft_iface, &cards)
        .err()
        .expect("unknown interface family must refuse");
    assert!(matches!(
        err,
        LoweringError::MissingInterfaceConstructor { .. }
    ));
    assert!(err.to_string().contains("mystery_gap"));
    assert!(err.to_string().contains("magneto-strictive-gap"));

    // 4. Explicit UnknownProperty refuses instead of defaulting.
    let mut catalog = EntityCatalog::new();
    let root = catalog
        .declare(EntityDeclaration::assembly("unknown_prop_rig"))
        .unwrap();
    let body = catalog
        .declare(EntityDeclaration::part(root, "unknown_body"))
        .unwrap();
    let mut mref = mat_ref(steel_hash);
    mref.unknown.push(UnknownProperty::Unknown {
        property: ParameterKey::new("density", DENSITY),
        bound: None,
    });
    let draft_unknown = DraftAssembly::admit(
        DraftAssemblyInput {
            name: "unknown_prop_rig".to_owned(),
            root,
            catalog,
            reference_state: Some(ReferenceState::Unloaded),
            environment: DraftEnvironment::earth_lab_air(),
            time_plan: fixed_time_plan(0.01, 1.0e-3),
            parts: vec![DraftPartDecl {
                name: "unknown_body".to_owned(),
                entity: body,
                kind: DraftPartKind::RigidBody,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Sphere {
                    radius: QtyAny::new(0.02, LENGTH),
                }),
                materials: vec![mref],
                physics: BTreeSet::from([PhysicsDomain::Mechanics]),
                fidelity: policy_rungs(vec![FidelityRung::Lumped]),
                uncertainty: Vec::new(),
            }],
            joints: Vec::new(),
            gear_meshes: Vec::new(),
            interfaces: Vec::new(),
            sources: Vec::new(),
            controllers: Vec::new(),
            probes: Vec::new(),
        },
        &scenario,
    )
    .unwrap();
    let err = registry
        .lower(&draft_unknown, &cards)
        .err()
        .expect("unknown property must refuse");
    assert!(matches!(
        err,
        LoweringError::UnknownPropertyDefaulted { .. }
    ));
    assert!(err.to_string().contains("density"));
}

#[test]
#[allow(clippy::too_many_lines)]
fn receipt_completeness_refuses_unsourced_parameters() {
    let mut cards = CardStore::new();
    let steel_hash = cards.insert_material(steel_card());

    // Register a buggy custom constructor that omits all parameter provenance.
    let mut registry = LoweringRegistry::new();
    registry.register_part(
        PartLoweringKey::new(
            DraftPartKind::RigidBody,
            [PhysicsDomain::Mechanics],
            FidelityRungKind::Lumped,
        ),
        |ctx, graph| {
            let node = graph.add_with_tangent(
                &ctx.part.name,
                Some(ctx.part.entity),
                Mass::new(1.0).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s: 0.0,
                },
            )?;
            let p = graph.port(node, "p")?;
            graph.bind_zero_boundary(p);
            Ok(LoweredPartOutput {
                created_nodes: vec![ctx.part.name.clone()],
                property_receipts: Vec::new(),
                parameter_sources: Vec::new(), // <-- missing provenance!
            })
        },
    );

    let mut catalog = EntityCatalog::new();
    let root = catalog
        .declare(EntityDeclaration::assembly("unsourced_rig"))
        .unwrap();
    let body = catalog
        .declare(EntityDeclaration::part(root, "sneaky_mass"))
        .unwrap();
    let scenario = base_scenario("unsourced_rig");
    let draft = DraftAssembly::admit(
        DraftAssemblyInput {
            name: "unsourced_rig".to_owned(),
            root,
            catalog,
            reference_state: Some(ReferenceState::Unloaded),
            environment: DraftEnvironment::earth_lab_air(),
            time_plan: fixed_time_plan(0.01, 1.0e-3),
            parts: vec![DraftPartDecl {
                name: "sneaky_mass".to_owned(),
                entity: body,
                kind: DraftPartKind::RigidBody,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Sphere {
                    radius: QtyAny::new(0.02, LENGTH),
                }),
                materials: vec![mat_ref(steel_hash)],
                physics: BTreeSet::from([PhysicsDomain::Mechanics]),
                fidelity: policy_rungs(vec![FidelityRung::Lumped]),
                uncertainty: Vec::new(),
            }],
            joints: Vec::new(),
            gear_meshes: Vec::new(),
            interfaces: Vec::new(),
            sources: Vec::new(),
            controllers: Vec::new(),
            probes: Vec::new(),
        },
        &scenario,
    )
    .unwrap();

    let err = registry
        .lower(&draft, &cards)
        .err()
        .expect("unsourced parameters must refuse");
    assert!(matches!(err, LoweringError::MissingParameterSource { .. }));

    // Also verify forged MatDbClaim without matching PropertyUsageReceipt refuses.
    let mut forged_registry = LoweringRegistry::new();
    forged_registry.register_part(
        PartLoweringKey::new(
            DraftPartKind::RigidBody,
            [PhysicsDomain::Mechanics],
            FidelityRungKind::Lumped,
        ),
        |ctx, graph| {
            let node = graph.add_with_tangent(
                &ctx.part.name,
                Some(ctx.part.entity),
                Mass::new(1.0).unwrap(),
                Motion {
                    position_m: 0.0,
                    velocity_m_s: 0.0,
                },
            )?;
            let p = graph.port(node, "p")?;
            graph.bind_zero_boundary(p);
            Ok(LoweredPartOutput {
                created_nodes: vec![ctx.part.name.clone()],
                property_receipts: Vec::new(), // <-- missing actual receipt!
                parameter_sources: vec![ParameterProvenance {
                    parameter: "density".to_owned(),
                    value: 8000.0,
                    dims: DENSITY,
                    source: ParameterSource::MatDbClaim {
                        card_hash: hash_bytes(b"fake-card"),
                        property: "density".to_owned(),
                        receipt_hash: hash_bytes(b"fake-receipt"),
                        selected_claim: hash_bytes(b"fake-claim"),
                        observation_backed: false,
                    },
                }],
            })
        },
    );
    let err2 = forged_registry
        .lower(&draft, &cards)
        .err()
        .expect("forged MatDbClaim without matching PropertyUsageReceipt must refuse");
    assert!(matches!(err2, LoweringError::MissingParameterSource { .. }));
    assert!(
        err2.to_string()
            .contains("absent from item property_receipts")
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn owner_extensibility_custom_em_constructors_register_without_kernel_edits() {
    let mut cards = CardStore::new();
    let mut cap_claims = ClaimSet::new();
    cap_claims
        .insert_claim(scalar_claim(
            "capacitance",
            1.0e-3,
            CAPACITANCE,
            "film-cap-spec",
        ))
        .unwrap();
    let cap_card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "polypropylene".to_owned(),
            phase: "film".to_owned(),
            process: "metallized".to_owned(),
            revision: 0,
        },
        cap_claims,
        Vec::new(),
    )
    .unwrap();
    let cap_hash = cards.insert_material(cap_card);

    let mut ohmic_claims = ClaimSet::new();
    ohmic_claims
        .insert_claim(scalar_claim(
            "electrical_conductance",
            0.25,
            ELEC_CONDUCTANCE,
            "ohmic-shunt",
        ))
        .unwrap();
    let ohmic_hash = cards
        .insert_interface("ohmic-contact", interface_card(ohmic_claims, "copper-bus"))
        .unwrap();

    // External physics owner registers its own (Solid, {Electromagnetics}, Lumped)
    // and (ElectricalContact, "ohmic-contact") constructors without modifying fs-system.
    let mut registry = LoweringRegistry::with_kernel_defaults();
    registry.register_part(
        PartLoweringKey::new(
            DraftPartKind::Solid,
            [PhysicsDomain::Electromagnetics],
            FidelityRungKind::Lumped,
        ),
        |ctx, graph| {
            let mut rec = ProvenanceRecorder::default();
            let (c, _) = ctx.query_material(&mut rec, 0, "capacitance", "capacitance")?;
            let node = graph.add_with_tangent(
                &ctx.part.name,
                Some(ctx.part.entity),
                Storage::in_domain(PortKind::ElectricalVoltageCurrent, c.value).unwrap(),
                5.0e-3, // initial charge 5 mC => 5 V across 1 mF capacitor
            )?;
            let p = graph.port(node, "p")?;
            graph.expose_port(ctx.part.entity, "electrical", p);
            let (property_receipts, parameter_sources) = rec.finish();
            Ok(LoweredPartOutput {
                created_nodes: vec![ctx.part.name.clone()],
                property_receipts,
                parameter_sources,
            })
        },
    );
    registry.register_interface(
        InterfaceLoweringKey::new(InterfaceKindTag::ElectricalContact, "ohmic-contact"),
        |ctx, graph| {
            let fs_system::ResolvedInterfaceCard::System(entry) = ctx.card else {
                unreachable!()
            };
            let mut rec = ProvenanceRecorder::default();
            let (g, _) = rec.query_interface_property(
                &format!("interface:{}", ctx.interface.name),
                entry,
                &[],
                "electrical_conductance",
                "electrical_conductance",
                ctx.query_env(),
            )?;
            let node_name = format!("interface.{}", ctx.interface.name);
            let node = graph.add_with_tangent(
                &node_name,
                Some(ctx.interface.entity),
                Admittance::in_domain(PortKind::ElectricalVoltageCurrent, g.value).unwrap(),
                (),
            )?;
            let (a, b, heat) = (
                graph.port(node, "a")?,
                graph.port(node, "b")?,
                graph.port(node, "heat")?,
            );
            let cap_port = graph
                .resolve_entity_port(
                    &ctx.draft.input().catalog,
                    ctx.interface.feature_a,
                    "electrical",
                )
                .unwrap();
            graph.attach_common_effort(cap_port, a, 1.0);
            graph.bind_zero_boundary(b);
            graph.add_heat_route(
                heat,
                Vec::new(),
                Some(heat_sink("ambient", ctx.ambient_temperature_k()).unwrap()),
                RouteTemperature::Fixed(ctx.ambient_temperature_k()),
            );
            let (property_receipts, parameter_sources) = rec.finish();
            Ok(LoweredInterfaceOutput {
                created_nodes: vec![node_name],
                property_receipts,
                parameter_sources,
            })
        },
    );

    let mut catalog = EntityCatalog::new();
    let root = catalog
        .declare(EntityDeclaration::assembly("rc_discharge"))
        .unwrap();
    let cap_entity = catalog
        .declare(EntityDeclaration::part(root, "bus_cap"))
        .unwrap();
    let chassis = catalog
        .declare(EntityDeclaration::part(root, "chassis"))
        .unwrap();
    let cap_term = catalog
        .declare(EntityDeclaration::surface(cap_entity, "bus_cap_terminal"))
        .unwrap();
    let gnd_term = catalog
        .declare(EntityDeclaration::surface(chassis, "chassis_ground"))
        .unwrap();
    let shunt_entity = catalog
        .declare(EntityDeclaration::interface(
            root,
            "bleed_resistor",
            InterfacePair::ordered(cap_term, gnd_term),
        ))
        .unwrap();
    let scenario = base_scenario("rc_discharge");
    let draft = DraftAssembly::admit(
        DraftAssemblyInput {
            name: "rc_discharge".to_owned(),
            root,
            catalog,
            reference_state: Some(ReferenceState::Unloaded),
            environment: DraftEnvironment::earth_lab_air(),
            time_plan: fixed_time_plan(0.01, 1.0e-4),
            parts: vec![
                DraftPartDecl {
                    name: "bus_cap".to_owned(),
                    entity: cap_entity,
                    kind: DraftPartKind::Solid,
                    geometry: DraftGeometrySource::Primitive(PrimitiveShape::Cylinder {
                        radius: QtyAny::new(0.01, LENGTH),
                        height: QtyAny::new(0.02, LENGTH),
                    }),
                    materials: vec![mat_ref(cap_hash)],
                    physics: BTreeSet::from([PhysicsDomain::Electromagnetics]),
                    fidelity: policy_rungs(vec![FidelityRung::Lumped]),
                    uncertainty: Vec::new(),
                },
                DraftPartDecl {
                    name: "chassis".to_owned(),
                    entity: chassis,
                    kind: DraftPartKind::DrivenBody,
                    geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                        lx: QtyAny::new(0.1, LENGTH),
                        ly: QtyAny::new(0.1, LENGTH),
                        lz: QtyAny::new(0.01, LENGTH),
                    }),
                    materials: vec![mat_ref(cap_hash)],
                    physics: BTreeSet::from([PhysicsDomain::Mechanics]),
                    fidelity: policy_rungs(vec![FidelityRung::Lumped]),
                    uncertainty: Vec::new(),
                },
            ],
            joints: Vec::new(),
            gear_meshes: Vec::new(),
            interfaces: vec![DraftInterfaceDecl {
                name: "bleed_resistor".to_owned(),
                entity: shunt_entity,
                feature_a: cap_term,
                feature_b: gnd_term,
                kind: DraftInterfaceKind::ElectricalContact {
                    law_card: ohmic_hash,
                    third_body: None,
                    state: Vec::new(),
                },
            }],
            sources: Vec::new(),
            controllers: Vec::new(),
            probes: Vec::new(),
        },
        &scenario,
    )
    .unwrap();

    let mut lowered = registry
        .lower(&draft, &cards)
        .expect("custom EM assembly lowers cleanly");
    assert_eq!(lowered.receipt.total_property_receipts(), 2);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for _ in 0..40 {
            let _ = lowered.system.step(1.0e-4, cx).unwrap();
        }
    });
    let total_dissipated = lowered.system.ledger().total_dissipated_j();
    assert!(
        total_dissipated > 0.0,
        "RC discharge must dissipate positive Joule heat"
    );
    assert!(lowered.system.ledger().energy_residual_j().abs() <= 1.0e-12);

    // A capacitor receives current, even though its output is voltage. Keep
    // a positive physical control: 2 A for 0.1 ms adds exactly 0.2 mC.
    let mut driven = draft.input().clone();
    driven.interfaces.clear();
    driven.sources.push(SourceDecl {
        name: "charging_current".to_owned(),
        target: SourceTarget::Port {
            entity: cap_entity,
            port: "electrical".to_owned(),
        },
        history: TimeHistory::Constant(QtyAny::new(2.0, Dims([0, 0, 0, 0, 1, 0]))),
    });
    let current_draft = DraftAssembly::admit(driven.clone(), &scenario).unwrap();
    let mut charged = registry
        .lower(&current_draft, &cards)
        .expect("current source lowers");
    with_cx(&gate, |cx| {
        charged.system.step(1.0e-4, cx).unwrap();
    });
    let charge: &f64 = charged.system.state(charged.node_map["bus_cap"]).unwrap();
    assert!((*charge - 5.2e-3).abs() <= 1.0e-15);
    assert!(charged.system.ledger().energy_residual_j().abs() <= 1.0e-12);

    driven.sources[0].history = TimeHistory::Constant(QtyAny::new(2.0, Dims([2, 1, -3, 0, -1, 0])));
    let voltage_draft = DraftAssembly::admit(driven, &scenario).unwrap();
    let error = registry
        .lower(&voltage_draft, &cards)
        .err()
        .expect("voltage is not current");
    assert!(
        matches!(error, LoweringError::Invalid { ref item, ref detail }
        if item == "source:charging_current" && detail.contains("dimensions"))
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn acoustic_assembly_families_lower_through_registry_and_refuse_unknown_family() {
    let cards = CardStore::new();
    let mut registry = LoweringRegistry::with_kernel_defaults();

    // 1. Coupled string + plate + Helmholtz cavity + duct lowers through registered families.
    let spec = AcousticAssemblySpec {
        name: "guitar_and_helmholtz_prototype".to_owned(),
        sample_rate_hz: 48_000,
        ambient_temperature_k: 293.15,
        ambient_pressure_pa: 101_325.0,
        members: vec![
            AcousticMemberSpec {
                name: "spruce_top".to_owned(),
                family: "plate".to_owned(),
                material_card: None,
                parameters: BTreeMap::from([
                    (
                        "modal_mass_kg".to_owned(),
                        (0.045, Dims([0, 1, 0, 0, 0, 0])),
                    ),
                    ("modal_frequency_hz".to_owned(), (220.0, FREQUENCY)),
                    ("modal_damping_n_s_m".to_owned(), (0.8, DAMPING)),
                ]),
                coupled_to: None,
            },
            AcousticMemberSpec {
                name: "e4_string".to_owned(),
                family: "string".to_owned(),
                material_card: None,
                parameters: BTreeMap::from([
                    ("length_m".to_owned(), (0.65, LENGTH)),
                    ("tension_n".to_owned(), (72.0, FORCE)),
                    (
                        "linear_density_kg_m".to_owned(),
                        (0.0004, Dims([-1, 1, 0, 0, 0, 0])),
                    ),
                    ("modal_damping_n_s_m".to_owned(), (0.015, DAMPING)),
                ]),
                coupled_to: Some("spruce_top".to_owned()),
            },
            AcousticMemberSpec {
                name: "body_cavity".to_owned(),
                family: "cavity".to_owned(),
                material_card: None,
                parameters: BTreeMap::from([(
                    "acoustic_compliance_m5_n".to_owned(),
                    (1.2e-7, Dims([4, -1, 2, 0, 0, 0])),
                )]),
                coupled_to: None,
            },
            AcousticMemberSpec {
                name: "soundhole_duct".to_owned(),
                family: "duct".to_owned(),
                material_card: None,
                parameters: BTreeMap::from([
                    (
                        "acoustic_inertance_kg_m4".to_owned(),
                        (18.0, Dims([-4, 1, 0, 0, 0, 0])),
                    ),
                    (
                        "viscothermal_resistance".to_owned(),
                        (500.0, Dims([-4, 1, -1, 0, 0, 0])),
                    ),
                ]),
                coupled_to: Some("body_cavity".to_owned()),
            },
        ],
    };

    let mut lowered = registry
        .lower_acoustic_spec(&spec, &cards)
        .expect("acoustic families lower cleanly");
    assert_eq!(lowered.receipt.acoustic_members.len(), 4);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for _ in 0..16 {
            lowered.system.step(1.0 / 48_000.0, cx).unwrap();
        }
    });

    // 2. Unregistered acoustic family ("reed") refuses by name until registered.
    let mut with_reed = spec.clone();
    with_reed.members.push(AcousticMemberSpec {
        name: "cane_reed".to_owned(),
        family: "reed".to_owned(),
        material_card: None,
        parameters: BTreeMap::from([("reed_inertance".to_owned(), (1.0e-4, INDUCTANCE))]),
        coupled_to: None,
    });
    let err = registry
        .lower_acoustic_spec(&with_reed, &cards)
        .err()
        .expect("unregistered reed family must refuse");
    assert!(matches!(
        err,
        LoweringError::MissingAcousticFamilyConstructor { .. }
    ));
    assert!(err.to_string().contains("cane_reed"));
    assert!(err.to_string().contains("reed"));

    // Registering "reed" family makes the exact same spec lower with zero kernel edits.
    registry.register_acoustic_family("reed", |ctx, graph| {
        let mut rec = ProvenanceRecorder::default();
        let (m, dims) = ctx.member.parameters["reed_inertance"];
        rec.record_declared_field(
            "reed_inertance",
            m,
            dims,
            format!("acoustic:{}", ctx.member.name),
            "reed_inertance",
        );
        let node_name = format!("{}.reed", ctx.member.name);
        let node = graph.add_with_tangent(
            &node_name,
            None,
            Mass::new(m).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.0,
            },
        )?;
        let p = graph.port(node, "p")?;
        graph.bind_zero_boundary(p);
        let (property_receipts, parameter_sources) = rec.finish();
        Ok(LoweredAcousticMemberOutput {
            created_nodes: vec![node_name],
            exposed_port: Some(p),
            property_receipts,
            parameter_sources,
        })
    });
    let lowered_reed = registry
        .lower_acoustic_spec(&with_reed, &cards)
        .expect("reed family lowers once registered");
    assert_eq!(lowered_reed.receipt.acoustic_members.len(), 5);
}
