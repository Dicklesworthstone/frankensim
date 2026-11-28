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
    DraftGeometrySource, DraftInitialState, DraftInterfaceDecl, DraftInterfaceKind, DraftPartDecl,
    DraftPartKind, FidelityErrorBudget, FidelityPolicy, FidelityRung, GearGenerator, GearMeshDecl,
    JointDecl, JointKind, PhysicsDomain, PrimitiveShape, SourceDecl, SourceTarget, StepPolicy,
    ThirdBodyDecl, TimeHistory, TimePlan, TimeScheme,
};
use fs_scenario::assembly::{MaterialRef, ParameterKey, ReferenceState, UnknownProperty};
use fs_scenario::{
    EntityCatalog, EntityDeclaration, EntityId, Environment, FrameId, InterfacePair, Scenario,
};
use fs_system::elements::{
    Admittance, Damper, Force, HeatCapacity, Mass, Motion, Slider, SliderEnd, Slip, SourceAccount,
    Spring, Storage, Waveform,
};
use fs_system::{
    AcousticAssemblySpec, AcousticMemberSpec, AnalysisKind, BlockSolveMode, BlockSolverPolicy,
    CardStore, FidelityRungKind, FidelitySelectionPolicy, GroupLoweringContext, GroupSelection,
    InterfaceKindTag, InterfaceLoweringKey, LoweredAcousticMemberOutput, LoweredGroupOutput,
    LoweredInterfaceOutput, LoweredPartOutput, LoweringError, LoweringGraph, LoweringRegistry,
    MatDbQueryEnv, MonolithicConfig, ParameterProvenance, ParameterSource, PartLoweringKey,
    ProvenanceRecorder, RouteTemperature, SystemBuilder, geometry_generator_label,
    geometry_volume_m3, heat_sink, resolve_geometry_fingerprint,
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
        scheme: TimeScheme::Explicit,
    }
}

fn base_scenario(name: &str) -> Scenario {
    Scenario::new(name, 112, Environment::earth_lab())
}

#[test]
#[allow(clippy::too_many_lines)]
fn step_sources_keep_exact_jumps_and_integrated_force_and_heat() {
    let mut cards = CardStore::new();
    let steel_hash = cards.insert_material(steel_card());
    let adjacent = f64::from_bits(0.5_f64.to_bits() + 1);
    let gate = CancelGate::new_clock_free();

    // G0/G1: each interval uses the authored constant force and heating power.
    // A jump at zero is already active; distinct adjacent jumps stay distinct.
    for (force_switch, heat_switch, expected_times) in [
        (0.0, 0.0, vec![1.0]),
        (0.375, 0.625, vec![0.375, 0.625, 1.0]),
        (0.5, adjacent, vec![0.5, adjacent, 1.0]),
        (adjacent, 0.5, vec![0.5, adjacent, 1.0]),
    ] {
        let mut catalog = EntityCatalog::new();
        let root = catalog
            .declare(EntityDeclaration::assembly("step_sources"))
            .unwrap();
        let body = catalog
            .declare(EntityDeclaration::part(root, "body"))
            .unwrap();
        let input = DraftAssemblyInput {
            name: "step_sources".to_owned(),
            root,
            catalog,
            reference_state: Some(ReferenceState::Unloaded),
            environment: DraftEnvironment::earth_lab_air(),
            time_plan: fixed_time_plan(1.0, 1.0),
            parts: vec![DraftPartDecl {
                name: "body".to_owned(),
                entity: body,
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
                initial: None,
            }],
            joints: Vec::new(),
            gear_meshes: Vec::new(),
            interfaces: Vec::new(),
            sources: vec![
                SourceDecl {
                    name: "force".to_owned(),
                    target: SourceTarget::Port {
                        entity: body,
                        port: "mechanics".to_owned(),
                    },
                    history: TimeHistory::Step {
                        t_switch: QtyAny::new(force_switch, TIME),
                        before: QtyAny::new(2.0, FORCE),
                        after: QtyAny::new(6.0, FORCE),
                    },
                },
                SourceDecl {
                    name: "heater".to_owned(),
                    target: SourceTarget::Port {
                        entity: body,
                        port: "thermal".to_owned(),
                    },
                    history: TimeHistory::Step {
                        t_switch: QtyAny::new(heat_switch, TIME),
                        before: QtyAny::new(100.0, POWER),
                        after: QtyAny::new(300.0, POWER),
                    },
                },
            ],
            controllers: Vec::new(),
            probes: Vec::new(),
        };
        let scenario = base_scenario("step_sources");
        let draft = DraftAssembly::admit(input, &scenario).unwrap();
        let initial_temperature = draft.input().environment.base.ambient_temperature.value;
        let mut lowered = LoweringRegistry::with_kernel_defaults()
            .lower(&draft, &cards)
            .unwrap();
        let mass = lowered.node_map["body.mechanics"];
        let capacity = lowered.node_map["body.thermal"];
        let force = lowered.node_map["source.force"];
        let heater = lowered.node_map["source.heater"];

        with_cx(&gate, |cx| {
            for end in expected_times {
                let report = lowered
                    .system
                    .step(1.0 - lowered.system.time_s(), cx)
                    .unwrap();
                assert_eq!(report.time_s, end, "no ramp or extra event before a jump");
                let since_force = (end - force_switch).max(0.0);
                let since_heat = (end - heat_switch).max(0.0);
                let velocity = (2.0 * end + 4.0 * since_force) / 2.0;
                let position = (end * end + 2.0 * since_force * since_force) / 2.0;
                let heat = 100.0 * end + 200.0 * since_heat;
                let motion = lowered.system.state::<Motion>(mass).unwrap();
                assert!((motion.velocity_m_s - velocity).abs() < 2.0e-14);
                assert!((motion.position_m - position).abs() < 2.0e-14);
                let temperature = *lowered.system.state::<f64>(capacity).unwrap();
                assert!((temperature - initial_temperature - heat / 1000.0).abs() < 2.0e-13);
                let force_account = lowered.system.state::<SourceAccount>(force).unwrap();
                let heat_account = lowered.system.state::<SourceAccount>(heater).unwrap();
                assert_eq!(force_account.time_s(), end);
                assert_eq!(heat_account.time_s(), end);
                assert!((force_account.energy_in_j + velocity * velocity).abs() < 1.0e-12);
                assert!((heat_account.energy_in_j + heat).abs() < 1.0e-12);
                assert!(lowered.system.ledger().energy_residual_j().abs() < 1.0e-11);
            }
        });
        assert_eq!(lowered.system.time_s(), 1.0);
    }
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
            initial: None,
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
    assert!(
        matches!(error, LoweringError::MissingPort { ref item, ref port, .. }
        if item == "controller:velocity_feedback" && port == "command")
    );
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
                initial: None,
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
                initial: None,
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
                initial: None,
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
                initial: None,
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
                initial: None,
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
    let pos = lowered
        .system
        .state::<Motion>(lowered.node_map["carriage.mechanics"])
        .unwrap()
        .position_m;
    let sx = *lowered
        .system
        .state::<f64>(lowered.node_map["interface.carriage_stop.spring"])
        .unwrap();
    assert!(sx.abs() > 1e-3 && (sx + pos).abs() <= 1e-12, "{sx} {pos}");
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
                initial: None,
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
                initial: None,
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
                    initial: None,
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
                    initial: None,
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
                initial: None,
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
                initial: None,
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
                    initial: None,
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
                    initial: None,
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
    lowered
        .system
        .check_analysis_capabilities(AnalysisKind::Stationary)
        .expect("every default acoustic owner offers its stage residual");
    assert!(
        lowered
            .system
            .block_modes()
            .contains(&BlockSolveMode::MonolithicNewtonKrylov)
    );
    lowered
        .system
        .set_solver_policy(BlockSolverPolicy::MonolithicNewtonKrylov(
            MonolithicConfig::default(),
        ));
    let mut reference = registry.lower_acoustic_spec(&spec, &cards).unwrap();
    reference
        .system
        .set_solver_policy(BlockSolverPolicy::InterfaceNewton);

    // Excite the Helmholtz pair: q is cavity integrated flow and u is duct flow.
    // The independent midpoint equations are q' = -u and L u' = q/C - R u.
    let cavity = lowered.node_map["body_cavity.compliance"];
    let duct = lowered.node_map["soundhole_duct.inertance"];
    let reference_cavity = reference.node_map["body_cavity.compliance"];
    let reference_duct = reference.node_map["soundhole_duct.inertance"];
    let compliance = 1.2e-7;
    let inertance = 18.0;
    let resistance = 500.0;
    let initial_q = 1.2e-5;
    lowered
        .system
        .set_accepted_state(cavity, initial_q)
        .unwrap();
    reference
        .system
        .set_accepted_state(reference_cavity, initial_q)
        .unwrap();
    let initial_energy = 0.5 * initial_q * initial_q / compliance;
    let dt = 1.0 / f64::from(spec.sample_rate_hz);
    let half_dt = 0.5 * dt;
    let denominator = inertance + half_dt * resistance + half_dt * half_dt / compliance;
    let mut expected_q = initial_q;
    let mut expected_u = 0.0;
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for _ in 0..32 {
            let midpoint_u =
                (inertance * expected_u + half_dt * expected_q / compliance) / denominator;
            expected_q -= dt * midpoint_u;
            expected_u = 2.0 * midpoint_u - expected_u;

            let report = lowered.system.step(dt, cx).unwrap();
            let reference_report = reference.system.step(dt, cx).unwrap();
            assert_eq!(report.dt_s, dt);
            assert_eq!(reference_report.dt_s, dt);
            let q = *lowered.system.state::<f64>(cavity).unwrap();
            let u = lowered.system.state::<Motion>(duct).unwrap().velocity_m_s;
            let reference_q = *reference.system.state::<f64>(reference_cavity).unwrap();
            let reference_u = reference
                .system
                .state::<Motion>(reference_duct)
                .unwrap()
                .velocity_m_s;
            assert!((q - expected_q).abs() < 2.0e-13);
            assert!((u - expected_u).abs() < 2.0e-11);
            assert!((q - reference_q).abs() < 2.0e-13);
            assert!((u - reference_u).abs() < 2.0e-11);
            let stored_energy = 0.5 * q * q / compliance + 0.5 * inertance * u * u;
            let dissipated_energy = lowered.system.ledger().total_dissipated_j();
            assert!((stored_energy + dissipated_energy - initial_energy).abs() < 2.0e-12);
        }
    });
    assert!(lowered.system.ledger().total_dissipated_j() > 0.0);

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

#[test]
fn typed_t_axis_matches_absolute_temperature_keeps_legacy_and_refuses_mismatch() {
    use fs_matdb::{MatDbError::AxisQuantityMismatch, SelectionPolicy::SingleClaimOnly};
    use fs_qty::semantic::{QuantityKind as Qk, SemanticType, ValueForm::Static};
    let env = fs_system::MatDbQueryEnv {
        ambient_temp_k: 293.15,
        ambient_pressure_pa: 101_325.0,
        policy: SingleClaimOnly,
    };
    let run = |validity: ValidityDomain, explicit_t: bool| {
        let (mut claims, mut cards) = (ClaimSet::new(), CardStore::new());
        let c = PropertyClaim {
            validity,
            ..scalar_claim("density", 7800.0, DENSITY, "typed-t")
        };
        claims.insert_claim(c).unwrap();
        let card = MaterialCard::assemble(steel_card().id().clone(), claims, vec![]).unwrap();
        let mut mr = mat_ref(cards.insert_material(card));
        if explicit_t {
            let k = ParameterKey::new("T", TEMPERATURE);
            mr.query_point.push((k, QtyAny::new(293.15, TEMPERATURE)));
        }
        let mut r = ProvenanceRecorder::default();
        r.query_material_property("part:m", &cards, &mr, "density", "density", env)
    };
    let tv = |k| {
        let q = fs_qty::QuantitySpec::semantic(SemanticType::new(k, Static));
        ValidityDomain::unconstrained().with_quantity("T", q, 293.15, 293.15)
    };
    assert!(run(tv(Qk::AbsoluteTemperature), false).is_ok());
    assert!(run(tv(Qk::AbsoluteTemperature), true).is_ok());
    let leg = ValidityDomain::unconstrained().with("T", 293.15, 293.15);
    assert!(run(leg, false).is_ok());
    assert!(matches!(
        run(tv(Qk::TemperatureDifference), true),
        Err(LoweringError::MatDb { source, .. }) if matches!(*source, AxisQuantityMismatch { .. })
    ));
}

#[test]
fn a_declared_key_takes_its_claims_axis_spec_or_is_refused_by_name() {
    let env = fs_system::MatDbQueryEnv {
        ambient_temp_k: 293.15,
        ambient_pressure_pa: 101_325.0,
        policy: fs_matdb::SelectionPolicy::SingleClaimOnly,
    };
    let one = fs_qty::QuantitySpec::dimensional(Dims::NONE);
    // `phase-liquid` bounded at 1, untyped or typed; one density claim per domain.
    let legacy = ValidityDomain::unconstrained().with("phase-liquid", 1.0, 1.0);
    let typed = legacy.clone().with_quantity("phase-liquid", one, 1.0, 1.0);
    let run = |domains: &[&ValidityDomain], dims: Dims| {
        let (mut claims, mut cards) = (ClaimSet::new(), CardStore::new());
        for (i, &validity) in domains.iter().enumerate() {
            let mut c = scalar_claim("density", 7800.0, DENSITY, &format!("claim {i}"));
            c.validity = validity.clone();
            claims.insert_claim(c).unwrap();
        }
        let card = MaterialCard::assemble(steel_card().id().clone(), claims, vec![]).unwrap();
        let mut mr = mat_ref(cards.insert_material(card));
        let key = ParameterKey::new("phase-liquid", dims);
        mr.query_point.push((key, QtyAny::new(1.0, dims)));
        let mut r = ProvenanceRecorder::default();
        let answer = r.query_material_property("part:m", &cards, &mr, "density", "density", env);
        answer.map(|_| ()).map_err(|e| e.to_string())
    };
    // Typed by its claim, where it went in untyped and was refused before.
    assert_eq!(run(&[&typed], Dims::NONE), Ok(()));
    assert_eq!(run(&[&legacy], Dims::NONE), Ok(()));
    let named = |why: &str| format!("query key `phase-liquid`: its {why}");
    let refusal = run(&[&typed], TEMPERATURE).unwrap_err();
    assert!(refusal.contains(&named("dims are not")), "{refusal}");
    let refusal = run(&[&typed, &legacy], Dims::NONE).unwrap_err();
    assert!(refusal.contains(&named("claims disagree")), "{refusal}");
}

/// A toy owner group: one heat capacity of the card's density and specific
/// heat over the part's box volume, each with its receipt.
fn toy_capacity(
    ctx: &GroupLoweringContext<'_>,
    graph: &mut LoweringGraph,
) -> Result<LoweredGroupOutput, LoweringError> {
    let input = ctx.draft.input();
    let part = &input.parts[ctx.selection.parts[0]];
    let scope = format!("part:{}", part.name);
    let env = MatDbQueryEnv {
        ambient_temp_k: input.environment.base.ambient_temperature.value,
        ambient_pressure_pa: input.environment.base.ambient_pressure.value,
        policy: ctx.selection_policy,
    };
    let mut rec = ProvenanceRecorder::default();
    let mut query = |property| {
        let material = &part.materials[0];
        let (sample, _) =
            rec.query_material_property(&scope, ctx.cards, material, property, property, env)?;
        Ok::<f64, LoweringError>(sample.value)
    };
    let density = query("density")?;
    let specific_heat = query("specific_heat_capacity")?;
    let (volume, formula) = geometry_volume_m3(&part.geometry).expect("a box");
    rec.record_geometry_derived(
        "volume",
        volume,
        Dims([3, 0, 0, 0, 0, 0]),
        resolve_geometry_fingerprint(&part.geometry),
        geometry_generator_label(&part.geometry),
        formula,
    );
    let capacity = HeatCapacity::new(density * specific_heat * volume).expect("a capacity");
    let node = graph.add("toy.block", Some(part.entity), capacity, env.ambient_temp_k)?;
    // Nothing in the draft reaches the block: it is insulated.
    let heat = graph.port(node, "heat")?;
    graph.bind_zero_boundary(heat);
    let (property_receipts, parameter_sources) = rec.finish();
    Ok(LoweredGroupOutput {
        parts: vec![LoweredPartOutput {
            created_nodes: vec!["toy.block".to_owned()],
            property_receipts,
            parameter_sources,
        }],
        ..LoweredGroupOutput::default()
    })
}

#[test]
fn an_owner_group_lowers_its_part_whole_and_two_groups_cannot_share_one() {
    let mut cards = CardStore::new();
    let steel_hash = cards.insert_material(steel_card());
    let mut catalog = EntityCatalog::new();
    let root = catalog
        .declare(EntityDeclaration::assembly("toy_group"))
        .unwrap();
    let block = catalog
        .declare(EntityDeclaration::part(root, "block"))
        .unwrap();
    let input = DraftAssemblyInput {
        name: "toy_group".to_owned(),
        root,
        catalog,
        reference_state: Some(ReferenceState::Unloaded),
        environment: DraftEnvironment::earth_lab_air(),
        time_plan: fixed_time_plan(1.0, 1.0),
        parts: vec![DraftPartDecl {
            name: "block".to_owned(),
            entity: block,
            kind: DraftPartKind::Solid,
            geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                lx: QtyAny::new(0.1, LENGTH),
                ly: QtyAny::new(0.05, LENGTH),
                lz: QtyAny::new(0.05, LENGTH),
            }),
            materials: vec![mat_ref(steel_hash)],
            physics: BTreeSet::from([PhysicsDomain::Thermal]),
            fidelity: policy_rungs(vec![FidelityRung::Lumped]),
            uncertainty: Vec::new(),
            initial: None,
        }],
        joints: Vec::new(),
        gear_meshes: Vec::new(),
        interfaces: Vec::new(),
        sources: Vec::new(),
        controllers: Vec::new(),
        probes: Vec::new(),
    };
    let draft = DraftAssembly::admit(input, &base_scenario("toy_group")).unwrap();
    let first_part = |_: &DraftAssembly, _: &CardStore| {
        Ok(vec![GroupSelection {
            parts: vec![0],
            ..GroupSelection::default()
        }])
    };
    let mut registry = LoweringRegistry::with_kernel_defaults();
    registry.register_group("toy", first_part, toy_capacity);
    // Registering the name again replaces the group: two groups named `toy`
    // would both take the block, and the lowering would refuse.
    registry.register_group("toy", first_part, toy_capacity);

    let lowered = registry.lower(&draft, &cards).unwrap();
    // The group's node, and not the kernel default's (named after the part).
    assert!(lowered.node_map.contains_key("toy.block"));
    assert!(!lowered.node_map.contains_key("block"));
    let receipt = &lowered.receipt.components[0];
    assert_eq!(
        receipt.fidelity_rationale,
        "lowered whole by owner group 'toy'"
    );
    assert_eq!(receipt.created_nodes, ["toy.block"]);
    assert_eq!(receipt.property_receipts.len(), 2);
    assert_eq!(receipt.parameter_sources.len(), 3);

    // A second group takes the same part (and creates nothing).
    registry.register_group("other", first_part, |_, _| {
        Ok(LoweredGroupOutput {
            parts: vec![LoweredPartOutput::default()],
            ..LoweredGroupOutput::default()
        })
    });
    let Err(LoweringError::Invalid { item, detail }) = registry.lower(&draft, &cards) else {
        panic!("two groups took one part");
    };
    assert_eq!(item, "part:block");
    assert_eq!(detail, "selected by owner groups 'toy' and 'other'");
}

#[test]
fn one_d_field_thermal_conductance_lowers_observes_and_refuses_invalid() {
    let (mut cards, mut catalog) = (CardStore::new(), EntityCatalog::new());
    let ch = cards.insert_material(tim_paste_card());
    let root = catalog.declare(EntityDeclaration::assembly("s")).unwrap();
    let slab = catalog
        .declare(EntityDeclaration::part(root, "slab"))
        .unwrap();
    let sources = [("h", "a", 307.15), ("c", "b", 287.15)]
        .map(|(n, p, t)| {
            port_src(
                n,
                slab,
                p,
                TimeHistory::Constant(QtyAny::new(t, TEMPERATURE)),
            )
        })
        .into();
    let slab_part = DraftPartDecl {
        name: "slab".into(),
        entity: slab,
        kind: DraftPartKind::Solid,
        geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
            lx: QtyAny::new(0.02, LENGTH),
            ly: QtyAny::new(0.5, LENGTH),
            lz: QtyAny::new(0.5, LENGTH),
        }),
        materials: vec![mat_ref(ch)],
        physics: BTreeSet::from([PhysicsDomain::Thermal]),
        fidelity: policy_rungs(vec![FidelityRung::OneDField { cells: 1 }]),
        uncertainty: Vec::new(),
        initial: None,
    };
    let mut input = empty_draft("s", root, catalog, vec![slab_part], sources);
    input.time_plan = fixed_time_plan(1.0, 0.1);
    let (reg, scen) = (LoweringRegistry::with_kernel_defaults(), base_scenario("s"));
    let mut lo = reg
        .lower(&DraftAssembly::admit(input.clone(), &scen).unwrap(), &cards)
        .unwrap();
    with_cx(&CancelGate::new(), |cx| lo.system.step(0.1, cx).unwrap());
    assert!(matches!(
        lo.observers[&slab](&lo.system),
        Some(fs_system::lowering::Observed::HeatRate(q)) if (q - 1000.0).abs() <= 1.0e-12
    ));
    let bad = |i| {
        reg.lower(&DraftAssembly::admit(i, &scen).unwrap(), &cards)
            .err()
            .unwrap()
            .to_string()
    };
    input.parts[0].fidelity = policy_rungs(vec![FidelityRung::OneDField { cells: 2 }]);
    assert!(bad(input.clone()).contains("cells: 2"));
    input.parts[0].fidelity = policy_rungs(vec![FidelityRung::OneDField { cells: 1 }]);
    input.parts[0].geometry = DraftGeometrySource::Primitive(PrimitiveShape::Cylinder {
        radius: QtyAny::new(0.1, LENGTH),
        height: QtyAny::new(0.02, LENGTH),
    });
    assert!(bad(input).contains("Cylinder"));
}

fn box_part(
    name: &str,
    entity: EntityId,
    kind: DraftPartKind,
    mat: ContentHash,
    phys: &[PhysicsDomain],
) -> DraftPartDecl {
    let q = QtyAny::new(0.1, LENGTH);
    DraftPartDecl {
        name: name.into(),
        entity,
        kind,
        geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
            lx: q,
            ly: q,
            lz: q,
        }),
        materials: vec![mat_ref(mat)],
        physics: phys.iter().copied().collect(),
        fidelity: policy_rungs(vec![FidelityRung::Lumped]),
        uncertainty: vec![],
        initial: None,
    }
}

fn empty_draft(
    name: &str,
    root: EntityId,
    catalog: EntityCatalog,
    parts: Vec<DraftPartDecl>,
    sources: Vec<SourceDecl>,
) -> DraftAssemblyInput {
    DraftAssemblyInput {
        name: name.into(),
        root,
        catalog,
        reference_state: Some(ReferenceState::Assembled),
        environment: DraftEnvironment::earth_lab_air(),
        time_plan: fixed_time_plan(0.2, 0.01),
        parts,
        joints: vec![],
        gear_meshes: vec![],
        interfaces: vec![],
        sources,
        controllers: vec![],
        probes: vec![],
    }
}

fn port_src(n: &str, e: EntityId, p: &str, history: TimeHistory) -> SourceDecl {
    SourceDecl {
        name: n.into(),
        target: SourceTarget::Port {
            entity: e,
            port: p.into(),
        },
        history,
    }
}

#[test]
fn lumped_initial_state_propagates_and_refuses_off_axis() {
    use PhysicsDomain as P;
    let (vel, q) = (Dims([1, 0, -1, 0, 0, 0]), QtyAny::new);
    let (mut cards, mut cat) = (CardStore::new(), EntityCatalog::new());
    let steel = cards.insert_material(steel_card());
    let root = cat.declare(EntityDeclaration::assembly("init")).unwrap();
    let body = cat.declare(EntityDeclaration::part(root, "body")).unwrap();
    let solid = cat.declare(EntityDeclaration::part(root, "solid")).unwrap();
    let mut ib = DraftInitialState {
        position: [q(0.25, LENGTH), q(0.0, LENGTH), q(0.0, LENGTH)],
        orientation: [1.0, 0.0, 0.0, 0.0],
        linear_velocity: [q(1.5, vel), q(0.0, vel), q(0.0, vel)],
        angular_velocity: [q(0.0, FREQUENCY); 3],
        temperature: Some(q(340.0, TEMPERATURE)),
    };
    let is = DraftInitialState {
        position: [q(0.0, LENGTH); 3],
        linear_velocity: [q(0.0, vel); 3],
        temperature: Some(q(365.0, TEMPERATURE)),
        ..ib.clone()
    };
    let mut bp = box_part(
        "body",
        body,
        DraftPartKind::RigidBody,
        steel,
        &[P::Mechanics, P::Thermal],
    );
    bp.initial = Some(ib.clone());
    let mut sp = box_part("solid", solid, DraftPartKind::Solid, steel, &[P::Thermal]);
    sp.initial = Some(is);
    let c0 = |d| TimeHistory::Constant(q(0.0, d));
    let srcs = vec![
        port_src("f0", body, "mechanics", c0(FORCE)),
        port_src("q0", body, "thermal", c0(POWER)),
        port_src("q1", solid, "thermal", c0(POWER)),
    ];
    let mut inp = empty_draft("init", root, cat, vec![bp, sp], srcs);
    let reg = LoweringRegistry::with_kernel_defaults();
    let scen = base_scenario("init");
    let go = |i| reg.lower(&DraftAssembly::admit(i, &scen).unwrap(), &cards);
    let tmp =
        |l: &fs_system::LoweredSystem, k: &str| *l.system.state::<f64>(l.node_map[k]).unwrap();
    let mut lo = go(inp.clone()).unwrap();
    let m0 = *lo
        .system
        .state::<Motion>(lo.node_map["body.mechanics"])
        .unwrap();
    assert!((m0.position_m - 0.25).abs() <= 1e-15 && (m0.velocity_m_s - 1.5).abs() <= 1e-15);
    assert!((tmp(&lo, "body.thermal") - 340.0).abs() <= 1e-15);
    assert!((tmp(&lo, "solid") - 365.0).abs() <= 1e-15);
    with_cx(&CancelGate::new(), |cx| lo.system.step(0.1, cx).unwrap());
    let m1 = *lo
        .system
        .state::<Motion>(lo.node_map["body.mechanics"])
        .unwrap();
    assert!((m1.position_m - 0.40).abs() <= 1e-14 && (m1.velocity_m_s - 1.5).abs() <= 1e-14);
    let mut check_bad = |ax: &str, val: &str, edit: fn(&mut DraftInitialState)| {
        let mut bi = ib.clone();
        edit(&mut bi);
        inp.parts[0].initial = Some(bi);
        let msg = go(inp.clone()).err().unwrap().to_string();
        assert!(msg.contains(ax) && msg.contains(val), "{msg}");
    };
    check_bad("initial.linear_velocity[1]", "0.75", |i| {
        i.linear_velocity[1].value = 0.75;
    });
    check_bad("initial.position[2]", "0.1", |i| i.position[2].value = 0.1);
    check_bad("initial.angular_velocity[0]", "2", |i| {
        i.angular_velocity[0].value = 2.0;
    });
    check_bad("initial.orientation", "0.707", |i| {
        let r = core::f64::consts::FRAC_1_SQRT_2;
        i.orientation = [r, r, 0.0, 0.0];
    });
    ib.temperature = None;
    inp.parts[0].initial = Some(ib);
    assert!((tmp(&go(inp).unwrap(), "body.thermal") - 293.15).abs() <= 1e-12);
}

/// The PID position servo of the two tests below: kp [N/m], ki [N/(m s)], kd [N s/m] and the
/// setpoint [m], sampled at 16 Hz, holding a free body against a constant push [N].
const SERVO: [f64; 4] = [72.0, 36.0, 38.4, 0.05];
const SERVO_PUSH_N: f64 = -12.0;

/// A free 8 kg steel box (a 0.1 m cube at 8000 kg/m^3) under the servo [`SERVO`], with a
/// thermal solid `warm` beside it; returns the warm solid's entity too.
fn servo_draft() -> (CardStore, DraftAssemblyInput, EntityId) {
    let ([kp, ki, kd, r], q) = (SERVO, QtyAny::new);
    let (mut cards, mut cat) = (CardStore::new(), EntityCatalog::new());
    let steel = cards.insert_material(steel_card());
    let root = cat.declare(EntityDeclaration::assembly("servo")).unwrap();
    let body = cat.declare(EntityDeclaration::part(root, "body")).unwrap();
    let warm = cat.declare(EntityDeclaration::part(root, "warm")).unwrap();
    let (mech, heat) = ([PhysicsDomain::Mechanics], [PhysicsDomain::Thermal]);
    let parts = vec![
        box_part("body", body, DraftPartKind::RigidBody, steel, &mech),
        box_part("warm", warm, DraftPartKind::Solid, steel, &heat),
    ];
    let (push, still) = (q(SERVO_PUSH_N, FORCE), q(0.0, POWER));
    let srcs = vec![
        port_src("push", body, "mechanics", TimeHistory::Constant(push)),
        port_src("q0", warm, "thermal", TimeHistory::Constant(still)),
    ];
    let mut inp = empty_draft("servo", root, cat, parts, srcs);
    inp.controllers.push(ControllerDecl {
        name: "servo".to_owned(),
        sensor_ports: vec![(body, "position".to_owned())],
        actuator_port: (body, "mechanics".to_owned()),
        sample_rate: q(16.0, FREQUENCY),
        law: ControllerLaw::Pid {
            kp: q(kp, STIFFNESS),
            ki: q(ki, Dims([0, 1, -3, 0, 0, 0])),
            kd: q(kd, DAMPING),
            setpoint: q(r, LENGTH),
        },
    });
    (cards, inp, warm)
}

#[test]
fn a_declared_pid_servo_holds_a_pushed_body_at_its_setpoint() {
    let ([kp, ki, kd, r], (m, t, dt)) = (SERVO, (8.0, 1.0 / 16.0, 1.0 / 64.0));
    let (cards, inp, _) = servo_draft();
    let draft = DraftAssembly::admit(inp, &base_scenario("servo")).unwrap();
    let mut lo = LoweringRegistry::with_kernel_defaults()
        .lower(&draft, &cards)
        .unwrap();
    let rows = &lo.receipt.sources;
    let row = rows.iter().find(|s| s.source_name == "controller:servo");
    let nodes = row.map(|s| s.created_nodes.clone()).unwrap_or_default();
    assert_eq!(nodes, ["controller.servo", "controller.servo.drive"]);

    // The same discrete loop written out: at each sample the error, its sum and its last
    // value give the command, held over four midpoint steps of the free body.
    let (gate, node) = (CancelGate::new(), lo.node_map["body"]);
    let (mut x, mut v, mut sum, mut e_prev, mut worst) = (0.0_f64, 0.0_f64, 0.0, 0.0, 0.0_f64);
    for _ in 0..192 {
        let e = r - x;
        sum += e;
        let u = kp * e + ki * t * sum + kd * (e - e_prev) / t;
        e_prev = e;
        for _ in 0..4 {
            let v1 = v + (u + SERVO_PUSH_N) * dt / m;
            x += 0.5 * (v + v1) * dt;
            v = v1;
            with_cx(&gate, |cx| lo.system.step(dt, cx).unwrap());
            let got = *lo.system.state::<Motion>(node).unwrap();
            worst = worst.max((got.position_m - x).abs().max((got.velocity_m_s - v).abs()));
        }
    }
    assert!(worst <= 1e-12, "{worst:e}");
    // 12 s on, the integral term holds the push: the body sits at its setpoint.
    assert!((x - r).abs() <= 1e-4 && v.abs() <= 1e-4, "{x} {v}");
}

#[test]
fn a_controller_that_cannot_lower_refuses_by_name() {
    type Edit = fn(&mut DraftAssemblyInput);
    let (cards, inp, warm) = servo_draft();
    let (reg, scen) = (
        LoweringRegistry::with_kernel_defaults(),
        base_scenario("servo"),
    );
    let cases: [(&str, Edit); 5] = [
        ("lookup-table", |i| {
            let (a, b) = (QtyAny::new(0.0, LENGTH), QtyAny::new(1.0, LENGTH));
            let (fa, fb) = (QtyAny::new(0.0, FORCE), QtyAny::new(1.0, FORCE));
            let breakpoints = vec![(a, fa), (b, fb)];
            i.controllers[0].law = ControllerLaw::LookupTable { breakpoints };
        }),
        ("card-backed", |i| {
            let (law_card, params) = (hash_bytes(b"servo law"), vec![]);
            i.controllers[0].law = ControllerLaw::CardBacked { law_card, params };
        }),
        ("one sensor, not 2", |i| {
            let s = i.controllers[0].sensor_ports[0].clone();
            i.controllers[0].sensor_ports.push(s);
        }),
        ("sensor 'acceleration'", |i| {
            i.controllers[0].sensor_ports[0].1 = "acceleration".to_owned();
        }),
        ("kp has dimensions [0, 1, -1, 0, 0, 0]", |i| {
            if let ControllerLaw::Pid { kp, .. } = &mut i.controllers[0].law {
                kp.dims = DAMPING;
            }
        }),
    ];
    let mut heat = inp.clone();
    heat.controllers[0].actuator_port = (warm, "thermal".to_owned());
    let mut refusals = vec![("need [1, 1, -3, 0, 0, 0]", heat)];
    for (want, edit) in cases {
        let mut i = inp.clone();
        edit(&mut i);
        refusals.push((want, i));
    }
    for (want, i) in refusals {
        let draft = DraftAssembly::admit(i, &scen).unwrap();
        let error = reg.lower(&draft, &cards).err();
        let named = matches!(&error, Some(LoweringError::Invalid { item, detail })
            if item == "controller:servo" && detail.contains(want));
        assert!(named, "{want}: {error:?}");
    }
}

/// A PI thermostat [W/K, W/(K s)] on the thermal solid `warm` of [`servo_draft`], which
/// loses no heat, raising it 5 K from ambient through a heater behind a limiter at 0.
/// The pushed body beside it shares no port with `warm`.
const THERMOSTAT: [f64; 3] = [8000.0, 16000.0, 5.0];

#[test]
fn a_declared_thermostat_heats_through_a_limiter_and_never_cools() {
    let ([kp, ki, rise], (t, dt), q) = (THERMOSTAT, (1.0 / 16.0, 1.0 / 64.0), QtyAny::new);
    let (cards, mut inp, warm) = servo_draft();
    let t0 = inp.environment.base.ambient_temperature.value;
    let r = t0 + rise;
    inp.controllers[0] = ControllerDecl {
        name: "thermostat".to_owned(),
        sensor_ports: vec![(warm, "temperature".to_owned())],
        actuator_port: (warm, "thermal".to_owned()),
        sample_rate: q(16.0, FREQUENCY),
        law: ControllerLaw::Pid {
            kp: q(kp, Dims([2, 1, -3, -1, 0, 0])),
            ki: q(ki, Dims([2, 1, -4, -1, 0, 0])),
            kd: q(0.0, Dims([2, 1, -2, -1, 0, 0])),
            setpoint: q(r, TEMPERATURE),
        },
    };
    let draft = DraftAssembly::admit(inp, &base_scenario("servo")).unwrap();
    let mut lo = LoweringRegistry::with_kernel_defaults()
        .lower(&draft, &cards)
        .unwrap();
    let row = lo
        .receipt
        .sources
        .iter()
        .find(|s| s.source_name == "controller:thermostat");
    let names = [
        "controller.thermostat",
        "controller.thermostat.limit",
        "controller.thermostat.drive",
    ];
    assert_eq!(row.unwrap().created_nodes, names);
    let part = lo
        .receipt
        .components
        .iter()
        .find(|c| c.part_name == "warm")
        .unwrap();
    let c = part
        .parameter_sources
        .iter()
        .find(|p| p.parameter == "heat_capacity")
        .unwrap()
        .value;

    // The loop written out: the heater takes max(u, 0) and the solid stores all of it.
    let (gate, node) = (CancelGate::new(), lo.node_map["warm"]);
    let (mut temp, mut sum, mut worst, mut work, mut cut) = (t0, 0.0, 0.0_f64, 0.0, 0);
    let mut last = t0;
    for _ in 0..96 {
        let e = r - temp;
        sum += e;
        let u: f64 = kp * e + ki * t * sum;
        cut += usize::from(u < 0.0);
        for _ in 0..4 {
            temp += u.max(0.0) * dt / c;
            work += u.max(0.0) * dt;
            with_cx(&gate, |cx| lo.system.step(dt, cx).unwrap());
            let got = *lo.system.state::<f64>(node).unwrap();
            assert!(got >= last, "the solid cooled: {last} K -> {got} K");
            (worst, last) = (worst.max((got - temp).abs()), got);
        }
    }
    assert!(worst <= 1e-10, "{worst:e} K");
    // It overshoots, the law then asks for cooling, and the heater rests at 0 W.
    assert!(cut > 0 && last > r + 0.1, "{cut} {last} {r}");
    assert!(
        (c * (last - t0) - work).abs() <= 1e-9 * work,
        "{} J {work} J",
        c * (last - t0)
    );
}

#[test]
fn moving_driven_body_and_load_line_graph_helpers() {
    let (vel, q) = (Dims([1, 0, -1, 0, 0, 0]), QtyAny::new);
    let (mut cards, mut cat) = (CardStore::new(), EntityCatalog::new());
    let steel = cards.insert_material(steel_card());
    let mut fc = ClaimSet::new();
    fc.insert_claim(scalar_claim("k", 1.0, Dims::NONE, "fixture"))
        .unwrap();
    let fcard = cards
        .insert_interface("toy-slider", interface_card(fc, "dry"))
        .unwrap();
    let root = cat.declare(EntityDeclaration::assembly("belt")).unwrap();
    let blk = cat.declare(EntityDeclaration::part(root, "block")).unwrap();
    let blt = cat.declare(EntityDeclaration::part(root, "belt")).unwrap();
    let sa = cat.declare(EntityDeclaration::surface(blk, "sa")).unwrap();
    let sb = cat.declare(EntityDeclaration::surface(blt, "sb")).unwrap();
    let pair = InterfacePair::ordered(sa, sb);
    let cid = cat
        .declare(EntityDeclaration::interface(root, "c", pair))
        .unwrap();
    let knots = vec![(q(0.0, TIME), q(0.0, vel)), (q(0.1, TIME), q(0.2, vel))];
    let src = port_src(
        "bv",
        blt,
        "mechanics",
        TimeHistory::PiecewiseLinear { knots },
    );
    let m = &[PhysicsDomain::Mechanics][..];
    let parts = vec![
        box_part("block", blk, DraftPartKind::RigidBody, steel, m),
        box_part("belt", blt, DraftPartKind::DrivenBody, steel, m),
    ];
    let mut inp = empty_draft("belt", root, cat, parts, vec![src]);
    inp.interfaces.push(DraftInterfaceDecl {
        name: "c".into(),
        entity: cid,
        feature_a: sa,
        feature_b: sb,
        kind: DraftInterfaceKind::MechanicalContact {
            law_card: fcard,
            third_body: None,
            state: vec![],
        },
    });
    let mut reg = LoweringRegistry::with_kernel_defaults();
    let key = InterfaceLoweringKey::new(InterfaceKindTag::MechanicalContact, "toy-slider");
    reg.register_interface(key, |ctx, g| {
        let sl = Slider::new(100.0, 0.5, 0.3).unwrap().told(SliderEnd::B);
        let n = g.add_on_line("s", Some(ctx.interface.entity), sl, Slip::stuck())?;
        let c = &ctx.draft.input().catalog;
        let (fa, fb) = (ctx.interface.feature_a, ctx.interface.feature_b);
        let pa = g.resolve_entity_port(c, fa, "mechanics").unwrap();
        let pb = g.resolve_entity_port(c, fb, "mechanics").unwrap();
        let sp = g.resolve_entity_port(c, fb, "span").unwrap();
        g.attach_common_flow(pa, g.port(n, "a")?, 1.0);
        g.attach_common_flow(pb, g.port(n, "b")?, 1.0);
        g.connect_signal(sp, g.port(n, "span_b")?)?;
        let hp = g.port(n, "heat")?;
        let amb = Some(heat_sink("amb", 300.0).unwrap());
        g.add_heat_route(hp, vec![], amb, RouteTemperature::Fixed(300.0));
        let mut r = ProvenanceRecorder::default();
        r.record_declared_field("normal_force", 10.0, FORCE, "toy", "n");
        let (property_receipts, parameter_sources) = r.finish();
        Ok(LoweredInterfaceOutput {
            created_nodes: vec!["s".into()],
            property_receipts,
            parameter_sources,
        })
    });
    let scen = base_scenario("belt");
    let go = |i| reg.lower(&DraftAssembly::admit(i, &scen).unwrap(), &cards);
    let err = |i| go(i).err().unwrap().to_string();
    let mut lo = go(inp.clone()).unwrap();
    assert!(lo.receipt.sources.len() == 1 && lo.part_ports.contains_key(&(blt, "span".into())));
    with_cx(&CancelGate::new(), |cx| lo.system.step(0.01, cx).unwrap());
    let bm = lo.system.state::<Motion>(lo.node_map["block"]).unwrap();
    assert!((bm.velocity_m_s - 0.02).abs() <= 1e-12, "stuck: {bm:?}");
    let mut bad = inp.clone();
    bad.sources[0].history = TimeHistory::Constant(q(5.0, FORCE));
    assert!(err(bad).contains("velocity source"));
    let mut dup = inp.sources[0].clone();
    dup.name = "bv2".into();
    inp.sources.push(dup);
    assert!(err(inp).contains("at most one"));
}

#[test]
fn part_surface_is_not_ground_and_non_thermal_contact_refuses() {
    use {DraftPartKind as K, EntityDeclaration as E, PhysicsDomain as P};
    let (mut cards, mut cat) = (CardStore::new(), EntityCatalog::new());
    let steel = cards.insert_material(steel_card());
    let mut tc = ClaimSet::new();
    let cl = scalar_claim("thermal_conductance", 5.0, CONDUCTANCE, "tc");
    tc.insert_claim(cl).unwrap();
    let card = interface_card(tc, "dry");
    let tcard = cards.insert_interface("thermal-conductance", card).unwrap();
    let root = cat.declare(E::assembly("g")).unwrap();
    let hot = cat.declare(E::part(root, "hot")).unwrap();
    let cold = cat.declare(E::part(root, "cold")).unwrap();
    let sh = cat.declare(E::surface(hot, "sh")).unwrap();
    let sc = cat.declare(E::surface(cold, "sc")).unwrap();
    let pair = InterfacePair::ordered(sh, sc);
    let cid = cat.declare(E::interface(root, "lk", pair)).unwrap();
    let parts = vec![
        box_part("hot", hot, K::Solid, steel, &[P::Thermal]),
        box_part("cold", cold, K::RigidBody, steel, &[P::Mechanics]),
    ];
    let mut inp = empty_draft("g", root, cat, parts, vec![]);
    let kind = DraftInterfaceKind::ThermalContact {
        law_card: tcard,
        third_body: None,
        state: vec![],
    };
    inp.interfaces.push(DraftInterfaceDecl {
        name: "lk".into(),
        entity: cid,
        feature_a: sh,
        feature_b: sc,
        kind,
    });
    let draft = DraftAssembly::admit(inp, &base_scenario("g")).unwrap();
    let err = LoweringRegistry::with_kernel_defaults()
        .lower(&draft, &cards)
        .err()
        .expect("must not ground");
    assert!(matches!(err, LoweringError::MissingPort { ref port, .. } if port == "thermal"));
}
