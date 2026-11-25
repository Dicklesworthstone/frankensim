//! Bead `frankensim-rc-1109-root-3gc12.11.1`: Canonical assembly draft validation.
//!
//! Proves that all eight target devices — the five kernel slices (A–E), the
//! heatsink project, a clarinet note, and a mechanical watch — are expressible
//! cleanly with zero escape hatches, round-trip identically in memory, and
//! refuse every invalid field by name.

#![allow(clippy::too_many_lines)]

use std::collections::BTreeSet;

use fs_blake3::{ContentHash, hash_bytes};
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::draft::{
    ControllerDecl, ControllerLaw, DraftAssembly, DraftAssemblyInput, DraftEnvironment,
    DraftGeometrySource, DraftInterfaceDecl, DraftInterfaceKind, DraftPartDecl, DraftPartKind,
    FidelityErrorBudget, FidelityPolicy, FidelityRung, GasComponent, GearGenerator, GearMeshDecl,
    JointDecl, JointKind, JointLimits, MicrophonePattern, PhysicsDomain, PrimitiveShape, ProbeDecl,
    ProbeKind, SourceDecl, SourceTarget, SpiralSpringGenerator, StepPolicy, ThirdBodyDecl,
    TimeHistory, TimePlan, TimingEdge,
};
use fs_scenario::assembly::{GeometrySource, MaterialRef, ParameterKey, ReferenceState};
use fs_scenario::{
    EntityCatalog, EntityDeclaration, EntityId, Environment, FrameId, InterfacePair, Scenario,
};

const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const FREQUENCY: Dims = Dims([0, 0, -1, 0, 0, 0]);
const FORCE: Dims = Dims([1, 1, -2, 0, 0, 0]);
const POWER: Dims = Dims([2, 1, -3, 0, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const TEMPERATURE: Dims = Dims([0, 0, 0, 1, 0, 0]);
const MAGNETIC_FIELD: Dims = Dims([0, 1, -2, 0, -1, 0]);

fn h(label: &str) -> ContentHash {
    hash_bytes(label.as_bytes())
}

fn mat(card_name: &str) -> MaterialRef {
    MaterialRef {
        card: h(card_name),
        query_point: vec![(
            ParameterKey::new("temperature", TEMPERATURE),
            QtyAny::new(293.15, TEMPERATURE),
        )],
        frame: FrameId(0),
        layer: None,
        unknown: Vec::new(),
    }
}

fn lumped_policy() -> FidelityPolicy {
    FidelityPolicy {
        rungs: vec![FidelityRung::Lumped],
        error_budget: FidelityErrorBudget {
            relative_tol: 1.0e-6,
            absolute_tol: QtyAny::dimensionless(1.0e-12),
            max_work_units: Some(1_000_000),
        },
    }
}

fn time_plan(horizon_s: f64, dt_s: f64, sample_hz: f64) -> TimePlan {
    TimePlan {
        horizon: QtyAny::new(horizon_s, TIME),
        step_policy: StepPolicy::Fixed {
            dt: QtyAny::new(dt_s, TIME),
        },
        default_sample_rate: QtyAny::new(sample_hz, FREQUENCY),
    }
}

fn add_part_with_surface(
    catalog: &mut EntityCatalog,
    root: EntityId,
    name: &str,
) -> (EntityId, EntityId) {
    let part = catalog
        .declare(EntityDeclaration::part(root, name))
        .unwrap();
    let surf = catalog
        .declare(EntityDeclaration::surface(part, &format!("{name}_surface")))
        .unwrap();
    (part, surf)
}

fn add_iface(
    catalog: &mut EntityCatalog,
    root: EntityId,
    name: &str,
    a: EntityId,
    b: EntityId,
) -> EntityId {
    catalog
        .declare(EntityDeclaration::interface(
            root,
            name,
            InterfacePair::unordered(a, b),
        ))
        .unwrap()
}

fn domain_set(domains: &[PhysicsDomain]) -> BTreeSet<PhysicsDomain> {
    domains.iter().copied().collect()
}

fn zero_pos() -> [QtyAny; 3] {
    [
        QtyAny::new(0.0, LENGTH),
        QtyAny::new(0.0, LENGTH),
        QtyAny::new(0.0, LENGTH),
    ]
}

fn build_slice_a_rigid_damper_heat(scenario: &Scenario) -> DraftAssembly {
    let mut cat = EntityCatalog::new();
    let root = cat.declare(EntityDeclaration::assembly("slice-a")).unwrap();
    let (body_a, surf_a) = add_part_with_surface(&mut cat, root, "rotor_body");
    let (sink_b, surf_b) = add_part_with_surface(&mut cat, root, "thermal_housing");
    let damp_iface = add_iface(&mut cat, root, "viscous_damper_route", surf_a, surf_b);

    let input = DraftAssemblyInput {
        name: "slice-a-rigid-damper-heat".to_owned(),
        root,
        catalog: cat,
        reference_state: Some(ReferenceState::Assembled),
        environment: DraftEnvironment::earth_lab_air(),
        time_plan: time_plan(1.0, 1.0e-3, 1_000.0),
        parts: vec![
            DraftPartDecl {
                name: "rotor_body".to_owned(),
                entity: body_a,
                kind: DraftPartKind::RigidBody,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Cylinder {
                    radius: QtyAny::new(0.05, LENGTH),
                    height: QtyAny::new(0.02, LENGTH),
                }),
                materials: vec![mat("mat:steel-4140")],
                physics: domain_set(&[PhysicsDomain::Mechanics]),
                fidelity: lumped_policy(),
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "thermal_housing".to_owned(),
                entity: sink_b,
                kind: DraftPartKind::Solid,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                    lx: QtyAny::new(0.15, LENGTH),
                    ly: QtyAny::new(0.15, LENGTH),
                    lz: QtyAny::new(0.08, LENGTH),
                }),
                materials: vec![mat("mat:al-6061")],
                physics: domain_set(&[PhysicsDomain::Thermal]),
                fidelity: lumped_policy(),
                uncertainty: Vec::new(),
            },
        ],
        joints: vec![JointDecl {
            name: "shaft_bearing".to_owned(),
            parent: root,
            child: body_a,
            frame: FrameId(0),
            kind: JointKind::Revolute,
            axis: [0.0, 0.0, 1.0],
            pitch: None,
            limits: None,
            damping: Some(QtyAny::dimensionless(0.25)),
            stiffness: None,
        }],
        gear_meshes: Vec::new(),
        interfaces: vec![DraftInterfaceDecl {
            name: "damper_heat_route".to_owned(),
            entity: damp_iface,
            feature_a: surf_a,
            feature_b: surf_b,
            kind: DraftInterfaceKind::ThermalContact {
                law_card: h("iface:viscous-dissipation-route"),
                third_body: Some(ThirdBodyDecl {
                    material_card: h("mat:silicone-damper-oil"),
                    film_thickness: QtyAny::new(2.0e-4, LENGTH),
                    query_point: Vec::new(),
                }),
                state: Vec::new(),
            },
        }],
        sources: vec![SourceDecl {
            name: "drive_torque".to_owned(),
            target: SourceTarget::Port {
                entity: body_a,
                port: "torque".to_owned(),
            },
            history: TimeHistory::Step {
                t_switch: QtyAny::new(0.1, TIME),
                before: QtyAny::new(0.0, FORCE),
                after: QtyAny::new(5.0, FORCE),
            },
        }],
        controllers: Vec::new(),
        probes: vec![ProbeDecl {
            name: "housing_tc".to_owned(),
            frame: FrameId(0),
            sample_rate: QtyAny::new(100.0, FREQUENCY),
            kind: ProbeKind::Thermocouple {
                location: sink_b,
                position: zero_pos(),
                time_constant: Some(QtyAny::new(0.02, TIME)),
            },
        }],
    };
    DraftAssembly::admit(input, scenario).unwrap()
}

fn build_slice_b_gas_piston(scenario: &Scenario) -> DraftAssembly {
    let mut cat = EntityCatalog::new();
    let root = cat.declare(EntityDeclaration::assembly("slice-b")).unwrap();
    let (cavity, surf_cav) = add_part_with_surface(&mut cat, root, "gas_cavity");
    let (piston, surf_pis) = add_part_with_surface(&mut cat, root, "piston");
    let (wall, surf_wall) = add_part_with_surface(&mut cat, root, "cylinder_wall");
    let p_if = add_iface(&mut cat, root, "piston_gas_face", surf_cav, surf_pis);
    let w_if = add_iface(&mut cat, root, "gas_wall_heat", surf_cav, surf_wall);

    let input = DraftAssemblyInput {
        name: "slice-b-gas-piston".to_owned(),
        root,
        catalog: cat,
        reference_state: Some(ReferenceState::Assembled),
        environment: DraftEnvironment::earth_lab_air(),
        time_plan: TimePlan {
            horizon: QtyAny::new(0.5, TIME),
            step_policy: StepPolicy::Adaptive {
                dt_initial: QtyAny::new(1.0e-4, TIME),
                dt_min: QtyAny::new(1.0e-7, TIME),
                dt_max: QtyAny::new(5.0e-3, TIME),
                relative_tol: 1.0e-6,
                order: 2,
            },
            default_sample_rate: QtyAny::new(2_000.0, FREQUENCY),
        },
        parts: vec![
            DraftPartDecl {
                name: "gas_cavity".to_owned(),
                entity: cavity,
                kind: DraftPartKind::Cavity,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Cylinder {
                    radius: QtyAny::new(0.04, LENGTH),
                    height: QtyAny::new(0.12, LENGTH),
                }),
                materials: vec![mat("mat:ideal-air-dry")],
                physics: domain_set(&[PhysicsDomain::Fluid, PhysicsDomain::Thermal]),
                fidelity: lumped_policy(),
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "piston".to_owned(),
                entity: piston,
                kind: DraftPartKind::RigidBody,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Cylinder {
                    radius: QtyAny::new(0.04, LENGTH),
                    height: QtyAny::new(0.03, LENGTH),
                }),
                materials: vec![mat("mat:al-6061")],
                physics: domain_set(&[PhysicsDomain::Mechanics]),
                fidelity: lumped_policy(),
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "cylinder_wall".to_owned(),
                entity: wall,
                kind: DraftPartKind::Solid,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Tube {
                    inner_radius: QtyAny::new(0.04, LENGTH),
                    outer_radius: QtyAny::new(0.05, LENGTH),
                    length: QtyAny::new(0.15, LENGTH),
                }),
                materials: vec![mat("mat:steel-304")],
                physics: domain_set(&[PhysicsDomain::Thermal]),
                fidelity: lumped_policy(),
                uncertainty: Vec::new(),
            },
        ],
        joints: vec![JointDecl {
            name: "piston_slide".to_owned(),
            parent: root,
            child: piston,
            frame: FrameId(0),
            kind: JointKind::Prismatic,
            axis: [0.0, 0.0, 1.0],
            pitch: None,
            limits: Some(JointLimits {
                lower: QtyAny::new(-0.08, LENGTH),
                upper: QtyAny::new(0.08, LENGTH),
            }),
            damping: None,
            stiffness: None,
        }],
        gear_meshes: Vec::new(),
        interfaces: vec![
            DraftInterfaceDecl {
                name: "piston_pressure_work".to_owned(),
                entity: p_if,
                feature_a: surf_cav,
                feature_b: surf_pis,
                kind: DraftInterfaceKind::PowerPort {
                    schema: h("schema:hydraulic-piston-work"),
                },
            },
            DraftInterfaceDecl {
                name: "gas_wall_convection".to_owned(),
                entity: w_if,
                feature_a: surf_cav,
                feature_b: surf_wall,
                kind: DraftInterfaceKind::ThermalContact {
                    law_card: h("iface:wall-convection-card"),
                    third_body: None,
                    state: Vec::new(),
                },
            },
        ],
        sources: vec![SourceDecl {
            name: "piston_load".to_owned(),
            target: SourceTarget::Port {
                entity: piston,
                port: "force".to_owned(),
            },
            history: TimeHistory::Constant(QtyAny::new(-500.0, FORCE)),
        }],
        controllers: Vec::new(),
        probes: vec![ProbeDecl {
            name: "cavity_pressure".to_owned(),
            frame: FrameId(0),
            sample_rate: QtyAny::new(2_000.0, FREQUENCY),
            kind: ProbeKind::FieldQuantity {
                target: cavity,
                position: zero_pos(),
                quantity: ParameterKey::new("pressure", PRESSURE),
            },
        }],
    };
    DraftAssembly::admit(input, scenario).unwrap()
}

fn build_slice_c_friction_heating(scenario: &Scenario) -> DraftAssembly {
    let mut cat = EntityCatalog::new();
    let root = cat.declare(EntityDeclaration::assembly("slice-c")).unwrap();
    let (slider, surf_s) = add_part_with_surface(&mut cat, root, "brake_pad");
    let (rail, surf_r) = add_part_with_surface(&mut cat, root, "steel_track");
    let tribo = add_iface(&mut cat, root, "pad_track_contact", surf_s, surf_r);

    let input = DraftAssemblyInput {
        name: "slice-c-friction-heating".to_owned(),
        root,
        catalog: cat,
        reference_state: Some(ReferenceState::Assembled),
        environment: DraftEnvironment::earth_lab_air(),
        time_plan: time_plan(0.25, 2.5e-4, 4_000.0),
        parts: vec![
            DraftPartDecl {
                name: "brake_pad".to_owned(),
                entity: slider,
                kind: DraftPartKind::RigidBody,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                    lx: QtyAny::new(0.05, LENGTH),
                    ly: QtyAny::new(0.03, LENGTH),
                    lz: QtyAny::new(0.015, LENGTH),
                }),
                materials: vec![mat("mat:ceramic-pad")],
                physics: domain_set(&[
                    PhysicsDomain::Mechanics,
                    PhysicsDomain::Thermal,
                    PhysicsDomain::Tribology,
                ]),
                fidelity: lumped_policy(),
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "steel_track".to_owned(),
                entity: rail,
                kind: DraftPartKind::Solid,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                    lx: QtyAny::new(0.5, LENGTH),
                    ly: QtyAny::new(0.06, LENGTH),
                    lz: QtyAny::new(0.02, LENGTH),
                }),
                materials: vec![mat("mat:cast-iron-rotor")],
                physics: domain_set(&[PhysicsDomain::Thermal]),
                fidelity: FidelityPolicy {
                    rungs: vec![FidelityRung::Lumped, FidelityRung::OneDField { cells: 32 }],
                    error_budget: lumped_policy().error_budget,
                },
                uncertainty: Vec::new(),
            },
        ],
        joints: vec![JointDecl {
            name: "slide_rail".to_owned(),
            parent: root,
            child: slider,
            frame: FrameId(0),
            kind: JointKind::Prismatic,
            axis: [1.0, 0.0, 0.0],
            pitch: None,
            limits: None,
            damping: None,
            stiffness: None,
        }],
        gear_meshes: Vec::new(),
        interfaces: vec![DraftInterfaceDecl {
            name: "stick_slip_friction".to_owned(),
            entity: tribo,
            feature_a: surf_s,
            feature_b: surf_r,
            kind: DraftInterfaceKind::MechanicalContact {
                law_card: h("iface:stribeck-coulomb-friction"),
                third_body: None,
                state: vec![(
                    ParameterKey::new("normal_load", FORCE),
                    QtyAny::new(200.0, FORCE),
                )],
            },
        }],
        sources: vec![SourceDecl {
            name: "oscillating_drive".to_owned(),
            target: SourceTarget::Port {
                entity: slider,
                port: "force".to_owned(),
            },
            history: TimeHistory::Harmonic {
                amplitude: QtyAny::new(150.0, FORCE),
                frequency: QtyAny::new(5.0, FREQUENCY),
                phase: QtyAny::dimensionless(0.0),
                offset: QtyAny::new(0.0, FORCE),
            },
        }],
        controllers: Vec::new(),
        probes: vec![
            ProbeDecl {
                name: "pad_tc".to_owned(),
                frame: FrameId(0),
                sample_rate: QtyAny::new(500.0, FREQUENCY),
                kind: ProbeKind::Thermocouple {
                    location: slider,
                    position: zero_pos(),
                    time_constant: None,
                },
            },
            ProbeDecl {
                name: "rail_strain".to_owned(),
                frame: FrameId(0),
                sample_rate: QtyAny::new(2_000.0, FREQUENCY),
                kind: ProbeKind::StrainGauge {
                    target: rail,
                    position: zero_pos(),
                    gauge_direction: [1.0, 0.0, 0.0],
                },
            },
        ],
    };
    DraftAssembly::admit(input, scenario).unwrap()
}

fn build_slice_d_modal_acoustics(scenario: &Scenario) -> DraftAssembly {
    let mut cat = EntityCatalog::new();
    let root = cat.declare(EntityDeclaration::assembly("slice-d")).unwrap();
    let (string, surf_st) = add_part_with_surface(&mut cat, root, "steel_string");
    let (plate, surf_pl) = add_part_with_surface(&mut cat, root, "soundboard");
    let bridge = add_iface(&mut cat, root, "bridge_coupling", surf_st, surf_pl);

    let input = DraftAssemblyInput {
        name: "slice-d-struck-string-plate".to_owned(),
        root,
        catalog: cat,
        reference_state: Some(ReferenceState::Strung),
        environment: DraftEnvironment::earth_lab_air(),
        time_plan: time_plan(0.5, 1.0 / 48_000.0, 48_000.0),
        parts: vec![
            DraftPartDecl {
                name: "steel_string".to_owned(),
                entity: string,
                kind: DraftPartKind::String,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Cylinder {
                    radius: QtyAny::new(4.5e-4, LENGTH),
                    height: QtyAny::new(0.65, LENGTH),
                }),
                materials: vec![mat("mat:music-wire-steel")],
                physics: domain_set(&[PhysicsDomain::Mechanics, PhysicsDomain::Acoustics]),
                fidelity: FidelityPolicy {
                    rungs: vec![FidelityRung::Modal { mode_count: 40 }],
                    error_budget: lumped_policy().error_budget,
                },
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "soundboard".to_owned(),
                entity: plate,
                kind: DraftPartKind::FlexibleBody,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                    lx: QtyAny::new(0.4, LENGTH),
                    ly: QtyAny::new(0.3, LENGTH),
                    lz: QtyAny::new(0.003, LENGTH),
                }),
                materials: vec![mat("mat:sitka-spruce")],
                physics: domain_set(&[PhysicsDomain::Mechanics, PhysicsDomain::Acoustics]),
                fidelity: FidelityPolicy {
                    rungs: vec![FidelityRung::Modal { mode_count: 64 }],
                    error_budget: lumped_policy().error_budget,
                },
                uncertainty: Vec::new(),
            },
        ],
        joints: Vec::new(),
        gear_meshes: Vec::new(),
        interfaces: vec![DraftInterfaceDecl {
            name: "bridge".to_owned(),
            entity: bridge,
            feature_a: surf_st,
            feature_b: surf_pl,
            kind: DraftInterfaceKind::PowerPort {
                schema: h("schema:string-bridge-force-vel"),
            },
        }],
        sources: vec![SourceDecl {
            name: "pluck_pulse".to_owned(),
            target: SourceTarget::Port {
                entity: string,
                port: "pluck".to_owned(),
            },
            history: TimeHistory::PiecewiseLinear {
                knots: vec![
                    (QtyAny::new(0.0, TIME), QtyAny::new(0.0, FORCE)),
                    (QtyAny::new(0.002, TIME), QtyAny::new(3.0, FORCE)),
                    (QtyAny::new(0.004, TIME), QtyAny::new(0.0, FORCE)),
                ],
            },
        }],
        controllers: Vec::new(),
        probes: vec![ProbeDecl {
            name: "studio_mic".to_owned(),
            frame: FrameId(0),
            sample_rate: QtyAny::new(48_000.0, FREQUENCY),
            kind: ProbeKind::Microphone {
                position: [
                    QtyAny::new(0.0, LENGTH),
                    QtyAny::new(0.0, LENGTH),
                    QtyAny::new(0.75, LENGTH),
                ],
                aim_direction: [0.0, 0.0, -1.0],
                pattern: MicrophonePattern::Cardioid,
            },
        }],
    };
    DraftAssembly::admit(input, scenario).unwrap()
}

fn build_slice_e_field_conduction(scenario: &Scenario) -> DraftAssembly {
    let mut cat = EntityCatalog::new();
    let root = cat.declare(EntityDeclaration::assembly("slice-e")).unwrap();
    let (slab_a, surf_a) = add_part_with_surface(&mut cat, root, "copper_slab");
    let (slab_b, surf_b) = add_part_with_surface(&mut cat, root, "steel_slab");
    let contact = add_iface(&mut cat, root, "mapped_thermal_trace", surf_a, surf_b);

    let fem_policy = FidelityPolicy {
        rungs: vec![FidelityRung::FiniteElement {
            order: 1,
            h_target: QtyAny::new(0.005, LENGTH),
        }],
        error_budget: FidelityErrorBudget {
            relative_tol: 1.0e-4,
            absolute_tol: QtyAny::new(1.0e-3, TEMPERATURE),
            max_work_units: Some(10_000_000),
        },
    };

    let input = DraftAssemblyInput {
        name: "slice-e-two-slab-fem".to_owned(),
        root,
        catalog: cat,
        reference_state: Some(ReferenceState::Assembled),
        environment: DraftEnvironment::earth_lab_air(),
        time_plan: time_plan(21.6, 0.1, 10.0),
        parts: vec![
            DraftPartDecl {
                name: "copper_slab".to_owned(),
                entity: slab_a,
                kind: DraftPartKind::Solid,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                    lx: QtyAny::new(0.05, LENGTH),
                    ly: QtyAny::new(0.05, LENGTH),
                    lz: QtyAny::new(0.02, LENGTH),
                }),
                materials: vec![mat("mat:copper-c110")],
                physics: domain_set(&[PhysicsDomain::Thermal]),
                fidelity: fem_policy.clone(),
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "steel_slab".to_owned(),
                entity: slab_b,
                kind: DraftPartKind::Solid,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                    lx: QtyAny::new(0.05, LENGTH),
                    ly: QtyAny::new(0.05, LENGTH),
                    lz: QtyAny::new(0.02, LENGTH),
                }),
                materials: vec![mat("mat:steel-304")],
                physics: domain_set(&[PhysicsDomain::Thermal]),
                fidelity: fem_policy,
                uncertainty: Vec::new(),
            },
        ],
        joints: vec![JointDecl {
            name: "bolted_stack".to_owned(),
            parent: slab_a,
            child: slab_b,
            frame: FrameId(0),
            kind: JointKind::Fixed,
            axis: [0.0, 0.0, 0.0],
            pitch: None,
            limits: None,
            damping: None,
            stiffness: None,
        }],
        gear_meshes: Vec::new(),
        interfaces: vec![DraftInterfaceDecl {
            name: "slab_contact_trace".to_owned(),
            entity: contact,
            feature_a: surf_a,
            feature_b: surf_b,
            kind: DraftInterfaceKind::ThermalContact {
                law_card: h("iface:cu-steel-contact-conductance"),
                third_body: None,
                state: Vec::new(),
            },
        }],
        sources: vec![SourceDecl {
            name: "heater_flux".to_owned(),
            target: SourceTarget::Port {
                entity: slab_a,
                port: "x_min".to_owned(),
            },
            history: TimeHistory::Constant(QtyAny::new(250.0, POWER)),
        }],
        controllers: Vec::new(),
        probes: vec![ProbeDecl {
            name: "cold_face_tc".to_owned(),
            frame: FrameId(0),
            sample_rate: QtyAny::new(10.0, FREQUENCY),
            kind: ProbeKind::Thermocouple {
                location: slab_b,
                position: zero_pos(),
                time_constant: None,
            },
        }],
    };
    DraftAssembly::admit(input, scenario).unwrap()
}

fn build_heatsink_project(scenario: &Scenario) -> DraftAssembly {
    let mut cat = EntityCatalog::new();
    let root = cat
        .declare(EntityDeclaration::assembly("heatsink"))
        .unwrap();
    let (die, surf_die) = add_part_with_surface(&mut cat, root, "silicon_die");
    let (sink, surf_sink) = add_part_with_surface(&mut cat, root, "finned_heatsink");
    let (duct, surf_duct) = add_part_with_surface(&mut cat, root, "cooling_channel");
    let tim_if = add_iface(&mut cat, root, "die_to_sink_tim", surf_die, surf_sink);
    let conv_if = add_iface(&mut cat, root, "sink_to_air", surf_sink, surf_duct);
    let elec_if = add_iface(&mut cat, root, "die_power_rail", surf_die, surf_sink);

    let input = DraftAssemblyInput {
        name: "electrothermal-heatsink".to_owned(),
        root,
        catalog: cat,
        reference_state: Some(ReferenceState::Assembled),
        environment: DraftEnvironment::earth_lab_air(),
        time_plan: time_plan(60.0, 0.25, 4.0),
        parts: vec![
            DraftPartDecl {
                name: "silicon_die".to_owned(),
                entity: die,
                kind: DraftPartKind::Solid,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                    lx: QtyAny::new(0.02, LENGTH),
                    ly: QtyAny::new(0.02, LENGTH),
                    lz: QtyAny::new(0.001, LENGTH),
                }),
                materials: vec![mat("mat:silicon-mono")],
                physics: domain_set(&[PhysicsDomain::Thermal, PhysicsDomain::Electromagnetics]),
                fidelity: FidelityPolicy {
                    rungs: vec![
                        FidelityRung::Lumped,
                        FidelityRung::FiniteElement {
                            order: 2,
                            h_target: QtyAny::new(0.001, LENGTH),
                        },
                    ],
                    error_budget: lumped_policy().error_budget,
                },
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "finned_heatsink".to_owned(),
                entity: sink,
                kind: DraftPartKind::Solid,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                    lx: QtyAny::new(0.10, LENGTH),
                    ly: QtyAny::new(0.10, LENGTH),
                    lz: QtyAny::new(0.04, LENGTH),
                }),
                materials: vec![mat("mat:al-6063")],
                physics: domain_set(&[PhysicsDomain::Thermal]),
                fidelity: FidelityPolicy {
                    rungs: vec![FidelityRung::FiniteElement {
                        order: 1,
                        h_target: QtyAny::new(0.003, LENGTH),
                    }],
                    error_budget: lumped_policy().error_budget,
                },
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "cooling_channel".to_owned(),
                entity: duct,
                kind: DraftPartKind::Duct,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                    lx: QtyAny::new(0.12, LENGTH),
                    ly: QtyAny::new(0.12, LENGTH),
                    lz: QtyAny::new(0.06, LENGTH),
                }),
                materials: vec![mat("mat:ambient-air")],
                physics: domain_set(&[PhysicsDomain::Fluid, PhysicsDomain::Thermal]),
                fidelity: FidelityPolicy {
                    rungs: vec![FidelityRung::LatticeBoltzmann {
                        dx: QtyAny::new(0.002, LENGTH),
                    }],
                    error_budget: lumped_policy().error_budget,
                },
                uncertainty: Vec::new(),
            },
        ],
        joints: Vec::new(),
        gear_meshes: Vec::new(),
        interfaces: vec![
            DraftInterfaceDecl {
                name: "tim_bondline".to_owned(),
                entity: tim_if,
                feature_a: surf_die,
                feature_b: surf_sink,
                kind: DraftInterfaceKind::ThermalContact {
                    law_card: h("iface:thermal-grease-tim"),
                    third_body: Some(ThirdBodyDecl {
                        material_card: h("mat:silicone-tim-paste"),
                        film_thickness: QtyAny::new(5.0e-5, LENGTH),
                        query_point: Vec::new(),
                    }),
                    state: Vec::new(),
                },
            },
            DraftInterfaceDecl {
                name: "fin_convection".to_owned(),
                entity: conv_if,
                feature_a: surf_sink,
                feature_b: surf_duct,
                kind: DraftInterfaceKind::ThermalContact {
                    law_card: h("iface:conjugate-convection"),
                    third_body: None,
                    state: Vec::new(),
                },
            },
            DraftInterfaceDecl {
                name: "die_contact_pad".to_owned(),
                entity: elec_if,
                feature_a: surf_die,
                feature_b: surf_sink,
                kind: DraftInterfaceKind::ElectricalContact {
                    law_card: h("iface:solder-contact-resistance"),
                    third_body: None,
                    state: Vec::new(),
                },
            },
        ],
        sources: vec![SourceDecl {
            name: "die_tdp".to_owned(),
            target: SourceTarget::Port {
                entity: die,
                port: "heat_gen".to_owned(),
            },
            history: TimeHistory::Constant(QtyAny::new(120.0, POWER)),
        }],
        controllers: vec![ControllerDecl {
            name: "fan_pid".to_owned(),
            sensor_ports: vec![(die, "temp_signal".to_owned())],
            actuator_port: (duct, "fan_speed_cmd".to_owned()),
            sample_rate: QtyAny::new(10.0, FREQUENCY),
            law: ControllerLaw::Pid {
                kp: QtyAny::dimensionless(0.2),
                ki: QtyAny::dimensionless(0.05),
                kd: QtyAny::dimensionless(0.0),
                setpoint: QtyAny::new(348.15, TEMPERATURE),
            },
        }],
        probes: vec![ProbeDecl {
            name: "junction_tc".to_owned(),
            frame: FrameId(0),
            sample_rate: QtyAny::new(4.0, FREQUENCY),
            kind: ProbeKind::Thermocouple {
                location: die,
                position: zero_pos(),
                time_constant: None,
            },
        }],
    };
    DraftAssembly::admit(input, scenario).unwrap()
}

fn build_clarinet_note(scenario: &Scenario) -> DraftAssembly {
    let mut cat = EntityCatalog::new();
    let root = cat
        .declare(EntityDeclaration::assembly("clarinet"))
        .unwrap();
    let (bore, surf_bore) = add_part_with_surface(&mut cat, root, "grenadilla_bore");
    let (reed, surf_reed) = add_part_with_surface(&mut cat, root, "cane_reed");
    let (keys, surf_keys) = add_part_with_surface(&mut cat, root, "tone_holes");
    let mouth_if = add_iface(&mut cat, root, "reed_mouthpiece", surf_reed, surf_bore);
    let hole_if = add_iface(&mut cat, root, "tone_hole_ports", surf_keys, surf_bore);

    let mut env = DraftEnvironment::earth_lab_air();
    env.relative_humidity = 0.85;
    env.gas_composition = vec![
        GasComponent {
            species: "N2".to_owned(),
            mole_fraction: 0.74,
        },
        GasComponent {
            species: "O2".to_owned(),
            mole_fraction: 0.16,
        },
        GasComponent {
            species: "CO2".to_owned(),
            mole_fraction: 0.04,
        },
        GasComponent {
            species: "H2O".to_owned(),
            mole_fraction: 0.06,
        },
    ];

    let input = DraftAssemblyInput {
        name: "clarinet-a4-note".to_owned(),
        root,
        catalog: cat,
        reference_state: Some(ReferenceState::Played),
        environment: env,
        time_plan: time_plan(0.5, 1.0 / 48_000.0, 48_000.0),
        parts: vec![
            DraftPartDecl {
                name: "grenadilla_bore".to_owned(),
                entity: bore,
                kind: DraftPartKind::Duct,
                geometry: DraftGeometrySource::Standard(GeometrySource::Parametric {
                    generator: "clarinet-cylindrical-bore".to_owned(),
                    params: vec![
                        (
                            ParameterKey::new("length", LENGTH),
                            QtyAny::new(0.66, LENGTH),
                        ),
                        (
                            ParameterKey::new("bore_radius", LENGTH),
                            QtyAny::new(0.0074, LENGTH),
                        ),
                    ],
                }),
                materials: vec![mat("mat:african-blackwood")],
                physics: domain_set(&[PhysicsDomain::Acoustics, PhysicsDomain::Fluid]),
                fidelity: FidelityPolicy {
                    rungs: vec![FidelityRung::OneDField { cells: 128 }],
                    error_budget: lumped_policy().error_budget,
                },
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "cane_reed".to_owned(),
                entity: reed,
                kind: DraftPartKind::SpringBeam,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Box {
                    lx: QtyAny::new(0.035, LENGTH),
                    ly: QtyAny::new(0.013, LENGTH),
                    lz: QtyAny::new(0.0005, LENGTH),
                }),
                materials: vec![mat("mat:arundo-donax-reed")],
                physics: domain_set(&[PhysicsDomain::Mechanics, PhysicsDomain::Acoustics]),
                fidelity: FidelityPolicy {
                    rungs: vec![FidelityRung::Modal { mode_count: 6 }],
                    error_budget: lumped_policy().error_budget,
                },
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "tone_holes".to_owned(),
                entity: keys,
                kind: DraftPartKind::ValvePattern,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Cylinder {
                    radius: QtyAny::new(0.004, LENGTH),
                    height: QtyAny::new(0.005, LENGTH),
                }),
                materials: vec![mat("mat:nickel-silver-pad")],
                physics: domain_set(&[PhysicsDomain::Acoustics]),
                fidelity: lumped_policy(),
                uncertainty: Vec::new(),
            },
        ],
        joints: Vec::new(),
        gear_meshes: Vec::new(),
        interfaces: vec![
            DraftInterfaceDecl {
                name: "reed_lay_contact".to_owned(),
                entity: mouth_if,
                feature_a: surf_reed,
                feature_b: surf_bore,
                kind: DraftInterfaceKind::MechanicalContact {
                    law_card: h("iface:beating-reed-lay"),
                    third_body: None,
                    state: Vec::new(),
                },
            },
            DraftInterfaceDecl {
                name: "tone_hole_radiation".to_owned(),
                entity: hole_if,
                feature_a: surf_keys,
                feature_b: surf_bore,
                kind: DraftInterfaceKind::PowerPort {
                    schema: h("schema:acoustic-tone-hole-junction"),
                },
            },
        ],
        sources: vec![SourceDecl {
            name: "mouth_pressure".to_owned(),
            target: SourceTarget::Port {
                entity: reed,
                port: "blowing_pressure".to_owned(),
            },
            history: TimeHistory::PiecewiseLinear {
                knots: vec![
                    (QtyAny::new(0.0, TIME), QtyAny::new(0.0, PRESSURE)),
                    (QtyAny::new(0.02, TIME), QtyAny::new(3_500.0, PRESSURE)),
                    (QtyAny::new(0.45, TIME), QtyAny::new(3_500.0, PRESSURE)),
                    (QtyAny::new(0.50, TIME), QtyAny::new(0.0, PRESSURE)),
                ],
            },
        }],
        controllers: Vec::new(),
        probes: vec![ProbeDecl {
            name: "bell_mic".to_owned(),
            frame: FrameId(0),
            sample_rate: QtyAny::new(48_000.0, FREQUENCY),
            kind: ProbeKind::Microphone {
                position: [
                    QtyAny::new(0.0, LENGTH),
                    QtyAny::new(0.0, LENGTH),
                    QtyAny::new(1.0, LENGTH),
                ],
                aim_direction: [0.0, 0.0, -1.0],
                pattern: MicrophonePattern::Omnidirectional,
            },
        }],
    };
    DraftAssembly::admit(input, scenario).unwrap()
}

fn build_mechanical_watch(scenario: &Scenario) -> DraftAssembly {
    let mut cat = EntityCatalog::new();
    let root = cat.declare(EntityDeclaration::assembly("watch")).unwrap();
    let (mainspring, surf_ms) = add_part_with_surface(&mut cat, root, "mainspring");
    let (barrel_gear, surf_bg) = add_part_with_surface(&mut cat, root, "barrel_wheel");
    let (escape_wheel, surf_ew) = add_part_with_surface(&mut cat, root, "escape_wheel");
    let (balance, surf_bal) = add_part_with_surface(&mut cat, root, "balance_wheel");
    let (hairspring, surf_hs) = add_part_with_surface(&mut cat, root, "hairspring");

    let ms_if = add_iface(&mut cat, root, "mainspring_barrel_drive", surf_ms, surf_bg);
    let esc_if = add_iface(&mut cat, root, "escapement_impulse", surf_ew, surf_bal);
    let hs_if = add_iface(&mut cat, root, "hairspring_collet", surf_hs, surf_bal);

    let mut env = DraftEnvironment::earth_lab_air();
    // 1.5 mT ambient stray magnetic field for antimagnetic chronometer testing.
    env.magnetic_field = [
        QtyAny::new(1.5e-3, MAGNETIC_FIELD),
        QtyAny::new(0.0, MAGNETIC_FIELD),
        QtyAny::new(0.0, MAGNETIC_FIELD),
    ];

    let input = DraftAssemblyInput {
        name: "mechanical-chronometer-watch".to_owned(),
        root,
        catalog: cat,
        reference_state: Some(ReferenceState::Assembled),
        environment: env,
        time_plan: time_plan(10.0, 1.0e-4, 1_000.0),
        parts: vec![
            DraftPartDecl {
                name: "mainspring".to_owned(),
                entity: mainspring,
                kind: DraftPartKind::SpringBeam,
                geometry: DraftGeometrySource::SpiralSpring(SpiralSpringGenerator {
                    inner_radius: QtyAny::new(0.002, LENGTH),
                    outer_radius: QtyAny::new(0.006, LENGTH),
                    turns: 12.0,
                    blade_thickness: QtyAny::new(0.00012, LENGTH),
                    blade_height: QtyAny::new(0.0014, LENGTH),
                }),
                materials: vec![mat("mat:nivaflex-mainspring")],
                physics: domain_set(&[PhysicsDomain::Mechanics]),
                fidelity: lumped_policy(),
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "barrel_wheel".to_owned(),
                entity: barrel_gear,
                kind: DraftPartKind::RigidBody,
                geometry: DraftGeometrySource::Gear(GearGenerator {
                    teeth: 80,
                    module: QtyAny::new(1.5e-4, LENGTH),
                    face_width: QtyAny::new(4.0e-4, LENGTH),
                    pressure_angle: QtyAny::dimensionless(20.0_f64.to_radians()),
                    profile_shift: 0.0,
                    bore_radius: QtyAny::new(1.0e-3, LENGTH),
                }),
                materials: vec![mat("mat:cu-zn-watch-brass")],
                physics: domain_set(&[PhysicsDomain::Mechanics]),
                fidelity: lumped_policy(),
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "escape_wheel".to_owned(),
                entity: escape_wheel,
                kind: DraftPartKind::RigidBody,
                geometry: DraftGeometrySource::Gear(GearGenerator {
                    teeth: 15,
                    module: QtyAny::new(1.5e-4, LENGTH),
                    face_width: QtyAny::new(2.5e-4, LENGTH),
                    pressure_angle: QtyAny::dimensionless(20.0_f64.to_radians()),
                    profile_shift: 0.0,
                    bore_radius: QtyAny::new(3.0e-4, LENGTH),
                }),
                materials: vec![mat("mat:hardened-steel-escape")],
                physics: domain_set(&[PhysicsDomain::Mechanics, PhysicsDomain::Tribology]),
                fidelity: lumped_policy(),
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "balance_wheel".to_owned(),
                entity: balance,
                kind: DraftPartKind::RigidBody,
                geometry: DraftGeometrySource::Primitive(PrimitiveShape::Tube {
                    inner_radius: QtyAny::new(0.004, LENGTH),
                    outer_radius: QtyAny::new(0.0052, LENGTH),
                    length: QtyAny::new(0.0006, LENGTH),
                }),
                materials: vec![mat("mat:glucydur-balance")],
                physics: domain_set(&[
                    PhysicsDomain::Mechanics,
                    PhysicsDomain::Electromagnetics,
                    PhysicsDomain::Thermal,
                ]),
                fidelity: lumped_policy(),
                uncertainty: Vec::new(),
            },
            DraftPartDecl {
                name: "hairspring".to_owned(),
                entity: hairspring,
                kind: DraftPartKind::SpringBeam,
                geometry: DraftGeometrySource::SpiralSpring(SpiralSpringGenerator {
                    inner_radius: QtyAny::new(0.0006, LENGTH),
                    outer_radius: QtyAny::new(0.0032, LENGTH),
                    turns: 14.0,
                    blade_thickness: QtyAny::new(0.000035, LENGTH),
                    blade_height: QtyAny::new(0.00018, LENGTH),
                }),
                materials: vec![mat("mat:nivarox-hairspring")],
                physics: domain_set(&[
                    PhysicsDomain::Mechanics,
                    PhysicsDomain::Thermal,
                    PhysicsDomain::Electromagnetics,
                ]),
                fidelity: FidelityPolicy {
                    rungs: vec![FidelityRung::Lumped, FidelityRung::Modal { mode_count: 12 }],
                    error_budget: lumped_policy().error_budget,
                },
                uncertainty: Vec::new(),
            },
        ],
        joints: vec![
            JointDecl {
                name: "barrel_arbor".to_owned(),
                parent: root,
                child: barrel_gear,
                frame: FrameId(0),
                kind: JointKind::Revolute,
                axis: [0.0, 0.0, 1.0],
                pitch: None,
                limits: None,
                damping: Some(QtyAny::dimensionless(1.0e-7)),
                stiffness: None,
            },
            JointDecl {
                name: "escape_pivot".to_owned(),
                parent: root,
                child: escape_wheel,
                frame: FrameId(0),
                kind: JointKind::Revolute,
                axis: [0.0, 0.0, 1.0],
                pitch: None,
                limits: None,
                damping: Some(QtyAny::dimensionless(2.0e-8)),
                stiffness: None,
            },
            JointDecl {
                name: "balance_staff".to_owned(),
                parent: root,
                child: balance,
                frame: FrameId(0),
                kind: JointKind::Revolute,
                axis: [0.0, 0.0, 1.0],
                pitch: None,
                limits: Some(JointLimits {
                    lower: QtyAny::dimensionless(-6.0),
                    upper: QtyAny::dimensionless(6.0),
                }),
                damping: Some(QtyAny::dimensionless(5.0e-9)),
                stiffness: None,
            },
        ],
        gear_meshes: vec![GearMeshDecl {
            name: "barrel_to_escape_train".to_owned(),
            driver: barrel_gear,
            driven: escape_wheel,
            teeth_driver: 80,
            teeth_driven: 15,
            pressure_angle: QtyAny::dimensionless(20.0_f64.to_radians()),
            center_distance: QtyAny::new(0.007125, LENGTH),
            backlash: Some(QtyAny::new(5.0e-6, LENGTH)),
            mesh_card: h("iface:involute-watch-gear-mesh"),
        }],
        interfaces: vec![
            DraftInterfaceDecl {
                name: "mainspring_torque".to_owned(),
                entity: ms_if,
                feature_a: surf_ms,
                feature_b: surf_bg,
                kind: DraftInterfaceKind::PowerPort {
                    schema: h("schema:rotational-torque-omega"),
                },
            },
            DraftInterfaceDecl {
                name: "swiss_lever_escapement".to_owned(),
                entity: esc_if,
                feature_a: surf_ew,
                feature_b: surf_bal,
                kind: DraftInterfaceKind::MechanicalContact {
                    law_card: h("iface:ruby-pallet-steel-impulse"),
                    third_body: Some(ThirdBodyDecl {
                        material_card: h("mat:moebius-9415-epilame-grease"),
                        film_thickness: QtyAny::new(2.0e-6, LENGTH),
                        query_point: Vec::new(),
                    }),
                    state: Vec::new(),
                },
            },
            DraftInterfaceDecl {
                name: "hairspring_collet_glue".to_owned(),
                entity: hs_if,
                feature_a: surf_hs,
                feature_b: surf_bal,
                kind: DraftInterfaceKind::Glue {
                    law_card: h("iface:rigid-collet-clamp"),
                },
            },
        ],
        sources: Vec::new(),
        controllers: Vec::new(),
        probes: vec![
            ProbeDecl {
                name: "beat_rate_timer".to_owned(),
                frame: FrameId(0),
                sample_rate: QtyAny::new(10_000.0, FREQUENCY),
                kind: ProbeKind::Timing {
                    target: balance,
                    threshold: QtyAny::dimensionless(0.0),
                    edge: TimingEdge::Both,
                },
            },
            ProbeDecl {
                name: "tick_acoustic_mic".to_owned(),
                frame: FrameId(0),
                sample_rate: QtyAny::new(48_000.0, FREQUENCY),
                kind: ProbeKind::Microphone {
                    position: [
                        QtyAny::new(0.0, LENGTH),
                        QtyAny::new(0.0, LENGTH),
                        QtyAny::new(0.05, LENGTH),
                    ],
                    aim_direction: [0.0, 0.0, -1.0],
                    pattern: MicrophonePattern::Supercardioid,
                },
            },
        ],
    };
    DraftAssembly::admit(input, scenario).unwrap()
}

#[test]
fn all_eight_target_devices_admit_and_round_trip_with_zero_escape_hatches() {
    let scenario = Scenario::new("draft-validation", 111, Environment::earth_lab());
    let devices = [
        build_slice_a_rigid_damper_heat(&scenario),
        build_slice_b_gas_piston(&scenario),
        build_slice_c_friction_heating(&scenario),
        build_slice_d_modal_acoustics(&scenario),
        build_slice_e_field_conduction(&scenario),
        build_heatsink_project(&scenario),
        build_clarinet_note(&scenario),
        build_mechanical_watch(&scenario),
    ];

    assert_eq!(devices.len(), 8);
    let mut hashes = BTreeSet::new();
    for dev in &devices {
        let rt = dev
            .round_trip(&scenario)
            .expect("admitted draft must round-trip cleanly");
        assert_eq!(rt, *dev);
        assert_eq!(rt.content_hash(), dev.content_hash());
        assert!(hashes.insert(dev.content_hash()));
    }

    // Parametric mechanical/primitive generators lower cleanly into v1 GeometrySource tables.
    let gear = GearGenerator {
        teeth: 24,
        module: QtyAny::new(1.0e-3, LENGTH),
        face_width: QtyAny::new(5.0e-3, LENGTH),
        pressure_angle: QtyAny::dimensionless(20.0_f64.to_radians()),
        profile_shift: 0.0,
        bore_radius: QtyAny::new(4.0e-3, LENGTH),
    };
    assert!((gear.pitch_radius_m() - 0.012).abs() < 1e-15);
    assert!(matches!(
        gear.to_geometry_source(),
        GeometrySource::Parametric { ref generator, .. } if generator == "Gear.draft.v1"
    ));
}

#[test]
fn draft_assembly_refuses_invalid_inputs_by_name() {
    let scenario = Scenario::new("draft-refusals", 111, Environment::earth_lab());
    let base = build_mechanical_watch(&scenario).input().clone();

    // 1. Missing reference state.
    let mut bad = base.clone();
    bad.reference_state = None;
    let err = DraftAssembly::admit(bad, &scenario).unwrap_err();
    assert_eq!(err.item, "assembly.reference_state");

    // 2. Gas composition mole fractions not summing to 1.
    let mut bad = base.clone();
    bad.environment.gas_composition[0].mole_fraction = 0.5;
    let err = DraftAssembly::admit(bad, &scenario).unwrap_err();
    assert_eq!(err.item, "environment.gas_composition");

    // 3. Invalid magnetic field dimensions.
    let mut bad = base.clone();
    bad.environment.magnetic_field[0] = QtyAny::new(1.0, LENGTH);
    let err = DraftAssembly::admit(bad, &scenario).unwrap_err();
    assert!(err.item.contains("magnetic_field"), "{err}");

    // 4. Relative humidity outside [0, 1].
    let mut bad = base.clone();
    bad.environment.relative_humidity = 1.25;
    let err = DraftAssembly::admit(bad, &scenario).unwrap_err();
    assert_eq!(err.item, "environment.relative_humidity");

    // 5. Overlapping spiral spring coils.
    let mut bad = base.clone();
    if let DraftGeometrySource::SpiralSpring(ref mut sp) = bad.parts[0].geometry {
        sp.turns = 100.0;
    }
    let err = DraftAssembly::admit(bad, &scenario).unwrap_err();
    assert!(
        err.item.contains("mainspring") && err.item.contains("blade_thickness"),
        "{err}"
    );

    // 6. Gear teeth < 3.
    let mut bad = base.clone();
    if let DraftGeometrySource::Gear(ref mut g) = bad.parts[1].geometry {
        g.teeth = 2;
    }
    let err = DraftAssembly::admit(bad, &scenario).unwrap_err();
    assert!(
        err.item.contains("barrel_wheel") && err.item.contains("teeth"),
        "{err}"
    );

    // 7. Zero joint axis on Revolute joint.
    let mut bad = base.clone();
    bad.joints[0].axis = [0.0, 0.0, 0.0];
    let err = DraftAssembly::admit(bad, &scenario).unwrap_err();
    assert!(err.item.contains("barrel_arbor"), "{err}");

    // 8. Helical joint without pitch.
    let mut bad = base.clone();
    bad.joints[0].kind = JointKind::Helical;
    bad.joints[0].pitch = None;
    let err = DraftAssembly::admit(bad, &scenario).unwrap_err();
    assert!(err.item.contains("barrel_arbor"), "{err}");

    // 9. Zero interface card hash.
    let mut bad = base.clone();
    if let DraftInterfaceKind::MechanicalContact {
        ref mut law_card, ..
    } = bad.interfaces[1].kind
    {
        *law_card = ContentHash([0; 32]);
    }
    let err = DraftAssembly::admit(bad, &scenario).unwrap_err();
    assert!(err.item.contains("swiss_lever_escapement"), "{err}");

    // 10. Part with empty physics set.
    let mut bad = base.clone();
    bad.parts[0].physics.clear();
    let err = DraftAssembly::admit(bad, &scenario).unwrap_err();
    assert!(err.item.contains("mainspring"), "{err}");

    // 11. Zero microphone aim vector.
    let mut bad = base.clone();
    if let ProbeKind::Microphone {
        ref mut aim_direction,
        ..
    } = bad.probes[1].kind
    {
        *aim_direction = [0.0, 0.0, 0.0];
    }
    let err = DraftAssembly::admit(bad, &scenario).unwrap_err();
    assert!(err.item.contains("tick_acoustic_mic"), "{err}");

    // 12. Time-step exceeding horizon.
    let mut bad = base.clone();
    bad.time_plan.step_policy = StepPolicy::Fixed {
        dt: QtyAny::new(20.0, TIME),
    };
    let err = DraftAssembly::admit(bad, &scenario).unwrap_err();
    assert!(err.item.contains("time_plan"), "{err}");

    // 13. Non-positive base ambient temperature.
    let mut bad = base.clone();
    bad.environment.base.ambient_temperature = QtyAny::new(-5.0, TEMPERATURE);
    let err = DraftAssembly::admit(bad, &scenario).unwrap_err();
    assert_eq!(err.item, "environment.base.ambient_temperature");

    // 14. Joint referencing a surface instead of root/declared part.
    let mut bad = base.clone();
    bad.joints[0].parent = bad.interfaces[0].feature_a;
    let err = DraftAssembly::admit(bad, &scenario).unwrap_err();
    assert!(err.item.contains("barrel_arbor"), "{err}");

    // 15. Parameter-targeted source whose history dimensions disagree with parameter dims.
    let mut bad = base;
    bad.sources.push(SourceDecl {
        name: "bad_param_src".to_owned(),
        target: SourceTarget::Parameter(fs_scenario::assembly::ParameterPath {
            entity: bad.parts[0].entity,
            parameter: ParameterKey::new("temp_bias", TEMPERATURE),
        }),
        history: TimeHistory::Constant(QtyAny::new(1.0, FORCE)),
    });
    let err = DraftAssembly::admit(bad, &scenario).unwrap_err();
    assert_eq!(err.item, "source:bad_param_src");
}
