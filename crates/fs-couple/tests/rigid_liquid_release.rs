//! G1/G3/G4: native rigid matter and motion enter the original thermal VOF.
//! Synthetic properties test physical accounting, not measured melt behavior.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::liquid_transfer::rigid::{
    RigidLiquidReleaseControls, RigidLiquidReleaseError, RigidTetrahedralLiquidState,
    release_rigid_liquid,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ga::{Se3, So3, So3Tangent, Twist, Vec3 as GaVec3};
use fs_lbm::d3q19::freesurface3::ContactModel3;
use fs_lbm::material::free_surface::enthalpy::FreeSurfaceEnthalpyControls;
use fs_lbm::material::free_surface::transfer::TetLiquidTransferLimits;
use fs_lbm::material::free_surface::transfer::enthalpy::TetEnthalpyTransferLimits;
use fs_lbm::material::free_surface::{MaterialFreeSurfaceDomain3, plan_free_surface_from_material};
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::surface::{LiquidCapillarity, resolve_liquid_surface_state};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd::RigidBodyState;
use fs_mbd_material::density_field::{DensityVertexBinding, TetrahedralDensityBody};
use fs_mbd_material::principal::PrincipalFrame;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};

const DX: f64 = 0.001;
const DT: f64 = 0.0002;
const RHO: f64 = 1000.0;
const TOPOLOGY: [[usize; 4]; 6] = [
    [0, 1, 3, 7],
    [0, 3, 2, 7],
    [0, 2, 6, 7],
    [0, 6, 4, 7],
    [0, 4, 5, 7],
    [0, 5, 1, 7],
];

fn with_cx<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 41,
                kernel_id: 28,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", temperature)
        .unwrap()
        .with("pressure", 100000.0)
        .unwrap()
}
fn card(label: &str, rho: [f64; 2]) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        ("density", Dims([-3, 1, 0, 0, 0, 0]), rho),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0]), [0.1, 0.02]),
        (
            "specific_heat_capacity",
            Dims([2, 0, -2, -1, 0, 0]),
            [2000.0; 2],
        ),
        ("thermal_conductivity", Dims([1, 1, -3, -1, 0, 0]), [0.5; 2]),
    ] {
        let quantity = if name == "thermal_conductivity" {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (500.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 300.0, 500.0)
                    .with("pressure", 100000.0, 100000.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic rigid-release conservation fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: label.into(),
            phase: "authored phase chart".into(),
            process: "analytic fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn chart(card: &MaterialCard, phase: SolidLiquidPhase) -> EquilibriumEnthalpyPhaseCurve {
    let fraction = if phase == SolidLiquidPhase::Liquid {
        1.0
    } else {
        0.0
    };
    EquilibriumEnthalpyPhaseCurve::try_single_phase(
        card.content_hash(),
        phase,
        vec![
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 0.0,
                temperature_k: 300.0,
                liquid_mass_fraction: fraction,
                bulk_density_kg_m3: RHO,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 400000.0,
                temperature_k: 500.0,
                liquid_mass_fraction: fraction,
                bulk_density_kg_m3: RHO,
            },
        ],
    )
    .unwrap()
}
fn controls() -> RigidLiquidReleaseControls {
    RigidLiquidReleaseControls {
        transfer: TetEnthalpyTransferLimits {
            mechanical: TetLiquidTransferLimits {
                max_tetrahedra: 16,
                max_overlap_tests: 30000,
                max_wet_cells: 4096,
                relative_volume_tolerance: 1e-10,
                mass_tolerance_kg: 1e-14,
                momentum_tolerance_kg_m_s: 1e-16,
                angular_momentum_tolerance_kg_m2_s: 1e-18,
                maximum_centroid_shift_m: 1e-10,
                maximum_projection_energy_change_j: 1e-6,
            },
            max_source_queries: 2048,
            enthalpy_tolerance_j: 1e-10,
        },
        thermal: FreeSurfaceEnthalpyControls::default(),
        source_kinetic_tolerance_j: 1e-15,
        source_centroid_tolerance_m: 1e-11,
        rheology: None,
        buoyancy: None,
    }
}
fn domain(gravity: [f64; 3]) -> MaterialFreeSurfaceDomain3<'static> {
    MaterialFreeSurfaceDomain3 {
        cells: [16; 3],
        gravity_m_s2: gravity,
        contact: ContactModel3::Neutral,
        solid_tiles: &[],
        max_cells: 4096,
    }
}
fn rotate(vector: GaVec3) -> So3 {
    So3::exp(So3Tangent::new(vector)).unwrap()
}
fn array(v: GaVec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vector(v: [f64; 3]) -> GaVec3 {
    GaVec3::new(v[0], v[1], v[2])
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}
fn close_vector(a: [f64; 3], b: GaVec3, tolerance: f64) {
    for (a, b) in a.into_iter().zip(array(b)) {
        close(a, b, tolerance);
    }
}
fn points(center: GaVec3, lengths: GaVec3, rotation: So3) -> Vec<GaVec3> {
    (0..8)
        .map(|i| {
            let local = GaVec3::new(
                (if i & 1 == 0 { -0.5 } else { 0.5 }) * lengths.x,
                (if i & 2 == 0 { -0.5 } else { 0.5 }) * lengths.y,
                (if i & 4 == 0 { -0.5 } else { 0.5 }) * lengths.z,
            );
            center + rotation.rotate(local).unwrap()
        })
        .collect()
}
fn geometry(
    card: &MaterialCard,
    points: &[GaVec3],
    topology: &[[usize; 4]],
    temperatures: &[f64],
) -> TetrahedralDensityBody {
    let query: Vec<_> = temperatures.iter().map(|&t| point(t)).collect();
    let bindings: Vec<_> = points
        .iter()
        .zip(&query)
        .map(|(&position_m, point)| DensityVertexBinding {
            position_m,
            material: card,
            point,
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
        .collect();
    TetrahedralDensityBody::resolve(&bindings, topology, points.len(), || false).unwrap()
}
fn motion(
    frame: &PrincipalFrame,
    world_center: GaVec3,
    rotation: So3,
    world_velocity: GaVec3,
    world_spin: GaVec3,
) -> RigidBodyState {
    let center = frame.inertia().center_of_mass();
    let pose = Se3::from_parts(rotation, world_center - rotation.rotate(center).unwrap()).unwrap();
    let omega = rotation.inverse().rotate(world_spin).unwrap();
    let v = rotation.inverse().rotate(world_velocity).unwrap() - omega.cross(center);
    frame
        .state_from_kinematics(pose, Twist::new(omega, v))
        .unwrap()
}
fn state<'a>(
    geometry: &'a TetrahedralDensityBody,
    frame: &'a PrincipalFrame,
    motion: RigidBodyState,
    h: &'a [f64],
) -> RigidTetrahedralLiquidState<'a> {
    RigidTetrahedralLiquidState {
        geometry,
        frame,
        motion,
        specific_enthalpy_j_kg: h,
        time_s: 1.25,
    }
}

#[test]
fn g1_full_inertia_world_spin_heat_and_clock_survive_release_and_native_continuation() {
    let card = card("rigid cuboid liquid", [RHO; 2]);
    let query = point(400.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &query,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = resolve_liquid_surface_state(
        &card,
        &query,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = chart(&card, SolidLiquidPhase::Liquid);
    let lengths = GaVec3::new(0.004, 0.006, 0.002);
    let center = GaVec3::new(0.003, -0.006, 0.005);
    let authored_rotation = rotate(GaVec3::new(0.3, -0.5, 0.2));
    let p = points(center, lengths, authored_rotation);
    let body = geometry(&card, &p, &TOPOLOGY, &[400.0; 8]);
    let frame = body.principal_frame(1e-12, || false).unwrap();
    assert!(body.inertia().inertia_com().m[1].abs() > 1e-13);
    let h: Vec<_> = p
        .iter()
        .map(|p| 200000.0 + 1e6 * (p.x - center.x))
        .collect();
    let untouched = (body.clone(), h.clone());
    let world_rotation = rotate(GaVec3::new(-0.2, 0.4, 0.1));
    let world_spin = GaVec3::new(1.7, -2.4, 3.1);
    let mass = RHO * lengths.x * lengths.y * lengths.z;
    // Independent rectangular-volume integrals, before either principal lowering
    // or native tetrahedral transfer is consulted.
    let moments = GaVec3::new(
        mass * (lengths.y.powi(2) + lengths.z.powi(2)) / 12.0,
        mass * (lengths.x.powi(2) + lengths.z.powi(2)) / 12.0,
        mass * (lengths.x.powi(2) + lengths.y.powi(2)) / 12.0,
    );
    let omega_box = authored_rotation
        .inverse()
        .rotate(world_rotation.inverse().rotate(world_spin).unwrap())
        .unwrap();
    let spin = world_rotation
        .rotate(
            authored_rotation
                .rotate(GaVec3::new(
                    moments.x * omega_box.x,
                    moments.y * omega_box.y,
                    moments.z * omega_box.z,
                ))
                .unwrap(),
        )
        .unwrap();
    let gravity = [0.2, -0.5, -2.0];
    let gate = CancelGate::new_clock_free();
    for (origin, boost, buoyant) in [
        ([0.0; 3], [0.0; 3], false),
        ([-0.01, 0.02, -0.03], [0.02, 0.01, -0.01], true),
    ] {
        let com = vector(origin) + GaVec3::new(0.008, 0.008, 0.008);
        let velocity = GaVec3::new(0.01, -0.006, 0.003) + vector(boost);
        let motion = motion(&frame, com, world_rotation, velocity, world_spin);
        let source = state(&body, &frame, motion, &h);
        let mut controls = controls();
        if buoyant {
            controls.buoyancy = Some(BoussinesqConfig {
                gravity_m_s2: gravity,
                max_relative_density_departure: 0.0,
                max_relative_specific_heat_departure: 0.0,
            });
        }
        with_cx(&gate, |cx| {
            let scaling = plan_free_surface_from_material(&surface, DX, DT, 0.1).unwrap();
            let (mut liquid, report) = release_rigid_liquid(
                source,
                scaling,
                domain(gravity),
                origin,
                &law,
                &caloric,
                controls,
                cx,
                &gate,
            )
            .unwrap();
            let expected_p = velocity.scale(mass);
            let expected_l = (com - vector(origin)).cross(expected_p) + spin;
            let expected_k = 0.5 * mass * velocity.dot(velocity) + 0.5 * world_spin.dot(spin);
            for actual in [report.rigid, report.transfer.mechanical.source] {
                close(actual.mass_kg, mass, 2e-18);
                close_vector(actual.momentum_kg_m_s, expected_p, 2e-19);
                close_vector(actual.angular_momentum_kg_m2_s, expected_l, 2e-21);
                close(actual.kinetic_energy_j, expected_k, 2e-20);
            }
            let mut actual_p = GaVec3::new(0.0, 0.0, 0.0);
            let mut actual_l = actual_p;
            let (mut actual_mass, mut actual_h) = (0.0, 0.0);
            for piece in &report.transfer.mechanical.components {
                for &coordinates in &piece.target_cells {
                    let [x, y, z] = coordinates;
                    let sample = liquid.flow().sample_at_si(x, y, z).unwrap().unwrap();
                    let mass = sample.density_kg_m3 * sample.fill_fraction * DX.powi(3);
                    let momentum = vector(sample.velocity_m_s).scale(mass);
                    let r = vector(coordinates.map(|i| (i as f64 + 0.5) * DX));
                    actual_mass += mass;
                    actual_p = actual_p + momentum;
                    actual_l = actual_l + r.cross(momentum);
                    actual_h += mass * liquid.specific_enthalpy_at(coordinates).unwrap();
                }
            }
            close(actual_mass, mass, 2e-18);
            close_vector(array(actual_p), expected_p, 2e-19);
            close_vector(array(actual_l), expected_l, 2e-21);
            close(actual_h, mass * 200000.0, 2e-11);
            close(report.rigid_enthalpy_j, mass * 200000.0, 2e-13);
            close(report.enthalpy_error_j, 0.0, 2e-11);
            assert!(report.kinetic_projection_change_j.abs() > 1e-14);
            close(
                report.kinetic_projection_change_j,
                report.transfer.mechanical.target.kinetic_energy_j - expected_k,
                2e-20,
            );
            close(liquid.physical_time_s(), source.time_s, 0.0);
            assert_eq!(liquid.accepted_steps(), 0);
            let step = liquid.step(&gate).unwrap();
            close(step.enthalpy_before_j, mass * 200000.0, 2e-11);
            close(step.external_heat_j, 0.0, 0.0);
            close(step.viscous_heat_j, 0.0, 0.0);
            close(liquid.physical_time_s(), source.time_s + DT, 0.0);
            assert_eq!(source.motion, motion);
        });
    }
    assert_eq!((body, h), untouched);
}

#[test]
fn g1_releasing_detached_rigidly_carried_pieces_preserves_each_pieces_motion_and_heat() {
    let card = card("two constrained liquid parcels", [RHO; 2]);
    let query = point(400.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &query,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = resolve_liquid_surface_state(
        &card,
        &query,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = chart(&card, SolidLiquidPhase::Liquid);
    let side = 0.0024;
    let mut p = Vec::new();
    for x in [-0.003, 0.003] {
        p.extend(points(
            GaVec3::new(x, 0.0, 0.0),
            GaVec3::new(side, side, side),
            So3::identity(),
        ));
    }
    let mut topology = TOPOLOGY.to_vec();
    topology.extend(TOPOLOGY.map(|tet| tet.map(|i| i + 8)));
    let body = geometry(&card, &p, &topology, &[400.0; 16]);
    let frame = body.principal_frame(1e-12, || false).unwrap();
    let com = GaVec3::new(0.008, 0.008, 0.008);
    let velocity = GaVec3::new(0.01, -0.005, 0.002);
    let omega = GaVec3::new(0.0, 0.0, 3.0);
    let motion = motion(&frame, com, So3::identity(), velocity, omega);
    let h: Vec<_> = (0..16)
        .map(|i| if i < 8 { 200000.0 } else { 220000.0 })
        .collect();
    let source = state(&body, &frame, motion, &h);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, DX, DT, 0.1).unwrap();
        let (_, report) = release_rigid_liquid(
            source,
            scaling,
            domain([0.0; 3]),
            [0.0; 3],
            &law,
            &caloric,
            controls(),
            cx,
            &gate,
        )
        .unwrap();
        assert_eq!(report.transfer.mechanical.components.len(), 2);
        for (index, piece) in report.transfer.mechanical.components.iter().enumerate() {
            let offset = GaVec3::new(if index == 0 { -0.003 } else { 0.003 }, 0.0, 0.0);
            let mass = RHO * side.powi(3);
            let momentum = (velocity + omega.cross(offset)).scale(mass);
            let angular = (com + offset).cross(momentum) + omega.scale(mass * side.powi(2) / 6.0);
            close_vector(piece.target.momentum_kg_m_s, momentum, 2e-20);
            close_vector(piece.target.angular_momentum_kg_m2_s, angular, 2e-22);
            close(
                report.transfer.components[index].target_enthalpy_j,
                mass * h[index * 8],
                2e-12,
            );
            assert_eq!(piece.source_tetrahedra.len(), 6);
            assert!(!piece.target_cells.is_empty());
        }
        assert_ne!(
            report.transfer.mechanical.components[0]
                .target
                .momentum_kg_m_s,
            report.transfer.mechanical.components[1]
                .target
                .momentum_kg_m_s
        );
    });
}

#[test]
fn g4_foreign_matter_frame_phase_budgets_and_cancellation_preserve_the_rigid_source() {
    let card = card("retry source", [RHO; 2]);
    let query = point(400.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &query,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = resolve_liquid_surface_state(
        &card,
        &query,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = chart(&card, SolidLiquidPhase::Liquid);
    let frozen = chart(&card, SolidLiquidPhase::Solid);
    let p = points(
        GaVec3::new(0.0, 0.0, 0.0),
        GaVec3::new(0.004, 0.004, 0.004),
        So3::identity(),
    );
    let body = geometry(&card, &p, &TOPOLOGY, &[400.0; 8]);
    let frame = body.principal_frame(1e-12, || false).unwrap();
    let motion = motion(
        &frame,
        GaVec3::new(0.008, 0.008, 0.008),
        So3::identity(),
        GaVec3::new(0.01, -0.005, 0.002),
        GaVec3::new(0.0, 0.0, 3.0),
    );
    let h = [200000.0; 8];
    let source = state(&body, &frame, motion, &h);
    let untouched = (body.clone(), motion, h);
    let foreign_card = self::card("same values, different material identity", [RHO; 2]);
    let foreign_body = geometry(&foreign_card, &p, &TOPOLOGY, &[400.0; 8]);
    let foreign_frame = foreign_body.principal_frame(1e-12, || false).unwrap();
    let resized = geometry(
        &card,
        &p.iter().map(|p| p.scale(1.1)).collect::<Vec<_>>(),
        &TOPOLOGY,
        &[400.0; 8],
    );
    let wrong_frame = resized.principal_frame(1e-12, || false).unwrap();
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, DX, DT, 0.1).unwrap();
        for bad in [
            RigidTetrahedralLiquidState {
                frame: &wrong_frame,
                ..source
            },
            state(&foreign_body, &foreign_frame, motion, &h),
        ] {
            assert!(matches!(
                release_rigid_liquid(
                    bad,
                    scaling,
                    domain([0.0; 3]),
                    [0.0; 3],
                    &law,
                    &caloric,
                    controls(),
                    cx,
                    &gate
                ),
                Err(RigidLiquidReleaseError::Invalid(_))
            ));
        }
        assert!(
            release_rigid_liquid(
                source,
                scaling,
                domain([0.0; 3]),
                [0.0; 3],
                &law,
                &frozen,
                controls(),
                cx,
                &gate
            )
            .is_err()
        );
        for bad in [
            {
                let mut c = controls();
                c.transfer.mechanical.max_overlap_tests = 1;
                c
            },
            {
                let mut c = controls();
                c.transfer.mechanical.maximum_projection_energy_change_j = 0.0;
                c
            },
            {
                let mut c = controls();
                c.transfer.max_source_queries = 1;
                c
            },
        ] {
            assert!(
                release_rigid_liquid(
                    source,
                    scaling,
                    domain([0.0; 3]),
                    [0.0; 3],
                    &law,
                    &caloric,
                    bad,
                    cx,
                    &gate
                )
                .is_err()
            );
        }
        let (_, a) = release_rigid_liquid(
            source,
            scaling,
            domain([0.0; 3]),
            [0.0; 3],
            &law,
            &caloric,
            controls(),
            cx,
            &gate,
        )
        .unwrap();
        let (_, b) = release_rigid_liquid(
            source,
            scaling,
            domain([0.0; 3]),
            [0.0; 3],
            &law,
            &caloric,
            controls(),
            cx,
            &gate,
        )
        .unwrap();
        assert_eq!(a, b);
    });
    let gate = CancelGate::new_clock_free();
    gate.request();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, DX, DT, 0.1).unwrap();
        assert!(matches!(
            release_rigid_liquid(
                source,
                scaling,
                domain([0.0; 3]),
                [0.0; 3],
                &law,
                &caloric,
                controls(),
                cx,
                &gate
            ),
            Err(RigidLiquidReleaseError::Cancelled)
        ));
    });
    assert_eq!((body, motion, h), untouched);
}

#[test]
fn g3_equal_total_mass_does_not_authorize_erasing_the_original_density_distribution() {
    let card = card("graded reference matter", [900.0, 1100.0]);
    let query = point(400.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &query,
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let surface = resolve_liquid_surface_state(
        &card,
        &query,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let caloric = chart(&card, SolidLiquidPhase::Liquid);
    let side: f64 = 0.004;
    let p = points(
        GaVec3::new(0.0, 0.0, 0.0),
        GaVec3::new(side, side, side),
        So3::identity(),
    );
    let temperatures: Vec<_> = (0..8)
        .map(|i| if i & 1 == 0 { 300.0 } else { 500.0 })
        .collect();
    let body = geometry(&card, &p, &TOPOLOGY, &temperatures);
    close(body.inertia().mass(), RHO * side.powi(3), 2e-18);
    let frame = body.principal_frame(1e-12, || false).unwrap();
    let motion = motion(
        &frame,
        GaVec3::new(0.008, 0.008, 0.008),
        So3::identity(),
        GaVec3::new(0.01, 0.0, 0.0),
        GaVec3::new(0.0, 0.0, 3.0),
    );
    let h = [200000.0; 8];
    let source = state(&body, &frame, motion, &h);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let scaling = plan_free_surface_from_material(&surface, DX, DT, 0.1).unwrap();
        assert!(matches!(
            release_rigid_liquid(
                source,
                scaling,
                domain([0.0; 3]),
                [0.0; 3],
                &law,
                &caloric,
                controls(),
                cx,
                &gate
            ),
            Err(RigidLiquidReleaseError::Invalid(_))
        ));
    });
}
