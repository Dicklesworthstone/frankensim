//! Resolved fields, native body momentum and independent physical work checks.
//! All material cards are synthetic numerical fixtures, not measured devices.

use fs_couple::electrostatic::motion::{
    ElectricalControl, Error, MotionControls, MotionInput, PlanarCapacitorMotion,
    PlanarCapacitorPath,
};
use fs_em::electrostatic::{ElectrostaticControls, PlanarElectrostatics};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_ga::Vec3;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY},
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_mbd_material::{RegionBinding, ResolvedBody, VolumeMoments, resolve_body};
use fs_qty::{Density, Dims};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::CheckedCgConfig;

fn card(name: &str, dims: Dims, value: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic resolved capacitor motion".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "numerical fixture".into(),
            phase: "solid".into(),
            process: "uniform".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn body(mass: f64) -> ResolvedBody {
    let source = card("density", Density::DIMS, mass);
    resolve_body(
        &[RegionBinding {
            key: "moving-electrode",
            material: &source,
            point: &QueryPoint::new(),
            selection: MaterialPropertySelection::SingleClaimOnly,
            geometry: VolumeMoments::cuboid(Vec3::new(1.0, 1.0, 1.0)).unwrap(),
        }],
        1,
        || false,
    )
    .unwrap()
}
fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 64,
        max_dofs: 128,
        max_triplets: 2048,
    }
}
fn field_controls() -> ElectrostaticControls {
    ElectrostaticControls {
        linear: CheckedCgConfig {
            tolerance: 1e-12,
            max_iterations: 256,
            max_corrections: 2,
        },
        charge_tolerance_c: 1e-10,
        energy_tolerance_j: 1e-10,
    }
}
fn controls() -> MotionControls {
    MotionControls {
        displacement_bracket_m: [-0.08, 0.08],
        displacement_tolerance_m: 1e-13,
        energy_tolerance_j: 1e-10,
        maximum_iterations: 96,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}
fn reference(left: f64, right: f64, tapered: bool) -> PlanarElectrostatics {
    let mut positions = vec![
        [0.0, 0.0],
        [0.5, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [0.5, 1.0],
        [1.0, 1.0],
    ];
    if tapered {
        positions[2][1] = 0.2;
        positions[5][1] = 0.8;
    }
    let mesh = TriComplex2::from_indexed_triangles(
        "resolved-moving-capacitor",
        positions,
        vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]],
        Metric2::Planar { thickness: 1.0 },
    )
    .unwrap();
    let materials = [left, left, right, right]
        .into_iter()
        .map(|epsilon| {
            let source = card(
                ELECTRIC_PERMITTIVITY_PROPERTY,
                ELECTRIC_PERMITTIVITY_DIMS,
                epsilon,
            );
            resolve_material_state_point(
                &source,
                &QueryPoint::new(),
                &[ScalarPropertyRequirement::try_new(
                    ELECTRIC_PERMITTIVITY_PROPERTY,
                    ELECTRIC_PERMITTIVITY_DIMS,
                    ScalarAdmissibility::StrictlyPositive,
                )
                .unwrap()],
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap()
        })
        .collect();
    PlanarElectrostatics::new(mesh, materials, limits(), &mut || false).unwrap()
}
fn path(epsilon: f64) -> PlanarCapacitorPath {
    PlanarCapacitorPath::new(
        reference(epsilon, epsilon, false),
        [vec![2, 5], vec![0, 3]],
        vec![
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
        ],
        [-0.4, 0.4],
        limits(),
        field_controls(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn isolated_charge_drives_native_material_body_with_exact_uniform_gap_motion() {
    for (epsilon, mass, charge) in [
        (2.0, 1.0, 0.4),
        (4.0, 1.0, -0.4),
        (2.0, 3.0, 0.4),
        (2.0, 1.0, 0.0),
    ] {
        let state = PlanarCapacitorMotion::new(
            path(epsilon),
            body(mass),
            ElectricalControl::IsolatedCharge { charge_c: charge },
            0.0,
            0.02,
            &mut || false,
        )
        .unwrap();
        let before = (
            state.coordinate_m(),
            state.velocity_m_s(),
            state.field_energy_j(),
        );
        let input = MotionInput {
            duration_s: 0.1,
            external_force_n: 0.01,
        };
        let force = -charge * charge / (2.0 * epsilon);
        let acceleration = (force + input.external_force_n) / mass;
        let expected =
            0.02 * input.duration_s + 0.5 * acceleration * input.duration_s * input.duration_s;
        let candidate = state
            .step_candidate(input, controls(), &mut || false)
            .unwrap();
        let next = candidate.state();
        let step = candidate.step();
        near(state.field_force_n(), force, 1e-11);
        near(next.coordinate_m(), expected, 1e-11);
        near(
            next.velocity_m_s(),
            0.02 + input.duration_s * acceleration,
            1e-11,
        );
        near(
            next.capacitance().capacitance_f(),
            epsilon / (1.0 + expected),
            1e-10,
        );
        near(next.voltage_v(), charge * (1.0 + expected) / epsilon, 1e-11);
        near(step.discrete_field_force_n, force, 1e-10);
        near(step.field_work_on_body_j, force * expected, 1e-11);
        near(
            step.native.kinetic_energy_change,
            0.5 * mass * (next.velocity_m_s().powi(2) - 0.02_f64.powi(2)),
            1e-11,
        );
        near(
            step.native.linear_momentum_change_world.x,
            mass * input.duration_s * acceleration,
            1e-11,
        );
        assert_eq!(step.electrical_supply_j, 0.0);
        near(
            step.total_energy_residual_j,
            0.0,
            controls().energy_tolerance_j,
        );
        assert_eq!(
            before,
            (
                state.coordinate_m(),
                state.velocity_m_s(),
                state.field_energy_j()
            )
        );
    }
}

#[test]
fn held_voltage_accounts_for_source_charge_work_in_attraction_and_energy_return() {
    let epsilon = 2.0;
    let voltage = 0.7;
    let dt = 0.04;
    let state = PlanarCapacitorMotion::new(
        path(epsilon),
        body(1.0),
        ElectricalControl::HeldVoltage { voltage_v: voltage },
        0.0,
        0.0,
        &mut || false,
    )
    .unwrap();
    let candidate = state
        .step_candidate(
            MotionInput {
                duration_s: dt,
                external_force_n: 0.0,
            },
            controls(),
            &mut || false,
        )
        .unwrap();
    // Exact root of the midpoint/discrete-gradient equations for C=epsilon/g.
    let k = 0.25 * dt * dt * epsilon * voltage * voltage;
    let displacement = -2.0 * k / (1.0 + (1.0 - 4.0 * k).sqrt());
    near(candidate.state().coordinate_m(), displacement, 1e-11);
    near(
        candidate.state().velocity_m_s(),
        2.0 * displacement / dt,
        1e-11,
    );
    let step = candidate.step();
    assert!(step.electrical_supply_j > 0.0 && step.field_work_on_body_j > 0.0);
    near(
        step.electrical_supply_j,
        2.0 * step.field_energy_change_j,
        1e-11,
    );
    near(step.field_work_on_body_j, step.field_energy_change_j, 1e-11);
    near(
        step.electrical_supply_j,
        voltage * (candidate.state().charge_c() - state.charge_c()),
        1e-11,
    );

    let receding = PlanarCapacitorMotion::new(
        path(epsilon),
        body(1.0),
        ElectricalControl::HeldVoltage { voltage_v: voltage },
        0.0,
        0.1,
        &mut || false,
    )
    .unwrap();
    let next = receding
        .step_candidate(
            MotionInput {
                duration_s: 0.01,
                external_force_n: 0.0,
            },
            controls(),
            &mut || false,
        )
        .unwrap();
    assert!(next.state().coordinate_m() > 0.0);
    assert!(next.step().electrical_supply_j < 0.0);
    assert!(next.step().field_work_on_body_j < 0.0);
    near(
        next.step().total_energy_residual_j,
        0.0,
        controls().energy_tolerance_j,
    );
}

#[test]
fn real_layer_geometry_and_nonparallel_fields_supply_the_shape_derivative() {
    let geometry = PlanarCapacitorPath::new(
        reference(2.0, 5.0, false),
        [vec![2, 5], vec![0, 3]],
        vec![
            [0.0, 0.0],
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [1.0, 0.0],
        ],
        [-0.2, 0.2],
        limits(),
        field_controls(),
        &mut || false,
    )
    .unwrap();
    for coordinate in [-0.1, 0.0, 0.13] {
        let sample = geometry.sample(coordinate, &mut || false).unwrap();
        let denominator = 0.5 / 2.0 + (0.5 + coordinate) / 5.0;
        near(sample.capacitance_f(), 1.0 / denominator, 1e-10);
        near(
            sample.derivative_f_per_m(),
            -1.0 / (5.0 * denominator * denominator),
            1e-10,
        );
        let rebuilt = geometry.problem_at(coordinate, &mut || false).unwrap();
        assert_eq!(
            rebuilt.mesh().lineage(),
            geometry.reference().mesh().lineage()
        );
        assert_eq!(rebuilt.materials(), geometry.reference().materials());
    }
    let tapered = PlanarCapacitorPath::new(
        reference(2.0, 5.0, true),
        [vec![2, 5], vec![0, 3]],
        vec![
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
        ],
        [-0.3, 0.3],
        limits(),
        field_controls(),
        &mut || false,
    )
    .unwrap();
    let coordinate = 0.12;
    let h = 1e-5;
    let sample = tapered.sample(coordinate, &mut || false).unwrap();
    let a = tapered
        .sample(coordinate - h, &mut || false)
        .unwrap()
        .capacitance_f();
    let b = tapered
        .sample(coordinate + h, &mut || false)
        .unwrap()
        .capacitance_f();
    near(sample.derivative_f_per_m(), (b - a) / (2.0 * h), 1e-7);
}

#[test]
fn failed_geometry_feedback_and_late_cancellation_leave_reusable_state() {
    let state = PlanarCapacitorMotion::new(
        path(2.0),
        body(1.0),
        ElectricalControl::IsolatedCharge { charge_c: 0.4 },
        0.0,
        0.02,
        &mut || false,
    )
    .unwrap();
    let input = MotionInput {
        duration_s: 0.1,
        external_force_n: 0.01,
    };
    let accepted = (state.coordinate_m(), state.velocity_m_s(), state.charge_c());
    let mut polls = 0;
    assert!(matches!(
        state.step_candidate(input, controls(), &mut || {
            polls += 1;
            polls >= 20
        }),
        Err(Error::Cancelled)
    ));
    let mut too_narrow = controls();
    too_narrow.displacement_bracket_m = [-1e-8, 1e-8];
    assert!(matches!(
        state.step_candidate(input, too_narrow, &mut || false),
        Err(Error::OutsideBracket { .. })
    ));
    assert!(state.path().sample(0.5, &mut || false).is_err());
    assert_eq!(
        accepted,
        (state.coordinate_m(), state.velocity_m_s(), state.charge_c())
    );
    let first = state
        .step_candidate(input, controls(), &mut || false)
        .unwrap();
    let replay = state
        .step_candidate(input, controls(), &mut || false)
        .unwrap();
    assert_eq!(first.state().coordinate_m(), replay.state().coordinate_m());
    assert_eq!(first.state().velocity_m_s(), replay.state().velocity_m_s());
    let (next, _) = first.into_parts();
    let continued = next
        .step_candidate(input, controls(), &mut || false)
        .unwrap();
    assert!(continued.state().coordinate_m() > next.coordinate_m());
    near(
        continued.step().total_energy_residual_j,
        0.0,
        controls().energy_tolerance_j,
    );
}

#[test]
fn path_rejects_folding_and_inconsistent_rigid_electrode_motion() {
    let mut motion = vec![
        [0.0, 0.0],
        [10.0, 0.0],
        [1.0, 0.0],
        [0.0, 0.0],
        [10.0, 0.0],
        [1.0, 0.0],
    ];
    assert!(
        PlanarCapacitorPath::new(
            reference(2.0, 2.0, false),
            [vec![2, 5], vec![0, 3]],
            motion.clone(),
            [-0.2, 0.2],
            limits(),
            field_controls(),
            &mut || false
        )
        .is_err()
    );
    motion[1] = [0.5, 0.0];
    motion[4] = [0.5, 0.0];
    motion[5] = [0.9, 0.0];
    assert!(
        PlanarCapacitorPath::new(
            reference(2.0, 2.0, false),
            [vec![2, 5], vec![0, 3]],
            motion,
            [-0.2, 0.2],
            limits(),
            field_controls(),
            &mut || false
        )
        .is_err()
    );
}
