//! G0/G3/G4: actual material-resolved field, native rotor momentum and paired acceptance.
//! Synthetic coefficients are numerical fixtures, not measured machine parameters.

use fs_couple::planar_induction::{
    CircularConductor, PlanarInductionRotor, RotorControls, RotorError,
};
use fs_em::planar::{
    PlanarMagnetostatics,
    eddy::{
        EddyConductor, EddyControls, PlanarEddy,
        moving::{MovingEddyControls, mechanical::PlanarLorentzLoadLimits},
    },
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_ga::{Mat3, Se3, So3, Vec3};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::{MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY},
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_mbd::articulated::{ArticulatedModel, JointModel, Link, SpatialInertia};
use fs_qty::{QuantitySpec, Temperature};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::CheckedCgConfig;

const DEPTH: f64 = 0.3;
const SIGMA: f64 = 6.0;
const MU: f64 = 0.2;
const INERTIA: f64 = 0.2;

fn source(sigma: f64) -> ResolvedMaterialStatePoint {
    let basis = ElectricalPropertyBasis::Conductivity;
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (basis.property_name(), basis.dimensions(), sigma),
        (
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            MU,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    200.0,
                    1000.0,
                ),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic native induction rotor".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "rotor-fixture".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300.0)
            .unwrap(),
        &[
            basis.requirement().unwrap(),
            ScalarPropertyRequirement::try_new(
                MAGNETIC_PERMEABILITY_PROPERTY,
                MAGNETIC_PERMEABILITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn linear() -> CheckedCgConfig {
    CheckedCgConfig {
        tolerance: 1e-11,
        max_iterations: 1000,
        max_corrections: 2,
    }
}
fn field_controls() -> MovingEddyControls {
    MovingEddyControls {
        relative_tolerance: 1e-11,
        absolute_residual_a_m: 1e-12,
        restart: 32,
        max_cycles: 40,
        max_krylov_scalars: 100000,
        max_cell_magnetic_peclet: 30.0,
        energy_tolerance_j: 1e-9,
        current_balance_tolerance_a: 1e-9,
    }
}
fn controls() -> RotorControls {
    RotorControls {
        maximum_iterations: 64,
        relaxation: 0.7,
        angular_velocity_tolerance_rad_s: 1e-10,
        energy_tolerance_j: 2e-9,
        angular_impulse_tolerance_nms: 1e-10,
        linear_impulse_tolerance_ns: 1e-10,
    }
}
fn limits() -> PlanarLorentzLoadLimits {
    PlanarLorentzLoadLimits {
        max_vertices: 256,
        max_cells: 512,
    }
}
fn section(inner: f64) -> CircularConductor {
    CircularConductor {
        inner_radius_m: inner,
        outer_radius_m: 1.0,
        maximum_boundary_error_m: 0.025,
    }
}
fn model(inertia: f64, axis_x: f64) -> ArticulatedModel {
    let inertia = SpatialInertia::new(
        2.0,
        Vec3::default(),
        Mat3 {
            m: [inertia, 0.0, 0.0, 0.0, inertia, 0.0, 0.0, 0.0, inertia],
        },
    )
    .unwrap();
    ArticulatedModel::new(vec![Link::new(
        "native electromagnetic rotor",
        None,
        Se3::from_parts(So3::identity(), Vec3::new(axis_x, 0.0, DEPTH / 2.0)).unwrap(),
        JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
        inertia,
    )])
    .unwrap()
}
fn problem(inner: f64, sectors: usize, nonuniform: bool) -> PlanarMagnetostatics {
    let mut vertices = Vec::new();
    let radii = if inner == 0.0 {
        vertices.push([0.0; 2]);
        vec![0.5, 1.0]
    } else {
        vec![inner, (inner + 1.0) / 2.0, 1.0]
    };
    let offset = vertices.len();
    for radius in &radii {
        for i in 0..sectors {
            let angle = std::f64::consts::TAU * i as f64 / sectors as f64;
            vertices.push([radius * angle.cos(), radius * angle.sin()]);
        }
    }
    let mut faces = Vec::new();
    if inner == 0.0 {
        for i in 0..sectors {
            faces.push([0, (1 + i) as u32, (1 + (i + 1) % sectors) as u32]);
        }
    }
    for ring in 0..radii.len() - 1 {
        for i in 0..sectors {
            let a = (offset + ring * sectors + i) as u32;
            let b = (offset + (ring + 1) * sectors + i) as u32;
            let c = (offset + (ring + 1) * sectors + (i + 1) % sectors) as u32;
            let d = (offset + ring * sectors + (i + 1) % sectors) as u32;
            faces.extend([[a, b, c], [a, c, d]]);
        }
    }
    let mesh = TriComplex2::from_indexed_triangles(
        "fixed circular conductor approximation",
        vertices,
        faces,
        Metric2::planar(DEPTH).unwrap(),
    )
    .unwrap();
    let mut materials = vec![source(SIGMA); mesh.faces().len()];
    if nonuniform {
        materials[0] = source(SIGMA * 2.0);
    }
    PlanarMagnetostatics::new(
        mesh,
        materials,
        WeightedAssemblyLimits {
            max_cells: 512,
            max_dofs: 1024,
            max_triplets: 10000,
        },
        &mut || false,
    )
    .unwrap()
}
fn history(problem: &PlanarMagnetostatics, inner: f64) -> PlanarEddy<'_> {
    let boundary: Vec<_> = problem
        .mesh()
        .vertices()
        .iter()
        .enumerate()
        .filter(|(_, p)| {
            let r = p[0].hypot(p[1]);
            (r - 1.0).abs() < 1e-12 || inner > 0.0 && (r - inner).abs() < 1e-12
        })
        .map(|(i, p)| (i, p[1]))
        .collect();
    let zero = vec![0.0; problem.mesh().faces().len()];
    let field = problem
        .solve(&zero, &boundary, linear(), &mut || false)
        .unwrap();
    let conductors = problem
        .materials()
        .iter()
        .enumerate()
        .map(|(triangle, source)| EddyConductor {
            triangle,
            region: 7,
            conductor: ResolvedConductorState::try_from_material_with_basis(
                source.clone(),
                ElectricalPropertyBasis::Conductivity,
            )
            .unwrap(),
        })
        .collect();
    PlanarEddy::new(
        &field,
        conductors,
        &[],
        &boundary,
        0.0,
        EddyControls {
            linear: linear(),
            energy_tolerance_j: 1e-9,
            current_balance_tolerance_a: 1e-9,
        },
        &mut || false,
    )
    .unwrap()
}
fn rotor(problem: &PlanarMagnetostatics, inner: f64, speed: f64) -> PlanarInductionRotor<'_> {
    PlanarInductionRotor::new(
        history(problem, inner),
        model(INERTIA, 0.0),
        section(inner),
        0.37,
        speed,
        limits(),
        &mut || false,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() < tolerance,
        "{a} != {b}, tolerance {tolerance}"
    );
}

#[test]
fn finite_inertia_braking_uses_same_native_field_and_mechanical_midpoint() {
    let problem = problem(0.0, 16, false);
    let zero = vec![0.0; problem.mesh().faces().len()];
    for initial_speed in [2.0, -2.0] {
        let original = rotor(&problem, 0.0, initial_speed);
        let duration = 0.02;
        let trial = original
            .step_candidate(
                duration,
                &zero,
                0.0,
                field_controls(),
                controls(),
                &mut || false,
            )
            .unwrap();
        assert!(trial.induction.report.as_ref().unwrap().iters > 0);
        assert!(trial.balance.iterations > 1);
        assert!(trial.state().angular_velocity_rad_s().abs() < initial_speed.abs());
        assert!(trial.state().angular_velocity_rad_s() * initial_speed > 0.0);
        assert!(trial.balance.electromagnetic_torque_nm * initial_speed < 0.0);
        assert!(trial.induction.heat_j > 0.0);
        assert!(trial.induction.stored_after_j > trial.induction.stored_before_j + 1e-10);
        assert!(
            trial
                .state()
                .field()
                .potential_wb_m()
                .iter()
                .zip(original.field().potential_wb_m())
                .any(|(a, b)| (a - b).abs() > 1e-7)
        );
        near(
            trial.balance.field_angular_velocity_rad_s,
            0.5 * (initial_speed + trial.state().angular_velocity_rad_s()),
            1.1e-10,
        );
        near(
            trial.mechanics.angular_momentum_change_world.z,
            duration * trial.balance.electromagnetic_torque_nm,
            1e-11,
        );
        near(
            trial.mechanics.angular_momentum_change_world.z,
            INERTIA * (trial.state().angular_velocity_rad_s() - initial_speed),
            1e-11,
        );
        near(
            trial.mechanics.kinetic_energy_change + trial.induction.stored_after_j
                - trial.induction.stored_before_j
                + trial.induction.heat_j,
            0.0,
            3e-9,
        );
        near(trial.balance.total_energy_residual_j, 0.0, 2e-9);
        near(trial.mechanics.support_work, 0.0, 1e-14);
        // An actual native solve held at the OLD speed gives a different finite
        // inertia response: merely exporting a torque cannot pass this check.
        let velocity: Vec<_> = problem
            .mesh()
            .vertices()
            .iter()
            .map(|p| [-initial_speed * p[1], initial_speed * p[0]])
            .collect();
        let old_speed_field = original
            .field()
            .step_with_eulerian_velocity_candidate(
                duration,
                &velocity,
                &zero,
                field_controls(),
                &mut || false,
            )
            .unwrap();
        let selected: Vec<_> = (0..zero.len()).collect();
        let held_torque = old_speed_field
            .lorentz_loads(&selected, limits(), &mut || false)
            .unwrap()
            .wrench(&Se3::identity(), [0.0, 0.0, DEPTH / 2.0])
            .unwrap()
            .torque
            .z;
        let one_way_speed = initial_speed + duration * held_torque / INERTIA;
        assert!((trial.state().angular_velocity_rad_s() - one_way_speed).abs() > 1e-4);
        let state = trial.into_state();
        let continued = state
            .step_candidate(0.04, &zero, 0.0, field_controls(), controls(), &mut || {
                false
            })
            .unwrap();
        assert_eq!(continued.state().field().time_s(), 0.04);
        assert!(continued.induction.heat_j > 0.0);
        assert_eq!(original.field().time_s(), 0.0);
        assert_eq!(original.angular_velocity_rad_s(), initial_speed);
    }
}

#[test]
fn annular_rotor_drive_retains_shaft_work_support_impulse_and_field_history() {
    let problem = problem(0.35, 16, false);
    let original = rotor(&problem, 0.35, 0.0);
    let zero = vec![0.0; problem.mesh().faces().len()];
    let trial = original
        .step_candidate(0.02, &zero, 1.0, field_controls(), controls(), &mut || {
            false
        })
        .unwrap();
    assert!(trial.state().angular_velocity_rad_s() > 0.0);
    assert!(trial.state().angular_velocity_rad_s() < 0.02 / INERTIA);
    assert!(trial.induction.heat_j > 0.0 && trial.balance.shaft_work_j > 0.0);
    near(trial.mechanics.support_impulse_world.torque.z, 0.02, 1e-12);
    near(
        trial.mechanics.angular_momentum_change_world.z,
        0.02 * (trial.balance.electromagnetic_torque_nm + 1.0),
        1e-12,
    );
    near(
        trial.mechanics.kinetic_energy_change + trial.induction.stored_after_j
            - trial.induction.stored_before_j
            + trial.induction.heat_j,
        trial.balance.shaft_work_j,
        3e-9,
    );
    near(trial.induction.supplied_work_j, 0.0, 1e-14);
    assert!(original.boundary_error_m() < section(0.35).maximum_boundary_error_m);
}

#[test]
fn noninvariant_geometry_and_failed_coupling_publish_no_field_or_rotor_state() {
    let square = problem(0.0, 4, false);
    assert!(
        PlanarInductionRotor::new(
            history(&square, 0.0),
            model(INERTIA, 0.0),
            section(0.0),
            0.0,
            2.0,
            limits(),
            &mut || false
        )
        .is_err()
    );
    let mixed = problem(0.0, 16, true);
    assert!(
        PlanarInductionRotor::new(
            history(&mixed, 0.0),
            model(INERTIA, 0.0),
            section(0.0),
            0.0,
            2.0,
            limits(),
            &mut || false
        )
        .is_err()
    );
    let problem = problem(0.0, 16, false);
    assert!(
        PlanarInductionRotor::new(
            history(&problem, 0.0),
            model(INERTIA, 0.1),
            section(0.0),
            0.0,
            2.0,
            limits(),
            &mut || false
        )
        .is_err()
    );
    let original = rotor(&problem, 0.0, 2.0);
    let potential = original.field().potential_wb_m().to_vec();
    let zero = vec![0.0; problem.mesh().faces().len()];
    let exhausted = RotorControls {
        maximum_iterations: 1,
        ..controls()
    };
    assert!(matches!(
        original.step_candidate(0.02, &zero, 0.0, field_controls(), exhausted, &mut || false),
        Err(RotorError::Nonconvergence { .. })
    ));
    let mut internal_source = zero.clone();
    internal_source[0] = 1.0;
    assert!(
        original
            .step_candidate(
                0.02,
                &internal_source,
                0.0,
                field_controls(),
                controls(),
                &mut || false
            )
            .is_err()
    );
    let mut successful_polls = 0;
    original
        .step_candidate(0.02, &zero, 0.0, field_controls(), controls(), &mut || {
            successful_polls += 1;
            false
        })
        .unwrap();
    let mut polls = 0;
    let cancelled =
        original.step_candidate(0.02, &zero, 0.0, field_controls(), controls(), &mut || {
            polls += 1;
            polls == successful_polls
        });
    assert!(matches!(cancelled, Err(RotorError::Cancelled)));
    assert_eq!(polls, successful_polls);
    assert_eq!(original.field().potential_wb_m(), potential);
    assert_eq!(original.field().time_s(), 0.0);
    assert_eq!(original.angular_velocity_rad_s(), 2.0);
}
