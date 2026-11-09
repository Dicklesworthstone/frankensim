//! Original constrained M/C/K, virtual work and modal limiting cases.
//! Analytic mechanics is the oracle; every exported operator uses native owners.
use super::*;
use fs_math::c64::C64;
use fs_mbd::articulated::motion::reconfigure::constraints::{
    ConstraintCoordinate, JointConstraint,
};
use fs_mbd::articulated::motion::spatial_loads::prescribed::equilibrium::{
    constrained::{self as static_owner, JointEquilibriumControls},
    linearized::{self as unconstrained, CentralDifferenceControls, LinearizationControls},
};
use static_owner::linearized::{self as response, ConstrainedLinearizationControls};

fn response_controls<'a>(q: &'a [f64], v: &'a [f64], f: &'a [f64])
    -> LinearizationControls<'a> {
    LinearizationControls {
        coordinate_scales: q, velocity_scales: v, effort_scales: f,
        equilibrium_tolerance: 1e-10,
        maximum_stencil_distortion: 1e-8,
        difference: CentralDifferenceControls {
            initial_step: 1e-3, absolute_tolerance: 1e-9,
            relative_tolerance: 1e-7, maximum_refinements: 8,
            maximum_evaluations: 18,
        },
    }
}
fn subspace_controls() -> ConstrainedLinearizationControls {
    ConstrainedLinearizationControls {
        maximum_physical_dofs: 8,
        minimum_reciprocal_condition: 1e-10,
        minimum_scaled_qr_diagonal: 1e-10,
        maximum_constraint_residual: 1e-12,
        maximum_orthogonality_error: 1e-12,
    }
}
fn row(coefficients: Vec<f64>, position: f64) -> JointConstraint {
    JointConstraint {
        coordinate: ConstraintCoordinate::Length, coefficients, position,
        position_tolerance: 1e-10, velocity_tolerance: 1e-10,
    }
}
fn joint_controls<'a>(reaction: &'a [f64], position: &'a [f64])
    -> JointEquilibriumControls<'a> {
    JointEquilibriumControls {
        reaction_scales: reaction, position_scales: position,
        minimum_reciprocal_condition: 1e-10,
    }
}
fn two_sliders() -> (ArticulatedModel, [AxialConnection; 2]) {
    let mut links = vec![Link::new("support", None, Se3::identity(), JointModel::FIXED,
        inertia(3.0, Vec3::default()))];
    let connections = core::array::from_fn(|i| {
        links.push(Link::new(format!("slider-{i}"), Some(0),
            Se3::from_parts(So3::identity(), Vec3::new(1.0, i as f64, 0.0)).unwrap(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(2.0+i as f64, Vec3::default())));
        AxialConnection {
            a: MaterialPoint { link: 0, point_body: Vec3::new(0.0, i as f64, 0.0) },
            b: MaterialPoint { link: i+1, point_body: Vec3::default() },
            rest_length_m: 1.0, stiffness_n_per_m: 100.0*(i+1) as f64,
            damping_ns_per_m: 7.0,
        }
    });
    (ArticulatedModel::new(links).unwrap(), connections)
}

#[test]
fn geared_reaction_balanced_root_has_correct_reduced_mechanics_and_work_ports() {
    let (model, connections) = two_sliders();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[10.0, 40.0],
            external_body_wrench: &[Wrench::default(); 3], gravity_world: Vec3::default(),
        }, world_points: &[], connections: &connections,
    };
    let initial = seed(vec![0.12, -0.04]);
    let before = initial.clone();
    let rows = [row(vec![-2.0, 1.0], 0.1)];
    let root = static_owner::equilibrate(&model, &initial, loads,
        EquilibriumBushings::Diagonal(&[]), &rows,
        controls(&[0.1; 2], &[40.0; 2]), joint_controls(&[20.0], &[0.1]), || false).unwrap();
    // A reaction-balanced point is NOT an unconstrained force-balanced point.
    assert!(unconstrained::prepare(&model, &root.state, loads,
        EquilibriumBushings::Diagonal(&[]),
        response_controls(&[1.0; 2], &[1.0; 2], &[40.0; 2]), || false).is_err());
    for scale in [1.0_f64, -1.0, 1e-6, 1e6] {
        let mut scaled_row = rows[0].clone();
        for coefficient in &mut scaled_row.coefficients { *coefficient *= scale; }
        scaled_row.position *= scale;
        scaled_row.position_tolerance *= scale.abs();
        scaled_row.velocity_tolerance *= scale.abs();
        let mut scaled_root = root.clone();
        scaled_root.multipliers[0] /= scale;
        // All cached reports, including reaction_effort, deliberately stay unchanged.
        let prepared = response::prepare(&model, &scaled_root, loads,
            EquilibriumBushings::Diagonal(&[]), &[scaled_row],
            response_controls(&[1.0; 2], &[1.0; 2], &[40.0; 2]),
            subspace_controls(), || false).unwrap();
        assert_eq!(prepared.subspace().dimension(), 1);
        let n = prepared.subspace().basis();
        close(n[1], 2.0*n[0], 1e-14);
        close(n[0]*n[0]+n[1]*n[1], 1.0, 1e-14);
        close(prepared.project_effort(prepared.mean_reaction_effort(), || false).unwrap()[0],
            0.0, 2e-13);
        let force = [3.0, -5.0];
        let z = [0.7];
        let physical = prepared.lift_perturbation(&z, || false).unwrap();
        let drive = prepared.project_effort(&force, || false).unwrap();
        close(z[0]*drive[0], physical[0]*force[0]+physical[1]*force[1], 2e-14);
        close(prepared.inertia_action(&[2.0], || false).unwrap()[0], 5.6, 1e-12);
        let pencil = prepared.quadratic_model(1, 64, || false).unwrap();
        close(pencil.mass()[0], 2.8, 1e-12);
        close(pencil.stiffness()[0], 180.0, 2e-6);
        close(pencil.velocity_operator()[0], 7.0, 2e-6);
        let poles = pencil.poles(128, || false).unwrap();
        assert_eq!(poles.len(), 2);
        let decay = 7.0/(2.0*2.8);
        let frequency = (180.0_f64/2.8-decay*decay).sqrt();
        for pole in poles {
            close(pole.re, -decay, 2e-7);
            close(pole.im.abs(), frequency, 2e-7);
        }
        let expanded = pencil.expand_mode(&[C64::new(2.0, -0.5)], || false).unwrap();
        close(expanded[0].re, 2.0*n[0], 1e-14);
        close(expanded[1].im, -0.5*n[1], 1e-14);
        assert!(prepared.scaled_equilibrium_residual() <= 1e-10);
        assert_eq!(prepared.state(), &root.state);
    }
    assert_eq!(initial, before);
}

#[test]
fn fully_locked_and_zero_constraint_limits_do_not_create_massless_modes() {
    let (model, connection) = slider();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[10.0], external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        }, world_points: &[], connections: &[connection],
    };
    let rows = [row(vec![1.0], 0.2)];
    let locked = static_owner::equilibrate(&model, &seed(vec![0.0]), loads,
        EquilibriumBushings::Diagonal(&[]), &rows,
        controls(&[0.1], &[20.0]), joint_controls(&[20.0], &[0.1]), || false).unwrap();
    let prepared = response::prepare(&model, &locked, loads,
        EquilibriumBushings::Diagonal(&[]), &rows,
        response_controls(&[0.1], &[0.2], &[20.0]), subspace_controls(), || false).unwrap();
    assert_eq!(prepared.subspace().dimension(), 0);
    assert_eq!(prepared.lift_perturbation(&[], || false).unwrap(), vec![0.0]);
    assert!(prepared.project_effort(&[12.0], || false).unwrap().is_empty());
    assert!(prepared.inertia_action(&[], || false).unwrap().is_empty());
    let action = prepared.stiffness_action(&[], || false).unwrap();
    assert_eq!(action.original.residual_evaluations, 0);
    assert_eq!(action.original.effort, vec![0.0]);
    let pencil = prepared.quadratic_model(0, 0, || false).unwrap();
    assert_eq!(pencil.dimension(), 0);
    assert_eq!(pencil.source_evaluations(), 0);
    assert!(pencil.mass().is_empty());
    assert!(pencil.poles(1, || false).unwrap().is_empty());
    assert!(pencil.prepare_modes(fs_modal::quadratic_modes::PoleOptions {
        spectral_scale: 1.0, step_tolerance: 1e-10, residual_tolerance: 1e-10,
        max_iterations: 16, max_dimension: 8, max_dense_bytes: 100_000, max_work: 100_000,
    }, || false).unwrap().is_none());
    assert_eq!(pencil.expand_mode(&[], || false).unwrap(), vec![C64::ZERO]);

    let unlocked = static_owner::equilibrate(&model, &seed(vec![0.0]), loads,
        EquilibriumBushings::Diagonal(&[]), &[],
        controls(&[0.1], &[20.0]), joint_controls(&[], &[]), || false).unwrap();
    let prepared = response::prepare(&model, &unlocked, loads,
        EquilibriumBushings::Diagonal(&[]), &[],
        response_controls(&[0.1], &[0.2], &[20.0]), subspace_controls(), || false).unwrap();
    let pencil = prepared.quadratic_model(1, 64, || false).unwrap();
    close(prepared.subspace().basis()[0], 0.1, 0.0);
    close(pencil.mass()[0], 0.02, 1e-15);
    close(pencil.stiffness()[0], 1.0, 1e-8);
    close(pencil.velocity_operator()[0], 0.07, 1e-8);

    let fixed = ArticulatedModel::new(vec![Link::new("fixed", None, Se3::identity(),
        JointModel::FIXED, inertia(2.0, Vec3::default()))]).unwrap();
    let fixed_loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[], external_body_wrench: &[Wrench::default()],
            gravity_world: Vec3::new(0.0, -9.81, 0.0),
        }, world_points: &[], connections: &[],
    };
    let fixed_root = static_owner::equilibrate(&fixed, &seed(vec![]), fixed_loads,
        EquilibriumBushings::Diagonal(&[]), &[], controls(&[], &[]),
        joint_controls(&[], &[]), || false).unwrap();
    let fixed_response = response::prepare(&fixed, &fixed_root, fixed_loads,
        EquilibriumBushings::Diagonal(&[]), &[], response_controls(&[], &[], &[]),
        subspace_controls(), || false).unwrap();
    assert_eq!(fixed_response.subspace().physical_dimension(), 0);
    assert_eq!(fixed_response.quadratic_model(0, 0, || false).unwrap().source_evaluations(), 0);
}

#[test]
fn locked_branch_preserves_two_free_axes_gyroscopic_operator() {
    let model = ArticulatedModel::new(vec![
        Link::new("support", None, Se3::identity(), JointModel::FIXED,
            inertia(3.0, Vec3::default())),
        Link::new("x carriage", Some(0), Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(1.0, Vec3::default())),
        Link::new("y carriage", Some(1), Se3::identity(),
            JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None).unwrap(),
            inertia(2.0, Vec3::default())),
        Link::new("locked z", Some(0), Se3::identity(),
            JointModel::prismatic(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(4.0, Vec3::default())),
    ]).unwrap();
    let bushings = [SpatialBushing {
        a: MaterialFrame { link: 0, link_from_attachment: Se3::identity() },
        b: MaterialFrame { link: 2, link_from_attachment: Se3::identity() },
        rest_translation_a: Vec3::default(),
        translation_stiffness: Vec3::new(100.0, 120.0, 0.0),
        translation_damping: Vec3::new(4.0, 5.0, 0.0),
        rotation_stiffness: Vec3::default(), rotation_damping: Vec3::default(),
    }];
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0, 0.0, 2.0],
            external_body_wrench: &[Wrench::default(); 4], gravity_world: Vec3::default(),
        }, world_points: &[], connections: &[],
    };
    let rows = [row(vec![0.0, 0.0, 1.0], 0.0)];
    let mut initial = seed(vec![0.0; 3]);
    initial.base.twist_body.angular = Vec3::new(0.0, 0.0, 3.0);
    let root = static_owner::equilibrate_rotating(&model, &initial, loads,
        EquilibriumBushings::Diagonal(&bushings), &rows,
        controls(&[1.0; 3], &[100.0; 3]), joint_controls(&[10.0], &[1.0]), || false).unwrap();
    let prepared = response::prepare(&model, &root, loads,
        EquilibriumBushings::Diagonal(&bushings), &rows,
        response_controls(&[1.0; 3], &[1.0; 3], &[100.0; 3]),
        subspace_controls(), || false).unwrap();
    let pencil = prepared.quadratic_model(2, 128, || false).unwrap();
    assert_eq!(pencil.dimension(), 2);
    let n = prepared.subspace().basis();
    for i in 0..2 {
        close(n[4+i], 0.0, 1e-14);
        for j in 0..2 {
            let (xi, yi, xj, yj) = (n[i], n[2+i], n[j], n[2+j]);
            close(pencil.mass()[2*i+j], 3.0*xi*xj+2.0*yi*yj, 1e-12);
            close(pencil.stiffness()[2*i+j], 73.0*xi*xj+102.0*yi*yj, 2e-6);
            close(pencil.velocity_operator()[2*i+j],
                4.0*xi*xj-12.0*xi*yj+12.0*yi*xj+5.0*yi*yj, 2e-6);
        }
    }
    let c = pencil.velocity_operator();
    assert!((c[1]-c[2]).abs() > 23.9);
    let z = [0.7, -0.4];
    let v = prepared.lift_perturbation(&z, || false).unwrap();
    let power = z[0]*(c[0]*z[0]+c[1]*z[1])+z[1]*(c[2]*z[0]+c[3]*z[1]);
    close(power, 4.0*v[0]*v[0]+5.0*v[1]*v[1], 2e-6);
    assert_eq!(prepared.state(), &root.state);
}

#[test]
fn constrained_tangent_retains_an_inverted_pendulums_growth_rate() {
    let model = ArticulatedModel::new(vec![
        Link::new("support", None, Se3::identity(), JointModel::FIXED,
            inertia(3.0, Vec3::default())),
        Link::new("inverted pendulum", Some(0), Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(2.0, Vec3::new(0.0, 1.0, 0.0))),
        Link::new("locked slider", Some(0), Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(4.0, Vec3::default())),
    ]).unwrap();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0, 3.0], external_body_wrench: &[Wrench::default(); 3],
            gravity_world: Vec3::new(0.0, -9.81, 0.0),
        }, world_points: &[], connections: &[],
    };
    let rows = [row(vec![0.0, 1.0], 0.2)];
    let root = static_owner::equilibrate(&model, &seed(vec![0.0, 0.2]), loads,
        EquilibriumBushings::Diagonal(&[]), &rows, controls(&[1.0; 2], &[20.0; 2]),
        joint_controls(&[10.0], &[1.0]), || false).unwrap();
    let prepared = response::prepare(&model, &root, loads, EquilibriumBushings::Diagonal(&[]),
        &rows, response_controls(&[1.0; 2], &[1.0; 2], &[20.0; 2]),
        subspace_controls(), || false).unwrap();
    let pencil = prepared.quadratic_model(1, 64, || false).unwrap();
    close(pencil.mass()[0], 3.0, 1e-12);
    close(pencil.stiffness()[0], -19.62, 3e-6);
    close(pencil.velocity_operator()[0], 0.0, 1e-12);
    let poles = pencil.poles(128, || false).unwrap();
    let growth = (19.62_f64/3.0).sqrt();
    assert_eq!(poles.len(), 2);
    for pole in &poles {
        close(pole.im, 0.0, 1e-12);
        close(pole.re.abs(), growth, 3e-7);
    }
    assert!(poles.iter().any(|p| p.re > 0.0));
    assert!(poles.iter().any(|p| p.re < 0.0));
}

#[test]
fn forged_public_equilibrium_fields_and_original_manifold_violations_refuse() {
    let (model, connections) = two_sliders();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[10.0, 40.0],
            external_body_wrench: &[Wrench::default(); 3], gravity_world: Vec3::default(),
        }, world_points: &[], connections: &connections,
    };
    let rows = [row(vec![-2.0, 1.0], 0.1)];
    let root = static_owner::equilibrate(&model, &seed(vec![0.0; 2]), loads,
        EquilibriumBushings::Diagonal(&[]), &rows, controls(&[0.1; 2], &[40.0; 2]),
        joint_controls(&[20.0], &[0.1]), || false).unwrap();
    let original = root.clone();
    let mut forged = root.clone();
    forged.multipliers[0] += 1.0;
    forged.scaled_force_residual_norm = 0.0;
    forged.joint_effort_residual.fill(0.0);
    forged.required_effort.fill(0.0);
    forged.applied_effort.fill(0.0);
    forged.reaction_effort.fill(0.0);
    assert!(response::prepare(&model, &forged, loads, EquilibriumBushings::Diagonal(&[]),
        &rows, response_controls(&[1.0; 2], &[1.0; 2], &[40.0; 2]),
        subspace_controls(), || false).is_err());
    let mut off_manifold = root.clone();
    off_manifold.state.joint_position[1] += 0.01;
    off_manifold.position_residuals.fill(0.0);
    assert!(response::prepare(&model, &off_manifold, loads, EquilibriumBushings::Diagonal(&[]),
        &rows, response_controls(&[1.0; 2], &[1.0; 2], &[40.0; 2]),
        subspace_controls(), || false).is_err());
    assert!(response::prepare(&model, &root, loads, EquilibriumBushings::Diagonal(&[]),
        &rows, response_controls(&[1.0; 2], &[1.0; 2], &[40.0; 2]),
        subspace_controls(), || true).is_err());
    let prepared = response::prepare(&model, &root, loads, EquilibriumBushings::Diagonal(&[]),
        &rows, response_controls(&[1.0; 2], &[1.0; 2], &[40.0; 2]),
        subspace_controls(), || false).unwrap();
    assert!(prepared.quadratic_model(1, 8, || false).is_err());
    let mut checkpoints = 0;
    let expected = prepared.quadratic_model(1, 64, || {
        checkpoints += 1;
        false
    }).unwrap();
    let mut calls = 0;
    assert!(prepared.quadratic_model(1, 64, || {
        calls += 1;
        calls == checkpoints
    }).is_err());
    let retry = prepared.quadratic_model(1, 64, || false).unwrap();
    assert_eq!(retry.mass(), expected.mass());
    assert_eq!(retry.stiffness(), expected.stiffness());
    assert_eq!(retry.velocity_operator(), expected.velocity_operator());
    assert_eq!(root, original);
}

#[path = "constrained_linearized/forced.rs"]
mod forced;
