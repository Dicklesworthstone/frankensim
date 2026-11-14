//! G1 original M/C/K about stationary and rotating equilibria.
//!
//! Analytic mechanics below is an independent oracle. Every action, matrix
//! and pole is obtained from the original articulated/material/modal owners.
use super::*;
use fs_mbd::articulated::motion::spatial_loads::prescribed::equilibrium::linearized::{
    CentralDifferenceControls, LinearizationControls, LinearizationError, prepare,
};
use fs_solver::finite_difference::central::CentralDifferenceError;

fn response_controls<'a>(
    coordinates: &'a [f64],
    velocities: &'a [f64],
    efforts: &'a [f64],
) -> LinearizationControls<'a> {
    LinearizationControls {
        coordinate_scales: coordinates,
        velocity_scales: velocities,
        effort_scales: efforts,
        equilibrium_tolerance: 1e-10,
        maximum_stencil_distortion: 1e-8,
        difference: CentralDifferenceControls {
            initial_step: 1e-3,
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-7,
            maximum_refinements: 8,
            maximum_evaluations: 18,
        },
    }
}

fn close_rows(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (&actual, &expected) in actual.iter().zip(expected) {
        close(actual, expected, tolerance);
    }
}

#[test]
fn native_slider_actions_and_poles_retain_damping_and_centrifugal_softening() {
    let (model, connection) = slider();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[connection],
    };
    for omega in [0.0_f64, 3.0] {
        // Physical radius is 1+q. k*q = m*Omega^2*(1+q).
        let stiffness = 100.0 - 2.0 * omega * omega;
        let mut state = seed(vec![2.0 * omega * omega / stiffness]);
        state.base.twist_body.angular = Vec3::new(0.0, 0.0, omega);
        let before = state.clone();
        let response = prepare(
            &model, &state, loads, EquilibriumBushings::Diagonal(&[]),
            response_controls(&[0.1], &[0.2], &[100.0]), || false,
        ).unwrap();
        close_rows(&response.inertia_action(&[2.5], || false).unwrap(), &[5.0], 1e-12);
        let k = response.stiffness_action(&[0.3], || false).unwrap();
        let c = response.velocity_action(&[-0.7], || false).unwrap();
        close_rows(&k.effort, &[0.3 * stiffness], 2e-7);
        close_rows(&c.effort, &[-4.9], 2e-7);
        assert!(k.residual_evaluations >= 4 && c.residual_evaluations >= 4);
        assert!(k.maximum_stencil_distortion <= 1e-8);
        assert!(c.maximum_stencil_distortion <= 1e-8);
        let pencil = response.quadratic_model(1, 64, || false).unwrap();
        close_rows(pencil.mass(), &[2.0], 1e-12);
        close_rows(pencil.stiffness(), &[stiffness], 2e-7);
        close_rows(pencil.velocity_operator(), &[7.0], 2e-7);
        assert_eq!(pencil.dimension(), 1);
        let poles = pencil.poles(128, || false).unwrap();
        assert_eq!(poles.len(), 2);
        let damped_frequency = (stiffness / 2.0 - 1.75_f64.powi(2)).sqrt();
        for pole in &poles {
            close(pole.re, -1.75, 1e-7);
            close(pole.im.abs(), damped_frequency, 1e-7);
        }
        assert!(poles[0].im * poles[1].im < 0.0);
        assert_eq!(response.state(), &before);
        assert_eq!(state, before);
    }
}

#[test]
fn transverse_preload_supplies_geometric_stiffness_without_false_axial_damping() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "support", None, Se3::identity(), JointModel::FIXED,
            inertia(3.0, Vec3::default()),
        ),
        Link::new(
            "transverse slider", Some(0),
            Se3::from_parts(So3::identity(), Vec3::new(0.0, 1.0, 0.0)).unwrap(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(2.0, Vec3::default()),
        ),
    ]).unwrap();
    let connection = AxialConnection {
        a: MaterialPoint { link: 0, point_body: Vec3::default() },
        b: MaterialPoint { link: 1, point_body: Vec3::default() },
        rest_length_m: 0.5,
        stiffness_n_per_m: 100.0,
        damping_ns_per_m: 7.0,
    };
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[connection],
    };
    let state = seed(vec![0.0]);
    let response = prepare(
        &model, &state, loads, EquilibriumBushings::Diagonal(&[]),
        response_controls(&[0.1], &[0.2], &[100.0]), || false,
    ).unwrap();
    // At length 1, tension is 50 N. The transverse stiffness is tension/length,
    // not the 100 N/m axial material coefficient. Transverse axial strain rate
    // vanishes to first order even though the original dashpot is nonzero.
    let k = response.stiffness_action(&[1.0], || false).unwrap();
    close_rows(&k.effort, &[50.0], 5e-7);
    close_rows(&response.velocity_action(&[1.0], || false).unwrap().effort, &[0.0], 1e-12);
    let pencil = response.quadratic_model(1, 64, || false).unwrap();
    close_rows(pencil.mass(), &[2.0], 1e-12);
    close_rows(pencil.stiffness(), &[50.0], 5e-7);
    close_rows(pencil.velocity_operator(), &[0.0], 1e-12);
    let poles = pencil.poles(128, || false).unwrap();
    assert_eq!(poles.len(), 2);
    for pole in &poles {
        close(pole.re, 0.0, 1e-12);
        close(pole.im.abs(), 5.0, 1e-6);
    }
    assert!(poles[0].im * poles[1].im < 0.0);
}

#[test]
fn inverted_pendulum_reports_its_negative_stiffness_and_growing_mode() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "support", None, Se3::identity(), JointModel::FIXED,
            inertia(3.0, Vec3::default()),
        ),
        Link::new(
            "inverted pendulum", Some(0), Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(2.0, Vec3::new(0.0, 1.0, 0.0)),
        ),
    ]).unwrap();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::new(0.0, -9.81, 0.0),
        },
        world_points: &[],
        connections: &[],
    };
    let state = seed(vec![0.0]);
    let response = prepare(
        &model, &state, loads, EquilibriumBushings::Diagonal(&[]),
        response_controls(&[0.1], &[0.2], &[20.0]), || false,
    ).unwrap();
    // I_hinge = I_COM + m*l^2 = 3; d^2(m*g*l*cos(q))/dq^2 = -19.62.
    close_rows(&response.inertia_action(&[1.0], || false).unwrap(), &[3.0], 1e-12);
    let pencil = response.quadratic_model(1, 64, || false).unwrap();
    close_rows(pencil.mass(), &[3.0], 1e-12);
    close_rows(pencil.stiffness(), &[-19.62], 5e-7);
    close_rows(pencil.velocity_operator(), &[0.0], 1e-12);
    let poles = pencil.poles(128, || false).unwrap();
    assert_eq!(poles.len(), 2);
    let growth = (19.62_f64 / 3.0).sqrt();
    for pole in &poles {
        close(pole.im, 0.0, 1e-12);
        close(pole.re.abs(), growth, 1e-7);
    }
    assert!(poles.iter().any(|pole| pole.re > 0.0));
    assert!(poles.iter().any(|pole| pole.re < 0.0));
    assert_eq!(response.state(), &state);
}

#[test]
fn two_axis_rotating_tree_keeps_skew_gyroscopic_and_material_velocity_terms() {
    let model = ArticulatedModel::new(vec![
        Link::new(
            "support", None, Se3::identity(), JointModel::FIXED,
            inertia(3.0, Vec3::default()),
        ),
        Link::new(
            "x carriage", Some(0), Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(1.0, Vec3::default()),
        ),
        Link::new(
            "y carriage", Some(1), Se3::identity(),
            JointModel::prismatic(Vec3::new(0.0, 1.0, 0.0), None).unwrap(),
            inertia(2.0, Vec3::default()),
        ),
    ]).unwrap();
    let bushings = [SpatialBushing {
        a: MaterialFrame { link: 0, link_from_attachment: Se3::identity() },
        b: MaterialFrame { link: 2, link_from_attachment: Se3::identity() },
        rest_translation_a: Vec3::default(),
        translation_stiffness: Vec3::new(100.0, 120.0, 0.0),
        translation_damping: Vec3::new(4.0, 5.0, 0.0),
        rotation_stiffness: Vec3::default(),
        rotation_damping: Vec3::default(),
    }];
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0; 2],
            external_body_wrench: &[Wrench::default(); 3],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[],
    };
    let mut state = seed(vec![0.0, 0.0]);
    state.base.twist_body.angular = Vec3::new(0.0, 0.0, 3.0);
    let before = state.clone();
    let response = prepare(
        &model, &state, loads, EquilibriumBushings::Diagonal(&bushings),
        response_controls(&[0.1, 0.2], &[0.3, 0.4], &[100.0, 120.0]), || false,
    ).unwrap();
    // x accelerates both carriages; y accelerates only the 2 kg tip.
    // Centrifugal K removes diag(3,2)*Omega^2. The tip supplies
    // C_xy=-2*m_tip*Omega and C_yx=+2*m_tip*Omega; these must not be symmetrized.
    close_rows(&response.inertia_action(&[2.0, -3.0], || false).unwrap(), &[6.0, -6.0], 1e-12);
    close_rows(&response.stiffness_action(&[0.2, -0.3], || false).unwrap().effort, &[14.6, -30.6], 2e-7);
    close_rows(&response.velocity_action(&[0.2, -0.3], || false).unwrap().effort, &[4.4, 0.9], 2e-7);
    let pencil = response.quadratic_model(2, 128, || false).unwrap();
    close_rows(pencil.mass(), &[3.0, 0.0, 0.0, 2.0], 1e-12);
    close_rows(pencil.stiffness(), &[73.0, 0.0, 0.0, 102.0], 2e-7);
    close_rows(pencil.velocity_operator(), &[4.0, -12.0, 12.0, 5.0], 2e-7);
    let c = pencil.velocity_operator();
    let (vx, vy) = (0.7, -0.4);
    let power = vx * (c[0] * vx + c[1] * vy) + vy * (c[2] * vx + c[3] * vy);
    close(power, 4.0 * vx * vx + 5.0 * vy * vy, 2e-7);
    assert_eq!(response.state(), &before);
    assert_eq!(state, before);
}

#[test]
fn work_limits_final_cancellation_and_unbalanced_inputs_publish_no_partial_response() {
    let (model, connection) = slider();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[connection],
    };
    let state = seed(vec![0.0]);
    let before = state.clone();
    let response = prepare(
        &model, &state, loads, EquilibriumBushings::Diagonal(&[]),
        response_controls(&[0.1], &[0.2], &[100.0]), || false,
    ).unwrap();
    // One native mass extraction plus two central pairs for K and C = 9.
    assert!(matches!(
        response.quadratic_model(1, 8, || false),
        Err(LinearizationError::Budget { used: 5, maximum: 8 })
    ));
    let mut checkpoints = 0;
    let accepted = response.quadratic_model(1, 9, || {
        checkpoints += 1;
        false
    }).unwrap();
    assert_eq!(accepted.source_evaluations(), 9);
    assert!(checkpoints > 9);
    let mut attempted = 0;
    assert!(matches!(
        response.quadratic_model(1, 9, || {
            attempted += 1;
            attempted == checkpoints
        }),
        Err(LinearizationError::Equilibrium(EquilibriumError::Motion(
            PrescribedBaseMotionError::Cancelled
        )))
    ));
    assert_eq!(attempted, checkpoints);
    assert!(matches!(
        response.stiffness_action(&[f64::NAN], || false),
        Err(LinearizationError::Invalid(_))
    ));
    let unbalanced = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[1.0],
            ..loads.held
        },
        ..loads
    };
    assert!(matches!(
        prepare(
            &model, &state, unbalanced, EquilibriumBushings::Diagonal(&[]),
            response_controls(&[0.1], &[0.2], &[100.0]), || false,
        ),
        Err(LinearizationError::Unbalanced { .. })
    ));
    // Refusals and cancellation do not modify the prepared root or its seed.
    assert_eq!(response.state(), &before);
    assert_eq!(state, before);
    let retried = response.quadratic_model(1, 9, || false).unwrap();
    assert_eq!(retried.mass(), accepted.mass());
    assert_eq!(retried.stiffness(), accepted.stiffness());
    assert_eq!(retried.velocity_operator(), accepted.velocity_operator());
}

#[test]
fn quantized_physical_stencils_and_erased_scaled_directions_are_refused() {
    let (model, _) = slider();
    let loads = SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.0],
            external_body_wrench: &[Wrench::default(); 2],
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[],
    };
    // The slider sits where its locked inertia is still representable. At
    // q = 1e15 the composite inertia about the system centre holds 2 kg m^2
    // against 1.2e30, below what f64 can certify, and prepare refuses it as
    // NonPhysicalInertia before any stencil is formed.
    let state = seed(vec![1.5e6]);
    let response = prepare(
        &model, &state, loads, EquilibriumBushings::Diagonal(&[]),
        response_controls(&[2.0], &[1.0], &[1.0]), || false,
    ).unwrap();
    // Requested displacement 2*0.001 = 0.002 is represented at q = 1.5e6, whose
    // spacing is 2.3e-10, as 0.002 + 9.5e-11: a relative distortion of 4.75e-8,
    // beyond the relative allowance of 1e-8.
    // Both stencils could still return identical zero forces for this neutral
    // slider. That agreement cannot admit a materially distorted physical chart.
    assert!(matches!(
        response.stiffness_action(&[2.0], || false),
        Err(LinearizationError::Difference(CentralDifferenceError::Source(
            EquilibriumError::Invalid("physical tangent stencil exceeds representation allowance")
        )))
    ));
    // This positive physical direction vanishes when divided by its scale.
    assert!(matches!(
        response.stiffness_action(&[f64::from_bits(1)], || false),
        Err(LinearizationError::Invalid("scaled tangent direction is not representable"))
    ));
    let zero = response.stiffness_action(&[0.0], || false).unwrap();
    assert_eq!(zero.effort, [0.0]);
    assert_eq!(zero.residual_evaluations, 0);
    assert_eq!(response.state(), &state);
}
