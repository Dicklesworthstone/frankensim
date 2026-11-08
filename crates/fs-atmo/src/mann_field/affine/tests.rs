use super::*;
use crate::mann::MannParams;
use crate::mann_field::fixed_step::{FixedStepMannField, MannOuBridge};
use crate::mann_field::{HomogeneousMannField, MannQuadrature, SpectralMode};

const MODES: usize = 16;
const BYTES: usize = 88 * MODES;
const LOOSE: MannAffineAccuracy = MannAffineAccuracy {
    maximum_velocity_error_mps: 1e6,
    maximum_gradient_error_per_s: 1e6,
};

fn source(seed: u64) -> HomogeneousMannField {
    let model = MannSpectralModel::build(
        31,
        MannParams {
            alpha_eps23: 0.2,
            length_m: 3.0,
            gamma: 3.9,
        },
        MannQuadrature {
            n_modes: MODES,
            min_kl: 0.03,
            max_kl: 16.0,
        },
        [2.0, -0.4, 0.3],
    )
    .unwrap();
    HomogeneousMannField::at_anchor(model, seed, vec![0.3; MODES]).unwrap()
}

fn check_errors(exact: FlowSample, approximate: MannAffineSample) {
    let velocity = norm((0..3).map(|i| exact.u[i] - approximate.flow.u[i])).unwrap();
    let gradient = norm(
        (0..3).flat_map(|i| (0..3).map(move |j| exact.grad[i][j] - approximate.flow.grad[i][j])),
    )
    .unwrap();
    let roundoff =
        2e-13 * (1.0 + norm(exact.u).unwrap() + norm(exact.grad.into_iter().flatten()).unwrap());
    assert!(
        velocity <= approximate.error.velocity_mps + roundoff,
        "velocity {velocity} exceeds {:?}",
        approximate.error
    );
    assert!(
        gradient <= approximate.error.gradient_per_s + roundoff,
        "gradient {gradient} exceeds {:?}",
        approximate.error
    );
}

#[test]
fn homogeneous_realizations_bound_both_errors_and_retain_the_original_anchor() {
    let anchor = [0.7, -0.4, 1.3];
    for seed in [7, 19, 81] {
        let mut field = source(seed);
        for step in [0, 3, 11] {
            field
                .advance_to(crate::ou::STATIONARY_ANCHOR_TICK + step)
                .unwrap();
            let before = field.checkpoint();
            let observation = field.affine_observation(anchor, None, || false).unwrap();
            assert_eq!(
                observation.sample(anchor, LOOSE).unwrap().flow,
                field.sample(anchor).unwrap()
            );
            let fixed = FixedStepMannField::from_snapshot(field.clone(), 43, 0.05, BYTES).unwrap();
            assert_eq!(
                observation,
                fixed.affine_observation(anchor, None, || false).unwrap()
            );
            let trace = observation.anchor.grad[0][0]
                + observation.anchor.grad[1][1]
                + observation.anchor.grad[2][2];
            assert!(trace.abs() < 2e-13);
            for radius in [1e-4, 0.01, 0.1, 0.7] {
                for direction in [[1.0, 0.2, -0.3], [-0.4, 1.0, 0.7], [0.1, -0.8, 1.0]] {
                    let point = core::array::from_fn(|i| anchor[i] + radius * direction[i]);
                    let approximate = observation.sample(point, LOOSE).unwrap();
                    check_errors(field.sample(point).unwrap(), approximate);
                    assert_eq!(approximate.flow.grad, observation.anchor.grad);
                    assert!(approximate.error.velocity_mps > 0.0);
                    assert!(approximate.error.gradient_per_s > 0.0);
                }
            }
            assert_eq!(field.checkpoint(), before);
        }
    }
}

#[test]
fn wall_remainder_includes_exponential_derivatives_and_preserves_slip() {
    let field = source(19);
    let wall = MannPlaneWall::new(-0.7).unwrap();
    let anchor = [0.7, -0.4, wall.height_m()];
    let observation = field
        .affine_observation(anchor, Some(wall), || false)
        .unwrap();
    let homogeneous = field.affine_observation(anchor, None, || false).unwrap();
    assert!(observation.velocity_curvature > homogeneous.velocity_curvature);
    assert!(observation.gradient_curvature > observation.velocity_curvature);
    for height in [0.0, 1e-4, 0.1, 0.7] {
        for shift in [[0.0, 0.0], [0.3, -0.4], [-0.7, 0.2]] {
            let point = [
                anchor[0] + shift[0],
                anchor[1] + shift[1],
                anchor[2] + height,
            ];
            let approximate = observation.sample(point, LOOSE).unwrap();
            check_errors(field.sample_wall_blocked(point, wall).unwrap(), approximate);
            if height == 0.0 {
                assert_eq!(approximate.flow.u[2], 0.0);
                assert_eq!(approximate.flow.grad[2][0], 0.0);
                assert_eq!(approximate.flow.grad[2][1], 0.0);
            }
            let g = approximate.flow.grad;
            assert!((g[0][0] + g[1][1] + g[2][2]).abs() < 2e-13);
        }
    }

    // Manufactured solenoidal Fourier mode w=cos(x), with the ORIGINAL
    // harmonic correction (-sin(x)e^-z, 0, -cos(x)e^-z). At x=0 the
    // gradient error approaches sqrt(2) dz, exceeding the homogeneous-only
    // Lipschitz constant 1. This detects an omitted wall-envelope derivative.
    let mut model = field.model().clone();
    model.quadrature.n_modes = 1;
    model.modes = vec![SpectralMode {
        k: [1.0, 0.0, 0.0],
        weight: 1.0,
        polarization: [[0.0, 0.0, 1.0], [0.0, 1.0, 0.0]],
    }];
    let amplitudes: [&[f64]; 4] = [&[1.0], &[0.0], &[0.0], &[0.0]];
    let wall = MannPlaneWall::new(0.0).unwrap();
    let affine = prepare_modes(&model, 0.0, amplitudes, [0.0; 3], Some(wall), || false).unwrap();
    assert_eq!(affine.velocity_curvature, 2.0);
    assert_eq!(affine.gradient_curvature, 1.0 + core::f64::consts::SQRT_2);
    let point = [0.0, 0.0, 1e-4];
    let exact =
        crate::mann_field::wall::sample_modes(&model, 0.0, amplitudes, point, wall).unwrap();
    let approximate = affine.sample(point, LOOSE).unwrap();
    let gradient_error = norm(
        (0..3).flat_map(|i| (0..3).map(move |j| exact.grad[i][j] - approximate.flow.grad[i][j])),
    )
    .unwrap();
    assert!(gradient_error > 1.4e-4);
    check_errors(exact, approximate);
}

#[test]
fn both_temporal_owners_and_reordered_bridge_stages_use_their_actual_realization() {
    let mut field = FixedStepMannField::from_snapshot(source(11), 43, 0.05, BYTES).unwrap();
    let before = field.checkpoint();
    let anchor = [0.0, 0.0, 1.0];
    let old = field.affine_observation(anchor, None, || false).unwrap();
    let bridge = MannOuBridge::prepare(&field, 18, 128, BYTES, 4 * MODES, || false).unwrap();
    let mut frame = bridge.frame(BYTES, || false).unwrap();
    let at = |fraction: f64| {
        bridge.start_time_s() + fraction * (bridge.end_time_s() - bridge.start_time_s())
    };
    let mut retained = None;
    for fraction in [0.0, 0.75, 0.25, 0.5, 0.75, 1.0] {
        bridge
            .sample_frame(at(fraction), &mut frame, 100_000, || false)
            .unwrap();
        let observation = frame.affine_observation(anchor, None, || false).unwrap();
        assert_eq!(observation.time_s(), at(fraction));
        assert_eq!(
            observation.sample(anchor, LOOSE).unwrap().flow,
            frame.sample(anchor).unwrap()
        );
        check_errors(
            frame.sample([0.03, -0.02, 1.01]).unwrap(),
            observation.sample([0.03, -0.02, 1.01], LOOSE).unwrap(),
        );
        if fraction == 0.0 {
            assert_eq!(observation, old);
        }
        if fraction == 0.75 {
            if let Some(previous) = retained {
                assert_eq!(observation, previous);
            }
            retained = Some(observation);
        }
    }
    assert_eq!(field.checkpoint(), before);
    field.advance_to(1, 4 * MODES, || false).unwrap();
    let next = field.affine_observation(anchor, None, || false).unwrap();
    assert_eq!(
        next,
        frame.affine_observation(anchor, None, || false).unwrap()
    );
    assert_ne!(next.anchor, old.anchor);
    assert_ne!(next.velocity_curvature, old.velocity_curvature);
}

#[test]
fn spatial_budget_and_unrepresentable_nonzero_remainders_refuse_explicitly() {
    let field = source(7);
    let observation = field.affine_observation([0.0; 3], None, || false).unwrap();
    let zero = MannAffineAccuracy {
        maximum_velocity_error_mps: 0.0,
        maximum_gradient_error_per_s: 0.0,
    };
    assert_eq!(
        observation.sample([0.0; 3], zero).unwrap().error,
        MannAffineError::default()
    );
    let point = [0.1, 0.2, -0.3];
    let bound = observation.sample(point, LOOSE).unwrap().error;
    for accuracy in [
        MannAffineAccuracy {
            maximum_velocity_error_mps: 0.5 * bound.velocity_mps,
            ..LOOSE
        },
        MannAffineAccuracy {
            maximum_gradient_error_per_s: 0.5 * bound.gradient_per_s,
            ..LOOSE
        },
        MannAffineAccuracy {
            maximum_gradient_error_per_s: f64::NAN,
            ..LOOSE
        },
    ] {
        assert_eq!(
            observation.sample(point, accuracy).unwrap_err().code,
            "mann-affine-accuracy"
        );
    }
    assert_eq!(
        observation
            .sample([1e-180, 0.0, 0.0], LOOSE)
            .unwrap_err()
            .code,
        "mann-affine-numerical"
    );
    assert_eq!(
        observation
            .sample([f64::MAX, f64::MAX, 0.0], LOOSE)
            .unwrap_err()
            .code,
        "mann-affine-numerical"
    );
    let wall = MannPlaneWall::new(0.0).unwrap();
    assert_eq!(
        field
            .affine_observation([0.0, 0.0, 0.1], Some(wall), || false)
            .unwrap_err()
            .code,
        "mann-affine-anchor"
    );
    let blocked = field
        .affine_observation([0.0; 3], Some(wall), || false)
        .unwrap();
    assert_eq!(
        blocked.sample([0.0, 0.0, -0.1], LOOSE).unwrap_err().code,
        "mann-affine-domain"
    );
}

#[test]
fn cancellation_preserves_source_and_zero_realizations_do_not_invent_error() {
    let field = FixedStepMannField::from_snapshot(source(11), 43, 0.05, BYTES).unwrap();
    let before = field.checkpoint();
    for stop in [1, MODES / 2, MODES + 2] {
        let mut polls = 0;
        assert_eq!(
            field
                .affine_observation([0.0; 3], None, || {
                    polls += 1;
                    polls == stop
                })
                .unwrap_err()
                .code,
            "mann-affine-cancelled"
        );
        assert_eq!(field.checkpoint(), before);
    }
    let mut checkpoint = before;
    for amplitudes in &mut checkpoint.amplitudes {
        amplitudes.fill(0.0);
    }
    let zero = FixedStepMannField::resume(checkpoint, BYTES).unwrap();
    let observation = zero.affine_observation([0.0; 3], None, || false).unwrap();
    let sample = observation
        .sample(
            [10.0, -20.0, 30.0],
            MannAffineAccuracy {
                maximum_velocity_error_mps: 0.0,
                maximum_gradient_error_per_s: 0.0,
            },
        )
        .unwrap();
    assert_eq!(sample.flow, zero.sample([10.0, -20.0, 30.0]).unwrap());
    assert_eq!(sample.error, MannAffineError::default());
}
