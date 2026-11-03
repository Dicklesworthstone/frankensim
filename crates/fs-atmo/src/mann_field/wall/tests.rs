use super::*;
use crate::mann::MannParams;

fn model() -> MannSpectralModel {
    MannSpectralModel::build(
        41,
        MannParams {
            alpha_eps23: 0.32,
            length_m: 1.0,
            gamma: 3.9,
        },
        MannQuadrature {
            n_modes: 32,
            min_kl: 0.03,
            max_kl: 12.0,
        },
        [1.7, -0.2, 0.4],
    )
    .unwrap()
}

fn field() -> HomogeneousMannField {
    HomogeneousMannField::at_anchor(model(), 79, vec![0.7; 32]).unwrap()
}

// Independent tensor oracle: multiply the complex harmonic transfer matrices
// by the original Phi, without using polarization factors or ModePoint.
fn tensor_covariance(
    model: &MannSpectralModel,
    first: [f64; 3],
    second: [f64; 3],
    lag: f64,
    wall: MannPlaneWall,
    tau: f64,
) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for mode in &model.modes {
        let q = mode.k[0].hypot(mode.k[1]);
        let transfer = |point: [f64; 3], time: f64| {
            let position = core::array::from_fn(|i| point[i] - model.advection_mps[i] * time);
            let phase = dot(mode.k, position);
            let trace = dot(
                mode.k,
                [
                    position[0],
                    position[1],
                    wall.height_m() - model.advection_mps[2] * time,
                ],
            );
            let decay = (-q * (point[2] - wall.height_m())).exp();
            let mut real = [[0.0; 3]; 3];
            let mut imaginary = [[0.0; 3]; 3];
            for i in 0..3 {
                real[i][i] = phase.cos();
                imaginary[i][i] = phase.sin();
                if i < 2 {
                    real[i][2] -= mode.k[i] / q * decay * trace.sin();
                    imaginary[i][2] += mode.k[i] / q * decay * trace.cos();
                } else {
                    real[i][2] -= decay * trace.cos();
                    imaginary[i][2] -= decay * trace.sin();
                }
            }
            (real, imaginary)
        };
        let (ar, ai) = transfer(first, 0.0);
        let (br, bi) = transfer(second, lag);
        let phi = mann_tensor(&model.parameters, mode.k).unwrap();
        let weight = mode.weight * (-lag.abs() / tau).exp();
        for i in 0..3 {
            for j in 0..3 {
                for a in 0..3 {
                    for b in 0..3 {
                        result[i][j] +=
                            weight * phi[a][b] * (ar[i][a] * br[j][b] + ai[i][a] * bi[j][b]);
                    }
                }
            }
        }
    }
    result
}

#[test]
fn wall_normal_velocity_and_tangent_derivatives_vanish_without_losing_divergence() {
    let field = field();
    let wall = MannPlaneWall::new(0.0).unwrap();
    let point = [0.2, -0.3, 0.0];
    let boundary = field.sample_wall_blocked(point, wall).unwrap();
    assert_eq!(boundary.u[2], 0.0);
    assert_eq!(boundary.grad[2][0], 0.0);
    assert_eq!(boundary.grad[2][1], 0.0);
    let small = 1e-12;
    let near = field
        .sample_wall_blocked([point[0], point[1], small], wall)
        .unwrap();
    assert!(
        (near.u[2] / small - boundary.grad[2][2]).abs()
            < 1e-8 * boundary.grad[2][2].abs().max(1.0)
    );
    for height in [0.0, small, 0.13, 2.3] {
        let sample = field
            .sample_wall_blocked([point[0], point[1], height], wall)
            .unwrap();
        let trace: f64 = (0..3).map(|i| sample.grad[i][i]).sum();
        let scale = sample.grad.iter().flatten().map(|v| v.abs()).sum::<f64>();
        assert!(trace.abs() < 2e-12 * scale.max(1.0));
    }
}

#[test]
fn all_nine_wall_gradients_match_independent_spatial_differences() {
    let field = field();
    let wall = MannPlaneWall::new(-0.4).unwrap();
    for point in [[0.2, -0.3, -0.27], [0.1, 0.7, 1.9]] {
        let sample = field.sample_wall_blocked(point, wall).unwrap();
        let step = 1e-5;
        for j in 0..3 {
            let mut plus = point;
            let mut minus = point;
            plus[j] += step;
            minus[j] -= step;
            let plus = field.sample_wall_blocked(plus, wall).unwrap();
            let minus = field.sample_wall_blocked(minus, wall).unwrap();
            for i in 0..3 {
                let difference = (plus.u[i] - minus.u[i]) / (2.0 * step);
                assert!(
                    (difference - sample.grad[i][j]).abs()
                        < 3e-6 * sample.grad[i][j].abs().max(1.0),
                    "[{i},{j}] {difference} != {}",
                    sample.grad[i][j]
                );
            }
        }
    }
}

#[test]
fn covariance_transforms_the_complete_original_tensor_at_both_heights_and_times() {
    let field = field();
    let wall = MannPlaneWall::new(-0.4).unwrap();
    let first = [0.2, -0.3, -0.27];
    let second = [0.6, 0.1, 0.8];
    let variance = field.model().covariance([0.0; 3]).unwrap();
    assert!(variance[0][2].abs() > 0.01 * variance[0][0]);
    for lag in [-0.13, 0.0, 0.21] {
        let actual = field.covariance_wall_blocked(first, second, lag, wall).unwrap();
        let expected = tensor_covariance(field.model(), first, second, lag, wall, 0.7);
        let reverse = field.covariance_wall_blocked(second, first, -lag, wall).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let scale = (variance[i][i] * variance[j][j]).sqrt();
                assert!((actual[i][j] - expected[i][j]).abs() < 2e-10 * scale);
                assert!((actual[i][j] - reverse[j][i]).abs() < 2e-12 * scale);
            }
        }
        if lag == 0.0 {
            assert_eq!(
                actual,
                field.model().covariance_wall_blocked(first, second, wall).unwrap()
            );
        }
    }
    let on_wall = [0.2, -0.3, wall.height_m()];
    let cross = field.covariance_wall_blocked(on_wall, second, 0.21, wall).unwrap();
    assert_eq!(cross[2], [0.0; 3]);
    let cross = field.covariance_wall_blocked(first, on_wall, -0.13, wall).unwrap();
    assert!((0..3).all(|i| cross[i][2] == 0.0));
}

#[test]
fn far_field_recovers_the_original_realization_and_full_covariance() {
    let field = field();
    let min_q = field.model().modes.iter()
        .map(|mode| mode.k[0].hypot(mode.k[1]))
        .fold(f64::INFINITY, f64::min);
    let first = [0.2, -0.3, 0.1];
    let second = [0.6, 0.1, 0.8];
    for wall_height in [-1000.0 / min_q, -1e308] {
        let wall = MannPlaneWall::new(wall_height).unwrap();
        let original = field.sample(first).unwrap();
        let blocked = field.sample_wall_blocked(first, wall).unwrap();
        for (a, b) in original
            .u
            .iter()
            .chain(original.grad.iter().flatten())
            .zip(blocked.u.iter().chain(blocked.grad.iter().flatten()))
        {
            assert!((a - b).abs() < 2e-12 * a.abs().max(1.0));
        }
        let separation = core::array::from_fn(|i| second[i] - first[i]);
        let expected = field.covariance(separation, 0.21).unwrap();
        let actual = field.covariance_wall_blocked(first, second, 0.21, wall).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                assert!((actual[i][j] - expected[i][j]).abs() < 2e-12);
            }
        }
    }
}

#[test]
fn vertical_and_tiny_horizontal_modes_need_no_epsilon_wavenumber() {
    let mut model = model();
    model.modes = vec![SpectralMode {
        k: [0.0, 0.0, 2.0],
        weight: 1.0,
        polarization: [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
    }];
    model.advection_mps = [0.0; 3];
    let wall = MannPlaneWall::new(0.0).unwrap();
    let amplitudes = [&[0.3][..], &[-0.4][..], &[0.7][..], &[0.8][..]];
    let point = [0.1, 0.2, 0.4];
    assert_eq!(
        sample_modes(&model, 0.0, amplitudes, point, wall).unwrap(),
        super::super::sample_modes(&model, 0.0, amplitudes, point, None).unwrap()
    );
    let mode = SpectralMode {
        k: [1e-300, 0.0, 1.0],
        weight: 1.0,
        polarization: [[-1.0, 0.0, 1e-300], [0.0, 1.0, 0.0]],
    };
    let geometry = ModePoint::new(&mode, [0.0; 3], 0.0, [0.0, 0.0, 1e-10], wall, 1e-10)
        .unwrap();
    assert_eq!(geometry.q, 1e-300);
    assert_eq!(geometry.horizontal, [1.0, 0.0]);
    let value = geometry.values([0.0; 3], mode.polarization[0]);
    assert!((value[2] / 1e-310 - 1.0).abs() < 1e-12);
}

#[test]
fn half_space_refusals_and_wall_queries_preserve_both_original_clock_owners() {
    let source = field();
    let original = source.checkpoint();
    let wall = MannPlaneWall::new(0.0).unwrap();
    assert!(MannPlaneWall::new(f64::NAN).is_err());
    for point in [[0.0, 0.0, -f64::MIN_POSITIVE], [f64::INFINITY, 0.0, 1.0]] {
        assert!(source.sample_wall_blocked(point, wall).is_err());
        assert!(source.covariance_wall_blocked([0.0; 3], point, 0.0, wall).is_err());
    }
    assert!(source.covariance_wall_blocked([0.0; 3], [0.0; 3], f64::NAN, wall).is_err());
    let point = [0.2, -0.3, 0.1];
    let mut fixed = fixed_step::FixedStepMannField::from_snapshot(
        source.clone(), 123, 0.0079, usize::MAX,
    ).unwrap();
    assert_eq!(fixed.sample_wall_blocked(point, wall), source.sample_wall_blocked(point, wall));
    assert_eq!(
        fixed.covariance_wall_blocked(point, [0.0, 0.0, 0.3], 0.21, wall),
        source.covariance_wall_blocked(point, [0.0, 0.0, 0.3], 0.21, wall)
    );
    let before = fixed.checkpoint();
    fixed.sample_wall_blocked(point, wall).unwrap();
    assert_eq!(fixed.checkpoint(), before);
    fixed.advance_to(3, 4 * 32 * 3, || false).unwrap();
    let restored = fixed_step::FixedStepMannField::resume(fixed.checkpoint(), usize::MAX).unwrap();
    assert_eq!(fixed.sample_wall_blocked(point, wall), restored.sample_wall_blocked(point, wall));
    assert_eq!(source.checkpoint(), original);
}
