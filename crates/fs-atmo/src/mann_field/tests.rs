use super::*;
use crate::mann::MANN_TARGET_V1;
use crate::ou::STATIONARY_ANCHOR_TICK;

fn model(gamma: f64, count: usize) -> MannSpectralModel {
    MannSpectralModel::build(
        41,
        MannParams {
            gamma,
            ..MANN_TARGET_V1
        },
        MannQuadrature {
            n_modes: count,
            min_kl: 0.03,
            max_kl: 30.0,
        },
        [7.0, -0.5, 0.2],
    )
    .unwrap()
}

// Independent covariance oracle: integrate the original tensor directly,
// without using the polarization factors or production covariance evaluator.
fn tensor_covariance(
    model: &MannSpectralModel,
    separation: [f64; 3],
    lag: f64,
    times: &[f64],
) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for (mode, &tau) in model.modes.iter().zip(times) {
        let phi = mann_tensor(&model.parameters, mode.k).unwrap();
        let phase: f64 = (0..3)
            .map(|j| mode.k[j] * (separation[j] - model.advection_mps[j] * lag))
            .sum();
        let scale = mode.weight * phase.cos() * (-lag.abs() / tau).exp();
        for i in 0..3 {
            for j in 0..3 {
                result[i][j] += scale * phi[i][j];
            }
        }
    }
    result
}

#[test]
fn all_modal_auto_and_cross_spectra_match_the_mann_tensor() {
    for gamma in [0.0, 3.9, 10.0] {
        let model = model(gamma, MAX_MODES);
        for mode in &model.modes {
            let phi = mann_tensor(&model.parameters, mode.k).unwrap();
            let scale = phi[0][0].max(phi[1][1]).max(phi[2][2]) * mode.weight;
            for i in 0..3 {
                for j in 0..3 {
                    let represented = mode.polarization[0][i] * mode.polarization[0][j]
                        + mode.polarization[1][i] * mode.polarization[1][j];
                    assert!(
                        (represented - mode.weight * phi[i][j]).abs() < 2e-10 * scale,
                        "gamma={gamma}, k={:?}, [{i},{j}]: {represented} versus {}",
                        mode.k,
                        mode.weight * phi[i][j]
                    );
                }
            }
            for polarization in mode.polarization {
                let norm = dot(mode.k, mode.k).sqrt() * dot(polarization, polarization).sqrt();
                assert!(dot(mode.k, polarization).abs() < 3e-15 * norm);
            }
        }
    }
}

#[test]
fn isotropic_energy_matches_an_independent_radial_integral_without_renormalization() {
    let model = model(0.0, MAX_MODES);
    let p = model.parameters();
    // Incompressible isotropic Phi has trace E/(2 pi k^2), hence total
    // represented kinetic variance is 2 integral E(k) dk. Integrate this
    // independent scalar expression with a dense midpoint rule in log(kL).
    let q = model.quadrature();
    let bins = 20_000;
    let dlog = (q.max_kl / q.min_kl).ln() / bins as f64;
    let mut expected_trace = 0.0;
    for bin in 0..bins {
        let kl = q.min_kl * ((bin as f64 + 0.5) * dlog).exp();
        let energy = p.alpha_eps23 * p.length_m.powf(5.0 / 3.0) * kl.powi(4)
            / (1.0 + kl * kl).powf(17.0 / 6.0);
        expected_trace += 2.0 * energy * kl / p.length_m * dlog;
    }
    let covariance = model.covariance([0.0; 3]).unwrap();
    let trace: f64 = (0..3).map(|i| covariance[i][i]).sum();
    assert!(
        (trace / expected_trace - 1.0).abs() < 0.005,
        "trace={trace}, independent integral={expected_trace}"
    );
    for i in 0..3 {
        assert!((covariance[i][i] / (expected_trace / 3.0) - 1.0).abs() < 0.08);
    }
    let narrow = MannSpectralModel::build(
        41,
        p,
        MannQuadrature {
            min_kl: 1.0,
            max_kl: 2.0,
            ..q
        },
        [0.0; 3],
    )
    .unwrap();
    let narrow_cov = narrow.covariance([0.0; 3]).unwrap();
    let narrow_trace: f64 = (0..3).map(|i| narrow_cov[i][i]).sum();
    assert!(
        narrow_trace < 0.4 * trace,
        "cutting the band must remove energy, not rescale it back"
    );
}

#[test]
fn independent_realizations_recover_sheared_oblique_two_point_covariance() {
    let model = model(3.9, 64);
    let times = vec![0.8; 64];
    let zero = tensor_covariance(&model, [0.0; 3], 0.0, &times);
    let separations = [[0.31, -0.24, 0.17], [1.1, 0.7, -0.9]];
    let expected = separations.map(|r| tensor_covariance(&model, r, 0.0, &times));
    assert!(
        zero[0][2] < -0.1 * (zero[0][0] * zero[2][2]).sqrt(),
        "nonzero shear is essential to this oracle"
    );
    for (r, target) in separations.iter().zip(expected) {
        let represented = model.covariance(*r).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (represented[i][j] - target[i][j]).abs()
                        < 1e-11 * (zero[i][i] * zero[j][j]).sqrt()
                );
            }
        }
    }
    let origin = [0.6, -0.2, 0.4];
    let translated = [3.1, 1.4, -2.0];
    let mut observed = [[[0.0; 3]; 3]; 2];
    let mut translated_stress = [[0.0; 3]; 3];
    let samples = 4096;
    for seed in 0..samples {
        let field =
            HomogeneousMannField::at_anchor(model.clone(), seed as u64, times.clone()).unwrap();
        let a = field.sample(origin).unwrap().u;
        let c = field.sample(translated).unwrap().u;
        for (index, r) in separations.iter().enumerate() {
            let b = field
                .sample(core::array::from_fn(|j| origin[j] + r[j]))
                .unwrap()
                .u;
            for i in 0..3 {
                for j in 0..3 {
                    observed[index][i][j] += a[i] * b[j] / samples as f64;
                }
            }
        }
        for i in 0..3 {
            for j in 0..3 {
                translated_stress[i][j] += c[i] * c[j] / samples as f64;
            }
        }
    }
    for index in 0..2 {
        for i in 0..3 {
            for j in 0..3 {
                let scale = (zero[i][i] * zero[j][j]).sqrt();
                assert!(
                    (observed[index][i][j] - expected[index][i][j]).abs() < 0.085 * scale,
                    "r={:?}, [{i},{j}]: {} vs {}",
                    separations[index],
                    observed[index][i][j],
                    expected[index][i][j]
                );
                assert!(
                    (translated_stress[i][j] - zero[i][j]).abs() < 0.085 * scale,
                    "translation must not change stationary stress [{i},{j}]"
                );
            }
        }
    }
}

#[test]
fn temporal_covariance_retains_vector_advection_and_ou_decay() {
    let model = model(3.9, 4);
    let times = vec![0.15, 0.2, 0.3, 0.4];
    let lag_ticks = 24;
    let lag = lag_ticks as f64 / TICK_HZ;
    let moving: [f64; 3] = model.advection_mps.map(|v| v * lag);
    let zero = tensor_covariance(&model, [0.0; 3], 0.0, &times);
    let expected = [[0.0; 3], moving].map(|r| tensor_covariance(&model, r, lag, &times));
    let mut observed = [[[0.0; 3]; 3]; 2];
    let origin = [-0.3, 0.1, 0.5];
    let samples = 2048;
    for seed in 0..samples {
        let mut field =
            HomogeneousMannField::at_anchor(model.clone(), seed as u64, times.clone()).unwrap();
        let a = field.sample(origin).unwrap().u;
        if seed == 0 {
            for (index, r) in [[0.0; 3], moving].iter().enumerate() {
                let represented = field.covariance(*r, lag).unwrap();
                for i in 0..3 {
                    for j in 0..3 {
                        assert!(
                            (represented[i][j] - expected[index][i][j]).abs()
                                < 1e-11 * (zero[i][i] * zero[j][j]).sqrt()
                        );
                    }
                }
            }
        }
        field
            .advance_to(STATIONARY_ANCHOR_TICK + lag_ticks)
            .unwrap();
        for (index, r) in [[0.0; 3], moving].iter().enumerate() {
            let b = field
                .sample(core::array::from_fn(|j| origin[j] + r[j]))
                .unwrap()
                .u;
            for i in 0..3 {
                for j in 0..3 {
                    observed[index][i][j] += a[i] * b[j] / samples as f64;
                }
            }
        }
    }
    for index in 0..2 {
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (observed[index][i][j] - expected[index][i][j]).abs()
                        < 0.1 * (zero[i][i] * zero[j][j]).sqrt(),
                    "temporal covariance {index}[{i},{j}]: {} vs {}",
                    observed[index][i][j],
                    expected[index][i][j]
                );
            }
        }
    }
    assert!(
        (expected[0][0][0] - expected[1][0][0]).abs() > 0.01 * zero[0][0],
        "advection must change the fixed-probe covariance"
    );
}

#[test]
fn gradients_are_solenoidal_and_this_homogeneous_field_has_no_wall() {
    let model = model(3.9, 32);
    let field = HomogeneousMannField::at_anchor(model, 27, vec![1.0; 32]).unwrap();
    let checkpoint = field.checkpoint();
    for point in [[0.1, 0.2, 0.0], [1.2, -0.7, -0.4]] {
        let sample = field.sample(point).unwrap();
        let trace: f64 = (0..3).map(|i| sample.grad[i][i]).sum();
        let scale: f64 = sample.grad.iter().flatten().map(|v| v.abs()).sum();
        assert!(trace.abs() < 1e-14 * scale);
        assert!(
            sample.u[2].abs() > 1e-8,
            "z=0 must not silently acquire a wall projection"
        );
        let h = 1e-5;
        for j in 0..3 {
            let mut left = point;
            let mut right = point;
            left[j] -= h;
            right[j] += h;
            let a = field.sample(left).unwrap();
            let b = field.sample(right).unwrap();
            for i in 0..3 {
                let finite_difference = (b.u[i] - a.u[i]) / (2.0 * h);
                assert!((finite_difference - sample.grad[i][j]).abs() < 2e-7 * scale.max(1.0));
            }
        }
    }
    assert_eq!(
        field.checkpoint(),
        checkpoint,
        "observation must not move any random path"
    );
}

#[test]
fn checkpoints_resume_exactly_and_invalid_models_or_advances_refuse() {
    let model = model(3.9, 8);
    let mut field = HomogeneousMannField::at_anchor(model.clone(), 22, vec![0.4; 8]).unwrap();
    field.advance_to(STATIONARY_ANCHOR_TICK + 7).unwrap();
    let checkpoint = field.checkpoint();
    let mut resumed = HomogeneousMannField::resume(checkpoint.clone()).unwrap();
    assert_eq!(resumed, field);
    for tick in [STATIONARY_ANCHOR_TICK + 11, STATIONARY_ANCHOR_TICK + 29] {
        field.advance_to(tick).unwrap();
        resumed.advance_to(tick).unwrap();
        assert_eq!(resumed, field);
        assert_eq!(
            resumed.sample([0.2, -0.6, 0.1]).unwrap(),
            field.sample([0.2, -0.6, 0.1]).unwrap()
        );
    }
    let before = field.clone();
    assert_eq!(
        field.advance_to(field.tick() - 1).unwrap_err().code,
        "ou-advance-backwards"
    );
    assert_eq!(
        field.advance_to(i64::MAX).unwrap_err().code,
        "ou-advance-span-exceeded"
    );
    assert_eq!(field, before);
    for component in 0..4 {
        let mut bad = checkpoint.clone();
        bad.paths[component].seed ^= 1;
        assert!(HomogeneousMannField::resume(bad).is_err());
        let mut bad = checkpoint.clone();
        bad.paths[component].tick += 1;
        assert!(HomogeneousMannField::resume(bad).is_err());
        let mut bad = checkpoint.clone();
        bad.paths[component].amplitudes[0] = f64::NAN;
        assert!(HomogeneousMannField::resume(bad).is_err());
    }
    assert!(HomogeneousMannField::at_anchor(model.clone(), 1, vec![0.0; 8]).is_err());
    assert!(HomogeneousMannField::at_anchor(model.clone(), 1, vec![1.0; 7]).is_err());
    assert!(field.sample([0.0, f64::NAN, 0.0]).is_err());
    for q in [
        MannQuadrature {
            n_modes: 0,
            ..model.quadrature()
        },
        MannQuadrature {
            n_modes: MAX_MODES + 1,
            ..model.quadrature()
        },
        MannQuadrature {
            min_kl: 0.0,
            ..model.quadrature()
        },
        MannQuadrature {
            max_kl: f64::INFINITY,
            ..model.quadrature()
        },
        MannQuadrature {
            max_kl: model.quadrature().min_kl,
            ..model.quadrature()
        },
    ] {
        assert!(MannSpectralModel::build(1, MANN_TARGET_V1, q, [0.0; 3]).is_err());
    }
}
