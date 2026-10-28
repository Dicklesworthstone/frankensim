//! G0/G3/G5: physical covariance, exact curl, and evolving-state replay.

use super::*;
use crate::ou::STATIONARY_ANCHOR_TICK;
use crate::{DEC17_AIR, Mode};

fn basis() -> TurbulenceField {
    TurbulenceField {
        modes: vec![Mode {
            kx: 0.4,
            ky: 0.7,
            kh: 0.3,
            a: [2.0, -1.0, 0.6],
            phi: 0.2,
        }],
        u_adv_mps: 6.0,
    }
}

fn field(seed: u64) -> EvolvingTurbulenceField {
    EvolvingTurbulenceField::at_anchor(basis(), seed, vec![0.6]).unwrap()
}

fn covariance_at_height(h: f64) -> [[f64; 3]; 3] {
    // Independent curl calculation for psi = [2 sin(kh h), -sin(kh h),
    // 0.6 cos(kh h)] cos(theta). These are the cosine/sine velocity
    // coefficients, so the phase-average tensor is (p p^T + q q^T)/2.
    let ch = (0.3 * h).cos();
    let sh = (0.3 * h).sin();
    let p = [0.3 * ch, 0.6 * ch, 0.0];
    let q = [-0.42 * ch, 0.24 * ch, 1.8 * sh];
    core::array::from_fn(|i| core::array::from_fn(|j| 0.5 * (p[i] * p[j] + q[i] * q[j])))
}

#[test]
fn paired_quadratures_preserve_full_stress_and_horizontal_homogeneity() {
    let checkpoint = field(7).checkpoint();
    let expected = covariance_at_height(4.0);
    for [x, y] in [[0.0, 0.0], [1.7, -4.3], [30.0, 12.0]] {
        let mut covariance = [[0.0; 3]; 3];
        // Four-point quadrature integrates this quadratic moment exactly:
        // E[c²] = E[s²] = 1/2, E[c s] = E[c] = E[s] = 0.
        for [c, s] in [[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0]] {
            let mut cp = checkpoint.clone();
            cp.cosine.amplitudes[0] = c;
            cp.sine.amplitudes[0] = s;
            let f = EvolvingTurbulenceField::resume(cp).unwrap();
            let u = f.sample(x, y, 4.0).unwrap().u;
            for i in 0..3 {
                for j in 0..3 {
                    covariance[i][j] += u[i] * u[j] * 0.25;
                }
            }
        }
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (covariance[i][j] - expected[i][j]).abs() < 2e-14,
                    "stress ({i},{j}) at ({x},{y}) changed: {} vs {}",
                    covariance[i][j],
                    expected[i][j]
                );
            }
        }
    }
    assert!(
        expected[0][2] < -0.05,
        "nonzero wall-normal shear-stress oracle"
    );
}

#[test]
fn ou_evolution_has_the_advected_stationary_covariance() {
    // Independent seeds sample the stationary ensemble at two times.
    // The per-component covariance is Cii exp(-dt/tau) cos(kx Uadv dt).
    let expected = covariance_at_height(4.0);
    let lag = 36;
    let dt = lag as f64 / TICK_HZ;
    let correlation = (-dt / 0.6).exp() * (0.4 * 6.0 * dt).cos();
    let mut variance0 = [0.0; 3];
    let mut variance1 = [0.0; 3];
    let mut lagged = [0.0; 3];
    let n = 4096;
    for seed in 0..n {
        let mut f = field(seed);
        let u0 = f.sample(1.7, -4.3, 4.0).unwrap().u;
        f.advance_to(STATIONARY_ANCHOR_TICK + lag).unwrap();
        let u1 = f.sample(1.7, -4.3, 4.0).unwrap().u;
        for i in 0..3 {
            variance0[i] += u0[i] * u0[i] / n as f64;
            variance1[i] += u1[i] * u1[i] / n as f64;
            lagged[i] += u0[i] * u1[i] / n as f64;
        }
    }
    for i in 0..3 {
        assert!((variance0[i] / expected[i][i] - 1.0).abs() < 0.07);
        assert!((variance1[i] / expected[i][i] - 1.0).abs() < 0.07);
        assert!(
            (lagged[i] / expected[i][i] - correlation).abs() < 0.07,
            "component {i}: correlation {} vs {correlation}",
            lagged[i] / expected[i][i]
        );
    }
}

#[test]
fn evolving_curl_keeps_wall_parity_and_exact_derivatives() {
    let mut f = field(19);
    for tick in [STATIONARY_ANCHOR_TICK, -37, 0, 240] {
        f.advance_to(tick).unwrap();
        for [x, y, h] in [[0.3, 1.4, 0.0], [1.7, -4.3, 4.0], [-7.0, 3.0, 0.1]] {
            let sample = f.sample(x, y, h).unwrap();
            let div = sample.grad[0][0] + sample.grad[1][1] + sample.grad[2][2];
            assert!(div.abs() < 1e-14);
            assert_eq!(f.sample(x, y, 0.0).unwrap().u[2], 0.0);
        }
    }
    let p = [1.7, -4.3, 4.0];
    let exact = f.sample(p[0], p[1], p[2]).unwrap();
    let error = |epsilon: f64| {
        let mut worst = 0.0_f64;
        for j in 0..3 {
            let mut lo = p;
            let mut hi = p;
            lo[j] -= epsilon;
            hi[j] += epsilon;
            let left = f.sample(lo[0], lo[1], lo[2]).unwrap();
            let right = f.sample(hi[0], hi[1], hi[2]).unwrap();
            for i in 0..3 {
                let fd = (right.u[i] - left.u[i]) / (2.0 * epsilon);
                worst = worst.max((fd - exact.grad[i][j]).abs());
            }
        }
        worst
    };
    let errors = [error(0.02), error(0.01), error(0.005)];
    for pair in errors.windows(2) {
        assert!(((pair[0] / pair[1]).log2() - 2.0).abs() < 0.05);
    }
    assert!(errors[2] < 2e-6);
}

#[test]
fn queries_are_read_only_and_checkpoint_resume_is_bitwise() {
    let spatial = TurbulenceField::build(42, 16, 1.4, 30.0, 10.7).unwrap();
    let mut straight = EvolvingTurbulenceField::at_anchor(spatial, 55, vec![1.0; 16]).unwrap();
    let mut paused = straight.clone();
    straight.advance_to(480).unwrap();
    paused.advance_to(-20).unwrap();
    paused = EvolvingTurbulenceField::resume(paused.checkpoint()).unwrap();
    paused.advance_to(0).unwrap();
    let state = paused.checkpoint();
    let points = [[1.0, 2.0, 3.0], [-1.0, 7.0, 0.0], [5.0, 8.0, 2.0]];
    let mut output = [FlowSample {
        u: [0.0; 3],
        grad: [[0.0; 3]; 3],
    }; 3];
    paused.sample_batch(&points, &mut output).unwrap();
    for i in (0..points.len()).rev() {
        let [x, y, h] = points[i];
        assert_eq!(output[i], paused.sample(x, y, h).unwrap());
    }
    assert_eq!(state, paused.checkpoint());
    paused.advance_to(480).unwrap();
    assert_eq!(straight, paused);
    for [x, y, h] in points {
        assert_eq!(
            straight.sample(x, y, h).unwrap(),
            paused.sample(x, y, h).unwrap()
        );
    }
}

#[test]
fn admission_refuses_mismatched_state_without_partial_mutation() {
    let mut f = field(3);
    let before = f.clone();
    for target in [STATIONARY_ANCHOR_TICK - 1, i64::MAX] {
        assert!(f.advance_to(target).is_err());
        assert_eq!(f, before);
    }
    let mut cp = f.checkpoint();
    cp.sine.tick += 1;
    assert_eq!(
        EvolvingTurbulenceField::resume(cp).unwrap_err().code,
        "evolving-checkpoint-invalid"
    );
    let mut cp = f.checkpoint();
    cp.sine.seed = cp.cosine.seed;
    assert_eq!(
        EvolvingTurbulenceField::resume(cp).unwrap_err().code,
        "evolving-checkpoint-invalid"
    );
    let mut cp = f.checkpoint();
    cp.cosine.amplitudes[0] = f64::NAN;
    assert_eq!(
        EvolvingTurbulenceField::resume(cp).unwrap_err().code,
        "ou-checkpoint-invalid"
    );
    assert!(EvolvingTurbulenceField::at_anchor(basis(), 3, vec![]).is_err());
    for tau in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert_eq!(
            EvolvingTurbulenceField::at_anchor(basis(), 3, vec![tau])
                .unwrap_err()
                .code,
            "ou-params-invalid"
        );
    }
    let sentinel = FlowSample {
        u: [19.0; 3],
        grad: [[19.0; 3]; 3],
    };
    let mut output = [sentinel; 2];
    assert!(
        f.sample_batch(&[[1.0; 3], [0.0, 0.0, -1.0]], &mut output)
            .is_err()
    );
    assert_eq!(output, [sentinel; 2]);
    assert!(f.sample(f64::NAN, 0.0, 0.0).is_err());
    assert!(f.sample_batch(&[], &mut output).is_err());
    assert!(f.sample_batch(&[], &mut []).is_ok());
    let largest = TurbulenceField::build(1, MAX_MODES, 1.0, 20.0, 5.0).unwrap();
    assert!(EvolvingTurbulenceField::at_anchor(largest, 2, vec![1.0; MAX_MODES]).is_ok());
}

#[test]
fn evolving_air_state_drives_velocity_gradient_pressure_and_reynolds_together() {
    let mut atmo = Atmosphere {
        mean: FlatSiteLogLaw {
            scenario_effective_z0_m: 0.003,
            displacement_height_m: 0.02,
            reference_height_m: 6.0,
            reference_speed_mps: 10.73,
        },
        turbulence: basis(),
        air: DEC17_AIR,
    }
    .into_evolving(17, vec![0.6])
    .unwrap();
    atmo.turbulence.advance_to(0).unwrap();
    let turbulent = atmo.turbulence.sample(1.0, 2.0, 4.0).unwrap();
    let air = atmo.sample_air_state(1.0, 2.0, 4.0).unwrap();
    assert_eq!(air.velocity_mps[0], turbulent.u[0] + atmo.mean.speed(4.0));
    assert_eq!(air.velocity_mps[1..], turbulent.u[1..]);
    assert_eq!(
        air.grad[0][2],
        turbulent.grad[0][2] + atmo.mean.dspeed_dh(4.0)
    );
    let speed = air.velocity_mps.iter().map(|u| u * u).sum::<f64>().sqrt();
    assert!((air.dynamic_pressure_pa() - 0.5 * DEC17_AIR.rho_kg_m3 * speed * speed).abs() < 1e-12);
    assert!(
        (air.reynolds(2.0) - DEC17_AIR.rho_kg_m3 * speed * 2.0 / DEC17_AIR.mu_kg_m_s).abs() < 1e-9
    );
    assert_eq!(air.provenance, DEC17_AIR.provenance);
    atmo.turbulence.advance_to(1).unwrap();
    assert_ne!(
        air.velocity_mps,
        atmo.sample_air_state(1.0, 2.0, 4.0).unwrap().velocity_mps
    );
}
