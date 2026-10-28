//! Full aperture radiation stays coupled to the retained bore modes.
use super::*;
use crate::flanged::{FLANGED_MAX_KA, FLANGED_MAX_MODES, flanged_modal_impedance};
use crate::profile::{RadiusProfile, RiccatiSteps};
use fs_material::gas::GasSpec;

fn air(temperature: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 101_325.0).unwrap()
}

fn matvec(matrix: &[C64], flow: &[C64]) -> Vec<C64> {
    let n = flow.len();
    (0..n)
        .map(|i| (0..n).fold(C64::ZERO, |sum, j| sum + matrix[i * n + j] * flow[j]))
        .collect()
}

fn power(pressure: &[C64], flow: &[C64]) -> f64 {
    pressure
        .iter()
        .zip(flow)
        .map(|(p, u)| (*p * u.conj()).re * 0.5)
        .sum()
}

// Independent physical pressure/flow transfer of a uniform cylinder:
// Pin = C Pout - i Zc S Uout; Uin = -i S/Zc Pout + C Uout.
fn cylinder_reference(
    gas: &GasState,
    radius: f64,
    length: f64,
    omega: f64,
    boundary: &[C64],
    n: usize,
) -> Vec<C64> {
    let mut roots = vec![0.0];
    roots.extend(j1_roots(n - 1).unwrap());
    let k0 = omega / gas.sound_speed;
    let z0 = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    let mut a = vec![C64::ZERO; n * n];
    let mut b = a.clone();
    for i in 0..n {
        let k = C64::from_re(k0 * k0 - (roots[i] / radius).powi(2)).sqrt();
        let zc = C64::from_re(z0 * k0) / k;
        let forward = crate::exp_i(k.scale(length));
        let backward = crate::exp_i(k.scale(-length));
        let cosine = (forward + backward).scale(0.5);
        let sine = (forward - backward) * C64::new(0.0, -0.5);
        let shunt = C64::new(0.0, -1.0) * sine / zc;
        for j in 0..n {
            a[i * n + j] = cosine * boundary[i * n + j];
            b[i * n + j] = shunt * boundary[i * n + j];
        }
        a[i * n + i] = a[i * n + i] + C64::new(0.0, -1.0) * zc * sine;
        b[i * n + i] = b[i * n + i] + cosine;
    }
    mat_mul(&a, &mat_inv(&b, n).unwrap(), n)
}

#[test]
fn full_flanged_matrix_preserves_signed_aperture_power() {
    let gas = air(293.15);
    let radius = 0.015;
    let omega = 1.2 * gas.sound_speed / radius;
    let n = 4;
    let z = flanged_modal_impedance(&gas, radius, omega, n).unwrap();
    let one = flanged_modal_impedance(&gas, radius, omega, 1).unwrap();
    let (r, x) =
        fs_phs::baffled_piston_impedance(gas.density, gas.sound_speed, radius, omega).unwrap();
    assert_eq!(one[0].re.to_bits(), r.to_bits());
    assert_eq!(one[0].im.to_bits(), x.to_bits());
    for i in 0..n {
        for j in 0..n {
            assert_eq!(z[i * n + j], z[j * n + i], "reciprocal aperture");
        }
    }
    let flow = [
        C64::new(0.8, 0.1),
        C64::new(-0.3, 0.7),
        C64::new(0.2, -0.6),
        C64::new(-0.4, -0.2),
    ];
    let actual = power(&matvec(&z, &flow), &flow);
    // Independent Rayleigh far-field power, using physical radial integrals,
    // not the production closed-form Hankel transforms or packed matrix.
    let nodes = crate::levine_schwinger::gauss_legendre(48);
    let mut roots = vec![0.0];
    roots.extend(j1_roots(n - 1).unwrap());
    let ka = omega * radius / gas.sound_speed;
    let mut hemispherical = 0.0;
    for &(cos_theta, angular_weight) in &nodes {
        let transverse = ka * det::sqrt(1.0 - cos_theta * cos_theta);
        let mut transform = C64::ZERO;
        for i in 0..n {
            let norm = bessel_j0_certified(roots[i]).abs();
            let radial: f64 = nodes
                .iter()
                .map(|&(eta, weight)| {
                    weight * eta * bessel_j0_certified(roots[i] * eta) / norm
                        * bessel_j0_certified(transverse * eta)
                })
                .sum();
            transform = transform + flow[i].scale(radial);
        }
        hemispherical += angular_weight * transform.abs().powi(2);
    }
    let expected =
        gas.density * omega * omega / (core::f64::consts::PI * gas.sound_speed) * hemispherical;
    assert!(actual > 0.0);
    assert!(
        (actual - expected).abs() < 2.0e-11 * expected,
        "{actual:e} vs {expected:e}"
    );
    let diagonal: f64 = (0..n)
        .map(|i| flow[i].abs().powi(2) * z[i * n + i].re * 0.5)
        .sum();
    assert!(
        (actual - diagonal).abs() > 1.0e-3 * actual,
        "mutual terms must do real work"
    );
}

#[test]
fn full_flanged_cylinder_matches_physical_transfer_and_work() {
    let gas = air(343.15);
    let exterior = air(293.15);
    let (radius, length, n) = (0.04, 0.013, 4);
    let omega = 4.8 * gas.sound_speed / radius;
    let z = flanged_modal_impedance(&exterior, radius, omega, n).unwrap();
    let expected = cylinder_reference(&gas, radius, length, omega, &z, n);
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let profile = RadiusProfile::exponential(radius, 0.0, length).unwrap();
    let stepped =
        mm_input_impedance_flanged(&duct, &gas, &exterior, omega, LossModel::Lossless, n, 1)
            .unwrap();
    let continuous = mm_profile_impedance_flanged(
        &profile,
        &gas,
        &exterior,
        omega,
        LossModel::Lossless,
        ModeCount(n),
        RiccatiSteps(32),
    )
    .unwrap();
    let scale = expected.iter().map(|z| z.abs()).fold(0.0, f64::max);
    for result in [&stepped, &continuous] {
        for (a, b) in result.impedance_matrix.iter().zip(&expected) {
            assert!((*a - *b).abs() <= 2.0e-10 * scale, "{a:?} vs {b:?}");
        }
        assert_eq!(result.mouth_ka, omega * radius / exterior.sound_speed);
    }
    let flow = [
        C64::new(1.0e-6, 0.2e-6),
        C64::new(-0.3e-6, 0.4e-6),
        C64::new(0.1e-6, -0.2e-6),
        C64::new(0.2e-6, 0.1e-6),
    ];
    let field = mm_profile_field_flanged(
        &profile,
        &gas,
        &exterior,
        omega,
        LossModel::Lossless,
        ModeCount(n),
        RiccatiSteps(32),
        &flow,
    )
    .unwrap();
    let mouth_p = matvec(&z, field.flow.last().unwrap());
    for (actual, expected) in field.pressure.last().unwrap().iter().zip(mouth_p) {
        assert!((*actual - expected).abs() <= 2.0e-10 * scale * 1.0e-6);
    }
    let radiated = power(field.pressure.last().unwrap(), field.flow.last().unwrap());
    assert!(radiated > 0.0);
    for (p, u) in field.pressure.iter().zip(&field.flow) {
        assert!((power(p, u) - radiated).abs() <= 2.0e-10 * radiated);
    }
    let in_p = matvec(&expected, &flow);
    for (actual, expected) in field.pressure[0].iter().zip(in_p) {
        assert!((*actual - expected).abs() <= 2.0e-10 * scale * 1.0e-6);
    }
}

#[test]
fn full_flanged_continuous_field_is_regular_at_transverse_cutoff() {
    let gas = air(293.15);
    let radius = 0.03;
    let profile = RadiusProfile::exponential(radius, 0.0, 0.02).unwrap();
    let cutoff = j1_roots(1).unwrap()[0] * gas.sound_speed / radius;
    let responses: Vec<_> = [1.0 - 1.0e-7, 1.0, 1.0 + 1.0e-7]
        .into_iter()
        .map(|factor| {
            mm_profile_impedance_flanged(
                &profile,
                &gas,
                &gas,
                cutoff * factor,
                LossModel::Lossless,
                ModeCount(2),
                RiccatiSteps(32),
            )
            .unwrap()
        })
        .collect();
    let scale = responses[1]
        .impedance_matrix
        .iter()
        .map(|z| z.abs())
        .fold(0.0, f64::max);
    for edge in [&responses[0], &responses[2]] {
        for (a, b) in edge
            .impedance_matrix
            .iter()
            .zip(&responses[1].impedance_matrix)
        {
            assert!((*a - *b).abs() < 2.0e-5 * scale);
        }
    }
}

fn peak_width(values: &[(f64, f64)]) -> (f64, f64) {
    let peak = (1..values.len() - 1)
        .max_by(|&i, &j| values[i].1.total_cmp(&values[j].1))
        .unwrap();
    let (left, centre, right) = (values[peak - 1].1, values[peak].1, values[peak + 1].1);
    assert!(centre > left && centre > right, "resolved interior peak");
    let shift = 0.5 * (left - right) / (left - 2.0 * centre + right);
    let frequency = values[peak].0 + shift * (values[peak + 1].0 - values[peak].0);
    let half = (centre - 0.25 * (left - right) * shift) * 0.5;
    let crossing = |i: usize| {
        let (f0, y0) = values[i];
        let (f1, y1) = values[i + 1];
        f0 + (half - y0) / (y1 - y0) * (f1 - f0)
    };
    let lo = (0..peak)
        .rev()
        .find(|&i| values[i].1 <= half)
        .expect("low half-power crossing");
    let hi = (peak..values.len() - 1)
        .find(|&i| values[i + 1].1 <= half)
        .expect("high half-power crossing");
    (frequency, crossing(hi) - crossing(lo))
}

#[test]
fn complete_mouth_changes_resonance_and_linewidth_on_a_mode_ladder() {
    let gas = air(293.15);
    let (radius, length) = (0.025, 0.01);
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let mut summaries = Vec::new();
    for n in [1, 2, 4] {
        let mut spectrum = Vec::new();
        for i in 0..=400 {
            let ka = 0.15 + 2.5 * i as f64 / 400.0;
            let omega = ka * gas.sound_speed / radius;
            let result =
                mm_input_impedance_flanged(&duct, &gas, &gas, omega, LossModel::Lossless, n, 1)
                    .unwrap();
            if n == 4 {
                let boundary = flanged_modal_impedance(&gas, radius, omega, n).unwrap();
                let expected = cylinder_reference(&gas, radius, length, omega, &boundary, n)[0];
                assert!((result.plane_impedance - expected).abs() < 2.0e-10 * expected.abs());
            }
            spectrum.push((ka, result.plane_impedance.abs().powi(2)));
        }
        let (peak, width) = peak_width(&spectrum);
        summaries.push((peak, peak / width));
        eprintln!(
            "full modal mouth n={n}: ka_peak={peak:e}, Q={:e}",
            peak / width
        );
    }
    let movement = |i: usize, component: usize| {
        let a = if component == 0 {
            summaries[i].0
        } else {
            summaries[i].1
        };
        let b = if component == 0 {
            summaries[i + 1].0
        } else {
            summaries[i + 1].1
        };
        (a - b).abs()
    };
    for component in 0..2 {
        assert!(
            movement(0, component) > 1.0e-4,
            "mode coupling must change the resonance and Q"
        );
        assert!(
            movement(1, component) < movement(0, component),
            "refined mode ladder: {summaries:?}"
        );
    }
}

#[test]
fn full_flanged_mouth_refuses_invalid_media_and_unresolved_domains() {
    let gas = air(293.15);
    let omega = 4000.0;
    for n in [0, FLANGED_MAX_MODES + 1] {
        assert!(matches!(
            flanged_modal_impedance(&gas, 0.01, omega, n),
            Err(DuctError::BadParameter { .. })
        ));
    }
    assert!(matches!(
        flanged_modal_impedance(&gas, 0.01, 0.0, 2),
        Err(DuctError::BadParameter { .. })
    ));
    assert!(matches!(
        flanged_modal_impedance(
            &gas,
            0.01,
            (FLANGED_MAX_KA + 1.0) * gas.sound_speed / 0.01,
            2
        ),
        Err(DuctError::RadiationKaTooLarge { .. })
    ));
    let exterior = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 90_000.0).unwrap();
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.01,
            length: 0.1,
        }],
    };
    let profile = RadiusProfile::exponential(0.01, 0.0, 0.1).unwrap();
    assert!(
        mm_input_impedance_flanged(&duct, &gas, &exterior, omega, LossModel::Lossless, 2, 1)
            .is_err()
    );
    assert!(
        mm_profile_impedance_flanged(
            &profile,
            &gas,
            &exterior,
            omega,
            LossModel::Lossless,
            ModeCount(2),
            RiccatiSteps(32)
        )
        .is_err()
    );
}
