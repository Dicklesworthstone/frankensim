use super::*;
use crate::modal::*;
use crate::profile::{RadiusProfile, RiccatiSteps};
use crate::{Duct, LossModel, Segment};
use fs_material::gas::GasSpec;

fn air(temperature: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), temperature, 101_325.0).unwrap()
}

fn matvec(matrix: &[C64], vector: &[C64]) -> Vec<C64> {
    let n = vector.len();
    (0..n)
        .map(|i| (0..n).fold(C64::ZERO, |sum, j| sum + matrix[i * n + j] * vector[j]))
        .collect()
}

fn power(pressure: &[C64], flow: &[C64]) -> f64 {
    pressure
        .iter()
        .zip(flow)
        .map(|(p, u)| (*p * u.conj()).re * 0.5)
        .sum()
}

fn flanged(gas: GasState, radius: f64, omega: f64, modes: usize) -> ModalMouthImpedance {
    let z = crate::flanged::flanged_modal_impedance(&gas, radius, omega, modes).unwrap();
    ModalMouthImpedance::new(omega, radius, gas, ModeCount(modes), z).unwrap()
}

#[test]
fn original_signed_radial_basis_has_unit_area_gram_and_zero_higher_net_flow() {
    // More modes than the staircase cap: this is also the continuous/BEM basis.
    let basis = AxisymmetricMouthBasis::new(ModeCount(12)).unwrap();
    let quadrature = crate::levine_schwinger::gauss_legendre(96);
    assert_eq!(basis.radial_eigenvalues()[0], 0.0);
    for i in 0..12 {
        let mean: f64 = quadrature
            .iter()
            .map(|&(r, w)| 2.0 * r * w * basis.shape(i, r).unwrap())
            .sum();
        assert!((mean - if i == 0 { 1.0 } else { 0.0 }).abs() < 3.0e-12);
        for j in 0..12 {
            let gram: f64 = quadrature
                .iter()
                .map(|&(r, w)| {
                    2.0 * r * w * basis.shape(i, r).unwrap() * basis.shape(j, r).unwrap()
                })
                .sum();
            assert!(
                (gram - if i == j { 1.0 } else { 0.0 }).abs() < 3.0e-12,
                "modal Gram ({i},{j})={gram:e}"
            );
        }
    }
    assert_eq!(basis.shape(0, 0.73).unwrap(), 1.0);
    assert!(basis.shape(1, 0.0).unwrap() > 0.0);
    assert_eq!(basis.shape(1, 1.0).unwrap(), -1.0);
    let maximum = AxisymmetricMouthBasis::new(ModeCount(certified_mode_cap())).unwrap();
    assert!(
        maximum
            .shape(maximum.modes().0 - 1, 0.73)
            .unwrap()
            .is_finite()
    );
    for modes in [0, certified_mode_cap() + 1] {
        assert!(AxisymmetricMouthBasis::new(ModeCount(modes)).is_err());
    }
    for radial in [-0.1, 1.1, f64::NAN] {
        assert!(basis.shape(0, radial).is_err());
    }
    assert!(basis.shape(12, 0.5).is_err());
}

#[test]
fn complete_hermitian_passivity_preserves_entries_and_resolves_small_loss() {
    let gas = air(293.15);
    let samples = [
        vec![
            C64::new(2.0e-8, 1.0e9),
            C64::new(0.3e-8, 0.5e9),
            C64::new(-0.1e-8, 0.5e9),
            C64::new(3.0e-8, 2.0e9),
        ],
        vec![C64::ONE; 4], // Singular positive semidefinite resistance.
        vec![
            C64::new(0.0, -2.0),
            C64::new(0.0, 0.5),
            C64::new(0.0, 0.5),
            C64::new(0.0, -1.0),
        ],
        vec![C64::ZERO, C64::ONE, -C64::ONE, C64::ZERO], // Lossless nonreciprocal.
        vec![
            C64::from_re(2.0),
            C64::new(0.0, 0.5),
            C64::new(0.0, -0.5),
            C64::ONE,
        ], // Complex Hermitian, passive.
    ];
    for z in samples {
        let load = ModalMouthImpedance::new(3000.0, 0.01, gas, ModeCount(2), z.clone()).unwrap();
        assert_eq!(load.impedance(), z);
        assert_eq!(load.basis().modes(), ModeCount(2));
        assert_eq!(load.exterior(), &gas);
    }
    // Positive diagonal resistance alone cannot establish passivity.
    let active = vec![C64::ONE, C64::new(0.0, 2.0), C64::new(0.0, -2.0), C64::ONE];
    assert!(ModalMouthImpedance::new(3000.0, 0.01, gas, ModeCount(2), active).is_err());
    // Enormous reactive energy must not mask even very small active loss.
    for resistance in [-1.0e-30, -f64::from_bits(1)] {
        assert!(
            ModalMouthImpedance::new(
                3000.0,
                0.01,
                gas,
                ModeCount(1),
                vec![C64::new(resistance, 1.0e30)]
            )
            .is_err()
        );
    }
    for z in [vec![C64::ONE; 3], vec![C64::new(f64::NAN, 0.0); 4]] {
        assert!(ModalMouthImpedance::new(3000.0, 0.01, gas, ModeCount(2), z).is_err());
    }
    for (omega, radius) in [(0.0, 0.01), (3000.0, 0.0), (f64::INFINITY, 0.01)] {
        assert!(
            ModalMouthImpedance::new(omega, radius, gas, ModeCount(1), vec![C64::ONE]).is_err()
        );
    }
}

// Independent physical transfer on a uniform cylinder, using pressure/flow
// blocks rather than either production reflection or Magnus propagation.
fn cylinder_reference(
    gas: &GasState,
    radius: f64,
    length: f64,
    omega: f64,
    boundary: &[C64],
    n: usize,
) -> Vec<C64> {
    let basis = AxisymmetricMouthBasis::new(ModeCount(n)).unwrap();
    let k0 = omega / gas.sound_speed;
    let z0 = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
    let mut a = vec![C64::ZERO; n * n];
    let mut b = a.clone();
    for i in 0..n {
        let k = C64::from_re(k0 * k0 - (basis.radial_eigenvalues()[i] / radius).powi(2)).sqrt();
        let zc = C64::from_re(z0 * k0) / k;
        let forward = crate::exp_i(k.scale(length));
        let backward = crate::exp_i(k.scale(-length));
        let cosine = (forward + backward).scale(0.5);
        let sine = (forward - backward) * C64::new(0.0, -0.5);
        for j in 0..n {
            a[i * n + j] = cosine * boundary[i * n + j];
            b[i * n + j] = C64::new(0.0, -1.0) * sine / zc * boundary[i * n + j];
        }
        a[i * n + i] = a[i * n + i] + C64::new(0.0, -1.0) * zc * sine;
        b[i * n + i] = b[i * n + i] + cosine;
    }
    mat_mul(&a, &mat_inv(&b, n).unwrap(), n)
}

#[test]
fn full_matrix_matches_independent_line_and_preserves_work_with_evanescent_modes() {
    let gas = air(343.15);
    let exterior = air(293.15);
    let (radius, length, n) = (0.04, 0.013, 4);
    let duct = Duct {
        segments: vec![Segment::Cylinder { radius, length }],
    };
    let profile = RadiusProfile::exponential(radius, 0.0, length).unwrap();
    let prepared = PreparedModalDuct::new(&duct, n, 1).unwrap();
    let flow = [
        C64::new(1.0e-6, 0.2e-6),
        C64::new(-0.3e-6, 0.4e-6),
        C64::new(0.1e-6, -0.2e-6),
        C64::new(0.2e-6, 0.1e-6),
    ];
    for ka in [0.8, 4.8] {
        let omega = ka * gas.sound_speed / radius;
        let load = flanged(exterior, radius, omega, n);
        let expected = cylinder_reference(&gas, radius, length, omega, load.impedance(), n);
        let scale = expected.iter().map(|z| z.abs()).fold(0.0, f64::max);
        let stepped =
            mm_input_impedance_matrix(&duct, &gas, omega, LossModel::Lossless, &load, n, 1)
                .unwrap();
        let cached = prepared
            .input_impedance_matrix(&gas, omega, LossModel::Lossless, &load)
            .unwrap();
        let continuous = mm_profile_impedance_matrix(
            &profile,
            &gas,
            omega,
            LossModel::Lossless,
            &load,
            ModeCount(n),
            RiccatiSteps(32),
        )
        .unwrap();
        assert_eq!(stepped.impedance_matrix, cached.impedance_matrix);
        for actual in [&stepped, &continuous] {
            for (a, b) in actual.impedance_matrix.iter().zip(&expected) {
                assert!((*a - *b).abs() < 3.0e-10 * scale);
            }
            assert_eq!(actual.mouth_ka, omega * radius / exterior.sound_speed);
        }
        let waves = prepared
            .input_reflection_matrix(&gas, omega, LossModel::Lossless, &load)
            .unwrap();
        let reference = gas.density * gas.sound_speed / (core::f64::consts::PI * radius * radius);
        let reflection = (stepped.plane_impedance - C64::from_re(reference))
            / (stepped.plane_impedance + C64::from_re(reference));
        assert!((waves.plane_reflection(reference).unwrap() - reflection).abs() < 3.0e-12);
        let mouth =
            outlet::mouth_response_matrix(&duct, &gas, omega, LossModel::Lossless, &load, n, 1)
                .unwrap();
        assert_eq!(mouth.input.impedance_matrix, stepped.impedance_matrix);
        let mouth_flow = matvec(&mouth.flow_map, &flow);
        let field = mm_profile_field_matrix(
            &profile,
            &gas,
            omega,
            LossModel::Lossless,
            &load,
            ModeCount(n),
            RiccatiSteps(32),
            &flow,
        )
        .unwrap();
        let mouth_pressure = matvec(load.impedance(), &mouth_flow);
        for (a, b) in mouth_flow.iter().zip(field.flow.last().unwrap()) {
            assert!((*a - *b).abs() < 3.0e-15);
        }
        for (a, b) in mouth_pressure.iter().zip(field.pressure.last().unwrap()) {
            assert!((*a - *b).abs() < 3.0e-10 * scale * 1.0e-6);
        }
        let outward = power(&mouth_pressure, &mouth_flow);
        assert!(outward > 0.0);
        for (p, u) in field.pressure.iter().zip(&field.flow) {
            assert!((power(p, u) - outward).abs() < 3.0e-10 * outward);
        }
    }
}

#[test]
fn generic_matrix_retains_native_flanged_paths_and_principal_tail_comparison() {
    let gas = air(293.15);
    let exterior = air(313.15);
    let profile = RadiusProfile::exponential(0.018, 8.0, 0.02).unwrap();
    let radius = profile.radius(profile.length()).unwrap();
    let omega = 1.3 * exterior.sound_speed / radius;
    let load = flanged(exterior, radius, omega, 4);
    let matrix = mm_profile_impedance_matrix(
        &profile,
        &gas,
        omega,
        LossModel::Bessel,
        &load,
        ModeCount(4),
        RiccatiSteps(40),
    )
    .unwrap();
    let native = mm_profile_impedance_flanged(
        &profile,
        &gas,
        &exterior,
        omega,
        LossModel::Bessel,
        ModeCount(4),
        RiccatiSteps(40),
    )
    .unwrap();
    assert_eq!(matrix.impedance_matrix, native.impedance_matrix);
    assert_eq!(matrix.tail_bound_estimate, native.tail_bound_estimate);
    let flow = [C64::from_re(1.0e-6), C64::ZERO, C64::ZERO, C64::ZERO];
    let matrix = mm_profile_field_matrix(
        &profile,
        &gas,
        omega,
        LossModel::Bessel,
        &load,
        ModeCount(4),
        RiccatiSteps(40),
        &flow,
    )
    .unwrap();
    let native = mm_profile_field_flanged(
        &profile,
        &gas,
        &exterior,
        omega,
        LossModel::Bessel,
        ModeCount(4),
        RiccatiSteps(40),
        &flow,
    )
    .unwrap();
    assert_eq!(matrix.pressure, native.pressure);
    assert_eq!(matrix.flow, native.flow);
    let duct = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: 0.018,
            outlet_radius: radius,
            length: 0.02,
        }],
    };
    let prepared = PreparedModalDuct::new(&duct, 4, 8).unwrap();
    assert!(prepared.mouth_radius_m() < radius);
    let load = flanged(exterior, prepared.mouth_radius_m(), omega, 4);
    let matrix =
        outlet::mouth_response_matrix(&duct, &gas, omega, LossModel::Bessel, &load, 4, 8).unwrap();
    let native =
        outlet::mouth_response_flanged(&duct, &gas, &exterior, omega, LossModel::Bessel, 4, 8)
            .unwrap();
    assert_eq!(matrix.input.impedance_matrix, native.input.impedance_matrix);
    assert_eq!(matrix.flow_map, native.flow_map);
}

#[test]
fn matrix_stays_frequency_geometry_basis_and_medium_bound() {
    let gas = air(293.15);
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.02,
            length: 0.03,
        }],
    };
    let profile = RadiusProfile::exponential(0.02, 0.0, 0.03).unwrap();
    let omega = 7000.0;
    let load = flanged(gas, 0.02, omega, 2);
    for wrong_omega in [f64::from_bits(omega.to_bits() + 1), omega * 1.01] {
        assert!(
            mm_input_impedance_matrix(&duct, &gas, wrong_omega, LossModel::Lossless, &load, 2, 1)
                .is_err()
        );
        assert!(
            mm_profile_impedance_matrix(
                &profile,
                &gas,
                wrong_omega,
                LossModel::Lossless,
                &load,
                ModeCount(2),
                RiccatiSteps(16)
            )
            .is_err()
        );
    }
    for n in [1, 3] {
        assert!(
            mm_input_impedance_matrix(&duct, &gas, omega, LossModel::Lossless, &load, n, 1)
                .is_err()
        );
        assert!(
            mm_profile_impedance_matrix(
                &profile,
                &gas,
                omega,
                LossModel::Lossless,
                &load,
                ModeCount(n),
                RiccatiSteps(16)
            )
            .is_err()
        );
    }
    let wrong_radius = flanged(gas, 0.021, omega, 2);
    assert!(
        mm_input_impedance_matrix(&duct, &gas, omega, LossModel::Lossless, &wrong_radius, 2, 1)
            .is_err()
    );
    let other_pressure = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 90_000.0).unwrap();
    assert!(
        mm_input_impedance_matrix(
            &duct,
            &other_pressure,
            omega,
            LossModel::Lossless,
            &load,
            2,
            1
        )
        .is_err()
    );
    let one = flanged(gas, 0.02, omega, 1);
    let matrix =
        mm_input_impedance_matrix(&duct, &gas, omega, LossModel::Lossless, &one, 1, 1).unwrap();
    let scalar = mm_input_impedance_load(
        &duct,
        &gas,
        omega,
        LossModel::Lossless,
        one.impedance()[0],
        1,
        1,
    )
    .unwrap();
    assert!(
        (matrix.plane_impedance - scalar.plane_impedance).abs()
            < 1.0e-13 * scalar.plane_impedance.abs()
    );
}

#[test]
fn full_matrix_continuous_field_is_regular_at_bare_transverse_cutoff() {
    let gas = air(293.15);
    let radius = 0.03;
    let profile = RadiusProfile::exponential(radius, 0.0, 0.02).unwrap();
    let omega = j1_roots(1).unwrap()[0] * gas.sound_speed / radius;
    let load = flanged(gas, radius, omega, 2);
    let flow = [C64::from_re(1.0e-6), C64::new(0.1e-6, 0.2e-6)];
    let field = mm_profile_field_matrix(
        &profile,
        &gas,
        omega,
        LossModel::Lossless,
        &load,
        ModeCount(2),
        RiccatiSteps(32),
        &flow,
    )
    .unwrap();
    let expected = matvec(load.impedance(), field.flow.last().unwrap());
    for (a, b) in field.pressure.last().unwrap().iter().zip(expected) {
        assert!((*a - b).abs() < 2.0e-11);
    }
    let outward = power(field.pressure.last().unwrap(), field.flow.last().unwrap());
    for (p, u) in field.pressure.iter().zip(&field.flow) {
        assert!((power(p, u) - outward).abs() < 2.0e-11 * outward);
    }
}

#[test]
fn retained_mutual_radiation_changes_resonance_frequency_and_bandwidth() {
    use crate::resonance::{ResonanceOptions, refine_response_resonance};
    let gas = air(293.15);
    let radius = 0.025;
    let duct = Duct {
        segments: vec![Segment::Cylinder {
            radius,
            length: 0.01,
        }],
    };
    let prepared = PreparedModalDuct::new(&duct, 4, 1).unwrap();
    let mut results = Vec::new();
    for keep_mutual in [true, false] {
        let result = refine_response_resonance(
            [
                0.15 * gas.sound_speed / radius,
                2.65 * gas.sound_speed / radius,
            ],
            ResonanceOptions::default(),
            |omega| {
                let mut z = crate::flanged::flanged_modal_impedance(&gas, radius, omega, 4)?;
                if !keep_mutual {
                    for i in 0..4 {
                        for j in 0..4 {
                            if i != j {
                                z[i * 4 + j] = C64::ZERO;
                            }
                        }
                    }
                }
                let load = ModalMouthImpedance::new(omega, radius, gas, ModeCount(4), z)?;
                Ok(prepared
                    .input_impedance_matrix(&gas, omega, LossModel::Lossless, &load)?
                    .plane_impedance)
            },
            || false,
        )
        .unwrap();
        eprintln!(
            "mutual={keep_mutual}: ka_peak={:e}, Q_width={:e}, evaluations={}",
            result.angular_frequency_rad_s * radius / gas.sound_speed,
            result.bandwidth_quality_factor,
            result.evaluations
        );
        results.push(result);
    }
    assert!(
        (results[0].angular_frequency_rad_s - results[1].angular_frequency_rad_s).abs()
            > 1.0e-3 * results[0].angular_frequency_rad_s
    );
    assert!(
        (results[0].bandwidth_quality_factor - results[1].bandwidth_quality_factor).abs() > 1.0e-3
    );
}
