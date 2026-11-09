use super::super::tests::{medium, near, near_complex};
use super::*;

#[test]
fn g1_coherent_observer_reduces_to_each_direction_with_signed_power() {
    let prepared = medium(0.04, 0.003, 0.002)
        .prepare_wave_absorption(100.0)
        .unwrap();
    let drive = C64::new(0.2, -0.3);
    let length = 0.07;
    for reverse in [false, true] {
        let (f, b) = if reverse {
            (C64::ZERO, drive)
        } else {
            (drive, C64::ZERO)
        };
        let segment = prepared.coherent_segment(length, f, b).unwrap();
        assert_eq!(segment.absorption(), &prepared);
        assert_eq!(segment.length_m(), length);
        let sign = if reverse { -1.0 } else { 1.0 };
        for x in [0.0, 0.013, length] {
            let distance = if reverse { length - x } else { x };
            let reference = prepared.at_distance(drive, distance).unwrap();
            let actual = segment.at_distance(x).unwrap();
            for (a, r) in [
                (actual.velocity_m_s, reference.velocity_m_s),
                (actual.pressure_pa, reference.pressure_pa.scale(sign)),
                (
                    actual.normal_traction_pa,
                    reference.normal_traction_pa.scale(sign),
                ),
                (
                    actual.temperature_perturbation_k,
                    reference.temperature_perturbation_k.scale(sign),
                ),
                (
                    actual.conductive_heat_flux_w_m2,
                    reference.conductive_heat_flux_w_m2,
                ),
            ] {
                near_complex(a, r, 2e-14);
            }
            near(
                actual.acoustic_power_flux_w_m2(),
                sign * reference.acoustic_power_flux_w_m2(),
                2e-14,
            );
            near(
                actual.viscous_heating_w_m3,
                reference.viscous_heating_w_m3,
                2e-14,
            );
            near(
                actual.conductive_entropy_w_m3_k,
                reference.conductive_entropy_w_m3_k,
                2e-14,
            );
        }
        let actual = segment.absorption_between(0.0, length, 0.7).unwrap();
        let reference = prepared
            .absorption_between(drive, 0.0, length, 0.7)
            .unwrap();
        near(actual.viscous_heat_w, reference.viscous_heat_w, 2e-14);
        near(
            actual.conductive_entropy_w_k,
            reference.conductive_entropy_w_k,
            2e-14,
        );
        near(
            actual.left_acoustic_power_w - actual.right_acoustic_power_w,
            actual.acoustic_power_loss_w,
            3e-12,
        );
        if reverse {
            assert!(actual.left_acoustic_power_w < 0.0 && actual.right_acoustic_power_w < 0.0);
        }
    }
}

#[test]
fn g1_standing_wave_nodes_separate_viscous_and_conductive_loss() {
    let prepared = medium(0.08, 0.003, 0.002)
        .prepare_wave_absorption(100.0)
        .unwrap();
    let length = 0.02;
    let pressure_node = prepared
        .coherent_segment(length, C64::ONE, C64::ONE)
        .unwrap();
    let a = pressure_node.at_distance(0.5 * length).unwrap();
    assert_eq!(a.pressure_pa, C64::ZERO);
    assert_eq!(a.temperature_perturbation_k, C64::ZERO);
    assert_eq!(a.viscous_heating_w_m3, 0.0);
    assert!(a.conductive_entropy_w_m3_k > 0.0);
    assert_eq!(a.acoustic_power_flux_w_m2(), 0.0);
    let velocity_node = prepared
        .coherent_segment(length, C64::ONE, -C64::ONE)
        .unwrap();
    let b = velocity_node.at_distance(0.5 * length).unwrap();
    assert_eq!(b.velocity_m_s, C64::ZERO);
    assert_eq!(b.conductive_heat_flux_w_m2, C64::ZERO);
    assert_eq!(b.conductive_entropy_w_m3_k, 0.0);
    assert!(b.viscous_heating_w_m3 > 0.0);
    assert_eq!(b.acoustic_power_flux_w_m2(), 0.0);
    // Equal incident powers with different relative phase have different
    // spatially integrated physical losses; incoherent power addition fails.
    let ca = pressure_node.absorption_between(0.0, length, 1.0).unwrap();
    let cb = velocity_node.absorption_between(0.0, length, 1.0).unwrap();
    assert!(ca.conductive_entropy_w_k > 100.0 * cb.conductive_entropy_w_k);
    assert!(cb.viscous_heat_w > 100.0 * ca.viscous_heat_w);
}

#[test]
fn g1_coherent_losses_match_original_field_quadrature_and_flux_balance() {
    // Independent heat-equation temperature coefficient and Simpson integration
    // of the original positive constitutive forms, not production Gram weights.
    for expansion in [-0.08, 0.0, 0.04] {
        for omega in [100.0, 10_000.0] {
            let equilibrium = medium(expansion, 0.003, 0.002);
            let prepared = equilibrium.prepare_wave_absorption(omega).unwrap();
            let k = prepared.wave().wavenumber_rad_m();
            let ik = C64::new(-k.im, k.re);
            let cv = equilibrium.absorption().state().specific_heat_cv_j_kg_k();
            let temperature = -ik
                .scale(expansion * 12.5 / equilibrium.isothermal_compressibility_pa_inv())
                / ((k * k).scale(0.6) + C64::new(0.0, -omega * 2.0 * cv));
            // Exercise the small-cell series and both signs of sinc(k_re L).
            for phase_span in [0.2, 1.2, 4.4] {
                let length = phase_span / k.re;
                let (front, back, area) = (C64::new(0.4, 0.2), C64::new(-0.3, 0.5), 0.7);
                let segment = prepared.coherent_segment(length, front, back).unwrap();
                let cell = segment.absorption_between(0.0, length, area).unwrap();
                let phasor = |distance: f64| {
                    let (s, c) = det::sin_cos(k.re * distance);
                    C64::new(c, s).scale(det::exp(-k.im * distance))
                };
                let mut viscous = 0.0;
                let mut entropy = 0.0;
                const STEPS: u32 = 2048;
                for i in 0..=STEPS {
                    let x = length * f64::from(i) / f64::from(STEPS);
                    let f = front * phasor(x);
                    let b = back * phasor(length - x);
                    let weight = if i == 0 || i == STEPS {
                        1.0
                    } else if i % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    viscous += weight * 0.003 * (ik * (f - b)).norm_sq();
                    entropy +=
                        weight * 0.6 / (2.0 * 12.5 * 12.5) * (ik * temperature * (f + b)).norm_sq();
                }
                let scale = area * length / (3.0 * f64::from(STEPS));
                near(cell.viscous_heat_w, viscous * scale, 4e-12);
                if expansion == 0.0 {
                    assert_eq!(cell.conductive_entropy_w_k, 0.0);
                } else {
                    near(cell.conductive_entropy_w_k, entropy * scale, 4e-12);
                }
                near(
                    cell.acoustic_power_loss_w,
                    cell.viscous_heat_w + 12.5 * cell.conductive_entropy_w_k,
                    3e-15,
                );
                assert!(
                    cell.balance_residual_w.abs()
                        < 5e-13
                            * (cell.left_acoustic_power_w.abs()
                                + cell.right_acoustic_power_w.abs()
                                + cell.acoustic_power_loss_w)
                );
            }
        }
    }
}

#[test]
fn g1_thin_node_integrals_retain_the_positive_cubic_length_limit() {
    let prepared = medium(0.08, 0.003, 0.002)
        .prepare_wave_absorption(100.0)
        .unwrap();
    let k = prepared.wave().wavenumber_rad_m();
    let length = 1e-12 / k.re;
    // Equal endpoint amplitudes: uf-ub ~ ik(2x-L). Integrated viscosity
    // is eta |k|^4 L^3/6, well below direct d-c subtraction precision.
    let a = prepared
        .coherent_segment(length, C64::ONE, C64::ONE)
        .unwrap()
        .absorption_between(0.0, length, 1.0)
        .unwrap();
    let expected = 0.006 * k.norm_sq().powi(2) * length.powi(3) / 6.0;
    assert!(a.viscous_heat_w > 0.0);
    near(a.viscous_heat_w, expected, 3e-12);
    let b = prepared
        .coherent_segment(length, C64::ONE, -C64::ONE)
        .unwrap()
        .absorption_between(0.0, length, 1.0)
        .unwrap();
    let temperature = prepared
        .at_distance(C64::ONE, 0.0)
        .unwrap()
        .temperature_perturbation_k;
    let thermal_expected =
        (0.6 / (6.0 * 12.5)) * temperature.norm_sq() * k.norm_sq().powi(2) * length.powi(3);
    assert!(b.thermal_availability_loss_w > 0.0);
    near(b.thermal_availability_loss_w, thermal_expected, 3e-12);
}

#[test]
fn g3_coherent_partition_phase_and_direction_reversal_preserve_loss() {
    let prepared = medium(-0.04, 0.003, 0.002)
        .prepare_wave_absorption(100.0)
        .unwrap();
    let length = 0.3;
    let (f, b) = (C64::new(0.4, 0.2), C64::new(-0.3, 0.5));
    let segment = prepared.coherent_segment(length, f, b).unwrap();
    let full = segment.absorption_between(0.0, length, 0.7).unwrap();
    let mut sums = [0.0; 3];
    for i in 0..31 {
        let cell = segment
            .absorption_between(
                length * f64::from(i) / 31.0,
                length * f64::from(i + 1) / 31.0,
                0.7,
            )
            .unwrap();
        sums[0] += cell.viscous_heat_w;
        sums[1] += cell.conductive_entropy_w_k;
        sums[2] += cell.acoustic_power_loss_w;
    }
    for (actual, expected) in sums.into_iter().zip([
        full.viscous_heat_w,
        full.conductive_entropy_w_k,
        full.acoustic_power_loss_w,
    ]) {
        near(actual, expected, 2e-14);
    }
    let phase = C64::new(0.6, 0.8);
    let shifted = prepared
        .coherent_segment(length, phase * f, phase * b)
        .unwrap();
    let shifted_cell = shifted.absorption_between(0.0, length, 0.7).unwrap();
    near(
        shifted_cell.acoustic_power_loss_w,
        full.acoustic_power_loss_w,
        2e-14,
    );
    let reversed = prepared.coherent_segment(length, -b, -f).unwrap();
    for x in [0.0, 0.07, length] {
        let original = segment.at_distance(length - x).unwrap();
        let opposite = reversed.at_distance(x).unwrap();
        near_complex(opposite.velocity_m_s, -original.velocity_m_s, 2e-14);
        near_complex(opposite.pressure_pa, original.pressure_pa, 2e-14);
        near_complex(
            opposite.conductive_heat_flux_w_m2,
            -original.conductive_heat_flux_w_m2,
            2e-14,
        );
        near(
            opposite.acoustic_power_flux_w_m2(),
            -original.acoustic_power_flux_w_m2(),
            2e-14,
        );
    }
    let mirrored = reversed.absorption_between(0.0, length, 0.7).unwrap();
    near(
        mirrored.acoustic_power_loss_w,
        full.acoustic_power_loss_w,
        2e-14,
    );
}

#[test]
fn g1_opaque_two_sided_segment_uses_only_decaying_modes() {
    let prepared = medium(0.08, 0.003, 0.002)
        .prepare_wave_absorption(10_000.0)
        .unwrap();
    let length = f64::MAX;
    let (f, b) = (C64::new(0.4, 0.2), C64::new(-0.3, 0.5));
    let segment = prepared.coherent_segment(length, f, b).unwrap();
    assert_eq!(
        segment.at_distance(0.5 * length).unwrap().velocity_m_s,
        C64::ZERO
    );
    let both = segment.absorption_between(0.0, length, 0.7).unwrap();
    let first = prepared.absorption_between(f, 0.0, length, 0.7).unwrap();
    let second = prepared.absorption_between(b, 0.0, length, 0.7).unwrap();
    near(
        both.acoustic_power_loss_w,
        first.acoustic_power_loss_w + second.acoustic_power_loss_w,
        2e-14,
    );
    near(
        both.left_acoustic_power_w,
        first.incoming_acoustic_power_w,
        2e-14,
    );
    near(
        both.right_acoustic_power_w,
        -second.incoming_acoustic_power_w,
        2e-14,
    );
    assert!(both.balance_residual_w.abs() < 2e-14 * both.acoustic_power_loss_w);
}

#[test]
fn g0_coherent_segment_refuses_invalid_domains_and_preserves_zero_cells() {
    let prepared = medium(0.04, 0.003, 0.002)
        .prepare_wave_absorption(100.0)
        .unwrap();
    for length in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            prepared
                .coherent_segment(length, C64::ONE, C64::ONE)
                .is_err()
        );
    }
    for bad in [C64::new(f64::NAN, 0.0), C64::new(0.0, f64::INFINITY)] {
        assert!(prepared.coherent_segment(1.0, bad, C64::ONE).is_err());
        assert!(prepared.coherent_segment(1.0, C64::ONE, bad).is_err());
    }
    let segment = prepared.coherent_segment(1.0, C64::ONE, C64::ONE).unwrap();
    for x in [-0.01, 1.01, f64::NAN, f64::INFINITY] {
        assert!(segment.at_distance(x).is_err());
    }
    for (a, b, area) in [
        (0.1, 0.0, 1.0),
        (0.0, 1.1, 1.0),
        (0.0, 1.0, 0.0),
        (0.0, 1.0, f64::NAN),
        (f64::NAN, 1.0, 1.0),
    ] {
        assert!(segment.absorption_between(a, b, area).is_err());
    }
    let empty = segment.absorption_between(0.1, 0.1, 0.7).unwrap();
    assert_eq!(empty.acoustic_power_loss_w, 0.0);
    assert_eq!(empty.left_acoustic_power_w, empty.right_acoustic_power_w);
    assert_eq!(empty.balance_residual_w, 0.0);
    let zero = prepared
        .coherent_segment(0.0, C64::ZERO, C64::ZERO)
        .unwrap();
    assert_eq!(
        zero.at_distance(0.0).unwrap().acoustic_power_flux_w_m2(),
        0.0
    );
    assert_eq!(
        zero.absorption_between(0.0, 0.0, 1.0)
            .unwrap()
            .acoustic_power_loss_w,
        0.0
    );
}
