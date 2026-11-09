//! Original-source controls for coherent spatial frequency derivatives.
use super::super::*;
use crate::plane_wave::{
    dispersion::IncidenceFrequencyConstraint,
    uniaxial::oriented::{
        OrientedCrystalField, OrientedCrystalFieldSide, OrientedUniaxialInterface,
    },
};

fn components(field: OrientedCrystalField) -> [[C64; 3]; 2] {
    [field.electric_v_m, field.magnetic_a_m]
}

#[test]
fn spatial_jets_follow_original_frequency_sources_and_both_incidence_constraints() {
    use OrientedCrystalFieldSide::*;
    let outside = incident(2.4, MU);
    let (ordinary, extraordinary) = axes(
        [1.4, 1.6],
        [1.8, 2.1],
        [0.0; 2],
        Some(1.2 * MU),
        [400.0, 800.0],
        2,
    );
    let source = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap(),
        [1.0, 2.0, -0.7],
    )
    .unwrap();
    let omega = frequency(550.0);
    let electric = [C64::new(0.6, -0.2), C64::new(0.1, 0.7)];
    // The first input is signed and propagating; the second exceeds both
    // principal indices and has two decaying transmitted Maxwell modes.
    for p in [-0.4, 2.25] {
        let prepared = source
            .prepare_field_at_tangential_index(omega, p, electric)
            .unwrap();
        for mode in &prepared.boundary_response().transmitted {
            assert_eq!(
                mode.wave_kind,
                if p < 0.0 {
                    UniaxialWaveKind::Propagating
                } else {
                    UniaxialWaveKind::Evanescent
                }
            );
        }
        for constraint in [
            IncidenceFrequencyConstraint::FixedAngle,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
        ] {
            let boundary = prepared.boundary_frequency_dispersion(constraint).unwrap();
            assert!(core::ptr::eq(boundary.source(), &source));
            for side in [
                Incident,
                Reflected,
                TotalIncidentSide,
                Ordinary,
                Extraordinary,
                Transmitted,
            ] {
                let crystal_side = matches!(side, Ordinary | Extraordinary | Transmitted);
                let point = [0.17e-6, -3.0, if crystal_side { 0.23e-6 } else { -0.23e-6 }];
                let native = prepared.field(side, point).unwrap();
                let actual = boundary.field_frequency_dispersion(side, point).unwrap();
                assert_eq!(actual.carrier_field(), native);
                let jet = actual.vector_frequency_jet();
                assert_eq!(jet.electric_v_m(), native.electric_v_m);
                assert_eq!(jet.magnetic_a_m(), native.magnetic_a_m);

                // Independent frequency samples use the original material
                // tables and the original coherent field solver. This stencil
                // differentiates u=(omega_sample-omega)/omega, so its two
                // outputs are exactly the public scaled ordinary derivatives.
                let h = 1.0 / 1024.0;
                let values = [-2.0, -1.0, 0.0, 1.0, 2.0].map(|step| {
                    let ratio = 1.0 + step * h;
                    let sampled_p = match constraint {
                        // The original incident material is nondispersive.
                        IncidenceFrequencyConstraint::FixedAngle => p,
                        IncidenceFrequencyConstraint::FixedTransverseWavenumber => p / ratio,
                    };
                    components(
                        source
                            .prepare_field_at_tangential_index(omega * ratio, sampled_p, electric)
                            .unwrap()
                            .field(side, point)
                            .unwrap(),
                    )
                });
                assert_eq!(values[2], components(native));
                let first = [
                    jet.scaled_electric_derivatives_v_m()[0],
                    jet.scaled_magnetic_derivatives_a_m()[0],
                ];
                let second = [
                    jet.scaled_electric_derivatives_v_m()[1],
                    jet.scaled_magnetic_derivatives_a_m()[1],
                ];
                let physical_first = [
                    jet.electric_derivative_v_s_m(),
                    jet.magnetic_derivative_a_s_m(),
                ];
                let physical_second = [
                    jet.electric_second_derivative_v_s2_m(),
                    jet.magnetic_second_derivative_a_s2_m(),
                ];
                for family in 0..2 {
                    for component in 0..3 {
                        let y = values.map(|value| value[family][component]);
                        let expected_first = (y[0] - y[1].scale(8.0) + y[3].scale(8.0) - y[4])
                            .scale(1.0 / (12.0 * h));
                        let expected_second =
                            (-y[4] + y[3].scale(16.0) - y[2].scale(30.0) + y[1].scale(16.0) - y[0])
                                .scale(1.0 / (12.0 * h * h));
                        let scale = y[2]
                            .abs()
                            .max(expected_first.abs())
                            .max(expected_second.abs())
                            .max(1e-5);
                        complex_near(first[family][component], expected_first, 3e-8 * scale);
                        complex_near(second[family][component], expected_second, 3e-7 * scale);
                        complex_near(
                            physical_first[family][component].scale(omega),
                            first[family][component],
                            2e-14 * scale,
                        );
                        complex_near(
                            physical_second[family][component].scale(omega).scale(omega),
                            second[family][component],
                            2e-14 * scale,
                        );
                    }
                }
                // Point queries retain the prepared native carrier and its
                // coherent Poynting calculation, including modal interference.
                assert_eq!(prepared.field(side, point).unwrap(), native);
            }
        }
    }
}

#[test]
fn spatial_jets_retain_deep_derivative_tails_and_zero_carrier_frequency_fields() {
    use OrientedCrystalFieldSide::*;
    let outside = incident(2.0, MU);
    let (ordinary, extraordinary) = axes([1.0; 2], [1.7; 2], [0.0; 2], Some(MU), [400.0, 800.0], 0);
    let source = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap(),
        [1.0, 2.0, 0.7],
    )
    .unwrap();
    let omega = frequency(550.0);
    for (amplitude, target_exposure) in [(1.0, 2.0), (1.0, 750.0), (1e100, 800.0)] {
        let electric = [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)].map(|z| z.scale(amplitude));
        let prepared = source.prepare_field(omega, 0.1, electric).unwrap();
        let boundary = prepared
            .boundary_frequency_dispersion(IncidenceFrequencyConstraint::FixedAngle)
            .unwrap();
        let mode = prepared.boundary_response().transmitted[0];
        assert_eq!(mode.wave_kind, UniaxialWaveKind::Evanescent);
        assert_eq!(mode.normal_wavenumber_rad_m.re, 0.0);
        let depth = target_exposure / -mode.normal_wavenumber_rad_m.im;
        let exposure = -mode.normal_wavenumber_rad_m.im * depth;
        let point = [0.0, 0.0, depth];
        let native = prepared.field(Ordinary, point).unwrap();
        let actual = boundary
            .field_frequency_dispersion(Ordinary, point)
            .unwrap();
        assert_eq!(actual.carrier_field(), native);
        let jet = actual.vector_frequency_jet();
        let derivatives = [
            jet.scaled_electric_derivatives_v_m(),
            jet.scaled_magnetic_derivatives_a_m(),
        ];
        let original = [mode.electric_v_m, mode.magnetic_a_m];
        // Every index and boundary amplitude is nondispersive. At fixed angle
        // the ordinary phase is -i*exposure, hence F1=-exposure*F and
        // F2=exposure^2*F. Evaluate the independent reference in log magnitude
        // so exp(-exposure) need not be representable on its own.
        for family in 0..2 {
            for order in 0..2 {
                let multiplier = if order == 0 {
                    -exposure
                } else {
                    exposure * exposure
                };
                for component in 0..3 {
                    let a = original[family][component];
                    let weighted = [a.re * multiplier, a.im * multiplier];
                    let expected = weighted.map(|value| {
                        if value == 0.0 {
                            0.0
                        } else {
                            value.signum() * (value.abs().ln() - exposure).exp()
                        }
                    });
                    let got = derivatives[family][order][component];
                    for (got, expected) in [(got.re, expected[0]), (got.im, expected[1])] {
                        near(
                            got,
                            expected,
                            3e-11 * expected.abs() + 8.0 * f64::from_bits(1),
                        );
                    }
                }
            }
        }
        if target_exposure == 750.0 {
            assert_eq!(native.electric_v_m, [C64::ZERO; 3]);
            assert_eq!(native.magnetic_a_m, [C64::ZERO; 3]);
            assert!(
                jet.scaled_electric_derivatives_v_m()[1]
                    .iter()
                    .any(|value| *value != C64::ZERO)
            );
            assert_eq!(jet.electric_component_logarithmic_dispersion(), [None; 3]);
        }
        if amplitude == 1e100 {
            assert!(native.electric_v_m.iter().any(|value| *value != C64::ZERO));
        }
        assert!(
            boundary
                .field_frequency_dispersion(Ordinary, [0.0, 0.0, -1.0])
                .is_err()
        );
        assert!(
            boundary
                .field_frequency_dispersion(Transmitted, [f64::NAN, 0.0, 0.0])
                .is_err()
        );
    }

    // A matched carrier can have zero reflection but a nonzero frequency
    // derivative because the original crystal source has a different slope.
    // Such a wave must remain active in the coherent frequency accumulator.
    let outside = incident(1.5, MU);
    let omega = frequency(600.0);
    let epsilon = DielectricResponse::Lossless(&outside)
        .lossless_permittivity(omega)
        .unwrap();
    let ni = (C0 * fs_math::det::sqrt(MU)) * fs_math::det::sqrt(epsilon);
    let (ordinary, extraordinary) = axes(
        [ni - 0.25, ni + 0.25],
        [ni - 0.25, ni + 0.25],
        [0.0; 2],
        Some(MU),
        [400.0, 800.0],
        0,
    );
    let crystal = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap();
    assert_eq!(crystal.phase_indices(omega).unwrap(), [ni; 2]);
    let source = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [1.0, 2.0, 3.0],
    )
    .unwrap();
    let prepared = source
        .prepare_field(omega, 1.0, [C64::ZERO, C64::ONE])
        .unwrap();
    let boundary = prepared
        .boundary_frequency_dispersion(IncidenceFrequencyConstraint::FixedAngle)
        .unwrap();
    assert_eq!(boundary.reflected_field().electric_v_m(), [C64::ZERO; 3]);
    let point = [0.17e-6, -3.0, -0.23e-6];
    let actual = boundary
        .field_frequency_dispersion(Reflected, point)
        .unwrap();
    assert_eq!(
        actual.carrier_field(),
        prepared.field(Reflected, point).unwrap()
    );
    let jet = actual.vector_frequency_jet();
    assert_eq!(jet.electric_v_m(), [C64::ZERO; 3]);
    assert_eq!(jet.electric_component_logarithmic_dispersion(), [None; 3]);
    let phase = -prepared
        .boundary_response()
        .incident_normal_wavenumber_rad_m
        * point[2];
    let (sine, cosine) = phase.sin_cos();
    let rotation = C64::new(cosine, -sine);
    let a = boundary.reflected_field().scaled_electric_derivatives_v_m();
    assert!(a[0][1].abs() > 1e-4);
    for component in 0..3 {
        complex_near(
            jet.scaled_electric_derivatives_v_m()[0][component],
            a[0][component] * rotation,
            2e-12,
        );
        complex_near(
            jet.scaled_electric_derivatives_v_m()[1][component],
            (a[1][component] + a[0][component] * C64::new(0.0, -2.0 * phase)) * rotation,
            2e-12,
        );
    }
    assert!(
        boundary
            .field_frequency_dispersion(Reflected, [0.0, 0.0, 1.0])
            .is_err()
    );
    // Retaining a point jet does not expand either original material domain.
    for wavelength in [400.0, 800.0] {
        let prepared = source
            .prepare_field(frequency(wavelength), 1.0, [C64::ZERO, C64::ONE])
            .unwrap();
        assert!(prepared.field(Transmitted, [0.0, 0.0, 1e-7]).is_ok());
        assert!(
            prepared
                .boundary_frequency_dispersion(IncidenceFrequencyConstraint::FixedAngle)
                .is_err()
        );
    }
}
