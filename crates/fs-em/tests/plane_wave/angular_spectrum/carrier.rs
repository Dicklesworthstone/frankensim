//! Arbitrary-angle Snell/Fresnel and vector-field controls on the native beam.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumCoherence, AngularSpectrumCoherenceLimits, AngularSpectrumExposure,
    AngularSpectrumExposureLimits, AngularSpectrumIlluminationLimits, AngularSpectrumLayerLimits,
};

fn grid() -> AngularSpectrumGrid {
    AngularSpectrumGrid {
        width: 4,
        height: 4,
        spacing_m: [40.0, 35.0],
    }
}

#[test]
fn arbitrary_carrier_refracts_with_exact_vector_fields_and_fresnel_power() {
    let air = material(1.0, 1.0);
    let glass = material(2.25, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&air, DielectricResponse::Lossless(&glass), 1.0, &[])
            .unwrap();
    // This direction is not any integer mode of the authored periodic cell.
    let carrier = [0.36, -0.27];
    let p = [0.8, -0.6];
    let s = [0.6, 0.8];
    let qi = (1.0_f64 - 0.45 * 0.45).sqrt();
    let qt = (2.25_f64 - 0.45 * 0.45).sqrt();
    let area = 16.0 * 40.0 * 35.0;
    for (polarization, axis, reflectance, admittance) in [
        (Polarization::Te, s, (qi - qt) / (qi + qt), qi),
        (
            Polarization::Tm,
            p,
            (qt - 2.25 * qi) / (qt + 2.25 * qi),
            1.0 / qi,
        ),
    ] {
        let input = vec![[C64::from_re(axis[0]), C64::from_re(axis[1])]; 16];
        let beam = PlanarAngularSpectrum::prepare_with_carrier(
            &stack,
            grid(),
            carrier,
            &input,
            limits(),
            || false,
        )
        .unwrap();
        assert_eq!(beam.carrier_transverse_wavenumber_rad_m(), carrier);
        relative(beam.power().incident_w, 0.5 * area * admittance);
        relative(
            beam.power().reflected_w,
            0.5 * area * admittance * reflectance.powi(2),
        );
        relative(
            beam.power().transmitted_w,
            0.5 * area * admittance * (1.0 - reflectance.powi(2)),
        );
        let field = beam.field(Side::Transmitted, 0.73, || false).unwrap();
        let t = 1.0 + reflectance;
        let e = [
            axis[0] * t,
            axis[1] * t,
            -(carrier[0] * axis[0] + carrier[1] * axis[1]) * t / qt,
        ];
        let h = [
            carrier[1] * e[2] - qt * e[1],
            qt * e[0] - carrier[0] * e[2],
            carrier[0] * e[1] - carrier[1] * e[0],
        ];
        for row in 0..4 {
            for column in 0..4 {
                let expected_phase = phase(
                    -carrier[0] * column as f64 * 40.0 - carrier[1] * row as f64 * 35.0 - qt * 0.73,
                );
                let index = row * 4 + column;
                for c in 0..3 {
                    complex_near(
                        field.electric_v_per_m()[index][c],
                        expected_phase.scale(e[c]),
                        2e-13,
                    );
                    complex_near(
                        field.magnetic_a_per_m()[index][c],
                        expected_phase.scale(h[c]),
                        2e-13,
                    );
                }
            }
        }
        relative(field.normal_power_w(), beam.power().transmitted_w);
        let response = stack.response(qi, polarization).unwrap();
        relative(response.exit_normal_wavenumber_rad_m.re, qt);
    }
}

#[test]
fn modulated_carrier_matches_independent_coherent_snell_fresnel_superposition() {
    let air = material(1.0, 1.0);
    let glass = material(2.25, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&air, DielectricResponse::Lossless(&glass), 1.0, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 8,
        height: 4,
        spacing_m: [10.0, 11.0],
    };
    let carrier = [0.34, -0.23];
    // Two signed, nonzero envelope bins with independent complex Ex/Ey.
    // Neither physical transverse wavevector lies on the unshifted FFT grid.
    let envelope_modes = [
        ([1.0, -1.0], [C64::new(0.6, -0.2), C64::new(0.35, 0.5)]),
        ([-2.0, 1.0], [C64::new(-0.4, 0.7), C64::new(0.2, -0.3)]),
    ]
    .map(|(bin, amplitude)| {
        (
            [
                core::f64::consts::TAU * bin[0] / (grid.width as f64 * grid.spacing_m[0]),
                core::f64::consts::TAU * bin[1] / (grid.height as f64 * grid.spacing_m[1]),
            ],
            amplitude,
        )
    });
    let mut input = vec![[C64::ZERO; 2]; grid.width * grid.height];
    for row in 0..grid.height {
        for column in 0..grid.width {
            let x = column as f64 * grid.spacing_m[0];
            let y = row as f64 * grid.spacing_m[1];
            for (kt, amplitude) in envelope_modes {
                let weight = phase(-(kt[0] * x + kt[1] * y));
                for c in 0..2 {
                    let sample = &mut input[row * grid.width + column][c];
                    *sample = *sample + amplitude[c] * weight;
                }
            }
        }
    }
    let beam = PlanarAngularSpectrum::prepare_with_carrier(
        &stack,
        grid,
        carrier,
        &input,
        limits(),
        || false,
    )
    .unwrap();
    // Independent half-space Fresnel amplitudes for the tangential electric
    // field, followed by k.E=0 and H=k cross E/(omega*mu). No native response,
    // modal decomposition or FFT is used by this reference superposition.
    let transmitted_modes = envelope_modes.map(|(offset, amplitude)| {
        let [kx, ky] = [carrier[0] + offset[0], carrier[1] + offset[1]];
        let kt_squared = kx * kx + ky * ky;
        let kt = kt_squared.sqrt();
        let qi = (1.0 - kt_squared).sqrt();
        let qt = (2.25 - kt_squared).sqrt();
        let p = [kx / kt, ky / kt];
        let s = [-p[1], p[0]];
        let te = (amplitude[0].scale(s[0]) + amplitude[1].scale(s[1])).scale(2.0 * qi / (qi + qt));
        let tm = (amplitude[0].scale(p[0]) + amplitude[1].scale(p[1]))
            .scale(2.0 * qt / (qt + 2.25 * qi));
        let ex = te.scale(s[0]) + tm.scale(p[0]);
        let ey = te.scale(s[1]) + tm.scale(p[1]);
        let ez = -(ex.scale(kx) + ey.scale(ky)).scale(1.0 / qt);
        let e = [ex, ey, ez];
        let h = [
            ez.scale(ky) - ey.scale(qt),
            ex.scale(qt) - ez.scale(kx),
            ey.scale(kx) - ex.scale(ky),
        ];
        ([kx, ky], qt, e, h)
    });
    for distance in [0.0, 0.73, 7.5] {
        let field = beam.field(Side::Transmitted, distance, || false).unwrap();
        for row in 0..grid.height {
            for column in 0..grid.width {
                let x = column as f64 * grid.spacing_m[0];
                let y = row as f64 * grid.spacing_m[1];
                let mut expected_e = [C64::ZERO; 3];
                let mut expected_h = [C64::ZERO; 3];
                for (kt, qt, e, h) in transmitted_modes {
                    let weight = phase(-(kt[0] * x + kt[1] * y + qt * distance));
                    for c in 0..3 {
                        expected_e[c] = expected_e[c] + e[c] * weight;
                        expected_h[c] = expected_h[c] + h[c] * weight;
                    }
                }
                let index = row * grid.width + column;
                for c in 0..3 {
                    complex_near(field.electric_v_per_m()[index][c], expected_e[c], 2e-13);
                    complex_near(field.magnetic_a_per_m()[index][c], expected_h[c], 2e-13);
                }
            }
        }
    }

    let legacy = PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
    let zero = PlanarAngularSpectrum::prepare_with_carrier(
        &stack,
        grid,
        [0.0; 2],
        &input,
        limits(),
        || false,
    )
    .unwrap();
    let power_bits = |beam: &PlanarAngularSpectrum<'_, '_>| {
        let power = beam.power();
        [
            power.incident_w,
            power.reflected_w,
            power.transmitted_w,
            power.film_absorbed_w,
        ]
        .map(f64::to_bits)
    };
    assert_eq!(power_bits(&zero), power_bits(&legacy));
    let complex_bits = |z: C64| [z.re.to_bits(), z.im.to_bits()];
    for side in [
        Side::Incident,
        Side::Reflected,
        Side::TotalIncidentSide,
        Side::Transmitted,
    ] {
        let a = zero.field(side, 0.73, || false).unwrap();
        let b = legacy.field(side, 0.73, || false).unwrap();
        assert_eq!(a.normal_power_w().to_bits(), b.normal_power_w().to_bits());
        for (left, right) in a
            .electric_v_per_m()
            .iter()
            .chain(a.magnetic_a_per_m())
            .zip(b.electric_v_per_m().iter().chain(b.magnetic_a_per_m()))
        {
            assert_eq!(left.map(complex_bits), right.map(complex_bits));
        }
    }
}

#[test]
fn carrier_resolves_brewster_and_total_internal_refraction_without_angle_quantization() {
    let low = material(1.0, 1.0);
    let high = material(2.25, 1.0);
    let brewster =
        PlanarDielectricLayers::try_new(&low, DielectricResponse::Lossless(&high), 1.0, &[])
            .unwrap();
    let input = vec![[C64::ONE, C64::ZERO]; 16];
    let beam = PlanarAngularSpectrum::prepare_with_carrier(
        &brewster,
        grid(),
        [1.5 / 3.25_f64.sqrt(), 0.0],
        &input,
        limits(),
        || false,
    )
    .unwrap();
    assert!(beam.power().reflected_w < 1e-28 * beam.power().incident_w);

    let tir = PlanarDielectricLayers::try_new(&high, DielectricResponse::Lossless(&low), 1.0, &[])
        .unwrap();
    let beam = PlanarAngularSpectrum::prepare_with_carrier(
        &tir,
        grid(),
        [1.2, 0.0],
        &input,
        limits(),
        || false,
    )
    .unwrap();
    assert_eq!(beam.power().transmitted_w, 0.0);
    relative(beam.power().reflected_w, beam.power().incident_w);
    let a = beam.field(Side::Transmitted, 0.0, || false).unwrap();
    let b = beam.field(Side::Transmitted, 0.8, || false).unwrap();
    let decay = (-0.8 * 0.44_f64.sqrt()).exp();
    for (left, right) in a.electric_v_per_m().iter().zip(b.electric_v_per_m()) {
        for c in 0..3 {
            complex_near(right[c], left[c].scale(decay), 2e-13);
        }
    }
}

#[test]
fn carrier_reaches_internal_fields_coherence_and_both_illumination_faces() {
    let medium = material(1.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&medium),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&medium),
            thickness_m: 0.6,
        }],
    )
    .unwrap();
    let input = vec![[C64::ZERO, C64::ONE]; 16];
    let beam = PlanarAngularSpectrum::prepare_with_carrier(
        &stack,
        grid(),
        [0.4, 0.0],
        &input,
        limits(),
        || false,
    )
    .unwrap();
    let layer = beam
        .prepare_layer(
            0,
            AngularSpectrumLayerLimits {
                maximum_layer_mode_pairs: 128,
                maximum_workspace_bytes: 4 << 20,
                maximum_fft_work: 1 << 20,
            },
            || false,
        )
        .unwrap();
    let internal = layer.field(0.3, || false).unwrap();
    let q = 0.84_f64.sqrt();
    for row in 0..4 {
        for column in 0..4 {
            complex_near(
                internal.field().electric_v_per_m()[row * 4 + column][1],
                phase(-0.4 * column as f64 * 40.0 - q * 0.3),
                2e-13,
            );
        }
    }
    let illumination = beam
        .prepare_illumination(
            &input,
            AngularSpectrumIlluminationLimits {
                maximum_layer_mode_pairs: 1024,
                maximum_workspace_bytes: 4 << 20,
                maximum_fft_work: 1 << 20,
            },
            || false,
        )
        .unwrap();
    relative(illumination.power().incoming_w[0], beam.power().incident_w);
    relative(illumination.power().incoming_w[1], beam.power().incident_w);
    let field = illumination.layer_field(0, 0.3, || false).unwrap();
    for (both, left) in field
        .field()
        .electric_v_per_m()
        .iter()
        .zip(internal.field().electric_v_per_m())
    {
        complex_near(both[1], left[1].scale(2.0), 2e-13);
    }
    let factors = [beam];
    let bank = AngularSpectrumCoherence::try_new(
        &factors[0],
        &factors,
        AngularSpectrumCoherenceLimits {
            maximum_modes: 2,
            maximum_point_mode_pairs: 128,
            maximum_workspace_bytes: 4 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap();
    relative(bank.power().incident_w, factors[0].power().incident_w);
}

#[test]
fn carrier_admission_and_exposure_boundaries_are_explicit() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let input = vec![[C64::ZERO, C64::ONE]; 16];
    let prepare = |carrier, budget| {
        PlanarAngularSpectrum::prepare_with_carrier(&stack, grid(), carrier, &input, budget, || {
            false
        })
    };
    let beam = prepare([0.4, 0.0], limits()).unwrap();
    assert!(
        prepare(
            [0.4, 0.0],
            AngularSpectrumLimits {
                maximum_workspace_bytes: beam.workspace_bytes() - 1,
                ..limits()
            }
        )
        .is_err()
    );
    assert!(matches!(
        prepare([0.99, 0.0], limits()),
        Err(AngularSpectrumError::IncidentLightCone { .. })
    ));
    assert!(prepare([f64::NAN, 0.0], limits()).is_err());
    let exposure_limits = AngularSpectrumExposureLimits {
        maximum_lines: 2,
        maximum_point_line_pairs: 128,
        maximum_point_frequency_pairs: 128,
        maximum_workspace_bytes: 4 << 20,
        maximum_fft_work: 1 << 20,
    };
    assert!(AngularSpectrumExposure::try_new(&beam, &[], exposure_limits, || false).is_err());
    let zero = prepare([0.0; 2], limits()).unwrap();
    let factors = [beam, prepare([0.2, 0.0], limits()).unwrap()];
    assert!(
        AngularSpectrumCoherence::try_new(
            &factors[0],
            &factors,
            AngularSpectrumCoherenceLimits {
                maximum_modes: 2,
                maximum_point_mode_pairs: 128,
                maximum_workspace_bytes: 4 << 20,
                maximum_fft_work: 1 << 20,
            },
            || false,
        )
        .is_err()
    );
    assert!(AngularSpectrumExposure::try_new(&zero, &factors, exposure_limits, || false).is_err());
    assert!(
        PlanarAngularSpectrum::prepare_with_carrier(
            &stack,
            grid(),
            [0.4, 0.0],
            &input,
            limits(),
            || true,
        )
        .is_err()
    );
}

#[test]
fn carrier_aperture_preserves_beat_wavenumbers_below_carrier_roundoff() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let dx = 2.0_f64.powi(60);
    let grid = AngularSpectrumGrid {
        width: 2,
        height: 1,
        spacing_m: [dx, 1.0 / dx],
    };
    // Both physical kx values round to 0.4, but the two original envelope
    // bins remain orthogonal over the full cell. Subtract the shared carrier
    // analytically before preparing any cycle-mean interference kernel.
    assert_eq!(0.4 - core::f64::consts::PI / dx, 0.4);
    let input = [[C64::ZERO, C64::from_re(2.0)], [C64::ZERO; 2]];
    let beam = PlanarAngularSpectrum::prepare_with_carrier(
        &stack,
        grid,
        [0.4, 0.0],
        &input,
        limits(),
        || false,
    )
    .unwrap();
    let detector = AngularSpectrumAperture::try_new([0.0; 2], [2.0 * dx, 1.0 / dx]).unwrap();
    let bounds = AngularSpectrumApertureLimits {
        maximum_mode_pairs: 6,
        maximum_workspace_bytes: 1 << 20,
    };
    let expected = 2.0 * 0.84_f64.sqrt();
    relative(beam.power().incident_w, expected);
    for side in [Side::Incident, Side::Transmitted] {
        let measured = beam
            .measure_aperture(side, 0.0, detector, bounds, || false)
            .unwrap();
        relative(measured.normal_power_w, expected);
        relative(
            measured.normal_power_w,
            beam.field(side, 0.0, || false).unwrap().normal_power_w(),
        );
    }
    let bank = AngularSpectrumCoherence::try_new(
        &beam,
        core::slice::from_ref(&beam),
        AngularSpectrumCoherenceLimits {
            maximum_modes: 2,
            maximum_point_mode_pairs: 128,
            maximum_workspace_bytes: 4 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap();
    relative(
        bank.measure_aperture(Side::Transmitted, 0.0, detector, bounds, || false)
            .unwrap()
            .normal_power_w,
        expected,
    );
}
