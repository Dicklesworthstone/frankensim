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
fn g1_authored_oblique_modes_resolve_fine_grid_fields_power_and_detector_fringes() {
    let air = material(1.0, 1.0);
    let glass = material(2.25, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&air, DielectricResponse::Lossless(&glass), 1.0, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 16,
        height: 1,
        spacing_m: [1.0; 2],
    };
    let carrier = [0.23, 0.0];
    let amplitudes = [C64::new(0.7, 0.2), C64::new(-0.4, 0.6)];
    let mut modes = [[C64::ZERO; 2]; 16];
    modes[1][1] = amplitudes[0];
    modes[15][1] = amplitudes[1];
    // Most bins of this fine grid lie outside the incident light cone. The
    // authored support is exact: analyzing synthesized samples is unnecessary.
    let beam = PlanarAngularSpectrum::prepare_modes_with_carrier(
        &stack,
        grid,
        carrier,
        &modes,
        limits(),
        || false,
    )
    .unwrap();
    assert_eq!(beam.carrier_transverse_wavenumber_rad_m(), carrier);

    // Independent TE Fresnel amplitudes and Maxwell reconstruction at two
    // distinct physical angles; these references use no native film response.
    let transverse = [1.0, -1.0].map(|bin| carrier[0] + core::f64::consts::TAU * bin / 16.0);
    let incident_q = transverse.map(|kx| (1.0 - kx * kx).sqrt());
    let exit_q = transverse.map(|kx| (2.25 - kx * kx).sqrt());
    let reflection: [f64; 2] =
        std::array::from_fn(|i| (incident_q[i] - exit_q[i]) / (incident_q[i] + exit_q[i]));
    let transmission: [f64; 2] =
        std::array::from_fn(|i| 2.0 * incident_q[i] / (incident_q[i] + exit_q[i]));
    let incident_w: f64 = (0..2)
        .map(|i| 8.0 * incident_q[i] * amplitudes[i].norm_sq())
        .sum();
    let reflected_w: f64 = (0..2)
        .map(|i| 8.0 * incident_q[i] * (amplitudes[i].scale(reflection[i])).norm_sq())
        .sum();
    let transmitted_w: f64 = (0..2)
        .map(|i| 8.0 * exit_q[i] * (amplitudes[i].scale(transmission[i])).norm_sq())
        .sum();
    relative(beam.power().incident_w, incident_w);
    relative(beam.power().reflected_w, reflected_w);
    relative(beam.power().transmitted_w, transmitted_w);
    near(beam.power().balance_defect_w(), 0.0, 5e-14);
    for distance in [0.0, 0.73, 7.5] {
        let field = beam.field(Side::Transmitted, distance, || false).unwrap();
        for x in 0..16 {
            let mut electric = [C64::ZERO; 3];
            let mut magnetic = [C64::ZERO; 3];
            for i in 0..2 {
                let ey = amplitudes[i].scale(transmission[i])
                    * phase(-transverse[i] * x as f64 - exit_q[i] * distance);
                electric[1] = electric[1] + ey;
                magnetic[0] = magnetic[0] - ey.scale(exit_q[i]);
                magnetic[2] = magnetic[2] + ey.scale(transverse[i]);
            }
            for c in 0..3 {
                complex_near(field.electric_v_per_m()[x][c], electric[c], 2e-13);
                complex_near(field.magnetic_a_per_m()[x][c], magnetic[c], 2e-13);
            }
        }
        relative(field.normal_power_w(), transmitted_w);
    }

    let detector = AngularSpectrumAperture::try_new([2.3, 0.2], [0.83, 0.7]).unwrap();
    let distance = 0.73;
    let mut expected_flux = 0.0;
    for i in 0..2 {
        for j in 0..2 {
            let delta = transverse[i] - transverse[j];
            let argument = 0.5 * delta * detector.extent_m()[0];
            let sinc = if argument == 0.0 {
                1.0
            } else {
                argument.sin() / argument
            };
            let center = detector.origin_m()[0] + 0.5 * detector.extent_m()[0];
            let cross = amplitudes[i]
                * amplitudes[j].conj()
                * phase(-delta * center - (exit_q[i] - exit_q[j]) * distance);
            expected_flux += 0.5 * exit_q[j] * transmission[i] * transmission[j] * sinc * cross.re;
        }
    }
    let reading = beam
        .measure_aperture(
            Side::Transmitted,
            distance,
            detector,
            AngularSpectrumApertureLimits {
                maximum_mode_pairs: 512,
                maximum_workspace_bytes: 4 << 20,
            },
            || false,
        )
        .unwrap();
    relative(reading.normal_power_w, expected_flux * detector.area_m2());
}

#[test]
fn g0_g4_modal_carriers_preserve_exact_support_zero_carrier_and_cancellation() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 16,
        height: 1,
        spacing_m: [1.0; 2],
    };
    let mut modes = [[C64::ZERO; 2]; 16];
    modes[1] = [C64::new(0.2, -0.1), C64::new(0.7, 0.2)];
    modes[15] = [C64::new(-0.3, 0.4), C64::new(-0.4, 0.6)];
    let original =
        PlanarAngularSpectrum::prepare_modes(&stack, grid, &modes, limits(), || false).unwrap();
    let zero = PlanarAngularSpectrum::prepare_modes_with_carrier(
        &stack,
        grid,
        [0.0; 2],
        &modes,
        limits(),
        || false,
    )
    .unwrap();
    assert_eq!(original.power(), zero.power());
    assert_eq!(original.workspace_bytes(), zero.workspace_bytes());
    for side in [
        Side::Incident,
        Side::Reflected,
        Side::TotalIncidentSide,
        Side::Transmitted,
    ] {
        let a = original.field(side, 0.4, || false).unwrap();
        let b = zero.field(side, 0.4, || false).unwrap();
        assert_eq!(a.normal_power_w().to_bits(), b.normal_power_w().to_bits());
        for (left, right) in a
            .electric_v_per_m()
            .iter()
            .chain(a.magnetic_a_per_m())
            .zip(b.electric_v_per_m().iter().chain(b.magnetic_a_per_m()))
        {
            assert_eq!(
                left.map(|z| [z.re.to_bits(), z.im.to_bits()]),
                right.map(|z| [z.re.to_bits(), z.im.to_bits()])
            );
        }
    }
    for tiny in [f64::from_bits(1), 1e-20, 1.0] {
        let mut unsupported = modes;
        unsupported[8][0] = C64::from_re(tiny);
        assert!(matches!(
            PlanarAngularSpectrum::prepare_modes_with_carrier(
                &stack,
                grid,
                [0.23, 0.0],
                &unsupported,
                limits(),
                || false,
            ),
            Err(AngularSpectrumError::IncidentLightCone { column: 8, row: 0 })
        ));
    }
    assert!(matches!(
        PlanarAngularSpectrum::prepare_modes_with_carrier(
            &stack,
            grid,
            [0.23, 0.0],
            &modes,
            limits(),
            || true,
        ),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(
        PlanarAngularSpectrum::prepare_modes_with_carrier(
            &stack,
            grid,
            [f64::NAN, 0.0],
            &modes,
            limits(),
            || false,
        )
        .is_err()
    );
}

#[test]
fn adjacent_light_cone_carriers_keep_their_axial_wavelength_and_poynting_flux() {
    let medium = material(1.0, 1.0);
    // G1/G3: independent high-precision axial cosines for the exact authored floats.
    for (x, y, cosine) in [
        (0.8, 0.6_f64.next_down(), 9.424_321_830_774_483e-9),
        (0.6, 0.8_f64.next_down(), 1.154_238_982_858_484e-8),
    ] {
        for exponent in [-400, 0, 400] {
            let scale = 2.0_f64.powi(exponent);
            let stack = PlanarDielectricLayers::try_new(
                &medium,
                DielectricResponse::Lossless(&medium),
                scale,
                &[],
            )
            .unwrap();
            let grid = AngularSpectrumGrid {
                width: 1,
                height: 1,
                spacing_m: [scale.recip(); 2],
            };
            let input = [[C64::ONE, C64::ZERO]];
            let beam = PlanarAngularSpectrum::prepare_with_carrier(
                &stack,
                grid,
                [x * scale, y * scale],
                &input,
                limits(),
                || false,
            )
            .unwrap();
            // Maxwell: Ez=-kx/kz and Hy=(kz+kx²/kz)/(omega*mu).
            let expected_ez = -x / cosine;
            let expected_power =
                0.5 * grid.spacing_m[0] * grid.spacing_m[1] * (cosine + x * x / cosine);
            let check_power = |actual: f64| {
                assert!((actual / expected_power - 1.0).abs() < 2e-14);
            };
            check_power(beam.power().incident_w);
            check_power(beam.power().transmitted_w);
            assert_eq!(beam.power().reflected_w, 0.0);
            let wavelength_z_m = core::f64::consts::TAU / (scale * cosine);
            for (fraction, gain) in [(0.5, -1.0), (1.0, 1.0)] {
                let field = beam
                    .field(Side::Transmitted, fraction * wavelength_z_m, || false)
                    .unwrap();
                for c in 0..2 {
                    complex_near(
                        field.electric_v_per_m()[0][c],
                        input[0][c].scale(gain),
                        2e-14,
                    );
                }
                complex_near(
                    field.electric_v_per_m()[0][2],
                    C64::from_re(gain * expected_ez),
                    2e-14 * expected_ez.abs(),
                );
                check_power(field.normal_power_w());
            }
            // One adjacent authored direction lies outside the physical cone.
            assert!(matches!(
                PlanarAngularSpectrum::prepare_with_carrier(
                    &stack,
                    grid,
                    [0.8 * scale, 0.6 * scale],
                    &input,
                    limits(),
                    || false,
                ),
                Err(AngularSpectrumError::IncidentLightCone { column: 0, row: 0 })
            ));
        }
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
    // Only the constant envelope's DC bin is occupied; shifted dark bins
    // outside the cone do not forbid this physical near-grazing plane wave.
    assert!(prepare([0.99, 0.0], limits()).is_ok());
    assert!(matches!(
        prepare([1.01, 0.0], limits()),
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
    assert!(AngularSpectrumExposure::try_new(&beam, &[], exposure_limits, || false).is_ok());
    let matched = [prepare([0.4, 0.0], limits()).unwrap()];
    assert!(
        AngularSpectrumExposure::try_new(
            &matched[0],
            &matched,
            AngularSpectrumExposureLimits {
                maximum_point_frequency_pairs: 1024,
                ..exposure_limits
            },
            || false,
        )
        .is_ok()
    );
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
