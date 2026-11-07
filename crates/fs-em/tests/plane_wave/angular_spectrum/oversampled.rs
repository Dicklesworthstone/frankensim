//! Propagating source support and observation resolution are independent.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumCoherence, AngularSpectrumCoherenceLimits, AngularSpectrumCorrelationLimits,
    AngularSpectrumExposure, AngularSpectrumExposureLimits, AngularSpectrumIlluminationLimits,
    AngularSpectrumIlluminationSide as DualSide, AngularSpectrumLayerLimits,
    AngularSpectrumObservation,
};
use fs_fft::harmonic_window::HarmonicWindow;

fn layer_limits() -> AngularSpectrumLayerLimits {
    AngularSpectrumLayerLimits {
        maximum_layer_mode_pairs: 4096,
        maximum_workspace_bytes: 4 << 20,
        maximum_fft_work: 1 << 20,
    }
}
fn dual_limits() -> AngularSpectrumIlluminationLimits {
    AngularSpectrumIlluminationLimits {
        maximum_layer_mode_pairs: 4096,
        maximum_workspace_bytes: 4 << 20,
        maximum_fft_work: 1 << 20,
    }
}
fn fringe_modes(sign: f64) -> [[C64; 2]; 16] {
    let mut modes = [[C64::ZERO; 2]; 16];
    modes[1][1] = C64::ONE;
    modes[15][1] = C64::from_re(sign);
    modes
}
fn fine_grid() -> AngularSpectrumGrid {
    AngularSpectrumGrid {
        width: 16,
        height: 1,
        spacing_m: [1.0; 2],
    }
}

#[test]
fn g1_empty_modes_skip_unrepresentable_phases_and_still_accept_right_illumination() {
    let medium = material(1.0, 1.0);
    let source =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 2.0, &[])
            .unwrap();
    let empty = [[C64::ZERO; 2]; 16];
    let beam =
        PlanarAngularSpectrum::prepare_modes(&source, fine_grid(), &empty, limits(), || false)
            .unwrap();
    for side in [
        Side::Incident,
        Side::Reflected,
        Side::TotalIncidentSide,
        Side::Transmitted,
    ] {
        // The nonzero wave number times this finite distance overflows, but
        // an exactly empty mode has no propagation phase to evaluate.
        for distance in [0.0, 1.0, f64::MAX] {
            let field = beam.field(side, distance, || false).unwrap();
            assert_eq!(field.normal_power_w(), 0.0);
            assert!(
                field
                    .electric_v_per_m()
                    .iter()
                    .chain(field.magnetic_a_per_m())
                    .all(|&value| value == [C64::ZERO; 3])
            );
        }
        for invalid in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(beam.field(side, invalid, || false).is_err());
        }
        assert!(matches!(
            beam.field(side, 1.0, || true),
            Err(AngularSpectrumError::Cancelled)
        ));
    }

    // Empty in-cone bins retain their physical port geometry: a later source
    // from the right may populate them, independently of the original left beam.
    let mut right = empty;
    right[0][1] = C64::ONE;
    let dual = beam
        .prepare_illumination_modes(&right, dual_limits(), || false)
        .unwrap();
    assert_eq!(dual.power().incoming_w, [0.0, 8.0]);
    assert_eq!(dual.power().outgoing_w, [8.0, 0.0]);
    assert_eq!(dual.power().balance_defect_w(), 0.0);
    for side in [DualSide::TotalLeft, DualSide::TotalRight] {
        let field = dual.field(side, 0.0, || false).unwrap();
        relative(field.normal_power_w(), -8.0);
        for (&electric, &magnetic) in field
            .electric_v_per_m()
            .iter()
            .zip(field.magnetic_a_per_m())
        {
            assert_eq!(electric, [C64::ZERO, C64::ONE, C64::ZERO]);
            assert_eq!(magnetic, [C64::ONE, C64::ZERO, C64::ZERO]);
        }
    }
}

#[test]
fn g1_oversampled_uniform_samples_retain_original_films_and_vector_fields() {
    let medium = material(1.0, 1.0);
    let film = material(4.0, 1.0);
    let source = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&medium),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&film),
            thickness_m: 0.37,
        }],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 8,
        height: 4,
        spacing_m: [0.5, 0.75],
    };
    let amplitude = [C64::new(0.7, 0.2), C64::new(-0.1, 0.9)];
    let samples = [amplitude; 32];
    let beam = PlanarAngularSpectrum::prepare(&source, grid, &samples, limits(), || false).unwrap();
    let mut modes = [[C64::ZERO; 2]; 32];
    modes[0] = amplitude;
    let modal =
        PlanarAngularSpectrum::prepare_modes(&source, grid, &modes, limits(), || false).unwrap();
    let area = 12.0;
    let power = area * 0.5 * (amplitude[0].norm_sq() + amplitude[1].norm_sq());
    relative(beam.power().incident_w, power);
    relative(modal.power().incident_w, power);
    let native = source.response(1.0, Polarization::Te).unwrap();
    relative(
        beam.power().transmitted_w,
        power * native.transmitted_power_fraction,
    );
    relative(
        beam.power().reflected_w,
        power * native.reflected_power_fraction,
    );
    let expected_force = 2.0 * beam.power().reflected_w;
    near(
        beam.radiation_force(|| false).unwrap()[2],
        expected_force,
        3e-13,
    );
    for side in [
        Side::Incident,
        Side::Reflected,
        Side::TotalIncidentSide,
        Side::Transmitted,
    ] {
        let actual = beam.field(side, 0.7, || false).unwrap();
        let authored = modal.field(side, 0.7, || false).unwrap();
        assert_eq!(actual.electric_v_per_m(), authored.electric_v_per_m());
        assert_eq!(actual.magnetic_a_per_m(), authored.magnetic_a_per_m());
    }
    let layer = beam.prepare_layer(0, layer_limits(), || false).unwrap();
    let native_layer = source
        .prepare_layer_field(0, 1.0, Polarization::Te)
        .unwrap();
    let at = native_layer.field(0.2).unwrap();
    let field = layer.field(0.2, || false).unwrap();
    for e in field.field().electric_v_per_m() {
        complex_near(e[0], amplitude[0] * at.electric_v_per_m, 2e-13);
        complex_near(e[1], amplitude[1] * at.electric_v_per_m, 2e-13);
    }
    assert_eq!(field.absorbed_power_density_w_m3(), &[0.0; 32]);
    let dual = beam
        .prepare_illumination(&[[C64::ZERO; 2]; 32], dual_limits(), || false)
        .unwrap();
    near(
        dual.radiation_force(|| false).unwrap()[2],
        expected_force,
        3e-13,
    );
    let transmitted = dual.field(DualSide::TotalRight, 0.7, || false).unwrap();
    let one_sided = beam.field(Side::Transmitted, 0.7, || false).unwrap();
    for (a, b) in transmitted
        .electric_v_per_m()
        .iter()
        .zip(one_sided.electric_v_per_m())
    {
        for c in 0..3 {
            complex_near(a[c], b[c], 3e-13);
        }
    }
}

#[test]
fn g1_oversampled_fringe_modes_preserve_local_aperture_and_ensemble_interference() {
    let medium = material(1.0, 1.0);
    let thickness = 0.8;
    let source = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&medium),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&medium),
            thickness_m: thickness,
        }],
    )
    .unwrap();
    let grid = fine_grid();
    let beams = [1.0, -1.0].map(|sign| {
        PlanarAngularSpectrum::prepare_modes(&source, grid, &fringe_modes(sign), limits(), || false)
            .unwrap()
    });
    let q = core::f64::consts::TAU / 16.0;
    let kz = (1.0 - q * q).sqrt();
    for (beam, sign) in beams.iter().zip([1.0, -1.0]) {
        let field = beam.field(Side::Transmitted, 0.6, || false).unwrap();
        let z_phase = phase(-kz * (thickness + 0.6));
        for x in 0..16 {
            let a = phase(-q * x as f64);
            let b = phase(q * x as f64).scale(sign);
            let ey = (a + b) * z_phase;
            complex_near(field.electric_v_per_m()[x][1], ey, 5e-14);
            complex_near(field.magnetic_a_per_m()[x][0], ey.scale(-kz), 5e-14);
            complex_near(
                field.magnetic_a_per_m()[x][2],
                (a - b).scale(q) * z_phase,
                5e-14,
            );
            near(
                field.normal_flux_w_m2(x).unwrap(),
                kz * (1.0 + sign * (2.0 * q * x as f64).cos()),
                5e-14,
            );
        }
        relative(field.normal_power_w(), 16.0 * kz);
        let aperture = AngularSpectrumAperture::try_new([-0.15, 0.0], [0.3, 1.0]).unwrap();
        let measured = beam
            .measure_aperture(
                Side::Transmitted,
                0.6,
                aperture,
                AngularSpectrumApertureLimits {
                    maximum_mode_pairs: 272,
                    maximum_workspace_bytes: 4 << 20,
                },
                || false,
            )
            .unwrap();
        near(
            measured.mean_poynting_w_m2[2],
            kz * (1.0 + sign * (0.3 * q).sin() / (0.3 * q)),
            5e-14,
        );
    }
    let ensemble = AngularSpectrumCoherence::try_new(
        &beams[0],
        &beams,
        AngularSpectrumCoherenceLimits {
            maximum_modes: 2,
            maximum_point_mode_pairs: 32,
            maximum_workspace_bytes: 4 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap();
    // Direct observations traverse the same dormant-bin-safe Maxwell owner as
    // FFT planes. The two independent fringe factors have Wyy(x,y)=2 cos(q(x-y)).
    let observe = |x| AngularSpectrumObservation {
        side: Side::Transmitted,
        distance_from_face_m: 0.6,
        position_m: [x, 0.0],
    };
    let correlation = ensemble
        .mutual_coherence(
            observe(0.3),
            observe(1.7),
            AngularSpectrumCorrelationLimits {
                maximum_mode_evaluations: 64,
                maximum_workspace_bytes: 4 << 20,
            },
            || false,
        )
        .unwrap();
    complex_near(
        correlation.electric_cross_v2_m2()[1][1],
        C64::from_re(2.0 * (q * (0.3 - 1.7)).cos()),
        5e-14,
    );
    complex_near(
        correlation.electric_at_a_v2_m2()[1][1],
        C64::from_re(2.0),
        5e-14,
    );
    let ensemble_aperture = ensemble
        .measure_aperture(
            Side::Transmitted,
            0.6,
            AngularSpectrumAperture::try_new([-0.15, 0.0], [0.3, 1.0]).unwrap(),
            AngularSpectrumApertureLimits {
                maximum_mode_pairs: 544,
                maximum_workspace_bytes: 4 << 20,
            },
            || false,
        )
        .unwrap();
    near(ensemble_aperture.mean_poynting_w_m2[2], 2.0 * kz, 5e-14);
    near(ensemble_aperture.normal_power_w, 0.6 * kz, 5e-14);
    // Opposite fringe phases add as independent second moments, not fields.
    let ensemble_layer = ensemble.prepare_layer(0, layer_limits(), || false).unwrap();
    for p in ensemble_layer
        .field(0.4, || false)
        .unwrap()
        .field()
        .points()
    {
        near(p.complex_poynting_w_m2[2].re, 2.0 * kz, 5e-14);
    }
    let bank = AngularSpectrumExposure::try_new(
        &beams[0],
        &beams,
        AngularSpectrumExposureLimits {
            maximum_lines: 2,
            maximum_point_line_pairs: 32,
            maximum_point_frequency_pairs: 16,
            maximum_workspace_bytes: 4 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap();
    // Coherently summing those same fields removes the negative spatial mode.
    let exposure = bank
        .field(
            Side::Transmitted,
            0.6,
            HarmonicWindow::new(0.0, core::f64::consts::TAU).unwrap(),
            || false,
        )
        .unwrap();
    for p in exposure.points() {
        near(p.poynting_w_m2[2], 2.0 * kz, 8e-14);
    }
}

#[test]
fn g1_oversampled_two_sided_modes_retain_interior_electric_nodes() {
    let medium = material(1.0, 1.0);
    let source = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&medium),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&medium),
            thickness_m: 0.8,
        }],
    )
    .unwrap();
    let modes = fringe_modes(1.0);
    let beam =
        PlanarAngularSpectrum::prepare_modes(&source, fine_grid(), &modes, limits(), || false)
            .unwrap();
    let right = modes.map(|e| e.map(|z| -z));
    let dual = beam
        .prepare_illumination_modes(&right, dual_limits(), || false)
        .unwrap();
    let mid = dual.layer_field(0, 0.4, || false).unwrap();
    for e in mid.field().electric_v_per_m() {
        for &component in e {
            complex_near(component, C64::ZERO, 1e-13);
        }
    }
    let q = core::f64::consts::TAU / 16.0;
    let kz = (1.0 - q * q).sqrt();
    for (x, h) in mid.field().magnetic_a_per_m().iter().enumerate() {
        complex_near(
            h[0],
            phase(-kz * 0.4).scale(-4.0 * kz * (q * x as f64).cos()),
            1e-13,
        );
    }
    relative(dual.power().incoming_w[0], 16.0 * kz);
    relative(dual.power().incoming_w[1], 16.0 * kz);
    near(dual.power().balance_defect_w(), 0.0, 1e-12);
    for force in dual.radiation_force(|| false).unwrap() {
        near(force, 0.0, 1e-12);
    }
}

#[test]
fn g0_g4_oversampling_never_filters_nonzero_unsupported_sources() {
    let medium = material(1.0, 1.0);
    let source =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let grid = fine_grid();
    let modes = fringe_modes(1.0);
    let beam =
        PlanarAngularSpectrum::prepare_modes(&source, grid, &modes, limits(), || false).unwrap();
    for tiny in [f64::from_bits(1), 1e-300, 1e-20, 1.0] {
        let mut unsupported = modes;
        unsupported[8][0] = C64::from_re(tiny);
        assert!(matches!(
            PlanarAngularSpectrum::prepare_modes(&source, grid, &unsupported, limits(), || false),
            Err(AngularSpectrumError::IncidentLightCone { column: 8, row: 0 })
        ));
        assert!(
            beam.prepare_illumination_modes(&unsupported, dual_limits(), || false)
                .is_err()
        );
    }
    let critical = AngularSpectrumGrid {
        width: 4,
        height: 1,
        spacing_m: [core::f64::consts::PI, 1.0],
    };
    let mut grazing = [[C64::ZERO; 2]; 4];
    grazing[2][1] = C64::ONE;
    assert!(matches!(
        PlanarAngularSpectrum::prepare_modes(&source, critical, &grazing, limits(), || false),
        Err(AngularSpectrumError::IncidentLightCone { column: 2, row: 0 })
    ));
    assert!(
        PlanarAngularSpectrum::prepare_modes(&source, grid, &modes[..15], limits(), || false)
            .is_err()
    );
    for bounded in [
        AngularSpectrumLimits {
            maximum_grid_points: 15,
            ..limits()
        },
        AngularSpectrumLimits {
            maximum_workspace_bytes: beam.workspace_bytes() - 1,
            ..limits()
        },
        AngularSpectrumLimits {
            maximum_fft_work: beam.field_fft_work() - 1,
            ..limits()
        },
    ] {
        assert!(
            PlanarAngularSpectrum::prepare_modes(&source, grid, &modes, bounded, || false).is_err()
        );
    }
    let mut polls = 0;
    assert!(matches!(
        PlanarAngularSpectrum::prepare_modes(&source, grid, &modes, limits(), || {
            polls += 1;
            polls > 4
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert!(matches!(
        beam.prepare_illumination_modes(&modes, dual_limits(), || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    let repeat =
        PlanarAngularSpectrum::prepare_modes(&source, grid, &modes, limits(), || false).unwrap();
    assert_eq!(
        beam.field(Side::Transmitted, 0.4, || false)
            .unwrap()
            .electric_v_per_m(),
        repeat
            .field(Side::Transmitted, 0.4, || false)
            .unwrap()
            .electric_v_per_m(),
    );
}
