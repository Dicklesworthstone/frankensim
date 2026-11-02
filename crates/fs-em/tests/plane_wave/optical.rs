use super::*;
use fs_em::plane_wave::{
    OpticalConductorResponse,
    angular_spectrum::{
        AngularSpectrumExposure, AngularSpectrumExposureLimits, AngularSpectrumGrid,
        AngularSpectrumLayerLimits, AngularSpectrumLimits, AngularSpectrumSide,
        PlanarAngularSpectrum,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_fft::harmonic_window::HarmonicWindow;
use fs_material::state_point::{
    VISIBLE_COMPLEX_IOR_ETA_PROPERTIES, VISIBLE_COMPLEX_IOR_K_PROPERTIES,
    VISIBLE_COMPLEX_IOR_WAVELENGTHS_NM, VisibleConductorStatePoint,
};
use fs_math::det;

const C0: f64 = 299_792_458.0;
const WAVELENGTH_FREQUENCY: f64 = core::f64::consts::TAU * C0 * 1e9;

fn index(wavelength_nm: f64) -> C64 {
    C64::new(
        0.25 + 0.0005 * (wavelength_nm - 380.0),
        -(3.0 + 0.001 * (wavelength_nm - 380.0)),
    )
}

fn metal(mu: Option<f64>) -> VisibleConductorStatePoint {
    let mut properties = Vec::new();
    for i in 0..9 {
        let n = index(VISIBLE_COMPLEX_IOR_WAVELENGTHS_NM[i]);
        properties.push((VISIBLE_COMPLEX_IOR_ETA_PROPERTIES[i], Dims::NONE, n.re));
        properties.push((VISIBLE_COMPLEX_IOR_K_PROPERTIES[i], Dims::NONE, -n.im));
    }
    if let Some(mu) = mu {
        properties.push((MU_NAME, MU, mu));
    }
    VisibleConductorStatePoint::try_from_material(source(
        &properties,
        ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        QueryPoint::new().with("T", 300.0).unwrap(),
    ))
    .unwrap()
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    near(
        actual,
        expected,
        tolerance * actual.abs().max(expected.abs()).max(1e-300),
    );
}

fn host() -> ResolvedMaterialStatePoint {
    // Explicit physical coefficients: this incident medium has n=1. The
    // implementation cannot substitute either coefficient from this fixture.
    material(1.0 / (M * C0 * C0), M)
}

#[test]
fn spectral_metal_reflection_retains_negative_epsilon_phase_and_heat() {
    let incident = host();
    let source = metal(Some(1.4 * M));
    let optical = OpticalConductorResponse::try_new(&source).unwrap();
    let response = DielectricResponse::OpticalConductor(optical);
    let interface = PlanarDielectricInterface::try_new(&incident, response).unwrap();
    assert!(std::ptr::eq(optical.source(), &source));
    assert!(std::ptr::eq(
        interface.transmitted_material(),
        source.resolved()
    ));
    assert_eq!(optical.permeability_h_m(), 1.4 * M);
    // No electric-permittivity or DC-conductivity property is supplied: all
    // electric response comes from the source's spectral complex index.
    assert!(source.resolved().property(EPS_NAME).is_none());
    for wavelength in [380.0, 405.0, 550.0, 725.0, 780.0] {
        let omega = WAVELENGTH_FREQUENCY / wavelength;
        let n = index(wavelength);
        complex_near(optical.complex_index(omega).unwrap(), n, 2e-15);
        let epsilon = (n * n).scale(1.0 / (1.4 * M * C0 * C0));
        assert!(epsilon.re < 0.0 && epsilon.im < 0.0);
        complex_near(
            optical.permittivity_f_m(omega).unwrap(),
            epsilon,
            3e-15 * epsilon.abs(),
        );
        for cosine in [1.0, 0.8, 0.1, 1e-5] {
            let ki = omega / C0;
            let kt = ki * det::sqrt((1.0 - cosine) * (1.0 + cosine));
            // Independent dimensional wave relation and Maxwell admittances.
            let kz = ((n * n).scale(ki * ki) - C64::from_re(kt * kt)).sqrt();
            for pol in [Polarization::Te, Polarization::Tm] {
                let (yi, yt) = match pol {
                    Polarization::Te => {
                        (ki * cosine / (omega * M), kz.scale(1.0 / (omega * 1.4 * M)))
                    }
                    Polarization::Tm => (
                        omega / (M * C0 * C0 * ki * cosine),
                        epsilon.scale(omega) / kz,
                    ),
                };
                let expected = (C64::from_re(yi) - yt) / (C64::from_re(yi) + yt);
                let actual = interface.response(omega, cosine, pol).unwrap();
                complex_near(actual.reflected_tangential_e, expected, 3e-14);
                complex_near(
                    actual.transmitted_normal_wavenumber_rad_m,
                    kz,
                    4e-15 * kz.abs(),
                );
                let et = actual.transmitted_tangential_e;
                let ez = if pol == Polarization::Tm {
                    -et.scale(kt) / kz
                } else {
                    C64::ZERO
                };
                let integrated_heat =
                    -omega * epsilon.im * (et.norm_sq() + ez.norm_sq()) / (-4.0 * kz.im);
                relative(
                    integrated_heat / (0.5 * yi),
                    actual.transmitted_power_fraction,
                    4e-13,
                );
                near(
                    actual.reflected_power_fraction + actual.transmitted_power_fraction,
                    1.0,
                    3e-14,
                );
            }
        }
    }
}

#[test]
fn spectral_metal_films_retain_coherent_transmission_and_opaque_limits() {
    let incident = host();
    let source = metal(Some(M));
    let optical = OpticalConductorResponse::try_new(&source).unwrap();
    let response = DielectricResponse::OpticalConductor(optical);
    let wavelength = 550.0;
    let omega = WAVELENGTH_FREQUENCY / wavelength;
    let n = index(wavelength);
    let epsilon = (n * n).scale(1.0 / (M * C0 * C0));
    for cosine in [1.0, 0.6, 0.02] {
        let ki = omega / C0;
        let kt2 = ki * ki * (1.0 - cosine) * (1.0 + cosine);
        let kz = ((n * n).scale(ki * ki) - C64::from_re(kt2)).sqrt();
        for pol in [Polarization::Te, Polarization::Tm] {
            let (yi, yt) = match pol {
                Polarization::Te => (ki * cosine / (omega * M), kz.scale(1.0 / (omega * M))),
                Polarization::Tm => (
                    omega / (M * C0 * C0 * ki * cosine),
                    epsilon.scale(omega) / kz,
                ),
            };
            let r01 = (C64::from_re(yi) - yt) / (C64::from_re(yi) + yt);
            for thickness in [1e-11, 5e-9, 70e-9, 1000.0 / -kz.im] {
                let stack = PlanarDielectricLayers::try_new(
                    &incident,
                    DielectricResponse::Lossless(&incident),
                    omega,
                    &[DielectricLayer {
                        response,
                        thickness_m: thickness,
                    }],
                )
                .unwrap();
                let actual = stack.response(cosine, pol).unwrap();
                let (s, c) = det::sin_cos(kz.re * thickness);
                let p = C64::new(c, -s).scale(det::exp(kz.im * thickness));
                // Independent two-interface multiple-reflection sum. This
                // moderate-contrast oracle is not the production E/H sweep.
                let denominator = C64::ONE - r01 * r01 * p * p;
                let r = r01 * (C64::ONE - p * p) / denominator;
                let t = (C64::ONE - r01 * r01) * p / denominator;
                complex_near(actual.reflected_tangential_e, r, 3e-12);
                complex_near(actual.transmitted_tangential_e, t, 3e-12);
                near(
                    actual.reflected_power_fraction
                        + actual.transmitted_power_fraction
                        + actual.absorbed_power_fraction,
                    1.0,
                    4e-13,
                );
                assert!(actual.absorbed_power_fraction > 0.0);
                if p == C64::ZERO {
                    complex_near(actual.reflected_tangential_e, r01, 2e-14);
                    assert_eq!(actual.transmitted_power_fraction, 0.0);
                }
                let split = PlanarDielectricLayers::try_new(
                    &incident,
                    DielectricResponse::Lossless(&incident),
                    omega,
                    &[DielectricLayer {
                        response,
                        thickness_m: thickness * 0.5,
                    }; 2],
                )
                .unwrap()
                .response(cosine, pol)
                .unwrap();
                complex_near(
                    split.reflected_tangential_e,
                    actual.reflected_tangential_e,
                    2e-13,
                );
                relative(
                    split.absorbed_power_fraction,
                    actual.absorbed_power_fraction,
                    2e-11,
                );
            }
        }
    }
}

#[test]
fn spectral_metal_beams_and_internal_heat_share_the_same_source_loss() {
    let incident = host();
    let source = metal(Some(M));
    let optical = OpticalConductorResponse::try_new(&source).unwrap();
    let response = DielectricResponse::OpticalConductor(optical);
    let omega = WAVELENGTH_FREQUENCY / 550.0;
    let thickness = 35e-9;
    let stack = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&incident),
        omega,
        &[DielectricLayer {
            response,
            thickness_m: thickness,
        }],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 4,
        height: 4,
        spacing_m: [700e-9, 800e-9],
    };
    let mut input = vec![[C64::ZERO; 2]; 16];
    input[2] = [C64::ONE, C64::new(0.2, -0.3)];
    input[11] = [C64::new(0.1, 0.8), C64::new(-0.4, 0.2)];
    let limits = AngularSpectrumLimits {
        maximum_grid_points: 16,
        maximum_workspace_bytes: 1 << 20,
        maximum_fft_work: 1 << 16,
    };
    let beam = PlanarAngularSpectrum::prepare(&stack, grid, &input, limits, || false).unwrap();
    let layer = beam
        .prepare_layer(
            0,
            AngularSpectrumLayerLimits {
                maximum_layer_mode_pairs: 32,
                maximum_workspace_bytes: 1 << 20,
                maximum_fft_work: 1 << 16,
            },
            || false,
        )
        .unwrap();
    relative(
        layer.absorbed_power_w(),
        beam.power().film_absorbed_w,
        3e-12,
    );
    let epsilon = optical.permittivity_f_m(omega).unwrap();
    let intervals = 128;
    let dz = thickness / intervals as f64;
    let mut heat = 0.0;
    for i in 0..=intervals {
        let plane = layer.field(i as f64 * dz, || false).unwrap();
        for (e, &density) in plane
            .field()
            .electric_v_per_m()
            .iter()
            .zip(plane.absorbed_power_density_w_m3())
        {
            relative(
                density,
                -0.5 * omega * epsilon.im * e.iter().map(|v| v.norm_sq()).sum::<f64>(),
                3e-15,
            );
        }
        let weight = if i == 0 || i == intervals {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        heat += weight * plane.absorbed_power_per_depth_w_m();
    }
    relative(heat * dz / 3.0, layer.absorbed_power_w(), 2e-8);
    assert!(beam.power().balance_defect_w().abs() < 3e-12 * beam.power().incident_w);

    // A metal EXIT also needs negative-real-epsilon admission in the original
    // vector-beam owner, not just a film response that happens to accept it.
    let half = PlanarDielectricLayers::try_new(&incident, response, omega, &[]).unwrap();
    let transmitted =
        PlanarAngularSpectrum::prepare(&half, grid, &input, limits, || false).unwrap();
    assert!(transmitted.power().transmitted_w > 0.0);
    assert_eq!(transmitted.power().film_absorbed_w, 0.0);
    assert!(transmitted.power().balance_defect_w().abs() < 3e-12 * transmitted.power().incident_w);
}

#[test]
fn spectral_wave_admission_preserves_source_band_and_required_magnetic_data() {
    let incomplete = metal(None);
    assert!(OpticalConductorResponse::try_new(&incomplete).is_err());
    let source = metal(Some(M));
    let optical = OpticalConductorResponse::try_new(&source).unwrap();
    let minimum = WAVELENGTH_FREQUENCY / 780.0;
    let maximum = WAVELENGTH_FREQUENCY / 380.0;
    for omega in [minimum, maximum] {
        assert!(optical.permittivity_f_m(omega).is_ok());
    }
    let incident = host();
    for omega in [
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        f64::from_bits(minimum.to_bits() - 1),
        f64::from_bits(maximum.to_bits() + 1),
    ] {
        assert!(optical.permittivity_f_m(omega).is_err());
        assert!(
            PlanarDielectricLayers::try_new(
                &incident,
                DielectricResponse::Lossless(&incident),
                omega,
                &[DielectricLayer {
                    response: DielectricResponse::OpticalConductor(optical),
                    thickness_m: 0.0
                }],
            )
            .is_err()
        );
    }
}

#[test]
fn spectral_metal_exposure_retains_one_original_source_across_frequencies() {
    let incident = host();
    let source = metal(Some(M));
    let grid = AngularSpectrumGrid {
        width: 1,
        height: 1,
        spacing_m: [1e-6, 2e-6],
    };
    let limits = AngularSpectrumLimits {
        maximum_grid_points: 1,
        maximum_workspace_bytes: 1 << 20,
        maximum_fft_work: 1 << 16,
    };
    let stacks = [550.0, 650.0].map(|wavelength| {
        // Independently prepared adapters keep the identical original source.
        PlanarDielectricLayers::try_new(
            &incident,
            DielectricResponse::OpticalConductor(
                OpticalConductorResponse::try_new(&source).unwrap(),
            ),
            WAVELENGTH_FREQUENCY / wavelength,
            &[],
        )
        .unwrap()
    });
    let lines = stacks
        .iter()
        .map(|stack| {
            PlanarAngularSpectrum::prepare(
                stack,
                grid,
                &[[C64::new(1.0, 0.2), C64::ZERO]],
                limits,
                || false,
            )
            .unwrap()
        })
        .collect::<Vec<_>>();
    let exposure_limits = AngularSpectrumExposureLimits {
        maximum_lines: 2,
        maximum_point_line_pairs: 2,
        maximum_point_frequency_pairs: 4,
        maximum_workspace_bytes: 2 << 20,
        maximum_fft_work: 1 << 18,
    };
    let exposure =
        AngularSpectrumExposure::try_new(&lines[0], &lines, exposure_limits, || false).unwrap();
    // These two optical frequencies are11 and13 times a common fundamental.
    // One full fundamental period averages every cross/carrier term to zero.
    let window = HarmonicWindow::new(0.0, 7150.0e-9 / C0).unwrap();
    let reflected = exposure
        .field(AngularSpectrumSide::Reflected, 0.0, window, || false)
        .unwrap();
    relative(
        reflected.normal_power_w(),
        -lines
            .iter()
            .map(|line| line.power().reflected_w)
            .sum::<f64>(),
        2e-12,
    );
    // Numerically identical constants from another owner cannot silently
    // replace the original spectral line bank's physical source identity.
    let foreign_source = source.clone();
    let foreign = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::OpticalConductor(
            OpticalConductorResponse::try_new(&foreign_source).unwrap(),
        ),
        WAVELENGTH_FREQUENCY / 550.0,
        &[],
    )
    .unwrap();
    let foreign_beam =
        PlanarAngularSpectrum::prepare(&foreign, grid, &[[C64::ONE, C64::ZERO]], limits, || false)
            .unwrap();
    assert!(
        AngularSpectrumExposure::try_new(&foreign_beam, &lines, exposure_limits, || false).is_err()
    );
}
