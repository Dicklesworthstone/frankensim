use super::*;
use crate::{
    conductor::{
        CoatingIorSample, CoatingLayer, CoherentAbsorbingCoating, ConductorDataStatus,
        ConductorSource,
    },
    dielectric::{DielectricIorSample, TabulatedIor},
    thermal_surface::{
        CoatingThermalEmission,
        crystal_emission::{
            UniaxialSlabThermalEmission,
            tests::{Fixture, close, planck_si},
        },
    },
};
use fs_blake3::ContentHash;
use fs_em::plane_wave::{
    DielectricResponse, UniaxialOpticalResponse,
    uniaxial::oriented::{
        slab::OrientedUniaxialSlab,
        stack::{UniaxialStackLayer, UniaxialStackLimits},
    },
};

const MU0: f64 = 4.0e-7 * core::f64::consts::PI;
const LAMBDA: f64 = 10_000.0;

fn layer(fixture: &Fixture, axis: [f64; 3], depth: f64) -> UniaxialStackLayer<'_> {
    UniaxialStackLayer::try_new(
        UniaxialOpticalResponse::try_new(&fixture.ordinary, &fixture.extraordinary, 6.1e-6)
            .unwrap(),
        axis,
        depth,
    )
    .unwrap()
}

fn stack<'s, 'm>(
    exterior: &'m Fixture,
    layers: &'s [UniaxialStackLayer<'m>],
) -> OrientedUniaxialStack<'s, 'm> {
    OrientedUniaxialStack::try_new(
        DielectricResponse::Lossless(&exterior.exterior),
        DielectricResponse::Lossless(&exterior.exterior),
        layers,
        UniaxialStackLimits {
            maximum_layers: 8,
            maximum_workspace_bytes: 1024 * 1024,
        },
    )
    .unwrap()
}

#[test]
fn one_original_layer_matches_slab_for_complex_analyzers_and_stokes() {
    let fixture = Fixture::new([[1.55, 0.08], [2.1, 0.25]], 800.0, 1.0, MU0);
    let entrance = fixture.interface([1.0, 2.0, 3.0]);
    let slab = OrientedUniaxialSlab::try_new(
        &entrance,
        DielectricResponse::Lossless(&fixture.exterior),
        4.0e-6,
    )
    .unwrap();
    let layers = [layer(&fixture, [1.0, 2.0, 3.0], 4.0e-6)];
    let source = stack(&fixture, &layers);
    let emitter =
        UniaxialStackThermalEmission::try_from_material_temperatures(&source, &["T"], || false)
            .unwrap();
    assert!(core::ptr::eq(emitter.source(), &source));
    assert_eq!(emitter.layer_temperatures_k(), &[800.0]);
    for cosine in [1.0, 0.63] {
        let response = emitter
            .prepare(LAMBDA, cosine, [3.0, 4.0], || false)
            .unwrap();
        let single = UniaxialSlabThermalEmission::try_from_material_temperature(&slab, "T")
            .unwrap()
            .prepare(LAMBDA, cosine, [3.0, 4.0])
            .unwrap();
        assert!(core::ptr::eq(response.source(), &source));
        assert!(core::ptr::eq(
            response.source().layers()[0]
                .crystal_response()
                .ordinary_source(),
            &fixture.ordinary
        ));
        assert_eq!(response.layer_count(), 1);
        assert_eq!(response.layer_temperature_k(0), Some(800.0));
        assert_eq!(response.layer_temperature_k(1), None);
        assert_eq!(response.wavelength_nm(), LAMBDA);
        assert_eq!(response.outgoing_direction(), single.outgoing_direction());
        assert_eq!(response.polarization_axes(), single.polarization_axes());
        for analyzer in [
            [C64::ONE, C64::ZERO],
            [C64::ZERO, C64::ONE],
            [C64::ONE, C64::new(0.0, 1.0)],
            [C64::new(0.3, -0.2), C64::new(-0.1, 0.7)],
        ] {
            let expected = single.analyzer_radiance_w_m2_sr_nm(analyzer).unwrap();
            close(
                response
                    .analyzer_radiance_w_m2_sr_nm(analyzer, || false)
                    .unwrap(),
                expected,
            );
            close(
                response
                    .layer_analyzer_radiance_w_m2_sr_nm(0, analyzer, || false)
                    .unwrap(),
                expected,
            );
        }
        let actual = response.stokes(|| false).unwrap();
        let expected = single.stokes().unwrap();
        for (a, b) in [actual.i, actual.q, actual.u, actual.v]
            .into_iter()
            .zip([expected.i, expected.q, expected.u, expected.v])
        {
            assert!((a - b).abs() < 3.0e-10 * expected.i);
        }
        close(response.radiance_w_m2_sr_nm(|| false).unwrap(), expected.i);
    }
}

#[test]
fn distinct_layer_temperatures_weight_original_reciprocal_work_without_averaging() {
    let cold = Fixture::new([[1.45, 0.05], [2.0, 0.12]], 350.0, 1.0, MU0);
    let hot = Fixture::new([[1.7, 0.18], [2.3, 0.07]], 1100.0, 1.0, MU0);
    let layers = [
        layer(&cold, [1.0, 2.0, 1.0], 1.2e-6),
        layer(&hot, [-2.0, 1.0, 3.0], 3.1e-6),
    ];
    let source = stack(&cold, &layers);
    let emitter =
        UniaxialStackThermalEmission::try_from_material_temperatures(&source, &["T", "T"], || {
            false
        })
        .unwrap();
    assert_eq!(emitter.layer_temperatures_k(), &[350.0, 1100.0]);
    let c = 0.65_f64;
    let u = [0.8, 0.6];
    let mut differs_from_average = false;
    for wavelength in [8500.0, 12_000.0] {
        let response = emitter.prepare(wavelength, c, u, || false).unwrap();
        let omega = core::f64::consts::TAU * LIGHT_M_S * 1.0e9 / wavelength;
        for analyzer in [
            [C64::ONE, C64::ZERO],
            [C64::ZERO, C64::ONE],
            [C64::ONE, C64::new(0.2, 0.7)],
        ] {
            let norm = (analyzer[0].norm_sq() + analyzer[1].norm_sq()).sqrt();
            let a = analyzer.map(|v| v.scale(1.0 / norm));
            // Independent whole-stack illumination, not the polarization
            // cache or isolated-film absorption. Keep the reversed direction
            // and conjugated full global transverse analyzer explicitly.
            let electric = [
                (a[0].scale(-u[1]) + a[1].scale(c * u[0])).conj(),
                (a[0].scale(u[0]) + a[1].scale(c * u[1])).conj(),
            ];
            let field = source
                .prepare_absorbing_field_in_plane(omega, c, [-u[0], -u[1]], electric, || false)
                .unwrap();
            let incident = field.response().incident_normal_power_w_m2;
            let mut expected = 0.0;
            let mut absorbed = 0.0;
            let mut contributions = [0.0; 2];
            for (index, temperature) in [350.0, 1100.0].into_iter().enumerate() {
                assert_eq!(response.layer_temperature_k(index), Some(temperature));
                let work = field
                    .layer_absorption(index, 1.0, 0.0, layers[index].thickness_m())
                    .unwrap();
                let fraction = work.absorbed_power_w() / incident;
                let radiance = 0.5 * planck_si(temperature, wavelength) * fraction;
                close(
                    response
                        .layer_analyzer_radiance_w_m2_sr_nm(index, analyzer, || false)
                        .unwrap(),
                    radiance,
                );
                contributions[index] = radiance;
                expected += radiance;
                absorbed += fraction;
            }
            let actual = response
                .analyzer_radiance_w_m2_sr_nm(analyzer, || false)
                .unwrap();
            close(actual, expected);
            let average = (350.0 * 1.2 + 1100.0 * 3.1) / 4.3;
            let averaged_emission = 0.5 * planck_si(average, wavelength) * absorbed;
            differs_from_average |= (actual - averaged_emission).abs() > 0.01 * actual;
            // The positive sources are independent even though their
            // reciprocal optical fields retain all coherent reflections.
            assert!(actual < (contributions[0].sqrt() + contributions[1].sqrt()).powi(2));
        }
    }
    assert!(differs_from_average);
}

#[test]
fn isotropic_multilayer_emission_matches_existing_per_film_thermal_owner() {
    let first = Fixture::new([[1.55, 0.03]; 2], 400.0, 1.0, MU0);
    let second = Fixture::new([[2.05, 0.17]; 2], 900.0, 1.0, MU0);
    let layers = [
        layer(&first, [1.0, 2.0, 3.0], 1.3e-6),
        layer(&second, [3.0, 1.0, -2.0], 2.4e-6),
    ];
    let source = stack(&first, &layers);
    let emitter =
        UniaxialStackThermalEmission::try_from_material_temperatures(&source, &["T", "T"], || {
            false
        })
        .unwrap();
    let film = |eta, extinction, depth, tag| {
        CoherentAbsorbingCoating::try_new_in_band(
            core::array::from_fn(|i| {
                CoatingIorSample::try_new_in_band(8000.0 + 750.0 * i as f64, eta, extinction)
                    .unwrap()
            }),
            depth,
            750.0,
            ConductorSource::try_new(ContentHash([tag; 32]), ConductorDataStatus::Representative)
                .unwrap(),
        )
        .unwrap()
    };
    let films = [film(1.55, 0.03, 1.3e-6, 93), film(2.05, 0.17, 2.4e-6, 94)];
    let vacuum = TabulatedIor::try_new_in_band(
        core::array::from_fn(|i| {
            DielectricIorSample::try_new_in_band(8000.0 + 750.0 * i as f64, 1.0).unwrap()
        }),
        750.0,
    )
    .unwrap()
    .into();
    let scalar = CoatingThermalEmission::try_transmitting(
        &[
            CoatingLayer::Absorbing(&films[0]),
            CoatingLayer::Absorbing(&films[1]),
        ],
        &[400.0, 900.0],
        vacuum,
        false,
        ContentHash([95; 32]),
    )
    .unwrap();
    for wavelength in [9000.0, 12_000.0] {
        for cosine in [1.0, 0.67] {
            let response = emitter
                .prepare(wavelength, cosine, [3.0, 4.0], || false)
                .unwrap();
            let expected = scalar.evaluate(wavelength, cosine).unwrap();
            for (pol, analyzer) in [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]]
                .into_iter()
                .enumerate()
            {
                close(
                    response
                        .analyzer_radiance_w_m2_sr_nm(analyzer, || false)
                        .unwrap(),
                    expected.polarized_radiance_w_m2_sr_nm()[pol],
                );
                for layer in 0..2 {
                    close(
                        response
                            .layer_analyzer_radiance_w_m2_sr_nm(layer, analyzer, || false)
                            .unwrap(),
                        expected.film_radiance_w_m2_sr_nm()[layer][pol],
                    );
                }
            }
        }
    }
}

#[test]
fn hot_lossless_layer_is_dark_and_weak_absorbing_polarization_remains_positive() {
    let transparent = Fixture::new([[1.4, 0.0], [1.9, 0.0]], 1500.0, 1.0, MU0);
    let weak = Fixture::new([[1.45, 0.12], [2.1, 1.0e-24]], 800.0, 1.0, MU0);
    let layers = [
        layer(&transparent, [1.0, 0.0, 0.0], 2.0e-6),
        layer(&weak, [1.0, 0.0, 0.0], 4.0e-6),
    ];
    let source = stack(&transparent, &layers);
    let response =
        UniaxialStackThermalEmission::try_from_material_temperatures(&source, &["T", "T"], || {
            false
        })
        .unwrap()
        .prepare(LAMBDA, 1.0, [1.0, 0.0], || false)
        .unwrap();
    for analyzer in [
        [C64::ONE, C64::ZERO],
        [C64::ZERO, C64::ONE],
        [C64::ONE, C64::new(0.0, 1.0)],
    ] {
        assert_eq!(
            response
                .layer_analyzer_radiance_w_m2_sr_nm(0, analyzer, || false)
                .unwrap(),
            0.0
        );
    }
    let strong = response
        .analyzer_radiance_w_m2_sr_nm([C64::ONE, C64::ZERO], || false)
        .unwrap();
    let faint = response
        .analyzer_radiance_w_m2_sr_nm([C64::ZERO, C64::ONE], || false)
        .unwrap();
    assert!(faint > 0.0 && faint < 1.0e-20 * strong);
    assert_eq!(
        faint,
        response
            .layer_analyzer_radiance_w_m2_sr_nm(1, [C64::ZERO, C64::ONE], || false)
            .unwrap()
    );
    let stokes = response.stokes(|| false).unwrap();
    assert_eq!(stokes.i - stokes.q, 0.0);
}

#[test]
fn original_layer_temperatures_counts_and_cancellation_are_enforced() {
    let fixture = Fixture::new([[1.55, 0.08], [2.1, 0.25]], 800.0, 1.0, MU0);
    let layers = [
        layer(&fixture, [1.0, 2.0, 3.0], 1.2e-6),
        layer(&fixture, [3.0, 2.0, 1.0], 3.1e-6),
    ];
    let source = stack(&fixture, &layers);
    assert!(matches!(
        UniaxialStackThermalEmission::try_from_material_temperatures(&source, &["T"], || false),
        Err(CrystalStackEmissionError::LayerTemperatures {
            expected: 2,
            observed: 1
        })
    ));
    assert!(matches!(
        UniaxialStackThermalEmission::try_from_material_temperatures(
            &source,
            &["T", "missing"],
            || false
        ),
        Err(CrystalStackEmissionError::LayerSource {
            layer: 1,
            source: CrystalEmissionError::TemperatureAxis
        })
    ));
    assert!(matches!(
        UniaxialStackThermalEmission::try_layer_temperatures(&source, &[800.0, f64::NAN], || false),
        Err(CrystalStackEmissionError::LayerSource { layer: 1, .. })
    ));
    assert!(matches!(
        UniaxialStackThermalEmission::try_layer_temperatures(&source, &[800.0, 900.0], || true),
        Err(CrystalStackEmissionError::Cancelled)
    ));
    let emitter =
        UniaxialStackThermalEmission::try_layer_temperatures(&source, &[800.0, 900.0], || false)
            .unwrap();
    let mut polls = 0;
    assert!(matches!(
        emitter.prepare(LAMBDA, 0.7, [1.0, 0.0], || {
            polls += 1;
            polls >= 4
        }),
        Err(CrystalStackEmissionError::Cancelled)
    ));
    let response = emitter.prepare(LAMBDA, 0.7, [1.0, 0.0], || false).unwrap();
    assert_eq!(response.layer_temperature_k(1), Some(900.0));
    let analyzer = [C64::ONE, C64::new(0.2, 0.7)];
    let expected = response
        .analyzer_radiance_w_m2_sr_nm(analyzer, || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        response.analyzer_radiance_w_m2_sr_nm(analyzer, || {
            polls += 1;
            polls == 3
        }),
        Err(CrystalStackEmissionError::Cancelled)
    ));
    assert_eq!(
        response
            .analyzer_radiance_w_m2_sr_nm(analyzer, || false)
            .unwrap(),
        expected
    );
    assert!(matches!(
        response.layer_analyzer_radiance_w_m2_sr_nm(2, analyzer, || false),
        Err(CrystalStackEmissionError::LayerIndex { index: 2, count: 2 })
    ));
    assert!(matches!(
        response.stokes(|| true),
        Err(CrystalStackEmissionError::Cancelled)
    ));
    assert!(emitter.prepare(7000.0, 0.7, [1.0, 0.0], || false).is_err());
}
