//! Original finite-band source owners through Maxwell interfaces, films and beams.
use super::*;
use fs_em::plane_wave::{
    SpectralOpticalResponse, SpectralOpticalSource,
    angular_spectrum::{
        AngularSpectrumExposure, AngularSpectrumExposureLayerLimits, AngularSpectrumExposureLimits,
        AngularSpectrumGrid, AngularSpectrumLimits, AngularSpectrumSide, PlanarAngularSpectrum,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_material::state_point::{
    SPECTRAL_CONDUCTOR_ETA_PROPERTIES, SPECTRAL_CONDUCTOR_K_PROPERTIES,
    SPECTRAL_CONDUCTOR_K_QUANTITY, SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
    SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY, SPECTRAL_DIELECTRIC_IOR_PROPERTIES,
    SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES, SpectralComplexIndexKnot,
    SpectralComplexIndexStatePoint, SpectralConductorStatePoint,
    SpectralTabulatedDielectricStatePoint,
};
use fs_qty::{
    QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const C0: f64 = 299_792_458.0;
const WL_FREQUENCY_M: f64 = core::f64::consts::TAU * C0;
const WL_FREQUENCY_NM: f64 = WL_FREQUENCY_M * 1e9;

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature(), t)
        .unwrap()
}
fn claim(key: PropertyKey, cold: f64, hot: f64) -> PropertyClaim {
    let dims = key.dims();
    PropertyClaim {
        key,
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: fs_qty::Temperature::DIMS,
            knots: vec![(300.0, cold), (600.0, hot)],
            dims,
        },
        validity: ValidityDomain::unconstrained().with_quantity("T", temperature(), 300.0, 600.0),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "authored finite-band Maxwell fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn card(claims: ClaimSet) -> MaterialCard {
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored infrared optical law".into(),
            phase: "solid".into(),
            process: "source-selected wavelength grid".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn resolve(card: &MaterialCard) -> ResolvedMaterialStatePoint {
    let requirements: Vec<_> = card
        .claims()
        .claims_ordered()
        .map(|(_, claim)| {
            ScalarPropertyRequirement::try_with_quantity(
                claim.key.name(),
                claim.key.quantity(),
                ScalarAdmissibility::NonNegative,
            )
            .unwrap()
        })
        .collect();
    resolve_material_state_point(
        card,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn generic_source(lossy: bool, with_mu: bool) -> (MaterialCard, SpectralComplexIndexStatePoint) {
    let mut claims = ClaimSet::new();
    let mut bindings = Vec::new();
    for (i, wavelength_nm) in [8000.0, 10_000.0, 14_000.0].into_iter().enumerate() {
        let n_key = PropertyKey::new(format!("original_n_{i}"), Dims::NONE);
        let k_key =
            PropertyKey::with_quantity(format!("original_k_{i}"), SPECTRAL_CONDUCTOR_K_QUANTITY);
        let n = 1.4 + 1e-5 * (wavelength_nm - 8000.0);
        let k = if lossy {
            0.05 + 1e-5 * (wavelength_nm - 8000.0)
        } else {
            0.0
        };
        claims
            .insert_claim(claim(n_key.clone(), n, n + 0.2))
            .unwrap();
        claims
            .insert_claim(claim(k_key.clone(), k, if lossy { k + 0.03 } else { 0.0 }))
            .unwrap();
        bindings.push(SpectralComplexIndexKnot {
            wavelength_nm,
            real_index_property: n_key,
            extinction_index_property: k_key,
        });
    }
    if with_mu {
        claims
            .insert_claim(claim(PropertyKey::new(MU_NAME, MU), M, 2.0 * M))
            .unwrap();
    }
    let card = card(claims);
    let source =
        SpectralComplexIndexStatePoint::try_from_resolved(resolve(&card), &bindings).unwrap();
    (card, source)
}
fn si_grid_card(dielectric: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (i, micrometres) in [8.0, 8.5, 9.0, 10.0, 11.0, 11.5, 12.0, 13.0, 14.0]
        .into_iter()
        .enumerate()
    {
        let wavelength = micrometres * 1e-6;
        claims
            .insert_claim(claim(
                PropertyKey::new(
                    SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
                    Dims([1, 0, 0, 0, 0, 0]),
                ),
                wavelength,
                1.1 * wavelength,
            ))
            .unwrap();
        let n_name = if dielectric {
            SPECTRAL_DIELECTRIC_IOR_PROPERTIES[i]
        } else {
            SPECTRAL_CONDUCTOR_ETA_PROPERTIES[i]
        };
        let n = 1.5 + 0.01 * (micrometres - 8.0);
        claims
            .insert_claim(claim(PropertyKey::new(n_name, Dims::NONE), n, n + 0.2))
            .unwrap();
        let (loss_name, quantity, loss, hot) = if dielectric {
            let alpha = 3000.0 + 1000.0 * (micrometres - 8.0);
            (
                SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
                SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
                alpha,
                alpha + 2000.0,
            )
        } else {
            let k = 0.05 + 0.01 * (micrometres - 8.0);
            (
                SPECTRAL_CONDUCTOR_K_PROPERTIES[i],
                SPECTRAL_CONDUCTOR_K_QUANTITY,
                k,
                k + 0.03,
            )
        };
        claims
            .insert_claim(claim(
                PropertyKey::with_quantity(loss_name, quantity),
                loss,
                hot,
            ))
            .unwrap();
    }
    claims
        .insert_claim(claim(PropertyKey::new(MU_NAME, MU), M, 2.0 * M))
        .unwrap();
    card(claims)
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 8e-13 * (1.0 + b.abs()), "{a:e} != {b:e}");
}
fn close_complex(a: C64, b: C64) {
    assert!((a - b).abs() <= 8e-13 * (1.0 + b.abs()), "{a:?} != {b:?}");
}
fn exposure_limits() -> AngularSpectrumExposureLimits {
    AngularSpectrumExposureLimits {
        maximum_lines: 4,
        maximum_point_line_pairs: 16,
        maximum_point_frequency_pairs: 16,
        maximum_workspace_bytes: 1 << 20,
        maximum_fft_work: 1024,
    }
}
fn beam<'s, 'm>(source: &'s PlanarDielectricLayers<'m>) -> PlanarAngularSpectrum<'s, 'm> {
    PlanarAngularSpectrum::prepare(
        source,
        AngularSpectrumGrid {
            width: 1,
            height: 1,
            spacing_m: [20e-6; 2],
        },
        &[[C64::ONE, C64::ZERO]],
        AngularSpectrumLimits {
            maximum_grid_points: 1,
            maximum_workspace_bytes: 1 << 20,
            maximum_fft_work: 1024,
        },
        || false,
    )
    .unwrap()
}

#[test]
fn finite_band_generic_complex_source_drives_original_fresnel_and_refreshes_mu() {
    let (card, source) = generic_source(true, true);
    assert_eq!(source.samples().len(), 3); // Arbitrary source knots, no nine-point resampling.
    let response =
        SpectralOpticalResponse::try_new(SpectralOpticalSource::ComplexIndex(&source), 4.1e-6)
            .unwrap();
    assert!(core::ptr::eq(response.material(), source.resolved()));
    assert!(source.resolved().property(EPS_NAME).is_none());
    let vacuum = material(1.0 / (M * C0 * C0), M);
    let interface =
        PlanarDielectricInterface::try_new(&vacuum, DielectricResponse::SpectralOptical(response))
            .unwrap();
    for wavelength in [8000.0, 9300.0, 12_000.0, 14_000.0] {
        let omega = WL_FREQUENCY_NM / wavelength;
        let index = C64::new(
            1.4 + 1e-5 * (wavelength - 8000.0),
            -(0.05 + 1e-5 * (wavelength - 8000.0)),
        );
        close_complex(response.complex_index(omega).unwrap(), index);
        let reflected = (C64::ONE - index) / (C64::ONE + index);
        for polarization in [Polarization::Te, Polarization::Tm] {
            let observed = interface.response(omega, 1.0, polarization).unwrap();
            close_complex(observed.reflected_tangential_e, reflected);
            close(observed.reflected_power_fraction, reflected.norm_sq());
        }
        let epsilon = response.permittivity_f_m(omega).unwrap();
        close_complex(epsilon.scale(M * C0 * C0), index * index);
    }
    for omega in [
        (WL_FREQUENCY_NM / 14_000.0).next_down(),
        (WL_FREQUENCY_NM / 8000.0).next_up(),
        WL_FREQUENCY_NM / 550.0,
    ] {
        assert!(response.complex_index(omega).is_err());
    }
    for gap in [0.0, f64::NAN, 2e-6] {
        assert!(SpectralOpticalResponse::try_new(response.source(), gap).is_err());
    }
    let (_, missing_mu) = generic_source(true, false);
    assert!(
        SpectralOpticalResponse::try_new(SpectralOpticalSource::ComplexIndex(&missing_mu), 4.1e-6)
            .is_err()
    );
    let refreshed = response.at_query_point(&card, &point(600.0)).unwrap();
    let warm = refreshed.response().unwrap();
    close(warm.permeability_h_m() / M, 2.0);
    let index = C64::new(1.62, -0.10);
    close_complex(
        warm.complex_index(WL_FREQUENCY_NM / 10_000.0).unwrap(),
        index,
    );
    close_complex(
        warm.permittivity_f_m(WL_FREQUENCY_NM / 10_000.0)
            .unwrap()
            .scale(2.0 * M * C0 * C0),
        index * index,
    );
    close(response.permeability_h_m() / M, 1.0);
    assert!(response.at_query_point(&card, &point(601.0)).is_err());
}

#[test]
fn transparent_infrared_incidence_retains_snell_film_phase_and_beam_source_identity() {
    let (_, source) = generic_source(false, true);
    let optical =
        SpectralOpticalResponse::try_new(SpectralOpticalSource::ComplexIndex(&source), 4.1e-6)
            .unwrap();
    let response = DielectricResponse::SpectralOptical(optical);
    let vacuum = material(1.0 / (M * C0 * C0), M);
    let interface = PlanarDielectricInterface::try_new_with_incident_response(
        response,
        DielectricResponse::Lossless(&vacuum),
    )
    .unwrap();
    let thickness = 1e-6;
    let wavelengths = [9000.0, 12_000.0];
    let stacks = wavelengths.map(|wavelength| {
        PlanarDielectricLayers::try_new_with_incident_response(
            response,
            response,
            WL_FREQUENCY_NM / wavelength,
            &[DielectricLayer {
                response,
                thickness_m: thickness,
            }],
        )
        .unwrap()
    });
    let beams = stacks.each_ref().map(beam);
    for (i, wavelength) in wavelengths.into_iter().enumerate() {
        let n = 1.4 + 1e-5 * (wavelength - 8000.0);
        let cosine: f64 = 0.8;
        let exit_cosine = (1.0 - n * n * (1.0 - cosine * cosine)).sqrt();
        let expected = (n * cosine - exit_cosine) / (n * cosine + exit_cosine);
        close(
            interface
                .response(WL_FREQUENCY_NM / wavelength, cosine, Polarization::Te)
                .unwrap()
                .reflected_power_fraction,
            expected * expected,
        );
        let phase = core::f64::consts::TAU * n * thickness / (wavelength * 1e-9);
        let expected = C64::new(phase.cos(), -phase.sin());
        let transmitted = stacks[i].response(1.0, Polarization::Te).unwrap();
        close(transmitted.reflected_power_fraction, 0.0);
        close(transmitted.absorbed_power_fraction, 0.0);
        close_complex(transmitted.transmitted_tangential_e, expected);
        let field = beams[i]
            .field(AngularSpectrumSide::Transmitted, 0.0, || false)
            .unwrap();
        close_complex(field.electric_v_per_m()[0][0], expected);
    }
    let bank =
        AngularSpectrumExposure::try_new(&beams[0], &beams, exposure_limits(), || false).unwrap();
    assert!(
        bank.prepare_layer(
            0,
            AngularSpectrumExposureLayerLimits {
                maximum_frequency_groups: 4,
                maximum_frequency_layer_mode_pairs: 16,
                maximum_point_line_pairs: 16,
                maximum_point_reservoir_frequency_pairs: 64,
                maximum_workspace_bytes: 2 << 20,
                maximum_fft_work: 1024,
            },
            || false
        )
        .is_err()
    ); // Optical table data do not supply causal heat reservoirs.
    let foreign_source = source.clone();
    let foreign_response = DielectricResponse::SpectralOptical(
        SpectralOpticalResponse::try_new(
            SpectralOpticalSource::ComplexIndex(&foreign_source),
            4.1e-6,
        )
        .unwrap(),
    );
    let foreign_stack = PlanarDielectricLayers::try_new_with_incident_response(
        foreign_response,
        foreign_response,
        WL_FREQUENCY_NM / 12_000.0,
        &[DielectricLayer {
            response: foreign_response,
            thickness_m: thickness,
        }],
    )
    .unwrap();
    let foreign = [beam(&stacks[0]), beam(&foreign_stack)];
    assert!(
        AngularSpectrumExposure::try_new(&foreign[0], &foreign, exposure_limits(), || false)
            .is_err()
    );
}

#[test]
fn dielectric_interpolates_original_power_absorption_and_refreshes_the_si_grid() {
    let card = si_grid_card(true);
    let source = SpectralTabulatedDielectricStatePoint::try_from_material(resolve(&card)).unwrap();
    let response = SpectralOpticalResponse::try_new(
        SpectralOpticalSource::TabulatedDielectric(&source),
        1.2e-6,
    )
    .unwrap();
    let wavelength = 9.25e-6;
    let alpha = 4250.0;
    let k = alpha * wavelength / (4.0 * core::f64::consts::PI);
    close_complex(
        response.complex_index(WL_FREQUENCY_M / wavelength).unwrap(),
        C64::new(1.5125, -k),
    );
    let interpolated_k =
        (0.75 * 4000.0 * 9e-6 + 0.25 * 5000.0 * 10e-6) / (4.0 * core::f64::consts::PI);
    assert!(
        (interpolated_k - k).abs() > 1e-6,
        "interpolating k would change source alpha"
    );
    let limited = SpectralOpticalResponse::try_new(response.source(), 1.05e-6).unwrap();
    assert!(limited.at_query_point(&card, &point(600.0)).is_err()); // Warm gap exceeds the SAME policy.
    let refreshed = response.at_query_point(&card, &point(600.0)).unwrap();
    let warm = refreshed.response().unwrap();
    assert!(warm.complex_index(WL_FREQUENCY_M / 8.4e-6).is_err());
    let wavelength = 1.1 * 9.25e-6;
    let alpha = 6250.0;
    let index = C64::new(1.7125, -alpha * wavelength / (4.0 * core::f64::consts::PI));
    close_complex(
        warm.complex_index(WL_FREQUENCY_M / wavelength).unwrap(),
        index,
    );
    close_complex(
        warm.permittivity_f_m(WL_FREQUENCY_M / wavelength)
            .unwrap()
            .scale(2.0 * M * C0 * C0),
        index * index,
    );
    assert!(response.complex_index(WL_FREQUENCY_M / 8.4e-6).is_ok());
    for omega in [
        (WL_FREQUENCY_M / (14.0 * 1e-6)).next_down(),
        (WL_FREQUENCY_M / (8.0 * 1e-6)).next_up(),
    ] {
        assert!(response.complex_index(omega).is_err());
    }
}

#[test]
fn si_conductor_refresh_keeps_source_wavelength_roles_distinct_from_fixed_bindings() {
    let card = si_grid_card(false);
    let source = SpectralConductorStatePoint::try_from_material(resolve(&card)).unwrap();
    let response =
        SpectralOpticalResponse::try_new(SpectralOpticalSource::Conductor(&source), 1.2e-6)
            .unwrap();
    let fixed = SpectralOpticalResponse::try_new(
        SpectralOpticalSource::ComplexIndex(source.complex_index()),
        1.2e-6,
    )
    .unwrap();
    let omega = WL_FREQUENCY_NM / 9250.0;
    close_complex(
        response.complex_index(omega).unwrap(),
        fixed.complex_index(omega).unwrap(),
    );
    let refreshed = response.at_query_point(&card, &point(600.0)).unwrap();
    let fixed_refreshed = fixed.at_query_point(&card, &point(600.0)).unwrap();
    assert!(
        refreshed
            .response()
            .unwrap()
            .complex_index(WL_FREQUENCY_NM / 8400.0)
            .is_err()
    );
    assert!(
        fixed_refreshed
            .response()
            .unwrap()
            .complex_index(WL_FREQUENCY_NM / 8400.0)
            .is_ok()
    );
    close_complex(
        refreshed
            .response()
            .unwrap()
            .complex_index(WL_FREQUENCY_NM / (1.1 * 9250.0))
            .unwrap(),
        C64::new(1.7125, -0.0925),
    );
    assert!(response.complex_index(WL_FREQUENCY_NM / 8400.0).is_ok());
}
