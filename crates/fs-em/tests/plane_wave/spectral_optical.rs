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

#[test]
fn absorbing_crystal_slab_retains_original_sources_and_both_faces() {
    use fs_em::plane_wave::{
        UniaxialOpticalResponse,
        uniaxial::oriented::{
            OrientedUniaxialInterface,
            slab::{OrientedUniaxialSlab, SlabRegion},
        },
    };

    let mut claims = ClaimSet::new();
    for (prefix, indices, extinction) in [
        ("ordinary", [1.4, 1.6], 0.08),
        ("extraordinary", [1.8, 2.0], 0.2),
    ] {
        for (i, index) in indices.into_iter().enumerate() {
            claims
                .insert_claim(claim(
                    PropertyKey::new(format!("{prefix}_n_{i}"), Dims::NONE),
                    index,
                    index,
                ))
                .unwrap();
            claims
                .insert_claim(claim(
                    PropertyKey::with_quantity(
                        format!("{prefix}_k_{i}"),
                        SPECTRAL_CONDUCTOR_K_QUANTITY,
                    ),
                    extinction,
                    extinction,
                ))
                .unwrap();
        }
    }
    claims
        .insert_claim(claim(PropertyKey::new(MU_NAME, MU), M, M))
        .unwrap();
    let original_card = card(claims);
    let original = resolve(&original_card);
    let bindings = |prefix: &str| {
        [0usize, 1].map(|i| SpectralComplexIndexKnot {
            wavelength_nm: [8000.0, 14_000.0][i],
            real_index_property: PropertyKey::new(format!("{prefix}_n_{i}"), Dims::NONE),
            extinction_index_property: PropertyKey::with_quantity(
                format!("{prefix}_k_{i}"),
                SPECTRAL_CONDUCTOR_K_QUANTITY,
            ),
        })
    };
    let ordinary =
        SpectralComplexIndexStatePoint::try_from_resolved(original.clone(), &bindings("ordinary"))
            .unwrap();
    let extraordinary =
        SpectralComplexIndexStatePoint::try_from_resolved(original, &bindings("extraordinary"))
            .unwrap();
    let crystal = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, 6.1e-6).unwrap();
    let outside = material(1.0 / (M * C0 * C0), M);
    let exit = material(1.2 * 1.2 / (0.8 * M * C0 * C0), 0.8 * M);
    let entrance = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [1.0, 2.0, 3.0],
    )
    .unwrap();
    let slab = OrientedUniaxialSlab::try_new(&entrance, DielectricResponse::Lossless(&exit), 5e-6)
        .unwrap();
    let omega = WL_FREQUENCY_NM / 10_000.0;
    let electric = [C64::new(0.6, -0.2), C64::new(0.1, 0.7)];
    let expected_loss = crystal
        .principal_permittivities_f_m(omega)
        .unwrap()
        .map(|epsilon| -epsilon.im);
    assert!(expected_loss.iter().all(|loss| *loss > 0.0));

    // Exercise both public entry points, including a non-axis-aligned
    // incidence plane with global Cartesian incident polarization.
    let prepared_fields = [
        slab.prepare_absorbing_field(omega, 0.8, electric).unwrap(),
        slab.prepare_absorbing_field_at_tangential_index_in_plane(omega, 0.4, [3.0, 4.0], electric)
            .unwrap(),
    ];
    for prepared in &prepared_fields {
        assert!(core::ptr::eq(prepared.source(), &slab));
        assert!(core::ptr::eq(
            prepared.source().entrance_interface(),
            &entrance
        ));
        assert!(core::ptr::eq(
            prepared.source().exit_response().material(),
            &exit
        ));
        let retained_crystal = prepared.source().entrance_interface().crystal_response();
        assert!(core::ptr::eq(retained_crystal.ordinary_source(), &ordinary));
        assert!(core::ptr::eq(
            retained_crystal.extraordinary_source(),
            &extraordinary
        ));
        assert_eq!(prepared.response().omega_rad_s(), omega);
        assert_eq!(
            prepared.response().principal_loss_permittivities_f_m(),
            expected_loss
        );
        assert!(prepared.response().absorbed_power_fraction > 0.0);
        assert!(
            prepared
                .field_at(SlabRegion::Crystal, [0.0, 0.0, 0.5 * slab.thickness_m()])
                .unwrap()
                .absorbed_power_density_w_m3
                > 0.0
        );

        let area = 0.025;
        let depth = slab.thickness_m();
        let heat = prepared.slab_absorption(area, 0.0, depth).unwrap();
        assert!(core::ptr::eq(heat.prepared_field(), prepared));
        assert!(core::ptr::eq(heat.source(), &slab));
        assert_eq!(heat.area_m2(), area);
        assert_eq!(heat.start_depth_m(), 0.0);
        assert_eq!(heat.end_depth_m(), depth);
        assert_eq!(heat.volume_m3(), area * depth);
        assert!(heat.absorbed_power_w() > 0.0);
        let split = 0.37 * depth;
        let left = prepared.slab_absorption(area, 0.0, split).unwrap();
        let right = prepared.slab_absorption(area, split, depth).unwrap();
        near(
            left.absorbed_power_w() + right.absorbed_power_w(),
            heat.absorbed_power_w(),
            5e-10 * heat.absorbed_power_w(),
        );
        assert!(
            prepared
                .slab_absorption(area, -f64::MIN_POSITIVE, depth)
                .is_err()
        );
        assert!(
            prepared
                .slab_absorption(area, 0.0, f64::from_bits(depth.to_bits() + 1))
                .is_err()
        );
        assert!(prepared.slab_absorption(area, depth, depth).is_err());

        for (z, left, right) in [
            (0.0, SlabRegion::TotalIncidentSide, SlabRegion::Crystal),
            (
                slab.thickness_m(),
                SlabRegion::Crystal,
                SlabRegion::Transmitted,
            ),
        ] {
            let point = [0.3e-6, -0.2e-6, z];
            let a = prepared.field_at(left, point).unwrap();
            let b = prepared.field_at(right, point).unwrap();
            for j in 0..2 {
                complex_near(a.field.electric_v_m[j], b.field.electric_v_m[j], 3e-12);
                complex_near(a.field.magnetic_a_m[j], b.field.magnetic_a_m[j], 3e-14);
            }
            near(a.field.poynting_w_m2[2], b.field.poynting_w_m2[2], 3e-14);
        }
    }
    // Opposite-face drives retain the same material snapshots and use authored
    // global Cartesian phasors at each physical face, including a rotated plane.
    let right_electric = [C64::new(0.2, 0.1), C64::new(-0.3, 0.4)];
    let both_sides = [
        slab.prepare_absorbing_field_from_both_sides(omega, 0.8, electric, right_electric)
            .unwrap(),
        slab.prepare_absorbing_field_from_both_sides_at_tangential_index(
            omega,
            0.4,
            electric,
            right_electric,
        )
        .unwrap(),
        slab.prepare_absorbing_field_from_both_sides_at_tangential_index_in_plane(
            omega,
            0.4,
            [3.0, 4.0],
            electric,
            right_electric,
        )
        .unwrap(),
    ];
    for prepared in &both_sides {
        assert!(core::ptr::eq(prepared.source(), &slab));
        assert!(core::ptr::eq(
            prepared.source().entrance_interface(),
            &entrance
        ));
        assert_eq!(
            prepared.response().principal_loss_permittivities_f_m(),
            expected_loss
        );
        assert!(prepared.response().incident_normal_power_w_m2 > 0.0);
        assert!(prepared.response().right_incident_normal_power_w_m2 > 0.0);
        let depth = slab.thickness_m();
        let incoming = prepared.exit_incident_field_at([0.0, 0.0, depth]).unwrap();
        for j in 0..2 {
            complex_near(
                prepared.response().incident_electric_v_m[j],
                electric[j],
                3e-14,
            );
            complex_near(incoming.field.electric_v_m[j], right_electric[j], 3e-14);
        }
        assert!(incoming.field.poynting_w_m2[2] < 0.0);
        let inside = prepared
            .field_at(SlabRegion::Crystal, [0.0, 0.0, depth])
            .unwrap();
        let outside = prepared.total_exit_field_at([0.0, 0.0, depth]).unwrap();
        for j in 0..2 {
            complex_near(
                inside.field.electric_v_m[j],
                outside.field.electric_v_m[j],
                3e-12,
            );
            complex_near(
                inside.field.magnetic_a_m[j],
                outside.field.magnetic_a_m[j],
                3e-14,
            );
        }
        let heat = prepared.slab_absorption(0.025, 0.0, depth).unwrap();
        assert!(core::ptr::eq(heat.prepared_field(), prepared));
        near(
            heat.absorbed_power_w(),
            0.025
                * prepared.response().total_incident_normal_power_w_m2()
                * prepared.response().absorbed_power_fraction,
            5e-10 * heat.absorbed_power_w(),
        );
    }
    let right_only = slab
        .prepare_absorbing_field_from_both_sides(omega, 0.8, [C64::ZERO; 2], right_electric)
        .unwrap();
    assert_eq!(right_only.response().incident_normal_power_w_m2, 0.0);
    assert!(right_only.response().right_incident_normal_power_w_m2 > 0.0);
    assert!(
        slab.prepare_absorbing_field_from_both_sides(omega, 0.8, [C64::ZERO; 2], [C64::ZERO; 2],)
            .is_err()
    );

    // Preparing a film must still obey the original finite spectral support.
    assert!(
        slab.prepare_absorbing_field(WL_FREQUENCY_NM / 7000.0, 0.8, electric)
            .is_err()
    );
    check_absorbing_film_illumination(&slab, omega, electric);
    check_absorbing_crystal_stack(&slab, omega, electric);
}

fn check_absorbing_crystal_stack(
    slab: &fs_em::plane_wave::uniaxial::oriented::slab::OrientedUniaxialSlab<'_, '_>,
    omega: f64,
    electric: [C64; 2],
) {
    use fs_em::plane_wave::uniaxial::oriented::{
        absorbing::stack::StackRegion,
        slab::SlabRegion,
        stack::{OrientedUniaxialStack, UniaxialStackLayer, UniaxialStackLimits},
    };
    let interface = slab.entrance_interface();
    let layers = [UniaxialStackLayer::try_new(
        interface.crystal_response(),
        interface.optic_axis(),
        slab.thickness_m() * 0.5,
    )
    .unwrap(); 2];
    let stack = OrientedUniaxialStack::try_new(
        interface.incident_response(),
        slab.exit_response(),
        &layers,
        UniaxialStackLimits {
            maximum_layers: 2,
            maximum_workspace_bytes: 1 << 20,
        },
    )
    .unwrap();
    let prepared = stack
        .prepare_absorbing_field(omega, 0.8, electric, || false)
        .unwrap();
    let single = slab.prepare_absorbing_field(omega, 0.8, electric).unwrap();
    assert!(core::ptr::eq(prepared.source(), &stack));
    assert!(core::ptr::eq(
        prepared.source().layers()[1].crystal_response().material(),
        interface.crystal_response().material(),
    ));
    close(
        prepared.response().reflected_power_fraction,
        single.response().reflected_power_fraction,
    );
    close(
        prepared.response().transmitted_power_fraction,
        single.response().transmitted_power_fraction,
    );
    for (layer, fraction) in [(0, 0.2), (1, 0.7)] {
        let point = [0.3e-6, -0.2e-6, fraction * slab.thickness_m()];
        let a = prepared.field_at(StackRegion::Layer(layer), point).unwrap();
        let b = single.field_at(SlabRegion::Crystal, point).unwrap();
        for j in 0..3 {
            close_complex(a.field.electric_v_m[j], b.field.electric_v_m[j]);
            close_complex(a.field.magnetic_a_m[j], b.field.magnetic_a_m[j]);
        }
        close(a.absorbed_power_density_w_m3, b.absorbed_power_density_w_m3);
    }
    let area = 0.7e-4;
    let mut watts = 0.0;
    for (index, layer) in layers.iter().enumerate() {
        let volume = prepared
            .layer_absorption(index, area, 0.0, layer.thickness_m())
            .unwrap();
        assert!(core::ptr::eq(volume.source(), &stack));
        assert_eq!(volume.layer_index(), index);
        watts += volume.absorbed_power_w();
    }
    close(
        watts,
        single
            .slab_absorption(area, 0.0, slab.thickness_m())
            .unwrap()
            .absorbed_power_w(),
    );

    // The same original global axes also survive an arbitrary incidence plane.
    let plane = [0.6, 0.8];
    let direct = stack
        .prepare_absorbing_field_in_plane(omega, 0.8, plane, electric, || false)
        .unwrap();
    let cached = stack
        .prepare_absorbing_polarization(omega, 0.8, plane, || false)
        .unwrap();
    assert!(core::ptr::eq(cached.source(), &stack));
    close(
        cached.incident_normal_power_w_m2(electric).unwrap(),
        direct.response().incident_normal_power_w_m2,
    );
    for (index, layer) in layers.iter().enumerate() {
        close(
            cached
                .layer_absorbed_power_w(index, electric, area, 0.0, layer.thickness_m())
                .unwrap(),
            direct
                .layer_absorption(index, area, 0.0, layer.thickness_m())
                .unwrap()
                .absorbed_power_w(),
        );
    }
    assert!(prepared.layer_absorption(2, area, 0.0, 1e-6).is_err());
    assert!(
        prepared
            .layer_absorption(0, area, 0.0, slab.thickness_m())
            .is_err()
    );
    assert!(
        stack
            .prepare_absorbing_field(omega, 0.8, electric, || true)
            .is_err()
    );
    assert!(
        stack
            .prepare_absorbing_field(WL_FREQUENCY_NM / 7000.0, 0.8, electric, || false)
            .is_err()
    );
}

// G0/G1/G3: source-bound positive factors preserve film interference, original
// tensor work and the distinction between coherent and independent light.
fn check_absorbing_film_illumination(
    slab: &fs_em::plane_wave::uniaxial::oriented::slab::OrientedUniaxialSlab<'_, '_>,
    omega: f64,
    electric: [C64; 2],
) {
    use fs_em::plane_wave::uniaxial::oriented::{
        absorbing::{
            CrystalIlluminationError,
            slab::{AbsorbingSlabCoherentIllumination as Group, AbsorbingSlabIllumination as Bank},
        },
        slab::{OrientedUniaxialSlab, SlabRegion},
    };
    let fields = [
        slab.prepare_absorbing_field(omega, 0.8, electric).unwrap(),
        slab.prepare_absorbing_field(omega, 0.8, electric).unwrap(),
        slab.prepare_absorbing_field(omega, 0.8, electric.map(|value| -value))
            .unwrap(),
    ];
    let constructive = Group::try_new(slab, &fields[..2], 2, || false).unwrap();
    let destructive = Group::try_new(slab, &fields[1..], 2, || false).unwrap();
    let independent = [
        Group::try_new(slab, &fields[1..2], 1, || false).unwrap(),
        Group::try_new(slab, &fields[2..], 1, || false).unwrap(),
    ];
    let bank = Bank::try_new(slab, &independent, 2, 2, || false).unwrap();
    assert!(core::ptr::eq(constructive.source(), slab));
    assert!(core::ptr::eq(bank.source(), slab));
    assert_eq!(constructive.omega_rad_s(), Some(omega));
    assert_eq!(constructive.component_count(), 2);
    assert_eq!(bank.group_count(), 2);
    assert_eq!(bank.component_count(), 2);

    for (region, z) in [
        (SlabRegion::Incident, -0.4e-6),
        (SlabRegion::Reflected, -0.4e-6),
        (SlabRegion::TotalIncidentSide, -0.4e-6),
        (SlabRegion::Crystal, 0.37 * slab.thickness_m()),
        (SlabRegion::Transmitted, 1.4 * slab.thickness_m()),
    ] {
        let point = [0.3e-6, -0.2e-6, z];
        let single = fields[0].field_at(region, point).unwrap();
        let bright = constructive.field_at(region, point, || false).unwrap();
        let dark = destructive.field_at(region, point, || false).unwrap();
        let mixed = bank.field_statistics(region, point, || false).unwrap();
        let flux = bank.poynting_w_m2(region, point, || false).unwrap();
        for i in 0..3 {
            complex_near(
                bright.field.electric_v_m[i],
                single.field.electric_v_m[i].scale(2.0),
                3e-12,
            );
            complex_near(
                bright.field.magnetic_a_m[i],
                single.field.magnetic_a_m[i].scale(2.0),
                3e-14,
            );
            assert!(dark.field.electric_v_m[i].abs() < 3e-12);
            assert!(dark.field.magnetic_a_m[i].abs() < 3e-14);
            near(
                bright.field.poynting_w_m2[i],
                4.0 * single.field.poynting_w_m2[i],
                3e-13,
            );
            near(
                mixed.poynting_w_m2[i],
                2.0 * single.field.poynting_w_m2[i],
                3e-13,
            );
            assert_eq!(flux[i], mixed.poynting_w_m2[i]);
            near(
                mixed.magnetic_rms_squared_a2_m2[i],
                single.field.magnetic_a_m[i].norm_sq(),
                3e-14,
            );
            for j in 0..3 {
                complex_near(
                    mixed.electric_coherency_v2_m2[i][j],
                    single.field.electric_v_m[i] * single.field.electric_v_m[j].conj(),
                    3e-12,
                );
            }
        }
        let heat = single.absorbed_power_density_w_m3;
        near(
            bright.absorbed_power_density_w_m3,
            4.0 * heat,
            3e-11 * heat.max(1.0),
        );
        near(
            mixed.absorbed_power_density_w_m3,
            2.0 * heat,
            3e-11 * heat.max(1.0),
        );
        assert!(dark.absorbed_power_density_w_m3 < 1e-20 * heat.max(1.0));
    }

    // Equal transverse p/s powers are unpolarized at oblique incidence. A
    // complex unitary factor change to circular components preserves all
    // statistics after the anisotropic film, without cohering the factors.
    let p = [C64::from_re(0.8), C64::ZERO];
    let s = [C64::ZERO, C64::ONE];
    let linear = [p, s].map(|e| slab.prepare_absorbing_field(omega, 0.8, e).unwrap());
    let circular = [[p[0], C64::new(0.0, 1.0)], [C64::new(0.0, 0.8), s[1]]].map(|e| {
        slab.prepare_absorbing_field(
            omega,
            0.8,
            e.map(|value| value.scale(core::f64::consts::FRAC_1_SQRT_2)),
        )
        .unwrap()
    });
    let linear_groups: [_; 2] =
        core::array::from_fn(|i| Group::try_new(slab, &linear[i..i + 1], 1, || false).unwrap());
    let circular_groups: [_; 2] =
        core::array::from_fn(|i| Group::try_new(slab, &circular[i..i + 1], 1, || false).unwrap());
    let linear_bank = Bank::try_new(slab, &linear_groups, 2, 2, || false).unwrap();
    let circular_bank = Bank::try_new(slab, &circular_groups, 2, 2, || false).unwrap();
    for z in [0.0, 0.3 * slab.thickness_m(), slab.thickness_m()] {
        let a = linear_bank
            .field_statistics(SlabRegion::Crystal, [0.0, 0.0, z], || false)
            .unwrap();
        let b = circular_bank
            .field_statistics(SlabRegion::Crystal, [0.0, 0.0, z], || false)
            .unwrap();
        for i in 0..3 {
            for j in 0..3 {
                complex_near(
                    a.electric_coherency_v2_m2[i][j],
                    b.electric_coherency_v2_m2[i][j],
                    4e-12,
                );
            }
            near(
                a.magnetic_rms_squared_a2_m2[i],
                b.magnetic_rms_squared_a2_m2[i],
                4e-14,
            );
            near(a.poynting_w_m2[i], b.poynting_w_m2[i], 4e-13);
        }
        near(
            a.absorbed_power_density_w_m3,
            b.absorbed_power_density_w_m3,
            4e-11 * a.absorbed_power_density_w_m3.max(1.0),
        );
    }

    // Different directions and incidence planes retain their global phase and
    // the same fixed optic axis. Contract the total E independently in that
    // original loss tensor, instead of adding independently positive heat.
    let crossed = [
        slab.prepare_absorbing_field_at_tangential_index_in_plane(omega, 0.4, [3.0, 4.0], electric)
            .unwrap(),
        slab.prepare_absorbing_field_at_tangential_index_in_plane(
            omega,
            -0.3,
            [4.0, -3.0],
            electric.map(|value| value * C64::new(0.3, 0.7)),
        )
        .unwrap(),
    ];
    let group = Group::try_new(slab, &crossed, 2, || false).unwrap();
    let axis = slab.entrance_interface().optic_axis();
    let loss = crossed[0].response().principal_loss_permittivities_f_m();
    let mut interferes = false;
    for x in [-2e-6, 0.0, 3e-6] {
        let point = [x, 0.7e-6, 0.6 * slab.thickness_m()];
        let a = crossed[0].field_at(SlabRegion::Crystal, point).unwrap();
        let b = crossed[1].field_at(SlabRegion::Crystal, point).unwrap();
        let actual = group
            .field_at(SlabRegion::Crystal, point, || false)
            .unwrap();
        let e: [C64; 3] =
            core::array::from_fn(|i| a.field.electric_v_m[i] + b.field.electric_v_m[i]);
        let dot = (0..3).fold(C64::ZERO, |sum, i| sum + e[i].scale(axis[i]));
        let norm = axis.iter().map(|v| v * v).sum::<f64>();
        let parallel = dot.norm_sq() / norm;
        let perpendicular = e.iter().map(|value| value.norm_sq()).sum::<f64>() - parallel;
        let expected_heat = 0.5 * omega * (loss[0] * perpendicular + loss[1] * parallel);
        near(
            actual.absorbed_power_density_w_m3,
            expected_heat,
            4e-11 * expected_heat.max(1.0),
        );
        for i in 0..3 {
            complex_near(actual.field.electric_v_m[i], e[i], 3e-12);
        }
        interferes |=
            (expected_heat - a.absorbed_power_density_w_m3 - b.absorbed_power_density_w_m3).abs()
                > 1e-3 * expected_heat;
    }
    assert!(interferes);

    let other_frequency = [slab
        .prepare_absorbing_field(f64::from_bits(omega.to_bits() + 1), 0.8, electric)
        .unwrap()];
    let mismatched = [
        slab.prepare_absorbing_field(omega, 0.8, electric).unwrap(),
        slab.prepare_absorbing_field(f64::from_bits(omega.to_bits() + 1), 0.8, electric)
            .unwrap(),
    ];
    assert!(Group::try_new(slab, &mismatched, 2, || false).is_err());
    let spectral_groups = [
        Group::try_new(slab, &fields[..1], 1, || false).unwrap(),
        Group::try_new(slab, &other_frequency, 1, || false).unwrap(),
    ];
    let spectral = Bank::try_new(slab, &spectral_groups, 2, 2, || false).unwrap();
    let point = [0.1e-6, 0.2e-6, 0.5 * slab.thickness_m()];
    let expected = fields[0]
        .field_at(SlabRegion::Crystal, point)
        .unwrap()
        .absorbed_power_density_w_m3
        + other_frequency[0]
            .field_at(SlabRegion::Crystal, point)
            .unwrap()
            .absorbed_power_density_w_m3;
    near(
        spectral
            .field_statistics(SlabRegion::Crystal, point, || false)
            .unwrap()
            .absorbed_power_density_w_m3,
        expected,
        3e-11 * expected,
    );

    let other_slab = OrientedUniaxialSlab::try_new(
        slab.entrance_interface(),
        slab.exit_response(),
        slab.thickness_m(),
    )
    .unwrap();
    assert!(Group::try_new(&other_slab, &fields[..1], 1, || false).is_err());
    assert!(Bank::try_new(&other_slab, &independent, 2, 2, || false).is_err());
    assert!(Group::try_new(slab, &fields[..2], 1, || false).is_err());
    assert!(Bank::try_new(slab, &independent, 1, 2, || false).is_err());
    assert!(Bank::try_new(slab, &independent, 2, 1, || false).is_err());
    assert!(matches!(
        bank.field_statistics(SlabRegion::Crystal, point, || true),
        Err(CrystalIlluminationError::Cancelled)
    ));
    assert!(matches!(
        constructive.field_at(SlabRegion::Crystal, point, || true),
        Err(CrystalIlluminationError::Cancelled)
    ));
    assert!(
        bank.field_statistics(SlabRegion::Crystal, point, || false)
            .is_ok()
    );
    let empty = Group::try_new(slab, &[], 0, || false).unwrap();
    let dark = Bank::try_new(slab, &[], 0, 0, || false).unwrap();
    assert_eq!(empty.omega_rad_s(), None);
    assert_eq!(
        empty
            .field_at(SlabRegion::Crystal, point, || false)
            .unwrap()
            .absorbed_power_density_w_m3,
        0.0
    );
    assert_eq!(
        dark.field_statistics(SlabRegion::Crystal, point, || false)
            .unwrap()
            .electric_coherency_v2_m2,
        [[C64::ZERO; 3]; 3]
    );
    for (region, invalid) in [
        (SlabRegion::Crystal, [0.0, 0.0, -1.0]),
        (SlabRegion::Transmitted, [0.0; 3]),
        (SlabRegion::Incident, [0.0, 0.0, 1.0]),
        (SlabRegion::Crystal, [f64::NAN, 0.0, 0.0]),
    ] {
        assert!(empty.field_at(region, invalid, || false).is_err());
        assert!(dark.field_statistics(region, invalid, || false).is_err());
        assert!(dark.poynting_w_m2(region, invalid, || false).is_err());
    }
}

#[test]
fn g1_bulk_spectral_complex_index_and_si_sources_keep_original_wavelength_derivatives() {
    let (_, generic) = generic_source(true, true);
    let conductor_card = si_grid_card(false);
    let conductor =
        SpectralConductorStatePoint::try_from_material(resolve(&conductor_card)).unwrap();
    let dielectric_card = si_grid_card(true);
    let dielectric =
        SpectralTabulatedDielectricStatePoint::try_from_material(resolve(&dielectric_card))
            .unwrap();
    for source in [
        SpectralOpticalSource::ComplexIndex(&generic),
        SpectralOpticalSource::Conductor(&conductor),
        SpectralOpticalSource::TabulatedDielectric(&dielectric),
    ] {
        let law = SpectralOpticalResponse::try_new(source, 6.1e-6).unwrap();
        let wave = DielectricResponse::SpectralOptical(law);
        for wavelength in [8.25e-6, 9.25e-6, 12.5e-6] {
            let omega = WL_FREQUENCY_M / wavelength;
            let jet = wave.bulk_dispersion(omega).unwrap();
            assert!(core::ptr::eq(jet.source().material(), law.material()));
            assert_eq!(jet.permittivity_f_m(), law.permittivity_f_m(omega).unwrap());
            let group = if matches!(source, SpectralOpticalSource::ComplexIndex(_)) {
                1.32
            } else {
                1.42
            };
            assert!((jet.group_slowness_s_m() / (group / C0) - 1.0).abs() < 8e-15);
            assert_eq!(jet.group_velocity_dispersion_s2_m(), 0.0);
            if matches!(source, SpectralOpticalSource::TabulatedDielectric(_)) {
                // Original alpha slope is 1000 per micrometre, not k slope.
                let first = wavelength * 1e9 / (2.0 * omega);
                let second = -(wavelength * 1e9 / omega) / omega;
                assert!((jet.wavenumber_derivative_s_m().im / first - 1.0).abs() < 8e-15);
                assert!((jet.wavenumber_second_derivative_s2_m().im / second - 1.0).abs() < 8e-15);
            } else {
                // The authored linear extinction has intercept -0.03.
                assert!((jet.wavenumber_derivative_s_m().im / (0.03 / C0) - 1.0).abs() < 2e-14);
                assert_eq!(jet.wavenumber_second_derivative_s2_m(), C64::ZERO);
            }
        }
        for omega in [WL_FREQUENCY_M / 8e-6, WL_FREQUENCY_M / 14e-6] {
            assert!(wave.bulk_dispersion(omega).is_err());
        }
    }
}
