//! Synthetic source controls: these test propagation and source binding, not
//! measured glass qualification or a fitted temperature-dependent loss model.
use super::*;
use crate::plane_wave::{
    DielectricResponse, PlanarDielectricInterface, Polarization,
    angular_spectrum::{
        AngularSpectrumExposure, AngularSpectrumExposureLayerLimits, AngularSpectrumExposureLimits,
        AngularSpectrumGrid, AngularSpectrumLimits, PlanarAngularSpectrum,
    },
    dispersion::IncidenceFrequencyConstraint,
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    air_refraction::CiddorAir,
    catalog_refraction::{AirRelativeSellmeier, CatalogWavelength, CiddorAirState},
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
        MAGNETIC_PERMEABILITY_PROPERTY,
    },
    state_point::{
        MaterialPropertySelection, SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
        SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY, SPECTRAL_DIELECTRIC_IOR_PROPERTIES,
        SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES,
        SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES, ScalarAdmissibility, ScalarPropertyRequirement,
        SpectralSellmeierStatePoint, SpectralTabulatedDielectricStatePoint,
        VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES, VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES,
        resolve_material_state_point,
    },
    thermo_optic::ThermoOpticCoefficients,
};
use fs_qty::Dims;

const MU: f64 = 1.256_637_062_12e-6;
const B: [f64; 3] = [1.1, 0.2, 0.8];
const C: [f64; 3] = [0.01, 0.04, 100.0];

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn complex_near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}"
    );
}
fn claim(name: &str, quantity: QuantitySpec, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::with_quantity(name, quantity),
        value: PropertyValue::Scalar {
            value,
            dims: quantity.dims(),
        },
        validity: ValidityDomain::unconstrained(),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic thermal Maxwell glass".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn assemble(properties: Vec<PropertyClaim>, revision: u32) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored thermal glass".into(),
            phase: "solid".into(),
            process: "equation and transport control".into(),
            revision,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new().with("T", temperature).unwrap()
}
fn resolve(card: &MaterialCard, point: &QueryPoint) -> ResolvedMaterialStatePoint {
    let requirements: Vec<_> = card
        .claims()
        .claims_ordered()
        .map(|(_, claim)| {
            ScalarPropertyRequirement::try_with_key(&claim.key, ScalarAdmissibility::Finite)
                .unwrap()
        })
        .collect();
    resolve_material_state_point(
        card,
        point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn reference_phase(wavelength_nm: f64) -> f64 {
    let square = (wavelength_nm * 1e-3).powi(2);
    (1.0 + B
        .into_iter()
        .zip(C)
        .map(|(b, c)| b * square / (square - c))
        .sum::<f64>())
    .sqrt()
}
fn thermal_phase(wavelength_nm: f64, temperature: f64) -> f64 {
    let n = reference_phase(wavelength_nm);
    let dt = temperature - 300.0;
    let square = (wavelength_nm * 1e-3).powi(2);
    n + (n * n - 1.0) / (2.0 * n)
        * (1.2e-5 * dt - 2e-8 * dt.powi(2)
            + 3e-11 * dt.powi(3)
            + (4e-6 * dt - 5e-9 * dt.powi(2)) / (square - 0.18_f64.powi(2)))
}
fn glass_card(alpha: f64, revision: u32) -> MaterialCard {
    glass_card_spectrum([alpha; 9], revision)
}
fn glass_card_spectrum(alpha: [f64; 9], revision: u32) -> MaterialCard {
    let dimensionless = QuantitySpec::dimensional(Dims::NONE);
    let length = QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0]));
    let square = QuantitySpec::dimensional(Dims([2, 0, 0, 0, 0, 0]));
    let mut properties = Vec::new();
    for i in 0..3 {
        properties.push(claim(
            VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[i],
            dimensionless,
            B[i],
        ));
        properties.push(claim(
            VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES[i],
            square,
            C[i] * 1e-12,
        ));
    }
    for (name, bound) in SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES
        .into_iter()
        .zip([800e-9, 2000e-9])
    {
        properties.push(claim(name, length, bound));
    }
    for i in 0..9 {
        let wavelength = 800.0 + 150.0 * i as f64;
        properties.push(claim(
            SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
            length,
            wavelength * 1e-9,
        ));
        properties.push(claim(
            SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
            SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
            alpha[i],
        ));
        properties.push(claim(
            SPECTRAL_DIELECTRIC_IOR_PROPERTIES[i],
            dimensionless,
            reference_phase(wavelength) + if i == 4 { 0.01 } else { 0.0 },
        ));
    }
    let mut magnetic = claim(
        MAGNETIC_PERMEABILITY_PROPERTY,
        QuantitySpec::dimensional(MAGNETIC_PERMEABILITY_DIMS),
        MU,
    );
    magnetic.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(200.0, MU), (500.0, 2.0 * MU)],
        dims: MAGNETIC_PERMEABILITY_DIMS,
    };
    magnetic.validity = ValidityDomain::unconstrained().with("T", 200.0, 500.0);
    magnetic.interpolation = InterpolationPolicy::LinearInside;
    properties.push(magnetic);
    assemble(properties, revision)
}
fn law() -> SchottThermoOpticLaw {
    SchottThermoOpticLaw::try_new(
        ThermoOpticCoefficients {
            d: [1.2e-5, -2e-8, 3e-11],
            e: [4e-6, -5e-9],
            resonance_um: 0.18,
        },
        300.0,
        [200.0, 500.0],
        [800.0, 2000.0],
        glass_card(0.0, 0).content_hash(),
    )
    .unwrap()
}
fn port() -> ResolvedMaterialStatePoint {
    let card = assemble(
        vec![
            claim(
                ELECTRIC_PERMITTIVITY_PROPERTY,
                QuantitySpec::dimensional(ELECTRIC_PERMITTIVITY_DIMS),
                1.0 / (MU * VACUUM_SPEED_M_S.powi(2)),
            ),
            claim(
                MAGNETIC_PERMEABILITY_PROPERTY,
                QuantitySpec::dimensional(MAGNETIC_PERMEABILITY_DIMS),
                MU,
            ),
        ],
        0,
    );
    resolve(&card, &QueryPoint::new())
}
fn reference<'a>(source: &'a SpectralSellmeierStatePoint) -> SpectralOpticalResponse<'a> {
    SpectralOpticalResponse::try_new(SpectralOpticalSource::DielectricSellmeier(source), 2e-7)
        .unwrap()
}

#[test]
fn actual_temperature_and_mu_reach_phase_group_gvd_without_refitting() {
    let card = glass_card(0.0, 0);
    let source =
        SpectralSellmeierStatePoint::try_from_material(resolve(&card, &point(300.0))).unwrap();
    let reference = reference(&source);
    let law = law();
    for temperature in [210.0, 337.5, 480.0] {
        let actual = resolve(&card, &point(temperature));
        let response =
            ThermoOpticOpticalResponse::try_lossless(reference, &law, &actual, "T").unwrap();
        assert!(core::ptr::eq(response.material(), &actual));
        assert!(core::ptr::eq(
            response.reference().material(),
            source.resolved()
        ));
        assert_ne!(
            response.material().identity(),
            response.reference().material().identity()
        );
        for wavelength in [850.0_f64, 1550.0, 1950.0] {
            let omega = FREQUENCY_WAVELENGTH_NM / wavelength;
            let n = thermal_phase(wavelength, temperature);
            near(response.complex_index(omega).unwrap().re, n, 1e-15);
            near(
                response.permittivity_f_m(omega).unwrap().re
                    * response.permeability_h_m()
                    * VACUUM_SPEED_M_S.powi(2),
                n * n,
                4e-15,
            );
            let h = wavelength * 1e-3;
            let f =
                [-2.0, -1.0, 0.0, 1.0, 2.0].map(|k| thermal_phase(wavelength + k * h, temperature));
            let first = wavelength * (f[0] - 8.0 * f[1] + 8.0 * f[3] - f[4]) / (12.0 * h);
            let second = wavelength.powi(2)
                * (-f[4] + 16.0 * f[3] - 30.0 * f[2] + 16.0 * f[1] - f[0])
                / (12.0 * h * h);
            let bulk = DielectricResponse::ThermoOpticOptical(response)
                .bulk_dispersion(omega)
                .unwrap();
            near(
                bulk.group_slowness_s_m() * VACUUM_SPEED_M_S,
                n - first,
                1e-9,
            );
            near(
                bulk.group_velocity_dispersion_s2_m() * omega * VACUUM_SPEED_M_S,
                second,
                2e-8,
            );
            assert_eq!(bulk.wavenumber_derivative_s_m().im, 0.0);
            assert_eq!(bulk.wavenumber_second_derivative_s2_m().im, 0.0);
            near(
                bulk.propagation(0.002).unwrap().group_delay_s,
                0.002 * (n - first) / VACUUM_SPEED_M_S,
                1e-20,
            );
        }
        let restored = response.at_state(source.resolved(), "T").unwrap();
        let omega = FREQUENCY_WAVELENGTH_NM / 1550.0;
        assert_eq!(
            restored.complex_index(omega).unwrap(),
            reference.complex_index(omega).unwrap()
        );
    }
}

#[test]
fn thermal_interfaces_frequency_jets_and_coherent_films_follow_original_equation() {
    let card = glass_card(0.0, 0);
    let source =
        SpectralSellmeierStatePoint::try_from_material(resolve(&card, &point(300.0))).unwrap();
    let law = law();
    let ambient = port();
    let omega = FREQUENCY_WAVELENGTH_NM / 1550.0;
    let ci = 0.8_f64;
    let thickness = 120e-6;
    let mut transmitted = Vec::new();
    for temperature in [250.0, 450.0] {
        let actual = resolve(&card, &point(temperature));
        let response =
            ThermoOpticOpticalResponse::try_lossless(reference(&source), &law, &actual, "T")
                .unwrap();
        let interface = PlanarDielectricInterface::try_new(
            &ambient,
            DielectricResponse::ThermoOpticOptical(response),
        )
        .unwrap();
        let stack = PlanarDielectricLayers::try_new(
            &ambient,
            DielectricResponse::Lossless(&ambient),
            omega,
            &[DielectricLayer {
                response: DielectricResponse::ThermoOpticOptical(response),
                thickness_m: thickness,
            }],
        )
        .unwrap();
        for polarization in [Polarization::Te, Polarization::Tm] {
            let fresnel = |frequency: f64| {
                let n = thermal_phase(FREQUENCY_WAVELENGTH_NM / frequency, temperature);
                let q = (n * n - (1.0 - ci * ci)).sqrt();
                let yi = match polarization {
                    Polarization::Te => ci,
                    Polarization::Tm => ci.recip(),
                } / MU;
                let yt = match polarization {
                    Polarization::Te => q,
                    Polarization::Tm => n * n / q,
                } / response.permeability_h_m();
                ((yi - yt) / (yi + yt), q)
            };
            let (r, q) = fresnel(omega);
            let carrier = interface.response(omega, ci, polarization).unwrap();
            complex_near(carrier.reflected_tangential_e, C64::from_re(r), 3e-15);
            near(
                carrier.transmitted_normal_wavenumber_rad_m.re,
                omega / VACUUM_SPEED_M_S * q,
                5e-9,
            );
            near(
                carrier.reflected_power_fraction + carrier.transmitted_power_fraction,
                1.0,
                4e-15,
            );
            let jet = interface
                .frequency_dispersion(
                    omega,
                    ci,
                    polarization,
                    IncidenceFrequencyConstraint::FixedAngle,
                )
                .unwrap();
            let h = 1e-3;
            let f = [-2.0, -1.0, 0.0, 1.0, 2.0].map(|k| fresnel(omega * (1.0 + k * h)).0);
            complex_near(
                jet.reflected_tangential_e().scaled_derivatives()[0],
                C64::from_re((f[0] - 8.0 * f[1] + 8.0 * f[3] - f[4]) / (12.0 * h)),
                1e-9,
            );
            complex_near(
                jet.reflected_tangential_e().scaled_derivatives()[1],
                C64::from_re(
                    (-f[4] + 16.0 * f[3] - 30.0 * f[2] + 16.0 * f[1] - f[0]) / (12.0 * h * h),
                ),
                2e-8,
            );
            let angle = -omega / VACUUM_SPEED_M_S * q * thickness;
            let p = C64::new(angle.cos(), angle.sin());
            let denominator = C64::ONE - (p * p).scale(r * r);
            let film = stack.response(ci, polarization).unwrap();
            complex_near(
                film.reflected_tangential_e,
                (C64::ONE - p * p).scale(r) / denominator,
                3e-13,
            );
            complex_near(
                film.transmitted_tangential_e,
                p.scale(1.0 - r * r) / denominator,
                3e-13,
            );
            near(
                film.reflected_power_fraction + film.transmitted_power_fraction,
                1.0,
                4e-15,
            );
            assert_eq!(film.absorbed_power_fraction, 0.0);
            transmitted.push(film.transmitted_power_fraction);
        }
    }
    assert!((transmitted[0] - transmitted[2]).abs() > 1e-3);
}

#[test]
fn source_temperature_loss_band_and_original_table_cusps_are_preserved() {
    let card = glass_card(0.0, 0);
    let base = resolve(&card, &point(300.0));
    let source = SpectralSellmeierStatePoint::try_from_material(base.clone()).unwrap();
    let law = law();
    let actual = resolve(&card, &point(337.5));
    let response =
        ThermoOpticOpticalResponse::try_lossless(reference(&source), &law, &actual, "T").unwrap();
    for edge in response.frequency_range_rad_s() {
        assert!(response.complex_index(edge).is_ok());
        assert!(response.refractive_index_jet(edge).is_err());
    }
    let [lo, hi] = response.frequency_range_rad_s();
    for omega in [lo.next_down(), hi.next_up(), 0.0, f64::NAN] {
        assert!(response.complex_index(omega).is_err());
    }
    let wrong_reference =
        SpectralSellmeierStatePoint::try_from_material(resolve(&card, &point(301.0))).unwrap();
    assert!(
        ThermoOpticOpticalResponse::try_lossless(reference(&wrong_reference), &law, &actual, "T")
            .is_err()
    );
    let other = resolve(&glass_card(0.0, 1), &point(337.5));
    assert!(response.at_state(&other, "T").is_err());
    assert!(response.at_state(&actual, "wrong-axis").is_err());
    let pressure = resolve(&card, &point(337.5).with("pressure", 101325.0).unwrap());
    assert!(response.at_state(&pressure, "T").is_err());
    let absorbing_card = glass_card(0.1, 0);
    let absorbing =
        SpectralSellmeierStatePoint::try_from_material(resolve(&absorbing_card, &point(300.0)))
            .unwrap();
    let absorbing_actual = resolve(&absorbing_card, &point(337.5));
    assert!(
        ThermoOpticOpticalResponse::try_lossless(
            reference(&absorbing),
            &law,
            &absorbing_actual,
            "T"
        )
        .is_err()
    );
    let table = SpectralTabulatedDielectricStatePoint::try_from_material(base).unwrap();
    let table_response =
        SpectralOpticalResponse::try_new(SpectralOpticalSource::TabulatedDielectric(&table), 2e-7)
            .unwrap();
    let thermal_table =
        ThermoOpticOpticalResponse::try_lossless(table_response, &law, &actual, "T").unwrap();
    let cusp = core::f64::consts::TAU * VACUUM_SPEED_M_S / table.samples()[4].wavelength_m;
    assert!(thermal_table.complex_index(cusp).is_ok());
    assert!(thermal_table.refractive_index_jet(cusp).is_err());
}

#[test]
fn catalog_air_stays_fixed_and_native_beams_retain_complete_thermal_sources() {
    let card = assemble(
        vec![claim(
            MAGNETIC_PERMEABILITY_PROPERTY,
            QuantitySpec::dimensional(MAGNETIC_PERMEABILITY_DIMS),
            MU,
        )],
        0,
    );
    let magnetic = resolve(&card, &QueryPoint::new());
    let actual = resolve(&card, &point(337.5));
    let changed = resolve(&card, &point(338.5));
    let air = CiddorAirState::try_new(CiddorAir::try_new(0.01, 650.0).unwrap(), 293.15, 101330.0)
        .unwrap();
    let coordinate_air =
        CiddorAirState::try_new(CiddorAir::try_new(0.0, 450.0).unwrap(), 288.15, 101325.0).unwrap();
    let catalog = AirRelativeSellmeier::try_new(
        B,
        C,
        (800.0, 2000.0),
        air,
        CatalogWavelength::Air(coordinate_air),
    )
    .unwrap();
    let reference = SpectralOpticalResponse::try_new(
        SpectralOpticalSource::LosslessAirRelativeSellmeier {
            law: &catalog,
            magnetic: &magnetic,
        },
        2e-7,
    )
    .unwrap();
    let law = law();
    let response = ThermoOpticOpticalResponse::try_lossless(reference, &law, &actual, "T").unwrap();
    let omega = FREQUENCY_WAVELENGTH_NM / 1550.0;
    let original = catalog.refractive_index_jet(1550.0).unwrap();
    let expected = law
        .refractive_index_jet(
            1550.0,
            337.5,
            original.phase_index,
            original.group_index,
            original.scaled_wavelength_derivatives[1],
        )
        .unwrap();
    let observed = response.refractive_index_jet(omega).unwrap();
    near(observed.phase_index, expected.phase_index, 1e-15);
    near(observed.group_index, expected.group_index, 2e-15);
    near(
        observed.scaled_wavelength_derivatives[1],
        expected.scaled_wavelength_derivatives[1],
        2e-15,
    );
    assert_eq!(catalog.reference_air(), air);
    assert_eq!(
        catalog.wavelength_coordinate(),
        CatalogWavelength::Air(coordinate_air)
    );
    let ambient = port();
    let make = |response, wavelength| {
        PlanarDielectricLayers::try_new(
            &ambient,
            DielectricResponse::Lossless(&ambient),
            FREQUENCY_WAVELENGTH_NM / wavelength,
            &[DielectricLayer {
                response: DielectricResponse::ThermoOpticOptical(response),
                thickness_m: 1e-6,
            }],
        )
        .unwrap()
    };
    let stacks = [make(response, 1100.0), make(response, 1550.0)];
    assert!(stacks[0].same_sources_and_geometry(&stacks[1]));
    assert!(
        !stacks[0]
            .same_sources_and_geometry(&make(response.at_state(&changed, "T").unwrap(), 1550.0))
    );
    let copied_law = law;
    let other =
        ThermoOpticOpticalResponse::try_lossless(reference, &copied_law, &actual, "T").unwrap();
    assert!(!stacks[0].same_sources_and_geometry(&make(other, 1550.0)));
    let grid = AngularSpectrumGrid {
        width: 1,
        height: 1,
        spacing_m: [1e-6; 2],
    };
    let beams: Vec<_> = stacks
        .iter()
        .map(|stack| {
            PlanarAngularSpectrum::prepare(
                stack,
                grid,
                &[[C64::ONE, C64::ZERO]],
                AngularSpectrumLimits {
                    maximum_grid_points: 1,
                    maximum_workspace_bytes: 2_000_000,
                    maximum_fft_work: 100_000,
                },
                || false,
            )
            .unwrap()
        })
        .collect();
    let exposure = AngularSpectrumExposure::try_new(
        &beams[0],
        &beams,
        AngularSpectrumExposureLimits {
            maximum_lines: 2,
            maximum_point_line_pairs: 8,
            maximum_point_frequency_pairs: 8,
            maximum_workspace_bytes: 2_000_000,
            maximum_fft_work: 100_000,
        },
        || false,
    )
    .unwrap();
    assert_eq!(exposure.frequency_count(), 2);
    match exposure.prepare_layer(
        0,
        AngularSpectrumExposureLayerLimits {
            maximum_frequency_groups: 2,
            maximum_frequency_layer_mode_pairs: 8,
            maximum_point_line_pairs: 8,
            maximum_point_reservoir_frequency_pairs: 100,
            maximum_workspace_bytes: 2_000_000,
            maximum_fft_work: 100_000,
        },
        || false,
    ) {
        Err(error) => assert!(
            error
                .to_string()
                .contains("original causal material reservoir"),
            "{error}"
        ),
        Ok(_) => panic!("a finite-band thermal glass fit does not supply a causal reservoir"),
    }
}

#[test]
fn retained_reference_loss_keeps_intensity_absorption_and_full_frequency_jets() {
    let wavelength_nm = 1525.0;
    let wavelength_m = wavelength_nm * 1e-9;
    let omega = FREQUENCY_WAVELENGTH_NM / wavelength_nm;
    for slope_per_m2 in [0.0, 2e9] {
        let alpha_at = |lambda_m: f64| 1000.0 + slope_per_m2 * lambda_m;
        let card = glass_card_spectrum(
            core::array::from_fn(|i| alpha_at((800.0 + 150.0 * i as f64) * 1e-9)),
            0,
        );
        let original =
            SpectralSellmeierStatePoint::try_from_material(resolve(&card, &point(300.0))).unwrap();
        let reference = reference(&original);
        let law = law();
        for temperature in [210.0, 450.0] {
            let actual = resolve(&card, &point(temperature));
            assert!(
                ThermoOpticOpticalResponse::try_lossless(reference, &law, &actual, "T").is_err()
            );
            let response = ThermoOpticOpticalResponse::try_with_reference_absorption(
                reference, &law, &actual, "T",
            )
            .unwrap();
            assert!(response.retains_reference_absorption());
            let index = response.complex_index(omega).unwrap();
            assert_eq!(index.im, reference.complex_index(omega).unwrap().im);
            near(index.re, thermal_phase(wavelength_nm, temperature), 1e-15);
            let bulk = DielectricResponse::ThermoOpticOptical(response)
                .bulk_dispersion(omega)
                .unwrap();
            near(
                -2.0 * bulk.wavenumber_rad_m().im,
                alpha_at(wavelength_m),
                3e-11,
            );
            let group_loss = slope_per_m2 * wavelength_m.powi(2) / (4.0 * core::f64::consts::PI);
            near(
                bulk.wavenumber_derivative_s_m().im * VACUUM_SPEED_M_S,
                group_loss,
                2e-17,
            );
            near(
                bulk.wavenumber_second_derivative_s2_m().im * omega * VACUUM_SPEED_M_S,
                -2.0 * group_loss,
                4e-17,
            );
            if slope_per_m2 == 0.0 {
                assert_eq!(bulk.wavenumber_derivative_s_m().im, 0.0);
                assert_eq!(bulk.wavenumber_second_derivative_s2_m().im, 0.0);
            }
            // Differentiate an independent N(T,lambda)^2 equation. Both real/
            // imaginary cross terms must reach the native permittivity jet.
            let normalized_epsilon = |frequency: f64| {
                let lambda_nm = FREQUENCY_WAVELENGTH_NM / frequency;
                let lambda_m = lambda_nm * 1e-9;
                let n = C64::new(
                    thermal_phase(lambda_nm, temperature),
                    -alpha_at(lambda_m) * lambda_m / (4.0 * core::f64::consts::PI),
                );
                n * n
            };
            let h = 1e-3;
            let f = [-2.0, -1.0, 0.0, 1.0, 2.0].map(|k| normalized_epsilon(omega * (1.0 + k * h)));
            let first = (f[0] - f[1].scale(8.0) + f[3].scale(8.0) - f[4]).scale(1.0 / (12.0 * h));
            let second = (-f[4] + f[3].scale(16.0) - f[2].scale(30.0) + f[1].scale(16.0) - f[0])
                .scale(1.0 / (12.0 * h * h));
            let jet = response.frequency_jet(omega).unwrap();
            let scale = response.permeability_h_m() * VACUUM_SPEED_M_S.powi(2);
            complex_near(jet.epsilon_derivatives[0].scale(scale), first, 2e-9);
            complex_near(jet.epsilon_derivatives[1].scale(scale), second, 2e-8);
            let distance = 2e-4;
            near(
                bulk.propagation(distance).unwrap().transfer.norm_sq(),
                (-alpha_at(wavelength_m) * distance).exp(),
                3e-15,
            );
        }
    }
}

#[test]
fn absorbing_thermal_films_conserve_power_and_heat_the_native_interior_and_beam() {
    let alpha = 80_000.0;
    let card = glass_card(alpha, 0);
    let original =
        SpectralSellmeierStatePoint::try_from_material(resolve(&card, &point(300.0))).unwrap();
    let law = law();
    let actual = resolve(&card, &point(337.5));
    let response = ThermoOpticOpticalResponse::try_with_reference_absorption(
        reference(&original),
        &law,
        &actual,
        "T",
    )
    .unwrap();
    let ambient = port();
    let wavelength_nm = 1550.0;
    let omega = FREQUENCY_WAVELENGTH_NM / wavelength_nm;
    let k0 = omega / VACUUM_SPEED_M_S;
    let ci = 0.8_f64;
    let si = (1.0 - ci * ci).sqrt();
    let thickness = 2e-6;
    let n = C64::new(
        thermal_phase(wavelength_nm, 337.5),
        -alpha * wavelength_nm * 1e-9 / (4.0 * core::f64::consts::PI),
    );
    let square = n * n - C64::from_re(si * si);
    let q_re = (0.5 * (square.abs() + square.re)).sqrt();
    let q = C64::new(q_re, square.im / (2.0 * q_re));
    let propagate = |depth: f64| {
        let angle = -k0 * q.re * depth;
        C64::new(angle.cos(), angle.sin()).scale((k0 * q.im * depth).exp())
    };
    let p = propagate(thickness);
    let stack = PlanarDielectricLayers::try_new(
        &ambient,
        DielectricResponse::Lossless(&ambient),
        omega,
        &[DielectricLayer {
            response: DielectricResponse::ThermoOpticOptical(response),
            thickness_m: thickness,
        }],
    )
    .unwrap();
    let interface = PlanarDielectricInterface::try_new(
        &ambient,
        DielectricResponse::ThermoOpticOptical(response),
    )
    .unwrap();
    for polarization in [Polarization::Te, Polarization::Tm] {
        let yi = match polarization {
            Polarization::Te => ci,
            Polarization::Tm => ci.recip(),
        } / (MU * VACUUM_SPEED_M_S);
        let yt = match polarization {
            Polarization::Te => q,
            Polarization::Tm => n * n / q,
        }
        .scale(1.0 / (response.permeability_h_m() * VACUUM_SPEED_M_S));
        let r = (C64::from_re(yi) - yt) / (C64::from_re(yi) + yt);
        let denominator = C64::ONE - r * r * p * p;
        let reflected = r * (C64::ONE - p * p) / denominator;
        let transmitted = p * (C64::ONE - r * r) / denominator;
        complex_near(
            interface
                .response(omega, ci, polarization)
                .unwrap()
                .reflected_tangential_e,
            r,
            5e-15,
        );
        let observed = stack.response(ci, polarization).unwrap();
        complex_near(observed.reflected_tangential_e, reflected, 3e-14);
        complex_near(observed.transmitted_tangential_e, transmitted, 3e-14);
        near(
            observed.reflected_power_fraction
                + observed.transmitted_power_fraction
                + observed.absorbed_power_fraction,
            1.0,
            4e-15,
        );
        let expected_heat = 0.5 * yi * (1.0 - reflected.norm_sq() - transmitted.norm_sq());
        assert!(expected_heat > 0.0);
        let field = stack.prepare_layer_field(0, ci, polarization).unwrap();
        near(field.absorbed_power_w_m2(), expected_heat, 3e-17);
        let amplitude = (C64::ONE + r) / denominator;
        for depth in [0.0, 0.3 * thickness, thickness] {
            let forward = amplitude * propagate(depth);
            let backward = -r * amplitude * p * p / propagate(depth);
            let sample = field.field(depth).unwrap();
            complex_near(sample.electric_v_per_m, forward + backward, 5e-14);
            complex_near(
                sample.conjugate_magnetic_a_per_m,
                yt * (forward - backward),
                3e-16,
            );
        }
        // Independent local epsilon loss, including longitudinal TM electric
        // field. Its depth integral must equal the native positive film work.
        let epsilon = response.permittivity_f_m(omega).unwrap();
        let mut heat_integral = 0.0;
        let intervals = 1024;
        for i in 0..=intervals {
            let sample = field
                .field(thickness * i as f64 / intervals as f64)
                .unwrap();
            let longitudinal = match polarization {
                Polarization::Te => C64::ZERO,
                Polarization::Tm => {
                    (sample.conjugate_magnetic_a_per_m / epsilon).scale(-k0 * si / omega)
                }
            };
            let heat = -0.5
                * omega
                * epsilon.im
                * (sample.electric_v_per_m.norm_sq() + longitudinal.norm_sq());
            let weight = if i == 0 || i == intervals {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            heat_integral += weight * heat;
        }
        heat_integral *= thickness / (3.0 * intervals as f64);
        near(
            field.absorbed_power_w_m2(),
            heat_integral,
            expected_heat * 5e-8,
        );
        let grid = AngularSpectrumGrid {
            width: 1,
            height: 1,
            spacing_m: [2e-6; 2],
        };
        let jones = match polarization {
            Polarization::Te => [C64::ZERO, C64::ONE],
            Polarization::Tm => [C64::ONE, C64::ZERO],
        };
        let beam = PlanarAngularSpectrum::prepare_with_carrier(
            &stack,
            grid,
            [k0 * si, 0.0],
            &[jones],
            AngularSpectrumLimits {
                maximum_grid_points: 1,
                maximum_workspace_bytes: 2_000_000,
                maximum_fft_work: 100_000,
            },
            || false,
        )
        .unwrap();
        let power = beam.power();
        near(power.film_absorbed_w, expected_heat * 4e-12, 2e-27);
        near(
            power.reflected_w / power.incident_w,
            reflected.norm_sq(),
            4e-14,
        );
        near(
            power.transmitted_w / power.incident_w,
            transmitted.norm_sq(),
            4e-14,
        );
        near(power.balance_defect_w(), 0.0, 2e-27);
    }
}

#[test]
fn absorption_policy_refresh_and_original_loss_cusps_remain_owned() {
    let card = glass_card(0.0, 0);
    let original =
        SpectralSellmeierStatePoint::try_from_material(resolve(&card, &point(300.0))).unwrap();
    let law = law();
    let actual = resolve(&card, &point(337.5));
    let zero =
        ThermoOpticOpticalResponse::try_lossless(reference(&original), &law, &actual, "T").unwrap();
    let held = ThermoOpticOpticalResponse::try_with_reference_absorption(
        reference(&original),
        &law,
        &actual,
        "T",
    )
    .unwrap();
    assert!(!zero.retains_reference_absorption());
    assert!(held.retains_reference_absorption());
    assert!(!zero.same_source(held));
    let mut absorption = [2000.0; 9];
    absorption[4] = 3000.0;
    let card = glass_card_spectrum(absorption, 0);
    let original =
        SpectralSellmeierStatePoint::try_from_material(resolve(&card, &point(300.0))).unwrap();
    let actual = resolve(&card, &point(337.5));
    let changed = resolve(&card, &point(450.0));
    let response = ThermoOpticOpticalResponse::try_with_reference_absorption(
        reference(&original),
        &law,
        &actual,
        "T",
    )
    .unwrap();
    let refreshed = response.at_state(&changed, "T").unwrap();
    assert!(refreshed.retains_reference_absorption());
    assert!(!response.same_source(refreshed));
    let cusp =
        core::f64::consts::TAU * VACUUM_SPEED_M_S / original.absorption_samples()[4].wavelength_m;
    let index = response.complex_index(cusp).unwrap();
    assert_eq!(
        index.im,
        reference(&original).complex_index(cusp).unwrap().im
    );
    assert_eq!(refreshed.complex_index(cusp).unwrap().im, index.im);
    assert_ne!(refreshed.complex_index(cusp).unwrap().re, index.re);
    assert!(response.frequency_jet(cusp).is_err());
    assert!(refreshed.frequency_jet(cusp).is_err());
    for edge in response.frequency_range_rad_s() {
        assert!(response.complex_index(edge).is_ok());
        assert!(response.frequency_jet(edge).is_err());
    }
    let omega = FREQUENCY_WAVELENGTH_NM / 1525.0;
    assert!(response.frequency_jet(omega).is_ok());
    assert!(refreshed.frequency_jet(omega).is_ok());
}
