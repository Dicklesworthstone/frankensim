//! G0/G1/G3 source-to-Maxwell controls; scalar material cards are test fixtures,
//! not measurements of permeability or empirical optical-device qualification.
use super::*;
use crate::plane_wave::{
    DielectricResponse, PlanarDielectricInterface, Polarization,
    angular_spectrum::{
        AngularSpectrumExposure, AngularSpectrumExposureLayerLimits, AngularSpectrumExposureLimits,
        AngularSpectrumGrid, AngularSpectrumLimits, PlanarAngularSpectrum,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS,
        MAGNETIC_PERMEABILITY_PROPERTY,
    },
    state_point::{
        MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
        resolve_material_state_point,
    },
};
use fs_qty::{Dims, semantic::SemanticType};

const MU: f64 = 1.256_637_062_12e-6;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}; tolerance {tolerance:.3e}"
    );
}
fn complex_near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}; tolerance {tolerance:.3e}"
    );
}
fn phase(angle: f64) -> C64 {
    let (sine, cosine) = angle.sin_cos();
    C64::new(cosine, sine)
}

fn snapshot(
    properties: &[(&str, Dims, f64)],
    point: QueryPoint,
    validity: ValidityDomain,
) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for &(name, dims, value) in properties {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: validity.clone(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored CHARMS-silica Maxwell control".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "caller-selected Corning 7980 source applicability".into(),
            phase: "homogeneous optical control".into(),
            process: "source-admission and analytic wave control".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &point,
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn silica_state(temperature: f64, mu: f64) -> ResolvedMaterialStatePoint {
    snapshot(
        &[(
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            mu,
        )],
        QueryPoint::new().with("T", temperature).unwrap(),
        ValidityDomain::unconstrained().with("T", 0.0, 1000.0),
    )
}
fn unit_index_port() -> ResolvedMaterialStatePoint {
    snapshot(
        &[
            (
                ELECTRIC_PERMITTIVITY_PROPERTY,
                ELECTRIC_PERMITTIVITY_DIMS,
                1.0 / (MU * VACUUM_SPEED_M_S.powi(2)),
            ),
            (
                MAGNETIC_PERMEABILITY_PROPERTY,
                MAGNETIC_PERMEABILITY_DIMS,
                MU,
            ),
        ],
        QueryPoint::new(),
        ValidityDomain::unconstrained(),
    )
}
fn admittance(polarization: Polarization, n: f64, q: f64, mu: f64) -> f64 {
    match polarization {
        Polarization::Te => q / (VACUUM_SPEED_M_S * mu),
        Polarization::Tm => n * n / (q * VACUUM_SPEED_M_S * mu),
    }
}

#[test]
fn source_indices_thermal_derivatives_and_si_permeability_reach_native_response() {
    // 70-digit Decimal controls of NASA Table 3, with explicit temperature
    // powers and numerical derivatives; these are fit checks, not measurements.
    let controls = [
        (
            30.0,
            400.0,
            1.468_993_543_706_87,
            1.512_181_464_956_940_7,
            3.381_490_402_126_881e-6,
        ),
        (
            222.25,
            1550.0,
            1.443_885_066_753_925_3,
            1.462_442_027_459_719_8,
            6.762_005_789_946_712e-6,
        ),
        (
            300.0,
            2600.0,
            1.428_277_318_914_038_4,
            1.479_727_471_218_672_5,
            9.457_631_451_406_03e-6,
        ),
    ];
    for (temperature, wavelength, phase_index, group_index, thermal_derivative) in controls {
        for mu in [MU, 2.3 * MU] {
            let source = silica_state(temperature, mu);
            let law = SilicaOpticalResponse::try_lossless(&source, "T").unwrap();
            assert!(core::ptr::eq(law.source(), &source));
            assert_eq!(law.material().identity(), source.identity());
            assert_eq!(law.temperature_k(), temperature);
            assert_eq!(law.permeability_h_m(), mu);
            assert_eq!(
                law.refractive_source().identity(),
                Corning7980Silica::at_temperature(temperature)
                    .unwrap()
                    .identity()
            );
            let omega = FREQUENCY_WAVELENGTH_NM / wavelength;
            let sample = law.refractive_sample(omega).unwrap();
            near(sample.phase_index, phase_index, 8e-15);
            near(sample.group_index, group_index, 2e-14);
            near(
                sample.derivative_temperature_per_k,
                thermal_derivative,
                2e-18,
            );
            let epsilon = law.permittivity_f_m(omega).unwrap();
            assert_eq!(epsilon.im, 0.0);
            near(
                epsilon.re * mu * VACUUM_SPEED_M_S.powi(2),
                phase_index.powi(2),
                3e-14,
            );
            assert_eq!(law.complex_index(omega).unwrap().im, 0.0);
        }
    }
}

#[test]
fn source_temperature_quantity_and_frequency_support_are_not_inferred() {
    for temperature in [0.0, 29.999, 300.001, 310.0] {
        assert!(SilicaOpticalResponse::try_lossless(&silica_state(temperature, MU), "T").is_err());
    }
    let properties = [(
        MAGNETIC_PERMEABILITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS,
        MU,
    )];
    for kind in [
        QuantityKind::AbsoluteTemperature,
        QuantityKind::TemperatureDifference,
    ] {
        let quantity = QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
        let source = snapshot(
            &properties,
            QueryPoint::new()
                .with_quantity("glass-temperature", quantity, 295.0)
                .unwrap(),
            ValidityDomain::unconstrained().with_quantity(
                "glass-temperature",
                quantity,
                290.0,
                300.0,
            ),
        );
        assert_eq!(
            SilicaOpticalResponse::try_lossless(&source, "glass-temperature").is_ok(),
            kind == QuantityKind::AbsoluteTemperature
        );
    }
    let no_temperature = snapshot(
        &properties,
        QueryPoint::new(),
        ValidityDomain::unconstrained(),
    );
    assert!(SilicaOpticalResponse::try_lossless(&no_temperature, "T").is_err());
    let missing_mu = snapshot(
        &[("other", Dims::NONE, 1.0)],
        QueryPoint::new().with("T", 295.0).unwrap(),
        ValidityDomain::unconstrained(),
    );
    assert!(SilicaOpticalResponse::try_lossless(&missing_mu, "T").is_err());
    for unknown in ["omega", "pressure", "stress"] {
        let source = snapshot(
            &properties,
            QueryPoint::new()
                .with("T", 295.0)
                .unwrap()
                .with(unknown, 1.0)
                .unwrap(),
            ValidityDomain::unconstrained()
                .with("T", 290.0, 300.0)
                .with(unknown, 0.0, 2.0),
        );
        assert!(SilicaOpticalResponse::try_lossless(&source, "T").is_err());
    }
    let source = silica_state(295.0, MU);
    assert!(SilicaOpticalResponse::try_lossless(&source, "different-axis").is_err());
    let law = SilicaOpticalResponse::try_lossless(&source, "T").unwrap();
    let low = FREQUENCY_WAVELENGTH_NM / WAVELENGTH_RANGE_NM[1];
    let high = FREQUENCY_WAVELENGTH_NM / WAVELENGTH_RANGE_NM[0];
    for omega in [low, high] {
        law.permittivity_f_m(omega).unwrap();
        assert!(DielectricResponse::SilicaOptical(law).bulk_dispersion(omega).is_err());
    }
    for omega in [0.0, -1.0, f64::NAN, f64::INFINITY, low.next_down(), high.next_up()] {
        assert!(law.permittivity_f_m(omega).is_err());
    }
}

#[test]
fn thermal_silica_interfaces_obey_snell_fresnel_brewster_and_total_reflection() {
    let ambient = unit_index_port();
    for temperature in [30.0, 222.25, 300.0] {
        let source = silica_state(temperature, MU);
        let silica = SilicaOpticalResponse::try_lossless(&source, "T").unwrap();
        let forward = PlanarDielectricInterface::try_new(
            &ambient,
            DielectricResponse::SilicaOptical(silica),
        )
        .unwrap();
        let reverse = PlanarDielectricInterface::try_new_with_incident_response(
            DielectricResponse::SilicaOptical(silica),
            DielectricResponse::Lossless(&ambient),
        )
        .unwrap();
        for wavelength in [400.0, 1064.0, 1550.0, 2600.0] {
            let omega = FREQUENCY_WAVELENGTH_NM / wavelength;
            let n = silica.complex_index(omega).unwrap().re;
            let k0 = omega / VACUUM_SPEED_M_S;
            for ci in [1.0_f64, 0.8, 0.35] {
                let qt = (n * n - (1.0 - ci * ci)).sqrt();
                for polarization in [Polarization::Te, Polarization::Tm] {
                    let yi = admittance(polarization, 1.0, ci, MU);
                    let yt = admittance(polarization, n, qt, MU);
                    let r = (yi - yt) / (yi + yt);
                    let t = 2.0 * yi / (yi + yt);
                    let actual = forward.response(omega, ci, polarization).unwrap();
                    complex_near(actual.reflected_tangential_e, C64::from_re(r), 3e-15);
                    complex_near(actual.transmitted_tangential_e, C64::from_re(t), 3e-15);
                    near(actual.transmitted_normal_wavenumber_rad_m.re / k0, qt, 3e-15);
                    near(actual.transmitted_power_fraction, yt / yi * t * t, 5e-15);
                    near(
                        actual.reflected_power_fraction + actual.transmitted_power_fraction,
                        1.0,
                        3e-15,
                    );
                    let back = reverse.response(omega, qt / n, polarization).unwrap();
                    complex_near(back.reflected_tangential_e, C64::from_re(-r), 5e-15);
                    near(back.transmitted_normal_wavenumber_rad_m.re / k0, ci, 5e-15);
                }
            }
            let brewster = forward
                .response(omega, 1.0 / (1.0 + n * n).sqrt(), Polarization::Tm)
                .unwrap();
            near(brewster.reflected_power_fraction, 0.0, 1e-29);
            let ci = 0.4_f64;
            let gamma = (n * n * (1.0 - ci * ci) - 1.0).sqrt();
            for polarization in [Polarization::Te, Polarization::Tm] {
                let reflected = reverse.response(omega, ci, polarization).unwrap();
                near(reflected.reflected_power_fraction, 1.0, 3e-15);
                assert_eq!(reflected.transmitted_power_fraction, 0.0);
                near(
                    reflected.transmitted_normal_wavenumber_rad_m.im / k0,
                    -gamma,
                    3e-15,
                );
            }
        }
    }
}

#[test]
fn original_temperature_changes_coherent_film_phase_and_preserves_power() {
    let ambient = unit_index_port();
    let omega = FREQUENCY_WAVELENGTH_NM / 1550.0;
    let ci = 0.8_f64;
    let thickness = 12e-6;
    let mut powers = Vec::new();
    for temperature in [30.0, 300.0] {
        let source = silica_state(temperature, MU);
        let silica = SilicaOpticalResponse::try_lossless(&source, "T").unwrap();
        let n = silica.complex_index(omega).unwrap().re;
        let q = (n * n - (1.0 - ci * ci)).sqrt();
        let kz = omega / VACUUM_SPEED_M_S * q;
        let p = phase(-kz * thickness);
        let stack = PlanarDielectricLayers::try_new(
            &ambient,
            DielectricResponse::Lossless(&ambient),
            omega,
            &[DielectricLayer {
                response: DielectricResponse::SilicaOptical(silica),
                thickness_m: thickness,
            }],
        )
        .unwrap();
        for polarization in [Polarization::Te, Polarization::Tm] {
            let yi = admittance(polarization, 1.0, ci, MU);
            let yt = admittance(polarization, n, q, MU);
            let r = (yi - yt) / (yi + yt);
            let denominator = C64::ONE - (p * p).scale(r * r);
            let expected_r = (C64::ONE - p * p).scale(r) / denominator;
            let expected_t = p.scale(1.0 - r * r) / denominator;
            let actual = stack.response(ci, polarization).unwrap();
            complex_near(actual.reflected_tangential_e, expected_r, 3e-14);
            complex_near(actual.transmitted_tangential_e, expected_t, 3e-14);
            near(
                actual.reflected_power_fraction + actual.transmitted_power_fraction,
                1.0,
                4e-15,
            );
            assert_eq!(actual.absorbed_power_fraction, 0.0);
            powers.push(actual.transmitted_power_fraction);
            let field = stack.prepare_layer_field(0, ci, polarization).unwrap();
            let amplitude = C64::from_re(1.0 + r) / denominator;
            for z in [0.0, thickness * 0.4, thickness] {
                let forward = amplitude * phase(-kz * z);
                let backward = (amplitude * p * p * phase(kz * z)).scale(-r);
                let observed = field.field(z).unwrap();
                complex_near(observed.electric_v_per_m, forward + backward, 8e-14);
                complex_near(
                    observed.conjugate_magnetic_a_per_m,
                    (forward - backward).scale(yt),
                    3e-16,
                );
            }
        }
    }
    assert!(
        (powers[0] - powers[2]).abs() > 1e-3,
        "original thermal index must move the TE fringe"
    );
}

#[test]
fn bulk_group_delay_and_gvd_retain_original_fit_at_fixed_temperature() {
    // Independent 80-digit Decimal, explicit Table-3 polynomial powers and
    // five-point wavelength differences at 1e-8 um; units cancel in lambda² n''.
    for (temperature, wavelength, group, curvature) in [
        (50.0, 450.0, 1.498_295_167_660_902, 0.102_502_342_050_541_27),
        (
            222.25,
            1550.0,
            1.462_442_027_459_719_8,
            -0.010_357_736_539_282_036,
        ),
        (
            295.0,
            2500.0,
            1.477_252_812_050_778_2,
            -0.058_970_901_178_658_13,
        ),
    ] {
        for mu in [MU, 2.3 * MU] {
            let source = silica_state(temperature, mu);
            let silica = SilicaOpticalResponse::try_lossless(&source, "T").unwrap();
            let omega = FREQUENCY_WAVELENGTH_NM / wavelength;
            let jet = DielectricResponse::SilicaOptical(silica)
                .bulk_dispersion(omega)
                .unwrap();
            assert!(core::ptr::eq(jet.source().material(), &source));
            near(jet.group_slowness_s_m() * VACUUM_SPEED_M_S, group, 3e-14);
            near(
                jet.group_velocity_dispersion_s2_m() * (VACUUM_SPEED_M_S * omega),
                curvature,
                3e-15,
            );
            assert_eq!(jet.wavenumber_derivative_s_m().im, 0.0);
            assert_eq!(jet.wavenumber_second_derivative_s2_m().im, 0.0);
            let distance = 0.003;
            let propagation = jet.propagation(distance).unwrap();
            near(
                propagation.group_delay_s,
                distance * group / VACUUM_SPEED_M_S,
                5e-26,
            );
            near(propagation.transfer.norm_sq(), 1.0, 5e-16);
            assert!(propagation.group_delay_s > propagation.phase_delay_s);
        }
    }
}

#[test]
fn multifrequency_films_keep_source_state_and_refuse_unprovided_causal_storage() {
    let ambient = unit_index_port();
    let source = silica_state(200.0, MU);
    let changed = silica_state(201.0, MU);
    let silica = SilicaOpticalResponse::try_lossless(&source, "T").unwrap();
    let make = |source, wavelength| {
        PlanarDielectricLayers::try_new(
            &ambient,
            DielectricResponse::Lossless(&ambient),
            FREQUENCY_WAVELENGTH_NM / wavelength,
            &[DielectricLayer {
                response: DielectricResponse::SilicaOptical(source),
                thickness_m: 1e-6,
            }],
        )
        .unwrap()
    };
    let stacks = [make(silica, 532.0), make(silica, 1550.0)];
    assert!(stacks[0].same_sources_and_geometry(&stacks[1]));
    assert!(!stacks[0].same_sources_and_geometry(&make(
        SilicaOpticalResponse::try_lossless(&changed, "T").unwrap(),
        1550.0,
    )));
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
            error.to_string().contains("original causal material reservoir"),
            "{error}"
        ),
        Ok(_) => panic!("finite-band CHARMS data do not supply a causal time-domain reservoir"),
    }
}
