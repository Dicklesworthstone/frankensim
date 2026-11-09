//! Authored alpha spectra test composition and Maxwell transport, not measured
//! water/silica loss or material qualification.
use super::*;
use crate::plane_wave::{
    DielectricResponse, Polarization, SilicaOpticalResponse, WaterOpticalResponse,
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS, ELECTRIC_PERMITTIVITY_PROPERTY,
        MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY,
    },
    state_point::{
        DENSITY_PROPERTY, MaterialPropertySelection, ResolvedMaterialStatePoint,
        SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
        SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY, SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES,
        ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_qty::{Density, Dims, QuantitySpec};

const MU: f64 = 1.256_637_062_12e-6;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:.17e} != {expected:.17e}");
}
fn complex_near(actual: C64, expected: C64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual:?} != {expected:?}");
}
fn claim(name: &str, quantity: QuantitySpec, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::with_quantity(name, quantity),
        value: PropertyValue::Scalar { value, dims: quantity.dims() },
        validity: ValidityDomain::unconstrained().with("T", 270.0, 310.0),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "authored analytic-phase/current-alpha control".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn resolve(properties: Vec<PropertyClaim>, temperature: f64) -> ResolvedMaterialStatePoint {
    let mut requirements = Vec::new();
    let mut claims = ClaimSet::new();
    for property in properties {
        requirements.push(
            ScalarPropertyRequirement::try_with_key(&property.key, ScalarAdmissibility::Finite)
                .unwrap(),
        );
        claims.insert_claim(property).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "caller selects original water or silica applicability".into(),
            phase: "homogeneous authored control".into(),
            process: "complex propagation and absorbed-work control".into(),
            revision: 0,
        },
        claims,
        vec![],
    ).unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new().with("T", temperature).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}
fn absorption(band_nm: [f64; 2], alpha: [f64; 9], temperature: f64)
    -> SpectralDielectricAbsorptionStatePoint
{
    let mut properties = vec![
        claim(MAGNETIC_PERMEABILITY_PROPERTY, QuantitySpec::dimensional(MAGNETIC_PERMEABILITY_DIMS), MU),
        claim(DENSITY_PROPERTY, QuantitySpec::dimensional(Density::DIMS), 997.0),
    ];
    for i in 0..9 {
        let nm = band_nm[0] + (band_nm[1] - band_nm[0]) * i as f64 / 8.0;
        properties.push(claim(
            SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
            nm * 1e-9,
        ));
        properties.push(claim(
            SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
            SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
            alpha[i],
        ));
    }
    SpectralDielectricAbsorptionStatePoint::try_from_material(resolve(properties, temperature)).unwrap()
}
fn port() -> ResolvedMaterialStatePoint {
    resolve(vec![
        claim(MAGNETIC_PERMEABILITY_PROPERTY, QuantitySpec::dimensional(MAGNETIC_PERMEABILITY_DIMS), MU),
        claim(ELECTRIC_PERMITTIVITY_PROPERTY, QuantitySpec::dimensional(ELECTRIC_PERMITTIVITY_DIMS),
            1.0 / (MU * VACUUM_SPEED_M_S.powi(2))),
    ], 295.0)
}
fn responses(source: &SpectralDielectricAbsorptionStatePoint) -> [DielectricResponse<'_>; 2] {
    [
        DielectricResponse::WaterOptical(
            WaterOpticalResponse::try_with_current_absorption(source, "T", 60e-9).unwrap(),
        ),
        DielectricResponse::SilicaOptical(
            SilicaOpticalResponse::try_with_current_absorption(source, "T", 60e-9).unwrap(),
        ),
    ]
}
fn lossless(source: &SpectralDielectricAbsorptionStatePoint) -> [DielectricResponse<'_>; 2] {
    [
        DielectricResponse::WaterOptical(WaterOpticalResponse::try_lossless(source.resolved(), "T").unwrap()),
        DielectricResponse::SilicaOptical(SilicaOpticalResponse::try_lossless(source.resolved(), "T").unwrap()),
    ]
}

#[test]
fn original_current_loss_reaches_both_coherent_films_and_local_work_once() {
    let source = absorption([600.0, 1000.0], [200_000.0; 9], 295.0);
    let ambient = port();
    let wavelength_nm = 775.0;
    let omega = FREQUENCY_WAVELENGTH_NM / wavelength_nm;
    let k0 = omega / VACUUM_SPEED_M_S;
    let cosine = 0.8_f64;
    let thickness = 2e-6;
    for (response, lossless) in responses(&source).into_iter().zip(lossless(&source)) {
        assert!(core::ptr::eq(response.material(), source.resolved()));
        let phase = lossless.bulk_dispersion(omega).unwrap().wavenumber_rad_m().re / k0;
        let n = C64::new(phase, -200_000.0 * wavelength_nm * 1e-9 / (4.0 * core::f64::consts::PI));
        let q = (n * n - C64::from_re(1.0 - cosine * cosine)).sqrt();
        let angle = -k0 * q.re * thickness;
        let p = C64::new(angle.cos(), angle.sin()).scale((k0 * q.im * thickness).exp());
        let stack = PlanarDielectricLayers::try_new(
            &ambient, DielectricResponse::Lossless(&ambient), omega,
            &[DielectricLayer { response, thickness_m: thickness }],
        ).unwrap();
        for polarization in [Polarization::Te, Polarization::Tm] {
            let yi = match polarization {
                Polarization::Te => cosine, Polarization::Tm => cosine.recip(),
            } / (MU * VACUUM_SPEED_M_S);
            let yt = match polarization {
                Polarization::Te => q, Polarization::Tm => n * n / q,
            }.scale(1.0 / (MU * VACUUM_SPEED_M_S));
            let r = (C64::from_re(yi) - yt) / (C64::from_re(yi) + yt);
            let denominator = C64::ONE - r * r * p * p;
            let reflected = r * (C64::ONE - p * p) / denominator;
            let transmitted = p * (C64::ONE - r * r) / denominator;
            let expected_heat_fraction = 1.0 - reflected.norm_sq() - transmitted.norm_sq();
            assert!(expected_heat_fraction > 0.2);
            let observed = stack.response(cosine, polarization).unwrap();
            complex_near(observed.reflected_tangential_e, reflected, 3e-14);
            complex_near(observed.transmitted_tangential_e, transmitted, 3e-14);
            near(observed.absorbed_power_fraction, expected_heat_fraction, 4e-14);
            let field = stack.prepare_layer_field(0, cosine, polarization).unwrap();
            near(field.absorbed_power_w_m2(), 0.5 * yi * expected_heat_fraction, 8e-17);
        }
        // Angle-based incidence still requires an exactly lossless exterior.
        let absorbing_incident = crate::plane_wave::PlanarDielectricInterface::try_new_with_incident_response(
            response, DielectricResponse::Lossless(&ambient),
        ).unwrap();
        assert!(absorbing_incident.response(omega, cosine, Polarization::Te).is_err());
    }
}

#[test]
fn original_analytic_phase_and_linear_alpha_keep_exact_complex_group_and_gvd() {
    let wavelength_m = 775e-9_f64;
    let omega = core::f64::consts::TAU * VACUUM_SPEED_M_S / wavelength_m;
    for slope in [0.0, 2e11] {
        let alpha = core::array::from_fn(|i| 100_000.0 + slope * (50e-9 * i as f64));
        let source = absorption([600.0, 1000.0], alpha, 295.0);
        for (response, lossless) in responses(&source).into_iter().zip(lossless(&source)) {
            let jet = response.bulk_dispersion(omega).unwrap();
            let reference = lossless.bulk_dispersion(omega).unwrap();
            let expected_alpha = 100_000.0 + slope * (wavelength_m - 600e-9);
            near(-2.0 * jet.wavenumber_rad_m().im, expected_alpha, 2e-9);
            near(jet.group_slowness_s_m(), reference.group_slowness_s_m(), 2e-23);
            near(jet.group_velocity_dispersion_s2_m(), reference.group_velocity_dispersion_s2_m(), 2e-38);
            let group_loss = slope * wavelength_m.powi(2) / (4.0 * core::f64::consts::PI);
            near(jet.wavenumber_derivative_s_m().im * VACUUM_SPEED_M_S, group_loss, 3e-16);
            near(jet.wavenumber_second_derivative_s2_m().im * omega * VACUUM_SPEED_M_S,
                -2.0 * group_loss, 6e-16);
            if slope == 0.0 {
                assert_eq!(jet.wavenumber_derivative_s_m().im, 0.0);
                assert_eq!(jet.wavenumber_second_derivative_s2_m().im, 0.0);
            }
            let travel = jet.propagation(3e-6).unwrap();
            near(travel.transfer.norm_sq(), (-expected_alpha * 3e-6).exp(), 3e-15);
        }
    }
}

#[test]
fn zero_loss_keeps_its_source_band_and_absorption_cusps_are_not_smoothed() {
    let transparent = absorption([600.0, 1000.0], [0.0; 9], 295.0);
    let samples = transparent.absorption_samples();
    let constant = core::f64::consts::TAU * VACUUM_SPEED_M_S;
    for response in responses(&transparent) {
        for index in [0, 8] {
            let edge = constant / samples[index].wavelength_m;
            assert_eq!(response.permittivity(edge).unwrap().im, 0.0);
            assert!(response.bulk_dispersion(edge).is_err());
        }
        assert!(response.permittivity(FREQUENCY_WAVELENGTH_NM / 550.0).is_err());
        assert!(response.permittivity(FREQUENCY_WAVELENGTH_NM / 1050.0).is_err());
    }
    let water = WaterOpticalResponse::try_with_current_absorption(&transparent, "T", 60e-9).unwrap();
    let silica = SilicaOpticalResponse::try_with_current_absorption(&transparent, "T", 60e-9).unwrap();
    assert!(core::ptr::eq(water.current_absorption().unwrap(), &transparent));
    assert!(core::ptr::eq(silica.current_absorption().unwrap(), &transparent));
    assert!(!water.same_source(WaterOpticalResponse::try_lossless(transparent.resolved(), "T").unwrap()));
    assert!(!silica.same_source(SilicaOpticalResponse::try_lossless(transparent.resolved(), "T").unwrap()));
    assert!(!water.same_source(WaterOpticalResponse::try_with_current_absorption(&transparent, "T", 70e-9).unwrap()));
    assert!(!silica.same_source(SilicaOpticalResponse::try_with_current_absorption(&transparent, "T", 70e-9).unwrap()));
    for gap in [0.0, -1.0, 40e-9, f64::INFINITY, f64::NAN] {
        assert!(WaterOpticalResponse::try_with_current_absorption(&transparent, "T", gap).is_err());
        assert!(SilicaOpticalResponse::try_with_current_absorption(&transparent, "T", gap).is_err());
    }
    let cusp = absorption([600.0, 1000.0], core::array::from_fn(|i| if i == 4 { 200_000.0 } else { 100_000.0 }), 295.0);
    let omega = constant / cusp.absorption_samples()[4].wavelength_m;
    for response in responses(&cusp) {
        assert!(response.permittivity(omega).is_ok());
        assert!(response.bulk_dispersion(omega).is_err());
        assert!(response.bulk_dispersion(omega * 1.001).is_ok());
    }
    // These independently valid loss spectra do not enlarge either phase law.
    let outside = absorption([3000.0, 3400.0], [0.0; 9], 295.0);
    assert!(WaterOpticalResponse::try_with_current_absorption(&outside, "T", 60e-9).is_err());
    assert!(SilicaOpticalResponse::try_with_current_absorption(&outside, "T", 60e-9).is_err());
    let hot_silica = absorption([600.0, 1000.0], [100_000.0; 9], 305.0);
    assert!(SilicaOpticalResponse::try_with_current_absorption(&hot_silica, "T", 60e-9).is_err());
}
