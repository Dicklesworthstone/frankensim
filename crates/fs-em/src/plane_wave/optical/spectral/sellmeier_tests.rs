use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_material::field_storage::{MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY};
use fs_material::state_point::{
    MaterialPropertySelection, SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES,
    SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
    SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES,
    SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES, ScalarAdmissibility, ScalarPropertyRequirement,
    VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES, VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES,
    resolve_material_state_point,
};
use fs_qty::{Dims, QuantitySpec, Temperature};

fn point(temperature_k: f64) -> QueryPoint {
    QueryPoint::new().with("T", temperature_k).unwrap()
}
fn claim(name: &str, quantity: QuantitySpec, cold: f64, hot: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::with_quantity(name, quantity),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, cold), (600.0, hot)],
            dims: quantity.dims(),
        },
        validity: ValidityDomain::unconstrained().with("T", 300.0, 600.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(),
        provenance: Provenance {
            source: format!("synthetic analytic IR wave fixture: {name}"),
            license: "internal-test-use".into(),
            artifact: None,
        },
    }
}
fn card(include_mu: bool, hot_pole: bool) -> MaterialCard {
    let mut sources = Vec::new();
    for i in 0..3 {
        sources.push(claim(
            VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[i],
            QuantitySpec::dimensional(Dims::NONE),
            [1.1, 0.2, 0.8][i],
            [1.2, 0.2, 0.8][i],
        ));
        let cold = [0.01, 0.04, 100.0][i] * 1e-12;
        let hot = if i == 0 && hot_pole { 1.21e-12 } else { cold };
        sources.push(claim(
            VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES[i],
            QuantitySpec::dimensional(Dims([2, 0, 0, 0, 0, 0])),
            cold,
            hot,
        ));
    }
    let length = QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0]));
    for (name, value) in SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES
        .into_iter()
        .zip([1e-6, 2e-6])
    {
        sources.push(claim(name, length, value, value));
    }
    for (i, wavelength) in [1.2, 1.25, 1.35, 1.4, 1.5, 1.6, 1.65, 1.7, 1.8]
        .into_iter()
        .enumerate()
    {
        sources.push(claim(
            SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
            length,
            wavelength * 1e-6,
            wavelength * 1e-6,
        ));
        sources.push(claim(
            SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
            SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
            0.3 + 0.05 * i as f64,
            0.5 + 0.05 * i as f64,
        ));
    }
    if include_mu {
        sources.push(claim(
            MAGNETIC_PERMEABILITY_PROPERTY,
            QuantitySpec::dimensional(MAGNETIC_PERMEABILITY_DIMS),
            1.3e-6,
            1.3e-6,
        ));
    }
    let mut claims = ClaimSet::new();
    for source in sources {
        claims.insert_claim(source).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic analytic infrared medium".into(),
            phase: "solid".into(),
            process: "equation verification, not measured source accuracy".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}
fn source(card: &MaterialCard, temperature: f64) -> SpectralOpticalStatePoint {
    let requirements: Vec<_> = card
        .claims()
        .claims_ordered()
        .map(|(_, claim)| {
            ScalarPropertyRequirement::try_with_key(&claim.key, ScalarAdmissibility::Finite)
                .unwrap()
        })
        .collect();
    let resolved = resolve_material_state_point(
        card,
        &point(temperature),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    SpectralOpticalStatePoint::DielectricSellmeier(
        SpectralSellmeierStatePoint::try_from_material(resolved).unwrap(),
    )
}
fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(f64::MIN_POSITIVE),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_original_analytic_ir_source_keeps_si_permittivity_and_power_absorption() {
    let card = card(true, false);
    let source = source(&card, 300.0);
    let response = SpectralOpticalResponse::try_new((&source).into(), 2e-7).unwrap();
    let wavelength = 1.55e-6;
    let omega = FREQUENCY_WAVELENGTH_M / wavelength;
    let index = response.complex_index(omega).unwrap();
    let x = 1.55_f64.powi(2);
    let expected =
        (1.0 + 1.1 * x / (x - 0.01) + 0.2 * x / (x - 0.04) + 0.8 * x / (x - 100.0)).sqrt();
    let alpha = 0.525;
    near(index.re, expected, 3e-15);
    near(
        index.im,
        -alpha * wavelength / (4.0 * core::f64::consts::PI),
        3e-15,
    );
    let wavenumber = index.scale(core::f64::consts::TAU / wavelength);
    near(-2.0 * wavenumber.im, alpha, 4e-15);
    near(
        det::exp(2.0 * wavenumber.im * 0.7),
        det::exp(-alpha * 0.7),
        4e-15,
    );
    let expected_epsilon = (index * index).scale(1.0 / (1.3e-6 * VACUUM_SPEED_M_S.powi(2)));
    let actual_epsilon = response.permittivity_f_m(omega).unwrap();
    near(actual_epsilon.re, expected_epsilon.re, 4e-15);
    near(actual_epsilon.im, expected_epsilon.im, 4e-15);
    assert!(core::ptr::eq(response.material(), source.resolved()));
    assert!(matches!(
        response.source(),
        SpectralOpticalSource::DielectricSellmeier(_)
    ));
    for property in response.material().properties() {
        card.claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
}

#[test]
fn g0_common_source_band_frequency_edges_gap_and_same_state_permeability_are_explicit() {
    let card = card(true, false);
    let source = source(&card, 300.0);
    assert!(SpectralOpticalResponse::try_new((&source).into(), 1e-8).is_err());
    let response = SpectralOpticalResponse::try_new((&source).into(), 2e-7).unwrap();
    let SpectralOpticalStatePoint::DielectricSellmeier(value) = &source else {
        unreachable!()
    };
    let (low, high) = value.wavelength_bounds_m();
    let low_frequency = FREQUENCY_WAVELENGTH_M / high;
    let high_frequency = FREQUENCY_WAVELENGTH_M / low;
    assert!(response.complex_index(low_frequency).is_ok());
    assert!(response.complex_index(high_frequency).is_ok());
    for omega in [
        low_frequency.next_down(),
        high_frequency.next_up(),
        FREQUENCY_WAVELENGTH_M / 1.1e-6,
        FREQUENCY_WAVELENGTH_M / 1.9e-6,
        f64::NAN,
    ] {
        assert!(response.complex_index(omega).is_err());
    }
    let missing_mu = self::source(&self::card(false, false), 300.0);
    assert!(SpectralOpticalResponse::try_new((&missing_mu).into(), 2e-7).is_err());
}

#[test]
fn g1_spectral_refresh_preserves_analytic_family_and_original_selected_claims() {
    let card = card(true, false);
    let original = source(&card, 300.0);
    let response = SpectralOpticalResponse::try_new((&original).into(), 2e-7).unwrap();
    let omega = FREQUENCY_WAVELENGTH_M / 1.55e-6;
    let cold = response.complex_index(omega).unwrap();
    let refreshed = response.at_query_point(&card, &point(600.0)).unwrap();
    let hot = refreshed.response().unwrap();
    assert!(matches!(
        hot.source(),
        SpectralOpticalSource::DielectricSellmeier(_)
    ));
    assert_eq!(hot.maximum_knot_span_m(), response.maximum_knot_span_m());
    assert!(hot.complex_index(omega).unwrap().re > cold.re);
    assert!(hot.complex_index(omega).unwrap().im < cold.im);
    for (a, b) in response
        .material()
        .properties()
        .iter()
        .zip(hot.material().properties())
    {
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
    assert_eq!(response.complex_index(omega).unwrap(), cold);
}

#[test]
fn g0_active_pole_in_full_fit_refuses_even_when_outside_absorption_support() {
    let card = card(true, true);
    let original = source(&card, 300.0);
    let response = SpectralOpticalResponse::try_new((&original).into(), 2e-7).unwrap();
    let cold = response
        .complex_index(FREQUENCY_WAVELENGTH_M / 1.55e-6)
        .unwrap();
    // The new pole is1.1um; absorption begins at1.2um, but the original fit
    // starts at1.0um. Checking only the common support would miss this refusal.
    assert!(original.at_query_point(&card, &point(600.0)).is_ok());
    assert!(response.at_query_point(&card, &point(600.0)).is_err());
    assert_eq!(
        response
            .complex_index(FREQUENCY_WAVELENGTH_M / 1.55e-6)
            .unwrap(),
        cold
    );
}
