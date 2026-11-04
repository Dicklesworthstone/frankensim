use super::*;
use crate::state_point::{
    SpectralOpticalStatePoint, VISIBLE_DIELECTRIC_REFERENCE_DISTANCE_M_PROPERTY,
    VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES, VisibleDielectricAbsorption,
    resolve_visible_sellmeier_state_point,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_qty::Temperature;
use fs_qty::semantic::{AttenuationConvention, QuantityKind, SemanticType, ValueForm};

const WAVELENGTHS_UM: [f64; 9] = [1.2, 1.25, 1.35, 1.4, 1.5, 1.6, 1.65, 1.7, 1.8];

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
            source: format!("synthetic analytic infrared fixture: {name}"),
            license: "internal-test-use".into(),
            artifact: None,
        },
    }
}

fn assemble(claims: ClaimSet) -> MaterialCard {
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic Sellmeier source".into(),
            phase: "solid".into(),
            process: "numerical equation fixture, not measured glass data".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn card(edit: impl Fn(PropertyClaim) -> Option<PropertyClaim>) -> MaterialCard {
    let mut sources = Vec::new();
    for i in 0..3 {
        sources.push(claim(
            VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[i],
            INDEX_QUANTITY,
            [1.1, 0.2, 0.8][i],
            [1.2, 0.2, 0.8][i],
        ));
        let c = [0.01, 0.04, 100.0][i] * 1e-12;
        sources.push(claim(
            VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES[i],
            QuantitySpec::dimensional(Dims([2, 0, 0, 0, 0, 0])),
            c,
            c,
        ));
    }
    for (name, value) in SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES
        .into_iter()
        .zip([1e-6, 2e-6])
    {
        sources.push(claim(name, WAVELENGTH_QUANTITY, value, value));
    }
    for (i, wavelength_um) in WAVELENGTHS_UM.into_iter().enumerate() {
        sources.push(claim(
            SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
            WAVELENGTH_QUANTITY,
            wavelength_um * 1e-6,
            wavelength_um * 1e-6,
        ));
        sources.push(claim(
            SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i],
            SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
            0.1 + 0.02 * i as f64,
            0.4 + 0.02 * i as f64,
        ));
    }
    let mut claims = ClaimSet::new();
    for source in sources {
        if let Some(source) = edit(source) {
            claims.insert_claim(source).unwrap();
        }
    }
    assemble(claims)
}

fn resolve(card: &MaterialCard, temperature: f64) -> SpectralSellmeierStatePoint {
    resolve_spectral_sellmeier_state_point(
        card,
        &point(temperature),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn constant(source: &mut PropertyClaim, value: f64) {
    source.value = PropertyValue::Scalar {
        value,
        dims: source.key.dims(),
    };
    source.interpolation = InterpolationPolicy::ConstantWithinValidity;
}

fn visible_card_with_fit_bounds(bounds: &[&str]) -> MaterialCard {
    let source = card(|source| {
        (VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES.contains(&source.key.name())
            || VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES.contains(&source.key.name())
            || bounds.contains(&source.key.name()))
        .then_some(source)
    });
    let mut claims = source.claims().clone();
    for name in VISIBLE_DIELECTRIC_TRANSMITTANCE_PROPERTIES {
        claims
            .insert_claim(claim(name, INDEX_QUANTITY, 0.9, 0.8))
            .unwrap();
    }
    claims
        .insert_claim(claim(
            VISIBLE_DIELECTRIC_REFERENCE_DISTANCE_M_PROPERTY,
            WAVELENGTH_QUANTITY,
            0.01,
            0.01,
        ))
        .unwrap();
    assemble(claims)
}

#[test]
fn g0_visible_resolution_cannot_discard_either_explicit_fit_bound() {
    let bounds = SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES;
    for present in [&bounds[..1], &bounds[1..], &bounds[..]] {
        let card = visible_card_with_fit_bounds(present);
        let pins = card
            .claims()
            .claims_ordered()
            .filter(|(_, claim)| !bounds.contains(&claim.key.name()))
            .map(|(id, claim)| (claim.key.name().to_owned(), id))
            .collect();
        for selection in [
            MaterialPropertySelection::SingleClaimOnly,
            MaterialPropertySelection::PinnedByProperty(pins),
        ] {
            assert!(matches!(
                resolve_visible_sellmeier_state_point(&card, &point(300.0), selection),
                Err(MaterialStatePointError::InvalidDerived {
                    quantity: "explicit_sellmeier_band_requires_spectral_source",
                })
            ));
        }
    }
}

#[test]
fn g0_legacy_visible_sources_without_explicit_fit_limits_keep_their_original_path() {
    let card = visible_card_with_fit_bounds(&[]);
    let source = resolve_visible_sellmeier_state_point(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(source.resolved().properties().len(), 10);
    assert_eq!(source.sellmeier_coefficients_si().0, [1.1, 0.2, 0.8]);
    assert_eq!(
        source.absorption(),
        VisibleDielectricAbsorption::ReferenceRgb {
            transmittance_linear_rgb: [0.9; 3],
            reference_distance_m: 0.01,
        }
    );
    let refreshed = source.at_query_point(&card, &point(600.0)).unwrap();
    assert_eq!(refreshed.sellmeier_coefficients_si().0, [1.2, 0.2, 0.8]);
    assert_eq!(
        refreshed.absorption(),
        VisibleDielectricAbsorption::ReferenceRgb {
            transmittance_linear_rgb: [0.8; 3],
            reference_distance_m: 0.01,
        }
    );
}

#[test]
fn g0_original_coefficients_fit_bounds_absorption_and_receipts_are_retained() {
    let card = card(Some);
    let state = resolve(&card, 450.0);
    assert_eq!(state.resolved().properties().len(), 26);
    let (b, c) = state.sellmeier_coefficients_si();
    assert!((b[0] - 1.15).abs() < 3e-16);
    assert_eq!([b[1], b[2]], [0.2, 0.8]);
    assert_eq!(c, [0.01, 0.04, 100.0].map(|v| v * 1e-12));
    assert_eq!(state.dispersion_wavelength_bounds_m(), (1e-6, 2e-6));
    assert_eq!(state.wavelength_bounds_m(), (1.2e-6, 1.8e-6));
    for (i, sample) in state.absorption_samples().iter().enumerate() {
        assert_eq!(sample.wavelength_m, WAVELENGTHS_UM[i] * 1e-6);
        assert!((sample.absorption_per_m - (0.25 + 0.02 * i as f64)).abs() < 1e-16);
        assert_eq!(
            state
                .absorption_per_m_at(sample.wavelength_m)
                .unwrap()
                .to_bits(),
            sample.absorption_per_m.to_bits()
        );
    }
    assert!((state.absorption_per_m_at(1.3e-6).unwrap() - 0.28).abs() < 1e-15);
    for property in state.resolved().properties() {
        card.claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    assert_eq!(state.resolved().card_identity(), card.content_hash());
    assert_ne!(state.identity(), state.resolved().identity());
    assert_eq!(
        SpectralOpticalStatePoint::DielectricSellmeier(state.clone()).identity(),
        state.identity()
    );
}

#[test]
fn g0_absorption_availability_never_widens_the_original_fit_domain() {
    let original = card(Some);
    let narrower = card(|mut source| {
        if source.key.name() == SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES[0] {
            constant(&mut source, 1.3e-6);
        }
        if source.key.name() == SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES[1] {
            constant(&mut source, 1.7e-6);
        }
        Some(source)
    });
    let a = resolve(&original, 300.0);
    let b = resolve(&narrower, 300.0);
    assert_eq!(a.sellmeier_coefficients_si(), b.sellmeier_coefficients_si());
    assert_eq!(a.absorption_samples(), b.absorption_samples());
    assert_eq!(b.dispersion_wavelength_bounds_m(), (1.3e-6, 1.7e-6));
    assert_eq!(b.wavelength_bounds_m(), (1.3e-6, 1.7e-6));
    assert_ne!(a.identity(), b.identity());
    for wavelength in [1.19e-6, 1.81e-6, f64::NAN] {
        assert!(b.absorption_per_m_at(wavelength).is_err());
    }
}

#[test]
fn g0_refresh_pins_original_coefficients_bounds_and_absorption_without_source_switching() {
    let original = card(|mut source| {
        if source.key.name() == VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[0] {
            source.validity = ValidityDomain::unconstrained().with("T", 300.0, 450.0);
        }
        Some(source)
    });
    let mut claims = original.claims().clone();
    let mut alternate = claim(
        VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[0],
        INDEX_QUANTITY,
        2.0,
        2.1,
    );
    alternate.validity = ValidityDomain::unconstrained().with("T", 400.0, 600.0);
    alternate.provenance.source = "independent competing source".into();
    claims.insert_claim(alternate).unwrap();
    let card = assemble(claims);
    let initial = resolve(&card, 300.0);
    let original = SpectralOpticalStatePoint::DielectricSellmeier(initial.clone());
    let refreshed = original.at_query_point(&card, &point(425.0)).unwrap();
    for (a, b) in original
        .resolved()
        .properties()
        .iter()
        .zip(refreshed.resolved().properties())
    {
        assert_eq!(a.answer().receipt.selected, b.answer().receipt.selected);
    }
    assert!(
        resolve_spectral_sellmeier_state_point(
            &card,
            &point(425.0),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    assert!(
        resolve_spectral_sellmeier_state_point(
            &card,
            &point(550.0),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_ok()
    );
    assert!(refreshed.at_query_point(&card, &point(550.0)).is_err());
    assert_eq!(original.resolved(), initial.resolved());
}

#[test]
fn g0_missing_wrong_quantity_dynamic_and_disjoint_source_domains_refuse() {
    for name in [
        VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[2],
        SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES[0],
        SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[5],
        SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[7],
    ] {
        let card = card(|source| (source.key.name() != name).then_some(source));
        assert!(
            resolve_spectral_sellmeier_state_point(
                &card,
                &point(300.0),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    for (name, value) in [
        (VISIBLE_DIELECTRIC_SELLMEIER_B_PROPERTIES[0], -0.1),
        (VISIBLE_DIELECTRIC_SELLMEIER_C_M2_PROPERTIES[0], -1e-12),
        (
            SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES[0],
            0.0,
        ),
        (
            SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES[0],
            2.1e-6,
        ),
        (
            SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES[1],
            1e-6,
        ),
        (SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[3], -0.1),
    ] {
        let card = card(|mut source| {
            if source.key.name() == name {
                constant(&mut source, value);
            }
            Some(source)
        });
        assert!(
            resolve_spectral_sellmeier_state_point(
                &card,
                &point(300.0),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err(),
            "{name}={value}"
        );
    }
    for quantity in [
        QuantitySpec::dimensional(Dims([-1, 0, 0, 0, 0, 0])),
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Attenuation(AttenuationConvention::AmplitudeNepersPerMetre),
            ValueForm::Static,
        )),
    ] {
        let card = card(|mut source| {
            if source.key.name() == SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[0] {
                source.key = PropertyKey::with_quantity(source.key.name(), quantity);
            }
            Some(source)
        });
        assert!(
            resolve_spectral_sellmeier_state_point(
                &card,
                &point(300.0),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    let dynamic = card(|mut source| {
        source.validity = source.validity.with("omega", 1.0, 2.0);
        Some(source)
    });
    assert!(
        resolve_spectral_sellmeier_state_point(
            &dynamic,
            &point(300.0).with("omega", 1.5).unwrap(),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}

#[test]
fn g0_adopted_complete_state_retains_auxiliary_properties_and_readmits_refreshed_bounds() {
    let original = card(|mut source| {
        if source.key.name() == SPECTRAL_DIELECTRIC_SELLMEIER_WAVELENGTH_BOUNDS_M_PROPERTIES[1] {
            source = claim(source.key.name(), WAVELENGTH_QUANTITY, 2e-6, 1e-6);
        }
        Some(source)
    });
    let mut claims = original.claims().clone();
    claims
        .insert_claim(claim(
            "auxiliary_density",
            QuantitySpec::dimensional(Dims([-3, 1, 0, 0, 0, 0])),
            2_000.0,
            1_900.0,
        ))
        .unwrap();
    let card = assemble(claims);
    let requirements: Vec<_> = card
        .claims()
        .claims_ordered()
        .map(|(_, claim)| {
            ScalarPropertyRequirement::try_with_key(&claim.key, ScalarAdmissibility::Finite)
                .unwrap()
        })
        .collect();
    let accepted = resolve_material_state_point(
        &card,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = SpectralSellmeierStatePoint::try_from_material(accepted.clone()).unwrap();
    assert_eq!(source.resolved(), &accepted);
    assert_eq!(source.resolved().properties().len(), 27);
    let hot = source.at_query_point(&card, &point(450.0)).unwrap();
    assert_eq!(hot.dispersion_wavelength_bounds_m().1, 1.5e-6);
    assert!(source.at_query_point(&card, &point(600.0)).is_err());
    assert_eq!(source.resolved(), &accepted);
}
