use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue, Provenance,
    UncertaintyModel,
};
use fs_qty::Temperature;
use fs_qty::semantic::{AttenuationConvention, QuantityKind, SemanticType, ValueForm};

const WAVELENGTHS_UM: [f64; 9] = [8.0, 8.25, 9.0, 9.5, 10.0, 10.25, 11.0, 12.0, 14.0];

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new().with("T", temperature).unwrap()
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
            source: format!("synthetic infrared source: {name}"),
            license: "internal-test-use".into(),
            artifact: None,
        },
    }
}

fn assemble(claims: ClaimSet) -> MaterialCard {
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic finite-band optical fixture".into(),
            phase: "solid".into(),
            process: "original wavelength table".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn card(dielectric: bool, edit: impl Fn(PropertyClaim) -> Option<PropertyClaim>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (index, wavelength_um) in WAVELENGTHS_UM.into_iter().enumerate() {
        let wavelength = wavelength_um * 1.0e-6;
        let (index_name, n, delta_n, loss_name, loss, loss_quantity) = if dielectric {
            (
                SPECTRAL_DIELECTRIC_IOR_PROPERTIES[index],
                1.75 - wavelength_um / 64.0,
                0.0625,
                SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[index],
                if index == 4 { 6.25 } else { 0.25 },
                SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
            )
        } else {
            (
                SPECTRAL_CONDUCTOR_ETA_PROPERTIES[index],
                1.25 + (wavelength_um - 8.0) / 8.0,
                0.25,
                SPECTRAL_CONDUCTOR_K_PROPERTIES[index],
                2.0 + (wavelength_um - 8.0) / 4.0,
                SPECTRAL_CONDUCTOR_K_QUANTITY,
            )
        };
        for source in [
            claim(
                SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[index],
                WAVELENGTH_QUANTITY,
                wavelength,
                wavelength,
            ),
            claim(index_name, INDEX_QUANTITY, n, n + delta_n),
            claim(loss_name, loss_quantity, loss, loss + 0.5),
        ] {
            if let Some(source) = edit(source) {
                claims.insert_claim(source).unwrap();
            }
        }
    }
    assemble(claims)
}

fn constant(source: &mut PropertyClaim, value: f64) {
    source.value = PropertyValue::Scalar {
        value,
        dims: source.key.dims(),
    };
    source.interpolation = InterpolationPolicy::ConstantWithinValidity;
}

#[test]
fn g0_infrared_original_grid_coefficients_and_receipts_survive_state_resolution() {
    let conductor_card = card(false, Some);
    let conductor = resolve_spectral_conductor_state_point(
        &conductor_card,
        &point(450.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(conductor.resolved().properties().len(), 27);
    assert_eq!(conductor.wavelength_bounds_m(), (8.0e-6, 14.0 * 1.0e-6));
    for (sample, wavelength_um) in conductor.samples().iter().zip(WAVELENGTHS_UM) {
        assert_eq!(sample.wavelength_nm, (wavelength_um * 1.0e-6) * 1.0e9);
        assert_eq!(sample.eta, 1.375 + (wavelength_um - 8.0) / 8.0);
        assert_eq!(sample.k, 2.25 + (wavelength_um - 8.0) / 4.0);
    }
    assert_eq!(conductor.identity(), conductor.complex_index().identity());
    assert_eq!(conductor.resolved(), conductor.complex_index().resolved());
    let interpolated = conductor.complex_index().at_wavelength(10_125.0).unwrap();
    assert!((interpolated.eta - (1.375 + 2.125 / 8.0)).abs() < 2e-14);
    assert!((interpolated.k - (2.25 + 2.125 / 4.0)).abs() < 2e-14);
    let dielectric_card = card(true, Some);
    let dielectric = resolve_spectral_tabulated_dielectric_state_point(
        &dielectric_card,
        &point(450.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(dielectric.resolved().properties().len(), 27);
    assert_eq!(
        dielectric.wavelength_bounds_m(),
        conductor.wavelength_bounds_m()
    );
    for (index, sample) in dielectric.samples().iter().enumerate() {
        assert_eq!(sample.phase_index, 1.78125 - WAVELENGTHS_UM[index] / 64.0);
        assert_eq!(sample.absorption_per_m, if index == 4 { 6.5 } else { 0.5 });
    }
    for (card, resolved) in [
        (&conductor_card, conductor.resolved()),
        (&dielectric_card, dielectric.resolved()),
    ] {
        assert_eq!(resolved.card_identity(), card.content_hash());
        assert_eq!(resolved.query_point(), &[("T".to_owned(), 450.0)]);
        for property in resolved.properties() {
            card.claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
}

#[test]
fn g0_si_wavelength_refresh_rebuilds_the_original_complex_index_bindings() {
    let card = card(false, |mut source| {
        if source.key.name() == SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[8] {
            let PropertyValue::Curve { knots, .. } = &mut source.value else {
                unreachable!()
            };
            knots[1].1 = 16.0e-6;
        }
        Some(source)
    });
    let original = resolve_spectral_conductor_state_point(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let hot = original.at_query_point(&card, &point(600.0)).unwrap();
    assert_eq!(hot.wavelength_bounds_m().1, 16.0e-6);
    assert_eq!(
        hot.complex_index().bindings()[8].wavelength_nm,
        16.0e-6 * 1.0e9
    );
    assert!(hot.complex_index().at_wavelength(15_000.0).is_ok());
    assert!(original.complex_index().at_wavelength(15_000.0).is_err());
    for (old, new) in original
        .resolved()
        .properties()
        .iter()
        .zip(hot.resolved().properties())
    {
        assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
    }
}

#[test]
fn g0_refresh_pins_original_claims_even_when_a_competing_source_covers_later_states() {
    let original = card(false, |mut source| {
        if source.key.name() == SPECTRAL_CONDUCTOR_ETA_PROPERTIES[0] {
            source.validity = ValidityDomain::unconstrained().with("T", 300.0, 450.0);
        }
        Some(source)
    });
    let mut claims = original.claims().clone();
    let mut alternate = claim(
        SPECTRAL_CONDUCTOR_ETA_PROPERTIES[0],
        INDEX_QUANTITY,
        3.0,
        4.0,
    );
    alternate.validity = ValidityDomain::unconstrained().with("T", 400.0, 600.0);
    alternate.provenance.source = "synthetic competing source".into();
    claims.insert_claim(alternate).unwrap();
    let card = assemble(claims);
    let initial = resolve_spectral_conductor_state_point(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let before = initial.clone();
    assert!(
        resolve_spectral_conductor_state_point(
            &card,
            &point(425.0),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let continued = initial.at_query_point(&card, &point(425.0)).unwrap();
    for (old, new) in initial
        .resolved()
        .properties()
        .iter()
        .zip(continued.resolved().properties())
    {
        assert_eq!(old.requirement(), new.requirement());
        assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
    }
    assert!(
        resolve_spectral_conductor_state_point(
            &card,
            &point(550.0),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_ok()
    );
    assert!(continued.at_query_point(&card, &point(550.0)).is_err());
    assert_eq!(initial, before);
}

#[test]
fn g0_source_wavelengths_and_optical_family_are_part_of_identity() {
    let conductor_card = card(false, Some);
    let dielectric_card = card(true, Some);
    let mut claims = conductor_card.claims().clone();
    for (_, source) in dielectric_card.claims().claims_ordered() {
        if !SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES.contains(&source.key.name()) {
            claims.insert_claim(source.clone()).unwrap();
        }
    }
    let card = assemble(claims);
    let mut requirements = spectral_requirements(false).unwrap();
    requirements.extend(
        spectral_requirements(true)
            .unwrap()
            .into_iter()
            .filter(|item| !SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES.contains(&item.name())),
    );
    let resolved = resolve_material_state_point(
        &card,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let conductor = SpectralConductorStatePoint::try_from_material(resolved.clone()).unwrap();
    let dielectric =
        SpectralTabulatedDielectricStatePoint::try_from_material(resolved.clone()).unwrap();
    assert_eq!(conductor.resolved(), &resolved);
    assert_eq!(dielectric.resolved(), &resolved);
    assert_eq!(resolved.properties().len(), 45);
    assert_ne!(conductor.identity(), dielectric.identity());
    assert_ne!(conductor.identity(), resolved.identity());
    let mut changed = ClaimSet::new();
    for (_, source) in conductor_card.claims().claims_ordered() {
        let mut source = source.clone();
        if source.key.name() == SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[0] {
            constant(&mut source, 7.5e-6);
        }
        changed.insert_claim(source).unwrap();
    }
    let changed = resolve_spectral_conductor_state_point(
        &assemble(changed),
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original = resolve_spectral_conductor_state_point(
        &conductor_card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_ne!(original.identity(), changed.identity());
    assert_ne!(
        original.wavelength_bounds_m(),
        changed.wavelength_bounds_m()
    );
}

#[test]
fn g0_missing_noncanonical_and_dynamic_sources_cannot_be_frozen_into_optical_knots() {
    for quantity in [
        QuantitySpec::dimensional(Dims([-1, 0, 0, 0, 0, 0])),
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Attenuation(AttenuationConvention::AmplitudeNepersPerMetre),
            ValueForm::Static,
        )),
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::Attenuation(AttenuationConvention::DecibelsPerMetre),
            ValueForm::Static,
        )),
    ] {
        let card = card(true, |mut source| {
            if source.key.name() == SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[0] {
                source.key = PropertyKey::with_quantity(source.key.name(), quantity);
            }
            Some(source)
        });
        assert!(
            resolve_spectral_tabulated_dielectric_state_point(
                &card,
                &point(300.0),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    for missing in [true, false] {
        let card = card(false, |mut source| {
            if source.key.name() == SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[0] {
                if missing {
                    return None;
                }
                source.key = PropertyKey::new(source.key.name(), Dims::NONE);
                constant(&mut source, 8.0e-6);
            }
            Some(source)
        });
        assert!(
            resolve_spectral_conductor_state_point(
                &card,
                &point(300.0),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    let untyped_k = card(false, |mut source| {
        if source.key.name() == SPECTRAL_CONDUCTOR_K_PROPERTIES[0] {
            source.key = PropertyKey::new(source.key.name(), Dims::NONE);
        }
        Some(source)
    });
    assert!(
        resolve_spectral_conductor_state_point(
            &untyped_k,
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let dynamic = card(false, |mut source| {
        source.validity = source.validity.with("omega", 1.0, 2.0);
        Some(source)
    });
    assert!(
        resolve_spectral_conductor_state_point(
            &dynamic,
            &point(300.0).with("omega", 1.5).unwrap(),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}

#[test]
fn g0_wavelength_order_passivity_and_native_dielectric_domain_are_readmitted() {
    for (property, value) in [
        (SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[0], 0.0),
        (
            SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[4],
            WAVELENGTHS_UM[3] * 1.0e-6,
        ),
        (SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[4], 9.0e-6),
        (SPECTRAL_DIELECTRIC_IOR_PROPERTIES[3], 0.99),
        (SPECTRAL_DIELECTRIC_IOR_PROPERTIES[3], f64::MAX),
        (SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[3], -0.1),
    ] {
        let card = card(true, |mut source| {
            if source.key.name() == property {
                constant(&mut source, value);
            }
            Some(source)
        });
        assert!(
            resolve_spectral_tabulated_dielectric_state_point(
                &card,
                &point(300.0),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err(),
            "{property} = {value}"
        );
    }
    let card = card(false, |mut source| {
        if source.key.name() == SPECTRAL_CONDUCTOR_K_PROPERTIES[2] {
            constant(&mut source, -1.0);
        }
        Some(source)
    });
    assert!(
        resolve_spectral_conductor_state_point(
            &card,
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}

#[test]
fn g0_explicit_spectral_source_admits_high_index_with_native_square_headroom() {
    for index in [4.0, (f64::MAX / 128.0).sqrt()] {
        let card = card(true, |mut source| {
            if SPECTRAL_DIELECTRIC_IOR_PROPERTIES.contains(&source.key.name()) {
                constant(&mut source, index);
            }
            Some(source)
        });
        let state = resolve_spectral_tabulated_dielectric_state_point(
            &card,
            &point(450.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert!(
            state
                .samples()
                .iter()
                .all(|sample| sample.phase_index == index)
        );
    }
}

#[test]
fn g0_refresh_rechecks_adopted_passivity_without_publishing_partial_state() {
    let card = card(false, |mut source| {
        if source.key.name() == SPECTRAL_CONDUCTOR_K_PROPERTIES[2] {
            let PropertyValue::Curve { knots, .. } = &mut source.value else {
                unreachable!()
            };
            knots[1].1 = -1.0;
        }
        Some(source)
    });
    let requirements: Vec<_> = spectral_requirements(false)
        .unwrap()
        .into_iter()
        .map(|item| {
            ScalarPropertyRequirement::try_with_quantity(
                item.name(),
                item.quantity(),
                ScalarAdmissibility::Finite,
            )
            .unwrap()
        })
        .collect();
    let original = resolve_material_state_point(
        &card,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = SpectralConductorStatePoint::try_from_material(original.clone()).unwrap();
    assert!(source.at_query_point(&card, &point(600.0)).is_err());
    assert_eq!(source.resolved(), &original);
    assert!(source.at_query_point(&card, &point(601.0)).is_err());
    assert!(source.at_query_point(&card, &QueryPoint::new()).is_err());
}
