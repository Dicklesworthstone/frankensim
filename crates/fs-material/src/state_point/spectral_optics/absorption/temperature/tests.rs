//! G1/G3/G4: original moving spectra on a prepared material temperature line.

use super::*;
use crate::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
    ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_spectral_dielectric_absorption_state_point,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::{Density, Dims, Pressure, Temperature};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "thermal",
            quantity(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
        .unwrap()
}

fn claim(name: &str, q: QuantitySpec, knots: Vec<(f64, f64)>) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::with_quantity(name, q),
        value: PropertyValue::Curve {
            abscissa: "thermal".into(),
            abscissa_dims: Temperature::DIMS,
            knots,
            dims: q.dims(),
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "thermal",
                quantity(QuantityKind::AbsoluteTemperature),
                280.0,
                400.0,
            )
            .with_quantity("pressure", quantity(QuantityKind::Pressure), 90_000.0, 110_000.0),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: Vec::new(),
        provenance: Provenance {
            source: format!("manufactured original thermal absorption: {name}"),
            license: "internal-test-use".into(),
            artifact: None,
        },
    }
}

fn fixture(
    edit: impl FnOnce(&mut Vec<PropertyClaim>),
) -> (MaterialCard, SpectralDielectricAbsorptionStatePoint) {
    let mut claims = Vec::new();
    for i in 0..SPECTRAL_OPTICAL_SAMPLE_COUNT {
        let w = 400.0 + 50.0 * i as f64;
        claims.push(claim(
            SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[i],
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
            [(280.0, w), (330.0, w + 15.0), (400.0, w + 40.0)]
                .map(|(t, wavelength)| (t, wavelength * 1e-9))
                .to_vec(),
        ));
        let i = i as f64;
        claims.push(claim(
            SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[i as usize],
            SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
            vec![
                (280.0, 1.0 + i),
                (320.0, 3.0 + 2.0 * i),
                (360.0, 2.0 + 1.5 * i),
                (400.0, 5.0 + i),
            ],
        ));
    }
    claims.push(claim(
        DENSITY_PROPERTY,
        QuantitySpec::dimensional(Density::DIMS),
        vec![(280.0, 2500.0), (340.0, 2400.0), (400.0, 2300.0)],
    ));
    edit(&mut claims);
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured moving-band absorber".into(),
            phase: "solid".into(),
            process: "original claim path regression".into(),
            revision: 0,
        },
        set,
        Vec::new(),
    )
    .unwrap();
    let source = resolve_spectral_dielectric_absorption_state_point(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let resolved = source
        .resolved()
        .with_additional_requirements(
            &card,
            &[ScalarPropertyRequirement::try_new(
                DENSITY_PROPERTY,
                Density::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
    (
        card,
        SpectralDielectricAbsorptionStatePoint::try_from_material(resolved).unwrap(),
    )
}

fn prepare(
    card: &MaterialCard,
    source: SpectralDielectricAbsorptionStatePoint,
) -> Result<SpectralDielectricAbsorptionTemperaturePath, AbsorptionTemperaturePathError> {
    SpectralDielectricAbsorptionTemperaturePath::try_new(
        card, source, "thermal", [280.0, 400.0], 128, || false,
    )
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-13 * expected.abs().max(1e-12),
        "{actual:.17e} != {expected:.17e}",
    );
}

#[test]
fn original_pins_and_interlaced_knots_match_current_source_without_hot_queries() {
    let (card, source) = fixture(|claims| {
        let mut competitor = claim(
            SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[0],
            SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
            vec![(280.0, 1000.0), (400.0, 1000.0)],
        );
        competitor.validity = competitor.validity.with_quantity(
            "thermal",
            quantity(QuantityKind::AbsoluteTemperature),
            350.0,
            400.0,
        );
        claims.push(competitor);
    });
    let path = prepare(&card, source.clone()).unwrap();
    assert_eq!(path.source(), &source);
    assert_eq!(path.source().resolved().properties().len(), 19);
    assert_eq!(path.source_knot_count(), 66);
    assert_eq!(path.temperature_axis(), "thermal");
    assert_eq!(path.temperature_range_k(), [280.0, 400.0]);
    assert!(Arc::ptr_eq(&path.prepared, &path.clone().prepared));
    let (lo, hi) = path.common_wavelength_bounds_m().unwrap();
    near(lo, 440e-9);
    near(hi, 800e-9);
    assert!(resolve_spectral_dielectric_absorption_state_point(
        &card, &point(360.0), MaterialPropertySelection::SingleClaimOnly,
    ).is_err());

    for t in [280.0, 293.0, 300.0, 320.0, 330.0, 340.0, 355.0, 360.0, 387.0, 400.0] {
        let current = source.at_query_point(&card, &point(t)).unwrap();
        let knots = current.absorption_samples();
        let (lo, hi) = path.wavelength_bounds_m_at(t).unwrap();
        near(lo, knots[0].wavelength_m);
        near(hi, knots[8].wavelength_m);
        for i in 0..8 {
            let wavelength = 0.3 * knots[i].wavelength_m + 0.7 * knots[i + 1].wavelength_m;
            let expected = 0.3 * knots[i].absorption_per_m + 0.7 * knots[i + 1].absorption_per_m;
            near(path.absorption_per_m_at(t, wavelength).unwrap(), expected);
            near(
                path.absorption_per_m_at(t, wavelength).unwrap(),
                current.absorption_per_m_at(wavelength).unwrap(),
            );
        }
    }
    let future = SpectralDielectricAbsorptionTemperaturePath::try_new(
        &card, source.clone(), "thermal", [350.0, 400.0], 128, || false,
    ).unwrap();
    assert_eq!(future.source().resolved().state_coordinate("thermal"), Some(300.0));
    assert_ne!(future.identity(), path.identity());
    assert_eq!(path.source(), &source);
}

#[test]
fn fixed_pressure_surface_cuts_keep_original_breakpoints_and_refuse_holes() {
    let (card, source) = fixture(|claims| {
        let loss = claims.iter_mut()
            .find(|claim| claim.key.name() == SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[0])
            .unwrap();
        loss.value = PropertyValue::Surface2D {
            axes: ["thermal".into(), "pressure".into()],
            axis_dims: [Temperature::DIMS, Pressure::DIMS],
            samples: vec![
                ([280.0, 90_000.0], 1.0),
                ([400.0, 90_000.0], 5.0),
                ([400.0, 110_000.0], 9.0),
                ([280.0, 110_000.0], 3.0),
            ],
            triangles: vec![[0, 1, 2], [0, 2, 3]],
            dims: SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY.dims(),
        };
    });
    let path = prepare(&card, source.clone()).unwrap();
    for (temperature, expected) in [(280.0, 2.0), (310.0, 3.5), (340.0, 5.0), (370.0, 6.0), (400.0, 7.0)] {
        let wavelength = path.wavelength_bounds_m_at(temperature).unwrap().0;
        near(path.absorption_per_m_at(temperature, wavelength).unwrap(), expected);
        near(
            path.absorption_per_m_at(temperature, wavelength).unwrap(),
            source.at_query_point(&card, &point(temperature)).unwrap()
                .absorption_per_m_at(wavelength).unwrap(),
        );
    }
    let (card, source) = fixture(|claims| {
        let loss = claims.iter_mut()
            .find(|claim| claim.key.name() == SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[0])
            .unwrap();
        loss.value = PropertyValue::Surface2D {
            axes: ["thermal".into(), "pressure".into()],
            axis_dims: [Temperature::DIMS, Pressure::DIMS],
            samples: vec![
                ([280.0, 90_000.0], 1.0), ([320.0, 90_000.0], 2.0),
                ([320.0, 110_000.0], 3.0), ([280.0, 110_000.0], 2.0),
                ([360.0, 90_000.0], 4.0), ([400.0, 90_000.0], 5.0),
                ([400.0, 110_000.0], 6.0), ([360.0, 110_000.0], 5.0),
            ],
            triangles: vec![[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]],
            dims: SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY.dims(),
        };
    });
    assert!(source.at_query_point(&card, &point(280.0)).is_ok());
    assert!(source.at_query_point(&card, &point(400.0)).is_ok());
    assert!(matches!(
        prepare(&card, source),
        Err(AbsorptionTemperaturePathError::Source(MaterialStatePointError::Support { .. }))
    ));
}

#[test]
fn interior_spectrum_auxiliary_domains_budget_and_cancellation_admit_atomically() {
    let (card, source) = fixture(|claims| {
        let wavelength = claims.iter_mut()
            .find(|claim| claim.key.name() == SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[4])
            .unwrap();
        let PropertyValue::Curve { knots, .. } = &mut wavelength.value else { unreachable!() };
        knots[1].1 = 730e-9;
    });
    // Both extreme spectra are ordered. The original interior wavelength knot
    // crosses its neighbour, so checking only thermal endpoints is insufficient.
    assert!(source.at_query_point(&card, &point(280.0)).is_ok());
    assert!(source.at_query_point(&card, &point(400.0)).is_ok());
    assert!(matches!(
        prepare(&card, source),
        Err(AbsorptionTemperaturePathError::Spectrum {
            source: AbsorptionTemperatureSampleError::InvalidSpectrum, ..
        })
    ));

    for property in [DENSITY_PROPERTY, SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[0]] {
        let (card, source) = fixture(|claims| {
            let claim = claims.iter_mut().find(|claim| claim.key.name() == property).unwrap();
            let PropertyValue::Curve { knots, .. } = &mut claim.value else { unreachable!() };
            let index = if property == DENSITY_PROPERTY { 1 } else { 2 };
            knots[index].1 = -1.0;
        });
        assert!(matches!(
            prepare(&card, source),
            Err(AbsorptionTemperaturePathError::Source(
                MaterialStatePointError::OutsideConsumerDomain { .. }
            ))
        ));
    }

    let (card, source) = fixture(|_| {});
    assert!(matches!(
        SpectralDielectricAbsorptionTemperaturePath::try_new(
            &card, source.clone(), "thermal", [280.0, 400.0], 65, || false,
        ),
        Err(AbsorptionTemperaturePathError::Budget { .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        SpectralDielectricAbsorptionTemperaturePath::try_new(
            &card, source.clone(), "thermal", [280.0, 400.0], 66, || {
                polls += 1;
                polls == 23
            },
        ),
        Err(AbsorptionTemperaturePathError::Cancelled)
    ));
    assert!(matches!(
        SpectralDielectricAbsorptionTemperaturePath::try_new(
            &card, source.clone(), "pressure", [280.0, 400.0], 66, || false,
        ),
        Err(AbsorptionTemperaturePathError::TemperatureAxis)
    ));
    for range in [[400.0, 280.0], [0.0, 400.0], [f64::NAN, 400.0]] {
        assert!(matches!(
            SpectralDielectricAbsorptionTemperaturePath::try_new(
                &card, source.clone(), "thermal", range, 66, || false,
            ),
            Err(AbsorptionTemperaturePathError::TemperatureRange)
        ));
    }
    let original = source.clone();
    assert!(prepare(&card, source.clone()).is_ok());
    assert_eq!(source, original);
}

#[test]
fn transparent_moving_bands_keep_local_support_when_no_global_band_exists() {
    let (card, source) = fixture(|claims| {
        for claim in claims {
            if let PropertyValue::Curve { knots, .. } = &mut claim.value {
                if SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES.contains(&claim.key.name()) {
                    for (_, value) in knots {
                        *value = 0.0;
                    }
                } else if SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES.contains(&claim.key.name()) {
                    let base = knots[0].1;
                    knots[1].1 = base + 250e-9;
                    knots[2].1 = base + 500e-9;
                }
            }
        }
    });
    let path = prepare(&card, source.clone()).unwrap();
    assert_eq!(path.common_wavelength_bounds_m(), None);
    for temperature in [280.0, 300.0, 330.0, 375.0, 400.0] {
        let (lo, hi) = path.wavelength_bounds_m_at(temperature).unwrap();
        assert_eq!(path.absorption_per_m_at(temperature, 0.5 * (lo + hi)), Ok(0.0));
        assert_eq!(
            path.absorption_per_m_at(temperature, lo * 0.99),
            Err(AbsorptionTemperatureSampleError::WavelengthOutsideRange),
        );
    }
    assert_eq!(path.absorption_per_m_at(401.0, 1000e-9),
        Err(AbsorptionTemperatureSampleError::TemperatureOutsideRange));
    assert_eq!(path.absorption_per_m_at(300.0, f64::NAN),
        Err(AbsorptionTemperatureSampleError::WavelengthOutsideRange));
    let point_path = SpectralDielectricAbsorptionTemperaturePath::try_new(
        &card, source.clone(), "thermal", [300.0, 300.0], 19, || false,
    ).unwrap();
    assert_eq!(point_path.source_knot_count(), 19);
    assert_eq!(point_path.common_wavelength_bounds_m(), Some(source.wavelength_bounds_m()));
    assert_eq!(point_path.absorption_per_m_at(300.1, 600e-9),
        Err(AbsorptionTemperatureSampleError::TemperatureOutsideRange));
}
