//! Original measured-table interpolation, thermal lens jets and refusal bounds.

use super::*;
use crate::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
    ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_spectral_tabulated_dielectric_state_point,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, UncertaintyModel,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("thermal", quantity(QuantityKind::AbsoluteTemperature), temperature).unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0).unwrap()
}

fn claim(name: &str, q: QuantitySpec, knots: Vec<(f64, f64)>) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::with_quantity(name, q),
        value: PropertyValue::Curve {
            abscissa: "thermal".into(), abscissa_dims: Temperature::DIMS,
            knots, dims: q.dims(),
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity("thermal", quantity(QuantityKind::AbsoluteTemperature), 280.0, 400.0)
            .with_quantity("pressure", quantity(QuantityKind::Pressure), 90_000.0, 110_000.0),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: Vec::new(),
        provenance: Provenance {
            source: format!("manufactured source thermal refraction: {name}"),
            license: "internal-test-use".into(), artifact: None,
        },
    }
}

fn fixture(edit: impl FnOnce(&mut Vec<PropertyClaim>)) -> (MaterialCard, SpectralTabulatedDielectricStatePoint) {
    let mut claims = Vec::new();
    for index in 0..SPECTRAL_OPTICAL_SAMPLE_COUNT {
        let i = index as f64;
        let wavelength_nm = 400.0 + 100.0 * i;
        claims.push(claim(
            SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[index],
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
            vec![
                (280.0, wavelength_nm * 1e-9),
                (400.0, (wavelength_nm + 120.0 * (0.1 + 0.02 * i)) * 1e-9),
            ],
        ));
        claims.push(claim(
            SPECTRAL_DIELECTRIC_IOR_PROPERTIES[index],
            QuantitySpec::dimensional(Dims::NONE),
            vec![
                (280.0, 1.6 - 0.015 * i),
                (400.0, 1.6 - 0.015 * i + 120.0 * (0.0002 + 0.00003 * i)),
            ],
        ));
        claims.push(claim(
            SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[index],
            SPECTRAL_DIELECTRIC_ABSORPTION_QUANTITY,
            vec![(280.0, 1.0 + 0.8 * i), (400.0, 1.0 + 0.8 * i + 120.0 * (0.01 + 0.001 * i))],
        ));
    }
    // A real auxiliary source cusp must remain admitted and retained; it does
    // not invent an optical cusp when the active optical claims stay affine.
    claims.push(claim(DENSITY_PROPERTY, QuantitySpec::dimensional(Density::DIMS),
        vec![(280.0, 2500.0), (340.0, 2400.0), (400.0, 2250.0)]));
    edit(&mut claims);
    let mut set = ClaimSet::new();
    for claim in claims { set.insert_claim(claim).unwrap(); }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "manufactured tabulated thermal lens".into(), phase: "solid".into(),
        process: "moving source coordinates".into(), revision: 0,
    }, set, Vec::new()).unwrap();
    let source = resolve_spectral_tabulated_dielectric_state_point(
        &card, &point(300.0), MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let resolved = source.resolved().with_additional_requirements(&card, &[
        ScalarPropertyRequirement::try_new(DENSITY_PROPERTY, Density::DIMS,
            ScalarAdmissibility::StrictlyPositive).unwrap(),
    ], MaterialPropertySelection::SingleClaimOnly).unwrap();
    (card, SpectralTabulatedDielectricStatePoint::try_from_material(resolved).unwrap())
}

fn path(card: &MaterialCard, source: SpectralTabulatedDielectricStatePoint)
    -> Result<SpectralDielectricTemperaturePath, AbsorptionTemperaturePathError> {
    SpectralDielectricTemperaturePath::try_new(card, source, "thermal", [280.0, 400.0], 128, || false)
}

fn knots<'a>(claims: &'a mut [PropertyClaim], property: &str) -> &'a mut Vec<(f64, f64)> {
    let claim = claims.iter_mut().find(|claim| claim.key.name() == property).unwrap();
    let PropertyValue::Curve { knots, .. } = &mut claim.value else { unreachable!() };
    knots
}

fn near(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 3e-11 * expected.abs().max(1e-10),
        "{actual:.17e} != {expected:.17e}");
}

#[test]
fn exact_original_pins_and_moving_grid_give_rational_phase_group_and_thermal_jets() {
    let (card, source) = fixture(|claims| {
        let mut competitor = claim(SPECTRAL_DIELECTRIC_IOR_PROPERTIES[2],
            QuantitySpec::dimensional(Dims::NONE), vec![(280.0, 3.0), (400.0, 3.0)]);
        competitor.validity = competitor.validity.with_quantity("thermal",
            quantity(QuantityKind::AbsoluteTemperature), 350.0, 400.0);
        claims.push(competitor);
    });
    let prepared = path(&card, source.clone()).unwrap();
    assert_eq!(prepared.source(), &source);
    assert_eq!(prepared.source().resolved().properties().len(), 28);
    assert_eq!(prepared.source_knot_count(), 57);
    near(prepared.maximum_wavelength_knot_span_m(), 102.4e-9);
    assert_eq!(prepared.temperature_axis(), "thermal");
    assert_eq!(prepared.temperature_range_k(), [280.0, 400.0]);
    assert_eq!(prepared.source().resolved().state_coordinate("thermal"), Some(300.0));
    assert!(resolve_spectral_tabulated_dielectric_state_point(
        &card, &point(360.0), MaterialPropertySelection::SingleClaimOnly,
    ).is_err());
    let cell = prepared.prepare(650e-9).unwrap();
    assert_eq!(cell.wavelength_m(), 650e-9);
    assert_eq!(cell.temperature_range_k(), [280.0, 400.0]);

    // Independent polynomial quotient oracle in nanometres: N(T)/W(T).
    // This differentiates the expanded numerator, rather than copying the
    // production barycentric moving-grid derivative formula.
    let (l0, l1, w0, w1) = (600.0, 0.14, 100.0, 0.02);
    let (a0, a1, d0, d1) = (1.57, 0.00026, -0.015, 0.00003);
    let p0 = a0 * w0 + (650.0 - l0) * d0;
    let p1 = a0 * w1 + a1 * w0 + (650.0 - l0) * d1 - l1 * d0;
    let p2 = a1 * w1 - l1 * d1;
    for temperature in [280.0, 293.0, 300.0, 340.0, 377.0, 400.0] {
        let u = temperature - 280.0;
        let width = w0 + w1 * u;
        let p = p0 + p1 * u + p2 * u * u;
        let p_prime = p1 + 2.0 * p2 * u;
        let expected_n = p / width;
        let expected_first = (p_prime * width - p * w1) / (width * width);
        let expected_second = (2.0 * p2 * width * width - 2.0 * p_prime * width * w1
            + 2.0 * p * w1 * w1) / (width * width * width);
        let sample = cell.sample(temperature).unwrap();
        near(sample.phase_index, expected_n);
        near(sample.derivative_temperature_per_k, expected_first);
        near(sample.second_derivative_temperature_per_k2, expected_second);
        near(sample.group_index.unwrap(), expected_n - 650.0 * (d0 + d1 * u) / width);
        assert_ne!(sample.second_derivative_temperature_per_k2, 0.0);
        let current = source.at_query_point(&card, &point(temperature)).unwrap();
        let left = current.samples()[2];
        let right = current.samples()[3];
        let q = (650e-9 - left.wavelength_m) / (right.wavelength_m - left.wavelength_m);
        near(sample.phase_index, left.phase_index * (1.0 - q) + right.phase_index * q);
        near(sample.absorption_per_m, current.absorption_per_m_at(650e-9).unwrap());
        near(prepared.phase_index_at(temperature, 650e-9).unwrap(), sample.phase_index);
        near(prepared.absorption_per_m_at(temperature, 650e-9).unwrap(), sample.absorption_per_m);
    }
    assert_eq!(cell.sample(279.0), Err(SpectralDielectricTemperatureError::TemperatureOutsideRange));
    assert_eq!(cell.sample(f64::NAN), Err(SpectralDielectricTemperatureError::TemperatureOutsideRange));
    assert_eq!(prepared.prepare(400e-9).unwrap_err(), SpectralDielectricTemperatureError::WavelengthOutsideRange);
}

#[test]
fn genuine_source_temperature_cusps_and_moving_spectral_crossings_refuse_smooth_cells() {
    for (property, offset) in [
        (SPECTRAL_DIELECTRIC_IOR_PROPERTIES[2], 0.01),
        (SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[2], 0.5e-9),
        (SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[2], 0.3),
    ] {
        let (card, source) = fixture(|claims| {
            let values = knots(claims, property);
            values.insert(1, (340.0, (values[0].1 + values[1].1) * 0.5 + offset));
        });
        let prepared = path(&card, source.clone()).unwrap();
        assert!(prepared.phase_index_at(340.0, 650e-9).is_ok());
        assert!(prepared.absorption_per_m_at(340.0, 650e-9).is_ok());
        assert_eq!(prepared.prepare(650e-9).unwrap_err(), SpectralDielectricTemperatureError::UndefinedDerivative);
        let cut = SpectralDielectricTemperaturePath::try_new(
            &card, source, "thermal", [280.0, 330.0], 128, || false,
        ).unwrap();
        assert!(cut.prepare(650e-9).is_ok());
        assert_ne!(cut.identity(), prepared.identity());
    }
    let (card, source) = fixture(|_| {});
    let prepared = path(&card, source).unwrap();
    assert!(prepared.phase_index_at(280.0, 610e-9).is_ok());
    assert!(prepared.phase_index_at(400.0, 610e-9).is_ok());
    assert_eq!(prepared.prepare(610e-9).unwrap_err(), SpectralDielectricTemperatureError::UndefinedDerivative);
}

#[test]
fn stationary_spectral_knots_keep_phase_jets_and_exact_group_derivative_availability() {
    for spectral_cusp in [false, true] {
        let (card, source) = fixture(|claims| {
            for index in 0..SPECTRAL_OPTICAL_SAMPLE_COUNT {
                let wavelengths = knots(claims, SPECTRAL_OPTICAL_WAVELENGTH_M_PROPERTIES[index]);
                wavelengths[1].1 = wavelengths[0].1;
                let index_values = knots(claims, SPECTRAL_DIELECTRIC_IOR_PROPERTIES[index]);
                let offset = if spectral_cusp && index == 3 { 0.125 } else { 0.0 };
                *index_values = vec![(280.0, 1.5 + offset), (340.0, 1.625 + offset), (400.0, 1.75 + offset)];
                if !spectral_cusp {
                    *knots(claims, SPECTRAL_DIELECTRIC_ABSORPTION_PER_M_PROPERTIES[index]) =
                        vec![(280.0, 0.0), (400.0, 0.0)];
                }
            }
        });
        let wavelength = source.samples()[2].wavelength_m;
        let prepared = path(&card, source).unwrap();
        let cell = prepared.prepare(wavelength).unwrap();
        for temperature in [280.0, 300.0, 340.0, 400.0] {
            let sample = cell.sample(temperature).unwrap();
            near(sample.phase_index, 1.5 + (temperature - 280.0) / 480.0);
            near(sample.derivative_temperature_per_k, 1.0 / 480.0);
            assert_eq!(sample.second_derivative_temperature_per_k2, 0.0);
            if spectral_cusp { assert_eq!(sample.group_index, None); }
            else {
                near(sample.group_index.unwrap(), sample.phase_index);
                assert_eq!(sample.absorption_per_m, 0.0);
            }
        }
    }
}

#[test]
fn original_phase_headroom_auxiliary_domains_source_budget_and_cancellation_are_atomic() {
    for bad_phase in [0.9, 1e154] {
        let (card, source) = fixture(|claims| {
            let values = knots(claims, SPECTRAL_DIELECTRIC_IOR_PROPERTIES[4]);
            *values = vec![(280.0, 1.54), (380.0, 1.55), (390.0, bad_phase), (400.0, 1.56)];
        });
        assert!(source.at_query_point(&card, &point(280.0)).is_ok());
        assert!(source.at_query_point(&card, &point(400.0)).is_ok());
        assert!(matches!(path(&card, source), Err(AbsorptionTemperaturePathError::Source(
            MaterialStatePointError::InvalidDerived { .. }))));
    }
    let (card, source) = fixture(|claims| knots(claims, DENSITY_PROPERTY)[1].1 = -1.0);
    assert!(matches!(path(&card, source), Err(AbsorptionTemperaturePathError::Source(
        MaterialStatePointError::OutsideConsumerDomain { .. }))));
    let (card, source) = fixture(|_| {});
    let original = source.identity();
    assert!(matches!(SpectralDielectricTemperaturePath::try_new(
        &card, source.clone(), "thermal", [280.0, 400.0], 56, || false,
    ), Err(AbsorptionTemperaturePathError::Budget { required: 57, maximum: 56 })));
    assert!(matches!(SpectralDielectricTemperaturePath::try_new(
        &card, source.clone(), "thermal", [280.0, 400.0], 128, || true,
    ), Err(AbsorptionTemperaturePathError::Cancelled)));
    assert!(matches!(SpectralDielectricTemperaturePath::try_new(
        &card, source.clone(), "pressure", [280.0, 400.0], 128, || false,
    ), Err(AbsorptionTemperaturePathError::TemperatureAxis)));
    let point_path = SpectralDielectricTemperaturePath::try_new(
        &card, source.clone(), "thermal", [300.0, 300.0], 28, || false,
    ).unwrap();
    assert!(point_path.phase_index_at(300.0, 650e-9).is_ok());
    assert_eq!(point_path.prepare(650e-9).unwrap_err(), SpectralDielectricTemperatureError::UndefinedDerivative);
    assert_eq!(source.identity(), original);
}
