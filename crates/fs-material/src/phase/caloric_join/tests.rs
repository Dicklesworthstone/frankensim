//! Manufactured Cp/rho laws, not measured specimen data.
use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyValue, Provenance,
    UncertaintyModel,
};
use fs_qty::{Density, Temperature};

const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);

/// The absolute temperature that every span query sends. The claims'
/// temperature validity declares the same quantity: fs-matdb compares the two
/// descriptors exactly, and an untyped axis is not a wildcard.
fn kelvin() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn card(revision: u32, cp_slope: f64, varying_density: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, centre, slope) in [
        ("specific_heat_capacity", CP, 400.0, cp_slope),
        ("alternate_capacity", CP, 400.0, cp_slope),
        ("density", Density::DIMS, 6000.0, if varying_density { -100.0 } else { 0.0 }),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, centre - 20.0 * slope), (320.0, centre + 20.0 * slope)],
                dims,
            },
            validity: ValidityDomain::unconstrained().with_quantity("T", kelvin(), 280.0, 320.0)
                .with("pressure", 100_000.0, 200_000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured single-phase caloric join".into(),
                license: "test-only".into(), artifact: None,
            },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {
        chemistry: "manufactured caloric material".into(), phase: "single phase fixture".into(),
        process: "analytical reference".into(), revision,
    }, claims, vec![]).unwrap()
}

#[allow(clippy::too_many_arguments)]
fn span(card: &MaterialCard, limits: [f64; 2], pressure: f64, cp_key: &str,
    phase: SolidLiquidPhase, datum: f64, error: f64) -> ResolvedHeatCapacityCurve
{
    let point = |t| QueryPoint::new()
        .with_quantity("T", kelvin(), t).unwrap()
        .with("pressure", pressure).unwrap();
    ResolvedHeatCapacityCurve::try_new(card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new(cp_key, CP),
        density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
        lower: &point(limits[0]), upper: &point(limits[1]),
        selection: SelectionPolicy::SingleClaimOnly,
    }, phase, datum, error).unwrap()
}

fn pair(card: &MaterialCard, phase: SolidLiquidPhase, error: f64)
    -> (ResolvedHeatCapacityCurve, ResolvedHeatCapacityCurve)
{
    (span(card, [280.0, 300.0], 100_000.0, "specific_heat_capacity", phase, 1_234_567.0, error),
     span(card, [300.0, 320.0], 100_000.0, "specific_heat_capacity", phase, -7654.0, error))
}

#[test]
fn g1_original_cp_wall_exchange_crosses_joined_spans_without_chord_error() {
    for cp_slope in [-10.0, 0.0, 10.0] {
        let material = card(0, cp_slope, false);
        let (lower, upper) = pair(&material, SolidLiquidPhase::Liquid, 0.001);
        let joined = lower.try_join_at_reference(&upper, 1e12).unwrap();
        let whole = span(&material, [280.0, 320.0], 100_000.0,
            "specific_heat_capacity", SolidLiquidPhase::Liquid, 0.0, 1e6);
        // The numerical owner traverses original source spans, even when the
        // optional public enthalpy chart contains hundreds of subdivisions.
        assert_eq!(joined.source_heat_capacity_knots().len(), 3);
        assert_eq!(whole.source_heat_capacity_knots().len(), 2);
        for (initial, expected, wall) in [(285.0_f64, 313.0_f64, 330.0_f64),
            (315.0, 287.0, 270.0)] {
            // Independent primitive for Cp(T)=400+s*(T-300). Select the target
            // and integrate analytically instead of reproducing the inverse.
            let exposure = (400.0 + cp_slope * (wall - 300.0))
                * ((wall - initial) / (wall - expected)).ln()
                - cp_slope * (expected - initial);
            let expected_heat = (expected - initial)
                * (400.0 + 0.5 * cp_slope * (expected + initial - 600.0));
            for law in [&joined, &whole] {
                let exchange = law.linear_heat_exchange(initial, wall, exposure).unwrap();
                assert!((exchange.temperature_k - expected).abs() < 2e-12);
                assert!((exchange.specific_enthalpy_change_j_kg - expected_heat).abs() < 2e-9);
                assert!((exposure * (wall - exchange.mean_temperature_k) - expected_heat).abs() < 2e-9);
                let gain = (wall - expected) / (wall - initial);
                assert!((exchange.inlet_enthalpy_gain - gain).abs() < 2e-14);
            }
        }
    }
}

#[test]
fn g0_original_cp_exchange_retains_support_and_weak_heat_independently_of_datum() {
    let material = card(0, 10.0, false);
    let ordinary = span(&material, [280.0, 320.0], 100_000.0,
        "specific_heat_capacity", SolidLiquidPhase::Solid, 0.0, 1e6);
    let shifted = span(&material, [280.0, 320.0], 100_000.0,
        "specific_heat_capacity", SolidLiquidPhase::Solid, 1e12, 0.001);
    let before = shifted.clone();
    let weak = ordinary.linear_heat_exchange(300.0, 330.0, 1e-8).unwrap();
    assert!(weak.specific_enthalpy_change_j_kg > 0.0);
    assert_eq!(shifted.linear_heat_exchange(300.0, 330.0, 1e-8).unwrap(), weak);
    for (temperature, wall, exposure) in [(279.0, 330.0, 1.0),
        (300.0, 330.0, 1e6), (300.0, 270.0, 1e6),
        (300.0, 330.0, -1.0), (300.0, f64::NAN, 1.0)] {
        assert!(shifted.linear_heat_exchange(temperature, wall, exposure).is_err());
    }
    for (wall, exposure) in [(330.0, 0.0), (300.0, 1000.0)] {
        let exchange = shifted.linear_heat_exchange(300.0, wall, exposure).unwrap();
        assert_eq!(exchange.temperature_k, 300.0);
        assert_eq!(exchange.specific_enthalpy_change_j_kg, 0.0);
    }
    assert_eq!(shifted, before);
}

#[test]
fn joined_calorimetry_preserves_exact_reference_and_source_sensible_density() {
    let material = card(0, 0.0, true);
    for phase in [SolidLiquidPhase::Solid, SolidLiquidPhase::Liquid] {
        let (lower, upper) = pair(&material, phase, 1e-6);
        let before = (lower.clone(), upper.clone());
        let joined = lower.try_join_at_reference(&upper, 0.0).unwrap();
        let chart = joined.curve();
        assert_eq!(chart.material_card_identity(), material.content_hash());
        assert_eq!(chart.knots().len(), 3);
        assert_eq!(chart.knots()[1].temperature_k, 300.0);
        assert_eq!(chart.knots()[1].specific_enthalpy_j_kg.to_bits(), 0.0_f64.to_bits());
        assert_eq!(chart.density_law_on_interval(0).unwrap(), PhaseIntervalDensityLaw::AffineDensity);
        assert_eq!(chart.density_law_on_interval(1).unwrap(), PhaseIntervalDensityLaw::AffineDensity);
        for (h, t, rho) in [(-8000.0, 280.0, 8000.0), (-4000.0, 290.0, 7000.0),
            (0.0, 300.0, 6000.0), (4000.0, 310.0, 5000.0), (8000.0, 320.0, 4000.0)] {
            let state = chart.state_at_specific_enthalpy(h).unwrap();
            assert_eq!(state.temperature_k(), t);
            assert_eq!(state.bulk_density_kg_m3(), rho);
            assert_eq!(state.phase(), phase);
        }
        assert_eq!(joined.receipts().len(), lower.receipts().len() + upper.receipts().len() - 2);
        assert_eq!(&joined.receipts()[..lower.receipts().len()], lower.receipts());
        assert_eq!(&joined.receipts()[lower.receipts().len()..], &upper.receipts()[2..]);
        assert_eq!(chart.source_receipts(), joined.receipts());
        assert_eq!((lower, upper), before);
    }
}

#[test]
fn affine_capacity_integrals_and_fixed_mass_work_survive_rebasing() {
    let material = card(0, 1.0, true);
    let (lower, upper) = pair(&material, SolidLiquidPhase::Solid, 1e6);
    let joined = lower.try_join_at_reference(&upper, 0.0).unwrap();
    assert_eq!(joined.curve().knots().len(), 3);
    for knot in joined.curve().knots() {
        let dt = knot.temperature_k - 300.0;
        let exact = 400.0 * dt + 0.5 * dt * dt;
        assert!((knot.specific_enthalpy_j_kg - exact).abs() < 3e-9);
    }
    for temperature in [287.5, 295.0, 305.0, 312.5] {
        let dt = temperature - 300.0;
        let h = 400.0 * dt + 0.5 * dt * dt;
        let state = joined.curve().values_at_specific_enthalpy(h).unwrap();
        assert!((state.temperature_k - temperature).abs() < 2e-12);
        assert!((state.bulk_density_kg_m3 - (6000.0 - 100.0 * dt)).abs() < 2e-10);
        assert!(
            (joined
                .curve()
                .temperature_derivative_at_specific_enthalpy(h)
                .unwrap()
                - 1.0 / (400.0 + dt))
                .abs()
                < 2e-17
        );
    }
    let initial = joined.curve().state_at_specific_enthalpy(0.0).unwrap();
    let step = joined.curve().advance_isobaric_heat(initial, 2.0, 100_000.0, 16_400.0).unwrap();
    assert!((step.state.temperature_k() - 320.0).abs() < 1e-10);
    let work = 100_000.0 * 2.0 * (1.0 / 4000.0 - 1.0 / 6000.0);
    assert!((step.boundary_work_j - work).abs() < 1e-12);
    assert!((step.internal_energy_change_j + step.boundary_work_j - 16_400.0).abs() < 1e-10);
    assert_eq!(step.state.material_card_identity(), material.content_hash());
    assert_eq!(step.state.phase_curve_identity(), joined.curve().identity());
}

#[test]
fn numerical_parity_does_not_admit_foreign_sources_contexts_or_phases() {
    let material = card(0, 0.0, false);
    // A later revision of the same state, with the same claims. assemble
    // starts a lineage at revision 0, so the successor comes from supersede.
    let foreign = MaterialCard::supersede(&material, material.claims().clone(), vec![]).unwrap();
    let (lower, upper) = pair(&material, SolidLiquidPhase::Solid, 1e-6);
    let before = lower.clone();
    for incompatible in [
        span(&foreign, [300.0, 320.0], 100_000.0, "specific_heat_capacity", SolidLiquidPhase::Solid, 0.0, 1e-6),
        span(&material, [300.0, 320.0], 200_000.0, "specific_heat_capacity", SolidLiquidPhase::Solid, 0.0, 1e-6),
        span(&material, [300.0, 320.0], 100_000.0, "alternate_capacity", SolidLiquidPhase::Solid, 0.0, 1e-6),
        span(&material, [300.0, 320.0], 100_000.0, "specific_heat_capacity", SolidLiquidPhase::Liquid, 0.0, 1e-6),
        span(&material, [301.0, 320.0], 100_000.0, "specific_heat_capacity", SolidLiquidPhase::Solid, 0.0, 1e-6),
        span(&material, [299.0, 320.0], 100_000.0, "specific_heat_capacity", SolidLiquidPhase::Solid, 0.0, 1e-6),
    ] {
        assert!(lower.try_join_at_reference(&incompatible, 0.0).is_err());
    }
    for datum in [f64::NAN, f64::INFINITY, -f64::INFINITY, f64::MAX] {
        assert!(lower.try_join_at_reference(&upper, datum).is_err());
    }
    assert_eq!(lower, before);
    let joined = lower.try_join_at_reference(&upper, 0.0).unwrap();
    assert!(joined.curve().advance_specific_energy(
        lower.curve().state_at_specific_enthalpy(lower.curve().knots()[0].specific_enthalpy_j_kg).unwrap(),
        0.0,
    ).is_err());
    assert!(joined.curve().state_at_specific_enthalpy(-8001.0).is_err());
    assert!(joined.curve().state_at_specific_enthalpy(8001.0).is_err());
}

#[test]
fn joining_keeps_the_combined_native_knot_budget() {
    let material = card(0, 1.0, false);
    let (lower, upper) = pair(&material, SolidLiquidPhase::Solid, 1e-5);
    assert!(lower.curve().knots().len() <= MAX_HEAT_CAPACITY_ENTHALPY_KNOTS);
    assert!(upper.curve().knots().len() <= MAX_HEAT_CAPACITY_ENTHALPY_KNOTS);
    assert!(lower.curve().knots().len() + upper.curve().knots().len() - 1 > MAX_HEAT_CAPACITY_ENTHALPY_KNOTS);
    assert!(lower.try_join_at_reference(&upper, 0.0).is_err());
}

#[test]
fn independent_old_datums_do_not_change_the_joined_physics_or_source_identity() {
    let material = card(0, 0.0, false);
    let (lower, upper) = pair(&material, SolidLiquidPhase::Solid, 1e-6);
    let first = lower.try_join_at_reference(&upper, 0.0).unwrap();
    let lower = span(&material, [280.0, 300.0], 100_000.0, "specific_heat_capacity", SolidLiquidPhase::Solid, 0.0, 1e-6);
    let upper = span(&material, [300.0, 320.0], 100_000.0, "specific_heat_capacity", SolidLiquidPhase::Solid, 0.0, 1e-6);
    let second = lower.try_join_at_reference(&upper, 0.0).unwrap();
    assert_eq!(first, second);
    assert_eq!(first.curve().identity(), second.curve().identity());
    let shifted = lower.try_join_at_reference(&upper, 123.0).unwrap();
    assert_ne!(first.curve().identity(), shifted.curve().identity());
    assert_eq!(shifted.curve().state_at_specific_enthalpy(123.0).unwrap().temperature_k(), 300.0);
    let signed_zero = lower.try_join_at_reference(&upper, -0.0).unwrap();
    assert_eq!(signed_zero.curve().knots()[1].specific_enthalpy_j_kg.to_bits(), (-0.0_f64).to_bits());
}
