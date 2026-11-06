//! Original affine Cp(T), conductivity breakpoints and their exact chain rule.

use super::*;
use fs_material::phase::{
    HeatCapacityKnot, HeatCapacitySourceQuery,
    melting::{MeltingSourceQuery, ResolvedMeltingCurve},
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

fn conductivity_card() -> MaterialCard {
    card(vec![claim(
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, 2.0), (500.0, 10.0), (800.0, 4.0)],
            dims: KD,
        },
        [0.0, 1.0],
    )])
}

fn source_curve(
    card: &MaterialCard,
    capacities: [f64; 2],
    datum: f64,
) -> EquilibriumEnthalpyPhaseCurve {
    let source = [300.0, 800.0].map(|temperature_k| HeatCapacityKnot {
        temperature_k,
        specific_heat_capacity_j_kg_k: if temperature_k == 300.0 {
            capacities[0]
        } else {
            capacities[1]
        },
        bulk_density_kg_m3: 1000.0,
    });
    let curve = EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
        card.content_hash(),
        SolidLiquidPhase::Solid,
        datum,
        &source,
        1e9,
    )
    .unwrap();
    // The conductivity source's 500 K breakpoint lies strictly inside this
    // single original capacity interval; enthalpy-chord placement is wrong.
    assert_eq!(curve.knots().len(), 2);
    curve
}

fn capacity(capacities: [f64; 2], temperature: f64) -> f64 {
    capacities[0] + (capacities[1] - capacities[0]) * (temperature - 300.0) / 500.0
}

fn original_enthalpy(capacities: [f64; 2], datum: f64, temperature: f64) -> f64 {
    let dt = temperature - 300.0;
    datum + capacities[0] * dt + 0.5 * (capacities[1] - capacities[0]) / 500.0 * dt * dt
}

fn source_conductivity(card: &MaterialCard, temperature: f64) -> f64 {
    let point = QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            temperature,
        )
        .unwrap()
        .with_quantity("f", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap();
    card.claims()
        .query_typed(
            &PropertyKey::new("k", KD),
            &point,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap()
        .evidence
        .value
        .value
}

fn close_derivative(actual: f64, expected: f64) {
    assert!(
        (actual / expected - 1.0).abs() < 2e-12,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn assert_retained_capacities(
    path: &PreparedPhaseConductivity<'_>,
    index: usize,
    expected: [f64; 2],
) {
    let retained = path
        .sensible_heat_capacities_on_interval(index)
        .expect("the compiled interval must retain original source capacities");
    for (actual, expected) in retained.into_iter().zip(expected) {
        close_derivative(actual, expected);
    }
}

fn assert_source_sample(
    card: &MaterialCard,
    path: &PreparedPhaseConductivity<'_>,
    capacities: [f64; 2],
    datum: f64,
    temperature: f64,
    conductivity_temperature_slope: f64,
) {
    let h = original_enthalpy(capacities, datum, temperature);
    let sample = path.state_at_specific_enthalpy(h).unwrap();
    let cp = capacity(capacities, temperature);
    close(sample.temperature_k(), temperature);
    close(
        sample.conductivity_w_m_k(),
        source_conductivity(card, temperature),
    );
    assert_eq!(sample.liquid_mass_fraction(), 0.0);
    assert_eq!(sample.phase(), SolidLiquidPhase::Solid);
    close_derivative(sample.temperature_enthalpy_derivative_kg_k_j(), cp.recip());
    close_derivative(
        sample.conductivity_enthalpy_derivative_kg_w_m_k_j(),
        conductivity_temperature_slope / cp,
    );
    let native = ImplicitConductionLaw::state_at_specific_enthalpy(path, h).unwrap();
    close(native.temperature_k, temperature);
    close(
        native.conductivity_w_m_k,
        source_conductivity(card, temperature),
    );
    close_derivative(native.temperature_slope_kg_k_j, cp.recip());
    close_derivative(
        native.conductivity_slope_kg_w_m_k_j,
        conductivity_temperature_slope / cp,
    );
}

#[test]
fn affine_capacity_retains_original_source_breakpoint_enthalpy_and_transport_tangents() {
    let card = conductivity_card();
    for capacities in [[1000.0, 3500.0], [3500.0, 1000.0]] {
        for datum in [0.0, 1e12] {
            let curve = source_curve(&card, capacities, datum);
            let h = |t| original_enthalpy(capacities, datum, t);
            let path = prepare(&card, &curve, [h(300.0), h(800.0)], || Ok(())).unwrap();
            assert_eq!(
                path.knots(),
                &[[h(300.0), 2.0], [h(500.0), 10.0], [h(800.0), 4.0]]
            );
            assert_eq!(path.source_paths().len(), 1);
            assert_retained_capacities(&path, 0, [capacities[0], capacity(capacities, 500.0)]);
            assert_retained_capacities(&path, 1, [capacity(capacities, 500.0), capacities[1]]);
            for temperature in [300.0, 350.0, 450.0, 500.0, 550.0, 650.0, 800.0] {
                // Right derivative at the original 500 K source knot;
                // the 800 K upper endpoint keeps its left derivative.
                let slope = if temperature < 500.0 { 0.04 } else { -0.02 };
                assert_source_sample(&card, &path, capacities, datum, temperature, slope);
            }
            assert!(std::ptr::eq(path.caloric_curve(), &curve));
        }
    }
}

#[test]
fn clipped_caloric_conductivity_keeps_original_breakpoint_and_one_sided_source_law() {
    let card = conductivity_card();
    for capacities in [[1000.0, 3500.0], [3500.0, 1000.0]] {
        let datum = 1e12;
        let curve = source_curve(&card, capacities, datum);
        let h = |t| original_enthalpy(capacities, datum, t);
        let path = prepare(&card, &curve, [h(350.0), h(650.0)], || Ok(())).unwrap();
        assert_eq!(
            path.knots(),
            &[[h(350.0), 4.0], [h(500.0), 10.0], [h(650.0), 7.0]]
        );
        assert_retained_capacities(
            &path,
            0,
            [capacity(capacities, 350.0), capacity(capacities, 500.0)],
        );
        assert_retained_capacities(
            &path,
            1,
            [capacity(capacities, 500.0), capacity(capacities, 650.0)],
        );
        for temperature in [350.0, 400.0, 450.0, 500.0, 600.0, 650.0] {
            let slope = if temperature < 500.0 { 0.04 } else { -0.02 };
            assert_source_sample(&card, &path, capacities, datum, temperature, slope);
        }
        assert!(
            path.state_at_specific_enthalpy(h(350.0).next_down())
                .is_err()
        );
        assert!(path.state_at_specific_enthalpy(h(650.0).next_up()).is_err());

        // Ending the requested domain at a conductivity corner selects the
        // admitted interval's left derivative, while the full path uses right.
        let left = prepare(&card, &curve, [h(350.0), h(500.0)], || Ok(())).unwrap();
        assert_source_sample(&card, &left, capacities, datum, 500.0, 0.04);
    }
}

#[test]
fn prepared_solid_prefix_keeps_its_sensible_derivative_at_the_melting_boundary() {
    let absolute_temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let cp_dims = Dims([2, 0, -2, -1, 0, 0]);
    let density_dims = Dims([-3, 1, 0, 0, 0, 0]);
    let energy_dims = Dims([2, 0, -2, 0, 0, 0]);
    let curve_value = |dims, knots| PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots,
        dims,
    };
    let properties = vec![
        (
            "cp_solid",
            QuantitySpec::dimensional(cp_dims),
            [300.0, 400.0],
            curve_value(cp_dims, vec![(300.0, 1000.0), (400.0, 2000.0)]),
        ),
        (
            "cp_liquid",
            QuantitySpec::dimensional(cp_dims),
            [400.0, 500.0],
            curve_value(cp_dims, vec![(400.0, 2000.0), (500.0, 2000.0)]),
        ),
        (
            "rho",
            QuantitySpec::dimensional(density_dims),
            [300.0, 500.0],
            PropertyValue::Scalar {
                value: 1000.0,
                dims: density_dims,
            },
        ),
        (
            "k",
            QuantitySpec::dimensional(KD),
            [300.0, 500.0],
            curve_value(KD, vec![(300.0, 2.0), (500.0, 6.0)]),
        ),
        (
            "melting_temperature",
            absolute_temperature,
            [300.0, 500.0],
            PropertyValue::Scalar {
                value: 400.0,
                dims: Temperature::DIMS,
            },
        ),
        (
            "fusion_heat",
            QuantitySpec::dimensional(energy_dims),
            [300.0, 500.0],
            PropertyValue::Scalar {
                value: 20_000.0,
                dims: energy_dims,
            },
        ),
    ];
    let card = card(
        properties
            .into_iter()
            .map(|(name, quantity, bounds, value)| {
                let scalar = matches!(&value, PropertyValue::Scalar { .. });
                let mut property = claim(value, [0.0, 1.0]);
                property.key = PropertyKey::with_quantity(name, quantity);
                property.validity = ValidityDomain::unconstrained().with_quantity(
                    "T",
                    absolute_temperature,
                    bounds[0],
                    bounds[1],
                );
                property.interpolation = if scalar {
                    InterpolationPolicy::ConstantWithinValidity
                } else {
                    InterpolationPolicy::LinearInside
                };
                property
            })
            .collect(),
    );
    let point = |temperature| {
        QueryPoint::new()
            .with_quantity("T", absolute_temperature, temperature)
            .unwrap()
    };
    let points = [point(300.0), point(400.0), point(500.0)];
    let solid_cp = PropertyKey::new("cp_solid", cp_dims);
    let liquid_cp = PropertyKey::new("cp_liquid", cp_dims);
    let density = PropertyKey::new("rho", density_dims);
    let melting_temperature =
        PropertyKey::with_quantity("melting_temperature", absolute_temperature);
    let fusion_heat = PropertyKey::new("fusion_heat", energy_dims);
    let melting = ResolvedMeltingCurve::try_new(
        &card,
        MeltingSourceQuery {
            solid: HeatCapacitySourceQuery {
                heat_capacity: &solid_cp,
                density: &density,
                temperature_axis: "T",
                lower: &points[0],
                upper: &points[1],
                selection: SelectionPolicy::SingleClaimOnly,
            },
            liquid: HeatCapacitySourceQuery {
                heat_capacity: &liquid_cp,
                density: &density,
                temperature_axis: "T",
                lower: &points[1],
                upper: &points[2],
                selection: SelectionPolicy::SingleClaimOnly,
            },
            melting_point: &melting_temperature,
            latent_heat_fusion: &fusion_heat,
            transition: &points[1],
            transition_selection: SelectionPolicy::SingleClaimOnly,
            phase_selector: None,
        },
        0.0,
        1e6,
    )
    .unwrap();
    let solid_end = melting.plateau_enthalpy_bounds_j_kg()[0];
    assert_eq!(solid_end, 150_000.0);
    // The full chart correctly selects the latent interval to its right.
    assert_eq!(
        melting
            .curve()
            .temperature_derivative_at_specific_enthalpy(solid_end)
            .unwrap(),
        0.0
    );
    let source = PreparedPhaseConductivity::new(
        melting.curve(),
        ThermalMaterialPoint {
            material: &card,
            context: &points[0],
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        &PropertyKey::new("k", KD),
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        [0.0, solid_end],
        || Ok(()),
    )
    .unwrap();
    // This prepared solid prefix instead owns the LEFT sensible derivative.
    let sample = source.state_at_specific_enthalpy(solid_end).unwrap();
    assert_eq!(sample.temperature_k(), 400.0);
    assert_eq!(sample.conductivity_w_m_k(), 4.0);
    close_derivative(
        sample.temperature_enthalpy_derivative_kg_k_j(),
        1.0 / 2000.0,
    );
    close_derivative(
        sample.conductivity_enthalpy_derivative_kg_w_m_k_j(),
        0.02 / 2000.0,
    );
    let native = ImplicitConductionLaw::state_at_specific_enthalpy(&source, solid_end).unwrap();
    close_derivative(native.temperature_slope_kg_k_j, 1.0 / 2000.0);
    close_derivative(native.conductivity_slope_kg_w_m_k_j, 0.02 / 2000.0);
}
