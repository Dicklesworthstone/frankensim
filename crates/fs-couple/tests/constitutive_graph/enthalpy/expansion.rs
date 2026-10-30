//! Source geometry, electrical response and finite heat share the existing owners.
use super::*;
use fs_matdb::ClaimSelection;
use fs_material::thermal::transport::expansion::{
    FreeExpansionGeometryQuery, FreeExpansionStrainSources, PreparedFreeExpansion,
};
use fs_qty::QuantitySpec;

fn typed_point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::dimensional(Temperature::DIMS),
            temperature,
        )
        .unwrap()
}

fn expanding_material(slope: f64) -> (MaterialCard, ResolvedHeatCapacityCurve, PropertyKey) {
    let (original, _) = material(2.0, 2.0, 1.0);
    let mut claims = original.claims().clone();
    let strain = PropertyKey::new("engineering_length_change", Dims::NONE);
    claims
        .insert_claim(PropertyClaim {
            key: strain.clone(),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(200.0, -100.0 * slope), (800.0, 500.0 * slope)],
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained().with_quantity(
                "T",
                QuantitySpec::dimensional(Temperature::DIMS),
                200.0,
                800.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured free strain for graph/enthalpy verification".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &typed_point(200.0),
            upper: &typed_point(800.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-7,
    )
    .unwrap();
    (card, curve, strain)
}

#[test]
fn strain_data_closes_expanded_resistance_current_and_fixed_mass_enthalpy() {
    with_cx(|cx, _| {
        let mut endpoints = Vec::new();
        // Only the source length-change table varies. The final table makes
        // isotropic stretch proportional to resistivity, so resistance is constant.
        for expansion_slope in [0.0, 0.002, 0.005] {
            let (card, curve, strain) = expanding_material(expansion_slope);
            let reference = typed_point(300.0);
            let expansion = PreparedFreeExpansion::try_new(
                &card,
                FreeExpansionGeometryQuery {
                    reference: &reference,
                    strains: FreeExpansionStrainSources::Isotropic {
                        strain: &strain,
                        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
                    },
                    temperature_axis: "T",
                    lower: &typed_point(200.0),
                    upper: &typed_point(800.0),
                },
            )
            .unwrap();
            let original = resolve_conductor_state(
                &card,
                &reference,
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            let node = MaterialOhmicNode::try_new(&card, original, "T", [200.0, 800.0], 1.0, 1.0)
                .unwrap()
                .with_free_expansion(expansion, 0)
                .unwrap();
            let mut graph = ConstitutiveGraph::new();
            graph.add_node("resistor", Box::new(node.clone())).unwrap();
            let feedback = voltage_feedback();
            let mut thermal = timeline(
                &[(MaterialId(0), curve.curve() as &dyn ThermalLaw)],
                vec![birth(HEATER, MaterialId(0), 200.0)],
            );
            let reference_geometry = thermal.cells()[&HEATER].geometry;
            let mut previous_change = 0.0;
            let dt = 0.5;
            for step in 1..=2 {
                let committed_cells = thermal.cells().clone();
                let committed_energy = thermal.energy_j();
                let candidate = step_candidate(
                    cx,
                    &thermal,
                    f64::from(step) * dt,
                    feedback.input(&graph, &[]),
                    thermal_controls(),
                    &[0.0],
                    no_additional_heat,
                )
                .unwrap();
                // With C=mass*cp, R=(R0+b*x)/(1+a*x), x=T-Tref:
                // C*(x-x_old)*(R0+b*x) = dt*V²*(1+a*x).
                // This independent quadratic includes the actual source geometry.
                let capacity = 2.0;
                let r0 = 2.0;
                let supply = dt * VOLTAGE * VOLTAGE;
                let a = capacity * RESISTANCE_SLOPE;
                let b =
                    capacity * (r0 - RESISTANCE_SLOPE * previous_change) - supply * expansion_slope;
                let c = capacity * r0 * previous_change + supply;
                let expected_change = 2.0 * c / (b + (b * b + 4.0 * a * c).sqrt());
                let temperature = candidate.thermal().temperature_k(HEATER).unwrap();
                let current = candidate.graph().external[&key("resistor", "current")];
                let stretch = 1.0 + expansion_slope * expected_change;
                let resistance = (r0 + RESISTANCE_SLOPE * expected_change) / stretch;
                close(temperature, 300.0 + expected_change);
                close(current, VOLTAGE / resistance);
                let native = node.resolve_at_temperature_k(temperature).unwrap();
                close(native.length_m(), stretch);
                close(native.area_m2(), stretch * stretch);
                close(native.resistance_ohm(), resistance);
                let endpoint_power = VOLTAGE * current;
                let stored = capacity * (expected_change - previous_change);
                close(
                    candidate.graph().graph.outputs[&key("resistor", "heat")],
                    endpoint_power,
                );
                close(candidate.thermal().energy_j() - committed_energy, stored);
                close(stored, dt * endpoint_power);
                close(candidate.report().graph_heat_j, stored);
                close(candidate.report().thermal.external_heat_j, stored);
                close(candidate.report().other_external_heat_j, 0.0);
                // No repeated density*deformed-volume mass creation or second
                // thermal geometry owner; expansion remains absolute from Tref.
                close(candidate.thermal().mass_kg(), 1.0);
                assert_eq!(
                    candidate.thermal().cells()[&HEATER].geometry,
                    reference_geometry
                );
                assert_eq!(thermal.cells(), &committed_cells);
                assert_eq!(thermal.time_s(), f64::from(step - 1) * dt);
                assert!(candidate.graph().graph.next_state.is_empty());
                previous_change = expected_change;
                thermal = candidate.into_parts().0;
            }
            endpoints.push((
                thermal.temperature_k(HEATER).unwrap(),
                VOLTAGE
                    / node
                        .resolve_at_temperature_k(thermal.temperature_k(HEATER).unwrap())
                        .unwrap()
                        .resistance_ohm(),
            ));
        }
        assert!(
            endpoints
                .windows(2)
                .all(|pair| pair[1].0 > pair[0].0 && pair[1].1 > pair[0].1)
        );
        close(endpoints[2].0, 325.0);
        close(endpoints[2].1, 5.0);
    });
}
