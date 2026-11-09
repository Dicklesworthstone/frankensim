//! The native caloric inverse owns heat; a source endpoint pin is admission only.
use super::*;

#[test]
fn native_wall_trace_outside_source_support_refuses_despite_admitted_enthalpy_pin() {
    let fluid_card = fluid_card();
    let (_, fluid) = source_pair(&fluid_card);
    let solid = solid(false);
    let fixture_wall = wall_caloric(&solid);
    let original = initial(&fluid, &fixture_wall, 300.0);
    // The large, valid original enthalpy datum makes this source endpoint's
    // inverse round trip differ from its exact source-temperature coordinate.
    let wall = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        solid.content_hash(),
        SolidLiquidPhase::Solid,
        [
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 1.0e7,
                temperature_k: 280.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            },
            EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: 1.0e7 + 40000.0,
                temperature_k: 320.0,
                liquid_mass_fraction: 0.0,
                bulk_density_kg_m3: 1000.0,
            },
        ]
        .to_vec(),
    )
    .unwrap();
    let lower = 305.23456789;
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    let layer = layer(&solid, &wall, &point, [lower, 315.0]);
    let endpoint_h = layer.enthalpy_range_j_kg()[0];
    assert_eq!(layer.temperature_at_enthalpy(endpoint_h).unwrap(), lower);

    // Reuse the graph fixture's exact cells/contact with the original wall law
    // and its own endpoint enthalpy. No alternate temperature field is stored.
    let mut events: Vec<_> = original
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            if cell.id == CellId(6) {
                cell.specific_h_j_kg = endpoint_h;
            }
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    events.extend(original.links().values().map(|&link| TimedEvent {
        time_s: 0.0,
        event: Event::Contact(link),
    }));
    let thermal = DepositionSimulation::new(
        events,
        &[
            (MaterialId(1), fluid.curve() as &dyn ThermalLaw),
            (MaterialId(2), &wall as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 16,
            max_links: 16,
            max_events: 32,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let native = thermal.temperature_k(CellId(6)).unwrap();
    assert!(native < lower, "native={native:.17e}, source={lower:.17e}");
    assert!(
        layer
            .conductivity_curve()
            .conductivity_w_m_k(native)
            .is_err()
    );
    let before_h = thermal.cells()[&CellId(6)].specific_h_j_kg.to_bits();
    for surface_feedback in [false, true] {
        let result = if surface_feedback {
            PreparedGraphConjugateLinks::new_with_surface_feedback(
                &thermal,
                vec![binding(&layer)],
                1,
                || Ok(()),
            )
        } else {
            PreparedGraphConjugateLinks::new(&thermal, vec![binding(&layer)], 1, || Ok(()))
        };
        assert!(matches!(
            result,
            Err(GraphConjugateError::Thermal(
                DepositionError::MaterialDomain
            ))
        ));
        assert_eq!(
            thermal.cells()[&CellId(6)].specific_h_j_kg.to_bits(),
            before_h
        );
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.links()[&9].conductance_w_k, 0.0);
    }
}
