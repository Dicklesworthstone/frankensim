//! Native solid/liquid latent storage remains part of the coupled heat graph.
use super::*;
use crate::evaporation::gas::implicit::adaptive::AdaptiveImplicitGasControls;
use crate::evaporation::gas::implicit::timeline::advance_deposition_with_gas_timeline;
use crate::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_blake3::ContentHash;

#[test]
fn unselected_native_pcm_melts_then_supplies_latent_heat_to_a_newborn_evaporating_liquid() {
    // The production phase curve and its original ThermalLaw implementation
    // own both the isothermal plateau and its liquid fraction throughout.
    let knot = |h, t, fraction| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: fraction,
        bulk_density_kg_m3: 1000.0,
    };
    let pcm = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x71; 32]),
        vec![
            knot(0.0, 290.0, 0.0),
            knot(20_000.0, 310.0, 0.0),
            knot(120_000.0, 310.0, 1.0),
            knot(140_000.0, 330.0, 1.0),
        ],
    )
    .unwrap();
    let mut wall = cell(1, 310.0, 10.0, 0.0012);
    wall.material_id = MaterialId(2);
    wall.mass_kg = 0.001;
    wall.specific_h_j_kg = 20_000.0;
    wall.geometry.volume_m3 = 1e-6;
    let off = Boundary {
        source_w: 0.0,
        ..wall.boundary
    };
    let born = cell(2, 300.0, 0.0, 0.0012);
    let born_energy = born.mass_kg * born.specific_h_j_kg;
    let initial = DepositionSimulation::new(
        vec![
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(wall),
            },
            TimedEvent {
                time_s: 2.0,
                event: Event::Birth(born),
            },
            TimedEvent {
                time_s: 2.0,
                event: Event::Contact(ThermalLink {
                    id: 9,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 0.0002,
                    conductance_w_k: 0.2,
                }),
            },
            TimedEvent {
                time_s: 2.0,
                event: Event::Boundary(CellId(1), off),
            },
        ],
        &[(MaterialId(1), &ORIGINAL_LIQUID), (MaterialId(2), &pcm)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 4,
            max_substeps: 256,
        },
    )
    .unwrap();
    let source = catalog(0.0);
    let initial_gas = gas(&source, 300.0, 1000.0, 0.001);
    let binding = binding(0.0);
    let surfaces = [surface(2, &binding, 0.02, 0.0)];
    let mut config = AdaptiveImplicitGasControls::new(controls());
    config.solve.thermal.maximum_iterations = 10_000;
    config.solve.thermal.maximum_evaluations = 100_000;
    config.solve.thermal.newton.absolute_tolerance = 1e-10;
    config.initial_dt_s = 1.0;
    config.maximum_dt_s = 1.0;

    let prefix = advance_deposition_with_gas_timeline(
        &initial,
        &surfaces,
        &initial_gas,
        &source,
        2.0,
        config,
        |thermal, liquids, gas| {
            assert!(liquids.is_empty());
            assert_eq!(gas, &initial_gas);
            assert_eq!(thermal.temperature_k(CellId(1))?, 310.0);
            Ok(coefficients(()))
        },
        || Ok(()),
    )
    .unwrap();
    let melted = pcm
        .values_at_specific_enthalpy(prefix.thermal().cells()[&CellId(1)].specific_h_j_kg)
        .unwrap();
    near(melted.specific_enthalpy_j_kg, 40_000.0, 1e-7);
    assert_eq!(melted.temperature_k, 310.0);
    near(melted.liquid_mass_fraction, 0.2, 1e-12);
    near(prefix.report().external_heat_j, 20.0, 1e-9);
    assert_eq!(prefix.gas(), &initial_gas);
    assert_eq!(prefix.report().gas_iterations, 0);
    assert!(prefix.report().surface_totals.is_empty());
    assert!(prefix.thermal().cells().contains_key(&CellId(2)));
    assert!(prefix.thermal().links().contains_key(&9));

    let mut idle_calls = 0usize;
    let mut active_calls = 0usize;
    let full = advance_deposition_with_gas_timeline(
        &initial,
        &surfaces,
        &initial_gas,
        &source,
        4.0,
        config,
        |thermal, liquids, gas| {
            assert_eq!(thermal.temperature_k(CellId(1))?, 310.0);
            if liquids.is_empty() {
                idle_calls += 1;
                assert!(thermal.time_s() <= 2.0);
                assert_eq!(gas, &initial_gas);
            } else {
                active_calls += 1;
                assert!(thermal.time_s() > 2.0);
                assert_eq!(liquids.len(), 1);
                assert!(thermal.links().contains_key(&9));
                near(liquids[&CellId(2)].area_m2, 0.001, 1e-15);
            }
            Ok(coefficients(thermal.time_s()))
        },
        || Ok(()),
    )
    .unwrap();
    assert!(idle_calls > 0 && active_calls > 0);
    assert_eq!(full.payload(), Some(&4.0));
    assert_eq!(full.report().events, 4);
    near(full.report().external_heat_j, 20.0, 1e-8);
    let final_wall = &full.thermal().cells()[&CellId(1)];
    let freezing = pcm
        .values_at_specific_enthalpy(final_wall.specific_h_j_kg)
        .unwrap();
    assert_eq!(freezing.temperature_k, 310.0);
    assert!(freezing.liquid_mass_fraction > 0.0 && freezing.liquid_mass_fraction < 0.2);
    assert_eq!(final_wall.mass_kg, 0.001);
    let released = 0.001 * (40_000.0 - final_wall.specific_h_j_kg);
    near(
        released,
        0.001 * 100_000.0 * (0.2 - freezing.liquid_mass_fraction),
        1e-10,
    );
    assert!(released > 0.0);

    let phase = &full.report().surface_totals[&CellId(2)];
    assert!(phase.transferred_mass_kg > 0.0 && phase.latent_energy_j > 0.0);
    assert_eq!(full.report().surface_totals.len(), 1);
    near(
        phase.latent_energy_j,
        LATENT * phase.transferred_mass_kg,
        1e-9,
    );
    near(
        full.gas().species_mass_kg(&vapor()).unwrap()
            - initial_gas.species_mass_kg(&vapor()).unwrap(),
        phase.transferred_mass_kg,
        1e-14,
    );
    // After the heater switches off, the original plateau is the only energy
    // supplier. Its freezing energy must reach the born liquid plus real gas.
    let final_liquid = &full.thermal().cells()[&CellId(2)];
    let gas_energy_change = full.gas().mass_kg()
        * full.gas().state().specific_internal_energy_j_kg()
        - initial_gas.mass_kg() * initial_gas.state().specific_internal_energy_j_kg();
    near(
        final_liquid.mass_kg * final_liquid.specific_h_j_kg - born_energy + gas_energy_change,
        released,
        1e-7,
    );
    assert!(
        full.report().absolute_energy_residual_j <= config.solve.thermal.total_energy_tolerance_j
    );
    assert!(full.report().mass_residual_kg.abs() <= config.solve.phase.gas.mass_tolerance_kg);
}
