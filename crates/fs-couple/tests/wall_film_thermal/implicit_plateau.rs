//! Native wall fusion enthalpy supplies evaporation while wall T stays fixed.
use super::fixture::*;
use fs_blake3::ContentHash;
use fs_couple::finite_gas_thermal::wall_film::implicit::WallFilmImplicitControls;
use fs_material::evaporation::gas::wall_film::GasWallFilmInventory;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, UniformEnthalpyStepInput,
};

#[test]
fn native_wall_freezing_supplies_complete_film_evaporation_at_constant_plateau_temperature() {
    const WALL_LATENT: f64 = 100_000.0;
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x84; 32]),
        [
            (280.0, -20.0 * CP_WALL, 0.0),
            (300.0, 0.0, 0.0),
            (300.0, WALL_LATENT, 1.0),
            (340.0, WALL_LATENT + 40.0 * CP_WALL, 1.0),
        ]
        .map(
            |(temperature_k, specific_enthalpy_j_kg, liquid_mass_fraction)| EnthalpyPhaseKnot {
                temperature_k,
                specific_enthalpy_j_kg,
                liquid_mass_fraction,
                bulk_density_kg_m3: 1000.0,
            },
        )
        .to_vec(),
    )
    .unwrap();
    let input = UniformEnthalpyStepInput {
        initial: curve.state_at_specific_enthalpy(0.5 * WALL_LATENT).unwrap(),
        ..input(&curve, 300.0, 0.1)
    };
    let initial_wall = input.initial;
    let source = source();
    let gas = gas(&source, 290.0, 0.0);
    let initial_gas = gas.clone();
    let binding = binding(&source);
    let film_mass = 1e-7;
    let film = GasWallFilmInventory::Wet {
        mass_kg: film_mass,
        specific_enthalpy_j_kg: liquid_h(300.0).unwrap(),
    };
    let body = body(0.0);
    let environment = environment(&gas, &source, film, &binding, &body);
    let next = scope(|cx, _| {
        environment
            .advance_implicit(cx, input, WallFilmImplicitControls::default())
            .unwrap()
    });
    let report = &next.report;

    // Fusion storage, not a large apparent Cp, supplies the entire 0.244 J.
    // No prescribed heating, convection, or radiation pays this latent debit.
    let released = film_mass * LATENT;
    assert_eq!(initial_wall.temperature_k(), 300.0);
    assert_eq!(next.state.temperature_k(), 300.0);
    near(
        next.state.specific_enthalpy_j_kg(),
        initial_wall.specific_enthalpy_j_kg() - released / input.mass_kg,
        1e-6,
    );
    near(
        next.state.liquid_mass_fraction(),
        initial_wall.liquid_mass_fraction() - released / (input.mass_kg * WALL_LATENT),
        1e-11,
    );
    assert!(next.state.liquid_mass_fraction() < initial_wall.liquid_mass_fraction());
    assert_eq!(next.state.phase_curve_identity(), curve.identity());
    assert_eq!(report.next_film, GasWallFilmInventory::Dry);
    assert_eq!(report.next_film.specific_enthalpy_j_kg(), None);
    assert_eq!(report.phase.transferred_mass_kg, film_mass);
    assert_eq!(report.phase.duration_s, input.duration_s);
    let exhaustion = report.phase.dry_out_time_s.unwrap();
    assert!(exhaustion > 0.0 && exhaustion < input.duration_s);
    near(report.wall_heat_into_phase_system_j, released, 1e-9);
    near(report.phase.latent_energy_j, released, 1e-10);
    near(next.external_heat_j, -released, 1e-9);
    assert_eq!(report.solid_to_gas_heat_j, 0.0);
    assert_eq!(report.external_gas_heat_j, 0.0);
    assert_eq!(report.radiation_heat_j, 0.0);

    // The original NASA9 fixture has Cp/R=3.5. At fixed volume the newly
    // present vapor adds its actual Cv, while its stream carries Cp and p*v.
    let gas_capacity = 2.5 * 100_000.0 * GAS_VOLUME / 290.0;
    let expected_gas_t = 290.0
        + film_mass * (R * 290.0 + 3.5 * R * (300.0 - 290.0))
            / (gas_capacity + film_mass * 2.5 * R);
    near(
        report.next_gas.state().state().temperature,
        expected_gas_t,
        2e-7,
    );
    near(
        report.next_gas.species_mass_kg(&water()).unwrap(),
        film_mass,
        1e-14,
    );
    near(
        report.next_gas.species_mass_kg(&carrier()).unwrap(),
        gas.mass_kg(),
        1e-14,
    );
    let vapor_energy = film_mass * (-13.4e6 + 3.5 * R * 300.0);
    near(report.phase.energy_into_gas_j, vapor_energy, 1e-9);
    let wall_change = input.mass_kg
        * (next.state.specific_enthalpy_j_kg() - initial_wall.specific_enthalpy_j_kg());
    let gas_change = report.next_gas.mass_kg()
        * report.next_gas.state().specific_internal_energy_j_kg()
        - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    let film_change = -film_energy(film, &binding);
    near(gas_change, vapor_energy, 1e-9);
    near(wall_change + film_change + gas_change, 0.0, 1e-8);
    near(report.combined_energy_residual_j, 0.0, 1e-8);
    assert!(report.sum_absolute_energy_residual_j <= environment.solve.wall.energy_tolerance_j);
    near(report.mass_residual_kg, 0.0, 1e-14);
    assert_eq!(report.next_gas.volume_m3(), GAS_VOLUME);
    assert_eq!(
        report.next_gas.point().axes().get("wall-film-lot"),
        Some(&13.0)
    );
    assert_eq!(gas, initial_gas);
    assert_eq!(input.initial, initial_wall);
}
