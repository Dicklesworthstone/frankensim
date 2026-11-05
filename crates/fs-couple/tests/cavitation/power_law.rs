//! G1/G3/G4 source-dependent nonlinear transport on original liquid inventory.
use super::*;
use fs_couple::cavitation::MaterialCavitatingPowerLawFilm;
use fs_tribo::resistive_film::cavitation::{CavitationLiquid, CavitationReservoir};
use fs_tribo::resistive_film::power_law::{PowerLawControls, RegularizedPowerLaw};

fn controls() -> PowerLawControls {
    PowerLawControls {
        maximum_iterations: 96,
        maximum_line_search_steps: 48,
        absolute_flow_tolerance_m3_s: 1e-28,
        relative_flow_tolerance: 2e-12,
        absolute_power_tolerance_w: 1e-22,
        relative_power_tolerance: 2e-12,
    }
}

fn reference_law() -> RegularizedPowerLaw {
    RegularizedPowerLaw {
        consistency_pa_s_n: 0.03,
        flow_index: 0.5,
        minimum_shear_rate_per_s: 5.0,
    }
}

fn nonlinear(
    source: &ResolvedCavitationState,
) -> MaterialCavitatingPowerLawFilm<'_, ResolvedCavitationState> {
    film(source)
        .with_power_law(vec![reference_law()], controls())
        .unwrap()
}

// Independent dimensional integration for n=1/2 and a resolved Newtonian core:
// Q = w*h^2/(8*K^2) * [tau_w^2 + tau_core^4/(3*tau_w^2)].
// The larger quadratic root is the branch above tau_core. No native pressure
// solver, material viscosity mean or regularized-flow helper is used here.
fn thinning_pressure_drop(consistency: f64, flow: f64) -> f64 {
    let width = 1e-3;
    let length = 1e-2;
    let core_squared = consistency * consistency * reference_law().minimum_shear_rate_per_s;
    let b = 8.0 * consistency * consistency * flow / (width * CHANNEL_GAP_M.powi(2));
    let wall_squared = 0.5 * (b + (b * b - 4.0 * core_squared.powi(2) / 3.0).sqrt());
    assert!(wall_squared > core_squared);
    2.0 * length * wall_squared.sqrt() / CHANNEL_GAP_M
}

#[test]
fn g1_source_temperature_drives_nonlinear_pressure_and_conserves_expanding_mass() {
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let fixture = source(basis, false, 0);
        let cold = fixture.resolve(300.0);
        let hot = cold.at_query_point(&fixture.card, &point(400.0)).unwrap();
        let original = nonlinear(&cold);
        let state = original.initial_state(&[0.0]).unwrap();
        let saved = state.clone();
        let next = original
            .advance_with_temperature_source(&hot, &state, &[0.0], &[0.0], &[0.0], 0.5, || false)
            .unwrap();
        let volume = AREA_M2 * INITIAL_GAP_M;
        let expected_flow = (1000.0 / 800.0 - 1.0) * volume / 0.5;
        let ratio = if basis == LiquidViscosityBasis::Dynamic {
            0.5
        } else {
            0.4
        };
        let consistency = reference_law().consistency_pa_s_n * ratio;
        near(next.film.film().laws()[0].consistency_pa_s_n, consistency);
        assert_eq!(next.film.film().laws()[0].flow_index, 0.5);
        assert_eq!(next.film.film().laws()[0].minimum_shear_rate_per_s, 5.0);
        near(next.step.channel_volume_flow_m3_s[0], expected_flow);
        near(
            next.step.state.pressure_pa()[0],
            AMBIENT_PA + thinning_pressure_drop(consistency, expected_flow),
        );
        assert_eq!(next.step.report.liquid_mass_before_kg, 1000.0 * volume);
        near(next.step.report.liquid_mass_after_kg, 800.0 * volume);
        near(next.step.report.boundary_mass_out_kg, 200.0 * volume);
        assert_eq!(next.film.film().liquid().cavitation_pressure_pa, 20_000.0);
        assert_eq!(next.film.reference_material(), cold.state());
        assert_eq!(next.film.reference_laws(), &[reference_law()]);
        assert_eq!(next.film.material(), hot.state());
        assert!(std::ptr::eq(next.film.source(), &hot));
        assert_eq!(state, saved);
        assert_eq!(original.material(), cold.state());
    }
}

#[test]
fn g3_original_reference_laws_survive_repeated_temperature_and_inventory_continuation() {
    let fixture = source(LiquidViscosityBasis::Dynamic, false, 0);
    let cold = fixture.resolve(300.0);
    let warm = cold.at_query_point(&fixture.card, &point(350.0)).unwrap();
    let hot = cold.at_query_point(&fixture.card, &point(400.0)).unwrap();
    let mut model = nonlinear(&cold);
    let mut state = model.initial_state(&[0.0]).unwrap();
    let reference = reference_law();
    for source in [&hot, &warm, &cold, &warm, &hot, &cold] {
        let previous_mass = model.material().liquid_density_kg_m3()
            * (AREA_M2 * INITIAL_GAP_M * state.liquid_fraction()[0]);
        let next = model
            .advance_with_temperature_source(source, &state, &[0.0], &[0.0], &[0.0], 0.5, || false)
            .unwrap();
        let expected = reference.consistency_pa_s_n
            * (source.state().dynamic_viscosity_pa_s() / cold.state().dynamic_viscosity_pa_s());
        assert_eq!(next.film.film().laws()[0].consistency_pa_s_n, expected);
        assert_eq!(next.film.reference_laws(), &[reference]);
        assert_eq!(next.film.reference_material(), cold.state());
        assert_eq!(next.step.report.liquid_mass_before_kg, previous_mass);
        near(
            next.step.report.liquid_mass_after_kg + next.step.report.boundary_mass_out_kg,
            previous_mass,
        );
        model = next.film;
        state = next.step.state;
    }
    assert_eq!(model.film().laws(), &[reference]);
}

#[test]
fn g1_starved_sliding_uses_the_same_shifted_law_endpoint_density_and_explicit_supply() {
    let fixture = source(LiquidViscosityBasis::Dynamic, false, 0);
    let cold = fixture.resolve(300.0);
    let hot = cold.at_query_point(&fixture.card, &point(400.0)).unwrap();
    let supplied = [Some(CavitationReservoir {
        pressure_pa: 30_000.0,
        liquid_fraction: 1.0,
    })];
    let model = nonlinear(&cold).with_reservoirs(&supplied).unwrap();
    let state = model.initialize(&[0.0], &[0.2]).unwrap();
    let material = RegularizedPowerLaw {
        consistency_pa_s_n: 0.015,
        ..reference_law()
    };
    let expected = model
        .film()
        .advance_with_liquid(
            &state,
            &[0.0],
            &[1e-10],
            &[0.02],
            0.01,
            CavitationLiquid {
                density_kg_m3: 800.0,
                cavitation_pressure_pa: 20_000.0,
                ..model.film().liquid()
            },
            vec![material],
            || false,
        )
        .unwrap();
    let actual = model
        .advance_with_temperature_source(&hot, &state, &[0.0], &[1e-10], &[0.02], 0.01, || false)
        .unwrap();
    assert_eq!(
        actual.step.state.liquid_fraction(),
        expected.step.state.liquid_fraction()
    );
    assert_eq!(actual.step.state.pressure_pa(), &[20_000.0]);
    assert!(
        actual.step.state.liquid_fraction()[0] > 0.25
            && actual.step.state.liquid_fraction()[0] < 1.0
    );
    assert_eq!(
        actual.step.channel_volume_flow_m3_s,
        expected.step.channel_volume_flow_m3_s
    );
    assert_eq!(
        actual.step.generalized_force,
        expected.step.generalized_force
    );
    assert_eq!(actual.step.report, expected.step.report);
    assert_eq!(
        actual
            .film
            .film()
            .reservoir_conditions()
            .collect::<Vec<_>>(),
        supplied
    );
    let stationary = model
        .advance_with_temperature_source(&hot, &state, &[0.0], &[1e-10], &[0.0], 0.01, || false)
        .unwrap();
    assert_ne!(
        actual.step.channel_volume_flow_m3_s,
        stationary.step.channel_volume_flow_m3_s
    );
    assert_eq!(state.liquid_fraction(), &[0.2]);
    assert_eq!(state.pressure_pa(), &[2000.0]);
}

#[test]
fn g3_source_newtonian_limit_retains_the_original_endpoint_trajectory_exactly() {
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let fixture = source(basis, false, 0);
        let cold = fixture.resolve(300.0);
        let hot = cold.at_query_point(&fixture.card, &point(400.0)).unwrap();
        let native = film(&cold);
        let state = native.initialize(&[0.0], &[0.2]).unwrap();
        let law = RegularizedPowerLaw {
            consistency_pa_s_n: cold.state().dynamic_viscosity_pa_s(),
            flow_index: 1.0,
            minimum_shear_rate_per_s: 5.0,
        };
        let model = native.with_power_law(vec![law], controls()).unwrap();
        let expected = native
            .advance_with_temperature_source(&hot, &state, &[0.0], &[1e-10], 0.01, || false)
            .unwrap();
        let actual = model
            .advance_with_temperature_source(&hot, &state, &[0.0], &[1e-10], &[0.02], 0.01, || {
                false
            })
            .unwrap();
        assert_eq!(
            actual.film.film().laws()[0].consistency_pa_s_n,
            hot.state().dynamic_viscosity_pa_s()
        );
        assert_eq!(
            actual.step.state.liquid_fraction(),
            expected.step.state.liquid_fraction()
        );
        assert_eq!(
            actual.step.state.pressure_pa(),
            expected.step.state.pressure_pa()
        );
        assert_eq!(
            actual.step.channel_volume_flow_m3_s,
            expected.step.channel_volume_flow_m3_s
        );
        assert_eq!(
            actual.step.generalized_force,
            expected.step.generalized_force
        );
        assert_eq!(actual.step.report, expected.step.report);
    }
}

#[test]
fn g4_source_lineage_supply_domain_and_late_cancellation_keep_the_accepted_pair() {
    let fixture = source(LiquidViscosityBasis::Dynamic, false, 0);
    let cold = fixture.resolve(300.0);
    let hot = cold.at_query_point(&fixture.card, &point(400.0)).unwrap();
    let model = nonlinear(&cold);
    let state = model.initialize(&[0.0], &[0.2]).unwrap();
    let saved = state.clone();
    let foreign = source(LiquidViscosityBasis::Dynamic, false, 1).resolve(400.0);
    let changed_context = cold
        .at_query_point(
            &fixture.card,
            &point(400.0).with("pressure", 110_000.0).unwrap(),
        )
        .unwrap();
    let mut pins = fixture.pins.clone();
    pins.iter_mut()
        .find(|(name, _)| name == "dynamic_viscosity")
        .unwrap()
        .1 = fixture.other_viscosity;
    let changed_claim = fixture.resolve_with(400.0, pins);
    for candidate in [&foreign, &changed_context, &changed_claim] {
        assert!(
            model
                .advance_with_temperature_source(
                    candidate,
                    &state,
                    &[0.0],
                    &[0.0],
                    &[0.02],
                    0.01,
                    || false,
                )
                .is_err()
        );
        assert_eq!(state, saved);
    }
    let supplied = model
        .with_reservoirs(&[Some(CavitationReservoir {
            pressure_pa: 10_000.0,
            liquid_fraction: 1.0,
        })])
        .unwrap();
    assert!(
        supplied
            .advance_with_temperature_source(&hot, &state, &[0.0], &[0.0], &[0.02], 0.01, || false,)
            .is_err()
    );
    let mut polls = 0;
    model
        .advance_with_temperature_source(&hot, &state, &[0.0], &[0.0], &[0.02], 0.001, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls - 1, polls] {
        let mut count = 0;
        assert!(
            model
                .advance_with_temperature_source(
                    &hot,
                    &state,
                    &[0.0],
                    &[0.0],
                    &[0.02],
                    0.001,
                    || {
                        count += 1;
                        count == stop
                    },
                )
                .is_err()
        );
        assert_eq!(model.film().laws(), &[reference_law()]);
        assert_eq!(model.reference_laws(), &[reference_law()]);
        assert!(std::ptr::eq(model.source(), &cold));
        assert_eq!(state, saved);
    }
    model
        .advance(&state, &[0.0], &[0.0], &[0.02], 0.001, || false)
        .unwrap();
}
