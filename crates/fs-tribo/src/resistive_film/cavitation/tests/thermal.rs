//! G1/G3 conservative endpoint material changes, independent of a heat model.
use super::*;

#[test]
fn g1_drained_thermal_expansion_uses_original_mass_and_new_viscosity() {
    for alpha in [0.0, 1e-6] {
        let geometry = geometry(false).with_barus_viscosity(alpha).unwrap();
        let film = CavitatingFilm::new(geometry, liquid(), config()).unwrap();
        let state = film.initial_state(&[0.0]).unwrap();
        let next_liquid = CavitationLiquid {
            density_kg_m3: 800.0,
            cavitation_pressure_pa: 6000.0,
            ..liquid()
        };
        let next = film
            .advance_with_liquid(&state, &[0.0], &[0.0], 0.1, next_liquid, &[0.005], || false)
            .unwrap();
        // Fixed gap: the density decrease must expel the actual excess liquid.
        // The reference Poiseuille conductance at the reservoir is 2e-10.
        let expected_flow = (1000.0 / 800.0 - 1.0) * 1e-6 / 0.1;
        let equivalent_gauge = expected_flow / 2e-10;
        let gauge = if alpha == 0.0 {
            equivalent_gauge
        } else {
            -(1.0 - alpha * equivalent_gauge).ln() / alpha
        };
        assert_eq!(next.step.state.configuration(), state.configuration());
        assert_eq!(next.step.state.liquid_fraction(), &[1.0]);
        close(next.step.channel_volume_flow_m3_s[0], expected_flow);
        close(next.step.state.pressure_pa()[0], 100_000.0 + gauge);
        close(next.step.generalized_force[0], 0.001 * gauge);
        assert_eq!(next.step.report.liquid_mass_before_kg, 1000.0 * 1e-6);
        close(next.step.report.liquid_mass_after_kg, 800.0 * 1e-6);
        close(next.step.report.boundary_mass_out_kg, 200.0 * 1e-6);
        assert!(next.step.report.mass_residual_kg.abs() < 2e-18);
        assert_eq!(next.film.liquid(), next_liquid);
        assert_eq!(next.film.geometry().viscosity_pa_s(), Some(0.005));
        assert_eq!(
            next.film.geometry().pressure_viscosity_coefficient_per_pa(),
            alpha
        );
        assert!(Arc::ptr_eq(&film.geometry.cells, &next.film.geometry.cells));
        assert!(Arc::ptr_eq(
            &film.geometry.channels,
            &next.film.geometry.channels
        ));
    }
}

#[test]
fn g1_sealed_cavity_heating_and_cooling_preserve_mass_without_resetting_pressure() {
    let mut film = CavitatingFilm::new(geometry(true), liquid(), config()).unwrap();
    let mut state = film.initialize(&[0.0], &[0.25]).unwrap();
    let mut accepted_mass = 1000.0 * (1e-6 * 0.25);
    for (density, pressure, viscosity) in [
        (800.0, 6000.0, 0.005),
        (1200.0, 500.0, 0.03),
        (1000.0, 2000.0, 0.01),
    ] {
        let saved = state.clone();
        let next = film
            .advance_with_liquid(
                &state,
                &[0.0],
                &[0.0],
                0.1,
                CavitationLiquid {
                    density_kg_m3: density,
                    cavitation_pressure_pa: pressure,
                    ..liquid()
                },
                &[viscosity],
                || false,
            )
            .unwrap();
        assert_eq!(state, saved);
        assert_eq!(next.step.state.pressure_pa(), &[pressure]);
        close(next.step.state.liquid_fraction()[0], 250.0 / density);
        // No recomputation of old mass using the new density, even by rounding.
        assert_eq!(next.step.report.liquid_mass_before_kg, accepted_mass);
        close(next.step.report.liquid_mass_after_kg, accepted_mass);
        assert_eq!(next.step.report.boundary_mass_out_kg, 0.0);
        assert!(
            film.advance(&next.step.state, &[0.0], &[0.0], 0.1, || false)
                .is_err()
        );
        assert!(
            next.film
                .advance(&state, &[0.0], &[0.0], 0.1, || false)
                .is_err()
        );
        accepted_mass = next.step.report.liquid_mass_after_kg;
        state = next.step.state;
        film = next.film;
    }
}

#[test]
fn g1_thermal_upwind_transport_closes_each_cell_at_the_new_density() {
    let film = two_cells();
    let state = film.initialize(&[0.0], &[0.2, 0.1]).unwrap();
    let next = film
        .advance_with_liquid(
            &state,
            &[0.0],
            &[-1e-6, 1e-6, 1e-6],
            1.0,
            CavitationLiquid {
                density_kg_m3: 800.0,
                cavitation_pressure_pa: 6000.0,
                ambient_pressure_pa: 6000.0,
                ambient_liquid_fraction: 0.4,
            },
            &[0.005, 0.01, 0.02],
            || false,
        )
        .unwrap();
    // dt*Q/V=1; old filling contributes rho_old/rho_new, incoming
    // boundary filling contributes its actual endpoint-density mass.
    let theta0 = (1.25 * 0.2 + 0.4) / 2.0;
    let theta1 = (1.25 * 0.1 + theta0) / 2.0;
    close(next.step.state.liquid_fraction()[0], theta0);
    close(next.step.state.liquid_fraction()[1], theta1);
    let flow = &next.step.channel_volume_flow_m3_s;
    close(flow[0], -0.4e-6);
    close(flow[1], theta0 * 1e-6);
    close(flow[2], theta1 * 1e-6);
    let cell0 = 800.0 * (1e-6 * theta0 + flow[0] + flow[1]) - 1000.0 * 1e-6 * 0.2;
    let cell1 = 800.0 * (1e-6 * theta1 - flow[1] + flow[2]) - 1000.0 * 1e-6 * 0.1;
    assert!(cell0.abs().max(cell1.abs()) < 2e-18);
    close(
        next.step.report.boundary_mass_out_kg,
        800.0 * (flow[0] + flow[2]),
    );
    assert!(next.step.report.mass_residual_kg.abs() < 2e-18);
}

#[test]
fn g1_barus_refill_references_new_rupture_pressure_and_reservoir_viscosity() {
    let alpha = 2e-6;
    let film = CavitatingFilm::new(
        geometry(false).with_barus_viscosity(alpha).unwrap(),
        liquid(),
        config(),
    )
    .unwrap();
    let state = film.initialize(&[0.0], &[0.25]).unwrap();
    let next = film
        .advance_with_liquid(
            &state,
            &[0.0],
            &[0.0],
            0.01,
            CavitationLiquid {
                density_kg_m3: 800.0,
                cavitation_pressure_pa: 6000.0,
                ..liquid()
            },
            &[0.02],
            || false,
        )
        .unwrap();
    // Independent integral of exp(-alpha*(p-p_ambient))/mu_ambient.
    let inflow = 5e-11 * ((alpha * 94_000.0).exp() - 1.0) / alpha;
    let mass_before = 1000.0 * (1e-6 * 0.25);
    let mass_after = mass_before + 800.0 * 0.01 * inflow;
    assert_eq!(next.step.state.pressure_pa(), &[6000.0]);
    close(next.step.channel_volume_flow_m3_s[0], -inflow);
    close(
        next.step.state.liquid_fraction()[0],
        mass_after / (800.0 * 1e-6),
    );
    assert_eq!(next.step.report.liquid_mass_before_kg, mass_before);
    close(next.step.report.liquid_mass_after_kg, mass_after);
    close(next.step.generalized_force[0], -94.0);
}

#[test]
fn g3_unchanged_liquid_and_viscosity_only_transitions_keep_the_original_history() {
    for alpha in [0.0, 1e-6] {
        let film = CavitatingFilm::new(
            geometry(false).with_barus_viscosity(alpha).unwrap(),
            liquid(),
            config(),
        )
        .unwrap();
        let initial = film.initial_state(&[0.0]).unwrap();
        let mut state = initial;
        for q in [-0.001, 0.0, 0.0001] {
            let expected = film.advance(&state, &[q], &[0.0], 0.01, || false).unwrap();
            let next = film
                .advance_with_liquid(&state, &[q], &[0.0], 0.01, liquid(), &[0.01], || false)
                .unwrap();
            assert_eq!(next.step, expected);
            state = next.step.state;
        }
        let expected = film
            .with_viscosity(0.02)
            .unwrap()
            .advance(&state, &[0.0002], &[0.0], 0.01, || false)
            .unwrap();
        let next = film
            .advance_with_liquid(&state, &[0.0002], &[0.0], 0.01, liquid(), &[0.02], || false)
            .unwrap();
        assert_eq!(next.step, expected);
    }
}

#[test]
fn g0_refused_thermal_transition_keeps_both_accepted_owners_unchanged() {
    let film = CavitatingFilm::new(geometry(false), liquid(), config()).unwrap();
    let state = film.initialize(&[0.0], &[0.25]).unwrap();
    let saved = state.clone();
    let next_liquid = CavitationLiquid {
        density_kg_m3: 800.0,
        cavitation_pressure_pa: 6000.0,
        ..liquid()
    };
    for invalid in [
        CavitationLiquid {
            density_kg_m3: 0.0,
            ..next_liquid
        },
        CavitationLiquid {
            cavitation_pressure_pa: 200_000.0,
            ..next_liquid
        },
    ] {
        assert!(
            film.advance_with_liquid(&state, &[0.0], &[0.0], 0.01, invalid, &[0.005], || false,)
                .is_err()
        );
    }
    assert!(
        film.advance_with_liquid(&state, &[0.0], &[0.0], 0.01, next_liquid, &[], || false,)
            .is_err()
    );
    let mut polls = 0;
    film.advance_with_liquid(&state, &[0.0], &[0.0], 0.01, next_liquid, &[0.005], || {
        polls += 1;
        false
    })
    .unwrap();
    let mut cancelled_polls = 0;
    assert!(
        film.advance_with_liquid(&state, &[0.0], &[0.0], 0.01, next_liquid, &[0.005], || {
            cancelled_polls += 1;
            cancelled_polls == polls
        },)
            .is_err()
    );
    assert_eq!(cancelled_polls, polls);
    assert_eq!(state, saved);
    assert_eq!(film.liquid(), liquid());
    assert_eq!(film.geometry().viscosity_pa_s(), Some(0.01));
    // A full sealed pocket cannot accept thermal expansion without a real
    // compressible-liquid owner. There is no pressure clip or hidden drainage.
    let sealed = CavitatingFilm::new(geometry(true), liquid(), config()).unwrap();
    let full = sealed.initial_state(&[0.0]).unwrap();
    let saved_full = full.clone();
    assert!(
        sealed
            .advance_with_liquid(&full, &[0.0], &[0.0], 0.01, next_liquid, &[0.005], || false,)
            .is_err()
    );
    assert_eq!(full, saved_full);
}
