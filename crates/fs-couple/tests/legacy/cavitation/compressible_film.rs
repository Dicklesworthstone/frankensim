//! G1/G3/G4 source-backed compressive storage and thermal expansion.
//! Manufactured reference curves do not establish measured-fluid accuracy.

use super::{
    AMBIENT_PA, AREA_M2, CHANNEL_GAP_M, INITIAL_GAP_M, film, near, point, quantity, source,
    try_film,
};
use fs_couple::cavitation::{
    CavitationFilmEnvironment, CavitationSolveConfig, MaterialCavitatingFilm,
};
use fs_material::cavitation::{ResolvedCavitationState, resolve_cavitation_state};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::semantic::QuantityKind;
use fs_tribo::resistive_film::cavitation::CavitationReservoir;
use fs_tribo::resistive_film::{FilmCell, FilmChannel, FilmLimits, GapPort};

fn sealed(source: &ResolvedCavitationState) -> MaterialCavitatingFilm<'_, ResolvedCavitationState> {
    MaterialCavitatingFilm::new(
        source,
        vec![FilmCell {
            area_m2: AREA_M2,
            gap: GapPort {
                reference_m: INITIAL_GAP_M,
                closure: vec![1.0],
            },
        }],
        vec![FilmChannel {
            from: 0,
            to: None,
            width_m: 1e-3,
            length_m: 1e-2,
            gap: GapPort {
                reference_m: 0.0,
                closure: vec![0.0],
            },
        }],
        1,
        FilmLimits {
            minimum_cell_gap_m: 1e-6,
            maximum_gap_m: 1e-3,
            maximum_pressure_pa: 1e9,
        },
        CavitationFilmEnvironment {
            ambient_pressure_pa: AMBIENT_PA,
            ambient_liquid_fraction: 1.0,
            maximum_vapor_density_ratio: 0.01,
            solve: CavitationSolveConfig {
                maximum_iterations: 32,
                absolute_mass_tolerance_kg: 1e-17,
                relative_mass_tolerance: 1e-11,
            },
        },
    )
    .unwrap()
}

#[test]
fn sealed_source_heating_raises_pressure_and_cooling_recovers_compressed_mass() {
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let source = source(basis, false, 0);
        let cold = source.resolve(300.0);
        let hot = cold.at_query_point(&source.card, &point(400.0)).unwrap();
        let bulk = 1e8;
        let original = sealed(&cold)
            .with_isothermal_bulk_modulus("pressure", bulk)
            .unwrap();
        let pressure = AMBIENT_PA + 2e6;
        let state = original
            .initialize_with_pressure(&[0.0], &[1.0], &[pressure])
            .unwrap();
        let saved = state.clone();
        let volume = AREA_M2 * INITIAL_GAP_M;
        let mass = 1000.0 * ((pressure - AMBIENT_PA) / bulk).exp() * volume;
        let next = original
            .advance_with_temperature_source(&hot, &state, &[0.0], &[0.0], 0.1, || false)
            .unwrap();
        near(
            next.step.state.pressure_pa()[0],
            pressure + bulk * (1000.0_f64 / 800.0).ln(),
        );
        assert_eq!(next.step.state.liquid_fraction(), &[1.0]);
        assert_eq!(next.step.channel_mass_flow_kg_s, [0.0]);
        near(next.step.report.liquid_mass_before_kg, mass);
        near(next.step.report.liquid_mass_after_kg, mass);
        assert_eq!(next.step.report.boundary_mass_out_kg, 0.0);
        assert!(next.step.report.mass_residual_kg.abs() < 2e-17);
        assert!(std::ptr::eq(next.film.source(), &hot));
        assert_eq!(next.film.material(), hot.state());
        assert_eq!(next.film.pressure_axis(), "pressure");
        assert_eq!(next.film.film().bulk_modulus_pa(), bulk);
        let mu = if basis == LiquidViscosityBasis::Dynamic {
            0.0005
        } else {
            800.0 * 0.5e-6
        };
        assert_eq!(next.film.film().geometry().viscosity_pa_s(), Some(mu));
        assert!(hot.material().property("specific_heat_capacity").is_none());
        let cooled = next
            .film
            .advance_with_temperature_source(&cold, &next.step.state, &[0.0], &[0.0], 0.1, || false)
            .unwrap();
        near(cooled.step.state.pressure_pa()[0], pressure);
        near(cooled.step.report.liquid_mass_before_kg, mass);
        near(cooled.step.report.liquid_mass_after_kg, mass);
        assert_eq!(cooled.step.state.liquid_fraction(), &[1.0]);
        assert!(std::ptr::eq(cooled.film.source(), &cold));
        assert_eq!(state, saved);
        assert_eq!(original.material(), cold.state());
        // Changed reference density requires the matching returned owner.
        assert!(
            original
                .advance(&next.step.state, &[0.0], &[0.0], 0.1, || false)
                .is_err()
        );
        assert!(
            next.film
                .advance(&next.step.state, &[0.0], &[0.0], 0.1, || false)
                .is_ok()
        );
    }
}

#[test]
fn source_heating_with_changed_stiffness_retains_original_compressed_inventory() {
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let source = source(basis, false, 0);
        let cold = source.resolve(300.0);
        let hot = cold.at_query_point(&source.card, &point(400.0)).unwrap();
        let old_bulk = 1e8;
        let new_bulk = 4e7;
        let pressure = AMBIENT_PA + 2e6;
        let original = sealed(&cold)
            .with_isothermal_bulk_modulus("pressure", old_bulk)
            .unwrap();
        let state = original
            .initialize_with_pressure(&[0.0], &[1.0], &[pressure])
            .unwrap();
        let mass = 1000.0 * ((pressure - AMBIENT_PA) / old_bulk).exp() * AREA_M2 * INITIAL_GAP_M;
        let next = original
            .advance_with_temperature_source_and_bulk_modulus(
                &hot,
                &state,
                &[0.0],
                &[0.0],
                0.1,
                new_bulk,
                || false,
            )
            .unwrap();
        near(
            next.step.state.pressure_pa()[0],
            AMBIENT_PA
                + new_bulk * ((1000.0_f64 / 800.0).ln() + (pressure - AMBIENT_PA) / old_bulk),
        );
        near(next.step.report.liquid_mass_before_kg, mass);
        near(next.step.report.liquid_mass_after_kg, mass);
        assert_eq!(next.step.state.liquid_fraction(), &[1.0]);
        assert_eq!(next.step.channel_mass_flow_kg_s, [0.0]);
        assert_eq!(next.film.film().bulk_modulus_pa(), new_bulk);
        assert_eq!(next.film.material(), hot.state());
        assert!(std::ptr::eq(next.film.source(), &hot));
        assert_eq!(next.film.pressure_axis(), "pressure");
        let cooled = next
            .film
            .advance_with_temperature_source_and_bulk_modulus(
                &cold,
                &next.step.state,
                &[0.0],
                &[0.0],
                0.1,
                old_bulk,
                || false,
            )
            .unwrap();
        near(cooled.step.state.pressure_pa()[0], pressure);
        near(cooled.step.report.liquid_mass_after_kg, mass);
        assert_eq!(cooled.film.film().bulk_modulus_pa(), old_bulk);
        assert!(std::ptr::eq(cooled.film.source(), &cold));
        assert!(
            original
                .advance(&next.step.state, &[0.0], &[0.0], 0.1, || false)
                .is_err()
        );
        assert_eq!(original.film().bulk_modulus_pa(), old_bulk);
        assert_eq!(original.material(), cold.state());
    }
}

#[test]
fn sourced_viscosity_and_endpoint_modulus_control_drained_expansion_with_retained_supply() {
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let source = source(basis, false, 0);
        let cold = source.resolve(300.0);
        let hot = cold.at_query_point(&source.card, &point(400.0)).unwrap();
        let bulk = 2e7;
        let endpoint_bulk = 7e6;
        let supply = CavitationReservoir {
            pressure_pa: 130_000.0,
            liquid_fraction: 1.0,
        };
        let original = film(&cold)
            .with_isothermal_bulk_modulus("pressure", bulk)
            .unwrap()
            .with_reservoirs(&[Some(supply)])
            .unwrap();
        let initial_pressure = 150_000.0;
        let state = original
            .initialize_with_pressure(&[0.0], &[1.0], &[initial_pressure])
            .unwrap();
        let dt = 1.0;
        let next = original
            .advance_with_temperature_source_and_bulk_modulus(
                &hot,
                &state,
                &[0.0],
                &[0.0],
                dt,
                endpoint_bulk,
                || false,
            )
            .unwrap();
        let volume = AREA_M2 * INITIAL_GAP_M;
        let mu = if basis == LiquidViscosityBasis::Dynamic {
            0.0005
        } else {
            800.0 * 0.5e-6
        };
        let conductance = 1e-3 * CHANNEL_GAP_M.powi(3) / (12.0 * mu * 1e-2);
        // Independent closed form: rho*V + dt*G*B*(rho-rho_supply) = m_old.
        // This integrates rho(p)/mu through the passage, not endpoint rho*Q.
        let mass_before = 1000.0 * ((initial_pressure - AMBIENT_PA) / bulk).exp() * volume;
        let rho_supply = 800.0 * ((supply.pressure_pa - AMBIENT_PA) / endpoint_bulk).exp();
        let rho = (mass_before + dt * conductance * endpoint_bulk * rho_supply)
            / (volume + dt * conductance * endpoint_bulk);
        let pressure = AMBIENT_PA + endpoint_bulk * (rho / 800.0).ln();
        let mass_flow = conductance * endpoint_bulk * (rho - rho_supply);
        near(next.step.state.pressure_pa()[0], pressure);
        near(next.step.channel_mass_flow_kg_s[0], mass_flow);
        near(next.step.report.liquid_mass_before_kg, mass_before);
        near(next.step.report.liquid_mass_after_kg, rho * volume);
        near(next.step.report.boundary_mass_out_kg, dt * mass_flow);
        assert!(next.step.report.mass_residual_kg.abs() < 2e-17);
        assert_eq!(
            next.film.film().reservoir_conditions().collect::<Vec<_>>(),
            [Some(supply)]
        );
        assert_eq!(next.film.film().liquid().ambient_pressure_pa, AMBIENT_PA);
        assert_eq!(next.film.film().liquid().cavitation_pressure_pa, 20_000.0);
        assert_eq!(next.film.film().bulk_modulus_pa(), endpoint_bulk);
    }
}

#[test]
fn source_pressure_reference_and_original_lineage_are_required_before_mass_transition() {
    let source = source(LiquidViscosityBasis::Dynamic, false, 0);
    let cold = source.resolve(300.0);
    let uncompressed = film(&cold);
    for axis in ["", "missing", "temperature"] {
        assert!(
            uncompressed
                .with_isothermal_bulk_modulus(axis, 1e8)
                .is_err()
        );
    }
    for bulk in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            uncompressed
                .with_isothermal_bulk_modulus("pressure", bulk)
                .is_err()
        );
    }
    let extra_nonpressure = resolve_cavitation_state(
        &source.card,
        &point(300.0)
            .with_quantity(
                "not_pressure",
                quantity(QuantityKind::AbsoluteTemperature),
                AMBIENT_PA,
            )
            .unwrap(),
        &source.requirements,
        "temperature",
        MaterialPropertySelection::PinnedByProperty(source.pins.clone()),
    )
    .unwrap();
    assert!(
        film(&extra_nonpressure)
            .with_isothermal_bulk_modulus("not_pressure", 1e8)
            .is_err()
    );
    let wrong_reference = cold
        .at_query_point(
            &source.card,
            &point(300.0).with("pressure", 110_000.0).unwrap(),
        )
        .unwrap();
    assert!(
        film(&wrong_reference)
            .with_isothermal_bulk_modulus("pressure", 1e8)
            .is_err()
    );
    let original = uncompressed
        .with_isothermal_bulk_modulus("pressure", 1e8)
        .unwrap();
    let state = original.initialize(&[0.0], &[0.2]).unwrap();
    let saved = state.clone();
    let foreign = super::source(LiquidViscosityBasis::Dynamic, false, 1).resolve(400.0);
    let changed_context = cold
        .at_query_point(
            &source.card,
            &point(400.0).with("pressure", 110_000.0).unwrap(),
        )
        .unwrap();
    let mut pins = source.pins.clone();
    pins.iter_mut()
        .find(|(name, _)| name == "dynamic_viscosity")
        .unwrap()
        .1 = source.other_viscosity;
    let changed_claim = source.resolve_with(400.0, pins);
    for candidate in [&foreign, &changed_context, &changed_claim] {
        assert!(
            original
                .advance_with_temperature_source(candidate, &state, &[0.0], &[0.0], 0.01, || false,)
                .is_err()
        );
        assert!(
            original
                .advance_with_temperature_source_and_bulk_modulus(
                    candidate,
                    &state,
                    &[0.0],
                    &[0.0],
                    0.01,
                    4e7,
                    || false,
                )
                .is_err()
        );
    }
    for bulk in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            original
                .advance_with_temperature_source_and_bulk_modulus(
                    &cold,
                    &state,
                    &[0.0],
                    &[0.0],
                    0.01,
                    bulk,
                    || false,
                )
                .is_err()
        );
    }
    let hot = cold.at_query_point(&source.card, &point(400.0)).unwrap();
    // The source reference-density cap passes, but rho at rupture is lower.
    assert!(
        try_film(&cold, 1e-4)
            .unwrap()
            .with_isothermal_bulk_modulus("pressure", 5e4)
            .is_err()
    );
    let rupture_dilute = try_film(&cold, 2.7e-4)
        .unwrap()
        .with_isothermal_bulk_modulus("pressure", 1e6)
        .unwrap();
    let rupture_state = rupture_dilute.initialize(&[0.0], &[0.2]).unwrap();
    assert!(
        hot.state().saturated_vapor().unwrap().density_kg_m3() / hot.state().liquid_density_kg_m3()
            < 2.7e-4
    );
    assert!(
        rupture_dilute
            .advance_with_temperature_source(&hot, &rupture_state, &[0.0], &[0.0], 0.01, || false,)
            .is_err()
    );
    assert!(
        rupture_dilute
            .advance(&rupture_state, &[0.0], &[0.0], 0.01, || false)
            .is_ok()
    );
    let dilute_only = try_film(&cold, 1e-4)
        .unwrap()
        .with_isothermal_bulk_modulus("pressure", 1e8)
        .unwrap();
    let dilute_state = dilute_only.initialize(&[0.0], &[0.2]).unwrap();
    // The reference density remains admissible while the new rupture density
    // can violate the explicitly declared negligible-vapor-mass approximation.
    assert!(
        dilute_only
            .advance_with_temperature_source_and_bulk_modulus(
                &cold,
                &dilute_state,
                &[0.0],
                &[0.0],
                0.01,
                5e4,
                || false,
            )
            .is_err()
    );
    assert!(
        dilute_only
            .advance_with_temperature_source(&hot, &dilute_state, &[0.0], &[0.0], 0.01, || false,)
            .is_err()
    );
    let incompressible_state = uncompressed.initial_state(&[0.0]).unwrap();
    assert!(
        original
            .advance(&incompressible_state, &[0.0], &[0.0], 0.01, || false)
            .is_err()
    );
    assert_eq!(state, saved);
    assert_eq!(original.material(), cold.state());
    assert!(
        original
            .advance(&state, &[0.0], &[0.0], 0.01, || false)
            .is_ok()
    );
}

#[test]
fn constant_source_matches_native_bits_and_late_cancellation_preserves_the_complete_pair() {
    let source = source(LiquidViscosityBasis::Dynamic, true, 0);
    let cold = source.resolve(300.0);
    let hot = cold.at_query_point(&source.card, &point(400.0)).unwrap();
    let original = film(&cold)
        .with_isothermal_bulk_modulus("pressure", 1e8)
        .unwrap();
    let state = original.initial_state(&[0.0]).unwrap();
    let saved = state.clone();
    let expected = original
        .advance(&state, &[-INITIAL_GAP_M], &[0.0], 0.01, || false)
        .unwrap();
    let mut polls = 0;
    let next = original
        .advance_with_temperature_source_and_bulk_modulus(
            &hot,
            &state,
            &[-INITIAL_GAP_M],
            &[0.0],
            0.01,
            1e8,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(next.step, expected);
    assert!(std::ptr::eq(next.film.source(), &hot));
    assert_eq!(next.film.material().temperature_k(), 400.0);
    assert!(polls > 2);
    let mut calls = 0;
    assert!(
        original
            .advance_with_temperature_source_and_bulk_modulus(
                &hot,
                &state,
                &[-INITIAL_GAP_M],
                &[0.0],
                0.01,
                1e8,
                || {
                    calls += 1;
                    calls == polls
                },
            )
            .is_err()
    );
    assert_eq!(calls, polls);
    assert_eq!(state, saved);
    assert_eq!(original.material(), cold.state());
    let retry = original
        .advance_with_temperature_source(&hot, &state, &[-INITIAL_GAP_M], &[0.0], 0.01, || false)
        .unwrap();
    assert_eq!(retry.step, expected);
}
