//! G0/G1/G3/G4: supplied heat shares the original mass and latent-energy solve.
use super::*;

fn independent_rate(t: f64, bath: GasFilm) -> f64 {
    let saturation = 4000.0 + 100.0 * (t - 300.0);
    bath.mass_transfer_coefficient_m_s * bath.total_pressure_pa * 0.018
        / (MOLAR_GAS_CONSTANT * bath.temperature_k)
        * ((bath.total_pressure_pa - bath.vapor_partial_pressure_pa)
            / (bath.total_pressure_pa - saturation))
            .ln()
}

#[test]
fn supplied_heat_feeds_latent_transfer_without_a_sensible_only_overshoot() {
    for vapor_pressure in [0.0, 6000.0] {
        let bath = film(vapor_pressure);
        let duration = 2000.0;
        let dm = independent_rate(300.0, bath) * liquid().area_m2 * duration;
        let heat = dm * 2.4e6;
        if vapor_pressure == 0.0 {
            assert!(300.0 + heat / (liquid().mass_kg * 4000.0) > 330.0);
        }
        for datum in [0.0, 1e7, -1e7] {
            let initial = LiquidSurfaceState {
                specific_enthalpy_j_kg: datum,
                ..liquid()
            };
            let step = advance_surface_with_heat(
                &law(), bath, initial, duration, heat, controls(datum),
                |h| Ok(300.0 + (h - datum) / 4000.0), || Ok(()),
            ).unwrap();
            assert_eq!(step.surface_temperature_k, 300.0);
            assert_eq!(step.liquid.specific_enthalpy_j_kg, datum);
            assert_eq!(step.supplied_heat_j, heat);
            near(step.transferred_mass_kg, dm, 1e-16);
            near(step.latent_energy_j, heat, 1e-8);
            near(
                step.liquid.mass_kg * step.liquid.specific_enthalpy_j_kg
                    + step.vapor_energy_j
                    - initial.mass_kg * initial.specific_enthalpy_j_kg,
                heat, 1e-7,
            );
            assert_eq!(step.transferred_mass_kg.is_sign_positive(), heat.is_sign_positive());
        }
    }
}

#[test]
fn zero_heat_preserves_original_solver_and_sealed_heat_impulses_are_caloric() {
    for vapor_pressure in [0.0, 6000.0] {
        let original = advance_surface(
            &law(), film(vapor_pressure), liquid(), 20.0, controls(0.0),
            |h| Ok(300.0 + h / 4000.0), || Ok(()),
        ).unwrap();
        let explicit_zero = advance_surface_with_heat(
            &law(), film(vapor_pressure), liquid(), 20.0, 0.0, controls(0.0),
            |h| Ok(300.0 + h / 4000.0), || Ok(()),
        ).unwrap();
        assert_eq!(original, explicit_zero);
    }
    let sealed = GasFilm { mass_transfer_coefficient_m_s: 0.0, ..film(0.0) };
    for (bath, duration) in [(sealed, 0.0), (sealed, 20.0), (film(0.0), 0.0)] {
        for heat in [-400.0, 400.0] {
            let step = advance_surface_with_heat(
                &law(), bath, liquid(), duration, heat, controls(0.0),
                |h| Ok(300.0 + h / 4000.0), || Ok(()),
            ).unwrap();
            near(step.surface_temperature_k, 300.0 + heat / 400.0, 2e-10);
            assert_eq!(step.liquid.mass_kg, liquid().mass_kg);
            assert_eq!(step.transferred_mass_kg, 0.0);
            assert_eq!(step.vapor_energy_j, 0.0);
            assert_eq!(step.latent_energy_j, 0.0);
            near(step.liquid.mass_kg * step.liquid.specific_enthalpy_j_kg, heat, 1e-8);
        }
    }
}

#[test]
fn heating_and_cooling_reverse_phase_flux_with_endpoint_latent_heat() {
    struct VariableLatent;
    impl EvaporationLaw for VariableLatent {
        fn temperature_range_k(&self) -> [f64; 2] { [280.0, 330.0] }
        fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
            let mut p = law().properties(t)?;
            p.latent_heat_j_kg = 2.4e6 - 2000.0 * (t - 300.0);
            Ok(p)
        }
    }
    for (vapor_pressure, target, duration) in [
        (4500.0, 306.0, 20.0),
        (3500.0, 294.0, 20.0),
        (0.0, 315.0, 1000.0),
    ] {
        let bath = film(vapor_pressure);
        let dm = independent_rate(target, bath) * liquid().area_m2 * duration;
        let latent = 2.4e6 - 2000.0 * (target - 300.0);
        let heat = liquid().mass_kg * 4000.0 * (target - 300.0) + dm * latent;
        if vapor_pressure != 0.0 {
            assert!(independent_rate(300.0, bath) * dm < 0.0);
        }
        let step = advance_surface_with_heat(
            &VariableLatent, bath, liquid(), duration, heat, controls(0.0),
            |h| Ok(300.0 + h / 4000.0), || Ok(()),
        ).unwrap();
        near(step.surface_temperature_k, target, 2e-9);
        near(step.transferred_mass_kg, dm, 5e-13);
        near(step.latent_energy_j, dm * latent, 5e-7);
        near(
            step.liquid.mass_kg * step.liquid.specific_enthalpy_j_kg
                + step.vapor_energy_j, heat, 5e-7,
        );
    }
    // The upper source sample boils at 5000 Pa; the heated endpoint is 305 K.
    let bath = GasFilm {
        total_pressure_pa: 5000.0,
        vapor_partial_pressure_pa: 4000.0,
        ..film(0.0)
    };
    let dm = independent_rate(305.0, bath) * liquid().area_m2 * 20.0;
    let heat = 0.1 * 4000.0 * 5.0 + dm * 2.4e6;
    let step = advance_surface_with_heat(
        &law(), bath, liquid(), 20.0, heat, controls(0.0),
        |h| Ok(300.0 + h / 4000.0), || Ok(()),
    ).unwrap();
    near(step.surface_temperature_k, 305.0, 2e-9);
    near(step.transferred_mass_kg, dm, 2e-14);
}

#[test]
fn supplied_heat_inventory_retry_remains_finite_after_flux_reversal() {
    for vapor_pressure in [4000.0, 4500.0] {
        let bath = film(vapor_pressure);
        assert!(independent_rate(300.0, bath) <= 0.0);
        let duration = 1e6;
        let dm = independent_rate(310.0, bath) * liquid().area_m2 * duration;
        assert!(dm > liquid().mass_kg);
        let heat = liquid().mass_kg * 4000.0 * 10.0 + dm * 2.4e6;
        let Err(EvaporationError::InventoryLimit { maximum_duration_s }) =
            advance_surface_with_heat(
                &law(), bath, liquid(), duration, heat, controls(0.0),
                |h| Ok(300.0 + h / 4000.0), || Ok(()),
            )
        else { panic!("a heated exhausted endpoint must return an inventory retry"); };
        let expected = duration * liquid().mass_kg / dm;
        assert!(maximum_duration_s.is_finite() && maximum_duration_s > 0.0);
        assert!(maximum_duration_s < duration);
        near(maximum_duration_s, expected, expected * 1e-8);
    }
}

#[test]
fn supplied_heat_retains_source_budget_cancellation_and_finite_acceptance() {
    for heat in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            advance_surface_with_heat(
                &law(), film(0.0), liquid(), 20.0, heat, controls(0.0),
                |h| Ok(300.0 + h / 4000.0), || Ok(()),
            ),
            Err(EvaporationError::Invalid(_))
        ));
    }
    assert_eq!(
        advance_surface_with_heat(
            &law(), film(0.0), liquid(), 20.0, 2000.0,
            SolverConfig { max_iterations: 2, ..controls(0.0) },
            |h| Ok(300.0 + h / 4000.0), || Ok(()),
        ),
        Err(EvaporationError::NonConverged),
    );
    let mut polls = 0;
    assert_eq!(
        advance_surface_with_heat(
            &law(), film(0.0), liquid(), 20.0, 2000.0, controls(0.0),
            |h| Ok(300.0 + h / 4000.0),
            || {
                polls += 1;
                if polls == 3 { Err(EvaporationError::Cancelled) } else { Ok(()) }
            },
        ),
        Err(EvaporationError::Cancelled),
    );
    let source_error = EvaporationError::Invalid("original caloric source refused");
    assert_eq!(
        advance_surface_with_heat(
            &law(), film(0.0), liquid(), 20.0, 2000.0, controls(0.0),
            |h| if h > 0.0 { Err(source_error) } else { Ok(300.0 + h / 4000.0) },
            || Ok(()),
        ),
        Err(source_error),
    );
    // Each control is finite, but their combined energy allowance is not.
    assert!(matches!(
        advance_surface_with_heat(
            &law(), film(4000.0), liquid(), 0.0, f64::MAX,
            SolverConfig {
                relative_tolerance: 0.5,
                absolute_energy_tolerance_j: f64::MAX,
                ..controls(0.0)
            },
            |h| Ok(300.0 + h / 4000.0), || Ok(()),
        ),
        Err(EvaporationError::Unrepresentable(_))
    ));
}

#[test]
fn supplied_power_refines_to_independent_continuous_mass_energy_balance() {
    let power = 20.0;
    let mut reference_mass = liquid().mass_kg;
    let mut reference_temperature = 300.0;
    let dt = 0.001;
    for _ in 0..20_000 {
        let rate = independent_rate(reference_temperature, film(0.0)) * liquid().area_m2;
        let mid_mass = reference_mass - 0.5 * dt * rate;
        let mid_temperature = reference_temperature
            + 0.5 * dt * (power - rate * 2.4e6) / (reference_mass * 4000.0);
        let mid_rate = independent_rate(mid_temperature, film(0.0)) * liquid().area_m2;
        reference_mass -= dt * mid_rate;
        reference_temperature += dt * (power - mid_rate * 2.4e6) / (mid_mass * 4000.0);
    }
    let run = |count: usize| {
        let mut state = liquid();
        let mut temperature = 300.0;
        let dt = 20.0 / count as f64;
        for _ in 0..count {
            let step = advance_surface_with_heat(
                &law(), film(0.0), state, dt, power * dt, controls(0.0),
                |h| Ok(300.0 + h / 4000.0), || Ok(()),
            ).unwrap();
            state = step.liquid;
            temperature = step.surface_temperature_k;
        }
        (state.mass_kg, temperature)
    };
    let coarse = run(10);
    let fine = run(20);
    assert!((fine.1 - reference_temperature).abs() < 0.55 * (coarse.1 - reference_temperature).abs());
    assert!((fine.0 - reference_mass).abs() < 0.55 * (coarse.0 - reference_mass).abs());
}
