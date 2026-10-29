use super::*;

fn controls(lo: f64, hi: f64) -> DenseIsobaricControls {
    DenseIsobaricControls {
        minimum_temperature_k: lo,
        maximum_temperature_k: hi,
        enthalpy_tolerance_j_kg: 2e-6,
        maximum_iterations: 96,
        maximum_entropy_seam_j_kg_k: 0.1,
        pressure: DensePressureSolve {
            pressure_tolerance_pa: 1e-4,
            density_tolerance_kg_m3: 1e-10,
            maximum_iterations: 256,
        },
    }
}
fn chart(p: f64, lo: f64, hi: f64) -> DenseIsobaricWater {
    DenseIsobaricWater::new(Pressure::new(p), controls(lo, hi), || false).unwrap()
}
fn at_t(chart: &DenseIsobaricWater, t: f64) -> DenseIsobaricState {
    chart
        .state_at_temperature(Temperature::new(t), || false)
        .unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}, tolerance {tolerance:.8e}"
    );
}
fn entropy_solve() -> WaterEntropySolve {
    WaterEntropySolve {
        entropy_tolerance_j_kg_k: 5e-9,
        maximum_iterations: 96,
    }
}

#[test]
fn supercritical_and_dense_liquid_enthalpy_entropy_round_trips_keep_raw_sources() {
    let chart = chart(25e6, 624.0, 675.0);
    assert!(chart.saturation().is_none());
    assert_eq!(chart.saturation_entropy_defect_j_kg_k(), 0.0);
    for t in [624.0, 630.0, 647.096, 650.0, 660.0, 670.0, 675.0] {
        let state = at_t(&chart, t);
        let raw = state.single_phase_source().unwrap();
        assert_eq!(
            state.specific_enthalpy_j_kg(),
            raw.specific_enthalpy().value()
        );
        assert_eq!(
            state.specific_internal_energy_j_kg(),
            raw.specific_internal_energy().value()
        );
        assert_eq!(
            state.specific_entropy_j_kg_k(),
            raw.specific_entropy().value()
        );
        assert_eq!(state.pressure_pa(), 25e6);
        assert!(
            state.source_pressure_residuals_pa()[0].abs()
                <= controls(624.0, 675.0).pressure.pressure_tolerance_pa
        );
        let h_inverse = chart
            .state_at_specific_enthalpy(state.specific_enthalpy_j_kg(), || false)
            .unwrap();
        let s_inverse = chart
            .state_at_specific_entropy(state.specific_entropy_j_kg_k(), entropy_solve(), || false)
            .unwrap();
        near(h_inverse.state().temperature_k(), t, 2e-8);
        near(s_inverse.state().temperature_k(), t, 2e-8);
        near(h_inverse.enthalpy_residual_j_kg(), 0.0, 2e-6);
        near(s_inverse.entropy_residual_j_kg_k(), 0.0, 5e-9);
        if t > CRITICAL_TEMPERATURE_K {
            assert_eq!(state.phase(), DenseIsobaricPhase::Supercritical);
            assert_eq!(state.vapor_mass_fraction(), None);
            assert_eq!(state.liquid_mass_fraction(), None);
            assert_eq!(state.vapor_volume_fraction(), None);
        }
        assert_eq!(state.vapor_quality(), None);
    }
}

#[test]
fn high_pressure_boiling_uses_mass_quality_actual_volumes_and_one_latent_account() {
    let chart = chart(20e6, 625.0, 649.0);
    let pair = chart.saturation().unwrap();
    let liquid = chart.state_at_quality(0.0, || false).unwrap();
    let vapor = chart.state_at_quality(1.0, || false).unwrap();
    assert!(pair.temperature().value() > 623.15);
    assert!(pair.temperature().value() < CRITICAL_TEMPERATURE_K);
    let latent = pair.latent_enthalpy().value();
    near(
        chart.saturation_entropy_defect_j_kg_k(),
        latent / pair.temperature().value() - pair.latent_entropy().value(),
        0.0,
    );
    let mass = 0.2;
    let quality = 0.4;
    let mixture = chart.state_at_quality(quality, || false).unwrap();
    let heat = mass * (mixture.specific_enthalpy_j_kg() - liquid.specific_enthalpy_j_kg());
    let step = chart
        .advance_heat(liquid, mass, heat, 2e-5, || false)
        .unwrap();
    near(step.state.vapor_quality().unwrap(), quality, 1e-15);
    assert_eq!(step.state.temperature_k(), pair.temperature().value());
    assert_eq!(step.state.phase(), DenseIsobaricPhase::LiquidVapor);
    near(
        step.phase_masses_kg.unwrap()[0],
        mass * (1.0 - quality),
        1e-16,
    );
    near(step.phase_masses_kg.unwrap()[1], mass * quality, 1e-16);
    let volume =
        (1.0 - quality) * liquid.specific_volume_m3_kg() + quality * vapor.specific_volume_m3_kg();
    near(mixture.specific_volume_m3_kg(), volume, 1e-18);
    near(mixture.bulk_density_kg_m3() * volume, 1.0, 2e-16);
    near(
        mixture.vapor_volume_fraction().unwrap(),
        quality * vapor.specific_volume_m3_kg() / volume,
        2e-16,
    );
    assert!(mixture.vapor_volume_fraction().unwrap() > quality);
    assert!(step.boundary_work_j > 1000.0);
    near(step.enthalpy_change_j, mass * quality * latent, 2e-5);
    near(
        step.internal_energy_change_j + step.boundary_work_j,
        heat,
        2e-5,
    );
    near(
        step.thermal_entropy_change_j_k,
        heat / pair.temperature().value(),
        5e-9,
    );
    near(
        step.source_entropy_change_j_k + step.entropy_seam_correction_change_j_k,
        step.thermal_entropy_change_j_k,
        1e-12,
    );
    let entropy_inverse = chart
        .state_at_specific_entropy(mixture.specific_entropy_j_kg_k(), entropy_solve(), || false)
        .unwrap();
    near(
        entropy_inverse.state().vapor_quality().unwrap(),
        quality,
        1e-14,
    );
    let cooled = chart
        .advance_heat(step.state, mass, -heat, 2e-5, || false)
        .unwrap();
    assert_eq!(cooled.state, liquid);
    near(step.boundary_work_j + cooled.boundary_work_j, 0.0, 1e-10);
    near(
        step.thermal_entropy_change_j_k + cooled.thermal_entropy_change_j_k,
        0.0,
        1e-10,
    );
}

#[test]
fn independent_cp_and_cp_over_t_integrals_follow_supercritical_enthalpy_and_entropy() {
    let chart = chart(25e6, 630.0, 670.0);
    let cold = at_t(&chart, 630.0);
    let hot = at_t(&chart, 670.0);
    // The pseudocritical Cp peak is resolved by independent composite Simpson
    // integration. This oracle never differences the caloric source values.
    let panels = 1024;
    let dt = 40.0 / f64::from(panels);
    let mut h_integral = 0.0;
    let mut s_integral = 0.0;
    for i in 0..=panels {
        let t = 630.0 + f64::from(i) * dt;
        let source = at_t(&chart, t).single_phase_source().unwrap();
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        h_integral += weight * source.cp().value();
        s_integral += weight * source.cp().value() / t;
    }
    h_integral *= dt / 3.0;
    s_integral *= dt / 3.0;
    near(
        hot.specific_enthalpy_j_kg() - cold.specific_enthalpy_j_kg(),
        h_integral,
        0.01,
    );
    near(
        hot.specific_entropy_j_kg_k() - cold.specific_entropy_j_kg_k(),
        s_integral,
        2e-5,
    );
    assert!(hot.bulk_density_kg_m3() < cold.bulk_density_kg_m3());
}

#[test]
fn heating_subdivision_small_heat_and_closed_cycle_preserve_h_and_pressure_work() {
    let chart = chart(25e6, 629.0, 671.0);
    let cold = at_t(&chart, 630.0);
    let hot = at_t(&chart, 670.0);
    let mass = 0.5;
    let heat = mass * (hot.specific_enthalpy_j_kg() - cold.specific_enthalpy_j_kg());
    let step = chart
        .advance_heat(cold, mass, heat, 2e-5, || false)
        .unwrap();
    near(step.state.temperature_k(), 670.0, 2e-8);
    near(step.enthalpy_change_j, heat, 2e-5);
    near(
        step.internal_energy_change_j + step.boundary_work_j,
        heat,
        2e-5,
    );
    assert!(step.boundary_work_j > 10_000.0);
    assert_eq!(step.phase_masses_kg, None);
    assert_eq!(step.initial_phase_masses_kg, Some([mass, 0.0]));
    let mut state = cold;
    let mut supplied = 0.0;
    let mut work = 0.0;
    let mut thermal_entropy = 0.0;
    for _ in 0..20 {
        let piece = chart
            .advance_heat(state, mass, heat / 20.0, 2e-5, || false)
            .unwrap();
        state = piece.state;
        supplied += piece.supplied_heat_j;
        work += piece.boundary_work_j;
        thermal_entropy += piece.thermal_entropy_change_j_k;
    }
    near(state.temperature_k(), step.state.temperature_k(), 2e-7);
    near(supplied, heat, 1e-8);
    near(work, step.boundary_work_j, 2e-4);
    near(thermal_entropy, step.thermal_entropy_change_j_k, 1e-7);
    let cool = chart
        .advance_heat(state, mass, -heat, 2e-5, || false)
        .unwrap();
    near(cool.state.temperature_k(), cold.temperature_k(), 2e-7);
    near(work + cool.boundary_work_j, 0.0, 2e-4);
    near(thermal_entropy + cool.thermal_entropy_change_j_k, 0.0, 1e-7);
    let small = chart
        .advance_heat(cold, mass, 0.001, 2e-6, || false)
        .unwrap();
    assert!(small.state.temperature_k() > cold.temperature_k());
    assert!(small.enthalpy_change_j > 0.0);
    near(small.enthalpy_change_j, 0.001, 2e-6);
    let zero = chart.advance_heat(cold, mass, 0.0, 0.0, || false).unwrap();
    assert_eq!(zero.state, cold);
    assert_eq!(zero.enthalpy_change_j, 0.0);
    assert_eq!(zero.boundary_work_j, 0.0);
    assert_eq!(zero.inversion_iterations, 0);
}

#[test]
fn pressure_closure_must_fit_the_physical_joule_budget_even_at_exact_enthalpy_knots() {
    let mut c = controls(630.0, 670.0);
    c.pressure.pressure_tolerance_pa = 1e4;
    c.pressure.density_tolerance_kg_m3 = 1e-3;
    let chart = DenseIsobaricWater::new(Pressure::new(25e6), c, || false).unwrap();
    let cold = at_t(&chart, 630.0);
    let hot = at_t(&chart, 670.0);
    let mass = 0.5;
    let defect_change =
        mass * (hot.pressure_closure_defect_j_kg() - cold.pressure_closure_defect_j_kg());
    assert!(defect_change.abs() > 1e-4);
    let heat = mass * (hot.specific_enthalpy_j_kg() - cold.specific_enthalpy_j_kg());
    let budget = 0.25 * defect_change.abs();
    assert!(
        matches!(chart.advance_heat(cold, mass, heat, budget, || false),
        Err(DenseIsobaricError::EnergyBalance { maximum_residual_j, .. }) if maximum_residual_j > budget)
    );
    let accepted = chart
        .advance_heat(cold, mass, heat, 2.0 * defect_change.abs(), || false)
        .unwrap();
    near(accepted.represented_enthalpy_residual_j, 0.0, 0.0);
    near(accepted.represented_energy_residual_j, -defect_change, 1e-8);
    assert!(accepted.energy_residual_j.abs() > budget);
}

#[test]
fn region_ambiguity_budget_representation_and_cancellation_refusals_are_atomic() {
    let chart = chart(20e6, 625.0, 649.0);
    let state = at_t(&chart, 630.0);
    let pair = chart.saturation().unwrap();
    assert_eq!(
        chart.state_at_temperature(pair.temperature(), || false),
        Err(DenseIsobaricError::AmbiguousSaturation)
    );
    assert_eq!(
        chart.state_at_specific_enthalpy(chart.enthalpy_bounds_j_kg()[0].next_down(), || false),
        Err(DenseIsobaricError::OutsideDomain)
    );
    assert!(matches!(
        chart.state_at_quality(-0.01, || false),
        Err(DenseIsobaricError::Invalid(_))
    ));
    assert_eq!(
        chart.state_at_quality(0.5, || true),
        Err(DenseIsobaricError::Cancelled)
    );
    assert_eq!(
        chart.advance_heat(state, 0.5, 1.0, 1e-4, || true),
        Err(DenseIsobaricError::Cancelled)
    );
    assert!(matches!(
        chart.advance_heat(state, 0.5, f64::MIN_POSITIVE, 1.0, || false),
        Err(DenseIsobaricError::Unrepresentable(_))
    ));
    assert!(matches!(
        chart.advance_heat(state, 0.0, 1.0, 1.0, || false),
        Err(DenseIsobaricError::Invalid(_))
    ));
    let another = self::chart(25e6, 630.0, 670.0);
    assert_eq!(
        another.advance_heat(state, 0.5, 1.0, 1e-4, || false),
        Err(DenseIsobaricError::ChartMismatch)
    );
    assert_eq!(
        another.state_at_quality(0.5, || false),
        Err(DenseIsobaricError::OutsideDomain)
    );
    assert_eq!(
        DenseIsobaricWater::new(Pressure::new(20e6), controls(625.0, 660.0), || false),
        Err(DenseIsobaricError::OutsideDomain)
    );
    assert_eq!(
        DenseIsobaricWater::new(
            Pressure::new(CRITICAL_PRESSURE_PA),
            controls(630.0, 650.0),
            || false
        ),
        Err(DenseIsobaricError::OutsideDomain)
    );
    let mut seam = controls(625.0, 649.0);
    seam.maximum_entropy_seam_j_kg_k = 0.0;
    assert!(matches!(
        DenseIsobaricWater::new(Pressure::new(20e6), seam, || false),
        Err(DenseIsobaricError::EntropySeam { .. })
    ));
    let mut short = controls(630.0, 670.0);
    short.maximum_iterations = 1;
    let short = DenseIsobaricWater::new(Pressure::new(25e6), short, || false).unwrap();
    let target = at_t(&another, 654.0).specific_enthalpy_j_kg();
    assert!(matches!(
        short.state_at_specific_enthalpy(target, || false),
        Err(DenseIsobaricError::NotConverged { iterations: 1, .. })
    ));
    let mut polls = 0;
    let before = state;
    let result = chart.advance_heat(state, 0.5, 0.0, 0.0, || {
        polls += 1;
        polls >= 2
    });
    assert_eq!(result, Err(DenseIsobaricError::Cancelled));
    assert_eq!(state, before);
}
