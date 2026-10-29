use super::*;

fn controls(t_min: f64, t_max: f64) -> IsochoricWaterControls {
    IsochoricWaterControls {
        minimum_temperature_k: t_min,
        maximum_temperature_k: t_max,
        minimum_pressure_pa: 100.0,
        maximum_pressure_pa: 100e6,
        volume_tolerance_m3_kg: 1e-14,
        internal_energy_tolerance_j_kg: 1e-6,
        maximum_pressure_iterations: 96,
        maximum_temperature_iterations: 96,
    }
}
fn chart(volume: f64, t_min: f64, t_max: f64) -> IsochoricWater {
    IsochoricWater::new(
        MassSpecificVolume::new(volume),
        controls(t_min, t_max),
        || false,
    )
    .unwrap()
}
fn at_t(law: &IsochoricWater, t: f64) -> IsochoricWaterState {
    law.state_at_temperature(Temperature::new(t), || false)
        .unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a} != {b}, tolerance={tol}");
}

#[test]
fn original_liquid_and_vapor_round_trip_from_volume_and_internal_energy() {
    for (region, p, t) in [
        (If97Region::Region1, 10e6, 300.0),
        (If97Region::Region1, 30e6, 330.0),
        (If97Region::Region2, 1e5, 500.0),
        (If97Region::Region2, 3e6, 600.0),
    ] {
        let source =
            If97State::at_pressure_temperature(Pressure::new(p), Temperature::new(t), region)
                .unwrap();
        let volume = source.specific_volume().value();
        let mut cfg = controls(t - 1.0, t + 1.0);
        cfg.volume_tolerance_m3_kg = volume * 1e-14;
        let law = IsochoricWater::new(MassSpecificVolume::new(volume), cfg, || false).unwrap();
        let state = at_t(&law, t);
        near(
            state.specific_volume_m3_kg(),
            volume,
            cfg.volume_tolerance_m3_kg,
        );
        near(state.pressure_pa(), p, 1e-4);
        near(
            state.specific_internal_energy_j_kg(),
            source.specific_internal_energy().value(),
            2e-7,
        );
        near(
            state.specific_enthalpy_j_kg(),
            source.specific_enthalpy().value(),
            2e-7,
        );
        near(
            state.specific_entropy_j_kg_k(),
            source.specific_entropy().value(),
            1e-9,
        );
        let inverse = law
            .state_at_specific_internal_energy(source.specific_internal_energy().value(), || false)
            .unwrap();
        near(inverse.state().temperature_k(), t, 1e-8);
        assert!(
            inverse.internal_energy_residual_j_kg().abs() <= cfg.internal_energy_tolerance_j_kg
        );
        assert_eq!(state.vapor_quality(), None);
        assert_eq!(state.law(), &law);
    }
}

#[test]
fn original_saturation_inventory_and_isochoric_heat_capacity_match_independent_derivatives() {
    let law = chart(0.1, 300.0, 623.15);
    for t in [320.0, 380.0, 450.0, 550.0, 600.0] {
        let state = at_t(&law, t);
        let delta = 0.002;
        let left = at_t(&law, t - delta);
        let right = at_t(&law, t + delta);
        let numerical_cv = (right.specific_internal_energy_j_kg()
            - left.specific_internal_energy_j_kg())
            / (2.0 * delta);
        near(state.heat_capacity_at_volume_j_kg_k(), numerical_cv, 0.003);
        if let Some(q) = state.vapor_quality() {
            let pair = If97SaturationPair::at_temperature(Temperature::new(t)).unwrap();
            let (l, v) = (pair.liquid(), pair.vapor());
            let expected = (0.1 - l.specific_volume().value())
                / (v.specific_volume().value() - l.specific_volume().value());
            near(q, expected, 1e-15);
            near(
                state.specific_internal_energy_j_kg(),
                (1.0 - q) * l.specific_internal_energy().value()
                    + q * v.specific_internal_energy().value(),
                1e-9,
            );
            near(
                state.specific_entropy_j_kg_k(),
                (1.0 - q) * l.specific_entropy().value() + q * v.specific_entropy().value(),
                1e-11,
            );
            assert_eq!(
                state.saturation_gibbs_gap_j_kg(),
                Some(pair.gibbs_energy_gap().value())
            );
            near(state.pressure_pa(), pair.pressure().value(), 0.0);
        } else {
            assert_eq!(state.phase(), WaterPhase::SuperheatedVapor);
            assert_eq!(state.saturation_gibbs_gap_j_kg(), None);
        }
    }
    // Exact authored saturated volumes retain complete phase depletion.
    let pair = If97SaturationPair::at_temperature(Temperature::new(400.0)).unwrap();
    let vapor = chart(pair.vapor().specific_volume().value(), 399.0, 401.0);
    assert_eq!(at_t(&vapor, 400.0).vapor_quality(), Some(1.0));
    let mut cfg = controls(399.99, 400.01);
    cfg.volume_tolerance_m3_kg = 1e-17;
    let liquid = IsochoricWater::new(pair.liquid().specific_volume(), cfg, || false).unwrap();
    assert_eq!(at_t(&liquid, 400.0).vapor_quality(), Some(0.0));
}

#[test]
fn sealed_vessel_heating_evaporates_and_cooling_condenses_with_zero_boundary_work() {
    let law = chart(0.1, 300.0, 623.15);
    let mass = 0.05;
    let initial = at_t(&law, 300.0);
    let hot = at_t(&law, 600.0);
    assert_eq!(initial.phase(), WaterPhase::LiquidVapor);
    assert_eq!(hot.phase(), WaterPhase::SuperheatedVapor);
    assert!(hot.pressure_pa() > 100.0 * initial.pressure_pa());
    let heat =
        mass * (hot.specific_internal_energy_j_kg() - initial.specific_internal_energy_j_kg());
    let heated = law
        .advance_heat(initial, mass, heat, 1e-7, || false)
        .unwrap();
    near(heated.state.temperature_k(), 600.0, 1e-8);
    assert_eq!(heated.boundary_work_j, 0.0);
    near(heated.fixed_volume_m3, 0.005, 1e-18);
    near(heated.internal_energy_change_j, heat, 1e-7);
    near(
        heated.enthalpy_change_j - heated.pressure_enthalpy_change_j,
        heat,
        1e-7,
    );
    // Using constant-pressure heat would overcount the pressure rise's H-U change.
    assert!(heated.enthalpy_change_j - heat > 10_000.0);
    assert_eq!(heated.phase_masses_kg, [0.0, mass]);
    let cooled = law
        .advance_heat(
            heated.state,
            mass,
            mass * (initial.specific_internal_energy_j_kg()
                - heated.state.specific_internal_energy_j_kg()),
            1e-7,
            || false,
        )
        .unwrap();
    near(cooled.state.temperature_k(), initial.temperature_k(), 1e-9);
    near(cooled.state.pressure_pa(), initial.pressure_pa(), 1e-6);
    near(
        cooled.state.vapor_mass_fraction(),
        initial.vapor_mass_fraction(),
        1e-12,
    );
    near(
        heated.internal_energy_change_j + cooled.internal_energy_change_j,
        0.0,
        1e-7,
    );
    near(
        heated.enthalpy_change_j + cooled.enthalpy_change_j,
        0.0,
        1e-7,
    );
    near(
        heated.source_entropy_change_j_k + cooled.source_entropy_change_j_k,
        0.0,
        1e-9,
    );
    for step in [heated, cooled] {
        near(step.phase_masses_kg.iter().sum(), mass, 1e-17);
        near(step.internal_energy_residual_j, 0.0, 1e-7);
        near(step.enthalpy_residual_j, 0.0, 1e-7);
        assert!(
            step.endpoint_volume_residuals_m3
                .iter()
                .all(|r| r.abs() <= mass * law.controls.volume_tolerance_m3_kg)
        );
    }
}

#[test]
fn compressed_liquid_pressurization_uses_cv_and_preserves_water_density_anomaly() {
    let source = If97State::at_pressure_temperature(
        Pressure::new(10e6),
        Temperature::new(300.0),
        If97Region::Region1,
    )
    .unwrap();
    let mut cfg = controls(290.0, 330.0);
    cfg.volume_tolerance_m3_kg = 1e-17;
    let law = IsochoricWater::new(source.specific_volume(), cfg, || false).unwrap();
    let cold = at_t(&law, 300.0);
    let hot = at_t(&law, 325.0);
    assert_eq!(cold.phase(), WaterPhase::CompressedLiquid);
    assert_eq!(hot.phase(), WaterPhase::CompressedLiquid);
    assert!(hot.pressure_pa() > cold.pressure_pa() + 10e6);
    let numerical_cv = (at_t(&law, 300.001).specific_internal_energy_j_kg()
        - at_t(&law, 299.999).specific_internal_energy_j_kg())
        / 0.002;
    near(numerical_cv, source.cv().value(), 0.003);
    let step = law
        .advance_heat(
            cold,
            1.0,
            hot.specific_internal_energy_j_kg() - cold.specific_internal_energy_j_kg(),
            2e-6,
            || false,
        )
        .unwrap();
    near(step.state.temperature_k(), hot.temperature_k(), 1e-8);
    assert!(step.pressure_enthalpy_change_j > 10_000.0);
    near(
        step.enthalpy_change_j - step.internal_energy_change_j,
        step.pressure_enthalpy_change_j,
        2e-6,
    );
    // Around cold water's density anomaly, heating at fixed volume can LOWER p.
    let anomaly = If97State::at_pressure_temperature(
        Pressure::new(1e6),
        Temperature::new(274.0),
        If97Region::Region1,
    )
    .unwrap();
    let mut cfg = controls(273.16, 275.0);
    cfg.volume_tolerance_m3_kg = 1e-17;
    let cold_water = IsochoricWater::new(anomaly.specific_volume(), cfg, || false).unwrap();
    let a = at_t(&cold_water, 273.8);
    let b = at_t(&cold_water, 274.2);
    assert!(b.specific_internal_energy_j_kg() > a.specific_internal_energy_j_kg());
    assert!(b.pressure_pa() < a.pressure_pa());
}

#[test]
fn heat_subdivision_keeps_original_energy_and_phase_inventory() {
    let law = chart(0.1, 300.0, 623.15);
    let initial = at_t(&law, 310.0);
    let final_state = at_t(&law, 590.0);
    let mass = 0.2;
    let total_heat = mass
        * (final_state.specific_internal_energy_j_kg() - initial.specific_internal_energy_j_kg());
    let direct = law
        .advance_heat(initial, mass, total_heat, 1e-6, || false)
        .unwrap();
    let mut current = initial;
    let mut actual_heat = 0.0;
    for i in 1..=20 {
        let target_u = initial.specific_internal_energy_j_kg()
            + (final_state.specific_internal_energy_j_kg()
                - initial.specific_internal_energy_j_kg())
                * f64::from(i)
                / 20.0;
        let q = mass * (target_u - current.specific_internal_energy_j_kg());
        let step = law.advance_heat(current, mass, q, 1e-6, || false).unwrap();
        current = step.state;
        actual_heat += step.internal_energy_change_j;
    }
    near(current.temperature_k(), direct.state.temperature_k(), 1e-8);
    near(actual_heat, direct.internal_energy_change_j, 1e-6);
    let tiny = law
        .advance_heat(final_state, mass, 1e-4, 1e-7, || false)
        .unwrap();
    assert!(tiny.state.temperature_k() > final_state.temperature_k());
    near(tiny.internal_energy_change_j, 1e-4, 1e-7);
    assert_eq!(
        law.advance_heat(initial, mass, 0.0, 0.0, || false)
            .unwrap()
            .state,
        initial
    );
    assert!(matches!(
        law.advance_heat(initial, mass, 1e-30, 1.0, || false),
        Err(IsochoricWaterError::Unrepresentable(_))
    ));
}

#[test]
fn domain_work_energy_and_cancellation_refusals_publish_no_partial_state() {
    let law = chart(0.1, 300.0, 623.15);
    let initial = at_t(&law, 350.0);
    let snapshot = initial;
    assert_eq!(
        law.state_at_temperature(Temperature::new(650.0), || false),
        Err(IsochoricWaterError::OutsideDomain)
    );
    assert!(
        IsochoricWater::new(MassSpecificVolume::new(0.0), controls(300.0, 500.0), || {
            false
        })
        .is_err()
    );
    let mut cfg = controls(300.0, 500.0);
    cfg.maximum_temperature_k = 624.0;
    assert!(matches!(
        IsochoricWater::new(MassSpecificVolume::new(0.1), cfg, || false),
        Err(IsochoricWaterError::Invalid(_))
    ));
    assert!(matches!(
        IsochoricWater::new(
            MassSpecificVolume::new(1e-5),
            controls(300.0, 500.0),
            || false
        ),
        Err(IsochoricWaterError::VolumeOutsideDomain { .. })
    ));
    let different = chart(0.2, 300.0, 600.0);
    assert_eq!(
        different.advance_heat(initial, 0.1, 100.0, 1e-6, || false),
        Err(IsochoricWaterError::ChartMismatch)
    );
    let mut cfg = controls(300.0, 623.15);
    cfg.maximum_temperature_iterations = 1;
    let limited = IsochoricWater::new(MassSpecificVolume::new(0.1), cfg, || false).unwrap();
    assert!(matches!(
        limited.state_at_specific_internal_energy(
            at_t(&law, 450.0).specific_internal_energy_j_kg(),
            || false
        ),
        Err(IsochoricWaterError::EnergyNotConverged { .. })
    ));
    let mut calls = 0;
    let result = law.advance_heat(initial, 0.1, 10_000.0, 1e-6, || {
        calls += 1;
        calls == 5
    });
    assert_eq!(result, Err(IsochoricWaterError::Cancelled));
    assert_eq!(calls, 5);
    assert_eq!(initial, snapshot);
    let high = at_t(&law, 600.0);
    assert!(law.advance_heat(high, 0.1, 1e20, 1e-6, || false).is_err());
    // Coarse pressure-volume admission cannot be hidden by an energy budget.
    let mut coarse = controls(500.0, 600.0);
    coarse.volume_tolerance_m3_kg = 1e-3;
    let coarse = IsochoricWater::new(MassSpecificVolume::new(0.1), coarse, || false).unwrap();
    let a = at_t(&coarse, 510.0);
    let b = at_t(&coarse, 590.0);
    assert!(
        coarse
            .advance_heat(
                a,
                1.0,
                b.specific_internal_energy_j_kg() - a.specific_internal_energy_j_kg(),
                1e-5,
                || false
            )
            .is_err()
    );
}
