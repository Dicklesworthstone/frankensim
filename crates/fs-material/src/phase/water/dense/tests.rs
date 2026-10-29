use super::*;

fn controls() -> DenseWaterControls {
    DenseWaterControls {
        minimum_temperature_k: 650.0,
        maximum_temperature_k: 750.0,
        internal_energy_tolerance_j_kg: 1e-6,
        maximum_iterations: 96,
        interval: DenseIsochoricIntervalSolve {
            maximum_intervals: 65_535,
        },
    }
}
fn chart() -> DenseIsochore {
    DenseIsochore::new(Density::new(500.0), controls(), || false).unwrap()
}
fn at_t(chart: &DenseIsochore, t: f64) -> DenseWaterState {
    chart
        .state_at_temperature(Temperature::new(t), || false)
        .unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:.15e} != {expected:.15e}, tolerance {tolerance:.8e}"
    );
}

#[test]
fn original_source_round_trips_and_proven_cv_survive_preparation() {
    let chart = chart();
    let bounds = chart.source_interval();
    assert!(bounds.minimum_heat_capacity_j_kg_k() > 0.0);
    for t in [650.0, 650.001, 673.15, 700.0, 749.999, 750.0] {
        let original =
            If97DenseState::at_density_temperature(Density::new(500.0), Temperature::new(t))
                .unwrap();
        let state = at_t(&chart, t);
        assert_eq!(state.source(), original);
        let inverse = chart
            .state_at_specific_internal_energy(state.specific_internal_energy_j_kg(), || false)
            .unwrap();
        near(inverse.state().temperature_k(), t, 1e-9);
        assert!(
            inverse.internal_energy_residual_j_kg().abs()
                <= controls().internal_energy_tolerance_j_kg
        );
        assert!(chart.maximum_temperature_slope() >= 1.0 / original.cv().value());
        assert!(original.cv().value() >= bounds.minimum_heat_capacity_j_kg_k());
        assert!(original.cv().value() <= bounds.maximum_heat_capacity_j_kg_k());
        assert_eq!(state.density_kg_m3(), 500.0);
        assert_eq!(state.specific_volume_m3_kg(), 0.002);
    }
    // Published Table 33 anchors, retained through the material owner.
    near(at_t(&chart, 650.0).pressure_pa(), 25.5837018e6, 0.05);
    near(at_t(&chart, 750.0).pressure_pa(), 78.3095639e6, 0.05);
}

#[test]
fn dense_rigid_heating_uses_internal_energy_and_keeps_pressure_enthalpy_separate() {
    let chart = chart();
    let initial = at_t(&chart, 650.0);
    let hot = at_t(&chart, 750.0);
    let mass = 0.5;
    let heat =
        mass * (hot.specific_internal_energy_j_kg() - initial.specific_internal_energy_j_kg());
    let step = chart
        .advance_heat(initial, mass, heat, 1e-5, || false)
        .unwrap();
    assert_eq!(step.state, hot);
    assert_eq!(step.volume_m3, 0.001);
    assert_eq!(step.boundary_work_j, 0.0);
    assert!(step.pressure_enthalpy_change_j > 50_000.0);
    near(step.internal_energy_change_j, heat, 1e-5);
    near(
        step.enthalpy_change_j - step.pressure_enthalpy_change_j,
        heat,
        1e-5,
    );
    near(
        step.pressure_enthalpy_change_j,
        step.volume_m3 * (hot.pressure_pa() - initial.pressure_pa()),
        1e-5,
    );
    for residual in [
        step.energy_residual_j,
        step.enthalpy_residual_j,
        step.represented_energy_residual_j,
        step.represented_enthalpy_residual_j,
    ] {
        near(residual, 0.0, 1e-5);
    }
    assert!(step.entropy_change_j_k > 0.0);
    // Exact authored endpoints still retain source evaluation roundoff. A zero
    // physical budget must not silently promote that numerical discrepancy.
    assert!(matches!(
        chart.advance_heat(initial, mass, heat, 0.0, || false),
        Err(DenseWaterError::EnergyBalance { maximum_residual_j, .. })
            if maximum_residual_j > 0.0
    ));
    // A constant-Cp or Q=Delta H shortcut would charge a large spurious heat.
    assert!((step.enthalpy_change_j - heat).abs() > 50_000.0);
    let cold = chart
        .advance_heat(hot, mass, -heat, 1e-5, || false)
        .unwrap();
    assert_eq!(cold.state, initial);
    near(
        step.entropy_change_j_k + cold.entropy_change_j_k,
        0.0,
        1e-10,
    );
    near(
        step.internal_energy_change_j + cold.internal_energy_change_j,
        0.0,
        1e-8,
    );
}

#[test]
fn original_cv_integral_and_finite_heat_subdivision_agree() {
    let chart = chart();
    let initial = at_t(&chart, 650.0);
    let hot = at_t(&chart, 750.0);
    let mass = 0.5;
    let heat =
        mass * (hot.specific_internal_energy_j_kg() - initial.specific_internal_energy_j_kg());
    // Independent Simpson integration of analytic source Cv; no energy differences.
    let panels = 256;
    let width = 100.0 / f64::from(panels);
    let mut cv = 0.0;
    let mut cv_over_t = 0.0;
    for i in 0..=panels {
        let t = 650.0 + f64::from(i) * width;
        let source = at_t(&chart, t).source();
        let weight = if i == 0 || i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        cv += weight * source.cv().value();
        cv_over_t += weight * source.cv().value() / t;
    }
    near(mass * width * cv / 3.0, heat, 2e-4);
    let whole = chart
        .advance_heat(initial, mass, heat, 1e-5, || false)
        .unwrap();
    near(
        mass * width * cv_over_t / 3.0,
        whole.entropy_change_j_k,
        5e-7,
    );

    let mut current = initial;
    let mut total_heat = 0.0;
    let mut total_entropy = 0.0;
    for i in 1..=20 {
        let next = at_t(&chart, 650.0 + 5.0 * f64::from(i));
        let q =
            mass * (next.specific_internal_energy_j_kg() - current.specific_internal_energy_j_kg());
        let step = chart
            .advance_heat(current, mass, q, 1e-5, || false)
            .unwrap();
        total_heat += step.internal_energy_change_j;
        total_entropy += step.entropy_change_j_k;
        current = step.state;
    }
    near(current.temperature_k(), hot.temperature_k(), 1e-9);
    near(total_heat, heat, 1e-5);
    near(total_entropy, whole.entropy_change_j_k, 1e-8);
    let small = chart
        .advance_heat(at_t(&chart, 700.0), mass, 0.001, 1e-6, || false)
        .unwrap();
    assert!(small.state.temperature_k() > 700.0);
    near(small.internal_energy_change_j, 0.001, 1e-6);
}

#[test]
fn bounds_work_cancel_and_physical_energy_failures_publish_no_endpoint() {
    let chart = chart();
    let initial = at_t(&chart, 700.0);
    assert_eq!(
        chart
            .advance_heat(initial, 0.5, 0.0, 0.0, || false)
            .unwrap()
            .state,
        initial
    );
    assert!(matches!(
        chart.state_at_temperature(Temperature::new(750.0_f64.next_up()), || false),
        Err(DenseWaterError::OutsideDomain)
    ));
    assert!(
        chart
            .advance_heat(initial, 0.0, 1.0, 1e-5, || false)
            .is_err()
    );
    assert!(
        chart
            .advance_heat(initial, 0.5, 1e-30, 1.0, || false)
            .is_err()
    );
    assert!(
        chart
            .advance_heat(initial, 0.5, 1e9, 1e-5, || false)
            .is_err()
    );
    let mut cfg = controls();
    cfg.internal_energy_tolerance_j_kg = f64::NAN;
    assert!(DenseIsochore::new(Density::new(500.0), cfg, || false).is_err());
    cfg = controls();
    cfg.maximum_iterations = 1;
    let limited = DenseIsochore::new(Density::new(500.0), cfg, || false).unwrap();
    assert!(matches!(
        limited
            .state_at_specific_internal_energy(initial.specific_internal_energy_j_kg(), || false),
        Err(DenseWaterError::NotConverged { iterations: 1, .. })
    ));
    assert!(matches!(
        limited.advance_heat(initial, 0.5, 1.0, 1e-5, || false),
        Err(DenseWaterError::ChartMismatch)
    ));
    assert!(matches!(
        DenseIsochore::new(Density::new(500.0), controls(), || true),
        Err(DenseWaterError::Cancelled)
    ));

    let mut calls = 0;
    chart
        .advance_heat(initial, 0.5, 100.0, 1e-5, || {
            calls += 1;
            false
        })
        .unwrap();
    assert!(calls > 3);
    let mut late = 0;
    assert_eq!(
        chart.advance_heat(initial, 0.5, 100.0, 1e-5, || {
            late += 1;
            late == calls
        }),
        Err(DenseWaterError::Cancelled)
    );
    assert_eq!(initial, at_t(&chart, 700.0));
}
