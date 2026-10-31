//! Complete intervals use actual source-owned gas steps and accepted work only.
use super::*;
use fs_flux::capillary::material::sealed::advance::{
    SealedCapillaryAdvanceControls, SealedCapillaryAdvanceError,
};

fn adaptive() -> SealedCapillaryAdvanceControls {
    SealedCapillaryAdvanceControls {
        initial_step_s: 8.0,
        minimum_step_s: 1e-7,
        maximum_step_s: 8.0,
        maximum_attempts: 4000,
        relative_tolerance: 1e-7,
        length_absolute_tolerance_m: 1e-9,
        temperature_absolute_tolerance_k: 1e-6,
        pressure_absolute_tolerance_pa: 1e-3,
    }
}

#[test]
fn adaptive_interval_keeps_rejected_heat_out_of_the_conservative_account() {
    let (liquid, wetting) = fixture(1000.0, 1.0, 0.03, PI / 3.0);
    let g = geometry();
    let original_gas = gas(volume(g), 3.5);
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .seal(original_gas.clone(), g)
        .unwrap();
    let input = SealedCapillaryInput {
        inlet_pressure_pa: 125000.0,
        duration_s: 8.0,
        gas_heat_power_w: 1e-7,
        opposing_gravity_m_s2: 9.81,
    };
    let result = owner
        .advance_candidate(input, controls(), adaptive(), || false)
        .unwrap();
    let s = result.step();
    assert_eq!(s.elapsed_time_s, input.duration_s);
    assert!(s.rejected_intervals > 0);
    assert!(s.accepted_intervals > 1);
    assert_eq!(s.attempts, s.accepted_intervals + s.rejected_intervals);
    assert!(s.step_evaluations <= 3 * s.attempts);
    assert!(s.maximum_normalized_local_error <= 1.0);
    near(
        s.gas_supplied_heat_j,
        input.gas_heat_power_w * input.duration_s,
        1e-19,
    );
    let area = PI * g.radius_m.powi(2);
    near(
        s.transferred_liquid_volume_m3,
        area * (s.wet_length_m - g.wet_length_m),
        1e-20,
    );
    near(
        s.transferred_liquid_mass_kg,
        1000.0 * s.transferred_liquid_volume_m3,
        1e-16,
    );
    near(
        original_gas.volume_m3() - result.endpoint().gas().volume_m3(),
        s.transferred_liquid_volume_m3,
        1e-20,
    );
    assert_eq!(result.endpoint().gas().mass_kg(), original_gas.mass_kg());
    assert_eq!(result.endpoint().gas().model(), original_gas.model());
    assert_eq!(owner.gas(), &original_gas);
    assert_eq!(owner.geometry(), g);
    assert!(s.liquid_viscous_dissipation_j > 0.0);
    assert_eq!(s.contact_line_viscous_dissipation_j, 0.0);
    assert_eq!(s.pinning_dissipation_j, 0.0);
    let energy_change = original_gas
        .model()
        .thermodynamics()
        .specific_internal_energy_change(
            fs_qty::Temperature::new(original_gas.state().state().temperature),
            fs_qty::Temperature::new(result.endpoint().gas().state().state().temperature),
        )
        .unwrap()
        .value()
        * original_gas.mass_kg();
    near(s.gas_stored_energy_change_j, energy_change, 2e-14);
    near(
        s.gas_stored_energy_change_j + s.gas_work_by_gas_j,
        s.gas_supplied_heat_j,
        2e-14,
    );
    near(
        s.inlet_pressure_work_j + s.gas_supplied_heat_j
            - s.gas_stored_energy_change_j
            - s.wall_surface_energy_change_j
            - s.gravitational_energy_change_j
            - s.liquid_viscous_dissipation_j,
        s.energy_balance_residual_j,
        1e-18,
    );
    assert!(
        s.energy_balance_residual_j.abs()
            <= controls().energy_absolute_tolerance_j
                + controls().energy_relative_tolerance * (s.inlet_pressure_work_j.abs() + 1e-5)
    );
}

#[test]
fn adaptive_tolerance_refines_against_independent_adiabatic_travel_time() {
    let (liquid, wetting) = fixture(1000.0, 1.0, 0.03, PI / 3.0);
    let g = geometry();
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .seal(gas(volume(g), 3.5), g)
        .unwrap();
    let target = 0.025;
    let inlet = 125000.0;
    let area = PI * g.radius_m.powi(2);
    let rate = g.radius_m.powi(2) / 8.0;
    let integrand = |x: f64| {
        let v = volume(g) - area * (x - g.wet_length_m);
        let pressure = 100000.0 * (volume(g) / v).powf(1.4);
        x / (rate * (inlet + 300.0 - pressure))
    };
    // Independent continuous adiabatic law, not the production midpoint formula.
    let n = 2048;
    let h = (target - g.wet_length_m) / n as f64;
    let mut integral = integrand(g.wet_length_m) + integrand(target);
    for k in 1..n {
        integral += (if k % 2 == 0 { 2.0 } else { 4.0 }) * integrand(g.wet_length_m + k as f64 * h);
    }
    let duration = h * integral / 3.0;
    let input = SealedCapillaryInput {
        inlet_pressure_pa: inlet,
        duration_s: duration,
        ..boundary()
    };
    let mut errors = [0.0; 2];
    let mut accepted = [0; 2];
    for (index, tolerance) in [1e-6, 1e-9].into_iter().enumerate() {
        let numerics = SealedCapillaryAdvanceControls {
            initial_step_s: duration,
            maximum_step_s: duration,
            relative_tolerance: 0.0,
            length_absolute_tolerance_m: tolerance,
            temperature_absolute_tolerance_k: 10.0,
            pressure_absolute_tolerance_pa: 1e6,
            ..adaptive()
        };
        let result = owner
            .advance_candidate(input, controls(), numerics, || false)
            .unwrap();
        errors[index] = (result.step().wet_length_m - target).abs();
        accepted[index] = result.step().accepted_intervals;
    }
    assert!(errors[1] < 2e-7, "{errors:?}");
    assert!(errors[1] < errors[0] / 10.0, "{errors:?}");
    assert!(accepted[1] > accepted[0]);
}

#[test]
fn adaptive_cancel_after_staged_progress_and_budget_failure_preserve_the_owner() {
    let (liquid, wetting) = fixture(1000.0, 1.0, 0.03, PI / 3.0);
    let g = geometry();
    let original_gas = gas(volume(g), 3.5);
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .seal(original_gas.clone(), g)
        .unwrap();
    let numerics = SealedCapillaryAdvanceControls {
        initial_step_s: 0.01,
        maximum_step_s: 0.01,
        ..adaptive()
    };
    let one_interval = SealedCapillaryInput {
        duration_s: 0.01,
        ..boundary()
    };
    let mut first_polls = 0;
    owner
        .advance_candidate(one_interval, controls(), numerics, || {
            first_polls += 1;
            false
        })
        .unwrap();
    let mut polls = 0;
    let long = SealedCapillaryInput {
        duration_s: 0.1,
        ..boundary()
    };
    let cancelled = owner.advance_candidate(long, controls(), numerics, || {
        polls += 1;
        polls > first_polls + 2
    });
    assert!(matches!(
        cancelled,
        Err(SealedCapillaryAdvanceError::Cancelled)
    ));
    assert!(polls > first_polls);
    let exhausted = owner.advance_candidate(
        long,
        controls(),
        SealedCapillaryAdvanceControls {
            maximum_attempts: 1,
            ..numerics
        },
        || false,
    );
    assert!(matches!(
        exhausted,
        Err(SealedCapillaryAdvanceError::StepBudget {
            attempts: 1,
            accepted_intervals: 1,
        })
    ));
    assert_eq!(owner.gas(), &original_gas);
    assert_eq!(owner.geometry(), g);
}

#[test]
fn adaptive_minimum_step_refuses_unresolved_error_without_partial_publication() {
    let (liquid, wetting) = fixture(1000.0, 1.0, 0.03, PI / 3.0);
    let g = geometry();
    let original_gas = gas(volume(g), 3.5);
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .seal(original_gas.clone(), g)
        .unwrap();
    let result = owner.advance_candidate(
        SealedCapillaryInput {
            inlet_pressure_pa: 125000.0,
            ..boundary()
        },
        controls(),
        SealedCapillaryAdvanceControls {
            initial_step_s: 1.0,
            minimum_step_s: 1.0,
            maximum_step_s: 1.0,
            relative_tolerance: 0.0,
            length_absolute_tolerance_m: 1e-18,
            temperature_absolute_tolerance_k: 1e-15,
            pressure_absolute_tolerance_pa: 1e-12,
            ..adaptive()
        },
        || false,
    );
    assert!(
        matches!(result, Err(SealedCapillaryAdvanceError::Tolerance { normalized_error, .. }) if normalized_error > 1.0)
    );
    assert_eq!(owner.gas(), &original_gas);
    assert_eq!(owner.geometry(), g);
}

#[test]
fn adaptive_zero_interval_is_identity_and_bad_controls_are_rejected() {
    let (liquid, wetting) = fixture(1000.0, 1.0, 0.03, PI / 3.0);
    let g = geometry();
    let original_gas = gas(volume(g), 3.5);
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .seal(original_gas.clone(), g)
        .unwrap();
    let result = owner
        .advance_candidate(
            SealedCapillaryInput {
                duration_s: 0.0,
                gas_heat_power_w: 10.0,
                ..boundary()
            },
            controls(),
            adaptive(),
            || false,
        )
        .unwrap();
    assert_eq!(result.endpoint().gas(), &original_gas);
    assert_eq!(result.endpoint().geometry(), g);
    assert_eq!(result.step().attempts, 0);
    assert_eq!(result.step().elapsed_time_s, 0.0);
    assert_eq!(result.step().gas_supplied_heat_j, 0.0);
    assert_eq!(result.step().energy_balance_residual_j, 0.0);
    for invalid in [
        SealedCapillaryAdvanceControls {
            maximum_attempts: 0,
            ..adaptive()
        },
        SealedCapillaryAdvanceControls {
            minimum_step_s: 9.0,
            ..adaptive()
        },
        SealedCapillaryAdvanceControls {
            maximum_step_s: 7.0,
            ..adaptive()
        },
    ] {
        assert!(matches!(
            owner.advance_candidate(boundary(), controls(), invalid, || false),
            Err(SealedCapillaryAdvanceError::Invalid(_))
        ));
    }
    assert!(matches!(
        owner.advance_candidate(
            SealedCapillaryInput {
                gas_heat_power_w: f64::MAX,
                duration_s: 2.0,
                ..boundary()
            },
            controls(),
            adaptive(),
            || false,
        ),
        Err(SealedCapillaryAdvanceError::Invalid(_)),
    ));
}

#[test]
fn adaptive_final_tail_below_minimum_preserves_the_requested_time_and_heat() {
    let (liquid, wetting) = fixture(1000.0, 1.0, 0.03, PI / 3.0);
    let g = geometry();
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .seal(gas(volume(g), 3.5), g)
        .unwrap();
    for (duration_s, accepted_intervals) in [(0.003, 1), (0.025, 3)] {
        let result = owner
            .advance_candidate(
                SealedCapillaryInput {
                    duration_s,
                    gas_heat_power_w: 1e-7,
                    ..boundary()
                },
                controls(),
                SealedCapillaryAdvanceControls {
                    initial_step_s: 0.01,
                    minimum_step_s: 0.01,
                    maximum_step_s: 0.01,
                    maximum_attempts: 3,
                    ..adaptive()
                },
                || false,
            )
            .unwrap();
        assert_eq!(result.step().elapsed_time_s, duration_s);
        assert_eq!(result.step().accepted_intervals, accepted_intervals);
        assert_eq!(result.step().attempts, accepted_intervals);
        near(result.step().gas_supplied_heat_j, 1e-7 * duration_s, 1e-22);
    }
    assert_eq!(owner.geometry(), g);
}

#[test]
fn adaptive_heating_resolves_depinning_and_keeps_distinct_contact_line_losses() {
    use super::super::contact_line::{limits, sourced_wetting};
    use fs_flux::capillary::ContactLineMotion;

    let young: f64 = 1.0;
    let receding: f64 = 0.8;
    let (liquid, _) = fixture(1000.0, 1.0, 0.03, young);
    let wetting = sourced_wetting(&liquid, young, 2.0, receding, 1.2);
    let g = geometry();
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .with_contact_line(limits())
        .unwrap()
        .seal(gas(volume(g), 3.5), g)
        .unwrap();
    let input = SealedCapillaryInput {
        inlet_pressure_pa: 100000.0 - 2.0 * 0.03 * young.cos() / g.radius_m,
        gas_heat_power_w: 1e-6,
        duration_s: 2.0,
        opposing_gravity_m_s2: 0.0,
    };
    let numerics = SealedCapillaryAdvanceControls {
        initial_step_s: 2.0,
        maximum_step_s: 2.0,
        relative_tolerance: 0.0,
        length_absolute_tolerance_m: 1e-10,
        ..adaptive()
    };
    let whole = owner
        .advance_candidate(input, controls(), numerics, || false)
        .unwrap();

    // During pinning the actual manufactured gas is isochoric and has gamma=1.4.
    // Its exact pressure rise reaches the receding threshold at this time.
    let threshold = input.inlet_pressure_pa + 2.0 * 0.03 * receding.cos() / g.radius_m;
    let depinning_time = (threshold - 100000.0) * volume(g) / (0.4 * input.gas_heat_power_w);
    assert!(depinning_time > 0.0 && depinning_time < input.duration_s);
    let pinned = owner
        .step_candidate(
            SealedCapillaryInput {
                duration_s: depinning_time,
                ..input
            },
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(
        pinned.step().contact_line_motion,
        Some(ContactLineMotion::Pinned)
    );
    assert_eq!(pinned.step().transferred_liquid_volume_m3, 0.0);
    let after = pinned
        .endpoint()
        .advance_candidate(
            SealedCapillaryInput {
                duration_s: input.duration_s - depinning_time,
                ..input
            },
            controls(),
            numerics,
            || false,
        )
        .unwrap();
    near(whole.step().wet_length_m, after.step().wet_length_m, 3e-9);
    near(
        whole.endpoint().gas().state().state().pressure,
        after.endpoint().gas().state().state().pressure,
        0.02,
    );
    assert!(whole.step().wet_length_m < g.wet_length_m);
    assert!(whole.step().rejected_intervals > 0);
    assert!(whole.step().contact_line_viscous_dissipation_j > 0.0);
    assert!(whole.step().pinning_dissipation_j > 0.0);
    near(
        whole.step().gas_supplied_heat_j,
        input.gas_heat_power_w * input.duration_s,
        1e-18,
    );
    let step = whole.step();
    near(
        step.inlet_pressure_work_j + step.gas_supplied_heat_j
            - step.gas_stored_energy_change_j
            - step.wall_surface_energy_change_j
            - step.gravitational_energy_change_j
            - step.liquid_viscous_dissipation_j
            - step.contact_line_viscous_dissipation_j
            - step.pinning_dissipation_j,
        step.energy_balance_residual_j,
        1e-18,
    );
    assert_eq!(owner.geometry(), g);
}

#[test]
fn adaptive_boundary_refusal_never_claims_a_clipped_full_interval() {
    let (liquid, wetting) = fixture(1000.0, 1.0, 0.03, PI / 3.0);
    let g = SealedCapillaryGeometry {
        wet_length_m: 0.1,
        dead_volume_m3: 1e-9,
        ..geometry()
    };
    let original_gas = gas(volume(g), 3.5);
    let owner = MaterialCapillaryTube::new(&liquid, &wetting)
        .unwrap()
        .seal(original_gas.clone(), g)
        .unwrap();
    let result = owner.advance_candidate(
        SealedCapillaryInput {
            inlet_pressure_pa: 200000.0,
            ..boundary()
        },
        controls(),
        SealedCapillaryAdvanceControls {
            initial_step_s: 0.1,
            minimum_step_s: 0.1,
            maximum_step_s: 0.1,
            ..adaptive()
        },
        || false,
    );
    assert!(matches!(
        result,
        Err(SealedCapillaryAdvanceError::Step {
            error: SealedCapillaryError::StepTooLarge(_),
            ..
        })
    ));
    assert_eq!(owner.gas(), &original_gas);
    assert_eq!(owner.geometry(), g);
}
