//! G0/G1/G3/G4 explicit two-pool transport and signed exterior heat.

use super::*;

fn law(times: [f64; 2], exterior_kev: [f64; 2]) -> FusionEnergyConfinement {
    FusionEnergyConfinement::new(
        Time::new(times[0]),
        Time::new(times[1]),
        Temperature::new(exterior_kev[0] * KELVIN_PER_KEV),
        Temperature::new(exterior_kev[1] * KELVIN_PER_KEV),
    )
    .unwrap()
}

#[test]
fn g0_confinement_admits_only_physical_parameters_and_retains_small_energy() {
    for invalid in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            FusionEnergyConfinement::new(
                Time::new(invalid),
                Time::new(1.0),
                Temperature::new(0.0),
                Temperature::new(0.0),
            )
            .is_err()
        );
        assert!(
            FusionEnergyConfinement::new(
                Time::new(1.0),
                Time::new(invalid),
                Temperature::new(0.0),
                Temperature::new(0.0),
            )
            .is_err()
        );
    }
    for invalid in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(
            FusionEnergyConfinement::new(
                Time::new(1.0),
                Time::new(1.0),
                Temperature::new(invalid),
                Temperature::new(1.0),
            )
            .is_err()
        );
        assert!(
            FusionEnergyConfinement::new(
                Time::new(1.0),
                Time::new(1.0),
                Temperature::new(1.0),
                Temperature::new(invalid),
            )
            .is_err()
        );
    }
    let weak = FusionEnergyConfinement::new(
        Time::new(1e20),
        Time::new(1e20),
        Temperature::new(0.0),
        Temperature::new(3.0),
    )
    .unwrap();
    let (_, outgoing) =
        crate::fusion::thermal::confined_energy([1.0, 2.0], [1.0; 2], 1.0, weak).unwrap();
    close(outgoing[0], 1e-20, 2e-14, 0.0);
    close(outgoing[1], -1e-20, 2e-14, 0.0);
    let deep = law([1.0, 0.04], [0.0; 2]);
    let (remaining, outgoing) =
        crate::fusion::thermal::confined_energy([1.0, 1e300], [1.0; 2], 40.0, deep).unwrap();
    close(remaining[0], (-40.0_f64).exp(), 2e-14, 0.0);
    close(remaining[1], (1e300_f64.ln() - 1000.0).exp(), 2e-13, 0.0);
    assert_eq!(outgoing, [1.0, 1e300]);
    assert!(remaining.iter().all(|u| *u > 0.0));
}

#[test]
fn g1_confinement_relaxes_to_each_exterior_with_signed_heat_and_restart() {
    let (model, initial) = source::plasma([0.0, 0.0, 5e19], 1.0, [10.0, 5.0]);
    let closure = law([0.2, 0.07], [2.0, 8.0]);
    let transport = DtFusionTransport {
        energy_confinement: Some(closure),
        ..DtFusionTransport::default()
    };
    let run = |state: &FrozenIonizedGasStateV1, duration| {
        advance_dt_fusion_with_transport(
            state,
            &model,
            heating(0.8),
            &transport,
            Time::new(duration),
            controls(),
            || false,
        )
        .unwrap()
    };
    let first = run(&initial, 0.15);
    let resumed = run(first.state(), 0.08);
    let whole = run(&initial, 0.23);
    let initial_temperatures = [initial.heavy_temperature(), initial.electron_temperature()];
    let final_temperatures = [
        whole.state().heavy_temperature(),
        whole.state().electron_temperature(),
    ];
    let capacities = [
        initial.heavy_heat_capacity_j_k(),
        initial.electron_heat_capacity_j_k(),
    ];
    let out = [
        whole.balance().heavy_transport_energy_out_j,
        whole.balance().electron_transport_energy_out_j,
    ];
    for i in 0..2 {
        let bath = closure.exterior_temperatures()[i].value();
        let expected = bath
            + (initial_temperatures[i].value() - bath)
                * (-0.23 / closure.confinement_times()[i].value()).exp();
        close(final_temperatures[i].value(), expected, 2e-13, 1e-5);
        close(
            out[i],
            capacities[i] * (initial_temperatures[i].value() - expected),
            2e-13,
            1e-7,
        );
    }
    assert!(out[0] > 0.0 && out[1] < 0.0);
    assert_eq!(whole.state().amounts(), initial.amounts());
    assert_eq!(whole.transport(), &transport);
    assert_eq!(whole.balance().reaction_extent_mol, 0.0);
    assert_eq!(whole.balance().photon_energy_escaped_j, 0.0);
    assert!(whole.balance().total_energy_residual_j.abs() < 1e-7);
    close(
        resumed.state().heavy_temperature().value(),
        final_temperatures[0].value(),
        2e-13,
        1e-5,
    );
    close(
        resumed.state().electron_temperature().value(),
        final_temperatures[1].value(),
        2e-13,
        1e-5,
    );
    close(
        first.balance().heavy_transport_energy_out_j
            + resumed.balance().heavy_transport_energy_out_j,
        out[0],
        2e-13,
        1e-7,
    );
    close(
        first.balance().electron_transport_energy_out_j
            + resumed.balance().electron_transport_energy_out_j,
        out[1],
        2e-13,
        1e-7,
    );
}

#[test]
fn g1_stiff_confinement_reaches_finite_exteriors_without_negative_energy() {
    let (model, initial) = source::plasma([0.0, 0.0, 5e19], 1.0, [10.0, 5.0]);
    let transport = DtFusionTransport {
        energy_confinement: Some(law([1e-9, 1e-8], [1.0, 2.0])),
        ..DtFusionTransport::default()
    };
    let mut c = controls();
    c.initial_step_s = 0.2;
    c.maximum_step_s = 0.2;
    c.minimum_step_s = 0.2;
    c.maximum_steps = 1;
    let result = advance_dt_fusion_with_transport(
        &initial,
        &model,
        heating(0.8),
        &transport,
        Time::new(0.2),
        c,
        || false,
    )
    .unwrap();
    close(
        result.state().heavy_temperature().value(),
        KELVIN_PER_KEV,
        1e-13,
        1e-6,
    );
    close(
        result.state().electron_temperature().value(),
        2.0 * KELVIN_PER_KEV,
        1e-13,
        1e-6,
    );
    assert_eq!(result.balance().accepted_steps, 1);
    assert_eq!(result.balance().rejected_steps, 0);
    assert!(result.state().heavy_internal_energy().value() > 0.0);
    assert!(result.state().electron_internal_energy().value() > 0.0);
    assert!(result.balance().total_energy_residual_j.abs() < 1e-7);
}

#[test]
fn g1_stiff_confinement_resolves_continuous_heating_in_each_pool() {
    let (model, initial) = source::plasma([0.0, 0.0, 5e19], 1.0, [5.0, 5.0]);
    let original = initial.clone();
    let duration = 0.1;
    let power = 10_000.0;
    let tau = 0.001;
    for pool in 0..2 {
        // The other pool is unforced and exactly relaxes to its original
        // temperature. Its much shorter time must not force tiny global steps.
        let mut times = [1e-12; 2];
        times[pool] = tau;
        let transport = DtFusionTransport {
            energy_confinement: Some(law(times, [5.0; 2])),
            ..DtFusionTransport::default()
        };
        let mut heat = heating(0.8);
        if pool == 0 {
            heat.heavy_power = Power::new(power);
        } else {
            heat.electron_power = Power::new(power);
        }
        let mut c = controls();
        c.initial_step_s = duration;
        c.maximum_step_s = duration;
        let result = advance_dt_fusion_with_transport(
            &initial,
            &model,
            heat,
            &transport,
            Time::new(duration),
            c,
            || false,
        )
        .unwrap();
        let before = [initial.heavy_internal_energy(), initial.electron_internal_energy()];
        let after = [
            result.state().heavy_internal_energy(),
            result.state().electron_internal_energy(),
        ];
        // Independent forced linear ODE, not a composition of solver maps:
        // U-Uext = P*tau*(1-exp(-t/tau)). Saturated coarse/fine split maps
        // previously agreed within 7e-9 J while losing this entire 10 J store.
        let expected = power * tau * -(-duration / tau).exp_m1();
        close(after[pool].value() - before[pool].value(), expected, 5e-5, 1e-7);
        close(after[1 - pool].value(), before[1 - pool].value(), 2e-13, 1e-7);
        let outgoing = [
            result.balance().heavy_transport_energy_out_j,
            result.balance().electron_transport_energy_out_j,
        ];
        close(outgoing[pool], power * duration - expected, 1e-6, 1e-7);
        assert!(result.balance().rejected_steps > 0);
        assert!(result.balance().total_energy_residual_j.abs() < 1e-6);

        // Insufficient resolution must refuse, never certify that false
        // equilibrium. The caller's state stays unchanged in both cases.
        c.minimum_step_s = duration;
        c.maximum_steps = 1;
        assert_eq!(
            advance_dt_fusion_with_transport(
                &initial,
                &model,
                heat,
                &transport,
                Time::new(duration),
                c,
                || false,
            ),
            Err(DtFusionError::StepBudget),
        );
    }
    assert_eq!(initial, original);
}

#[test]
fn g3_confinement_changes_burn_while_all_thermal_and_nuclear_accounts_close() {
    let (model, initial) = source::plasma([5e19, 5e19, 0.0], 1.0, [10.0, 5.0]);
    let mut transport = thermal_sources(&initial);
    let mut heat = heating(0.8);
    heat.heavy_power = Power::new(200.0);
    heat.electron_power = Power::new(500.0);
    let baseline = advance_dt_fusion_with_transport(
        &initial,
        &model,
        heat,
        &transport,
        Time::new(0.1),
        controls(),
        || false,
    )
    .unwrap();
    transport.energy_confinement = Some(law([0.3, 0.2], [0.5, 0.5]));
    let result = advance_dt_fusion_with_transport(
        &initial,
        &model,
        heat,
        &transport,
        Time::new(0.1),
        controls(),
        || false,
    )
    .unwrap();
    let b = result.balance();
    assert!(b.reaction_extent_mol < baseline.balance().reaction_extent_mol);
    assert!(b.heavy_transport_energy_out_j > 0.0);
    assert!(b.electron_transport_energy_out_j > 0.0);
    assert!(b.photon_energy_escaped_j > 0.0);
    assert!(b.electron_to_heavy_coulomb_energy_j != 0.0);
    close(
        result.state().heavy_internal_energy().value() - initial.heavy_internal_energy().value(),
        0.2 * b.alpha_energy_deposited_j - b.reacting_ion_energy_removed_j
            + b.electron_to_heavy_coulomb_energy_j
            + b.external_heavy_energy_j
            - b.heavy_transport_energy_out_j,
        3e-12,
        1e-7,
    );
    close(
        result.state().electron_internal_energy().value()
            - initial.electron_internal_energy().value(),
        0.8 * b.alpha_energy_deposited_j - b.electron_to_heavy_coulomb_energy_j
            + b.external_electron_energy_j
            - b.photon_energy_escaped_j
            - b.electron_transport_energy_out_j,
        3e-12,
        1e-7,
    );
    assert!(b.total_energy_residual_j.abs() < 1e-6);
    assert!(b.relative_baryon_residual.abs() < 1e-12);
    assert!(b.relative_charge_residual.abs() < 1e-12);
}

#[test]
fn g4_confinement_domain_refusal_and_cancellation_leave_inputs_unchanged() {
    let (model, initial) = source::plasma([0.0, 0.0, 5e19], 1.0, [10.0, 5.0]);
    let original = initial.clone();
    let transport = DtFusionTransport {
        energy_confinement: Some(law([0.005; 2], [0.0; 2])),
        ..DtFusionTransport::default()
    };
    let mut c = controls();
    c.minimum_step_s = 1e-8;
    let refused = advance_dt_fusion_with_transport(
        &initial,
        &model,
        heating(0.8),
        &transport,
        Time::new(0.1),
        c,
        || false,
    );
    assert!(matches!(refused, Err(DtFusionError::OutsideDomain)));
    for stop in [0, 75] {
        let mut polls = 0;
        let result = advance_dt_fusion_with_transport(
            &initial,
            &model,
            heating(0.8),
            &transport,
            Time::new(0.1),
            c,
            || {
                polls += 1;
                polls > stop
            },
        );
        assert_eq!(result, Err(DtFusionError::Cancelled));
    }
    assert_eq!(initial, original);
}
