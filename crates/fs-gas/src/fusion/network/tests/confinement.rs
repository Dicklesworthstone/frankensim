//! G1/G3 energy confinement with competing channels and pressure work.

use super::*;
use crate::fusion::FusionEnergyConfinement;
use crate::fusion::volume::FusionVolumeKnot;
use fs_qty::Volume;

fn confinement() -> FusionEnergyConfinement {
    FusionEnergyConfinement::new(
        Time::new(0.2),
        Time::new(0.13),
        Temperature::new(2.0 * KELVIN_PER_KEV),
        Temperature::new(KELVIN_PER_KEV),
    )
    .unwrap()
}

#[test]
fn g1_network_stiff_confinement_retains_the_forced_steady_energy() {
    let (plasma, initial) = fixture::plasma(MIXED_DENSITY, 0.025, [15.0, 4.0]);
    let duration = 0.1;
    let tau = 0.001;
    let power = 1000.0;
    let transport = FusionNetworkTransport {
        energy_confinement: Some(
            FusionEnergyConfinement::new(
                Time::new(tau),
                Time::new(1e-12),
                initial.heavy_temperature(),
                initial.electron_temperature(),
            )
            .unwrap(),
        ),
        ..FusionNetworkTransport::default()
    };
    let mut heat = heating();
    heat.heavy_power = Power::new(power);
    let mut c = controls();
    c.initial_step_s = duration;
    c.maximum_step_s = duration;
    let run = |c| {
        advance_fusion_network(
            &initial,
            &plasma,
            FusionNetworkChannels::NONE,
            heat,
            &transport,
            Time::new(duration),
            c,
            || false,
        )
    };
    let result = run(c).unwrap();
    let expected = power * tau * -(-duration / tau).exp_m1();
    close(
        result.state().heavy_internal_energy().value() - initial.heavy_internal_energy().value(),
        expected,
        1e-4,
        1e-7,
    );
    assert!(result.balance().rejected_steps > 0);
    assert_eq!(result.state().amounts(), initial.amounts());
    check_accounts(&result);
    c.minimum_step_s = duration;
    c.maximum_steps = 1;
    assert_eq!(run(c), Err(FusionNetworkError::StepBudget));
}

#[test]
fn g1_confinement_during_compression_matches_independent_linear_ode() {
    let (plasma, initial) = fixture::plasma(MIXED_DENSITY, 0.025, [15.0, 4.0]);
    let duration = 0.1;
    let final_volume = initial.volume().value() * 0.7;
    let geometry = FusionVolumeHistory::new(
        vec![
            FusionVolumeKnot {
                time: Time::new(0.0),
                volume: initial.volume(),
            },
            FusionVolumeKnot {
                time: Time::new(duration),
                volume: Volume::new(final_volume),
            },
        ],
        || false,
    )
    .unwrap();
    let law = confinement();
    let transport = FusionNetworkTransport {
        energy_confinement: Some(law),
        ..FusionNetworkTransport::default()
    };
    let initial_temperatures = [initial.heavy_temperature(), initial.electron_temperature()];
    let capacities = [
        initial.heavy_heat_capacity_j_k(),
        initial.electron_heat_capacity_j_k(),
    ];
    let mut expected_temperatures = [0.0; 2];
    let mut expected_out = [0.0; 2];
    for i in 0..2 {
        // Direct continuum equation: T' = a*T + T_ext/tau at fixed stock.
        // Work and transport are independent integrals, not solver receipts.
        let tau = law.confinement_times()[i].value();
        let exterior = law.exterior_temperatures()[i].value();
        let a = -(2.0 / 3.0) * 0.7_f64.ln() / duration - 1.0 / tau;
        let equilibrium = -exterior / (tau * a);
        let amplitude = initial_temperatures[i].value() - equilibrium;
        expected_temperatures[i] = equilibrium + amplitude * (a * duration).exp();
        let integral = equilibrium * duration + amplitude * (a * duration).exp_m1() / a;
        expected_out[i] = capacities[i] / tau * (integral - exterior * duration);
    }
    let run = |c| {
        advance_fusion_network_with_volume(
            &initial,
            &plasma,
            FusionNetworkChannels::NONE,
            heating(),
            &transport,
            &geometry,
            c,
            || false,
        )
        .unwrap()
    };
    let result = run(controls());
    let b = result.balance();
    let actual_temperatures = [
        result.state().heavy_temperature(),
        result.state().electron_temperature(),
    ];
    for i in 0..2 {
        close(
            actual_temperatures[i].value(),
            expected_temperatures[i],
            2e-8,
            1e-3,
        );
    }
    for (actual, expected) in [
        b.heavy_transport_energy_out_j,
        b.electron_transport_energy_out_j,
    ]
    .into_iter()
    .zip(expected_out)
    {
        close(actual, expected, 8e-8, 1e-6);
    }
    assert_eq!(result.state().amounts(), initial.amounts());
    assert_eq!(result.volume_history(), Some(&geometry));
    let energies_before = [
        initial.heavy_internal_energy().value(),
        initial.electron_internal_energy().value(),
    ];
    let work = [b.boundary_work_on_heavy_j, b.boundary_work_on_electrons_j];
    for i in 0..2 {
        close(
            work[i],
            capacities[i] * expected_temperatures[i] - energies_before[i] + expected_out[i],
            2e-7,
            1e-6,
        );
        assert!(work[i] > 0.0);
    }
    check_accounts(&result);
    // The nonzero exterior makes geometric and confinement maps noncommuting.
    // Fixed-step refinement therefore tests their coupled second-order order.
    let mut errors = Vec::new();
    for count in [4, 8, 16] {
        let mut c = controls();
        c.relative_tolerance = 0.1;
        c.absolute_temperature_tolerance_k = 1e6;
        c.energy_tolerance_j = 1e3;
        c.initial_step_s = duration / f64::from(count);
        c.minimum_step_s = c.initial_step_s;
        c.maximum_step_s = c.initial_step_s;
        let result = run(c);
        assert_eq!(result.balance().rejected_steps, 0);
        errors.push(
            (result.state().heavy_temperature().value() - expected_temperatures[0]).abs()
                + (result.state().electron_temperature().value() - expected_temperatures[1]).abs(),
        );
    }
    for pair in errors.windows(2) {
        assert!(
            (3.7..4.3).contains(&(pair[0] / pair[1])),
            "errors={errors:?}"
        );
    }
}

#[test]
fn g3_confinement_is_shared_by_all_light_ion_channels_and_original_thermal_sources() {
    let (plasma, initial) = fixture::plasma(MIXED_DENSITY, 0.025, [15.0, 4.0]);
    let mut transport = transport(&initial);
    transport.energy_confinement = Some(confinement());
    let result = advance_fusion_network(
        &initial,
        &plasma,
        FusionNetworkChannels::ALL,
        mixed_heating(),
        &transport,
        Time::new(0.05),
        controls(),
        || false,
    )
    .unwrap();
    let b = result.balance();
    assert!(b.reaction_extents_mol.iter().all(|extent| *extent > 0.0));
    assert!(b.heavy_transport_energy_out_j > 0.0);
    assert!(b.electron_transport_energy_out_j > 0.0);
    assert!(b.photon_energy_escaped_j > 0.0);
    assert!(b.electron_to_heavy_coulomb_energy_j != 0.0);
    let delta = result.state().heavy_internal_energy().value()
        - initial.heavy_internal_energy().value()
        + result.state().electron_internal_energy().value()
        - initial.electron_internal_energy().value();
    close(
        delta
            + b.neutron_energy_escaped_j
            + b.photon_energy_escaped_j
            + b.heavy_transport_energy_out_j
            + b.electron_transport_energy_out_j,
        b.nuclear_energy_released_j,
        3e-12,
        1e-7,
    );
    check_accounts(&result);
}
