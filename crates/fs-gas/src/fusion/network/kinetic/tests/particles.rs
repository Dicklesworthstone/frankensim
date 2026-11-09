//! Thermal ash ports must not export or prematurely heat kinetic alpha particles.

use super::*;
use crate::fusion::FusionParticleExchange;
use fs_qty::Frequency;

fn particle_transport() -> FusionNetworkTransport {
    FusionNetworkTransport {
        particles: Some(
            FusionParticleExchange::new(
                [0.0, 2e-5, 2e-5, 0.0, 1e-6],
                [0.0, 0.0, 0.0, 0.0, 80.0].map(Frequency::new),
                Temperature::new(2.0 * KELVIN_PER_KEV),
                Temperature::new(KELVIN_PER_KEV),
            )
            .unwrap(),
        ),
        ..FusionNetworkTransport::default()
    }
}

#[test]
fn g3_particles_exhaust_only_thermal_ash_while_fast_alphas_slow_and_restart() {
    let (plasma, reference) = fixture::plasma([0.0, 5e19, 5e19, 0.0, 1.01e18], 1.0, [10.0, 5.0]);
    let grid = grid(5000.0);
    let mut populations = vec![Amount::new(0.0); grid.len()];
    populations[90] = Amount::new(1e16 / AVOGADRO);
    let initial = KineticFusionState::new(
        &reference,
        &plasma,
        &grid,
        &populations,
        controls().plasma_limits,
        || false,
    )
    .unwrap();
    let original = initial.clone();
    let alpha = alpha_transport(collision_source(&reference, plasma.indices()[4]));
    let thermal = particle_transport();
    let run = |state, duration| {
        advance_kinetic_fusion(
            state,
            FusionNetworkChannels::NONE,
            heating(),
            &thermal,
            &alpha,
            Time::new(duration),
            controls(),
            || false,
        )
        .unwrap()
    };
    let duration = 0.01;
    let result = run(&initial, duration);
    let half = run(&initial, 0.5 * duration);
    let restarted = run(half.state(), 0.5 * duration);
    let ash_index = plasma.indices()[4];
    let ash_before = initial.thermal_amounts()[ash_index].value();
    let expected_ash =
        ash_before * (-80.0_f64 * duration).exp() + 1e-6 / 80.0 * -(-80.0_f64 * duration).exp_m1();
    let b = result.balance();
    let p = b.particles;
    close(
        result.state().thermal_amounts()[ash_index].value(),
        expected_ash,
        1e-6,
        1e-18,
    );
    close(
        result.state().alpha_amount().value(),
        initial.alpha_amount().value(),
        3e-13,
        1e-22,
    );
    assert!(result.state().alpha_energy().value() < initial.alpha_energy().value());
    assert!(result.alpha_balance().deposited_electron_energy_j > 0.0);
    assert_eq!(result.alpha_balance().born_alpha_amount_mol, 0.0);
    assert!(p.exhausted_ions_mol[4] > 0.3 * ash_before);
    close(
        p.exhausted_electrons_mol,
        2.0 * p.exhausted_ions_mol[4],
        3e-13,
        1e-20,
    );
    close(
        result.state().amounts()[ash_index].value(),
        expected_ash + initial.alpha_amount().value(),
        1e-6,
        1e-18,
    );
    close(
        result.state().amounts()[plasma.electron_index()].value(),
        initial.amounts()[plasma.electron_index()].value() + p.injected_electrons_mol
            - p.exhausted_electrons_mol,
        3e-13,
        1e-20,
    );
    for (actual, expected) in restarted
        .state()
        .thermal_amounts()
        .into_iter()
        .zip(result.state().thermal_amounts())
    {
        close(actual.value(), expected.value(), 2e-6, 1e-18);
    }
    close(
        restarted.state().background_heavy_temperature().value(),
        result.state().background_heavy_temperature().value(),
        2e-6,
        1e-3,
    );
    close(
        restarted.state().electron_temperature().value(),
        result.state().electron_temperature().value(),
        2e-6,
        1e-3,
    );
    close(
        restarted.state().alpha_energy().value(),
        result.state().alpha_energy().value(),
        2e-5,
        1e-6,
    );
    for completed in [&result, &half, &restarted] {
        check_accounts(completed);
    }
    assert_eq!(initial, original);
}

#[test]
fn g1_particles_incoming_helium_is_thermal_and_matches_the_common_network_owner() {
    let (plasma, reference) = fixture::plasma([0.0, 5e19, 5e19, 0.0, 1e18], 1.0, [10.0, 5.0]);
    let energy_grid = AlphaEnergyGrid::new(
        vec![Energy::new(KEV_J), Energy::new(100.0 * KEV_J)],
        || false,
    )
    .unwrap();
    let initial = KineticFusionState::from_thermal(
        &reference,
        &plasma,
        &energy_grid,
        controls().plasma_limits,
        || false,
    )
    .unwrap();
    let alpha = alpha_transport(collision_source(&reference, plasma.indices()[4]));
    let ash_index = plasma.indices()[4];
    let ash_before = reference.amounts()[ash_index].value();
    let hazard = 1e4;
    let thermal = FusionNetworkTransport {
        particles: Some(
            FusionParticleExchange::new(
                [0.0, 0.0, 0.0, 0.0, hazard * ash_before],
                [0.0, 0.0, 0.0, 0.0, hazard].map(Frequency::new),
                reference.heavy_temperature(),
                reference.electron_temperature(),
            )
            .unwrap(),
        ),
        ..FusionNetworkTransport::default()
    };
    let duration = Time::new(0.04);
    let mut c = controls();
    c.relative_tolerance = 1e-4;
    c.energy_tolerance_j = 1e-5;
    c.initial_step_s = duration.value();
    c.maximum_step_s = duration.value();
    c.minimum_step_s = 1e-9;
    c.maximum_steps = 50_000;
    let kinetic = advance_kinetic_fusion(
        &initial,
        FusionNetworkChannels::NONE,
        heating(),
        &thermal,
        &alpha,
        duration,
        c,
        || false,
    )
    .unwrap();
    let maxwellian = crate::fusion::network::advance_fusion_network(
        &reference,
        &plasma,
        FusionNetworkChannels::NONE,
        heating(),
        &thermal,
        duration,
        c,
        || false,
    )
    .unwrap();
    // An exact steady thermal ash stock must survive 400 removal lifetimes.
    // Comparing only the two owners would miss their shared saturated error.
    close(
        kinetic.state().thermal_amounts()[ash_index].value(),
        ash_before,
        1e-3,
        0.0,
    );
    assert!(kinetic.balance().accepted_steps > 1 && kinetic.balance().rejected_steps > 0);
    assert_eq!(kinetic.state().alpha_amount().value(), 0.0);
    assert_eq!(kinetic.state().alpha_energy().value(), 0.0);
    assert!(
        kinetic
            .state()
            .alpha_populations_mol()
            .iter()
            .all(|n| *n == 0.0)
    );
    assert!(kinetic.balance().particles.injected_ions_mol[4] > 0.0);
    for (actual, expected) in kinetic
        .state()
        .thermal_amounts()
        .iter()
        .zip(maxwellian.state().amounts())
    {
        close(actual.value(), expected.value(), 2e-12, 1e-20);
    }
    close(
        kinetic.state().background_heavy_temperature().value(),
        maxwellian.state().heavy_temperature().value(),
        3e-12,
        1e-4,
    );
    close(
        kinetic.state().electron_temperature().value(),
        maxwellian.state().electron_temperature().value(),
        3e-12,
        1e-4,
    );
    check_accounts(&kinetic);
}
