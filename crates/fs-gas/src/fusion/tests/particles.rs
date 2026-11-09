//! G0/G1/G3 explicit fuel/ash ports, calorics and conservative restarts.

use super::*;
use fs_qty::Frequency;

fn transport(
    injection: [f64; 5],
    removal: [f64; 5],
    temperatures_kev: [f64; 2],
) -> DtFusionTransport {
    DtFusionTransport {
        particles: Some(
            FusionParticleExchange::new(
                injection,
                removal.map(Frequency::new),
                Temperature::new(temperatures_kev[0] * KELVIN_PER_KEV),
                Temperature::new(temperatures_kev[1] * KELVIN_PER_KEV),
            )
            .unwrap(),
        ),
        ..DtFusionTransport::default()
    }
}

#[test]
fn g1_particles_injection_mixes_source_calorics_and_companion_electrons() {
    let (plasma, initial) = source::plasma([8e19, 0.0, 2e19], 1.5, [12.0, 3.0]);
    let duration = 0.2;
    let rate = [0.0, 2e-4, 0.0, 0.0, 1e-4];
    let transport = transport(rate, [0.0; 5], [2.0, 1.0]);
    let result = advance_dt_fusion_with_transport(
        &initial,
        &plasma,
        heating(0.8),
        &transport,
        Time::new(duration),
        controls(),
        || false,
    )
    .unwrap();
    let b = result.balance();
    let p = b.particles;
    let [d, _, he4] = plasma.indices();
    for (slot, index) in [(1, d), (4, he4)] {
        close(
            result.state().amounts()[index].value(),
            initial.amounts()[index].value() + rate[slot] * duration,
            3e-14,
            1e-20,
        );
    }
    let added = [
        duration * (rate[1] + rate[4]),
        duration * (rate[1] + 2.0 * rate[4]),
    ];
    let before = [
        initial.heavy_internal_energy().value(),
        initial.electron_internal_energy().value(),
    ];
    let after = [
        result.state().heavy_internal_energy().value(),
        result.state().electron_internal_energy().value(),
    ];
    let temperatures = [
        result.state().heavy_temperature().value(),
        result.state().electron_temperature().value(),
    ];
    let amounts = [
        initial.amounts()[d].value() + initial.amounts()[he4].value(),
        initial.amounts()[plasma.electron_index()].value(),
    ];
    for pool in 0..2 {
        let incoming = CV_MOLAR * [2.0, 1.0][pool] * KELVIN_PER_KEV * added[pool];
        close(p.injected_thermal_energy_j[pool], incoming, 3e-14, 1e-7);
        close(after[pool] - before[pool], incoming, 3e-13, 1e-7);
        close(
            temperatures[pool],
            (before[pool] + incoming) / (CV_MOLAR * (amounts[pool] + added[pool])),
            3e-14,
            1e-5,
        );
    }
    close(p.injected_electrons_mol, added[1], 3e-14, 1e-20);
    assert_eq!(p.exhausted_ions_mol, [0.0; 5]);
    assert_eq!(b.reaction_extent_mol, 0.0);
    assert!(b.total_energy_residual_j.abs() < 1e-6);
    assert_eq!(result.state().thermodynamics(), initial.thermodynamics());
}

#[test]
fn g1_particles_exhaust_keeps_exact_survivors_and_exports_current_thermal_energy() {
    let (plasma, initial) = source::plasma([8e19, 0.0, 2e19], 1.0, [12.0, 3.0]);
    let duration = 0.3;
    let hazards = [0.0, 3.0, 0.0, 0.0, 120.0];
    let transport = transport([0.0; 5], hazards, [2.0, 1.0]);
    let mut c = controls();
    c.initial_step_s = duration;
    c.maximum_step_s = duration;
    c.maximum_steps = 1;
    let result = advance_dt_fusion_with_transport(
        &initial,
        &plasma,
        heating(0.8),
        &transport,
        Time::new(duration),
        c,
        || false,
    )
    .unwrap();
    assert_eq!(result.balance().accepted_steps, 1);
    assert_eq!(result.balance().rejected_steps, 0);
    let p = result.balance().particles;
    let mut removed = [0.0; 2];
    for (slot, index, charge) in [(1, plasma.indices()[0], 1.0), (4, plasma.indices()[2], 2.0)] {
        let before = initial.amounts()[index].value();
        let expected = before * (-hazards[slot] * duration).exp();
        close(
            result.state().amounts()[index].value(),
            expected,
            2e-13,
            0.0,
        );
        let loss = before * -(-hazards[slot] * duration).exp_m1();
        close(p.exhausted_ions_mol[slot], loss, 2e-13, 1e-20);
        removed[0] += loss;
        removed[1] += charge * loss;
    }
    assert!(result.state().amounts()[plasma.indices()[2]].value() > 0.0);
    for (pool, temperature) in [initial.heavy_temperature(), initial.electron_temperature()]
        .into_iter()
        .enumerate()
    {
        close(
            p.exhausted_thermal_energy_j[pool],
            CV_MOLAR * temperature.value() * removed[pool],
            3e-13,
            1e-6,
        );
    }
    close(
        result.state().heavy_temperature().value(),
        initial.heavy_temperature().value(),
        3e-13,
        1e-5,
    );
    close(
        result.state().electron_temperature().value(),
        initial.electron_temperature().value(),
        3e-13,
        1e-5,
    );
    close(p.exhausted_electrons_mol, removed[1], 2e-13, 1e-20);
    assert_eq!(p.injected_ions_mol, [0.0; 5]);
    assert!(result.balance().relative_baryon_residual.abs() < 1e-13);
    assert!(result.balance().relative_charge_residual.abs() < 1e-13);
}

#[test]
fn g3_particles_replenish_an_absent_reactant_and_restart_with_ash_exhaust() {
    let (plasma, initial) = source::plasma([5e19, 0.0, 1e18], 1.0, [10.0, 5.0]);
    let transport = transport(
        [0.0, 2e-5, 5e-5, 0.0, 0.0],
        [0.0, 0.2, 0.3, 0.0, 10.0],
        [2.0, 1.0],
    );
    let run = |state, duration| {
        advance_dt_fusion_with_transport(
            state,
            &plasma,
            heating(0.8),
            &transport,
            Time::new(duration),
            controls(),
            || false,
        )
        .unwrap()
    };
    let full = run(&initial, 0.05);
    let half = run(&initial, 0.025);
    let resumed = run(half.state(), 0.025);
    assert!(full.balance().reaction_extent_mol > 0.0);
    assert!(full.balance().particles.exhausted_ions_mol[4] > 0.0);
    assert!(full.state().amounts()[plasma.indices()[1]].value() > 0.0);
    for (actual, expected) in resumed.state().amounts().iter().zip(full.state().amounts()) {
        close(actual.value(), expected.value(), 3e-8, 1e-18);
    }
    close(
        resumed.state().heavy_temperature().value(),
        full.state().heavy_temperature().value(),
        3e-8,
        1e-3,
    );
    close(
        resumed.state().electron_temperature().value(),
        full.state().electron_temperature().value(),
        3e-8,
        1e-3,
    );
    close(
        half.balance().reaction_extent_mol + resumed.balance().reaction_extent_mol,
        full.balance().reaction_extent_mol,
        2e-7,
        1e-18,
    );
    let b = full.balance();
    let net = b.particles.injected_thermal_energy_j.iter().sum::<f64>()
        - b.particles.exhausted_thermal_energy_j.iter().sum::<f64>();
    close(
        full.state().internal_energy().value() - initial.internal_energy().value()
            + b.neutron_energy_escaped_j,
        b.nuclear_energy_released_j + net,
        3e-12,
        1e-6,
    );
}

#[test]
fn g0_particles_refuse_invalid_or_unsupported_sources_and_cancel_atomically() {
    let (plasma, initial) = source::plasma([5e19, 5e19, 0.0], 1.0, [10.0, 5.0]);
    let original = initial.clone();
    for bad in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(
            FusionParticleExchange::new(
                [bad; 5],
                [Frequency::new(0.0); 5],
                Temperature::new(1e6),
                Temperature::new(1e6)
            )
            .is_err()
        );
        assert!(
            FusionParticleExchange::new(
                [0.0; 5],
                [Frequency::new(bad); 5],
                Temperature::new(1e6),
                Temperature::new(1e6)
            )
            .is_err()
        );
    }
    for (injection, temperature) in [
        ([1e-5, 0.0, 0.0, 0.0, 0.0], 2.0),
        ([0.0, 1e-5, 0.0, 0.0, 0.0], 1e-6),
    ] {
        let law = transport(injection, [0.0; 5], [temperature, 1.0]);
        assert!(
            advance_dt_fusion_with_transport(
                &initial,
                &plasma,
                heating(0.8),
                &law,
                Time::new(0.01),
                controls(),
                || false
            )
            .is_err()
        );
    }
    let hot_electrons = transport([0.0, 1e-5, 0.0, 0.0, 0.0], [0.0; 5], [2.0, 100.0]);
    assert_eq!(
        advance_dt_fusion_with_transport(
            &initial,
            &plasma,
            heating(0.8),
            &hot_electrons,
            Time::new(0.01),
            controls(),
            || false
        ),
        Err(DtFusionError::OutsideDomain)
    );
    let law = transport(
        [0.0, 1e-5, 1e-5, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 2.0],
        [2.0, 1.0],
    );
    let mut calls = 0;
    advance_dt_fusion_with_transport(
        &initial,
        &plasma,
        heating(0.8),
        &law,
        Time::new(0.01),
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut cancelled_calls = 0;
    assert_eq!(
        advance_dt_fusion_with_transport(
            &initial,
            &plasma,
            heating(0.8),
            &law,
            Time::new(0.01),
            controls(),
            || {
                cancelled_calls += 1;
                cancelled_calls == calls
            }
        ),
        Err(DtFusionError::Cancelled)
    );
    assert_eq!(initial, original);
}

#[test]
fn g0_particles_local_error_retains_a_transfer_below_cumulative_roundoff() {
    use crate::fusion::particles::Accounts;
    let mut before = Accounts::default();
    before.exhausted[4] = solve::Sum::new(1e20);
    before.energy_out[0] = solve::Sum::new(1e30);
    let mut coarse = before;
    let mut fine = before;
    coarse.exhausted[4].add(1.0);
    fine.exhausted[4].add(2.0);
    coarse.energy_out[0].add(3.0);
    fine.energy_out[0].add(5.0);
    assert_eq!(coarse.exhausted[4].value, fine.exhausted[4].value);
    assert_eq!(coarse.energy_out[0].value, fine.energy_out[0].value);
    let mut differences = Vec::new();
    before.compare(coarse, fine, 0.1, 0.2, |a, b, tolerance| {
        if a != b {
            differences.push((a, b, tolerance));
        }
    });
    assert_eq!(differences, [(1.0, 2.0, 0.1), (3.0, 5.0, 0.2)]);
}

#[test]
fn g1_particles_stiff_replenishment_retains_the_forced_steady_inventory() {
    // N' = S - lambda*N with S = lambda*N(0) has a constant solution.
    // At 400 lifetimes, unresolved split injection/exhaust can make both
    // step-doubling estimates agree on a nearly empty, incorrect endpoint.
    let (plasma, initial) = source::plasma([8e19, 0.0, 2e19], 1.0, [12.0, 3.0]);
    let [d, _, he4] = plasma.indices();
    let duration = 0.04;
    let hazard = 1e4;
    let stock = initial.amounts()[d].value();
    let supply = hazard * stock;
    let law = transport(
        [0.0, supply, 0.0, 0.0, 0.0],
        [0.0, hazard, 0.0, 0.0, 0.0],
        [12.0, 3.0],
    );
    let mut c = controls();
    c.relative_tolerance = 1e-4;
    c.energy_tolerance_j = 1e-5;
    c.initial_step_s = duration;
    c.maximum_step_s = duration;
    c.minimum_step_s = 1e-9;
    c.maximum_steps = 50_000;
    let result = advance_dt_fusion_with_transport(
        &initial,
        &plasma,
        heating(0.8),
        &law,
        Time::new(duration),
        c,
        || false,
    )
    .unwrap();
    let b = result.balance();
    let p = b.particles;
    assert!(b.accepted_steps > 1 && b.rejected_steps > 0);
    assert_eq!(b.reaction_extent_mol, 0.0);
    close(result.state().amounts()[d].value(), stock, 1e-3, 0.0);
    close(
        result.state().amounts()[he4].value(),
        initial.amounts()[he4].value(),
        1e-13,
        0.0,
    );
    close(
        result.state().amounts()[plasma.electron_index()].value(),
        initial.amounts()[plasma.electron_index()].value(),
        1e-3,
        0.0,
    );
    close(p.injected_ions_mol[1], supply * duration, 1e-12, 1e-18);
    close(p.exhausted_ions_mol[1], supply * duration, 3e-6, 1e-18);
    close(p.injected_electrons_mol, p.injected_ions_mol[1], 1e-13, 1e-18);
    close(p.exhausted_electrons_mol, p.exhausted_ions_mol[1], 1e-13, 1e-18);
    for (pool, (actual, expected)) in [
        (result.state().heavy_temperature(), initial.heavy_temperature()),
        (
            result.state().electron_temperature(),
            initial.electron_temperature(),
        ),
    ]
    .into_iter()
    .enumerate()
    {
        close(actual.value(), expected.value(), 3e-9, 1e-3);
        close(
            p.injected_thermal_energy_j[pool],
            CV_MOLAR * expected.value() * supply * duration,
            1e-12,
            1e-6,
        );
        close(
            p.exhausted_thermal_energy_j[pool],
            CV_MOLAR * expected.value() * supply * duration,
            3e-6,
            1e-6,
        );
    }
    // A budget that cannot resolve the driven loss must refuse the solve.
    c.maximum_steps = 1;
    assert_eq!(
        advance_dt_fusion_with_transport(
            &initial,
            &plasma,
            heating(0.8),
            &law,
            Time::new(duration),
            c,
            || false,
        ),
        Err(DtFusionError::StepBudget)
    );
}
