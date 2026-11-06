//! G0/G1/G3: native OLD/NEW force epochs, work and time refinement.
use super::*;

fn si(flow: &NsfGrid3<D3q27>, acceleration_lattice: [f64; 3]) -> [f64; 3] {
    let scale = flow.law.dt_s * (flow.law.dt_s / flow.law.dx_m);
    acceleration_lattice.map(|a| a / scale)
}

#[test]
fn changing_uniform_drive_supplies_its_impulse_and_kinetic_work_without_heat() {
    let mut flow = grid([3; 3]);
    let gravity = [0.0003, -0.0002, 0.0001];
    flow.set_gravity_si(si(&flow, gravity)).unwrap();
    let mut old = [0.0; 3];
    let mut sum = 0.0;
    for new in [
        [0.006, -0.004, 0.002],
        [-0.003, 0.001, -0.001],
        [0.0; 3],
        [0.002, 0.003, -0.001],
    ] {
        let before = flow.cell([1; 3]);
        let energy = flow.total_energy_j();
        let endpoint = si(&flow, new);
        let field = if new == [0.0; 3] { Vec::new() } else { vec![endpoint; flow.n] };
        let ledger = flow.step_with_body_acceleration_si(&field, || false).unwrap();
        let after = flow.cell([1; 3]);
        for axis in 0..3 {
            near(after.1[axis] - before.1[axis],
                0.5 * (old[axis] + new[axis]) + gravity[axis], 3e-13);
        }
        let expected = 0.25 * flow.n as f64 * flow.law.energy_unit_j()
            * (0..3).map(|a| (old[a] + new[a])
                * (before.0 * before.1[a] + after.0 * after.1[a])).sum::<f64>();
        near(ledger.body_force_work_j, expected, 3e-12 * expected.abs().max(1e-20));
        near(after.2, before.2, 2e-8);
        near(after.2, 300.0, 2e-8);
        near(flow.total_energy_j() - energy,
            ledger.body_force_work_j + ledger.gravity_work_j, 3e-12 * energy);
        near(ledger.viscous_heating_j, 0.0, 3e-12 * energy);
        assert!(ledger.residual_rel < 3e-12);
        assert_eq!(flow.body_acceleration_m_s2([1; 3]).unwrap(), endpoint);
        sum += ledger.body_force_work_j;
        old = new;
    }
    assert_eq!(flow.cumulative_ledger().body_force_work_j, sum);
    // Acceptance retains the endpoint field for an ordinary held native step.
    let before = flow.cell([1; 3]);
    flow.step().unwrap();
    let after = flow.cell([1; 3]);
    for a in 0..3 {
        near(after.1[a] - before.1[a], old[a] + gravity[a], 3e-13);
    }
    near(after.2, before.2, 2e-8);
}

#[test]
fn identical_endpoint_and_disabled_routes_retain_every_held_step_bit() {
    for enabled in [false, true] {
        let mut held = grid([4; 3]);
        let mut endpoint = grid([4; 3]);
        let field = if enabled {
            (0..held.n).map(|i| [
                300.0, 2000.0 * (core::f64::consts::TAU * held.at(i)[0] as f64 / 4.0).sin(), -100.0,
            ]).collect::<Vec<_>>()
        } else {
            vec![[0.0; 3]; held.n]
        };
        held.set_body_acceleration_si(&field, || false).unwrap();
        endpoint.set_body_acceleration_si(&field, || false).unwrap();
        for _ in 0..3 {
            let h = held.step().unwrap();
            let e = endpoint.step_with_body_acceleration_si(&field, || false).unwrap();
            assert_eq!(h, e);
            assert_eq!(held.state_hash(), endpoint.state_hash());
            assert_eq!(held.force, endpoint.force);
            assert_eq!(held.div, endpoint.div);
            assert_eq!(held.cumulative_ledger(), endpoint.cumulative_ledger());
        }
    }
}

#[test]
fn smooth_time_dependent_drive_converges_at_second_order() {
    let amplitude = 1.0e6;
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let mut specification = spec([3; 3], [true; 3]);
        specification.dx_m = 0.012 / steps as f64;
        let mut flow = NsfGrid3::<D3q27>::new(
            specification, |_| false, |_| (1.0, [0.01, -0.02, 0.015], 300.0),
        ).unwrap();
        let before = flow.cell([1; 3]);
        let end_time = steps as f64 * flow.law.dt_s;
        let mut work = 0.0;
        for step in 1..=steps {
            // A(t)=A*(t/T)^2: exact integral is A*T/3.
            let fraction = step as f64 / steps as f64;
            let field = vec![[0.0, amplitude * fraction * fraction, 0.0]; flow.n];
            let ledger = flow.step_with_body_acceleration_si(&field, || false).unwrap();
            work += ledger.body_force_work_j;
            near(flow.cell([1; 3]).2, 300.0, 3e-8);
            assert!(ledger.residual_rel < 3e-12);
        }
        let after = flow.cell([1; 3]);
        let change = (after.1[1] - before.1[1]) * flow.law.dx_m / flow.law.dt_s;
        let exact = amplitude * end_time / 3.0;
        let error = (change - exact).abs();
        // Independent composite trapezoid error for a quadratic signal.
        near(error, amplitude * end_time / (6.0 * (steps * steps) as f64), 2e-9);
        errors.push(error);
        let kinetic_change = 0.5 * flow.n as f64 * flow.law.energy_unit_j()
            * (0..3).map(|a| after.0 * after.1[a].powi(2)
                - before.0 * before.1[a].powi(2)).sum::<f64>();
        near(work, kinetic_change, 2e-11 * work.abs());
    }
    assert!(errors[1] < 0.251 * errors[0], "{errors:?}");
    assert!(errors[2] < 0.251 * errors[1], "{errors:?}");
}

#[test]
fn spatial_endpoint_forces_generate_native_vorticity_and_disabling_preserves_the_flow() {
    for n in [12, 24] {
        let mut specification = spec([n, 3, 3], [true; 3]);
        specification.dx_m = 0.012 / n as f64;
        let mut flow = NsfGrid3::<D3q27>::new(
            specification, |_| false, |_| (1.0, [0.0; 3], 300.0),
        ).unwrap();
        let k = core::f64::consts::TAU / n as f64;
        let field = |amplitude: f64| (0..n * 9).map(|i| [
            0.0, amplitude * (k * ((i % n) as f64 + 0.5)).sin(), 0.0,
        ]).collect::<Vec<_>>();
        let old_amplitude = 10_000.0;
        let new_amplitude = 20_000.0;
        flow.set_body_acceleration_si(&field(old_amplitude), || false).unwrap();
        let ledger = flow.step_with_body_acceleration_si(&field(new_amplitude), || false).unwrap();
        let stream_filter = 2.0 / 3.0 + k.cos() / 3.0;
        let amplitude = 0.5 * (stream_filter * old_amplitude + new_amplitude);
        let mut norm = 0.0;
        for x in 0..n {
            let phase = k * (x as f64 + 0.5);
            // Independent native D3Q27 collision/stream filter on OLD force;
            // NEW acceleration enters directly through endpoint half force.
            let expected = amplitude * flow.law.dt_s.powi(2) / flow.law.dx_m
                * k.sin() * phase.cos();
            let omega = flow.vorticity_budget([x, 1, 1]).unwrap().flow.vorticity()[2];
            near(omega, expected, 3e-13);
            let curl = flow.body_acceleration_curl_si([x, 1, 1]).unwrap();
            near(curl[2], new_amplitude * k.sin() / flow.law.dx_m * phase.cos(),
                2e-6);
            norm += omega * omega;
        }
        assert!(ledger.body_force_work_j > 0.0);
        assert!(ledger.residual_rel < 3e-12);
        flow.step_with_body_acceleration_si(&[], || false).unwrap();
        assert!(flow.body_acceleration.is_empty());
        assert_eq!(flow.body_acceleration_curl_si([1; 3]).unwrap(), [0.0; 3]);
        let held = flow.step().unwrap();
        assert_eq!(held.body_force_work_j, 0.0);
        let retained: f64 = (0..n).map(|x|
            flow.vorticity_budget([x, 1, 1]).unwrap().flow.vorticity()[2].powi(2)).sum();
        assert!(retained > 0.5 * norm);
    }
}

#[test]
fn unsteady_body_force_uses_actual_mass_and_samples_sources_only_on_the_native_clock() {
    use std::sync::{Arc, atomic::{AtomicUsize, Ordering}};
    let calls = Arc::new(AtomicUsize::new(0));
    let observed = calls.clone();
    let mass_source = 0.001;
    let momentum_source = [0.00002, -0.00001, 0.00003];
    let source_energy = 0.001;
    let mut flow = NsfGrid3::<D3q27>::new_with_sources(
        spec([3; 3], [true; 3]), |_| false,
        |_| (1.0, [0.01, -0.02, 0.015], 300.0),
        Some(Box::new(move |_, _| {
            observed.fetch_add(1, Ordering::SeqCst);
            (mass_source, momentum_source, source_energy)
        })),
    ).unwrap();
    let mut old = [0.0; 3];
    for new in [[0.002, -0.001, 0.003], [-0.001, 0.003, 0.0], [0.0; 3]] {
        let before = flow.cell([1; 3]);
        let sampled = calls.load(Ordering::SeqCst);
        let endpoint = vec![si(&flow, new); flow.n];
        let ledger = flow.step_with_body_acceleration_si(&endpoint, || false).unwrap();
        assert_eq!(calls.load(Ordering::SeqCst) - sampled, 2 * flow.n);
        let after = flow.cell([1; 3]);
        near(after.0 - before.0, mass_source, 3e-14);
        for a in 0..3 {
            near(after.0 * after.1[a] - before.0 * before.1[a],
                0.5 * (before.0 * old[a] + after.0 * new[a]) + momentum_source[a], 3e-13);
        }
        let expected = 0.25 * flow.n as f64 * flow.law.energy_unit_j()
            * (0..3).map(|a| (old[a] + new[a])
                * (before.0 * before.1[a] + after.0 * after.1[a])).sum::<f64>();
        near(ledger.body_force_work_j, expected, 3e-12 * expected.abs().max(1e-20));
        near(ledger.source_energy_j,
            flow.n as f64 * flow.law.energy_unit_j() * source_energy,
            3e-12 * ledger.source_energy_j.abs());
        assert_eq!(ledger.viscous_heating_j, 0.0);
        assert!(ledger.residual_rel < 3e-12);
        old = new;
    }
}

#[test]
fn endpoint_refusal_and_cancellation_leave_every_accepted_epoch_unchanged() {
    let dims = [9, 4, 4];
    let make = || {
        let mut flow = grid(dims);
        flow.set_body_acceleration_si(&vec![[3000.0, -2000.0, 1000.0]; flow.n], || false).unwrap();
        flow
    };
    let field = vec![[-2000.0, 4000.0, -3000.0]; dims.into_iter().product()];
    let mut reference = make();
    let mut polls = 0;
    reference.step_with_body_acceleration_si(&field, || {
        polls += 1;
        false
    }).unwrap();
    assert!(polls > 12);
    for stop in [1, polls / 2, polls] {
        let mut flow = make();
        let before = flow.state_hash();
        let force = flow.force.clone();
        let div = flow.div.clone();
        let ledger = flow.cumulative_ledger();
        let steps = flow.steps();
        let mut calls = 0;
        assert!(flow.step_with_body_acceleration_si(&field, || {
            calls += 1;
            calls == stop
        }).is_err());
        assert_eq!(flow.state_hash(), before);
        assert_eq!(flow.force, force);
        assert_eq!(flow.div, div);
        assert_eq!(flow.cumulative_ledger(), ledger);
        assert_eq!(flow.steps(), steps);
        flow.step_with_body_acceleration_si(&field, || false).unwrap();
        assert_eq!(flow.state_hash(), reference.state_hash());
        assert_eq!(flow.force, reference.force);
        assert_eq!(flow.div, reference.div);
        assert_eq!(flow.cumulative_ledger(), reference.cumulative_ledger());
    }
    let mut flow = make();
    let before = flow.state_hash();
    let accepted_force = flow.force.clone();
    let accepted_ledger = flow.cumulative_ledger();
    for invalid in [
        vec![[0.0; 3]; flow.n - 1],
        vec![[f64::NAN, 0.0, 0.0]; flow.n],
        vec![[f64::from_bits(1), 0.0, 0.0]; flow.n],
        vec![si(&flow, [10.0, 0.0, 0.0]); flow.n],
    ] {
        assert!(flow.step_with_body_acceleration_si(&invalid, || false).is_err());
        assert_eq!(flow.state_hash(), before);
        assert_eq!(flow.force, accepted_force);
        assert_eq!(flow.cumulative_ledger(), accepted_ledger);
    }
}
