//! Independent gamma-law, variable-Cp, storage and cancellation regressions.
use super::*;
use crate::multicomponent::characteristic_tests::caloric_gas;
use crate::multicomponent::{Boundaries, Boundary, Conserved};
use super::super::{BranchEnd, PlenumConnection};
use fs_qty::Length;

mod pressure;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "actual={a:.16e}, expected={b:.16e}, tolerance={tolerance:e}");
}
fn orifice(area: f64) -> Orifice {
    Orifice::new(Area::new(area), 1.0, 1.0, ValveDirection::Bidirectional).unwrap()
}
fn gamma_rate(area: f64, pressure: f64, temperature: f64, backpressure: f64) -> f64 {
    if backpressure >= pressure { return 0.0; }
    let gamma = 1.4_f64;
    let critical = (2.0 / (gamma + 1.0)).powf(gamma / (gamma - 1.0));
    let ratio = (backpressure / pressure).max(critical);
    area * pressure * (2.0 * gamma / ((gamma - 1.0) * temperature)
        * (ratio.powf(2.0 / gamma) - ratio.powf((gamma + 1.0) / gamma))).sqrt()
}

#[test]
fn gamma_nozzle_matches_independent_subcritical_and_choked_flow_with_reverse_donor() {
    let gas = caloric_gas(false, 0.0);
    let valve = Orifice::new(Area::new(0.1), 0.7, 0.6, ValveDirection::Bidirectional).unwrap();
    let donor = gas.conserved(&[1.0, 0.0], 0.0, 2.0).unwrap();
    for backpressure in [1.9, 1.5, 0.2, 0.02] {
        let receiver = gas.conserved(&[0.0, backpressure], 0.0, backpressure).unwrap();
        let flow = valve.flow(&gas, donor.as_cell(), receiver.as_cell(), || false).unwrap();
        let expected = gamma_rate(0.042, 2.0, 2.0, backpressure);
        close(flow.mass_kg_s, expected, 2e-11);
        assert_eq!(flow.species_kg_s[1], 0.0);
        assert_eq!(flow.donor, Some(0));
        close(flow.energy_w / flow.mass_kg_s, 7.0 - 2.5 * 0.001, 2e-12);
        assert_eq!(flow.choked, backpressure < 1.0);
        if flow.choked {
            let exit = flow.exit.unwrap();
            close(exit.velocity, exit.sound_speed, 2e-11);
        }
        let reverse = valve.flow(&gas, receiver.as_cell(), donor.as_cell(), || false).unwrap();
        assert_eq!(reverse.mass_kg_s, -flow.mass_kg_s);
        assert_eq!(reverse.energy_w, -flow.energy_w);
        assert_eq!(reverse.donor, Some(1));
    }
}

#[test]
fn variable_cp_nozzle_preserves_original_enthalpy_entropy_and_energy_gauge() {
    let mut baseline = Vec::new();
    for reference in [0.0, 1e18] {
        let gas = caloric_gas(true, reference);
        let donor = gas.conserved(&[1.0, 0.0], 0.0, 8.0).unwrap();
        for (i, backpressure) in [6.0, 0.01].into_iter().enumerate() {
            let receiver = gas.conserved(&[backpressure, 0.0], 0.0, backpressure).unwrap();
            let flow = orifice(0.1).flow(&gas, donor.as_cell(), receiver.as_cell(), || false).unwrap();
            let exit = flow.exit.unwrap();
            let t = exit.temperature;
            close(3.5 * (8.0 - t) + 0.125 * (64.0 - t * t), 0.5 * exit.velocity.powi(2), 2e-10);
            close(3.5 * (t / 8.0).ln() + 0.25 * (t - 8.0) - (exit.pressure / 8.0).ln(), 0.0, 2e-11);
            if reference == 0.0 { baseline.push(flow); } else { assert_eq!(flow, baseline[i]); }
        }
    }
}

#[test]
fn closed_and_check_valves_have_no_leakage_and_invalid_inputs_refuse() {
    let gas = caloric_gas(false, 0.0);
    let high = gas.conserved(&[1.0, 0.0], 0.0, 2.0).unwrap();
    let low = gas.conserved(&[0.0, 0.1], 0.0, 0.1).unwrap();
    for valve in [
        orifice(0.0),
        Orifice::new(Area::new(1.0), 0.0, 1.0, ValveDirection::Bidirectional).unwrap(),
        Orifice::new(Area::new(1.0), 1.0, 0.0, ValveDirection::Bidirectional).unwrap(),
        Orifice::new(Area::new(1.0), 1.0, 1.0, ValveDirection::SecondToFirst).unwrap(),
    ] {
        assert_eq!(valve.flow(&gas, high.as_cell(), low.as_cell(), || false).unwrap(), ValveFlow::zero(2));
    }
    assert_eq!(orifice(1.0).flow(&gas, high.as_cell(), high.as_cell(), || false).unwrap(), ValveFlow::zero(2));
    for (a, cd, opening) in [(-1.0, 1.0, 1.0), (1.0, 1.1, 1.0), (1.0, 1.0, f64::NAN)] {
        assert!(Orifice::new(Area::new(a), cd, opening, ValveDirection::Bidirectional).is_err());
    }
    let moving = gas.conserved(&[1.0, 0.0], 0.1, 2.0).unwrap();
    assert!(orifice(1.0).flow(&gas, moving.as_cell(), low.as_cell(), || false).is_err());
}

fn tanks(gas: &MixtureGas) -> [Plenum; 2] {
    [Plenum::new(0.2, gas.conserved(&[1.0, 0.0], 0.0, 2.0).unwrap()).unwrap(),
     Plenum::new(0.1, gas.conserved(&[0.1, 0.0], 0.0, 0.1).unwrap()).unwrap()]
}
fn apply(plenums: &mut [Plenum; 2], update: &pair::Update) {
    for (i, p) in plenums.iter_mut().enumerate() {
        p.state.partial_densities.copy_from_slice(&update.states[i][..2]);
        p.state.energy = update.states[i][3];
    }
}

#[test]
fn finite_transfer_matches_adiabatic_blowdown_and_closes_actual_species_and_energy() {
    let gas = caloric_gas(false, 0.0);
    let initial = tanks(&gas);
    for duration in [0.001, 0.1, 10.0] {
        let update = pair::solve(&gas, [&initial[0], &initial[1]], orifice(0.1), duration,
            ValveControl::default(), &mut || false).unwrap();
        let x = update.exchange.mass_to_second_kg();
        assert!(x > 0.0 && x < 0.2);
        let remaining = 1.0 - x / 0.2;
        let donor = gas.primitive(cell(&update.states[0][..4], 2), || false).unwrap();
        let receiver = gas.primitive(cell(&update.states[1][..4], 2), || false).unwrap();
        close(donor.temperature, 2.0 * remaining.powf(0.4), 3e-11);
        close(donor.pressure, 2.0 * remaining.powf(1.4), 3e-11);
        assert!(donor.pressure > receiver.pressure);
        // Physical donor U0=p0*V/(gamma-1)=1 J; its formation-free reference
        // contribution travels with the exact same species transfer.
        close(update.exchange.energy_to_second_j + 0.0025 * x, 1.0 - remaining.powf(1.4), 3e-11);
        let mut changes = [0.0; 2];
        let mut energy_change = 0.0;
        for side in 0..2 {
            let q = cell(&update.states[side][..4], 2);
            let old = initial[side].state();
            let volume = initial[side].volume_m3();
            for k in 0..2 { changes[k] += volume * (q.partial_densities[k] - old.partial_densities[k]); }
            energy_change += volume * (gas.absolute_energy_density(q).unwrap() - gas.absolute_energy_density(old).unwrap());
        }
        for change in changes { close(change, 0.0, 3e-14); }
        close(energy_change, 0.0, 3e-13);
        close(x, duration * gamma_rate(0.1, donor.pressure, donor.temperature, receiver.pressure), 3e-10);
    }
}

#[test]
fn finite_transfer_refines_toward_independent_continuous_two_tank_dynamics() {
    let rate = |x: f64| {
        let fraction = (1.0 - x / 0.2).powf(1.4);
        gamma_rate(0.1, 2.0 * fraction, 2.0 * (1.0 - x / 0.2).powf(0.4),
            0.1 + 4.0 * (1.0 - fraction))
    };
    let end = 0.2;
    let h = end / 10_000.0;
    let mut exact = 0.0;
    for _ in 0..10_000 {
        let a = rate(exact);
        let b = rate(exact + 0.5 * h * a);
        let c = rate(exact + 0.5 * h * b);
        let d = rate(exact + h * c);
        exact += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let gas = caloric_gas(false, 0.0);
    let mut errors = Vec::new();
    for count in [8, 16] {
        let mut p = tanks(&gas);
        for _ in 0..count {
            let update = pair::solve(&gas, [&p[0], &p[1]], orifice(0.1), end / count as f64,
                ValveControl::default(), &mut || false).unwrap();
            apply(&mut p, &update);
        }
        errors.push((0.2 - p[0].state().density() * 0.2 - exact).abs());
    }
    assert!(errors[1] < 0.56 * errors[0], "first-order errors={errors:?}");
}

fn network(gas: &MixtureGas) -> MixtureNetwork1d {
    let branches = [([1.0, 0.0], 2.0), ([0.0, 1.0], 1.0)].into_iter().map(|(rho, p)| {
        MixtureEuler1d::new(gas.clone(), Length::new(0.1), 0.1,
            vec![gas.conserved(&rho, 0.0, p).unwrap(); 4],
            Boundaries::Ends { left: Boundary::Reflecting, right: Boundary::Reflecting },
            || false).unwrap()
    }).collect();
    MixtureNetwork1d::new(branches,
        vec![Plenum::new(0.02, gas.conserved(&[1.0, 0.0], 0.0, 2.0).unwrap()).unwrap(),
             Plenum::new(0.03, gas.conserved(&[0.0, 1.0], 0.0, 1.0).unwrap()).unwrap()],
        vec![PlenumConnection { branch: 0, end: BranchEnd::Right, plenum: 0 },
             PlenumConnection { branch: 1, end: BranchEnd::Left, plenum: 1 }], || false).unwrap()
}
fn owner() -> ValvedNetwork1d {
    let gas = caloric_gas(true, 0.0);
    ValvedNetwork1d::new(network(&gas),
        vec![ValveLink { first: 0, second: 1, orifice: orifice(0.01) }], || false).unwrap()
}
fn fields(flow: &ValvedNetwork1d) -> Vec<Conserved> {
    let owned = |q: Cell<'_>| Conserved { partial_densities: q.partial_densities.to_vec(), momentum: q.momentum, energy: q.energy };
    flow.branches().iter().flat_map(|b| (0..b.cell_count()).map(move |i| owned(b.cell(i).unwrap())))
        .chain(flow.plenums().iter().map(|p| owned(p.state()))).collect()
}

#[test]
fn valve_drives_connected_waves_and_independent_plenum_budgets_close() {
    let mut flow = owner();
    let initial = flow.totals();
    for _ in 0..6 {
        let old: Vec<_> = flow.plenums().iter().map(|p| p.state.clone()).collect();
        let report = flow.step(Time::new(0.002), ValveControl::default(), || false).unwrap();
        let valve = flow.accepted_valve_exchanges().unwrap()[0];
        let ports = flow.accepted_port_exchanges().unwrap();
        assert!(valve.mass_to_second_kg() > 0.0);
        for side in 0..2 {
            let sign = if side == 0 { -1.0 } else { 1.0 };
            let p = &flow.plenums()[side];
            for k in 0..2 {
                close(p.volume_m3() * (p.state().partial_densities[k] - old[side].partial_densities[k]),
                    sign * valve.species_to_second_kg[k] + ports[side].species_masses_to_plenum[k], 3e-13);
            }
            close(p.volume_m3() * (p.state().energy - old[side].energy),
                sign * valve.energy_to_second_j + ports[side].energy_to_plenum, 3e-12);
            assert_eq!(p.state().momentum, 0.0);
        }
        assert_eq!(report.elapsed, report.transport.elapsed);
        close(report.absolute_energy_residual_j, 0.0, 4e-12);
        for b in flow.branches() { assert_eq!(b.time(), flow.time()); }
    }
    assert!(flow.plenums()[1].state().partial_densities[0] > 0.0);
    assert!(flow.branches()[1].totals().species_masses[0] > 0.0);
    for k in 0..2 { close(flow.totals().species_masses[k], initial.species_masses[k], 3e-12); }
    close(flow.totals().absolute_energy(flow.gas()).unwrap(), initial.absolute_energy(flow.gas()).unwrap(), 3e-11);
    flow.set_orifice(0, orifice(0.0)).unwrap();
    flow.step(Time::new(0.002), ValveControl::default(), || false).unwrap();
    assert_eq!(flow.accepted_valve_exchanges().unwrap()[0].mass_to_second_kg(), 0.0);
    flow.set_orifice(0, orifice(0.01)).unwrap();
    flow.step(Time::new(0.002), ValveControl::default(), || false).unwrap();
    assert!(flow.accepted_valve_exchanges().unwrap()[0].mass_to_second_kg() > 0.0);
    let raw = flow.into_network();
    assert!(raw.last_report().is_none());
}

#[test]
fn cancellation_and_exhausted_roots_preserve_all_accepted_fields_and_reports() {
    let make = || { let mut f = owner(); f.step(Time::new(0.002), ValveControl::default(), || false).unwrap(); f };
    let mut reference = make();
    let mut polls = 0;
    reference.step(Time::new(0.002), ValveControl::default(), || { polls += 1; false }).unwrap();
    for stop in [1, polls / 2, polls] {
        let mut f = make();
        let before = fields(&f);
        let time = f.time();
        let report = f.last_report();
        let exchanges = f.accepted_valve_exchanges().unwrap().to_vec();
        let ports = f.accepted_port_exchanges().unwrap().to_vec();
        let branches = f.accepted_branch_transport_reports().unwrap().to_vec();
        let mut count = 0;
        assert_eq!(f.step(Time::new(0.002), ValveControl::default(), || { count += 1; count >= stop }), Err(ValveError::Cancelled));
        assert_eq!(fields(&f), before);
        assert_eq!(f.time(), time);
        assert_eq!(f.last_report(), report);
        assert_eq!(f.accepted_valve_exchanges().unwrap(), exchanges);
        assert_eq!(f.accepted_port_exchanges().unwrap(), ports);
        assert_eq!(f.accepted_branch_transport_reports().unwrap(), branches);
        let limited = ValveControl { maximum_iterations: 1, maximum_coupling_retries: 0, ..ValveControl::default() };
        assert!(f.step(Time::new(0.002), limited, || false).is_err());
        assert_eq!(fields(&f), before);
        assert_eq!(f.last_report(), report);
        f.step(Time::new(0.002), ValveControl::default(), || false).unwrap();
        for (a, b) in fields(&f).iter().zip(fields(&reference)) {
            for k in 0..2 { close(a.partial_densities[k], b.partial_densities[k], 3e-11); }
            close(a.momentum, b.momentum, 3e-11);
            close(a.energy, b.energy, 3e-10);
        }
    }
}

#[test]
fn topology_and_setting_refusals_leave_the_existing_owner_unchanged() {
    let gas = caloric_gas(false, 0.0);
    for (first, second) in [(0, 0), (0, 2), (2, 0)] {
        assert!(ValvedNetwork1d::new(network(&gas),
            vec![ValveLink { first, second, orifice: orifice(0.1) }], || false).is_err());
    }
    let mut flow = owner();
    let old = flow.links().to_vec();
    assert!(flow.set_orifice(1, orifice(0.0)).is_err());
    assert_eq!(flow.links(), old);
    assert_eq!(flow.time().value(), 0.0);
}
