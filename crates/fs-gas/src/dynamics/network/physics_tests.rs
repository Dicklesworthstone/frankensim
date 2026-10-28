//! Independent physical checks for the compact-plenum / resolved-passage model.
use super::{BranchEnd, GasNetwork1d, NetworkTotals, Plenum, PlenumConnection};
use crate::dynamics::tests::gas;
use crate::dynamics::{Boundaries, Boundary, Conserved, Euler1d, FrozenGas, Mesh1d, StepControl};
use fs_qty::{Length, Time};

fn walls() -> Boundaries {
    Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    }
}

fn uniform_branch(
    gas: &FrozenGas,
    n: usize,
    length: f64,
    area: f64,
    boundaries: Boundaries,
    initial: impl Fn(f64) -> [f64; 3],
) -> Euler1d {
    let states = (0..n)
        .map(|i| {
            let [rho, u, p] = initial(length * (i as f64 + 0.5) / n as f64);
            gas.conserved(rho, u, p).unwrap()
        })
        .collect();
    Euler1d::new(
        gas.clone(),
        Length::new(length / n as f64),
        area,
        states,
        boundaries,
        || false,
    )
    .unwrap()
}

// Independent direct integration, deliberately not network.totals() or a
// reduction over the reported balances. All branch energies include kinetic
// energy; the unresolved rigid plenum stores internal energy only.
fn integrated_storage(network: &GasNetwork1d) -> NetworkTotals {
    let mut mass = 0.0;
    let mut energy = 0.0;
    for branch in network.branches() {
        for (state, volume) in branch.states().iter().zip(branch.mesh().volumes()) {
            mass += state.density * volume;
            energy += state.energy * volume;
        }
    }
    for plenum in network.plenums() {
        mass += plenum.state().density * plenum.volume_m3();
        energy += plenum.state().energy * plenum.volume_m3();
    }
    NetworkTotals { mass, energy }
}

fn run_to(network: &mut GasNetwork1d, end: f64) {
    for _ in 0..20_000 {
        if network.time().value() >= end {
            return;
        }
        let report = network
            .step(
                Time::new(end - network.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        assert!(report.maximum_cfl <= 0.400_000_000_000_01);
    }
    panic!("network did not reach bounded test time");
}

#[test]
fn compact_plenum_relaxes_at_the_independent_acoustic_impedance_rate() {
    // R=1, rho=p=1, gamma=1.4. For a semi-infinite outgoing tube,
    // p'=rho*c*u. A rigid adiabatic plenum obeys p_dot=-(gamma*p/V) A*u,
    // hence p'(t)/p'(0)=exp(-c*A*t/V). Stagnation/static differences and
    // finite-amplitude thermodynamics enter only beyond linear order.
    // A stationary ghost passed to a two-sided Riemann solver incorrectly
    // doubles the acoustic port impedance and fails this physical oracle.
    let gamma: f64 = 1.4;
    let amplitude: f64 = 1e-4;
    let area = 0.1;
    let volume = 0.2;
    let end = 0.4;
    let exact = (-gamma.sqrt() * area * end / volume).exp();
    let mut errors = Vec::new();
    for n in [64, 128] {
        let gas = gas(0.0, 0.0);
        let branch = uniform_branch(
            &gas,
            n,
            1.0,
            area,
            Boundaries::Ends {
                left: Boundary::Transmissive,
                right: Boundary::Reflecting,
            },
            |_| [1.0, 0.0, 1.0],
        );
        let plenum = Plenum::new(
            volume,
            gas.conserved((1.0 + amplitude).powf(1.0 / gamma), 0.0, 1.0 + amplitude)
                .unwrap(),
        )
        .unwrap();
        let mut network = GasNetwork1d::new(
            vec![branch],
            vec![plenum],
            vec![PlenumConnection {
                branch: 0,
                end: BranchEnd::Left,
                plenum: 0,
            }],
            || false,
        )
        .unwrap();
        run_to(&mut network, end);
        // The distant rigid end cannot return a signal before 2L/c~=1.69 s.
        let pressure = gas
            .primitive(network.plenums()[0].state(), || false)
            .unwrap()
            .pressure;
        let normalized = (pressure - 1.0) / amplitude;
        errors.push((normalized - exact).abs());
        eprintln!(
            "compact-plenum acoustics: n={n}, normalized={normalized:.9e}, exact={exact:.9e}, error={:.9e}",
            errors[errors.len() - 1]
        );
    }
    // dx*A/V is 0.0078125 at 64 cells and halves with refinement. The
    // admitted 1% of initial excess also exceeds the O(1e-4) nonlinear term.
    assert!(
        errors[1] < 0.01,
        "incorrect compact-plenum impedance: {errors:?}"
    );
    assert!(
        errors[1] < 0.8 * errors[0] || errors[1] < 3e-4,
        "acoustic port error did not decrease toward the independent limit: {errors:?}"
    );
}

fn splitter(gas: &FrozenGas, plenum_pressure: f64) -> GasNetwork1d {
    let inlet = uniform_branch(gas, 96, 1.0, 0.1, walls(), |x| {
        let bump = if (0.35..0.55).contains(&x) {
            (core::f64::consts::PI * (x - 0.45) / 0.2).cos().powi(2)
        } else {
            0.0
        };
        let pressure = 1.0 + 0.4 * bump;
        // Independent right-running gamma-law simple-wave initial data. For
        // variable Cp this remains admitted initial data, not an exact wave.
        let density = pressure.powf(1.0 / 1.4);
        let velocity = 2.0 * 1.4_f64.sqrt() / 0.4 * (pressure.powf(0.4 / 2.8) - 1.0);
        [density, velocity, pressure]
    });
    let daughter = |area| uniform_branch(gas, 48, 0.7, area, walls(), |_| [1.0, 0.0, 1.0]);
    GasNetwork1d::new(
        vec![inlet, daughter(0.04), daughter(0.06)],
        vec![Plenum::new(0.002, gas.conserved(1.0, 0.0, plenum_pressure).unwrap()).unwrap()],
        vec![
            PlenumConnection {
                branch: 0,
                end: BranchEnd::Right,
                plenum: 0,
            },
            PlenumConnection {
                branch: 1,
                end: BranchEnd::Left,
                plenum: 0,
            },
            PlenumConnection {
                branch: 2,
                end: BranchEnd::Left,
                plenum: 0,
            },
        ],
        || false,
    )
    .unwrap()
}

#[test]
fn finite_pulse_reaches_both_branches_and_closes_actual_storage_and_port_exchange() {
    let gas = gas(0.0, 0.0);
    let mut network = splitter(&gas, 1.0);
    let initial = integrated_storage(&network);
    run_to(&mut network, 0.1);
    // The compact initial pulse ends at x=.55. Even its largest initial
    // characteristic cannot traverse the remaining .45 m in .1 seconds.
    for branch in &network.branches()[1..] {
        for state in branch.states() {
            let pressure = gas.primitive(*state, || false).unwrap().pressure;
            assert!(
                (pressure - 1.0).abs() < 1e-8,
                "spurious immediate branch response"
            );
        }
    }
    let mut peak = [1.0_f64; 2];
    let mut hottest_plenum = 1.0_f64;
    let end = 0.65;
    for _ in 0..20_000 {
        if network.time().value() >= end {
            break;
        }
        let old_plenum = network.plenums()[0].state();
        let old_branch: Vec<_> = network.branches().iter().map(Euler1d::totals).collect();
        let before = integrated_storage(&network);
        let report = network
            .step(
                Time::new(end - network.time().value()),
                StepControl::default(),
                || false,
            )
            .unwrap();
        let after = integrated_storage(&network);
        let ports = network.accepted_port_exchanges().unwrap();
        let branch_reports = network.accepted_branch_reports().unwrap();
        let new_plenum = network.plenums()[0].state();
        let volume = network.plenums()[0].volume_m3();
        let mass_from_ports: f64 = ports.iter().map(|p| p.mass_to_plenum).sum();
        let energy_from_ports: f64 = ports.iter().map(|p| p.energy_to_plenum).sum();
        assert!(
            (volume * (new_plenum.density - old_plenum.density) - mass_from_ports).abs() < 2e-13
        );
        assert!(
            (volume * (new_plenum.energy - old_plenum.energy) - energy_from_ports).abs() < 2e-12
        );
        assert_eq!(new_plenum.momentum, 0.0);
        for (index, branch) in network.branches().iter().enumerate() {
            let current = branch.totals();
            // The other end is an impermeable wall: exactly one material port
            // contributes to this branch's mass and total-energy changes.
            assert!(
                (current.mass - old_branch[index].mass + ports[index].mass_to_plenum).abs() < 2e-13
            );
            assert!(
                (current.energy - old_branch[index].energy + ports[index].energy_to_plenum).abs()
                    < 2e-12
            );
            assert!(
                (current.momentum
                    - old_branch[index].momentum
                    - branch_reports[index].boundary_change.momentum
                    - branch_reports[index].area_momentum_impulse)
                    .abs()
                    < 2e-12
            );
        }
        assert_eq!(report.boundary_change.mass, 0.0);
        assert_eq!(report.boundary_change.energy, 0.0);
        assert!((after.mass - before.mass).abs() < 2e-13);
        assert!((after.energy - before.energy).abs() < 2e-12);
        for (index, branch) in network.branches()[1..].iter().enumerate() {
            for state in branch.states() {
                peak[index] = peak[index].max(gas.primitive(*state, || false).unwrap().pressure);
            }
        }
        hottest_plenum =
            hottest_plenum.max(gas.primitive(new_plenum, || false).unwrap().temperature);
    }
    assert_eq!(network.time().value(), end);
    let final_state = integrated_storage(&network);
    assert!((final_state.mass - initial.mass).abs() < 2e-11);
    assert!((final_state.energy - initial.energy).abs() < 2e-10);
    assert!(
        peak.into_iter().all(|p| p > 1.01),
        "pulse did not enter both daughter ducts: {peak:?}"
    );
    assert!(
        hottest_plenum > 1.001,
        "arriving enthalpy did not heat the stored gas"
    );
    eprintln!("pulse splitter: daughter peaks={peak:?}, plenum Tmax={hottest_plenum:.9e}");
}

#[test]
fn variable_cp_network_and_every_accepted_exchange_ignore_formation_gauge() {
    let mut ordinary = splitter(&gas(0.4, 0.0), 1.3);
    let mut shifted = splitter(&gas(0.4, 1e18), 1.3);
    for _ in 0..70 {
        let a = ordinary
            .step(Time::new(0.002), StepControl::default(), || false)
            .unwrap();
        let b = shifted
            .step(Time::new(0.002), StepControl::default(), || false)
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(
            ordinary.accepted_port_exchanges(),
            shifted.accepted_port_exchanges()
        );
        assert_eq!(
            ordinary.accepted_branch_reports(),
            shifted.accepted_branch_reports()
        );
        for (a, b) in ordinary.branches().iter().zip(shifted.branches()) {
            assert_eq!(a.states(), b.states());
            assert_eq!(a.time(), b.time());
        }
        for (a, b) in ordinary.plenums().iter().zip(shifted.plenums()) {
            assert_eq!(a.state(), b.state());
        }
    }
}

fn mirrored_taper(gas: &FrozenGas, mirror: bool) -> GasNetwork1d {
    let n = 48;
    let positions: Vec<_> = (0..=n).map(|i| i as f64 / n as f64).collect();
    let areas: Vec<_> = positions.iter().map(|&x| 0.1 + 0.05 * x).collect();
    let volumes: Vec<_> = positions
        .windows(2)
        .map(|x| 0.1 * (x[1] - x[0]) + 0.025 * (x[1] * x[1] - x[0] * x[0]))
        .collect();
    let states: Vec<_> = (0..n)
        .map(|i| {
            let x = (i as f64 + 0.5) / n as f64;
            let bump = (-((x - 0.6) / 0.12).powi(2)).exp();
            gas.conserved(1.0, 0.2 * bump, 1.0 + 0.3 * bump).unwrap()
        })
        .collect();
    let (areas, volumes, states) = if mirror {
        (
            areas.into_iter().rev().collect(),
            volumes.into_iter().rev().collect(),
            states
                .into_iter()
                .rev()
                .map(|q| Conserved {
                    momentum: -q.momentum,
                    ..q
                })
                .collect(),
        )
    } else {
        (areas, volumes, states)
    };
    let branch = Euler1d::new_on_mesh(
        gas.clone(),
        Mesh1d::new(positions, areas, volumes).unwrap(),
        states,
        walls(),
        || false,
    )
    .unwrap();
    GasNetwork1d::new(
        vec![branch],
        vec![Plenum::new(0.03, gas.conserved(1.0, 0.0, 1.2).unwrap()).unwrap()],
        vec![PlenumConnection {
            branch: 0,
            end: if mirror {
                BranchEnd::Left
            } else {
                BranchEnd::Right
            },
            plenum: 0,
        }],
        || false,
    )
    .unwrap()
}

#[test]
fn reversing_a_taper_and_port_preserves_physics_and_reverses_axial_impulse() {
    let gas = gas(0.4, 0.0);
    let mut forward = mirrored_taper(&gas, false);
    let mut reversed = mirrored_taper(&gas, true);
    for _ in 0..60 {
        forward
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        reversed
            .step(Time::new(0.001), StepControl::default(), || false)
            .unwrap();
        let a = forward.accepted_port_exchanges().unwrap()[0];
        let b = reversed.accepted_port_exchanges().unwrap()[0];
        assert!((a.mass_to_plenum - b.mass_to_plenum).abs() < 2e-12);
        assert!((a.energy_to_plenum - b.energy_to_plenum).abs() < 2e-12);
        assert!((a.momentum_to_plenum + b.momentum_to_plenum).abs() < 2e-12);
        for (a, b) in forward.branches()[0]
            .states()
            .iter()
            .zip(reversed.branches()[0].states().iter().rev())
        {
            assert!((a.density - b.density).abs() < 2e-11);
            assert!((a.energy - b.energy).abs() < 2e-10);
            assert!((a.momentum + b.momentum).abs() < 2e-11);
        }
        let a = forward.plenums()[0].state();
        let b = reversed.plenums()[0].state();
        assert!((a.density - b.density).abs() < 2e-11);
        assert!((a.energy - b.energy).abs() < 2e-10);
    }
}
