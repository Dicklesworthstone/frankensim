use super::*;

fn surface(cells: &[([usize; 2], f64, [f64; 2])], force: [f64; 2]) -> FreeSurface {
    let mut grid = Grid::uniform(16, 8, 0.8);
    grid.flags.fill(Cell::Gas);
    grid.periodic_x = false;
    grid.periodic_y = false;
    for &(position, rho, velocity) in cells {
        let i = grid.idx(position[0], position[1]);
        grid.flags[i] = Cell::Interface;
        grid.fext[i] = force;
        grid.f[i] = equilibrium(rho, velocity[0], velocity[1]);
        crate::thermal::forcing::reframe2(&mut grid.f[i], velocity, [0.0; 2], force);
    }
    let mut sim = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    for &(position, _, _) in cells {
        let i = sim.grid.idx(position[0], position[1]);
        sim.mass[i] = 0.4 * sim.grid.moments(i).rho;
    }
    sim
}

fn inventory(sim: &FreeSurface, region: impl Fn(usize) -> bool) -> (f64, [f64; 2], f64) {
    let mut mass = 0.0;
    let mut momentum = [0.0; 2];
    let mut kinetic = 0.0;
    for i in 0..sim.mass.len() {
        if !region(i % sim.grid.nx) || !matches!(sim.grid.flags[i], Cell::Fluid | Cell::Interface) {
            continue;
        }
        let state = sim.grid.moments(i);
        let m = if sim.grid.flags[i] == Cell::Fluid {
            state.rho
        } else {
            sim.mass[i]
        };
        mass += m;
        for d in 0..2 {
            momentum[d] += m * state.u[d];
        }
        kinetic += 0.5 * m * norm2(state.u);
    }
    (mass, momentum, kinetic)
}

fn assert_balance(
    sim: &FreeSurface,
    before: (f64, [f64; 2], f64),
    report: ConversionMomentumReport2,
) {
    let after = inventory(sim, |_| true);
    let scale = before.0.abs().max(after.0.abs()).max(1.0);
    assert!((after.0 - before.0).abs() < 2e-13 * scale);
    for d in 0..2 {
        assert!((after.1[d] - before.1[d]).abs() < 2e-14 * scale);
    }
    assert!((after.2 - before.2 - report.kinetic_energy_change).abs() < 2e-15 * scale);
    assert!(report.absolute_mass_defect < 2e-13 * scale);
    assert!(
        report
            .absolute_momentum_defect
            .iter()
            .all(|x| *x < 2e-14 * scale)
    );
    assert_eq!(sim.carry, 0.0);
    for (i, flag) in sim.grid.flags.iter().enumerate() {
        if *flag == Cell::Interface {
            assert!(sim.mass[i] >= 0.0 && sim.mass[i] <= sim.grid.moments(i).rho);
        }
    }
}

#[test]
fn g1_advancing_front_carries_density_weighted_momentum_under_boost_and_force() {
    let gate = CancelGate::new_clock_free();
    for scale in [0.4, 1.0, 2.5] {
        for boost in [[0.0; 2], [0.02, -0.01]] {
            let incoming = [0.06 + boost[0], 0.01 + boost[1]];
            let receiving = [-0.025 + boost[0], 0.03 + boost[1]];
            let mut sim = surface(
                &[
                    ([3, 3], 0.8 * scale, incoming),
                    ([4, 3], 1.4 * scale, receiving),
                ],
                [0.0002 * scale, -0.0001 * scale],
            );
            let source = sim.grid.idx(3, 3);
            let receiver = sim.grid.idx(4, 3);
            let density_bits = [source, receiver].map(|i| sim.grid.moments(i).rho.to_bits());
            sim.mass[source] = 1.2 * sim.grid.moments(source).rho;
            let before = inventory(&sim, |_| true);
            let report = sim.convert_momentum(&gate).unwrap();
            assert_balance(&sim, before, report);
            assert!(report.routes > 1 && report.kinetic_energy_change < 0.0);
            assert_eq!(
                [source, receiver].map(|i| sim.grid.moments(i).rho.to_bits()),
                density_bits
            );
            assert!((sim.grid.moments(receiver).u[0] - receiving[0]).abs() > 1e-6);
            assert_eq!(sim.grid.flags[source], Cell::Fluid);
        }
    }
}

#[test]
fn g1_receding_front_carries_signed_debt_without_borrowing_from_dry_ghosts() {
    let gate = CancelGate::new_clock_free();
    let mut sim = surface(
        &[
            ([3, 3], 0.8, [-0.06, 0.02]),
            ([4, 3], 1.7, [0.02, -0.01]),
            ([3, 4], 0.5, [-0.01, 0.04]),
        ],
        [0.0002, -0.0001],
    );
    let source = sim.grid.idx(3, 3);
    let recipient = sim.grid.idx(4, 3);
    let dry = sim.grid.idx(3, 4);
    sim.mass[source] = -0.02;
    sim.mass[dry] = 0.0;
    let before = inventory(&sim, |_| true);
    let old = sim.grid.moments(recipient);
    let donor = sim.grid.moments(source);
    let old_mass = sim.mass[recipient];
    let report = sim.convert_momentum(&gate).unwrap();
    assert_balance(&sim, before, report);
    assert_eq!(sim.grid.flags[source], Cell::Gas);
    assert_eq!(sim.mass[dry], 0.0);
    for d in 0..2 {
        let expected = (old_mass * old.u[d] - 0.02 * donor.u[d]) / (old_mass - 0.02);
        assert!((sim.grid.moments(recipient).u[d] - expected).abs() < 2e-15);
    }
    assert!(report.kinetic_energy_change > 0.0);
}

#[test]
fn g1_disconnected_fronts_retain_their_own_momentum() {
    let gate = CancelGate::new_clock_free();
    let mut sim = surface(
        &[
            ([3, 3], 0.8, [0.05, -0.01]),
            ([4, 3], 1.2, [-0.02, 0.03]),
            ([5, 3], 1.6, [-0.04, -0.02]),
            ([11, 3], 0.6, [0.01, -0.03]),
            ([12, 3], 1.1, [0.02, 0.04]),
        ],
        [0.0001, -0.0001],
    );
    for (x, excess) in [(3, 0.12), (5, 0.18), (11, 0.1)] {
        let i = sim.grid.idx(x, 3);
        sim.mass[i] = sim.grid.moments(i).rho + excess;
    }
    let before = inventory(&sim, |_| true);
    let left = inventory(&sim, |x| x < 8);
    let right = inventory(&sim, |x| x >= 8);
    let report = sim.convert_momentum(&gate).unwrap();
    assert_balance(&sim, before, report);
    for (actual, expected) in [
        (inventory(&sim, |x| x < 8), left),
        (inventory(&sim, |x| x >= 8), right),
    ] {
        assert!((actual.0 - expected.0).abs() < 2e-13);
        for d in 0..2 {
            assert!((actual.1[d] - expected.1[d]).abs() < 2e-14);
        }
    }
}

#[test]
fn g0_small_periodic_grid_deduplicates_recipient_capacity() {
    let gate = CancelGate::new_clock_free();
    let mut grid = Grid::uniform(2, 1, 0.8);
    grid.periodic_x = true;
    grid.periodic_y = true;
    grid.flags.fill(Cell::Interface);
    grid.f[0] = equilibrium(1.0, 0.05, -0.02);
    grid.f[1] = equilibrium(1.2, -0.01, 0.03);
    let mut sim = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    sim.mass[0] = sim.grid.moments(0).rho + 0.1;
    let before = inventory(&sim, |_| true);
    let report = sim.convert_momentum(&gate).unwrap();
    assert_balance(&sim, before, report);
    assert_eq!(report.routes, 1);
}

#[test]
fn g1_velocity_retarget_preserves_density_and_physical_viscous_stress() {
    let previous = [0.02, -0.01];
    let next = [-0.01, 0.025];
    let force = [0.0003, -0.0002];
    for rho in [0.4, 0.8, 1.3, 2.5] {
        let mut f = equilibrium(rho, previous[0], previous[1]);
        for q in 0..Q {
            f[q] += 0.001 * rho * crate::W[q] * f64::from(E[q].0 * E[q].1);
        }
        crate::thermal::forcing::reframe2(&mut f, previous, [0.0; 2], force);
        let density = f.iter().sum::<f64>();
        let after = retarget_velocity(f, previous, next, force).unwrap();
        assert_eq!(density.to_bits(), after.iter().sum::<f64>().to_bits());
        for a in 0..2 {
            let directions = E.map(|(x, y)| [f64::from(x), f64::from(y)]);
            let momentum = after
                .iter()
                .zip(directions)
                .map(|(f, e)| f * e[a])
                .sum::<f64>()
                + 0.5 * force[a];
            assert!((momentum - density * next[a]).abs() < 2e-14 * rho);
            for b in 0..2 {
                let moment = |populations: [f64; Q], u: [f64; 2]| {
                    populations
                        .iter()
                        .zip(directions)
                        .map(|(f, e)| f * e[a] * e[b])
                        .sum::<f64>()
                        + 0.5 * (u[a] * force[b] + force[a] * u[b])
                        - density * u[a] * u[b]
                };
                assert!((moment(after, next) - moment(f, previous)).abs() < 2e-14 * rho);
            }
        }
        assert_eq!(retarget_velocity(f, previous, previous, force).unwrap(), f);
    }
}

#[test]
fn g3_native_step_replays_unchanged_without_conversion_and_refuses_atomically() {
    let gate = CancelGate::new_clock_free();
    let mut grid = Grid::uniform(8, 8, 0.8);
    grid.f.fill(equilibrium(1.2, 0.02, -0.01));
    let mut sim = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    let mut legacy = sim.clone();
    legacy.step();
    assert_eq!(
        sim.step_with_momentum(&gate).unwrap(),
        ConversionMomentumReport2::default()
    );
    assert_eq!(sim.grid.f, legacy.grid.f);
    assert_eq!(sim.mass, legacy.mass);

    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    let bits = sim.grid.f.clone();
    let mass = sim.mass.clone();
    assert_eq!(
        sim.step_with_momentum(&cancelled),
        Err(FreeSurfaceMomentumError2::Cancelled)
    );
    assert_eq!(sim.grid.f, bits);
    assert_eq!(sim.mass, mass);

    let mut grid = Grid::uniform(3, 3, 0.8);
    grid.flags.fill(Cell::Wall);
    grid.flags[4] = Cell::Interface;
    let mut isolated = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    isolated.mass[4] = 2.0;
    let bits = isolated.grid.f.clone();
    let mass = isolated.mass.clone();
    let flags = isolated.grid.flags.clone();
    assert!(matches!(
        isolated.step_with_momentum(&gate),
        Err(FreeSurfaceMomentumError2::Invalid(_))
    ));
    assert_eq!(isolated.grid.f, bits);
    assert_eq!(isolated.mass, mass);
    assert_eq!(isolated.grid.flags, flags);
    assert_eq!(isolated.conversions.to_fluid, 0);
}

#[test]
fn g0_newly_wet_cells_require_valid_material_state_before_publication() {
    let gate = CancelGate::new_clock_free();
    for bad_force in [true, false] {
        let mut sim = surface(&[([3, 3], 1.0, [0.01, 0.0])], [0.0; 2]);
        let source = sim.grid.idx(3, 3);
        let dry = sim.grid.idx(4, 3);
        sim.mass[source] = 1.2 * sim.grid.moments(source).rho;
        if bad_force {
            sim.grid.fext[dry][0] = f64::NAN;
        } else {
            sim.grid.tau[dry] = 0.5;
        }
        let bits = sim.grid.f.clone();
        let mass = sim.mass.clone();
        let flags = sim.grid.flags.clone();
        assert!(matches!(
            sim.step_with_momentum(&gate),
            Err(FreeSurfaceMomentumError2::Invalid(_))
        ));
        assert_eq!(sim.grid.f, bits);
        assert_eq!(sim.mass, mass);
        assert_eq!(sim.grid.flags, flags);
        assert_eq!(sim.conversions.gas_to_interface, 0);
    }
    assert!(admit_initialized_cell(&[f64::INFINITY; Q], [0.0; 2]).is_err());
}
