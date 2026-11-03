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

#[test]
fn g3_material_capillary_wall_trajectories_admit_conservative_steps() {
    let gate = CancelGate::new_clock_free();
    for wall in [
        PlanarWall2::XMin,
        PlanarWall2::XMax,
        PlanarWall2::YMin,
        PlanarWall2::YMax,
    ] {
        for angle in [
            core::f64::consts::PI / 3.0,
            core::f64::consts::FRAC_PI_2,
            2.0 * core::f64::consts::PI / 3.0,
        ] {
            // Same rho=1000 kg/m³, mu=0.1 Pa s, dx=1 mm, dt=0.2 ms
            // and sigma=0.05 N/m as the sourced material wetting battery.
            let mut grid = Grid::uniform(20, 20, 0.56);
            grid.periodic_x = false;
            grid.periodic_y = false;
            for y in 0..20 {
                for x in 0..20 {
                    let (normal, tangent) = match wall {
                        PlanarWall2::XMin => (x, y),
                        PlanarWall2::XMax => (19 - x, y),
                        PlanarWall2::YMin => (y, x),
                        PlanarWall2::YMax => (19 - y, x),
                    };
                    let i = grid.idx(x, y);
                    grid.flags[i] = if x == 0 || y == 0 || x == 19 || y == 19 {
                        Cell::Wall
                    } else if (normal * normal) as f64 + (tangent as f64 - 10.0).powi(2) < 36.0 {
                        Cell::Fluid
                    } else {
                        Cell::Gas
                    };
                }
            }
            let mut sim = FreeSurface::new(grid, 0.002, ContactModel::Neutral);
            sim.set_planar_contact_angle(wall, angle).unwrap();
            let initial_mass = sim.ledger_mass();
            for _ in 0..12 {
                sim.step_with_momentum(&gate).unwrap();
            }
            assert!((sim.ledger_mass() - initial_mass).abs() < 1e-10 * initial_mass);
            assert_eq!(sim.carry, 0.0);
        }
    }
}

#[test]
fn g3_material_gravity_advects_a_drop_through_bounded_momentum_conversions() {
    // Exact physical scaling of the material dispatch regression: nu=1e-4
    // m²/s, dx=1 mm, dt=0.4 ms and acceleration -9.81 m/s².
    let mut grid = Grid::uniform(20, 20, 0.62);
    grid.periodic_x = false;
    grid.periodic_y = false;
    for y in 0..20 {
        for x in 0..20 {
            let i = grid.idx(x, y);
            grid.flags[i] = if x == 0 || y == 0 || x == 19 || y == 19 {
                Cell::Wall
            } else if (x as f64 - 10.0).powi(2) + (y as f64 - 12.0).powi(2) < 16.0 {
                Cell::Fluid
            } else {
                Cell::Gas
            };
        }
    }
    grid.rebind_gravity([0.0, -9.81 * 0.0004_f64.powi(2) / 0.001], &mut Vec::new())
        .unwrap();
    let mut sim = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    let mass = sim.ledger_mass();
    assert!(inventory(&sim, |_| true).1[1].abs() < 2e-14);
    let mut routes = 0;
    for _ in 0..41 {
        routes += sim
            .step_with_momentum(&CancelGate::new_clock_free())
            .unwrap()
            .routes;
    }
    assert!(routes > 0);
    assert!(inventory(&sim, |_| true).1[1] < 0.0);
    assert!((sim.ledger_mass() - mass).abs() < 1e-11 * mass);
    assert_eq!(sim.carry, 0.0);
}

// Independent signed cell-centre inventory. A different measurement origin
// must not alter closure when the linear momentum also closes.
fn angular_inventory(sim: &FreeSurface, origin: [f64; 2], region: impl Fn(usize) -> bool) -> f64 {
    let mut angular = 0.0;
    for (i, flag) in sim.grid.flags.iter().enumerate() {
        if !region(i % sim.grid.nx) || !matches!(flag, Cell::Fluid | Cell::Interface) {
            continue;
        }
        let state = sim.grid.moments(i);
        let mass = if *flag == Cell::Fluid {
            state.rho
        } else {
            sim.mass[i]
        };
        let x = (i % sim.grid.nx) as f64 - origin[0];
        let y = (i / sim.grid.nx) as f64 - origin[1];
        angular += mass * (x * state.u[1] - y * state.u[0]);
    }
    angular
}

#[test]
fn g1_angular_conversion_closes_off_axis_fronts_in_each_selected_frame() {
    let gate = CancelGate::new_clock_free();
    for density_scale in [0.4, 1.0, 2.5] {
        // Each boost defines its own conservation target. Relocated first
        // mass moments do not permit a general boost-equivariance claim.
        for boost in [[0.0; 2], [0.02, -0.01]] {
            let mut sim = surface(
                &[
                    (
                        [3, 3],
                        0.8 * density_scale,
                        [0.06 + boost[0], 0.01 + boost[1]],
                    ),
                    (
                        [4, 4],
                        1.4 * density_scale,
                        [-0.02 + boost[0], 0.03 + boost[1]],
                    ),
                ],
                [0.0002 * density_scale, -0.0001 * density_scale],
            );
            let source = sim.grid.idx(3, 3);
            sim.mass[source] = 1.2 * sim.grid.moments(source).rho;
            let before = inventory(&sim, |_| true);
            let origins = [[0.0, 0.0], [13.5, -7.25]];
            let angular = origins.map(|origin| angular_inventory(&sim, origin, |_| true));
            let mut linear = sim.clone();
            let linear_report = linear.convert_momentum(&gate).unwrap();
            assert_eq!(linear_report.absolute_angular_momentum_defect, None);
            assert!(
                (angular_inventory(&linear, origins[0], |_| true) - angular[0]).abs()
                    > 1e-6 * density_scale
            );
            let report = sim.convert_angular_momentum(&gate).unwrap();
            assert_balance(&sim, before, report);
            for (origin, expected) in origins.into_iter().zip(angular) {
                assert!(
                    (angular_inventory(&sim, origin, |_| true) - expected).abs()
                        < 5e-13 * density_scale
                );
            }
            assert!(report.absolute_angular_momentum_defect.unwrap() < 5e-13 * density_scale);
        }
    }
}

#[test]
fn g1_angular_conversion_keeps_disconnected_debts_and_dual_role_sources_local() {
    let gate = CancelGate::new_clock_free();
    let mut sim = surface(
        &[
            ([3, 3], 0.8, [0.05, -0.01]),
            ([4, 3], 1.2, [-0.02, 0.03]),
            ([4, 4], 0.7, [0.01, 0.04]),
            ([11, 3], 0.6, [0.03, -0.04]),
            ([12, 3], 1.1, [-0.03, 0.01]),
            ([11, 4], 0.8, [0.01, 0.03]),
        ],
        [0.0001, -0.0001],
    );
    // The positive source becomes interface again beside the signed debt
    // and receives a transfer. Its original overshoot must not be forgotten.
    let positive = sim.grid.idx(4, 3);
    sim.mass[positive] = sim.grid.moments(positive).rho + 0.12;
    let left_debt = sim.grid.idx(3, 3);
    let right_debt = sim.grid.idx(11, 3);
    sim.mass[left_debt] = -0.02;
    sim.mass[right_debt] = -0.03;
    let before = inventory(&sim, |_| true);
    let left = inventory(&sim, |x| x < 8);
    let right = inventory(&sim, |x| x >= 8);
    let angular = [
        angular_inventory(&sim, [3.0, 3.0], |x| x < 8),
        angular_inventory(&sim, [11.0, 3.0], |x| x >= 8),
    ];
    let report = sim.convert_angular_momentum(&gate).unwrap();
    assert_balance(&sim, before, report);
    for (actual, expected) in [
        (inventory(&sim, |x| x < 8), left),
        (inventory(&sim, |x| x >= 8), right),
    ] {
        assert!((actual.0 - expected.0).abs() < 3e-13);
        for d in 0..2 {
            assert!((actual.1[d] - expected.1[d]).abs() < 3e-14);
        }
    }
    assert!((angular_inventory(&sim, [3.0, 3.0], |x| x < 8) - angular[0]).abs() < 3e-14);
    assert!((angular_inventory(&sim, [11.0, 3.0], |x| x >= 8) - angular[1]).abs() < 3e-14);
    assert_eq!(sim.grid.flags[positive], Cell::Interface);
}

#[test]
fn g3_angular_conversion_is_coordinate_translation_invariant() {
    let gate = CancelGate::new_clock_free();
    let make = |offset| {
        let mut sim = surface(
            &[
                ([3 + offset, 3], 0.8, [0.06, 0.01]),
                ([4 + offset, 4], 1.4, [-0.02, 0.03]),
            ],
            [0.0002, -0.0001],
        );
        let source = sim.grid.idx(3 + offset, 3);
        sim.mass[source] = 1.2 * sim.grid.moments(source).rho;
        sim
    };
    let mut original = make(0);
    let mut shifted = make(7);
    let a = original.convert_angular_momentum(&gate).unwrap();
    let b = shifted.convert_angular_momentum(&gate).unwrap();
    assert_eq!(a, b);
    for y in 2..=4 {
        for x in 2..=5 {
            let i = original.grid.idx(x, y);
            let j = shifted.grid.idx(x + 7, y);
            assert_eq!(original.grid.flags[i], shifted.grid.flags[j]);
            assert_eq!(original.mass[i].to_bits(), shifted.mass[j].to_bits());
            assert_eq!(
                original.grid.f[i].map(f64::to_bits),
                shifted.grid.f[j].map(f64::to_bits)
            );
        }
    }
}

#[test]
fn g4_public_angular_step_refuses_unresolved_spin_periodicity_and_cancellation_atomically() {
    let gate = CancelGate::new_clock_free();
    let mut grid = Grid::uniform(4, 3, 0.8);
    grid.periodic_x = false;
    grid.periodic_y = false;
    grid.flags.fill(Cell::Wall);
    let source = grid.idx(1, 1);
    let recipient = grid.idx(2, 1);
    grid.flags[source] = Cell::Interface;
    grid.flags[recipient] = Cell::Interface;
    grid.f[source] = equilibrium(1.0, 0.0, 0.04);
    let mut sim = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    sim.mass[source] = -0.02;
    let before = sim.clone();
    let mut linear = sim.clone();
    linear.step_with_momentum(&gate).unwrap();
    assert_eq!(
        sim.step_with_angular_momentum(&gate),
        Err(invalid(
            "converted cells cannot resolve angular momentum; refine"
        )),
    );
    assert_eq!(sim.grid.f, before.grid.f);
    assert_eq!(sim.mass, before.mass);
    assert_eq!(sim.grid.flags, before.grid.flags);
    assert_eq!(sim.conversions.to_gas, before.conversions.to_gas);

    let mut grid = Grid::uniform(8, 8, 0.8);
    grid.periodic_x = false;
    grid.periodic_y = false;
    let mut no_conversion = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    let mut legacy = no_conversion.clone();
    legacy.step();
    let report = no_conversion.step_with_angular_momentum(&gate).unwrap();
    assert_eq!(report.absolute_angular_momentum_defect, Some(0.0));
    assert_eq!(no_conversion.grid.f, legacy.grid.f);
    assert_eq!(no_conversion.mass, legacy.mass);
    let before = no_conversion.clone();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        no_conversion.step_with_angular_momentum(&cancelled),
        Err(FreeSurfaceMomentumError2::Cancelled)
    );
    assert_eq!(no_conversion.grid.f, before.grid.f);
    assert_eq!(no_conversion.mass, before.mass);
    no_conversion.grid.periodic_x = true;
    assert!(matches!(
        no_conversion.step_with_angular_momentum(&gate),
        Err(FreeSurfaceMomentumError2::Invalid(_))
    ));
    assert_eq!(no_conversion.grid.f, before.grid.f);
    assert_eq!(no_conversion.mass, before.mass);
}
