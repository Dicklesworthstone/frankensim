//! G1/G3 local signed inventories and G0/G4 numerical-domain refusals.
use super::*;

type CellState = ([usize; 2], f64, f64, [f64; 2]);

fn surface(cells: &[CellState]) -> FreeSurface {
    let mut grid = Grid::uniform(18, 8, 0.8);
    grid.flags.fill(Cell::Gas);
    grid.periodic_x = false;
    grid.periodic_y = false;
    grid.g = [0.00007, -0.00011];
    for &(p, rho, _, velocity) in cells {
        let i = grid.idx(p[0], p[1]);
        grid.flags[i] = Cell::Interface;
        grid.fext[i] = [0.0002, -0.0001];
        grid.f[i] = equilibrium(rho, velocity[0], velocity[1]);
        let force = core::array::from_fn(|d| rho * grid.g[d] + grid.fext[i][d]);
        crate::thermal::forcing::reframe2(&mut grid.f[i], velocity, [0.0; 2], force);
    }
    let mut flow = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    for &(p, _, mass, _) in cells {
        let i = flow.grid.idx(p[0], p[1]);
        flow.mass[i] = mass;
    }
    flow
}

#[derive(Clone, Copy, Debug, Default)]
struct Inventory {
    mass: f64,
    momentum: [f64; 2],
    angular: f64,
    kinetic: f64,
}

fn inventory(flow: &FreeSurface, right: bool) -> Inventory {
    let mut total = Inventory::default();
    for i in 0..flow.mass.len() {
        let x = i % flow.grid.nx;
        if (x >= 8) != right || !matches!(flow.grid.flags[i], Cell::Fluid | Cell::Interface) {
            continue;
        }
        let rho = flow.grid.f[i].iter().sum::<f64>();
        // Independent physical readback, not the conservation report or its
        // local-origin arithmetic. Angular origin is the first domain centre.
        let velocity: [f64; 2] = core::array::from_fn(|d| {
            let raw = flow.grid.f[i]
                .iter()
                .zip(E)
                .map(|(f, e)| f * f64::from([e.0, e.1][d]))
                .sum::<f64>();
            (raw + 0.5 * flow.grid.g[d].mul_add(rho, flow.grid.fext[i][d])) / rho
        });
        let mass = if flow.grid.flags[i] == Cell::Fluid {
            rho
        } else {
            flow.mass[i]
        };
        let p = velocity.map(|u| mass * u);
        total.mass += mass;
        total.angular += x as f64 * p[1] - (i / flow.grid.nx) as f64 * p[0];
        for d in 0..2 {
            total.momentum[d] += p[d];
            total.kinetic += 0.5 * mass * velocity[d] * velocity[d];
        }
    }
    total
}

#[test]
fn g1_g3_angular_conversion_keeps_retained_source_debt_and_detached_fronts_local() {
    let gate = CancelGate::new_clock_free();
    for boost in [[0.0; 2], [0.02, -0.01]] {
        let mut cells = [
            ([3, 3], 0.8, 0.96, [0.04, 0.02]),
            ([4, 3], 1.2, -0.02, [-0.03, 0.01]),
            ([4, 4], 0.9, 0.36, [0.01, -0.02]),
            ([12, 3], 1.1, 1.22, [-0.035, 0.015]),
            ([13, 3], 1.2, 0.48, [0.025, -0.015]),
            ([12, 4], 0.7, 0.28, [0.005, 0.02]),
        ];
        for cell in &mut cells {
            for d in 0..2 {
                cell.3[d] += boost[d];
            }
        }
        let mut flow = surface(&cells);
        let before = [inventory(&flow, false), inventory(&flow, true)];
        let density_bits = cells.map(|(p, _, _, _)| {
            flow.grid.f[flow.grid.idx(p[0], p[1])]
                .iter()
                .sum::<f64>()
                .to_bits()
        });
        let mut linear = flow.clone();
        let linear_report = linear.convert_momentum(&gate).unwrap();
        assert_eq!(linear_report.absolute_angular_momentum_defect, None);
        assert!(
            [false, true].into_iter().enumerate().any(|(i, region)| {
                (inventory(&linear, region).angular - before[i].angular).abs() > 1e-6
            }),
            "the control must actually expose conversion torque"
        );

        flow.set_conversion_redistribution(ConversionRedistribution2::ConservativeLocalAngular)
            .unwrap();
        let report = flow.convert_momentum(&gate).unwrap();
        let after = [inventory(&flow, false), inventory(&flow, true)];
        for i in 0..2 {
            assert!((after[i].mass - before[i].mass).abs() < 3e-13);
            for d in 0..2 {
                assert!((after[i].momentum[d] - before[i].momentum[d]).abs() < 3e-14);
            }
            assert!((after[i].angular - before[i].angular).abs() < 5e-13);
        }
        let kinetic_change =
            (after[0].kinetic - before[0].kinetic) + (after[1].kinetic - before[1].kinetic);
        assert!((kinetic_change - report.kinetic_energy_change).abs() < 5e-15);
        assert!(report.routes > 4 && report.absolute_angular_momentum_defect.unwrap() < 5e-13);
        assert_eq!(flow.carry, 0.0);
        for (index, (p, _, _, _)) in cells.into_iter().enumerate() {
            assert_eq!(
                flow.grid.f[flow.grid.idx(p[0], p[1])]
                    .iter()
                    .sum::<f64>()
                    .to_bits(),
                density_bits[index]
            );
        }
        let retained = flow.grid.idx(3, 3);
        assert_eq!(flow.grid.flags[retained], Cell::Interface);
        assert!(
            flow.mass[retained] < flow.grid.f[retained].iter().sum::<f64>(),
            "the original overfilled source must later pay its neighbor's debt"
        );
        assert_eq!(flow.grid.flags[flow.grid.idx(4, 3)], Cell::Gas);
    }
}

#[test]
fn g0_g4_angular_policy_admits_identity_and_refuses_periodicity_or_unresolved_fast_spin() {
    let gate = CancelGate::new_clock_free();
    let mut identity = surface(&[([3, 3], 1.0, 0.4, [0.02, -0.01])]);
    identity.grid.periodic_x = true;
    assert!(
        identity
            .set_conversion_redistribution(ConversionRedistribution2::ConservativeLocalAngular)
            .is_err()
    );
    assert_eq!(
        identity.conversion_redistribution(),
        ConversionRedistribution2::ConservativeLocal
    );
    identity.grid.periodic_x = false;
    let mut explicit = identity.clone();
    let explicit_report = explicit.step_with_angular_momentum(&gate).unwrap();
    assert_eq!(explicit_report.absolute_angular_momentum_defect, Some(0.0));
    assert_eq!(
        explicit.conversion_redistribution(),
        ConversionRedistribution2::ConservativeLocal
    );
    identity
        .set_conversion_redistribution(ConversionRedistribution2::ConservativeLocalAngular)
        .unwrap();
    let populations = identity.grid.f.clone();
    let report = identity.convert_momentum(&gate).unwrap();
    assert_eq!(report.absolute_angular_momentum_defect, Some(0.0));
    assert_eq!(report.routes, 0);
    assert_eq!(identity.grid.f, populations);
    identity.grid.periodic_y = true; // Public geometry may change after selection.
    let before = identity.clone();
    assert!(matches!(
        identity.step_with_momentum(&gate),
        Err(FreeSurfaceMomentumError2::Invalid(
            "angular conversion requires a nonperiodic domain"
        ))
    ));
    assert_eq!(identity.grid.f, before.grid.f);
    assert_eq!(identity.mass, before.mass);
    assert_eq!(identity.grid.flags, before.grid.flags);
    assert_eq!(identity.conversions.to_fluid, before.conversions.to_fluid);
    assert_eq!(identity.conversions.to_gas, before.conversions.to_gas);

    let mut spin = surface(&[
        ([3, 3], 1.0, -0.01, [0.0, 0.1]),
        ([4, 3], 1.0, 0.4, [0.0; 2]),
        ([4, 4], 1.0, 1e-8, [0.0; 2]),
    ]);
    // Linear conversion is admissible, but preserving the lost moment arm
    // would require a huge velocity at the tiny remaining off-axis inventory.
    let mut linear = spin.clone();
    linear.convert_momentum(&gate).unwrap();
    spin.set_conversion_redistribution(ConversionRedistribution2::ConservativeLocalAngular)
        .unwrap();
    assert!(matches!(
        spin.convert_momentum(&gate),
        Err(FreeSurfaceMomentumError2::Invalid(
            "conversion exceeds the finite low-Mach lattice domain; refine"
        ))
    ));
}
