//! G1 pressure/traction and G3/G4 conservation, replay and atomic refusal.
use super::*;

fn slab(rho: f64) -> FreeSurface {
    let mut grid = Grid::uniform(12, 12, 0.8);
    grid.periodic_y = false;
    for y in 0..12 {
        for x in 0..12 {
            let i = grid.idx(x, y);
            grid.f[i] = equilibrium(rho, 0.0, 0.0);
            grid.flags[i] = if y == 0 || y == 11 {
                Cell::Wall
            } else if y < 7 {
                Cell::Fluid
            } else {
                Cell::Gas
            };
        }
    }
    FreeSurface::new(grid, 0.0, ContactModel::Neutral)
}

fn same(a: &FreeSurface, b: &FreeSurface) {
    assert_eq!(a.grid.f, b.grid.f);
    assert_eq!(a.grid.flags, b.grid.flags);
    assert_eq!(a.mass, b.mass);
    assert_eq!(a.carry, b.carry);
    assert_eq!(a.conversions.to_fluid, b.conversions.to_fluid);
    assert_eq!(a.conversions.to_gas, b.conversions.to_gas);
    assert_eq!(
        a.conversions.gas_to_interface,
        b.conversions.gas_to_interface
    );
    assert_eq!(
        a.conversions.fluid_to_interface,
        b.conversions.fluid_to_interface
    );
}

#[test]
fn g1_loaded_flat_surface_remains_at_rest_and_keeps_mass() {
    let gate = CancelGate::new_clock_free();
    for rho in [0.985, 1.015] {
        let mut flow = slab(rho);
        let mass = flow.ledger_mass();
        for _ in 0..10 {
            flow.step_with_gas_pressure(&gate, 0.05, |_| Ok((rho - 1.0) * CS2))
                .unwrap();
            assert!((flow.ledger_mass() - mass).abs() < 2e-12);
            for i in 0..flow.grid.f.len() {
                if matches!(flow.grid.flags[i], Cell::Fluid | Cell::Interface) {
                    let state = flow.grid.moments(i);
                    assert!((state.rho - rho).abs() < 3e-14);
                    assert!(state.u.iter().all(|u| u.abs() < 2e-14));
                }
            }
        }
    }
}

#[test]
fn g1_distinct_neighbor_loads_enter_the_actual_incoming_populations() {
    let mut grid = Grid::uniform(5, 5, 0.8);
    grid.flags.fill(Cell::Gas);
    let cell = grid.idx(2, 2);
    grid.flags[cell] = Cell::Fluid;
    let mut flow = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    let gate = CancelGate::new_clock_free();
    let mass = flow.ledger_mass();
    let visited = std::cell::RefCell::new(Vec::new());
    flow.step_with_gas_link_pressure(&gate, 0.05, |link| {
        assert_eq!(link.liquid_cell, [2, 2]);
        let q = link.incoming_direction;
        assert_eq!(
            link.gas_cell,
            [(2 - E[q].0) as usize, (2 - E[q].1) as usize]
        );
        visited.borrow_mut().push(q);
        Ok(CS2 * 0.001 * q as f64)
    })
    .unwrap();
    assert_eq!(*visited.borrow(), (1..Q).collect::<Vec<_>>());
    for q in 1..Q {
        // Static anti-population traction: f_in = 2*w*rho_boundary - f_out.
        let expected = crate::W[q] * (1.0 + 0.002 * q as f64);
        assert!((flow.grid.f[cell][q] - expected).abs() < 2e-16);
    }
    assert_eq!(flow.ledger_mass(), mass);
}

#[test]
fn g1_pressure_reconstruction_uses_the_physical_forced_stress_frame() {
    let velocity = [0.08, -0.03];
    let force = [0.004, -0.007];
    let rho = 1.025;
    let mut cache = GasEquilibrium2::default();
    let eq = cache
        .at_link(
            || rho,
            velocity,
            force,
            None,
            GasPressureLink2 {
                liquid_cell: [1, 1],
                gas_cell: [0, 1],
                incoming_direction: 1,
            },
        )
        .unwrap();
    let u2 = velocity[0] * velocity[0] + velocity[1] * velocity[1];
    let uf = velocity[0] * force[0] + velocity[1] * force[1];
    for q in 1..Q {
        let eu = f64::from(E[q].0) * velocity[0] + f64::from(E[q].1) * velocity[1];
        let ef = f64::from(E[q].0) * force[0] + f64::from(E[q].1) * force[1];
        // Independent paired stress from the second-order Hermite source:
        // eq_q+eq_opp - (S_q+S_opp)/2. Odd momentum terms cancel.
        let expected = 2.0 * crate::W[q] * rho * (1.0 + 4.5 * eu * eu - 1.5 * u2)
            + crate::W[q] * (3.0 * uf - 9.0 * eu * ef);
        assert!((eq[q] + eq[OPP[q]] - expected).abs() < 2e-16);
    }
}

#[test]
fn g3_zero_cell_and_link_pressure_replay_the_original_step() {
    let gate = CancelGate::new_clock_free();
    let mut native = slab(1.0);
    native.grid.g = [0.0001, -0.0002];
    let mut cell = native.clone();
    let mut link = native.clone();
    for _ in 0..5 {
        native.step();
        cell.step_with_gas_pressure(&gate, 0.1, |_| Ok(0.0))
            .unwrap();
        link.step_with_gas_link_pressure(&gate, 0.1, |_| Ok(0.0))
            .unwrap();
        same(&native, &cell);
        same(&native, &link);
    }
}

#[test]
fn g0_g4_late_pressure_refusal_and_cancellation_leave_the_complete_state() {
    let initial = slab(1.0);
    for fail in [f64::NAN, f64::INFINITY, -1.0, 0.2] {
        let gate = CancelGate::new_clock_free();
        let mut flow = initial.clone();
        assert!(
            flow.step_with_gas_pressure(&gate, 0.05, |x| {
                Ok(if x[0] == 10 { fail } else { 0.001 })
            })
            .is_err()
        );
        same(&flow, &initial);
    }
    for bound in [0.0, -1.0, 1.0, f64::NAN] {
        let mut flow = initial.clone();
        assert!(
            flow.step_with_gas_pressure(&CancelGate::new_clock_free(), bound, |_| Ok(0.0))
                .is_err()
        );
        same(&flow, &initial);
    }
    let gate = CancelGate::new_clock_free();
    let mut flow = initial.clone();
    assert_eq!(
        flow.step_with_gas_link_pressure(&gate, 0.05, |link| {
            if link.liquid_cell[0] == 10 {
                gate.request();
            }
            Ok(0.001)
        }),
        Err(FreeSurfacePressureError2::Cancelled)
    );
    same(&flow, &initial);
}
