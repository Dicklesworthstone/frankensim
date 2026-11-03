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

fn advancing_front() -> FreeSurface {
    let mut flow = slab(1.0);
    flow.grid.g = [0.0001, -0.0002];
    let donor = flow.grid.idx(5, 6);
    let velocity = [0.025, -0.008];
    let mut populations = equilibrium(1.0, velocity[0], velocity[1]);
    crate::thermal::forcing::reframe2(&mut populations, velocity, [0.0; 2], flow.grid.g);
    flow.grid.f[donor] = populations;
    flow.mass[donor] = 1.2 * populations.iter().sum::<f64>();
    flow
}

// Independently sum the actual signed liquid inventories. The pre-conversion
// snapshot includes excess/debt; gas populations carry no represented mass.
fn inventory(flow: &FreeSurface) -> (f64, [f64; 2], f64) {
    let mut mass = 0.0;
    let mut momentum = [0.0; 2];
    let mut kinetic = 0.0;
    for i in 0..flow.grid.f.len() {
        if !matches!(flow.grid.flags[i], Cell::Fluid | Cell::Interface) {
            continue;
        }
        let state = flow.grid.moments(i);
        let m = if flow.grid.flags[i] == Cell::Fluid {
            state.rho
        } else {
            flow.mass[i]
        };
        mass += m;
        for d in 0..2 {
            momentum[d] += m * state.u[d];
        }
        kinetic += 0.5 * m * (state.u[0] * state.u[0] + state.u[1] * state.u[1]);
    }
    (mass, momentum, kinetic)
}

#[test]
fn g1_combined_loaded_front_preserves_conversion_momentum_with_real_routes() {
    let gate = CancelGate::new_clock_free();
    let mut flow = advancing_front();
    let pressure = |_: [usize; 2]| Ok(0.002 * CS2);
    // Measure immediately after the ORIGINAL pressure/force streaming owner,
    // before conversion. Boundary and body-force impulses belong to this part
    // of the step and must not be mistaken for a conversion residual.
    let boundary = GasPressure2::new(Field::Cell(&pressure), 0.05).unwrap();
    let mut streamed = flow.clone();
    streamed
        .stream_step_with_pressure(Some(&gate), Some(&boundary))
        .unwrap();
    let before = inventory(&streamed);
    let report = flow
        .step_with_gas_pressure_and_momentum(&gate, 0.05, pressure)
        .unwrap();
    let after = inventory(&flow);
    assert!(report.routes > 0);
    assert!(report.redistributed_mass > 0.1);
    assert_eq!(flow.carry, 0.0);
    assert!((after.0 - before.0).abs() < 2e-12);
    for d in 0..2 {
        assert!((after.1[d] - before.1[d]).abs() < 5e-13);
        assert!(report.absolute_momentum_defect[d] < 5e-13);
    }
    assert!((after.2 - before.2 - report.kinetic_energy_change).abs() < 1e-14);
    for i in 0..flow.grid.f.len() {
        if flow.grid.flags[i] == Cell::Interface {
            assert!(flow.mass[i] >= 0.0);
            assert!(flow.mass[i] <= flow.grid.f[i].iter().sum::<f64>());
        }
    }
}

#[test]
fn g3_combined_ambient_cell_and_link_pressure_replay_momentum_conversion_bits() {
    let gate = CancelGate::new_clock_free();
    let mut ordinary = advancing_front();
    let mut cell = ordinary.clone();
    let mut link = ordinary.clone();
    let expected = ordinary.step_with_momentum(&gate).unwrap();
    let cell_report = cell
        .step_with_gas_pressure_and_momentum(&gate, 0.1, |_| Ok(0.0))
        .unwrap();
    let link_report = link
        .step_with_gas_link_pressure_and_momentum(&gate, 0.1, |_| Ok(0.0))
        .unwrap();
    assert!(expected.routes > 0);
    assert_eq!(cell_report, expected);
    assert_eq!(link_report, expected);
    same(&ordinary, &cell);
    same(&ordinary, &link);
}

#[test]
fn g1_combined_pressure_drives_actual_forced_populations_in_the_physical_frame() {
    let gate = CancelGate::new_clock_free();
    let mut flow = slab(1.02);
    flow.grid.g = [0.0002, -0.0003];
    let velocity = [0.025, -0.01];
    for i in 0..flow.grid.f.len() {
        flow.grid.fext[i] = [0.0001, -0.0001];
        if matches!(flow.grid.flags[i], Cell::Fluid | Cell::Interface) {
            let force =
                core::array::from_fn(|d| flow.grid.g[d].mul_add(1.02, flow.grid.fext[i][d]));
            let mut f = equilibrium(1.02, velocity[0], velocity[1]);
            crate::thermal::forcing::reframe2(&mut f, velocity, [0.0; 2], force);
            flow.grid.f[i] = f;
        }
    }
    let cell = flow.grid.idx(5, 6);
    let state = flow.grid.moments(cell);
    let force: [f64; 2] =
        core::array::from_fn(|d| flow.grid.g[d].mul_add(state.rho, flow.grid.fext[cell][d]));
    let mut post = Vec::new();
    flow.grid.collide_into(&mut post);
    let mut ambient = flow.clone();
    ambient.step_with_momentum(&gate).unwrap();
    flow.step_with_gas_pressure_and_momentum(&gate, 0.05, |_| Ok(0.015 * CS2))
        .unwrap();
    let u2 = state.u[0] * state.u[0] + state.u[1] * state.u[1];
    let uf = state.u[0] * force[0] + state.u[1] * force[1];
    let mut normal_impulse = 0.0;
    for q in 1..Q {
        if E[q].1 != -1 {
            continue;
        }
        let e = [f64::from(E[q].0), f64::from(E[q].1)];
        let eu = e[0] * state.u[0] + e[1] * state.u[1];
        let ef = e[0] * force[0] + e[1] * force[1];
        let restored_pair = flow.grid.f[cell][q]
            + post[cell][OPP[q]]
            + crate::W[q] * (-uf / CS2 + eu * ef / (CS2 * CS2));
        let expected = 2.0 * crate::W[q] * 1.015 * (1.0 + 4.5 * eu * eu - 1.5 * u2);
        assert!((restored_pair - expected).abs() < 5e-15);
        // A higher upper gas load must send real downward momentum into the
        // liquid. The physical source-frame term cancels between these loads.
        let difference = flow.grid.f[cell][q] - ambient.grid.f[cell][q];
        let expected_difference = 2.0 * crate::W[q] * 0.015 * (1.0 + 4.5 * eu * eu - 1.5 * u2);
        assert!((difference - expected_difference).abs() < 5e-15);
        normal_impulse += e[1] * difference;
    }
    assert!(normal_impulse < -0.004);
}

#[test]
fn g1_combined_callbacks_visit_only_actual_gas_links_at_the_streaming_time_level() {
    let gate = CancelGate::new_clock_free();
    let mut cell = advancing_front();
    let mut link = cell.clone();
    let cells = std::cell::RefCell::new(Vec::new());
    let links = std::cell::RefCell::new(Vec::new());
    let a = cell
        .step_with_gas_pressure_and_momentum(&gate, 0.05, |p| {
            cells.borrow_mut().push(p);
            Ok(0.002 * CS2)
        })
        .unwrap();
    let b = link
        .step_with_gas_link_pressure_and_momentum(&gate, 0.05, |p| {
            links.borrow_mut().push(p);
            Ok(0.002 * CS2)
        })
        .unwrap();
    assert_eq!(a, b);
    same(&cell, &link);
    assert_eq!(*cells.borrow(), (0..12).map(|x| [x, 6]).collect::<Vec<_>>());
    let mut expected = Vec::new();
    for x in 0..12 {
        for q in 1..Q {
            if E[q].1 == -1 {
                expected.push(GasPressureLink2 {
                    liquid_cell: [x, 6],
                    gas_cell: [(x as i32 - E[q].0).rem_euclid(12) as usize, 7],
                    incoming_direction: q,
                });
            }
        }
    }
    assert_eq!(*links.borrow(), expected);
    assert!(a.routes > 0);
    assert!(cell.conversions.gas_to_interface > 0);
    // The newly exposed y=7 front was created AFTER this field was sampled.
    // It must not receive a second collision or stale previous-step pressure.
    let mut full = FreeSurface::new(Grid::uniform(8, 8, 0.8), 0.0, ContactModel::Neutral);
    full.step_with_gas_link_pressure_and_momentum(&gate, 0.05, |_| {
        panic!("walls and fully flooded interiors must not query gas pressure")
    })
    .unwrap();
}

#[test]
fn g0_g4_combined_pressure_and_conversion_failures_keep_their_owner_and_roll_back() {
    let gate = CancelGate::new_clock_free();
    let initial = slab(1.0);
    for bound in [0.0, -1.0, 1.0, f64::NAN] {
        let mut flow = initial.clone();
        assert!(matches!(
            flow.step_with_gas_pressure_and_momentum(&gate, bound, |_| {
                panic!("invalid bound queried the pressure field")
            }),
            Err(FreeSurfacePressureMomentumError2::Pressure(
                FreeSurfacePressureError2::Invalid {
                    coordinates: None,
                    ..
                }
            ))
        ));
        same(&flow, &initial);
    }
    let authored = FreeSurfacePressureError2::Invalid {
        coordinates: Some([10, 6]),
        reason: "selected pressure source is outside its physical domain",
    };
    let mut flow = initial.clone();
    assert_eq!(
        flow.step_with_gas_pressure_and_momentum(&gate, 0.05, |p| {
            if p == [10, 6] {
                Err(authored.clone())
            } else {
                Ok(0.001)
            }
        }),
        Err(FreeSurfacePressureMomentumError2::Pressure(
            authored.clone()
        ))
    );
    same(&flow, &initial);

    // Streaming and the prescribed load are valid; only the unpayable local
    // conversion debt fails. Preserve that original momentum-owner error.
    let mut grid = Grid::uniform(5, 5, 0.8);
    grid.flags.fill(Cell::Gas);
    let i = grid.idx(2, 2);
    grid.flags[i] = Cell::Fluid;
    let mut debt = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    debt.mass[i] = -0.25;
    let before = debt.clone();
    let visits = std::cell::Cell::new(0);
    assert!(matches!(
        debt.step_with_gas_link_pressure_and_momentum(&gate, 0.05, |_| {
            visits.set(visits.get() + 1);
            Ok(0.0)
        }),
        Err(FreeSurfacePressureMomentumError2::Momentum(
            FreeSurfaceMomentumError2::Invalid(_)
        ))
    ));
    assert_eq!(visits.get(), Q - 1);
    same(&debt, &before);
}

#[test]
fn g4_combined_cancellation_before_or_during_pressure_sampling_is_atomic() {
    let initial = advancing_front();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    let mut flow = initial.clone();
    assert_eq!(
        flow.step_with_gas_pressure_and_momentum(&cancelled, 0.05, |_| {
            panic!("cancelled step queried the pressure field")
        }),
        Err(FreeSurfacePressureMomentumError2::Pressure(
            FreeSurfacePressureError2::Cancelled
        ))
    );
    same(&flow, &initial);
    let gate = CancelGate::new_clock_free();
    assert_eq!(
        flow.step_with_gas_link_pressure_and_momentum(&gate, 0.05, |link| {
            if link.liquid_cell == [10, 6] {
                gate.request();
            }
            Ok(0.002 * CS2)
        }),
        Err(FreeSurfacePressureMomentumError2::Pressure(
            FreeSurfacePressureError2::Cancelled
        ))
    );
    same(&flow, &initial);
    let open = CancelGate::new_clock_free();
    let mut reference = initial.clone();
    let a = flow
        .step_with_gas_pressure_and_momentum(&open, 0.05, |_| Ok(0.002 * CS2))
        .unwrap();
    let b = reference
        .step_with_gas_pressure_and_momentum(&open, 0.05, |_| Ok(0.002 * CS2))
        .unwrap();
    assert_eq!(a, b);
    same(&flow, &reference);
}
