//! G1 traction/conversion separation and G3/G4 replay and atomic publication.
use super::*;

const ORIGIN: [f64; 2] = [1.25, -2.75];

fn front(source_mass: f64) -> FreeSurface {
    let mut grid = Grid::uniform(9, 7, 0.8);
    grid.periodic_x = false;
    grid.periodic_y = false;
    grid.flags.fill(Cell::Gas);
    grid.g = [0.00007, -0.00011];
    // Negative debt needs two surviving recipients to represent the transferred
    // spin. The separate rank-refusal test removes the third occupied cell.
    for (x, y, rho, u) in [
        (3, 3, 1.01, [0.025, -0.012]),
        (4, 3, 0.99, [-0.018, 0.021]),
        (4, 4, 1.02, [0.007, -0.006]),
    ] {
        let i = grid.idx(x, y);
        grid.flags[i] = Cell::Interface;
        grid.fext[i] = [0.0002, -0.0001];
        grid.f[i] = equilibrium(rho, u[0], u[1]);
        let force = core::array::from_fn(|d| rho * grid.g[d] + grid.fext[i][d]);
        crate::thermal::forcing::reframe2(&mut grid.f[i], u, [0.0; 2], force);
    }
    let mut flow = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    let source = flow.grid.idx(3, 3);
    let receiver = flow.grid.idx(4, 3);
    flow.mass[source] = source_mass;
    flow.mass[receiver] = 0.55;
    let second_receiver = flow.grid.idx(4, 4);
    flow.mass[second_receiver] = 0.45;
    flow
}

fn physical_momentum(flow: &FreeSurface, i: usize) -> [f64; 2] {
    let rho = flow.grid.f[i].iter().sum::<f64>();
    core::array::from_fn(|d| {
        flow.grid.f[i]
            .iter()
            .zip(E)
            .map(|(f, e)| f * f64::from([e.0, e.1][d]))
            .sum::<f64>()
            + 0.5 * flow.grid.g[d].mul_add(rho, flow.grid.fext[i][d])
    })
}

fn inventory(flow: &FreeSurface) -> (f64, [f64; 2], f64, f64) {
    let mut mass = 0.0;
    let mut momentum = [0.0; 2];
    let mut angular = 0.0;
    let mut kinetic = 0.0;
    for i in 0..flow.mass.len() {
        if !matches!(flow.grid.flags[i], Cell::Fluid | Cell::Interface) {
            continue;
        }
        let rho = flow.grid.f[i].iter().sum::<f64>();
        let m = if flow.grid.flags[i] == Cell::Fluid {
            rho
        } else {
            flow.mass[i]
        };
        let u = physical_momentum(flow, i).map(|p| p / rho);
        mass += m;
        for d in 0..2 {
            momentum[d] += m * u[d];
            kinetic += 0.5 * m * u[d] * u[d];
        }
        let position = [
            (i % flow.grid.nx) as f64 + 0.5 - ORIGIN[0],
            (i / flow.grid.nx) as f64 + 0.5 - ORIGIN[1],
        ];
        angular += m * (position[0] * u[1] - position[1] * u[0]);
    }
    (mass, momentum, angular, kinetic)
}

fn same(a: &FreeSurface, b: &FreeSurface) {
    assert_eq!(a.conversion_redistribution(), b.conversion_redistribution());
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
    assert_eq!(a.grid.fext, b.grid.fext);
    assert_eq!(a.grid.tau, b.grid.tau);
    assert_eq!(a.grid.g, b.grid.g);
    assert_eq!(a.scratch, b.scratch);
    assert_eq!(a.fill_raw, b.fill_raw);
    assert_eq!(a.fill_smooth, b.fill_smooth);
}

#[test]
fn g1_g3_loaded_front_separates_traction_from_signed_conversion_momentum() {
    let gate = CancelGate::new_clock_free();
    for policy in [
        ConversionRedistribution2::ConservativeLocal,
        ConversionRedistribution2::ConservativeLocalAngular,
    ] {
        for source_mass in [1.25, -0.06] {
            let mut initial = front(source_mass);
            initial.set_conversion_redistribution(policy).unwrap();
            let pressure =
                |link: GasPressureLink2| Ok(CS2 * 0.0005 * link.incoming_direction as f64);
            let boundary = GasPressure2::new(Field::Link(&pressure), 0.05).unwrap();
            // Observe the original pressure/force step BEFORE any cell conversion.
            let mut transported = initial.clone();
            transported
                .stream_step_with_pressure(Some(&gate), Some(&boundary))
                .unwrap();
            let mut unloaded = initial.clone();
            unloaded
                .stream_step_with_pressure(Some(&gate), None)
                .unwrap();
            let mut traction = [0.0; 2];
            let mut expected_traction = [0.0; 2];
            let mut torque = 0.0;
            let mut expected_torque = 0.0;
            for (x, y) in [(3, 3), (4, 3), (4, 4)] {
                let i = initial.grid.idx(x, y);
                let arm = [x as f64 + 0.5 - ORIGIN[0], y as f64 + 0.5 - ORIGIN[1]];
                let rho = initial.grid.f[i].iter().sum::<f64>();
                let u = physical_momentum(&initial, i).map(|p| p / rho);
                let u2 = u[0] * u[0] + u[1] * u[1];
                let density_difference = transported.grid.f[i].iter().sum::<f64>()
                    - unloaded.grid.f[i].iter().sum::<f64>();
                let force_difference = initial.grid.g.map(|g| g * density_difference);
                let udf = u[0] * force_difference[0] + u[1] * force_difference[1];
                for q in 1..Q {
                    let Some(source) = initial.grid.source(x, y, q) else {
                        continue;
                    };
                    if initial.grid.flags[source] != Cell::Gas {
                        continue;
                    }
                    let eu = f64::from(E[q].0) * u[0] + f64::from(E[q].1) * u[1];
                    let edf = f64::from(E[q].0) * force_difference[0]
                        + f64::from(E[q].1) * force_difference[1];
                    // Restore the physical paired load at the two actual
                    // endpoint densities. Held fext cancels; rho*g does not.
                    let delta = 2.0 * crate::W[q] * (0.0005 * q as f64)
                        * (1.0 + 4.5 * eu * eu - 1.5 * u2)
                        - crate::W[q] * (-udf / CS2 + eu * edf / (CS2 * CS2));
                    assert!(
                        (transported.grid.f[i][q] - unloaded.grid.f[i][q] - delta).abs() < 3e-16
                    );
                    let imposed: [f64; 2] = core::array::from_fn(|d| {
                        delta * (f64::from([E[q].0, E[q].1][d]) + 0.5 * initial.grid.g[d])
                    });
                    for d in 0..2 {
                        expected_traction[d] += imposed[d];
                    }
                    expected_torque += arm[0] * imposed[1] - arm[1] * imposed[0];
                }
                let loaded_p = physical_momentum(&transported, i);
                let unloaded_p = physical_momentum(&unloaded, i);
                for d in 0..2 {
                    traction[d] += loaded_p[d] - unloaded_p[d];
                }
                torque +=
                    arm[0] * (loaded_p[1] - unloaded_p[1]) - arm[1] * (loaded_p[0] - unloaded_p[0]);
            }
            assert!(expected_traction[0].hypot(expected_traction[1]) > 1e-5);
            for d in 0..2 {
                assert!((traction[d] - expected_traction[d]).abs() < 2e-15);
            }
            assert!(
                expected_torque.abs() > 1e-7,
                "the pressure loading must exert a real moment"
            );
            assert!((torque - expected_torque).abs() < 2e-14);

            // Pressure traction above is measured from actual population moments.
            // This independent inventory uses signed tracked liquid mass. Its
            // conversion-only moment starts AFTER pressure/force streaming, so
            // physical pressure torque is never included in a remap defect.
            let before_conversion = inventory(&transported);
            let mut accepted = initial.clone();
            let report = accepted
                .step_with_gas_link_pressure_and_momentum(&gate, 0.05, pressure)
                .unwrap();
            let after = inventory(&accepted);
            assert!(report.routes > 0 && report.redistributed_mass > 0.0);
            assert!((after.0 - before_conversion.0).abs() < 3e-13);
            for d in 0..2 {
                assert!((after.1[d] - before_conversion.1[d]).abs() < 3e-14);
                assert!(report.absolute_momentum_defect[d] < 3e-14);
            }
            if policy == ConversionRedistribution2::ConservativeLocalAngular {
                assert!((after.2 - before_conversion.2).abs() < 3e-13);
                assert!(report.absolute_angular_momentum_defect.unwrap() < 3e-14);
            } else {
                assert_eq!(report.absolute_angular_momentum_defect, None);
            }
            assert!((after.3 - before_conversion.3 - report.kinetic_energy_change).abs() < 3e-15);
            assert_eq!(accepted.carry, 0.0);
            let source = accepted.grid.idx(3, 3);
            assert_eq!(
                accepted.grid.flags[source],
                if source_mass > 0.0 {
                    Cell::Fluid
                } else {
                    Cell::Gas
                }
            );

            let mut zero = initial.clone();
            let mut unforced = initial.clone();
            let zero_report = zero
                .step_with_gas_pressure_and_momentum(&gate, 0.05, |_| Ok(0.0))
                .unwrap();
            assert_eq!(zero_report, unforced.step_with_momentum(&gate).unwrap());
            same(&zero, &unforced);
            let mut cell = initial.clone();
            let mut link = initial;
            let a = cell
                .step_with_gas_pressure_and_momentum(&gate, 0.05, |_| Ok(0.001))
                .unwrap();
            let b = link
                .step_with_gas_link_pressure_and_momentum(&gate, 0.05, |_| Ok(0.001))
                .unwrap();
            assert_eq!(a, b);
            same(&cell, &link);
        }
    }
}

#[test]
fn g4_pressure_conversion_and_cancellation_refusals_preserve_the_complete_candidate() {
    let mut initial = front(1.25);
    initial
        .set_conversion_redistribution(ConversionRedistribution2::ConservativeLocalAngular)
        .unwrap();
    for cancel in [false, true] {
        let gate = CancelGate::new_clock_free();
        let mut flow = initial.clone();
        let result = flow.step_with_gas_link_pressure_and_momentum(&gate, 0.05, |link| {
            if link.liquid_cell[0] == 4 {
                if cancel {
                    gate.request();
                } else {
                    return Err(FreeSurfacePressureError2::Invalid {
                        coordinates: Some(link.liquid_cell),
                        reason: "prescribed gas source refused",
                    });
                }
            }
            Ok(0.001)
        });
        if cancel {
            assert_eq!(
                result,
                Err(FreeSurfacePressureMomentumError2::Pressure(
                    FreeSurfacePressureError2::Cancelled
                ))
            );
        } else {
            assert!(matches!(
                result,
                Err(FreeSurfacePressureMomentumError2::Pressure(
                    FreeSurfacePressureError2::Invalid {
                        coordinates: Some([4, 3]),
                        ..
                    }
                ))
            ));
        }
        same(&flow, &initial);
    }
    let mut grid = Grid::uniform(3, 3, 0.8);
    grid.periodic_x = false;
    grid.periodic_y = false;
    grid.flags.fill(Cell::Wall);
    grid.flags[4] = Cell::Interface;
    grid.flags[5] = Cell::Gas;
    let mut flow = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    flow.set_conversion_redistribution(ConversionRedistribution2::ConservativeLocalAngular)
        .unwrap();
    flow.mass[4] = 3.0; // Two units of excess cannot fit the one adjacent gas cell.
    let before = flow.clone();
    let visits = std::cell::Cell::new(0);
    let result =
        flow.step_with_gas_pressure_and_momentum(&CancelGate::new_clock_free(), 0.05, |_| {
            visits.set(visits.get() + 1);
            Ok(0.001)
        });
    assert!(
        visits.get() > 0,
        "the pressure stage must have actually completed"
    );
    assert!(matches!(
        result,
        Err(FreeSurfacePressureMomentumError2::Momentum(
            FreeSurfaceMomentumError2::Invalid(_)
        ))
    ));
    same(&flow, &before);
}

#[test]
fn g4_pressure_loaded_single_survivor_refuses_unrepresentable_conversion_spin() {
    let gate = CancelGate::new_clock_free();
    let mut initial = front(-0.06);
    let removed = initial.grid.idx(4, 4);
    initial.grid.flags[removed] = Cell::Gas;
    initial.mass[removed] = 0.0;
    // The linear policy can complete this pressure step, so this is a spin
    // representation refusal rather than absent capacity or invalid traction.
    let mut linear = initial.clone();
    let report = linear
        .step_with_gas_pressure_and_momentum(&gate, 0.05, |_| Ok(0.001))
        .unwrap();
    assert!(report.routes > 0);
    initial
        .set_conversion_redistribution(ConversionRedistribution2::ConservativeLocalAngular)
        .unwrap();
    let before = initial.clone();
    let visits = std::cell::Cell::new(0);
    let result = initial.step_with_gas_pressure_and_momentum(&gate, 0.05, |_| {
        visits.set(visits.get() + 1);
        Ok(0.001)
    });
    assert!(visits.get() > 0);
    assert!(matches!(
        result,
        Err(FreeSurfacePressureMomentumError2::Momentum(
            FreeSurfaceMomentumError2::Invalid(_)
        ))
    ));
    same(&initial, &before);
}
