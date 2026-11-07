//! G0/G1/G3/G4: original sparse mass, Guo-frame moments and atomic publication.

use super::*;
use crate::d3q19::equilibrium3;
use crate::d3q19::freesurface3::ContactModel3;
use crate::impulse::ImpulseControls;
use std::collections::BTreeMap;

fn controls() -> WrenchImpulseControls3 {
    WrenchImpulseControls3 {
        impulse: ImpulseControls {
            max_cells: 128,
            maximum_speed: 0.15,
            maximum_momentum_defect: 2e-12,
            maximum_work_defect: 2e-13,
        },
        maximum_angular_momentum_defect: 5e-12,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn coordinates(index: usize) -> [usize; 3] {
    [index % 4, index / 4 % 4, index / 16]
}

fn position(index: usize) -> [f64; 3] {
    coordinates(index).map(|v| v as f64 + 0.5)
}

// Seed a retained VOF inventory and a nonzero physical shear moment in the
// original sparse owner. This fixture performs no collision or transport.
fn fixture(boost: [f64; 3]) -> FreeSurface3 {
    let mut flow = FreeSurface3::new(
        4,
        4,
        4,
        0.9,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        |_, _, _| true,
    )
    .unwrap();
    let key = flow.grid.active_keys()[0];
    let mut forces = BTreeMap::new();
    for lane in 0..64 {
        let rho = 0.96 + 0.002 * lane as f64;
        let force = [0.0001 + lane as f64 * 1e-6, -0.0002, 0.0003];
        let velocity = core::array::from_fn(|a| [0.02, -0.01, 0.015][a] + boost[a]);
        let mut f = equilibrium3(rho, velocity);
        for q in 0..19 {
            let e = [f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2)];
            let eu: f64 = e.into_iter().zip(velocity).map(|(e, u)| e * u).sum();
            let ef: f64 = e.into_iter().zip(force).map(|(e, f)| e * f).sum();
            let uf: f64 = velocity.into_iter().zip(force).map(|(u, f)| u * f).sum();
            f[q] += W3[q] * (9.0 * 0.0004 * e[0] * e[1] - 0.5 * (3.0 * (ef - uf) + 9.0 * eu * ef));
        }
        flow.grid.set_populations(0, lane, f);
        let rho = f.iter().sum::<f64>();
        let fill = match lane {
            0 => 0.2,
            3 => 0.5,
            12 => 0.75,
            _ => 1.0,
        };
        let tile = flow.tiles.get_mut(&key).unwrap();
        tile.cells[lane] = if fill < 1.0 {
            Cell3::Interface
        } else {
            Cell3::Fluid
        };
        tile.mass[lane] = fill * rho;
        forces.insert((key, lane), force);
    }
    flow.local_force = Some(forces);
    flow
}

fn mass(flow: &FreeSurface3, index: usize) -> f64 {
    let [x, y, z] = coordinates(index);
    flow.region_mass(|a, b, c| [a, b, c] == [x, y, z])
}

fn velocity_and_stress(flow: &FreeSurface3, index: usize) -> ([f64; 3], [[f64; 3]; 3]) {
    let [x, y, z] = coordinates(index);
    let (slot, lane) = flow
        .grid
        .resolve_source(x as i64, y as i64, z as i64)
        .unwrap();
    let f = flow.grid.populations(slot, lane);
    let rho = f.iter().sum::<f64>();
    let key = flow.grid.active_keys()[slot];
    let force = flow.force_at_lane(key, lane);
    let e = E3.map(|(x, y, z)| [f64::from(x), f64::from(y), f64::from(z)]);
    let velocity = core::array::from_fn(|a| {
        (f.iter().zip(e).map(|(f, e)| f * e[a]).sum::<f64>() + 0.5 * force[a]) / rho
    });
    let stress = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            f.iter().zip(e).map(|(f, e)| f * e[a] * e[b]).sum::<f64>()
                + 0.5 * (velocity[a] * force[b] + velocity[b] * force[a])
                - rho * velocity[a] * velocity[b]
                - if a == b { rho / 3.0 } else { 0.0 }
        })
    });
    (velocity, stress)
}

fn invariants(before: &FreeSurface3, after: &FreeSurface3) {
    assert_eq!(
        before.ledger_mass().to_bits(),
        after.ledger_mass().to_bits()
    );
    assert_eq!(before.carry().to_bits(), after.carry().to_bits());
    assert_eq!(before.conversions(), after.conversions());
    assert_eq!(before.geometry_epoch(), after.geometry_epoch());
    assert_eq!(before.grid.steps(), after.grid.steps());
    assert_eq!(before.local_force, after.local_force);
    for i in 0..64 {
        let [x, y, z] = coordinates(i);
        let a = before.sample_at(x, y, z).unwrap();
        let b = after.sample_at(x, y, z).unwrap();
        assert_eq!(a.cell, b.cell);
        assert_eq!(a.density.to_bits(), b.density.to_bits());
        assert_eq!(a.fill_fraction.to_bits(), b.fill_fraction.to_bits());
    }
}

// Construct the requested impulse from an independently authored rigid
// velocity jump. Unequal masses and partial fill make density-as-inertia fail.
fn authored(
    flow: &FreeSurface3,
    indices: &[usize],
    origin: [f64; 3],
) -> (WrenchImpulse3, Vec<[f64; 3]>, f64, [f64; 3]) {
    let total: f64 = indices.iter().map(|&i| mass(flow, i)).sum();
    let center = core::array::from_fn(|a| {
        indices
            .iter()
            .map(|&i| mass(flow, i) * position(i)[a])
            .sum::<f64>()
            / total
    });
    let jumps: Vec<_> = indices
        .iter()
        .map(|&i| {
            let arm = core::array::from_fn(|a| position(i)[a] - center[a]);
            let rotation = cross([0.003, -0.002, 0.004], arm);
            core::array::from_fn(|a| [0.005, -0.003, 0.002][a] + rotation[a])
        })
        .collect();
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for (&index, &jump) in indices.iter().zip(&jumps) {
        let impulse = jump.map(|value| value * mass(flow, index));
        let arm = core::array::from_fn(|a| position(index)[a] - origin[a]);
        let moment = cross(arm, impulse);
        for a in 0..3 {
            linear[a] += impulse[a];
            angular[a] += moment[a];
        }
    }
    (
        WrenchImpulse3 {
            origin,
            linear,
            angular,
        },
        jumps,
        total,
        center,
    )
}

#[test]
fn partial_liquid_wrench_matches_original_moments_and_retains_physical_stress() {
    let mut flow = fixture([0.0; 3]);
    let before = flow.clone();
    let indices = [0, 3, 12, 48];
    let (wrench, jumps, total, center) = authored(&flow, &indices, [-0.75, 0.25, 1.5]);
    let report = flow
        .apply_wrench_impulse(&indices, wrench, controls(), &CancelGate::new_clock_free())
        .unwrap();
    invariants(&before, &flow);
    assert_eq!(report.impulse.affected_cells, indices.len());
    assert_eq!(report.impulse.mass_change, 0.0);
    close(report.stencil_mass, total, 1e-15);
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    let mut kinetic = 0.0;
    let mut work = 0.0;
    for index in 0..64 {
        let old = velocity_and_stress(&before, index);
        let new = velocity_and_stress(&flow, index);
        if let Some(selected) = indices.iter().position(|&i| i == index) {
            let m = mass(&before, index);
            let impulse = core::array::from_fn(|a| m * (new.0[a] - old.0[a]));
            let arm = core::array::from_fn(|a| position(index)[a] - wrench.origin[a]);
            let moment = cross(arm, impulse);
            for a in 0..3 {
                close(new.0[a] - old.0[a], jumps[selected][a], 4e-15);
                linear[a] += impulse[a];
                angular[a] += moment[a];
                work += m * jumps[selected][a] * (new.0[a] + old.0[a]) * 0.5;
                kinetic += 0.5 * m * (new.0[a].powi(2) - old.0[a].powi(2));
                for b in 0..3 {
                    close(new.1[a][b], old.1[a][b], 5e-16);
                }
            }
            assert!(old.1[0][1].abs() > 1e-4);
        } else {
            assert_eq!(
                before.grid.populations(0, index),
                flow.grid.populations(0, index)
            );
        }
    }
    for a in 0..3 {
        close(report.stencil_center_of_mass[a], center[a], 2e-15);
        close(linear[a], wrench.linear[a], 2e-15);
        close(angular[a], wrench.angular[a], 8e-15);
        close(report.angular_momentum_change[a], angular[a], 2e-17);
    }
    close(report.impulse.work, work, 2e-16);
    close(report.impulse.kinetic_energy_change, kinetic, 2e-16);
}

#[test]
fn origin_shift_and_galilean_boost_keep_the_same_liquid_velocity_jump() {
    let indices = [0, 3, 12, 48];
    let mut rest = fixture([0.0; 3]);
    let boost = [-0.015, 0.007, 0.011];
    let mut moving = fixture(boost);
    let (wrench, _, _, _) = authored(&rest, &indices, [0.0; 3]);
    let shifted_origin = [2.5, -1.0, 4.0];
    let moment = cross(shifted_origin, wrench.linear);
    let shifted = WrenchImpulse3 {
        origin: shifted_origin,
        angular: core::array::from_fn(|a| wrench.angular[a] - moment[a]),
        ..wrench
    };
    let gate = CancelGate::new_clock_free();
    let a = rest
        .apply_wrench_impulse(&indices, wrench, controls(), &gate)
        .unwrap();
    let b = moving
        .apply_wrench_impulse(&indices, shifted, controls(), &gate)
        .unwrap();
    for index in indices {
        let u = velocity_and_stress(&rest, index).0;
        let v = velocity_and_stress(&moving, index).0;
        for axis in 0..3 {
            close(v[axis] - u[axis], boost[axis], 4e-15);
        }
    }
    let boost_work: f64 = boost
        .into_iter()
        .zip(wrench.linear)
        .map(|(u, j)| u * j)
        .sum();
    close(b.impulse.work - a.impulse.work, boost_work, 2e-16);
    close(
        b.impulse.kinetic_energy_change - a.impulse.kinetic_energy_change,
        boost_work,
        2e-16,
    );
}

#[test]
fn sparse_line_supports_transverse_couple_but_refuses_unrepresented_spin() {
    let mut flow = fixture([0.0; 3]);
    let gate = CancelGate::new_clock_free();
    let wrench = WrenchImpulse3 {
        origin: [0.5; 3],
        linear: [0.0; 3],
        angular: [0.0, 0.004, -0.006],
    };
    let report = flow
        .apply_wrench_impulse(&[0, 1, 2], wrench, controls(), &gate)
        .unwrap();
    for a in 0..3 {
        close(report.angular_momentum_change[a], wrench.angular[a], 4e-15);
    }
    let bits = flow.grid.state_bits();
    assert!(
        flow.apply_wrench_impulse(
            &[0, 1, 2],
            WrenchImpulse3 {
                angular: [0.001, 0.0, 0.0],
                ..wrench
            },
            controls(),
            &gate
        )
        .is_err()
    );
    assert_eq!(flow.grid.state_bits(), bits);
}

#[test]
fn empty_dry_duplicate_unlocated_and_late_speed_refusals_publish_nothing() {
    let mut flow = fixture([0.0; 3]);
    let indices = [0, 3, 12, 48];
    let (wrench, _, _, _) = authored(&flow, &indices, [0.0; 3]);
    let gate = CancelGate::new_clock_free();
    let original = flow.grid.state_bits();
    for indices in [&[][..], &[0, 0], &[3, 0], &[0, 64]] {
        assert!(
            flow.apply_wrench_impulse(indices, wrench, controls(), &gate)
                .is_err()
        );
        assert_eq!(flow.grid.state_bits(), original);
    }
    let key = flow.grid.active_keys()[0];
    flow.tiles.get_mut(&key).unwrap().cells[48] = Cell3::Gas;
    assert!(
        flow.apply_wrench_impulse(&indices, wrench, controls(), &gate)
            .is_err()
    );
    assert_eq!(flow.grid.state_bits(), original);
    flow.tiles.get_mut(&key).unwrap().cells[48] = Cell3::Interface;
    flow.tiles.get_mut(&key).unwrap().mass[48] = 2.0;
    assert!(
        flow.apply_wrench_impulse(&indices, wrench, controls(), &gate)
            .is_err()
    );
    assert_eq!(flow.grid.state_bits(), original);
    flow.tiles.get_mut(&key).unwrap().mass[48] = 0.3;
    let enormous = WrenchImpulse3 {
        linear: [0.0, 1.0, 0.0],
        ..wrench
    };
    assert!(
        flow.apply_wrench_impulse(&indices, enormous, controls(), &gate)
            .is_err()
    );
    assert_eq!(flow.grid.state_bits(), original);
}

#[test]
fn dependent_receipt_and_cancellation_are_atomic_and_zero_wrench_preserves_bits() {
    let mut flow = fixture([0.0; 3]);
    let before = flow.clone();
    let indices = [0, 3, 12, 48];
    let (wrench, _, _, _) = authored(&flow, &indices, [0.0; 3]);
    let gate = CancelGate::new_clock_free();
    let zero = WrenchImpulse3 {
        linear: [0.0; 3],
        angular: [0.0; 3],
        ..wrench
    };
    let report = flow
        .apply_wrench_impulse(&indices, zero, controls(), &gate)
        .unwrap();
    assert_eq!(report.impulse.affected_cells, 0);
    assert_eq!(before.grid.state_bits(), flow.grid.state_bits());
    let refused: Result<(), _> =
        flow.apply_wrench_impulse_with_receipt(&indices, wrench, controls(), &gate, |_| {
            Err(ImpulseError::Balance)
        });
    assert_eq!(refused, Err(ImpulseError::Balance));
    assert_eq!(before.grid.state_bits(), flow.grid.state_bits());
    let cancelled =
        flow.apply_wrench_impulse_with_receipt(&indices, wrench, controls(), &gate, |report| {
            gate.request();
            Ok(report)
        });
    assert_eq!(cancelled, Err(ImpulseError::Cancelled));
    assert_eq!(before.grid.state_bits(), flow.grid.state_bits());
    invariants(&before, &flow);
}
