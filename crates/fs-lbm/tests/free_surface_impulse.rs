//! G0/G1/G3: instantaneous momentum acts on retained liquid mass, without remap.

use fs_exec::CancelGate;
use fs_lbm::impulse::{ImpulseControls, ImpulseError};
use fs_lbm::{CS2, Cell, ContactModel, FreeSurface, Grid, Q, equilibrium};

const E: [[f64; 2]; Q] = [
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
    [-1.0, 0.0],
    [0.0, -1.0],
    [1.0, 1.0],
    [-1.0, 1.0],
    [-1.0, -1.0],
    [1.0, -1.0],
];
fn controls() -> ImpulseControls {
    ImpulseControls {
        max_cells: 1024,
        maximum_speed: 0.15,
        maximum_momentum_defect: 1e-11,
        maximum_work_defect: 1e-12,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "actual={a:.17e}, expected={b:.17e}"
    );
}
fn observe(s: &FreeSurface, i: usize) -> (f64, f64, [f64; 2], [[f64; 2]; 2]) {
    let f = s.grid.f[i];
    let rho = f.iter().sum::<f64>();
    let m = match s.grid.flags[i] {
        Cell::Fluid => rho,
        Cell::Interface => s.mass[i],
        _ => 0.0,
    };
    let force: [f64; 2] = core::array::from_fn(|a| s.grid.g[a].mul_add(rho, s.grid.fext[i][a]));
    let u: [f64; 2] = core::array::from_fn(|a| {
        (f.iter().zip(E).map(|(f, e)| f * e[a]).sum::<f64>() + 0.5 * force[a]) / rho
    });
    let stress = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            f.iter().zip(E).map(|(f, e)| f * e[a] * e[b]).sum::<f64>()
                - rho * u[a] * u[b]
                - if a == b { rho * CS2 } else { 0.0 }
                + 0.5 * (u[a] * force[b] + u[b] * force[a])
        })
    });
    (rho, m, u, stress)
}
fn fixture() -> FreeSurface {
    let mut grid = Grid::uniform(6, 3, 0.9);
    for i in 0..18 {
        grid.f[i] = equilibrium(0.8 + 0.03 * i as f64, 0.03, -0.02);
        grid.f[i][1] += 0.0005;
        grid.f[i][3] += 0.0005;
        grid.f[i][0] -= 0.001;
    }
    grid.flags[0] = Cell::Wall;
    grid.flags[17] = Cell::Gas;
    let mut s = FreeSurface::new(grid, 0.01, ContactModel::Neutral);
    s.rebind_forcing([0.0002, -0.0001], &vec![[0.0001, 0.00005]; 18])
        .unwrap();
    for i in 0..18 {
        if s.grid.flags[i] == Cell::Interface {
            s.mass[i] = [0.125, 0.5, 1.0][i % 3] * s.grid.f[i].iter().sum::<f64>();
        }
    }
    // This lane is metadata only for Fluid; its inertia is the population sum.
    assert_eq!(s.grid.flags[8], Cell::Fluid);
    s.mass[8] = 123.0;
    s.carry = 0.0125;
    s
}
fn unchanged(s: &FreeSurface, old: &FreeSurface) {
    assert_eq!(s.grid.f, old.grid.f);
    assert_eq!(s.grid.flags, old.grid.flags);
    assert_eq!(s.mass, old.mass);
    assert_eq!(s.grid.g, old.grid.g);
    assert_eq!(s.grid.fext, old.grid.fext);
    assert_eq!(s.carry, old.carry);
    assert_eq!(s.ledger_mass().to_bits(), old.ledger_mass().to_bits());
}

#[test]
fn liquid_inertia_stress_signed_work_and_exact_fill_survive_a_kick() {
    let mut s = fixture();
    let old = s.clone();
    let inputs: Vec<_> = s
        .grid
        .flags
        .iter()
        .map(|cell| {
            if matches!(cell, Cell::Fluid | Cell::Interface) {
                [-0.001, 0.0004]
            } else {
                [0.0; 2]
            }
        })
        .collect();
    let report = s
        .apply_impulse_density(&inputs, controls(), &CancelGate::new_clock_free())
        .unwrap();
    let mut work = 0.0;
    let mut kinetic = 0.0;
    let mut total_j = [0.0; 2];
    for i in 0..18 {
        if inputs[i] == [0.0; 2] {
            assert_eq!(s.grid.f[i], old.grid.f[i]);
            continue;
        }
        let (rho, m, before, stress) = observe(&old, i);
        let (rho_after, m_after, after, stress_after) = observe(&s, i);
        assert_eq!(rho_after.to_bits(), rho.to_bits());
        assert_eq!(m_after.to_bits(), m.to_bits());
        assert_eq!(s.fill(i).to_bits(), old.fill(i).to_bits());
        for a in 0..2 {
            close(after[a] - before[a], inputs[i][a] / m, 2e-16);
            for b in 0..2 {
                close(stress_after[a][b], stress[a][b], 2e-15);
            }
            total_j[a] += inputs[i][a];
            work += inputs[i][a] * (before[a] + after[a]) / 2.0;
            kinetic += 0.5 * m * (after[a] * after[a] - before[a] * before[a]);
        }
    }
    assert!(report.work < 0.0);
    assert_eq!(report.mass_change, 0.0);
    close(report.work, work, 1e-17);
    close(report.kinetic_energy_change, kinetic, 1e-17);
    for a in 0..2 {
        close(report.momentum_change[a], total_j[a], 2e-15);
    }
    assert_eq!(s.mass, old.mass);
    assert_eq!(s.grid.flags, old.grid.flags);
    assert_eq!(s.grid.g, old.grid.g);
    assert_eq!(s.grid.fext, old.grid.fext);
    assert_eq!(s.grid.tau, old.grid.tau);
    assert_eq!(s.sigma, old.sigma);
    assert_eq!(s.contact, old.contact);
    assert_eq!(s.carry, old.carry);
    assert_eq!(s.ledger_mass().to_bits(), old.ledger_mass().to_bits());
    assert_eq!(s.conversions.to_fluid, old.conversions.to_fluid);
    assert_eq!(s.conversions.to_gas, old.conversions.to_gas);
}

#[test]
fn half_the_liquid_mass_gives_twice_the_velocity_jump() {
    let mut grid = Grid::uniform(2, 1, 0.9);
    grid.flags.fill(Cell::Interface);
    let mut s = FreeSurface::new(grid, 0.0, ContactModel::Neutral);
    let rho = s.grid.f[0].iter().sum::<f64>();
    s.mass = vec![0.25 * rho, 0.5 * rho];
    let report = s
        .apply_impulse_density(
            &[[0.005, 0.0]; 2],
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let u0 = observe(&s, 0).2[0];
    let u1 = observe(&s, 1).2[0];
    close(u0, 0.005 / s.mass[0], 2e-16);
    close(u1, 0.005 / s.mass[1], 2e-16);
    close(u0, 2.0 * u1, 2e-16);
    close(report.work, 0.5 * 0.005 * (u0 + u1), 2e-18);
}

#[test]
fn zero_dry_and_late_invalid_impulses_preserve_all_accepted_liquid_state() {
    let mut s = fixture();
    let old = s.clone();
    let gate = CancelGate::new_clock_free();
    let zero = s
        .apply_impulse_density(&[[0.0; 2]; 18], controls(), &gate)
        .unwrap();
    assert_eq!(zero.affected_cells, 0);
    unchanged(&s, &old);
    for invalid_cell in [0, 17] {
        let mut impulses = [[0.0; 2]; 18];
        impulses[8] = [0.001, 0.0];
        impulses[invalid_cell] = [0.001, 0.0];
        assert!(
            s.apply_impulse_density(&impulses, controls(), &gate)
                .is_err()
        );
        unchanged(&s, &old);
    }
    let mut impulses = [[0.0; 2]; 18];
    impulses[8] = [0.001, 0.0];
    for bad in [f64::NAN, f64::INFINITY, 1.0] {
        impulses[16] = [bad, 0.0];
        assert!(
            s.apply_impulse_density(&impulses, controls(), &gate)
                .is_err()
        );
        unchanged(&s, &old);
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        s.apply_impulse_density(&[[0.0; 2]; 18], controls(), &cancelled),
        Err(ImpulseError::Cancelled)
    );
    unchanged(&s, &old);
}

#[test]
fn invalid_interface_mass_is_neither_clamped_nor_published() {
    let gate = CancelGate::new_clock_free();
    for mass in [0.0, -0.1, f64::INFINITY, 2.0] {
        let mut s = fixture();
        assert_eq!(s.grid.flags[16], Cell::Interface);
        s.mass[16] = mass;
        let old = s.clone();
        let mut impulses = [[0.0; 2]; 18];
        impulses[8] = [0.001, 0.0];
        impulses[16] = [0.001, 0.0];
        assert!(
            s.apply_impulse_density(&impulses, controls(), &gate)
                .is_err()
        );
        unchanged(&s, &old);
    }
}
