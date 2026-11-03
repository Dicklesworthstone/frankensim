//! SI force histories on the original material-owned dense liquid state.
use super::*;
use fs_lbm::core2::{Cell, ForceRebindError2, Grid};
use fs_lbm::freesurface::{ContactModel, FreeSurface};
use fs_lbm::material::MaterialLbmError;
use fs_lbm::material::free_surface::MaterialFreeSurfaceDomain2;

const E: [[f64; 2]; 9] = [
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

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-14 * expected.abs().max(1.0),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn domain2() -> MaterialFreeSurfaceDomain2 {
    MaterialFreeSurfaceDomain2 {
        cells: [16, 14],
        gravity_m_s2: [0.0; 2],
        contact: ContactModel::Neutral,
        max_cells: 16 * 14,
    }
}

fn physical(grid: &Grid, i: usize) -> [f64; 6] {
    let f = grid.f[i];
    let rho = f.iter().sum::<f64>();
    let force: [f64; 2] = core::array::from_fn(|d| rho * grid.g[d] + grid.fext[i][d]);
    let momentum: [f64; 2] = core::array::from_fn(|d| {
        f.iter().zip(E).map(|(f, e)| f * e[d]).sum::<f64>() + 0.5 * force[d]
    });
    let velocity = momentum.map(|p| p / rho);
    let stress = |a: usize, b: usize| {
        f.iter().zip(E).map(|(f, e)| f * e[a] * e[b]).sum::<f64>()
            + 0.5 * (velocity[a] * force[b] + force[a] * velocity[b])
    };
    [
        rho,
        momentum[0],
        momentum[1],
        stress(0, 0),
        stress(0, 1),
        stress(1, 1),
    ]
}

fn unchanged_geometry(actual: &FreeSurface, before: &FreeSurface) {
    assert_eq!(actual.grid.flags, before.grid.flags);
    assert_eq!(actual.grid.tau, before.grid.tau);
    assert_eq!(actual.mass, before.mass);
    assert_eq!(actual.carry, before.carry);
    assert_eq!(actual.sigma, before.sigma);
    assert_eq!(actual.contact, before.contact);
    let counts = |s: &FreeSurface| {
        [
            s.conversions.to_fluid,
            s.conversions.to_gas,
            s.conversions.gas_to_interface,
            s.conversions.fluid_to_interface,
        ]
    };
    assert_eq!(counts(actual), counts(before));
    for i in 0..actual.grid.f.len() {
        close(actual.fill(i), before.fill(i));
    }
    close(actual.ledger_mass(), before.ledger_mass());
}

#[test]
fn si_load_histories_apply_g_plus_force_over_density_only_when_stepped() {
    for rho in [800.0, 1_600.0] {
        let source = liquid(rho, rho * 0.0001, None);
        for (dx, dt) in [(0.001, 0.001), (0.002, 0.0005)] {
            let scaling = plan_free_surface_from_material(&source, dx, dt, 0.1).unwrap();
            let mut flow = scaling.free_surface2(domain2(), |_, _| true).unwrap();
            let mut expected = [0.0; 2];
            for (gravity, external) in [
                ([0.1, -0.2], [80.0, -40.0]),
                ([-0.05, 0.15], [-24.0, 56.0]),
                ([0.0; 2], [0.0; 2]),
            ] {
                let field = vec![external; 16 * 14];
                flow.rebind_forcing(gravity, &field).unwrap();
                let prepared = flow.sample_at_si(8, 7).unwrap().unwrap();
                for d in 0..2 {
                    close(prepared.velocity_m_s[d], expected[d]);
                }
                flow.step();
                for d in 0..2 {
                    expected[d] += dt * (gravity[d] + external[d] / rho);
                }
                let stepped = flow.sample_at_si(8, 7).unwrap().unwrap();
                for d in 0..2 {
                    close(stepped.velocity_m_s[d], expected[d]);
                }
            }
            assert!(std::ptr::eq(flow.scaling().material(), &source));
            assert_eq!(flow.scaling().cell_size_m(), dx);
            assert_eq!(flow.scaling().timestep_s(), dt);
        }
    }
}

#[test]
fn changing_gravity_retains_si_force_and_evolved_interface_moments() {
    let gate = fs_exec::CancelGate::new_clock_free();
    for reference_density in [700.0, 1_400.0] {
        let source = liquid(reference_density, reference_density * 0.0001, None);
        let mut flow = plan(&source)
            .with_pressure_reference(100_000.0)
            .unwrap()
            .free_surface2(domain2(), |_, y| y < 8)
            .unwrap();
        // Real pressure evolution produces nonunit local density and physical
        // nonequilibrium stress; neither may be reset by a force-chart change.
        flow.step_with_gas_pressure(&gate, 0.05, |_| 100_000.0 + 0.002 * reference_density)
            .unwrap();
        let before = flow.lattice().clone();
        let mut density_departure = 0.0_f64;
        let mut stress_departure = 0.0_f64;
        for i in 0..before.grid.f.len() {
            if matches!(before.grid.flags[i], Cell::Fluid | Cell::Interface) {
                let m = physical(&before.grid, i);
                density_departure = density_departure.max((m[0] - 1.0).abs());
                stress_departure =
                    stress_departure.max((m[3] - m[0] * CS2 - m[1] * m[1] / m[0]).abs());
            }
        }
        assert!(density_departure > 1e-6);
        assert!(stress_departure > 1e-8);
        let external: Vec<_> = (0..16 * 14)
            .map(|i| [6.0 + (i % 3) as f64, -3.0 + (i % 2) as f64])
            .collect();
        flow.rebind_forcing([0.3, -0.2], &external).unwrap();
        let held_external = flow.lattice().grid.fext.clone();
        flow.rebind_gravity([-0.15, 0.25]).unwrap();
        assert_eq!(flow.lattice().grid.fext, held_external);
        unchanged_geometry(flow.lattice(), &before);
        let grid = &flow.lattice().grid;
        let mut post = Vec::new();
        grid.collide_into(&mut post);
        for i in 0..grid.f.len() {
            if !matches!(grid.flags[i], Cell::Fluid | Cell::Interface) {
                continue;
            }
            for (a, b) in physical(grid, i).into_iter().zip(physical(&before.grid, i)) {
                close(a, b);
            }
            let rho = grid.f[i].iter().sum::<f64>();
            for d in 0..2 {
                // Reference density converts N/m³; CURRENT local density
                // multiplies gravity. Substituting either one for the other fails.
                let converted_external = external[i][d] / reference_density * 0.001;
                close(grid.fext[i][d], converted_external);
                let expected = rho * [-0.15, 0.25][d] * 0.001 + converted_external;
                let actual = post[i]
                    .iter()
                    .zip(&grid.f[i])
                    .zip(E)
                    .map(|((after, before), e)| (after - before) * e[d])
                    .sum::<f64>();
                close(actual, expected);
            }
        }
        assert!(std::ptr::eq(flow.scaling().material(), &source));
    }
}

#[test]
fn invalid_si_fields_and_underflow_leave_liquid_and_held_loads_unchanged() {
    let source = liquid(1_000.0, 0.1, None);
    let mut flow = plan(&source)
        .free_surface2(domain2(), |_, y| y < 8)
        .unwrap();
    flow.rebind_forcing([0.1, -0.1], &vec![[2.0, -1.0]; 16 * 14])
        .unwrap();
    for case in 0..5 {
        let before = flow.lattice().clone();
        let mut external = vec![[3.0, -2.0]; 16 * 14];
        let mut gravity = [0.2, -0.3];
        match case {
            0 => {
                external.pop();
            }
            1 => external[16 * 14 - 1][0] = f64::NAN,
            2 => external[16 * 14 - 1][1] = f64::INFINITY,
            3 => external[16 * 14 - 1][0] = f64::from_bits(1),
            _ => gravity[0] = f64::NAN,
        }
        let rejected = flow.rebind_forcing(gravity, &external).unwrap_err();
        if case == 0 {
            assert!(matches!(
                rejected,
                MaterialLbmError::DenseForcing(ForceRebindError2::Invalid(_))
            ));
        }
        unchanged_geometry(flow.lattice(), &before);
        assert_eq!(flow.lattice().grid.f, before.grid.f);
        assert_eq!(flow.lattice().grid.g, before.grid.g);
        assert_eq!(flow.lattice().grid.fext, before.grid.fext);
    }
    for bad in [f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        let before = flow.lattice().clone();
        assert!(flow.rebind_gravity([bad, 0.0]).is_err());
        unchanged_geometry(flow.lattice(), &before);
        assert_eq!(flow.lattice().grid.f, before.grid.f);
        assert_eq!(flow.lattice().grid.g, before.grid.g);
        assert_eq!(flow.lattice().grid.fext, before.grid.fext);
    }
}
