//! Physical gravity on the actual density carried by sparse liquid populations.
use super::*;
use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::FreeSurfaceError3;
use fs_lbm::material::MaterialLbmError;
use fs_lbm::material::free_surface::MaterialFreeSurfaceInitialCell3;

fn close(actual: f64, expected: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= 4e-12 * expected.abs().max(1.0),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn filled(velocity_m_s: [f64; 3]) -> Vec<MaterialFreeSurfaceInitialCell3> {
    let mut cells = Vec::new();
    for z in 0..12 {
        for y in 0..12 {
            for x in 0..12 {
                cells.push(MaterialFreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: 1.0,
                    velocity_m_s,
                });
            }
        }
    }
    cells
}

#[test]
fn source_density_and_pressure_density_share_the_same_si_acceleration() {
    let gate = CancelGate::new_clock_free();
    let initial_velocity = [0.01, -0.015, 0.005];
    let initial = filled(initial_velocity);
    with_cx(&gate, |cx| {
        for reference_density in [700.0, 1400.0] {
            let source = liquid(reference_density, reference_density * 0.0001, None);
            for (dx, dt) in [(0.001, 0.001), (0.002, 0.0005)] {
                for numerical_density in [0.96, 1.04] {
                    let domain = MaterialFreeSurfaceDomain3 {
                        cells: [12; 3],
                        gravity_m_s2: [0.0; 3],
                        max_cells: 12usize.pow(3),
                        ..domain()
                    };
                    let velocity_unit: f64 = dx / dt;
                    let pressure = 100_000.0
                        + reference_density
                            * velocity_unit.powi(2)
                            * CS2
                            * (numerical_density - 1.0);
                    let mut flow = plan_free_surface_from_material(&source, dx, dt, 0.1)
                        .unwrap()
                        .with_density_dependent_gravity()
                        .with_pressure_reference(100_000.0)
                        .unwrap()
                        .free_surface3_from_cells_with_pressure(
                            domain,
                            &initial,
                            0.05,
                            |_| pressure,
                            cx,
                        )
                        .unwrap();
                    let mut expected = initial_velocity;
                    for gravity in [[0.2, -0.3, 0.1], [-0.1, 0.15, -0.05], [0.0; 3]] {
                        let before = flow.sample_at_si(6, 6, 6).unwrap().unwrap();
                        flow.rebind_gravity(gravity, &gate).unwrap();
                        let rebound = flow.sample_at_si(6, 6, 6).unwrap().unwrap();
                        assert_eq!(rebound.density_kg_m3, before.density_kg_m3);
                        for d in 0..3 {
                            close(rebound.velocity_m_s[d], expected[d]);
                        }
                        flow.step_with_gate(&gate).unwrap();
                        let sample = flow.sample_at_si(6, 6, 6).unwrap().unwrap();
                        for d in 0..3 {
                            expected[d] += gravity[d] * dt;
                            close(sample.velocity_m_s[d], expected[d]);
                        }
                        // This interior has not received a wall disturbance.
                        // Force/rho, not force/rho_reference, gives the result.
                        close(sample.density_kg_m3 / reference_density, numerical_density);
                    }
                    assert!(std::ptr::eq(flow.scaling().material(), &source));
                    assert_eq!(flow.scaling().cell_size_m(), dx);
                    assert_eq!(flow.scaling().timestep_s(), dt);
                }
            }
        }
    });
}

#[test]
fn source_pressure_and_interface_evolution_refresh_the_same_gravity_law() {
    let source = liquid(1000.0, 0.1, None);
    let gate = CancelGate::new_clock_free();
    let gravity = [0.02, -0.01, -0.03];
    let domain = MaterialFreeSurfaceDomain3 {
        cells: [12; 3],
        gravity_m_s2: gravity,
        max_cells: 12usize.pow(3),
        ..domain()
    };
    let mut flow = plan(&source)
        .with_density_dependent_gravity()
        .with_pressure_reference(100_000.0)
        .unwrap()
        .free_surface3(domain, |_, _, z| z < 6)
        .unwrap();
    let acceleration = gravity.map(|g| g * 0.001);
    assert_eq!(flow.lattice().body_acceleration(), Some(acceleration));
    for component in flow.sample_at_si(5, 5, 5).unwrap().unwrap().velocity_m_s {
        close(component, 0.0);
    }
    let mass = flow.ledger_mass_kg().unwrap();
    let mut density_departure = 0.0_f64;
    for pressure in [100_001.0, 100_002.0, 100_001.5] {
        flow.step_with_gas_pressure(&gate, 0.05, |_| pressure)
            .unwrap();
        for z in 0..12 {
            for y in 0..12 {
                for x in 0..12 {
                    let Some(sample) = flow.lattice().sample_at(x, y, z) else {
                        continue;
                    };
                    density_departure = density_departure.max((sample.density - 1.0).abs());
                    let force = flow.lattice().force_density_at(x, y, z).unwrap();
                    for d in 0..3 {
                        assert_eq!(force[d], sample.density * acceleration[d]);
                    }
                }
            }
        }
        near(flow.ledger_mass_kg().unwrap(), mass, 2e-12);
    }
    assert!(density_departure > 1e-4);
}

#[test]
fn replacing_si_gravity_preserves_accepted_state_and_refuses_atomically() {
    let source = liquid(1200.0, 0.12, None);
    let gate = CancelGate::new_clock_free();
    let mut flow = plan(&source).free_surface3(domain(), drop).unwrap();
    assert_eq!(flow.lattice().body_acceleration(), None);
    flow.step_with_gate(&gate).unwrap();
    let before = flow.clone();
    let gravity = [-0.2, 0.1, 0.3];
    flow.rebind_gravity(gravity, &gate).unwrap();
    assert_eq!(
        flow.ledger_mass_kg().unwrap(),
        before.ledger_mass_kg().unwrap()
    );
    assert_eq!(
        flow.carry_mass_kg().unwrap(),
        before.carry_mass_kg().unwrap()
    );
    for z in 0..16 {
        for y in 0..16 {
            for x in 0..16 {
                let old = before.sample_at_si(x, y, z).unwrap();
                let new = flow.sample_at_si(x, y, z).unwrap();
                match (old, new) {
                    (None, None) => {}
                    (Some(old), Some(new)) => {
                        assert_eq!(old.cell, new.cell);
                        assert_eq!(old.density_kg_m3, new.density_kg_m3);
                        assert_eq!(old.fill_fraction, new.fill_fraction);
                        assert_eq!(old.pressure_perturbation_pa, new.pressure_perturbation_pa);
                        for d in 0..3 {
                            close(old.velocity_m_s[d], new.velocity_m_s[d]);
                        }
                    }
                    _ => panic!("gravity rebinding changed liquid topology"),
                }
            }
        }
    }
    let bits = flow.lattice().grid().state_bits();
    let accepted_force = flow.lattice().body_acceleration();
    let mass = flow.ledger_mass_kg().unwrap();
    flow.rebind_gravity(gravity, &gate).unwrap();
    assert_eq!(flow.lattice().grid().state_bits(), bits);
    for refused in [
        [f64::NAN, 0.0, 0.0],
        [f64::MAX, f64::MAX, f64::MAX],
        [f64::from_bits(1), 0.0, 0.0],
    ] {
        assert!(flow.rebind_gravity(refused, &gate).is_err());
        assert_eq!(flow.lattice().grid().state_bits(), bits);
        assert_eq!(flow.lattice().body_acceleration(), accepted_force);
        assert_eq!(flow.ledger_mass_kg().unwrap(), mass);
    }
    gate.request();
    assert!(matches!(
        flow.rebind_gravity([0.0; 3], &gate),
        Err(MaterialLbmError::FreeSurface(FreeSurfaceError3::Cancelled))
    ));
    assert_eq!(flow.lattice().grid().state_bits(), bits);
    assert_eq!(flow.lattice().body_acceleration(), accepted_force);
    assert_eq!(flow.ledger_mass_kg().unwrap(), mass);
}

#[test]
fn zero_acceleration_retains_native_trajectory_and_source_scaling() {
    let source = liquid(1000.0, 0.1, None);
    let zero = MaterialFreeSurfaceDomain3 {
        gravity_m_s2: [0.0; 3],
        ..domain()
    };
    let mut reference = plan(&source).free_surface3(zero, drop).unwrap();
    let mut acceleration = plan(&source)
        .with_density_dependent_gravity()
        .free_surface3(zero, drop)
        .unwrap();
    assert_eq!(reference.lattice().body_acceleration(), None);
    assert_eq!(acceleration.lattice().body_acceleration(), Some([0.0; 3]));
    for _ in 0..3 {
        assert_eq!(
            reference.lattice().grid().state_bits(),
            acceleration.lattice().grid().state_bits()
        );
        assert_eq!(
            reference.ledger_mass_kg().unwrap(),
            acceleration.ledger_mass_kg().unwrap()
        );
        reference.step().unwrap();
        acceleration.step().unwrap();
    }
}
