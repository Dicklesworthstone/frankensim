//! Physical gravity on the actual density carried by sparse liquid populations.
use super::*;
use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::FreeSurfaceError3;
use fs_lbm::material::MaterialLbmError;
use fs_lbm::material::free_surface::{
    MaterialAffineAcceleration3, MaterialFreeSurfaceInitialCell3,
};

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

#[test]
fn acceleration_endpoints_integrate_si_velocity_across_material_density_and_units() {
    let gate = CancelGate::new_clock_free();
    let initial_velocity = [0.01, -0.015, 0.005];
    let initial = filled(initial_velocity);
    let initial_gravity = [0.2, -0.3, 0.1];
    with_cx(&gate, |cx| {
        for reference_density in [700.0, 1400.0] {
            let source = liquid(reference_density, reference_density * 0.0001, None);
            for (dx, dt) in [(0.001, 0.001), (0.002, 0.0005)] {
                let velocity_unit: f64 = dx / dt;
                let numerical_density = 1.04;
                let pressure = 100_000.0
                    + reference_density * velocity_unit.powi(2) * CS2 * (numerical_density - 1.0);
                let mut flow = plan_free_surface_from_material(&source, dx, dt, 0.1)
                    .unwrap()
                    .with_density_dependent_gravity()
                    .with_pressure_reference(100_000.0)
                    .unwrap()
                    .free_surface3_from_cells_with_pressure(
                        MaterialFreeSurfaceDomain3 {
                            cells: [12; 3],
                            gravity_m_s2: initial_gravity,
                            max_cells: 12usize.pow(3),
                            ..domain()
                        },
                        &initial,
                        0.05,
                        |_| pressure,
                        cx,
                    )
                    .unwrap();
                let mut expected = initial_velocity;
                let mut previous = initial_gravity;
                for next in [[-0.1, 0.15, -0.05], [0.0; 3], [0.3, -0.1, 0.2]] {
                    flow.advance_gravity(next, &gate).unwrap();
                    let sample = flow.sample_at_si(6, 6, 6).unwrap().unwrap();
                    for axis in 0..3 {
                        expected[axis] += 0.5 * (previous[axis] + next[axis]) * dt;
                        close(sample.velocity_m_s[axis], expected[axis]);
                    }
                    close(sample.density_kg_m3 / reference_density, numerical_density);
                    previous = next;
                }
                assert!(std::ptr::eq(flow.scaling().material(), &source));
                assert_eq!(flow.scaling().reference_pressure_pa(), Some(100_000.0));
                assert_eq!(flow.scaling().cell_size_m(), dx);
                assert_eq!(flow.scaling().timestep_s(), dt);
            }
        }
    });
}

#[test]
fn changing_si_acceleration_and_gas_pressure_publish_one_density_force_endpoint() {
    let source = liquid(1000.0, 0.1, None);
    let gate = CancelGate::new_clock_free();
    let mut flow = plan(&source)
        .with_density_dependent_gravity()
        .with_pressure_reference(100_000.0)
        .unwrap()
        .free_surface3(
            MaterialFreeSurfaceDomain3 {
                cells: [12; 3],
                gravity_m_s2: [0.02, -0.01, -0.03],
                max_cells: 12usize.pow(3),
                ..domain()
            },
            |_, _, z| z < 6,
        )
        .unwrap();
    let mass = flow.ledger_mass_kg().unwrap();
    for (next, pressure) in [
        ([-0.03, 0.02, 0.01], 100_001.0),
        ([0.01, -0.02, -0.04], 100_002.0),
    ] {
        flow.advance_gravity_with_gas_link_pressure(next, &gate, 0.05, |_| pressure)
            .unwrap();
        let native = next.map(|g| flow.scaling().gravity_lattice(g).unwrap());
        assert_eq!(flow.lattice().body_acceleration(), Some(native));
        for z in 0..12 {
            for y in 0..12 {
                for x in 0..12 {
                    if let Some(sample) = flow.lattice().sample_at(x, y, z) {
                        assert_eq!(
                            flow.lattice().force_density_at(x, y, z).unwrap(),
                            native.map(|a| sample.density * a)
                        );
                    }
                }
            }
        }
        near(flow.ledger_mass_kg().unwrap(), mass, 2e-12);
    }
    let before = flow.clone();
    assert!(
        flow.advance_gravity_with_gas_link_pressure([0.0; 3], &gate, 0.05, |_| f64::NAN)
            .is_err()
    );
    assert_eq!(
        flow.lattice().grid().state_bits(),
        before.lattice().grid().state_bits()
    );
    assert_eq!(
        flow.lattice().body_acceleration(),
        before.lattice().body_acceleration()
    );
    assert_eq!(
        flow.ledger_mass_kg().unwrap(),
        before.ledger_mass_kg().unwrap()
    );
}

fn acceleration_si(field: MaterialAffineAcceleration3, position_m: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|axis| {
        field.acceleration_m_s2[axis]
            + (0..3)
                .map(|coordinate| {
                    field.gradient_per_s2[axis][coordinate]
                        * (position_m[coordinate] - field.origin_m[coordinate])
                })
                .sum::<f64>()
    })
}

#[test]
fn affine_si_field_converts_gradients_and_integrates_density_weighted_endpoints() {
    let gate = CancelGate::new_clock_free();
    let old = MaterialAffineAcceleration3 {
        origin_m: [0.003, -0.002, 0.001],
        acceleration_m_s2: [0.1, -0.2, 0.05],
        gradient_per_s2: [[10.0, 2.0, 0.0], [-3.0, -4.0, 1.0], [1.0, 0.0, 2.0]],
    };
    let next = MaterialAffineAcceleration3 {
        origin_m: [-0.001, 0.004, 0.003],
        acceleration_m_s2: [-0.2, 0.1, -0.05],
        gradient_per_s2: [[-5.0, 3.0, 1.0], [2.0, 7.0, 0.0], [0.0, -2.0, 4.0]],
    };
    for reference_density in [700.0, 1400.0] {
        let source = liquid(reference_density, reference_density * 0.0001, None);
        for (dx, dt) in [(0.001, 0.001), (0.002, 0.0005)] {
            let mut flow = plan_free_surface_from_material(&source, dx, dt, 0.1)
                .unwrap()
                .free_surface3(
                    MaterialFreeSurfaceDomain3 {
                        cells: [12; 3],
                        gravity_m_s2: [0.0; 3],
                        max_cells: 12usize.pow(3),
                        ..domain()
                    },
                    |_, _, _| true,
                )
                .unwrap();
            flow.rebind_affine_acceleration(old, &gate).unwrap();
            let positions = [[5, 5, 5], [7, 6, 4]];
            let before = positions.map(|p| flow.sample_at_si(p[0], p[1], p[2]).unwrap().unwrap());
            flow.advance_affine_acceleration(next, &gate).unwrap();
            for (xyz, incoming) in positions.into_iter().zip(before) {
                let position_m = xyz.map(|i| (i as f64 + 0.5) * dx);
                let a0 = acceleration_si(old, position_m);
                let a1 = acceleration_si(next, position_m);
                let sample = flow.sample_at_si(xyz[0], xyz[1], xyz[2]).unwrap().unwrap();
                let force = flow
                    .lattice()
                    .force_density_at(xyz[0], xyz[1], xyz[2])
                    .unwrap();
                // Independently dimensionalized at physical cell centres:
                // div(a_old)=8/s², and Guo's first step from rest integrates
                // dt*(rho_old*a_old + rho_next*a_next)/2 in the undisturbed bulk.
                close(
                    sample.density_kg_m3,
                    incoming.density_kg_m3 * (1.0 - 0.5 * 8.0 * dt * dt),
                );
                for axis in 0..3 {
                    close(incoming.velocity_m_s[axis], 0.0);
                    close(
                        force[axis] * reference_density * dx / (dt * dt),
                        sample.density_kg_m3 * a1[axis],
                    );
                    close(
                        sample.density_kg_m3 * sample.velocity_m_s[axis],
                        0.5 * dt
                            * (incoming.density_kg_m3 * a0[axis] + sample.density_kg_m3 * a1[axis]),
                    );
                }
            }
            assert!(std::ptr::eq(flow.scaling().material(), &source));
            assert!(flow.lattice().body_acceleration().is_none());
            assert!(flow.lattice().body_acceleration_field().is_some());
        }
    }
}

#[test]
fn affine_si_gas_pressure_keeps_spatial_force_and_rollback_consistent() {
    let source = liquid(1000.0, 0.1, None);
    let gate = CancelGate::new_clock_free();
    let field = MaterialAffineAcceleration3 {
        origin_m: [0.002, 0.003, 0.004],
        acceleration_m_s2: [0.02, -0.01, -0.03],
        gradient_per_s2: [[2.0, 0.0, 0.0], [0.0, -3.0, 0.0], [0.0, 0.0, 4.0]],
    };
    let mut flow = plan(&source)
        .with_pressure_reference(100_000.0)
        .unwrap()
        .free_surface3(
            MaterialFreeSurfaceDomain3 {
                cells: [12; 3],
                gravity_m_s2: [0.0; 3],
                max_cells: 12usize.pow(3),
                ..domain()
            },
            |_, _, z| z < 6,
        )
        .unwrap();
    flow.rebind_affine_acceleration(field, &gate).unwrap();
    let mass = flow.ledger_mass_kg().unwrap();
    let next = MaterialAffineAcceleration3 {
        acceleration_m_s2: [-0.03, 0.02, 0.01],
        ..field
    };
    flow.advance_affine_acceleration_with_gas_link_pressure(next, &gate, 0.05, |_| 100_002.0)
        .unwrap();
    let dx = flow.scaling().cell_size_m();
    let dt = flow.scaling().timestep_s();
    let mut density_departure = 0.0_f64;
    for z in 0..12 {
        for y in 0..12 {
            for x in 0..12 {
                if let Some(sample) = flow.sample_at_si(x, y, z).unwrap() {
                    let position = [x, y, z].map(|i| (i as f64 + 0.5) * dx);
                    let expected = acceleration_si(next, position);
                    let force = flow.lattice().force_density_at(x, y, z).unwrap();
                    for axis in 0..3 {
                        close(
                            force[axis] * 1000.0 * dx / (dt * dt),
                            sample.density_kg_m3 * expected[axis],
                        );
                    }
                    density_departure =
                        density_departure.max((sample.density_kg_m3 / 1000.0 - 1.0).abs());
                }
            }
        }
    }
    assert!(density_departure > 1e-4);
    near(flow.ledger_mass_kg().unwrap(), mass, 2e-12);
    let before = flow.clone();
    assert!(
        flow.advance_affine_acceleration_with_gas_link_pressure(field, &gate, 0.05, |_| f64::NAN)
            .is_err()
    );
    assert_eq!(
        flow.lattice().grid().state_bits(),
        before.lattice().grid().state_bits()
    );
    assert_eq!(
        flow.lattice().body_acceleration_field(),
        before.lattice().body_acceleration_field()
    );
    assert_eq!(
        flow.ledger_mass_kg().unwrap(),
        before.ledger_mass_kg().unwrap()
    );
}
