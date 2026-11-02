//! G0/G3/G4 pressure load, gauge, material scaling and atomic-refusal checks.

use super::*;
use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::FreeSurfaceError3;
use fs_lbm::material::MaterialLbmError;
use fs_lbm::material::free_surface::MaterialFreeSurface3;
use fs_lbm::material::free_surface::MaterialFreeSurfaceInitialCell3;

fn initial_slab(axis: usize, positive_gravity: bool) -> Vec<MaterialFreeSurfaceInitialCell3> {
    let mut initial = Vec::new();
    for z in 0..12 {
        for y in 0..12 {
            for x in 0..12 {
                let coordinates = [x, y, z];
                let occupied = if positive_gravity {
                    coordinates[axis] >= 4
                } else {
                    coordinates[axis] < 8
                };
                if occupied {
                    initial.push(MaterialFreeSurfaceInitialCell3 {
                        coordinates,
                        fill_fraction: 1.0,
                        velocity_m_s: [0.0; 3],
                    });
                }
            }
        }
    }
    initial
}

#[test]
fn g1_initial_hydrostatic_pressure_balances_interior_gravity_in_all_directions() {
    let source = liquid(1_000.0, 0.1, None);
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for axis in 0..3 {
            for sign in [-1.0, 1.0] {
                let initial = initial_slab(axis, sign > 0.0);
                let mut gravity = [0.0; 3];
                gravity[axis] = sign * 0.5;
                let domain = MaterialFreeSurfaceDomain3 {
                    cells: [12; 3],
                    gravity_m_s2: gravity,
                    max_cells: 12usize.pow(3),
                    ..domain()
                };
                let surface = if sign > 0.0 { 0.004 } else { 0.008 };
                // Independent hydrostatics: grad(p) = rho_source * g.
                let pressure =
                    |x: [f64; 3]| 100_000.0 + 1_000.0 * gravity[axis] * (x[axis] - surface);
                let mut balanced = plan(&source)
                    .with_pressure_reference(100_000.0)
                    .unwrap()
                    .free_surface3_from_cells_with_pressure(domain, &initial, 0.05, pressure, cx)
                    .unwrap();
                let mut uniform = plan(&source)
                    .free_surface3_from_cells(domain, &initial, cx)
                    .unwrap();
                let mut expected_mass = 0.0;
                for cell in &initial {
                    let position = cell.coordinates.map(|c| (c as f64 + 0.5) * 0.001);
                    let sample = balanced
                        .sample_at_si(
                            cell.coordinates[0],
                            cell.coordinates[1],
                            cell.coordinates[2],
                        )
                        .unwrap()
                        .unwrap();
                    assert!(
                        (sample.pressure_perturbation_pa - (pressure(position) - 100_000.0)).abs()
                            < 2e-12
                    );
                    expected_mass +=
                        1e-6 * (1.0 + (pressure(position) - 100_000.0) / (1_000.0 * CS2));
                }
                near(balanced.ledger_mass_kg().unwrap(), expected_mass, 5e-14);
                balanced.step().unwrap();
                uniform.step().unwrap();
                let c = if sign > 0.0 { 6 } else { 5 };
                let prepared = balanced.sample_at_si(c, c, c).unwrap().unwrap();
                let unprepared = uniform.sample_at_si(c, c, c).unwrap().unwrap();
                assert!((unprepared.velocity_m_s[axis] - gravity[axis] * 0.001).abs() < 1e-12);
                // Equilibrium populations need not be an exact steady discrete
                // hydrostatic state, but the spurious first-order acceleration
                // must disappear. This checks cells before boundary waves arrive.
                assert!(
                    prepared.velocity_m_s[axis].abs() < 1e-4 * unprepared.velocity_m_s[axis].abs(),
                    "axis={axis}, sign={sign}: {:?} vs {:?}",
                    prepared.velocity_m_s,
                    unprepared.velocity_m_s
                );
            }
        }
    });
}

#[test]
fn g1_initial_pressure_matches_a_loaded_flat_surface_without_startup_flow() {
    let source = liquid(1_000.0, 0.1, None);
    let initial = initial_slab(2, false);
    let domain = MaterialFreeSurfaceDomain3 {
        cells: [12; 3],
        gravity_m_s2: [0.0; 3],
        max_cells: 12usize.pow(3),
        ..domain()
    };
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let make = |reference, pressure| {
            plan(&source)
                .with_pressure_reference(reference)
                .unwrap()
                .free_surface3_from_cells_with_pressure(domain, &initial, 0.05, |_| pressure, cx)
                .unwrap()
        };
        let mut flow = make(100_000.0, 100_002.0);
        let mut shifted = make(100_500.0, 100_502.0);
        let mass = flow.ledger_mass_kg().unwrap();
        near(mass, initial.len() as f64 * 1e-6 * 1.006, 2e-14);
        for _ in 0..5 {
            flow.step_with_gas_pressure(&gate, 0.05, |_| 100_002.0)
                .unwrap();
            shifted
                .step_with_gas_pressure(&gate, 0.05, |_| 100_502.0)
                .unwrap();
        }
        assert_eq!(
            flow.lattice().grid().state_bits(),
            shifted.lattice().grid().state_bits()
        );
        near(flow.ledger_mass_kg().unwrap(), mass, 2e-13);
        for cell in &initial {
            let [x, y, z] = cell.coordinates;
            let sample = flow.sample_at_si(x, y, z).unwrap().unwrap();
            assert!(sample.velocity_m_s.iter().all(|u| u.abs() < 2e-14));
            assert!((sample.pressure_perturbation_pa - 2.0).abs() < 2e-12);
        }
        let plain = plan(&source)
            .free_surface3_from_cells(domain, &initial, cx)
            .unwrap();
        let zero = make(100_000.0, 100_000.0);
        assert_eq!(
            plain.lattice().grid().state_bits(),
            zero.lattice().grid().state_bits()
        );
        assert_eq!(
            plain.ledger_mass_kg().unwrap().to_bits(),
            zero.ledger_mass_kg().unwrap().to_bits()
        );
    });
}

#[test]
fn g1_initial_spherical_laplace_pressure_uses_physical_surface_tension() {
    let source = liquid(1_000.0, 0.1, Some(0.005));
    let radius_m = 0.0035;
    let laplace_pa = 2.0 * 0.005 / radius_m;
    let mut initial = Vec::new();
    for z in 0..16 {
        for y in 0..16 {
            for x in 0..16 {
                if drop(x, y, z) {
                    initial.push(MaterialFreeSurfaceInitialCell3 {
                        coordinates: [x, y, z],
                        fill_fraction: 1.0,
                        velocity_m_s: [0.0; 3],
                    });
                }
            }
        }
    }
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        let flow = plan(&source)
            .with_pressure_reference(100_000.0)
            .unwrap()
            .free_surface3_from_cells_with_pressure(
                MaterialFreeSurfaceDomain3 {
                    gravity_m_s2: [0.0; 3],
                    ..domain()
                },
                &initial,
                0.05,
                |_| 100_000.0 + laplace_pa,
                cx,
            )
            .unwrap();
        let sample = flow.sample_at_si(7, 7, 7).unwrap().unwrap();
        assert!((sample.pressure_perturbation_pa - laplace_pa).abs() < 1e-11);
        assert!(sample.velocity_m_s.iter().all(|u| u.abs() < 1e-15));
        // The supplied continuum jump is not added a second time at ingress.
        // Staircase occupancy here is not an exact static-drop equilibrium.
    });
}

#[test]
fn g0_initial_pressure_refuses_unsupported_or_cancelled_candidates() {
    let source = liquid(1_000.0, 0.1, None);
    let initial = [MaterialFreeSurfaceInitialCell3 {
        coordinates: [6, 6, 6],
        fill_fraction: 0.7,
        velocity_m_s: [0.0; 3],
    }];
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        assert!(
            plan(&source)
                .free_surface3_from_cells_with_pressure(domain(), &initial, 0.1, |_| 100_000.0, cx,)
                .is_err()
        );
        for (bound, pressure) in [
            (0.0, 100_000.0),
            (1.0, 100_000.0),
            (f64::NAN, 100_000.0),
            (0.05, 0.0),
            (0.05, -1.0),
            (0.05, f64::NAN),
            (0.05, f64::INFINITY),
            (0.01, 100_004.0),
            (0.01, 99_996.0),
        ] {
            assert!(
                plan(&source)
                    .with_pressure_reference(100_000.0)
                    .unwrap()
                    .free_surface3_from_cells_with_pressure(
                        domain(),
                        &initial,
                        bound,
                        |_| pressure,
                        cx
                    )
                    .is_err()
            );
        }
    });
    let cancelled = CancelGate::new_clock_free();
    with_cx(&cancelled, |cx| {
        let result = plan(&source)
            .with_pressure_reference(100_000.0)
            .unwrap()
            .free_surface3_from_cells_with_pressure(
                domain(),
                &initial,
                0.05,
                |_| {
                    cancelled.request();
                    100_001.0
                },
                cx,
            );
        assert!(matches!(
            result,
            Err(MaterialLbmError::FreeSurface(FreeSurfaceError3::Cancelled))
        ));
    });
}

fn slab<'a>(
    source: &'a ResolvedLiquidSurfaceState,
    axis: usize,
    reference: f64,
) -> MaterialFreeSurface3<'a> {
    plan(source)
        .with_pressure_reference(reference)
        .unwrap()
        .free_surface3(
            MaterialFreeSurfaceDomain3 {
                cells: [12; 3],
                gravity_m_s2: [0.0; 3],
                max_cells: 12 * 12 * 12,
                ..domain()
            },
            |x, y, z| (4..8).contains(&[x, y, z][axis]),
        )
        .unwrap()
}

fn mean_velocity(flow: &MaterialFreeSurface3<'_>) -> [f64; 3] {
    let mut momentum = [0.0; 3];
    let mut mass = 0.0;
    for z in 0..12 {
        for y in 0..12 {
            for x in 0..12 {
                if let Some(sample) = flow.sample_at_si(x, y, z).unwrap() {
                    let weight = sample.fill_fraction * sample.density_kg_m3;
                    mass += weight;
                    for (p, u) in momentum.iter_mut().zip(sample.velocity_m_s) {
                        *p += weight * u;
                    }
                }
            }
        }
    }
    momentum.map(|p| p / mass)
}

#[test]
fn zero_gas_pressure_load_preserves_original_capillary_trajectory_bits() {
    let source = liquid(1_000.0, 0.1, Some(0.005));
    let mut original = plan(&source).free_surface3(domain(), drop).unwrap();
    let mut loaded = plan(&source)
        .with_pressure_reference(100_000.0)
        .unwrap()
        .free_surface3(domain(), drop)
        .unwrap();
    let gate = CancelGate::new_clock_free();
    for _ in 0..6 {
        original.step().unwrap();
        loaded
            .step_with_gas_pressure(&gate, 0.1, |_| 100_000.0)
            .unwrap();
    }
    assert_eq!(
        original.lattice().grid().state_bits(),
        loaded.lattice().grid().state_bits()
    );
    assert_eq!(
        original.ledger_mass_kg().unwrap().to_bits(),
        loaded.ledger_mass_kg().unwrap().to_bits()
    );
    assert_eq!(
        original.lattice().conversions(),
        loaded.lattice().conversions()
    );
}

#[test]
fn pressure_difference_moves_liquid_with_axis_and_sign_symmetry_and_conserves_mass() {
    let source = liquid(1_000.0, 0.1, None);
    let gate = CancelGate::new_clock_free();
    let mut reference_speed: Option<f64> = None;
    for axis in 0..3 {
        for sign in [-1.0, 1.0] {
            let mut flow = slab(&source, axis, 100_000.0);
            let initial_mass = flow.ledger_mass_kg().unwrap();
            for _ in 0..10 {
                flow.step_with_gas_pressure(&gate, 0.05, |p| {
                    100_000.0
                        + if p[axis] < 0.006 {
                            sign * 0.75
                        } else {
                            -sign * 0.75
                        }
                })
                .unwrap();
            }
            near(flow.ledger_mass_kg().unwrap(), initial_mass, 1e-12);
            let velocity = mean_velocity(&flow);
            let speed = sign * velocity[axis];
            assert!(
                speed > 1e-4,
                "pressure gradient must move liquid toward low pressure: {velocity:?}"
            );
            for (other, component) in velocity.into_iter().enumerate() {
                if other != axis {
                    assert!(
                        component.abs() < 2e-13,
                        "spurious transverse flow: {velocity:?}"
                    );
                }
            }
            if let Some(reference) = reference_speed {
                near(speed, reference, 1e-9);
            } else {
                reference_speed = Some(speed);
            }
        }
    }
}

#[test]
fn absolute_pressure_reference_shifts_leave_dynamics_exactly_unchanged() {
    let source = liquid(1_000.0, 0.1, None);
    let gate = CancelGate::new_clock_free();
    let mut first = slab(&source, 0, 100_000.0);
    let mut shifted = slab(&source, 0, 100_500.0);
    for _ in 0..6 {
        first
            .step_with_gas_pressure(&gate, 0.05, |p| {
                100_000.0 + if p[0] < 0.006 { 0.75 } else { -0.75 }
            })
            .unwrap();
        shifted
            .step_with_gas_pressure(&gate, 0.05, |p| {
                100_500.0 + if p[0] < 0.006 { 0.75 } else { -0.75 }
            })
            .unwrap();
    }
    assert_eq!(
        first.lattice().grid().state_bits(),
        shifted.lattice().grid().state_bits()
    );
    assert_eq!(
        first.ledger_mass_kg().unwrap().to_bits(),
        shifted.ledger_mass_kg().unwrap().to_bits()
    );
    let sample = first.lattice().sample_at(4, 6, 6).unwrap();
    let pressure = first
        .scaling()
        .absolute_pressure_si(sample.density)
        .unwrap();
    let shifted_pressure = shifted
        .scaling()
        .absolute_pressure_si(sample.density)
        .unwrap();
    assert!((shifted_pressure - pressure - 500.0).abs() < 2e-11);
}

#[test]
fn pressure_and_material_scaling_use_the_original_density_and_timestep() {
    let source = liquid(1_000.0, 0.1, None);
    let denser = liquid(2_000.0, 0.2, None);
    let gate = CancelGate::new_clock_free();
    let mut first = slab(&source, 0, 100_000.0);
    let mut second = slab(&denser, 0, 100_000.0);
    for _ in 0..5 {
        first
            .step_with_gas_pressure(&gate, 0.05, |p| {
                100_000.0 + if p[0] < 0.006 { 0.75 } else { -0.75 }
            })
            .unwrap();
        second
            .step_with_gas_pressure(&gate, 0.05, |p| {
                100_000.0 + if p[0] < 0.006 { 1.5 } else { -1.5 }
            })
            .unwrap();
    }
    // Same kinematic viscosity and twice the physical load at twice the
    // density produce the same dimensionless equations and populations.
    assert_eq!(
        first.lattice().grid().state_bits(),
        second.lattice().grid().state_bits()
    );
    near(
        second.ledger_mass_kg().unwrap(),
        2.0 * first.ledger_mass_kg().unwrap(),
        2e-14,
    );
    near(
        first.scaling().absolute_pressure_si(1.003).unwrap(),
        100_001.0,
        2e-15,
    );
    assert!(plan(&source).absolute_pressure_si(1.0).is_err());
    for pressure in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(plan(&source).with_pressure_reference(pressure).is_err());
    }
}

#[test]
fn invalid_or_cancelled_pressure_loading_preserves_state_and_retry() {
    let source = liquid(1_000.0, 0.1, Some(0.005));
    let mut flow = plan(&source)
        .with_pressure_reference(100_000.0)
        .unwrap()
        .free_surface3(domain(), drop)
        .unwrap();
    let initial_bits = flow.lattice().grid().state_bits();
    let initial_mass = flow.ledger_mass_kg().unwrap().to_bits();
    let initial_conversions = flow.lattice().conversions();
    let gate = CancelGate::new_clock_free();
    for (bound, pressure) in [
        (0.1, f64::NAN),
        (0.1, f64::INFINITY),
        (0.1, 0.0),
        (0.1, -1.0),
        (0.1, 200_000.0),
        (0.1, 1.0),
        (0.0, 100_000.0),
        (1.0, 100_000.0),
        (f64::NAN, 100_000.0),
        // Zero imposed perturbation still has to admit the capillary load.
        (1e-8, 100_000.0),
    ] {
        assert!(
            flow.step_with_gas_pressure(&gate, bound, |_| pressure)
                .is_err()
        );
        assert_eq!(flow.lattice().grid().state_bits(), initial_bits);
        assert_eq!(flow.ledger_mass_kg().unwrap().to_bits(), initial_mass);
        assert_eq!(flow.lattice().conversions(), initial_conversions);
    }
    let cancelled = CancelGate::new_clock_free();
    let result = flow.step_with_gas_pressure(&cancelled, 0.1, |_| {
        cancelled.request();
        100_001.0
    });
    assert!(matches!(
        result,
        Err(MaterialLbmError::FreeSurface(FreeSurfaceError3::Cancelled))
    ));
    assert_eq!(flow.lattice().grid().state_bits(), initial_bits);
    assert_eq!(flow.ledger_mass_kg().unwrap().to_bits(), initial_mass);
    assert_eq!(flow.lattice().conversions(), initial_conversions);
    let mut expected = flow.clone();
    expected.step().unwrap();
    flow.step_with_gas_pressure(&gate, 0.1, |_| 100_000.0)
        .unwrap();
    assert_eq!(
        flow.lattice().grid().state_bits(),
        expected.lattice().grid().state_bits()
    );
}
