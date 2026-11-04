use super::*;

fn surface(initial: &[FreeSurfaceInitialCell3], force: [f64; 3]) -> FreeSurface3 {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 1941,
                kernel_id: 19,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        FreeSurface3::from_cells(
            [8; 3],
            0.8,
            force,
            0.0,
            ContactModel3::Neutral,
            &[],
            initial,
            &cx,
        )
        .unwrap()
    })
}

fn cell(coordinates: [usize; 3], density: f64, velocity: [f64; 3]) -> FreeSurfaceInitialCell3 {
    FreeSurfaceInitialCell3 {
        coordinates,
        density,
        velocity,
        fill_fraction: 0.5,
    }
}

fn unchanged(before: &FreeSurface3, after: &FreeSurface3) {
    assert_eq!(before.body_acceleration, after.body_acceleration);
    assert_eq!(before.local_force, after.local_force);
    assert_eq!(before.grid.active_keys(), after.grid.active_keys());
    for (&key, tile) in &before.tiles {
        assert_eq!(tile.cells, after.tiles[&key].cells);
        assert_eq!(
            tile.mass.map(f64::to_bits),
            after.tiles[&key].mass.map(f64::to_bits)
        );
        let old = before.grid.slot_of(key).unwrap();
        let new = after.grid.slot_of(key).unwrap();
        for lane in 0..TILE_CELLS {
            assert_eq!(
                before.grid.populations(old, lane).map(f64::to_bits),
                after.grid.populations(new, lane).map(f64::to_bits)
            );
        }
    }
    assert_eq!(before.carry, after.carry);
    assert_eq!(before.conversions, after.conversions);
    assert_eq!(before.conversion_momentum, after.conversion_momentum);
}

fn moments(surface: &FreeSurface3, xyz: [usize; 3]) -> ([f64; 3], [[f64; 3]; 3]) {
    let sample = surface.sample_at(xyz[0], xyz[1], xyz[2]).unwrap();
    let force = surface.force_density_at(xyz[0], xyz[1], xyz[2]).unwrap();
    let (slot, lane) = surface
        .grid
        .resolve_source(xyz[0] as i64, xyz[1] as i64, xyz[2] as i64)
        .unwrap();
    let mut momentum = force.map(|value| 0.5 * value);
    let mut stress = core::array::from_fn(|a| {
        core::array::from_fn(|b| {
            0.5 * (sample.velocity[a] * force[b] + sample.velocity[b] * force[a])
        })
    });
    for (population, direction) in surface.grid.populations(slot, lane).into_iter().zip(E3) {
        let e = [direction.0, direction.1, direction.2].map(f64::from);
        for a in 0..3 {
            momentum[a] += population * e[a];
            for b in 0..3 {
                stress[a][b] += population * e[a] * e[b];
            }
        }
    }
    (momentum, stress)
}

#[test]
fn g1_body_rebind_preserves_density_mass_and_nonequilibrium_moments() {
    let cells = [
        cell([2, 2, 2], 0.6, [0.02, -0.01, 0.005]),
        cell([5, 5, 5], 1.4, [-0.01, 0.02, -0.005]),
    ];
    let gate = CancelGate::new_clock_free();
    let mut flow = surface(&cells, [0.0002, -0.0001, 0.0003]);
    for cell in cells {
        flow.initialize_strain_at(
            cell.coordinates,
            [0.001, -0.001, 0.0, 0.002, 0.0, -0.001],
            0.8,
        )
        .unwrap();
    }
    for acceleration in [[0.002, -0.003, 0.001], [-0.001, 0.002, -0.004], [0.0; 3]] {
        let before = flow.clone();
        let mass = flow.ledger_mass();
        flow.rebind_body_acceleration(acceleration, &gate).unwrap();
        assert_eq!(flow.body_acceleration(), Some(acceleration));
        assert_eq!(flow.ledger_mass().to_bits(), mass.to_bits());
        for cell in cells {
            let [x, y, z] = cell.coordinates;
            let old = before.sample_at(x, y, z).unwrap();
            let next = flow.sample_at(x, y, z).unwrap();
            assert_eq!(old.density.to_bits(), next.density.to_bits());
            assert_eq!(old.fill_fraction.to_bits(), next.fill_fraction.to_bits());
            assert_eq!(
                flow.force_density_at(x, y, z).unwrap(),
                acceleration.map(|a| next.density * a)
            );
            let (p0, s0) = moments(&before, cell.coordinates);
            let (p1, s1) = moments(&flow, cell.coordinates);
            for a in 0..3 {
                assert!((p1[a] - p0[a]).abs() < 2e-14 * next.density);
                for b in 0..3 {
                    assert!((s1[a][b] - s0[a][b]).abs() < 2e-14 * next.density);
                }
            }
        }
        let before = flow.clone();
        flow.rebind_body_acceleration(acceleration, &gate).unwrap();
        unchanged(&before, &flow);
    }
    assert!(
        flow.local_force.is_some(),
        "zero acceleration retains explicit ownership/history"
    );
}

#[test]
fn g0_body_force_ownership_and_late_cancellation_refuse_atomically() {
    let gate = CancelGate::new_clock_free();
    let mut flow = surface(&[cell([3, 3, 3], 1.2, [0.01, 0.0, 0.0])], [0.0; 3]);
    let before = flow.clone();
    for acceleration in [[f64::NAN, 0.0, 0.0], [10.0, 0.0, 0.0]] {
        assert!(flow.rebind_body_acceleration(acceleration, &gate).is_err());
        unchanged(&before, &flow);
    }
    let mut checkpoints = 0;
    let before_publication = 2 * flow.tiles.len() + 3;
    assert!(matches!(
        flow.rebind_body_acceleration_checked([0.001, 0.0, 0.0], || {
            checkpoints += 1;
            if checkpoints == before_publication {
                Err(FreeSurfaceError3::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(FreeSurfaceError3::Cancelled)
    ));
    assert_eq!(
        checkpoints, before_publication,
        "cancellation after all cell replacements were prepared"
    );
    unchanged(&before, &flow);

    let (slot, lane) = flow.grid.resolve_source(3, 3, 3).unwrap();
    let key = flow.grid.active_keys()[slot];
    let forces = BTreeMap::from([((key, lane), [0.002, 0.0, 0.0])]);
    flow.rebind_force_field(forces.clone(), &gate).unwrap();
    let independent = flow.clone();
    assert!(
        flow.rebind_body_acceleration([0.001, 0.0, 0.0], &gate)
            .is_err()
    );
    unchanged(&independent, &flow);

    let mut flow = before;
    flow.rebind_body_acceleration([0.001, 0.0, 0.0], &gate)
        .unwrap();
    let accepted = flow.clone();
    assert!(flow.rebind_force_field(forces, &gate).is_err());
    unchanged(&accepted, &flow);
    let mut legacy = flow.with_conversion_redistribution(ConversionRedistribution3::LegacyGlobal);
    let before = legacy.clone();
    assert!(legacy.step().is_err());
    unchanged(&before, &legacy);
}

#[test]
fn g1_new_interfaces_use_their_own_density_acceleration_frame() {
    let velocity = [0.02, -0.01, 0.03];
    let mut cells = Vec::new();
    for z in 2..5 {
        for y in 2..5 {
            for x in 2..5 {
                cells.push(FreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: 1.0,
                    density: 0.7 + 0.1 * ((x + 2 * y + z) % 7) as f64,
                    velocity,
                });
            }
        }
    }
    let mut flow = surface(&cells, [0.0; 3]);
    let acceleration = [0.002, -0.003, 0.001];
    flow.rebind_body_acceleration(acceleration, &CancelGate::new_clock_free())
        .unwrap();
    let (slot, lane) = flow.grid.resolve_source(4, 3, 3).unwrap();
    let key = flow.grid.active_keys()[slot];
    let density = flow.grid.populations(slot, lane).iter().sum::<f64>();
    flow.tiles.get_mut(&key).unwrap().mass[lane] = 1.1 * density;
    let mass = flow.ledger_mass();
    flow.apply_conversions(None).unwrap();
    assert!(flow.conversions.gas_to_interface > 0);
    let sample = flow.sample_at(5, 3, 3).unwrap();
    assert_eq!(sample.cell, Cell3::Interface);
    assert_eq!(
        flow.force_density_at(5, 3, 3).unwrap(),
        acceleration.map(|a| sample.density * a)
    );
    for a in 0..3 {
        assert!((sample.velocity[a] - velocity[a]).abs() < 2e-14);
    }
    let (_, stress) = moments(&flow, [5, 3, 3]);
    for a in 0..3 {
        for b in 0..3 {
            let expected =
                sample.density * (velocity[a] * velocity[b] + if a == b { CS2_3 } else { 0.0 });
            assert!((stress[a][b] - expected).abs() < 2e-14 * sample.density);
        }
    }
    assert!((flow.ledger_mass() - mass).abs() < 2e-13 * mass);
}

#[test]
fn g1_gas_pressure_pairs_close_in_the_actual_endpoint_force_frame() {
    let mut flow = surface(&[cell([3, 3, 3], 0.75, [0.02, -0.01, 0.015])], [0.0; 3]);
    let acceleration = [0.002, -0.003, 0.001];
    let gate = CancelGate::new_clock_free();
    flow.rebind_body_acceleration(acceleration, &gate).unwrap();
    let before = flow.sample_at(3, 3, 3).unwrap();
    let (slot, lane) = flow.grid.resolve_source(3, 3, 3).unwrap();
    let force = flow.force_density_at(3, 3, 3).unwrap();
    let post = crate::d3q19::collision_kernel3()(
        flow.grid.populations(slot, lane),
        flow.collision_model(),
        force,
    )
    .unwrap();
    let boundary_density = |q: usize| 1.05 + 0.01 * f64::from(E3[q].0 - E3[q].1);
    flow.step_with_gas_link_pressure(&gate, 0.2, |link| {
        Ok((boundary_density(link.incoming_direction) - 1.0) * CS2_3)
    })
    .unwrap();
    let after = flow.sample_at(3, 3, 3).unwrap();
    assert!((after.density - before.density).abs() > 0.1);
    let force = flow.force_density_at(3, 3, 3).unwrap();
    assert_eq!(force, acceleration.map(|a| after.density * a));
    let f = flow.grid.populations(slot, lane);
    let uu: f64 = before.velocity.iter().map(|u| u * u).sum();
    let uf: f64 = before
        .velocity
        .iter()
        .zip(force)
        .map(|(u, force)| u * force)
        .sum();
    for q in 1..Q3 {
        let e = [E3[q].0, E3[q].1, E3[q].2].map(f64::from);
        let eu: f64 = e.iter().zip(before.velocity).map(|(e, u)| e * u).sum();
        let ef: f64 = e.iter().zip(force).map(|(e, force)| e * force).sum();
        // Independent even moment, using the accepted endpoint force, not
        // the force used for collision at the old density.
        let physical_pair = f[q] + post[OPP3[q]] + W3[q] * (-3.0 * uf + 9.0 * eu * ef);
        let prescribed_pair = 2.0 * W3[q] * boundary_density(q) * (1.0 - 1.5 * uu + 4.5 * eu * eu);
        assert!(
            (physical_pair - prescribed_pair).abs() < 2e-14,
            "direction {q}"
        );
    }
}

#[test]
fn g1_body_history_survives_strain_inflow_and_cross_tile_piston_remap() {
    let gate = CancelGate::new_clock_free();
    let acceleration = [0.002, -0.003, 0.001];
    let mut flow = surface(&[cell([3, 3, 3], 0.9, [0.01, -0.01, 0.005])], [0.0; 3]);
    flow.rebind_body_acceleration(acceleration, &gate).unwrap();
    let density = flow.sample_at(3, 3, 3).unwrap().density;
    let force = flow.force_density_at(3, 3, 3).unwrap();
    flow.initialize_strain_at([3, 3, 3], [0.001, -0.001, 0.0, 0.002, 0.0, 0.0], 0.8)
        .unwrap();
    let (slot, lane) = flow.grid.resolve_source(3, 3, 3).unwrap();
    let location = MassLocation3::Cell(flow.grid.active_keys()[slot], lane);
    flow.transport_set_interface_mass(location, 0.6 * density, &gate)
        .unwrap();
    let velocity = [-0.005, 0.015, 0.01];
    flow.transport_set_interface_velocity(location, velocity, &gate)
        .unwrap();
    assert_eq!(
        flow.sample_at(3, 3, 3).unwrap().density.to_bits(),
        density.to_bits()
    );
    assert_eq!(flow.force_density_at(3, 3, 3).unwrap(), force);
    let populations = flow.grid.populations(slot, lane);
    let mass = flow.ledger_mass();
    flow.remap_planar_piston_one_cell(PlanarWall3::XMin, true, 1e-12, &gate)
        .unwrap();
    let (slot, lane) = flow.grid.resolve_source(4, 3, 3).unwrap();
    assert_eq!(
        flow.grid.populations(slot, lane).map(f64::to_bits),
        populations.map(f64::to_bits)
    );
    assert_eq!(flow.body_acceleration(), Some(acceleration));
    assert_eq!(flow.force_density_at(4, 3, 3).unwrap(), force);
    assert_eq!(
        flow.sample_at(4, 3, 3).unwrap().density.to_bits(),
        density.to_bits()
    );
    assert_eq!(flow.ledger_mass().to_bits(), mass.to_bits());
    assert!(flow.sample_at(3, 3, 3).is_none());
    for (actual, expected) in flow
        .sample_at(4, 3, 3)
        .unwrap()
        .velocity
        .into_iter()
        .zip(velocity)
    {
        assert!((actual - expected).abs() < 2e-14);
    }
}
