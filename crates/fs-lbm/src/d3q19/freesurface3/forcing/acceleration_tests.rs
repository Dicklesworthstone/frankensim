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
    assert_eq!(before.material_gravity, after.material_gravity);
    same_physical_state(before, after);
}

fn same_physical_state(before: &FreeSurface3, after: &FreeSurface3) {
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

fn flooded_tank() -> FreeSurface3 {
    let cells: Vec<_> = (0..8)
        .flat_map(|z| {
            (0..8).flat_map(move |y| {
                (0..8).map(move |x| FreeSurfaceInitialCell3 {
                    coordinates: [x, y, z],
                    fill_fraction: 1.0,
                    density: 1.0
                        + 0.02 * (x as f64 * core::f64::consts::TAU / 8.0).cos()
                        + 0.01 * (y as f64 * core::f64::consts::TAU / 8.0).sin(),
                    velocity: [0.01, -0.008, 0.002],
                })
            })
        })
        .collect();
    surface(&cells, [0.0; 3])
}

#[test]
fn g1_rigid_wall_acceleration_keeps_actual_density_impulse_and_wall_work_separate() {
    let gate = CancelGate::new_clock_free();
    let old = affine_field();
    let next = BodyAccelerationField3::Affine {
        origin: [4.0; 3],
        acceleration: [-0.002, 0.001, -0.001],
        gradient: [[0.0, 0.0003, 0.0], [-0.0004, 0.0, 0.0], [0.0; 3]],
    };
    for wall in [PlanarWall3::ZMin, PlanarWall3::XMax] {
        let mut initial = flooded_tank();
        initial.rebind_body_acceleration_field(old, &gate).unwrap();
        let mut translation = [0.004, -0.003, 0.002];
        translation[wall.axis()] = 0.0;
        let motion = PlanarRigidWallMotion3 {
            wall,
            origin: [4.0; 3],
            velocity_at_origin: translation,
            angular_speed: 0.001,
        };
        let (held, held_wall) = initial
            .planar_rigid_wall_motion_candidate(motion, 1e-12, &gate)
            .unwrap();
        let (evolved, wall_report) = initial
            .planar_rigid_wall_motion_with_body_acceleration_candidate(
                motion, 1e-12, next, &gate,
            )
            .unwrap();
        assert_eq!(initial.body_acceleration_field(), Some(old));
        assert_eq!(evolved.body_acceleration_field(), Some(next));
        // Without gas links, both candidates collide/stream precisely the same
        // OLD source and reflection. Only the physical endpoint half source
        // differs. A zero-time reframe, or collision with the future source,
        // would violate these population and physical-momentum identities.
        assert_eq!(held.grid.state_bits(), evolved.grid.state_bits());
        assert_eq!(held_wall.fluid_impulse, wall_report.fluid_impulse);
        assert_eq!(
            held_wall.fluid_angular_impulse,
            wall_report.fluid_angular_impulse
        );
        assert_eq!(held_wall.fluid_response, wall_report.fluid_response);
        assert_eq!(held_wall.fluid_work, wall_report.fluid_work);
        assert_eq!(wall_report.wall_work, -wall_report.fluid_work);
        let work = (0..3)
            .map(|d| translation[d] * wall_report.fluid_impulse[d])
            .sum::<f64>()
            + motion.angular_speed * wall_report.fluid_angular_impulse[wall.axis()];
        assert!((wall_report.fluid_work - work).abs() < 2e-13);
        assert!((evolved.ledger_mass() - initial.ledger_mass()).abs() < 2e-11);
        let mut density_change = 0.0_f64;
        let mut impulse = [0.0; 3];
        let mut expected_impulse = [0.0; 3];
        for (&key, tile) in &evolved.tiles {
            for (lane, &cell) in tile.cells.iter().enumerate() {
                if cell == Cell3::Gas {
                    continue;
                }
                let (x, y, z) = FreeSurface3::global_of(key, lane);
                let xyz = [x as usize, y as usize, z as usize];
                let before = initial.sample_at(xyz[0], xyz[1], xyz[2]).unwrap();
                let after = evolved.sample_at(xyz[0], xyz[1], xyz[2]).unwrap();
                let p0 = moments(&held, xyz).0;
                let p1 = moments(&evolved, xyz).0;
                let a0 = old.at_cell(key, lane).unwrap();
                let a1 = next.at_cell(key, lane).unwrap();
                assert_eq!(
                    evolved.force_density_at(xyz[0], xyz[1], xyz[2]).unwrap(),
                    a1.map(|a| after.density * a)
                );
                density_change = density_change.max((after.density - before.density).abs());
                for d in 0..3 {
                    let expected = 0.5 * after.density * (a1[d] - a0[d]);
                    assert!((p1[d] - p0[d] - expected).abs() < 3e-15);
                    impulse[d] += p1[d] - p0[d];
                    expected_impulse[d] += expected;
                }
            }
        }
        assert!(density_change > 1e-4, "exercise transported density");
        for d in 0..3 {
            assert!((impulse[d] - expected_impulse[d]).abs() < 3e-13);
        }
        let mut committed = initial.clone();
        let receipt = committed
            .advance_body_acceleration_with_planar_rigid_wall_motion(
                next, motion, 1e-12, &gate,
            )
            .unwrap();
        unchanged(&evolved, &committed);
        assert_eq!(wall_report, receipt);
        let (equal, same_receipt) = initial
            .planar_rigid_wall_motion_with_body_acceleration_candidate(
                motion, 1e-12, old, &gate,
            )
            .unwrap();
        unchanged(&held, &equal);
        assert_eq!(held_wall, same_receipt);
    }
}

#[test]
fn g1_rigid_wall_acceleration_uses_next_field_for_converted_interfaces() {
    let gate = CancelGate::new_clock_free();
    let mut flow = FreeSurface3::new(
        8,
        8,
        8,
        0.8,
        [0.0; 3],
        0.0,
        ContactModel3::Neutral,
        |_, _, z| z < 4,
    )
    .unwrap();
    flow.rebind_body_acceleration([0.0; 3], &gate).unwrap();
    // One accepted interface has conversion excess. Its redistribution runs
    // during the same step as wall reflection and the new force endpoint.
    let (slot, lane) = flow.grid.resolve_source(3, 3, 3).unwrap();
    let key = flow.grid.active_keys()[slot];
    assert_eq!(flow.tiles[&key].cells[lane], Cell3::Interface);
    let rho = flow.grid.populations(slot, lane).iter().sum::<f64>();
    flow.tiles.get_mut(&key).unwrap().mass[lane] = 1.01 * rho;
    let before = flow.clone();
    let motion = PlanarRigidWallMotion3 {
        wall: PlanarWall3::ZMin,
        origin: [4.0, 4.0, 0.0],
        velocity_at_origin: [0.003, -0.002, 0.0],
        angular_speed: 0.0005,
    };
    let next = affine_field();
    flow.advance_body_acceleration_with_planar_rigid_wall_motion(next, motion, 1e-12, &gate)
        .unwrap();
    assert!(flow.conversions.gas_to_interface > before.conversions.gas_to_interface);
    assert!((flow.ledger_mass() - before.ledger_mass()).abs() < 2e-11);
    for (&key, tile) in &flow.tiles {
        for (lane, &cell) in tile.cells.iter().enumerate() {
            if cell != Cell3::Gas {
                let slot = flow.grid.slot_of(key).unwrap();
                let rho = flow.grid.populations(slot, lane).iter().sum::<f64>();
                assert_eq!(
                    flow.force_at_lane(key, lane),
                    next.at_cell(key, lane).unwrap().map(|a| rho * a)
                );
            }
        }
    }
}

#[test]
fn g0_acceleration_endpoint_mach_and_rigid_wall_refusals_keep_the_accepted_state() {
    let gate = CancelGate::new_clock_free();
    let mut flow = flooded_tank();
    flow.rebind_body_acceleration([0.0; 3], &gate).unwrap();
    let before = flow.clone();
    let excessive = BodyAccelerationField3::Uniform([0.5, 0.0, 0.0]);
    // Collision still sees zero force and produces valid populations. Only
    // actual endpoint velocity reveals the high-Mach half-source increment.
    assert!(matches!(
        flow.advance_body_acceleration_field(excessive, &gate),
        Err(FreeSurfaceError3::Transport(
            "body-acceleration endpoint exceeds the physical lattice Mach domain"
        ))
    ));
    unchanged(&before, &flow);
    let motion = PlanarRigidWallMotion3 {
        wall: PlanarWall3::ZMin,
        origin: [4.0, 4.0, 0.0],
        velocity_at_origin: [0.003, -0.002, 0.0],
        angular_speed: 0.0005,
    };
    for next in [excessive, BodyAccelerationField3::Uniform([f64::NAN; 3])] {
        assert!(
            flow.advance_body_acceleration_with_planar_rigid_wall_motion(
                next, motion, 1e-12, &gate,
            )
            .is_err()
        );
        unchanged(&before, &flow);
    }
    let next = affine_field();
    gate.request();
    assert!(matches!(
        flow.advance_body_acceleration_with_planar_rigid_wall_motion(next, motion, 1e-12, &gate),
        Err(FreeSurfaceError3::Cancelled)
    ));
    unchanged(&before, &flow);
    let retry = CancelGate::new_clock_free();
    let (expected, expected_receipt) = before
        .planar_rigid_wall_motion_with_body_acceleration_candidate(motion, 1e-12, next, &retry)
        .unwrap();
    let receipt = flow
        .advance_body_acceleration_with_planar_rigid_wall_motion(next, motion, 1e-12, &retry)
        .unwrap();
    unchanged(&expected, &flow);
    assert_eq!(expected_receipt, receipt);
}

#[test]
fn g3_material_force_handoff_restores_actual_density_gravity_without_impulse() {
    let cells = [
        cell([2, 2, 2], 0.6, [0.02, -0.01, 0.005]),
        cell([5, 5, 5], 1.4, [-0.01, 0.02, -0.005]),
    ];
    let gate = CancelGate::new_clock_free();
    let gravity = [0.002, -0.003, 0.001];
    let mut flow = surface(&cells, [0.0; 3]);
    flow.rebind_material_gravity(gravity, &gate).unwrap();
    let before = flow.clone();
    let forces = flow.material_gravity_forces(&gate).unwrap();
    flow.rebind_material_force_field(forces, &gate).unwrap();
    unchanged(&before, &flow);

    let extra = [0.0003, 0.0004, -0.0002];
    let mut forces = flow.material_gravity_forces(&gate).unwrap();
    for force in forces.values_mut() {
        for d in 0..3 {
            force[d] += extra[d];
        }
    }
    flow.rebind_material_force_field(forces, &gate).unwrap();
    assert_eq!(flow.body_acceleration(), None);
    assert_eq!(flow.material_gravity(), Some(gravity));
    for c in cells {
        let sample = flow
            .sample_at(c.coordinates[0], c.coordinates[1], c.coordinates[2])
            .unwrap();
        let actual = flow
            .force_density_at(c.coordinates[0], c.coordinates[1], c.coordinates[2])
            .unwrap();
        for d in 0..3 {
            assert_eq!(actual[d], sample.density * gravity[d] + extra[d]);
        }
    }
    let forces = flow.material_gravity_forces(&gate).unwrap();
    flow.rebind_material_force_field(forces, &gate).unwrap();
    assert_eq!(flow.body_acceleration(), Some(gravity));
    assert_eq!(flow.ledger_mass(), before.ledger_mass());
    for c in cells {
        let (old_p, old_s) = moments(&before, c.coordinates);
        let (new_p, new_s) = moments(&flow, c.coordinates);
        for (old, new) in old_p
            .into_iter()
            .chain(old_s.into_iter().flatten())
            .zip(new_p.into_iter().chain(new_s.into_iter().flatten()))
        {
            assert!((old - new).abs() < 3e-15);
        }
    }

    let before = flow.clone();
    let mut incomplete = flow.material_gravity_forces(&gate).unwrap();
    incomplete.pop_first();
    assert!(flow.rebind_material_force_field(incomplete, &gate).is_err());
    unchanged(&before, &flow);
    let forces = flow.material_gravity_forces(&gate).unwrap();
    gate.request();
    assert!(matches!(
        flow.rebind_material_force_field(forces, &gate),
        Err(FreeSurfaceError3::Cancelled)
    ));
    unchanged(&before, &flow);
}

#[test]
fn g1_complete_frame_replaces_material_gravity_but_not_independent_acceleration() {
    let cells = [
        cell([2, 2, 2], 0.6, [0.01; 3]),
        cell([5, 5, 5], 1.4, [-0.01; 3]),
    ];
    let gate = CancelGate::new_clock_free();
    let gravity = [0.001, -0.002, 0.003];
    let next = [-0.002, 0.001, -0.001];
    let mut flow = surface(&cells, [0.0; 3]);
    flow.rebind_material_gravity(gravity, &gate).unwrap();
    RigidFrame3::default()
        .prepare_free_surface_step(&mut flow, next, &gate)
        .unwrap();
    assert_eq!(flow.body_acceleration(), None);
    for c in cells {
        let actual = flow
            .force_density_at(c.coordinates[0], c.coordinates[1], c.coordinates[2])
            .unwrap();
        for d in 0..3 {
            assert!((actual[d] - c.density * next[d]).abs() < 2e-17);
        }
    }
    let mut independent = surface(&cells, [0.0; 3]);
    independent
        .rebind_body_acceleration(gravity, &gate)
        .unwrap();
    let before = independent.clone();
    assert!(
        RigidFrame3::default()
            .prepare_free_surface_step(&mut independent, next, &gate)
            .is_err()
    );
    unchanged(&before, &independent);
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
    assert_gas_pressure_endpoint(None);
}

#[test]
fn g1_acceleration_endpoint_and_pressure_close_in_the_same_force_frame() {
    assert_gas_pressure_endpoint(Some([-0.003, 0.001, 0.004]));
}

fn assert_gas_pressure_endpoint(next_acceleration: Option<[f64; 3]>) {
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
    let pressure =
        |link: GasPressureLink3| Ok((boundary_density(link.incoming_direction) - 1.0) * CS2_3);
    if let Some(next) = next_acceleration {
        flow.advance_body_acceleration_with_gas_link_pressure(next, &gate, 0.2, pressure)
            .unwrap();
    } else {
        flow.step_with_gas_link_pressure(&gate, 0.2, pressure)
            .unwrap();
    }
    let after = flow.sample_at(3, 3, 3).unwrap();
    assert!((after.density - before.density).abs() > 0.1);
    let force = flow.force_density_at(3, 3, 3).unwrap();
    let endpoint = next_acceleration.unwrap_or(acceleration);
    assert_eq!(flow.body_acceleration(), Some(endpoint));
    assert_eq!(force, endpoint.map(|a| after.density * a));
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
fn g1_acceleration_endpoints_integrate_a_ramp_at_actual_density() {
    let gate = CancelGate::new_clock_free();
    let velocity = [0.02, -0.01, 0.005];
    let initial_acceleration = [0.001, -0.002, 0.003];
    for density in [0.4, 1.0, 2.5] {
        let mut initial = Vec::new();
        for z in 0..8 {
            for y in 0..8 {
                for x in 0..8 {
                    initial.push(FreeSurfaceInitialCell3 {
                        coordinates: [x, y, z],
                        fill_fraction: 1.0,
                        density,
                        velocity,
                    });
                }
            }
        }
        let mut flow = surface(&initial, [0.0; 3]);
        flow.rebind_body_acceleration(initial_acceleration, &gate)
            .unwrap();
        let mut expected_velocity = velocity;
        let mut previous = initial_acceleration;
        for next in [[-0.003, 0.002, 0.001], [0.004, -0.001, -0.002], [0.0; 3]] {
            let before = flow.sample_at(4, 4, 4).unwrap();
            let (momentum_before, _) = moments(&flow, [4, 4, 4]);
            flow.advance_body_acceleration(next, &gate).unwrap();
            let after = flow.sample_at(4, 4, 4).unwrap();
            let (momentum_after, _) = moments(&flow, [4, 4, 4]);
            assert_eq!(flow.body_acceleration(), Some(next));
            assert_eq!(
                flow.force_density_at(4, 4, 4).unwrap(),
                next.map(|a| after.density * a)
            );
            // The bulk cell has not received any wall disturbance. The
            // analytic time integral of a linear acceleration ramp is exact.
            for axis in 0..3 {
                expected_velocity[axis] += 0.5 * (previous[axis] + next[axis]);
                assert!((after.velocity[axis] - expected_velocity[axis]).abs() < 5e-13);
                let impulse = 0.5 * (before.density * previous[axis] + after.density * next[axis]);
                assert!(
                    (momentum_after[axis] - momentum_before[axis] - impulse).abs()
                        < 5e-13 * density
                );
            }
            assert!((after.density - density).abs() < 5e-13 * density);
            previous = next;
        }
    }
}

#[test]
fn g3_equal_acceleration_endpoints_reproduce_held_force_population_bits() {
    let gate = CancelGate::new_clock_free();
    let acceleration = [0.001, -0.002, 0.0005];
    let mut initial = surface(&[cell([3, 3, 3], 0.9, [0.01, -0.02, 0.005])], [0.0; 3]);
    initial
        .rebind_body_acceleration(acceleration, &gate)
        .unwrap();
    for with_pressure in [false, true] {
        let mut held = initial.clone();
        let mut evolved = initial.clone();
        for _ in 0..3 {
            if with_pressure {
                held.step_with_gas_link_pressure(&gate, 0.1, |_| Ok(0.01))
                    .unwrap();
                evolved
                    .advance_body_acceleration_with_gas_link_pressure(
                        acceleration,
                        &gate,
                        0.1,
                        |_| Ok(0.01),
                    )
                    .unwrap();
            } else {
                held.step_with_gate(&gate).unwrap();
                evolved
                    .advance_body_acceleration(acceleration, &gate)
                    .unwrap();
            }
            unchanged(&held, &evolved);
            assert_eq!(held.worst_step_violation, evolved.worst_step_violation);
        }
    }
}

#[test]
fn g0_acceleration_endpoint_refusal_and_late_cancellation_preserve_history() {
    let gate = CancelGate::new_clock_free();
    let next = [-0.002, 0.001, 0.003];
    let mut flow = surface(&[cell([3, 3, 3], 0.9, [0.01, -0.02, 0.005])], [0.0; 3]);
    let unowned = flow.clone();
    assert!(flow.advance_body_acceleration(next, &gate).is_err());
    unchanged(&unowned, &flow);
    flow.rebind_body_acceleration([0.001, 0.0, -0.001], &gate)
        .unwrap();
    let before = flow.clone();
    for invalid in [[f64::NAN, 0.0, 0.0], [0.0, f64::INFINITY, 0.0]] {
        assert!(flow.advance_body_acceleration(invalid, &gate).is_err());
        unchanged(&before, &flow);
    }
    assert!(
        flow.advance_body_acceleration_with_gas_link_pressure(next, &gate, 0.1, |_| Ok(f64::NAN),)
            .is_err()
    );
    unchanged(&before, &flow);
    let calls = std::cell::Cell::new(0);
    assert!(matches!(
        flow.advance_body_acceleration_with_gas_link_pressure(next, &gate, 0.1, |_| {
            calls.set(calls.get() + 1);
            gate.request();
            Ok(0.01)
        },),
        Err(FreeSurfaceError3::Cancelled)
    ));
    assert!(calls.get() > 0);
    unchanged(&before, &flow);
    let retry = CancelGate::new_clock_free();
    let mut reference = before;
    reference
        .advance_body_acceleration_with_gas_link_pressure(next, &retry, 0.1, |_| Ok(0.01))
        .unwrap();
    flow.advance_body_acceleration_with_gas_link_pressure(next, &retry, 0.1, |_| Ok(0.01))
        .unwrap();
    unchanged(&reference, &flow);
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

fn affine_field() -> BodyAccelerationField3 {
    BodyAccelerationField3::Affine {
        origin: [0.0; 3],
        acceleration: [0.001, -0.002, 0.003],
        gradient: [[0.0005, 0.0, 0.0], [0.0, -0.0002, 0.0], [0.0, 0.0, 0.0001]],
    }
}

#[test]
fn g1_affine_acceleration_integrates_local_momentum_with_evolving_density() {
    let gate = CancelGate::new_clock_free();
    let old = affine_field();
    let next = BodyAccelerationField3::Affine {
        origin: [1.0, 2.0, 3.0],
        acceleration: [-0.002, 0.001, -0.001],
        gradient: [[0.0, 0.0003, 0.0], [-0.0004, 0.0, 0.0], [0.0, 0.0, 0.0002]],
    };
    for density in [0.6, 1.4] {
        let initial: Vec<_> = (0..8)
            .flat_map(|z| {
                (0..8).flat_map(move |y| {
                    (0..8).map(move |x| FreeSurfaceInitialCell3 {
                        coordinates: [x, y, z],
                        fill_fraction: 1.0,
                        density,
                        velocity: [0.0; 3],
                    })
                })
            })
            .collect();
        let mut flow = surface(&initial, [0.0; 3]);
        flow.rebind_body_acceleration_field(old, &gate).unwrap();
        let before = flow.clone();
        flow.advance_body_acceleration_field(next, &gate).unwrap();
        for xyz in [[3, 3, 3], [4, 4, 4]] {
            let [x, y, z] = xyz.map(|coordinate| coordinate as f64 + 0.5);
            let a0 = [0.001 + 0.0005 * x, -0.002 - 0.0002 * y, 0.003 + 0.0001 * z];
            let a1 = [
                -0.002 + 0.0003 * (y - 2.0),
                0.001 - 0.0004 * (x - 1.0),
                -0.001 + 0.0002 * (z - 3.0),
            ];
            let incoming = before.sample_at(xyz[0], xyz[1], xyz[2]).unwrap();
            let accepted = flow.sample_at(xyz[0], xyz[1], xyz[2]).unwrap();
            let (momentum, _) = moments(&flow, xyz);
            // At physical rest, streaming the linear Guo source has zeroth
            // moment -rho*div(a0)/2 and first moment rho*a0(x)/2. The latter
            // then receives the actual endpoint half source, rho_next*a1/2.
            let expected_density = incoming.density * (1.0 - 0.5 * (0.0005 - 0.0002 + 0.0001));
            assert!((accepted.density - expected_density).abs() < 3e-13 * density);
            assert!((accepted.density - incoming.density).abs() > 1e-5 * density);
            assert_eq!(
                flow.force_density_at(xyz[0], xyz[1], xyz[2]).unwrap(),
                a1.map(|a| accepted.density * a)
            );
            for axis in 0..3 {
                let impulse = 0.5 * (incoming.density * a0[axis] + accepted.density * a1[axis]);
                assert!((momentum[axis] - impulse).abs() < 3e-13 * density);
            }
        }
        assert_eq!(flow.body_acceleration_field(), Some(next));
        assert_eq!(flow.body_acceleration(), None);
    }
}

#[test]
fn g3_zero_gradient_field_retains_uniform_population_arithmetic() {
    let gate = CancelGate::new_clock_free();
    let acceleration = [0.001, -0.002, 0.0005];
    let field = |acceleration| BodyAccelerationField3::Affine {
        origin: [11.0, -3.0, 2.0],
        acceleration,
        gradient: [[0.0; 3]; 3],
    };
    let mut uniform = surface(&[cell([3, 3, 3], 0.9, [0.01, -0.02, 0.005])], [0.0; 3]);
    let mut affine = uniform.clone();
    uniform
        .rebind_body_acceleration(acceleration, &gate)
        .unwrap();
    affine
        .rebind_body_acceleration_field(field(acceleration), &gate)
        .unwrap();
    same_physical_state(&uniform, &affine);
    uniform.step_with_gate(&gate).unwrap();
    affine.step_with_gate(&gate).unwrap();
    same_physical_state(&uniform, &affine);
    for next in [[-0.002, 0.001, 0.003], [0.0; 3]] {
        uniform
            .advance_body_acceleration_with_gas_link_pressure(next, &gate, 0.1, |_| Ok(0.01))
            .unwrap();
        affine
            .advance_body_acceleration_field_with_gas_link_pressure(field(next), &gate, 0.1, |_| {
                Ok(0.01)
            })
            .unwrap();
        same_physical_state(&uniform, &affine);
    }
}

#[test]
fn g1_affine_gas_pressure_uses_each_positions_next_density_force() {
    let gate = CancelGate::new_clock_free();
    let coordinates = [[2, 2, 2], [5, 5, 5]];
    let initial = [
        cell(coordinates[0], 0.75, [0.02, -0.01, 0.015]),
        cell(coordinates[1], 1.1, [-0.01, 0.015, 0.02]),
    ];
    let mut flow = surface(&initial, [0.0; 3]);
    flow.rebind_body_acceleration_field(affine_field(), &gate)
        .unwrap();
    let incoming: Vec<_> = coordinates
        .into_iter()
        .map(|xyz| {
            let (slot, lane) = flow
                .grid
                .resolve_source(xyz[0] as i64, xyz[1] as i64, xyz[2] as i64)
                .unwrap();
            let force = flow.force_density_at(xyz[0], xyz[1], xyz[2]).unwrap();
            let post = crate::d3q19::collision_kernel3()(
                flow.grid.populations(slot, lane),
                flow.collision_model(),
                force,
            )
            .unwrap();
            (xyz, flow.sample_at(xyz[0], xyz[1], xyz[2]).unwrap(), post)
        })
        .collect();
    let next = BodyAccelerationField3::Affine {
        origin: [1.0; 3],
        acceleration: [-0.002, 0.001, 0.003],
        gradient: [[0.0002, 0.0, 0.0], [0.0, -0.0003, 0.0], [0.0, 0.0, 0.0004]],
    };
    let boundary_density = |q: usize| 1.05 + 0.01 * f64::from(E3[q].0 - E3[q].1);
    flow.advance_body_acceleration_field_with_gas_link_pressure(next, &gate, 0.2, |link| {
        Ok((boundary_density(link.incoming_direction) - 1.0) * CS2_3)
    })
    .unwrap();
    for (xyz, before, post) in incoming {
        let after = flow.sample_at(xyz[0], xyz[1], xyz[2]).unwrap();
        let r = xyz.map(|coordinate| coordinate as f64 - 0.5);
        let acceleration = [
            -0.002 + 0.0002 * r[0],
            0.001 - 0.0003 * r[1],
            0.003 + 0.0004 * r[2],
        ];
        let force = flow.force_density_at(xyz[0], xyz[1], xyz[2]).unwrap();
        assert_eq!(force, acceleration.map(|a| after.density * a));
        let (slot, lane) = flow
            .grid
            .resolve_source(xyz[0] as i64, xyz[1] as i64, xyz[2] as i64)
            .unwrap();
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
            let ef: f64 = e.iter().zip(force).map(|(e, f)| e * f).sum();
            let physical_pair = f[q] + post[OPP3[q]] + W3[q] * (-3.0 * uf + 9.0 * eu * ef);
            let imposed_pair = 2.0 * W3[q] * boundary_density(q) * (1.0 - 1.5 * uu + 4.5 * eu * eu);
            assert!((physical_pair - imposed_pair).abs() < 3e-14);
        }
    }
}

#[test]
fn g1_affine_field_initializes_new_interfaces_at_their_own_position() {
    let gate = CancelGate::new_clock_free();
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
    flow.rebind_body_acceleration_field(affine_field(), &gate)
        .unwrap();
    let (slot, lane) = flow.grid.resolve_source(4, 3, 3).unwrap();
    let key = flow.grid.active_keys()[slot];
    let density = flow.grid.populations(slot, lane).iter().sum::<f64>();
    flow.tiles.get_mut(&key).unwrap().mass[lane] = 1.1 * density;
    let mass = flow.ledger_mass();
    flow.apply_conversions(None).unwrap();
    let sample = flow.sample_at(5, 3, 3).unwrap();
    let acceleration = [
        0.001 + 0.0005 * 5.5,
        -0.002 - 0.0002 * 3.5,
        0.003 + 0.0001 * 3.5,
    ];
    assert!(flow.conversions.gas_to_interface > 0);
    assert_eq!(sample.cell, Cell3::Interface);
    assert_eq!(
        flow.force_density_at(5, 3, 3).unwrap(),
        acceleration.map(|a| sample.density * a)
    );
    let (_, stress) = moments(&flow, [5, 3, 3]);
    for a in 0..3 {
        assert!((sample.velocity[a] - velocity[a]).abs() < 3e-14);
        for b in 0..3 {
            let expected =
                sample.density * (velocity[a] * velocity[b] + if a == b { CS2_3 } else { 0.0 });
            assert!((stress[a][b] - expected).abs() < 3e-14 * sample.density);
        }
    }
    assert!((flow.ledger_mass() - mass).abs() < 3e-13 * mass);
}

#[test]
fn g1_affine_field_stays_in_world_coordinates_across_piston_remap() {
    let gate = CancelGate::new_clock_free();
    let velocity = [0.02, -0.01, 0.015];
    let mut flow = surface(&[cell([3, 3, 3], 0.9, velocity)], [0.0; 3]);
    let field = affine_field();
    flow.rebind_body_acceleration_field(field, &gate).unwrap();
    let before = flow.sample_at(3, 3, 3).unwrap();
    let (p0, s0) = moments(&flow, [3, 3, 3]);
    let mass = flow.ledger_mass();
    let report = flow
        .remap_planar_piston_one_cell(PlanarWall3::XMin, true, 1e-12, &gate)
        .unwrap();
    let after = flow.sample_at(4, 3, 3).unwrap();
    let (p1, s1) = moments(&flow, [4, 3, 3]);
    let acceleration = [
        0.001 + 0.0005 * 4.5,
        -0.002 - 0.0002 * 3.5,
        0.003 + 0.0001 * 3.5,
    ];
    assert_eq!(flow.body_acceleration_field(), Some(field));
    assert_eq!(
        flow.force_density_at(4, 3, 3).unwrap(),
        acceleration.map(|a| after.density * a)
    );
    assert_eq!(before.density, after.density);
    assert_eq!(flow.ledger_mass(), mass);
    assert_eq!(report.fluid_impulse, [0.0; 3]);
    assert_eq!(report.fluid_work, 0.0);
    assert!(flow.sample_at(3, 3, 3).is_none());
    for a in 0..3 {
        assert!((after.velocity[a] - before.velocity[a]).abs() < 3e-14);
        assert!((p1[a] - p0[a]).abs() < 3e-14);
        for b in 0..3 {
            assert!((s1[a][b] - s0[a][b]).abs() < 3e-14);
        }
    }
}

#[test]
fn g0_affine_force_ownership_and_late_refusal_preserve_the_accepted_field() {
    let gate = CancelGate::new_clock_free();
    let mut flow = surface(&[cell([3, 3, 3], 0.9, [0.01, -0.02, 0.005])], [0.0; 3]);
    flow.rebind_body_acceleration_field(affine_field(), &gate)
        .unwrap();
    assert!(flow.body_acceleration().is_none());
    assert!(flow.body_acceleration_field().is_some());
    let before = flow.clone();
    assert!(
        RigidFrame3::default()
            .prepare_free_surface_step(&mut flow, [0.0; 3], &gate)
            .is_err()
    );
    let independent = flow.local_force.clone().unwrap();
    assert!(flow.rebind_force_field(independent, &gate).is_err());
    unchanged(&before, &flow);
    let invalid = BodyAccelerationField3::Affine {
        origin: [0.0; 3],
        acceleration: [0.0; 3],
        gradient: [[f64::NAN; 3]; 3],
    };
    assert!(
        flow.advance_body_acceleration_field(invalid, &gate)
            .is_err()
    );
    unchanged(&before, &flow);
    let next = BodyAccelerationField3::Affine {
        origin: [1.0; 3],
        acceleration: [0.0; 3],
        gradient: [[0.0001; 3]; 3],
    };
    assert!(matches!(
        flow.advance_body_acceleration_field_with_gas_link_pressure(next, &gate, 0.1, |_| {
            gate.request();
            Ok(0.01)
        },),
        Err(FreeSurfaceError3::Cancelled)
    ));
    unchanged(&before, &flow);
}
