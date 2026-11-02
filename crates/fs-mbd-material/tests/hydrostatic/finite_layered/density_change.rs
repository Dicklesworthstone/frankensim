//! Closed-layer source updates consumed by the native interfaces and mechanics.
use super::*;
use fs_mbd_material::hydrostatic::FiniteTankBody;

fn batch_pose(x: f64, z: f64) -> Se3 {
    Se3::from_parts(So3::identity(), Vec3::new(x, 0.0, z)).unwrap()
}

fn batch_stack<'a>(
    sources: [&'a dyn FluidTransport; 2],
    volumes: [f64; 2],
    floor: f64,
    height: f64,
) -> FiniteLayeredHydrostaticTank<'a> {
    let mut tank_geometry = geometry(height);
    tank_geometry.world_from_tank = pose(floor);
    FiniteLayeredHydrostaticTank::new(
        &[0, 1].map(|i| FiniteHydrostaticLayer {
            source: sources[i],
            liquid_volume_m3: volumes[i],
        }),
        gravity(),
        tank_geometry,
        controls(),
        || false,
    )
    .unwrap()
}

#[test]
fn source_updates_keep_each_mass_and_change_actual_interfaces_energy_and_motion() {
    let bottom = fluid(1000.0);
    let top = fluid(800.0);
    let changed_bottom = fluid(800.0);
    let changed_top = fluid(640.0);
    let tank = tank(&bottom, &top, 8.0);
    let mesh = cube();
    let held = pose(1.5);
    let (next, audit) = tank
        .rebind_sources_at_pose(&[&changed_bottom, &changed_top], &mesh, held, || false)
        .unwrap();
    assert_eq!(
        next.layer_masses_kg().collect::<Vec<_>>(),
        vec![12000.0, 6400.0]
    );
    assert_eq!(
        next.liquid_mass_kg().to_bits(),
        tank.liquid_mass_kg().to_bits()
    );
    assert_eq!(audit.layers[0].liquid_volume_m3, [12.0, 15.0]);
    assert_eq!(audit.layers[1].liquid_volume_m3, [8.0, 10.0]);
    close(audit.interfaces_before[0].surface_height_m, 5.0 / 6.0);
    close(audit.interfaces_before[1].surface_height_m, 1.5);
    close(audit.interfaces_after[0].surface_height_m, 13.0 / 12.0);
    close(audit.interfaces_after[1].surface_height_m, 23.0 / 12.0);

    // Independent first moments of liquid slabs excluding the upright prism:
    // lower G/g=rho0*(6*h0²+1/2), upper G/g=6*rho1*(h1²-h0²).
    let energy = |rho: [f64; 2], h: [f64; 2]| {
        10.0 * (rho[0] * (6.0 * h[0].powi(2) + 0.5) + 6.0 * rho[1] * (h[1].powi(2) - h[0].powi(2)))
    };
    let expected = energy([800.0, 640.0], [13.0 / 12.0, 23.0 / 12.0])
        - energy([1000.0, 800.0], [5.0 / 6.0, 1.5]);
    close(audit.empty_tank_gravitational_energy_change_j, 27250.0);
    close(audit.displacement_potential_change_j, 7750.0);
    close(audit.liquid_gravitational_energy_change_j, expected);
    close(expected, 35000.0);

    let source = body(600.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let state = body.state_from_kinematics(held, Twist::zero()).unwrap();
    let before = tank.load(&mesh, &body, state, || false).unwrap();
    let after = next.load(&mesh, &body, state, || false).unwrap();
    close(before.displaced_mass_kg(), 10400.0 / 3.0);
    close(after.displaced_mass_kg(), 4000.0);
    close(
        after.pressure_potential_j() - before.pressure_potential_j(),
        audit.displacement_potential_change_j,
    );
    let step = next
        .step(
            &mesh,
            &body,
            state,
            HydrostaticStepControls {
                duration_s: 0.001,
                maximum_energy_defect_j: 1e-6,
            },
            || false,
        )
        .unwrap();
    // New buoyancy is 40000 N; body mass 4800 kg gives a=-5/3 m/s².
    close(
        body.world_from_body(step.dynamics.state_after)
            .unwrap()
            .translation()
            .z,
        1.5 - (5.0 / 6.0) * 0.001_f64.powi(2),
    );

    let (restored, reverse) = next
        .rebind_sources_at_pose(&[&bottom, &top], &mesh, held, || false)
        .unwrap();
    assert_eq!(
        restored.layer_masses_kg().collect::<Vec<_>>(),
        tank.layer_masses_kg().collect::<Vec<_>>()
    );
    assert_eq!(reverse.interfaces_after, audit.interfaces_before);
    close(reverse.liquid_gravitational_energy_change_j, -expected);
}

#[test]
fn expanding_a_dry_lower_liquid_accounts_for_lifting_unchanged_upper_matter() {
    let bottom = fluid(1000.0);
    let top = fluid(800.0);
    let changed_bottom = fluid(800.0);
    let mesh = cube();
    for floor in [0.0, 10000.0] {
        let mut tank_geometry = geometry(6.0);
        tank_geometry.world_from_tank = pose(floor);
        let tank = FiniteLayeredHydrostaticTank::new(
            &[
                FiniteHydrostaticLayer {
                    source: &bottom,
                    liquid_volume_m3: 12.0,
                },
                FiniteHydrostaticLayer {
                    source: &top,
                    liquid_volume_m3: 8.0,
                },
            ],
            gravity(),
            tank_geometry,
            controls(),
            || false,
        )
        .unwrap();
        let (_, audit) = tank
            .rebind_sources_at_pose(&[&changed_bottom, &top], &mesh, pose(floor + 4.0), || false)
            .unwrap();
        assert!(
            audit
                .interfaces_after
                .iter()
                .all(|fill| fill.displaced_below.volume_m3() == 0.0)
        );
        assert_eq!(audit.layers[1].liquid_volume_m3, [8.0, 8.0]);
        assert_eq!(audit.displacement_potential_change_j, 0.0);
        // Lower liquid gains 11250 J; its 3/16 m rise lifts the unchanged
        // 6400 kg upper liquid by another 12000 J at g=10 m/s².
        close(
            audit.empty_tank_gravitational_energy_change_j,
            11250.0 + 12000.0,
        );
        close(audit.liquid_gravitational_energy_change_j, 23250.0);
        close(audit.interfaces_after[1].surface_height_m, 23.0 / 16.0);
    }
}

#[test]
fn repeated_updates_copy_original_masses_and_equal_density_refresh_freezes_once() {
    struct Source {
        rho: f64,
        mu: f64,
        calls: std::cell::Cell<usize>,
    }
    impl FluidTransport for Source {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            self.calls.set(self.calls.get() + 1);
            FluidTransportState::try_new(self.rho, 4000.0, 0.5, self.mu)
        }
    }
    let source = |rho, mu| Source {
        rho,
        mu,
        calls: std::cell::Cell::new(0),
    };
    let bottom = source(997.3, 0.001);
    let top = source(801.7, 0.001);
    let changed_bottom = source(654.321, 0.001);
    let changed_top = source(400.1, 0.001);
    let same_bottom = source(654.321, 0.002);
    let same_top = source(400.1, 0.002);
    let tank = FiniteLayeredHydrostaticTank::new(
        &[
            FiniteHydrostaticLayer {
                source: &bottom,
                liquid_volume_m3: 12.123,
            },
            FiniteHydrostaticLayer {
                source: &top,
                liquid_volume_m3: 8.456,
            },
        ],
        gravity(),
        geometry(6.0),
        controls(),
        || false,
    )
    .unwrap();
    let masses: Vec<_> = tank.layer_masses_kg().map(f64::to_bits).collect();
    let mesh = cube();
    let held = pose(4.0);
    let (mut current, _) = tank
        .rebind_sources_at_pose(&[&changed_bottom, &changed_top], &mesh, held, || false)
        .unwrap();
    assert_eq!((bottom.calls.get(), top.calls.get()), (1, 1));
    assert_eq!(
        (changed_bottom.calls.get(), changed_top.calls.get()),
        (1, 1)
    );
    for ((input, state), mass) in current.layers().zip(current.layer_masses_kg()) {
        // These source values deliberately expose a rho_new*V_new round trip.
        assert_ne!(
            (state.density_kg_m3() * input.liquid_volume_m3).to_bits(),
            mass.to_bits()
        );
    }
    let (same, audit) = current
        .rebind_sources_at_pose(&[&same_bottom, &same_top], &mesh, held, || false)
        .unwrap();
    assert_eq!((same_bottom.calls.get(), same_top.calls.get()), (1, 1));
    assert_eq!(
        same.layer_masses_kg().map(f64::to_bits).collect::<Vec<_>>(),
        masses
    );
    assert_eq!(audit.interfaces_before, audit.interfaces_after);
    assert_eq!(
        audit.liquid_gravitational_energy_change_j.to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!(
        audit.empty_tank_gravitational_energy_change_j.to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!(
        audit.displacement_potential_change_j.to_bits(),
        0.0_f64.to_bits()
    );
    for ((a, _), (b, state)) in current.layers().zip(same.layers()) {
        assert_eq!(a.liquid_volume_m3.to_bits(), b.liquid_volume_m3.to_bits());
        assert_eq!(state.dynamic_viscosity_pa_s(), 0.002);
    }
    same.interfaces_at_pose(&mesh, held, || false).unwrap();
    assert_eq!((same_bottom.calls.get(), same_top.calls.get()), (1, 1));
    for _ in 0..4 {
        let (back, _) = current
            .rebind_sources_at_pose(&[&bottom, &top], &mesh, held, || false)
            .unwrap();
        let (forward, _) = back
            .rebind_sources_at_pose(&[&changed_bottom, &changed_top], &mesh, held, || false)
            .unwrap();
        current = forward;
        assert_eq!(
            current
                .layer_masses_kg()
                .map(f64::to_bits)
                .collect::<Vec<_>>(),
            masses
        );
        assert_eq!(
            current.liquid_mass_kg().to_bits(),
            tank.liquid_mass_kg().to_bits()
        );
    }
}

#[test]
fn a_single_layer_update_matches_the_original_uniform_inventory_owner() {
    let original = fluid(1000.0);
    let changed = fluid(800.0);
    let layered = FiniteLayeredHydrostaticTank::new(
        &[FiniteHydrostaticLayer {
            source: &original,
            liquid_volume_m3: 12.123,
        }],
        gravity(),
        geometry(6.0),
        controls(),
        || false,
    )
    .unwrap();
    let uniform =
        FiniteHydrostaticTank::new(&original, gravity(), geometry(6.0), 12.123, controls())
            .unwrap();
    let mesh = cube();
    let held = pose(1.5);
    let (new_layered, a) = layered
        .rebind_sources_at_pose(&[&changed], &mesh, held, || false)
        .unwrap();
    let (new_uniform, b) = uniform
        .rebind_source_at_pose(&changed, &mesh, held, || false)
        .unwrap();
    assert_eq!(
        a.empty_tank_gravitational_energy_change_j,
        b.empty_tank_gravitational_energy_change_j
    );
    assert_eq!(
        a.displacement_potential_change_j,
        b.displacement_potential_change_j
    );
    assert_eq!(
        a.liquid_gravitational_energy_change_j,
        b.liquid_gravitational_energy_change_j
    );
    assert_eq!(
        a.interfaces_after[0].surface_height_m,
        b.fill_after.surface_height_m
    );
    let source = body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let state = body.state_from_kinematics(held, Twist::zero()).unwrap();
    assert_eq!(
        new_layered.load(&mesh, &body, state, || false).unwrap(),
        new_uniform.load(&mesh, &body, state, || false).unwrap()
    );
    let step = HydrostaticStepControls {
        duration_s: 0.001,
        maximum_energy_defect_j: 0.1,
    };
    assert_eq!(
        new_layered
            .step(&mesh, &body, state, step, || false)
            .unwrap(),
        new_uniform
            .step(&mesh, &body, state, step, || false)
            .unwrap()
    );
}

#[test]
fn failed_sources_ordering_overflow_and_late_cancellation_leave_all_original_liquids_intact() {
    let bottom = fluid(1000.0);
    let top = fluid(800.0);
    let changed_bottom = fluid(800.0);
    let changed_top = fluid(640.0);
    let tank = FiniteLayeredHydrostaticTank::new(
        &[
            FiniteHydrostaticLayer {
                source: &bottom,
                liquid_volume_m3: 12.0,
            },
            FiniteHydrostaticLayer {
                source: &top,
                liquid_volume_m3: 8.0,
            },
        ],
        gravity(),
        geometry(2.2),
        controls(),
        || false,
    )
    .unwrap();
    let mesh = cube();
    let held = pose(1.5);
    let original = tank.interfaces_at_pose(&mesh, held, || false).unwrap();
    assert!(
        tank.rebind_sources_at_pose(&[&changed_bottom], &mesh, held, || false)
            .is_err()
    );
    let unstable = fluid(700.0);
    assert!(
        tank.rebind_sources_at_pose(&[&unstable, &top], &mesh, held, || false)
            .is_err()
    );
    for (rho_bottom, rho_top) in [(700.0, 560.0), (500.0, 400.0)] {
        let a = fluid(rho_bottom);
        let b = fluid(rho_top);
        // First empty stack fits but body displacement overflows; the second
        // empty stack itself exceeds capacity. Neither spills or removes mass.
        assert!(matches!(
            tank.rebind_sources_at_pose(&[&a, &b], &mesh, held, || false),
            Err(HydrostaticError::TankOverflow { .. })
        ));
    }
    struct Refused;
    impl FluidTransport for Refused {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            FluidTransportState::try_new(-1.0, 4000.0, 0.5, 0.001)
        }
    }
    assert!(matches!(
        tank.rebind_sources_at_pose(&[&Refused, &top], &mesh, held, || false),
        Err(HydrostaticError::Material(_))
    ));
    let mut polls = 0;
    let (_, expected) = tank
        .rebind_sources_at_pose(&[&changed_bottom, &changed_top], &mesh, held, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            tank.rebind_sources_at_pose(&[&changed_bottom, &changed_top], &mesh, held, || {
                count += 1;
                count == stop
            },),
            Err(HydrostaticError::Cancelled)
        ));
    }
    assert_eq!(
        tank.interfaces_at_pose(&mesh, held, || false).unwrap(),
        original
    );
    assert_eq!(
        tank.layer_masses_kg().collect::<Vec<_>>(),
        vec![12000.0, 6400.0]
    );
    assert_eq!(tank.liquid_mass_kg(), 18400.0);
    let (_, replay) = tank
        .rebind_sources_at_pose(&[&changed_bottom, &changed_top], &mesh, held, || false)
        .unwrap();
    assert_eq!(replay, expected);
}

// G1/G3/G4: source-driven density/volume changes reach native coupled motion.
#[test]
fn shared_layer_source_updates_match_independent_stack_integrals_and_body_motion() {
    let bottom = fluid(1000.0);
    let top = fluid(800.0);
    let changed_bottom = fluid(900.0);
    let changed_top = fluid(700.0);
    let tank = batch_stack([&bottom, &top], [16.0; 2], 0.0, 4.0);
    let sides = Vec3::new(1.0, 1.0, 1.0);
    let source = body(900.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let bindings = [FiniteTankBody {
        envelope: &mesh,
        frame: body.frame(),
    }; 2];
    let z = [1.1_f64, 1.05];
    let poses = [batch_pose(-1.0, z[0]), batch_pose(1.0, z[1])];
    let states = poses.map(|pose| body.state_from_kinematics(pose, Twist::zero()).unwrap());
    let (next, audit) = tank
        .rebind_sources_at_bodies(&[&changed_bottom, &changed_top], &bindings, &states, || {
            false
        })
        .unwrap();
    assert_eq!(audit.world_from_body, poses);
    assert_eq!(
        next.layer_masses_kg().collect::<Vec<_>>(),
        [16000.0, 12800.0]
    );
    assert_eq!(
        next.liquid_mass_kg().to_bits(),
        tank.liquid_mass_kg().to_bits()
    );
    let densities = [[1000.0, 800.0], [900.0, 700.0]];
    let heights = densities.map(|rho| {
        let v = [16000.0 / rho[0], 12800.0 / rho[1]];
        [
            (v[0] + 1.0 - z[0] - z[1]) / 14.0,
            (v[0] + v[1] + 2.0) / 16.0,
        ]
    });
    for (interfaces, h) in [&audit.interfaces_before, &audit.interfaces_after]
        .into_iter()
        .zip(heights)
    {
        for i in 0..2 {
            close(interfaces[i].surface_height_m, h[i]);
            assert!(interfaces[i].cumulative_volume_defect_m3.abs() <= 1e-12);
        }
    }
    // Direct first-moment integrals of each liquid slab minus both body parts.
    let energy = |rho: [f64; 2], h: [f64; 2]| {
        let lower = 8.0 * h[0] * h[0]
            - 0.5
                * z.iter()
                    .map(|z| h[0] * h[0] - (z - 0.5).powi(2))
                    .sum::<f64>();
        let upper = 8.0 * (h[1] * h[1] - h[0] * h[0])
            - 0.5
                * z.iter()
                    .map(|z| (z + 0.5).powi(2) - h[0] * h[0])
                    .sum::<f64>();
        10.0 * (rho[0] * lower + rho[1] * upper)
    };
    let expected = energy(densities[1], heights[1]) - energy(densities[0], heights[0]);
    close(audit.liquid_gravitational_energy_change_j, expected);
    let before = tank.load_bodies(&bindings, &states, || false).unwrap();
    let after = next.load_bodies(&bindings, &states, || false).unwrap();
    close(
        audit.displacement_potential_change_j,
        after.surface.liquid_gravitational_energy_change_j
            - before.surface.liquid_gravitational_energy_change_j,
    );
    let dt = 0.001_f64;
    let step = next
        .step_bodies(
            &bindings,
            &states,
            HydrostaticStepControls {
                duration_s: dt,
                maximum_energy_defect_j: 1e-5,
            },
            || false,
        )
        .unwrap();
    for i in 0..2 {
        let displaced_mass = 700.0 + 200.0 * (heights[1][0] - z[i] + 0.5);
        close(after.bodies[i].displaced_mass_kg(), displaced_mass);
        close(
            body.world_from_body(step.dynamics[i].state_after)
                .unwrap()
                .translation()
                .z,
            z[i] + 0.5 * dt * dt * 10.0 * (displaced_mass / 900.0 - 1.0),
        );
    }
    let (restored, reverse) = next
        .rebind_sources_at_bodies(&[&bottom, &top], &bindings, &states, || false)
        .unwrap();
    assert_eq!(
        restored.layer_masses_kg().collect::<Vec<_>>(),
        [16000.0, 12800.0]
    );
    assert_eq!(reverse.interfaces_after, audit.interfaces_before);
    close(reverse.liquid_gravitational_energy_change_j, -expected);

    let (_, one) = tank
        .rebind_sources_at_pose(&[&changed_bottom, &changed_top], &mesh, poses[0], || false)
        .unwrap();
    let (_, batch) = tank
        .rebind_sources_at_bodies(
            &[&changed_bottom, &changed_top],
            &bindings[..1],
            &states[..1],
            || false,
        )
        .unwrap();
    assert_eq!(one.layers, batch.layers);
    assert_eq!(
        one.liquid_gravitational_energy_change_j,
        batch.liquid_gravitational_energy_change_j
    );
    for (a, b) in one.interfaces_after.iter().zip(batch.interfaces_after) {
        assert_eq!(a.surface_height_m, b.surface_height_m);
        assert_eq!(a.displaced_below, b.displaced_below[0]);
        assert_eq!(a.volume_evaluations, b.volume_evaluations);
    }
}

#[test]
fn shared_density_change_preserves_uniform_parity_and_dry_upper_liquid_lift() {
    let bottom = fluid(1000.0);
    let top = fluid(800.0);
    let sides = Vec3::new(1.0, 1.0, 1.0);
    let source = body(900.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let bindings = [FiniteTankBody {
        envelope: &mesh,
        frame: body.frame(),
    }; 2];
    for floor in [0.0, 10000.0] {
        let tank = batch_stack([&bottom, &top], [16.0; 2], floor, 6.0);
        let states = [-1.0, 1.0].map(|x| {
            body.state_from_kinematics(batch_pose(x, floor + 4.0), Twist::zero())
                .unwrap()
        });
        let (_, audit) = tank
            .rebind_sources_at_bodies(&[&top, &top], &bindings, &states, || false)
            .unwrap();
        assert_eq!(audit.layers[1].liquid_volume_m3, [16.0, 16.0]);
        assert!(audit.interfaces_after.iter().all(|interface| {
            interface
                .displaced_below
                .iter()
                .all(|cut| cut.volume_m3() == 0.0)
        }));
        assert_eq!(audit.displacement_potential_change_j, 0.0);
        // 16000 kg lower liquid rises 1/8 m; 12800 kg upper liquid rises 1/4 m.
        close(
            audit.empty_tank_gravitational_energy_change_j,
            20000.0 + 32000.0,
        );
        close(audit.liquid_gravitational_energy_change_j, 52000.0);
    }
    let layered = batch_stack([&bottom, &bottom], [16.0; 2], 0.0, 6.0);
    let uniform =
        FiniteHydrostaticTank::new(&bottom, gravity(), geometry(6.0), 32.0, controls()).unwrap();
    let states = [-1.0, 1.0].map(|x| {
        body.state_from_kinematics(batch_pose(x, 2.1), Twist::zero())
            .unwrap()
    });
    let (new_layered, a) = layered
        .rebind_sources_at_bodies(&[&top, &top], &bindings, &states, || false)
        .unwrap();
    let (new_uniform, b) = uniform
        .rebind_source_at_bodies(&top, &bindings, &states, || false)
        .unwrap();
    assert_eq!(
        a.empty_tank_gravitational_energy_change_j,
        b.empty_tank_gravitational_energy_change_j
    );
    assert_eq!(
        a.displacement_potential_change_j,
        b.displacement_potential_change_j
    );
    assert_eq!(
        a.liquid_gravitational_energy_change_j,
        b.liquid_gravitational_energy_change_j
    );
    assert_eq!(
        new_layered
            .load_bodies(&bindings, &states, || false)
            .unwrap(),
        new_uniform
            .load_bodies(&bindings, &states, || false)
            .unwrap()
    );
    let (_, same) = new_layered
        .rebind_sources_at_bodies(&[&top, &top], &bindings, &states, || false)
        .unwrap();
    assert_eq!(same.interfaces_before, same.interfaces_after);
    assert_eq!(
        same.liquid_gravitational_energy_change_j.to_bits(),
        0.0_f64.to_bits()
    );
}

#[test]
fn shared_layers_retain_exact_mass_and_frozen_sources_through_repeats_and_atomic_refusal() {
    struct Source {
        rho: f64,
        calls: std::cell::Cell<usize>,
    }
    impl FluidTransport for Source {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            self.calls.set(self.calls.get() + 1);
            FluidTransportState::try_new(self.rho, 4000.0, 0.5, 0.001)
        }
    }
    let source = |rho| Source {
        rho,
        calls: std::cell::Cell::new(0),
    };
    let bottom = source(997.3);
    let top = source(801.7);
    let changed_bottom = source(654.321);
    let changed_top = source(400.1);
    let refused = source(-1.0);
    let tank = batch_stack([&bottom, &top], [16.123, 16.456], 0.0, 6.0);
    let sides = Vec3::new(1.0, 1.0, 1.0);
    let material = body(900.0, sides);
    let body = material.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let bindings = [FiniteTankBody {
        envelope: &mesh,
        frame: body.frame(),
    }; 2];
    let states = [-1.0, 1.0].map(|x| {
        body.state_from_kinematics(batch_pose(x, 2.1), Twist::zero())
            .unwrap()
    });
    let originals = states;
    let masses = tank.layer_masses_kg().map(f64::to_bits).collect::<Vec<_>>();
    let before = tank.load_bodies(&bindings, &states, || false).unwrap();
    let (mut next, _) = tank
        .rebind_sources_at_bodies(&[&changed_bottom, &changed_top], &bindings, &states, || {
            false
        })
        .unwrap();
    assert_eq!(
        (
            bottom.calls.get(),
            top.calls.get(),
            changed_bottom.calls.get(),
            changed_top.calls.get()
        ),
        (1, 1, 1, 1)
    );
    next.load_bodies(&bindings, &states, || false).unwrap();
    next.interfaces_at_bodies(&bindings, &states, || false)
        .unwrap();
    assert_eq!(
        (changed_bottom.calls.get(), changed_top.calls.get()),
        (1, 1)
    );
    for _ in 0..4 {
        let (back, _) = next
            .rebind_sources_at_bodies(&[&bottom, &top], &bindings, &states, || false)
            .unwrap();
        next = back
            .rebind_sources_at_bodies(&[&changed_bottom, &changed_top], &bindings, &states, || {
                false
            })
            .unwrap()
            .0;
        assert_eq!(
            next.layer_masses_kg().map(f64::to_bits).collect::<Vec<_>>(),
            masses
        );
        assert_eq!(
            next.liquid_mass_kg().to_bits(),
            tank.liquid_mass_kg().to_bits()
        );
    }
    assert!(
        tank.rebind_sources_at_bodies(&[&changed_bottom], &bindings, &states, || false)
            .is_err()
    );
    assert!(
        tank.rebind_sources_at_bodies(&[&changed_top, &changed_bottom], &bindings, &states, || {
            false
        })
        .is_err()
    );
    assert!(matches!(
        tank.rebind_sources_at_bodies(&[&refused, &top], &bindings, &states, || false),
        Err(HydrostaticError::Material(_))
    ));
    assert!(
        tank.rebind_sources_at_bodies(&[&bottom, &top], &bindings, &states[..1], || false)
            .is_err()
    );
    let mut polls = 0;
    let (_, expected) = tank
        .rebind_sources_at_bodies(&[&changed_bottom, &changed_top], &bindings, &states, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut count = 0;
        assert!(matches!(
            tank.rebind_sources_at_bodies(
                &[&changed_bottom, &changed_top],
                &bindings,
                &states,
                || {
                    count += 1;
                    count == stop
                }
            ),
            Err(HydrostaticError::Cancelled)
        ));
    }
    assert_eq!(states, originals);
    assert_eq!(
        tank.load_bodies(&bindings, &states, || false).unwrap(),
        before
    );
    assert_eq!(
        tank.rebind_sources_at_bodies(&[&changed_bottom, &changed_top], &bindings, &states, || {
            false
        })
        .unwrap()
        .1,
        expected
    );

    let low = fluid(1000.0);
    let high = fluid(800.0);
    let expanded_low = fluid(14000.0 / 15.5);
    let expanded_high = fluid(11200.0 / 15.0);
    let shallow = batch_stack([&low, &high], [14.0; 2], 0.0, 2.0);
    let poses = [batch_pose(-1.0, 1.0), batch_pose(1.0, 1.0)];
    let submerged = poses.map(|pose| body.state_from_kinematics(pose, Twist::zero()).unwrap());
    for pose in poses {
        shallow
            .rebind_sources_at_pose(&[&expanded_low, &expanded_high], &mesh, pose, || false)
            .unwrap();
    }
    assert!(matches!(
        shallow.rebind_sources_at_bodies(
            &[&expanded_low, &expanded_high],
            &bindings,
            &submerged,
            || false
        ),
        Err(HydrostaticError::TankOverflow { .. })
    ));
    assert_eq!(
        shallow.layer_masses_kg().collect::<Vec<_>>(),
        [14000.0, 11200.0]
    );
}
