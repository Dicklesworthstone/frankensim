//! Closed-layer source updates consumed by the native interfaces and mechanics.
use super::*;

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
