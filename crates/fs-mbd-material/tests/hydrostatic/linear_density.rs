//! G1/G3/G4: continuous density is integrated through force, torque and motion.
use super::*;
use fs_mbd_material::hydrostatic::{HydrostaticStepControls, LinearDensityHydrostaticMedium};

#[test]
fn zero_gradient_is_the_exact_uniform_load_and_step() {
    let water = fluid(1000.0);
    let uniform = HydrostaticMedium::new(&water, gravity(), origin()).unwrap();
    let linear = LinearDensityHydrostaticMedium::new(&water, 0.0, gravity(), origin()).unwrap();
    let source = body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = cube();
    for z in [2.0, 0.2, -2.0] {
        let state = body.state_from_kinematics(pose(z), Twist::zero()).unwrap();
        assert_eq!(
            linear.load(&mesh, &body, state, || false).unwrap(),
            uniform.load(&mesh, &body, state, || false).unwrap()
        );
        let controls = HydrostaticStepControls {
            duration_s: 0.001,
            maximum_energy_defect_j: 0.1,
        };
        assert_eq!(
            linear
                .step(&mesh, &body, state, controls, || false)
                .unwrap(),
            uniform
                .step(&mesh, &body, state, controls, || false)
                .unwrap()
        );
    }
}

#[test]
fn dry_partial_and_full_immersion_match_independent_slab_integrals() {
    let water = fluid(1000.0);
    let medium = LinearDensityHydrostaticMedium::new(&water, 400.0, gravity(), origin()).unwrap();
    let source = body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = cube();
    for z in [2.0_f64, 1.0, 0.4, 0.0, -1.0, -2.0] {
        let (a, b) = ((-z - 1.0).max(0.0), (1.0 - z).max(0.0));
        let volume = 4.0 * (b - a);
        let first = 2.0 * (b * b - a * a);
        let second = (4.0 / 3.0) * (b * b * b - a * a * a);
        let mass = 1000.0 * volume + 400.0 * first;
        let state = body.state_from_kinematics(pose(z), Twist::zero()).unwrap();
        let load = medium.load(&mesh, &body, state, || false).unwrap();
        close(load.displaced().volume_m3(), volume);
        close(load.displaced_mass_kg(), mass);
        close(load.wrench().force_world.z, 10.0 * mass);
        close(
            load.pressure_potential_j(),
            10.0 * (1000.0 * first + 200.0 * second),
        );
        close(load.wrench().torque_body.norm_squared(), 0.0);
        if mass > 0.0 {
            vec_close(
                load.center_of_buoyancy_world().unwrap(),
                Vec3::new(0.0, 0.0, -(1000.0 * first + 400.0 * second) / mass),
            );
        } else {
            assert_eq!(load.center_of_buoyancy_world(), None);
        }
    }
    close(medium.density_at_depth_kg_m3(0.0).unwrap(), 1000.0);
    close(medium.density_at_depth_kg_m3(2.0).unwrap(), 1800.0);
}

#[test]
fn all_tetrahedron_orders_match_independent_stratified_slice_moments() {
    let water = fluid(1000.0);
    let vertices = [
        origin(),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];
    let card = density_card("synthetic-tet", 1500.0, 1500.0);
    let source = resolve_body(
        &[RegionBinding {
            key: "solid",
            material: &card,
            point: &point(300.0),
            geometry: VolumeMoments::tetrahedron(vertices).unwrap(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap();
    let body = source.principal_body(1e-12, || false).unwrap();
    let state = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    for h in [0.25_f64, 0.5, 0.75] {
        let medium =
            LinearDensityHydrostaticMedium::new(&water, 400.0, gravity(), Vec3::new(0.0, 0.0, h))
                .unwrap();
        // The horizontal cross section has area (1-z)^2/2 and x-mean (1-z)/3.
        let volume = (1.0 - (1.0 - h).powi(3)) / 6.0;
        let x = (1.0 - (1.0 - h).powi(4)) / 24.0;
        let z = h * h / 4.0 - h.powi(3) / 3.0 + h.powi(4) / 8.0;
        let zz = h.powi(3) / 6.0 - h.powi(4) / 4.0 + h.powi(5) / 10.0;
        let xz = h * h / 12.0 - h.powi(3) / 6.0 + h.powi(4) / 8.0 - h.powi(5) / 30.0;
        let depth = h * volume - z;
        let mass = 1000.0 * volume + 400.0 * depth;
        let cx = (1000.0 * x + 400.0 * (h * x - xz)) / mass;
        let cz = (1000.0 * z + 400.0 * (h * z - zz)) / mass;
        let potential = 10.0 * (1000.0 * depth + 200.0 * (h * h * volume - 2.0 * h * z + zz));
        for a in 0..4 {
            for b in 0..4 {
                for c in 0..4 {
                    for d in 0..4 {
                        if a == b || a == c || a == d || b == c || b == d || c == d {
                            continue;
                        }
                        let mesh = DisplacementMesh::new(
                            vec![[vertices[a], vertices[b], vertices[c], vertices[d]]],
                            1,
                            || false,
                        )
                        .unwrap();
                        let load = medium.load(&mesh, &body, state, || false).unwrap();
                        close(load.displaced_mass_kg(), mass);
                        vec_close(
                            load.center_of_buoyancy_world().unwrap(),
                            Vec3::new(cx, cx, cz),
                        );
                        close(load.pressure_potential_j(), potential);
                    }
                }
            }
        }
    }
}

#[test]
fn stratified_force_and_torque_are_the_actual_pressure_potential_gradient() {
    let water = fluid(1000.0);
    let medium = LinearDensityHydrostaticMedium::new(&water, 400.0, gravity(), origin()).unwrap();
    let sides = Vec3::new(1.2, 1.8, 0.8);
    let source = body(800.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.15, 0.3, 0.1))).unwrap();
    for z in [0.05, -1.0] {
        let state = body
            .state_from_kinematics(
                Se3::from_parts(rotation, Vec3::new(0.1, -0.2, z)).unwrap(),
                Twist::zero(),
            )
            .unwrap();
        let load = medium.load(&mesh, &body, state, || false).unwrap();
        let h = 1e-6;
        for axis in 0..3 {
            let unit = match axis {
                0 => fs_mbd::Vec3::new(1.0, 0.0, 0.0),
                1 => fs_mbd::Vec3::new(0.0, 1.0, 0.0),
                _ => fs_mbd::Vec3::new(0.0, 0.0, 1.0),
            };
            for rotational in [false, true] {
                let energy = |sign: f64| {
                    let old = state.pose();
                    let pos = if rotational {
                        old.position_world()
                    } else {
                        old.position_world().add(unit.scale(sign * h))
                    };
                    let orient = if rotational {
                        old.orientation().right_exp(unit.scale(sign * h)).unwrap()
                    } else {
                        old.orientation()
                    };
                    let changed = RigidBodyState::new(
                        Pose::new(pos, orient).unwrap(),
                        state.linear_momentum_world(),
                        state.angular_momentum_body(),
                    )
                    .unwrap();
                    medium
                        .load(&mesh, &body, changed, || false)
                        .unwrap()
                        .pressure_potential_j()
                };
                let force = if rotational {
                    load.wrench().torque_body
                } else {
                    load.wrench().force_world
                };
                let expected = force.dot(unit);
                let measured = -(energy(1.0) - energy(-1.0)) / (2.0 * h);
                assert!(
                    (measured - expected).abs() < 2e-5,
                    "{z}, {axis}, {rotational}: {measured} != {expected}"
                );
            }
        }
        // A nonuniform density shifts the actual center of buoyancy even when
        // the full volume is immersed; a center-density force is insufficient.
        if z == -1.0 {
            assert!(load.wrench().torque_body.norm_squared() > 1.0);
        }
    }
}

#[test]
fn fully_immersed_motion_restores_neutral_density_depth_with_second_order_convergence() {
    let water = fluid(1000.0);
    let medium = LinearDensityHydrostaticMedium::new(&water, 250.0, gravity(), origin()).unwrap();
    let sides = Vec3::new(1.0, 1.0, 0.5);
    let source = body(1500.0, sides);
    let body = source.principal_body(1e-12, || false).unwrap();
    let mesh = DisplacementMesh::cuboid(sides, || false).unwrap();
    // rho_body = rho_surface + beta*2. Thus z=-2 is the fully submerged
    // equilibrium and omega^2 = beta*g/rho_body, without added fluid mass.
    let exact = -2.0 + 0.05 * ((5.0_f64 / 3.0).sqrt() * 1.3).cos();
    let mut errors = Vec::new();
    for (dt, n) in [(0.02, 65), (0.01, 130), (0.005, 260)] {
        let mut state = body
            .state_from_kinematics(pose(-1.95), Twist::zero())
            .unwrap();
        for _ in 0..n {
            let result = medium
                .step(
                    &mesh,
                    &body,
                    state,
                    HydrostaticStepControls {
                        duration_s: dt,
                        maximum_energy_defect_j: 0.01,
                    },
                    || false,
                )
                .unwrap();
            state = result.dynamics.state_after;
            close(result.load_after.displaced().volume_m3(), 0.5);
            close(body.properties().mass(), 750.0);
        }
        errors.push((body.world_from_body(state).unwrap().translation().z - exact).abs());
    }
    for pair in errors.windows(2) {
        assert!((3.8..4.2).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn source_is_frozen_and_invalid_profiles_or_cancelled_loads_publish_nothing() {
    struct Counted(std::cell::Cell<usize>);
    impl FluidTransport for Counted {
        fn transport_state(&self) -> Result<FluidTransportState, fs_material::MaterialError> {
            self.0.set(self.0.get() + 1);
            Ok(fluid(1000.0))
        }
    }
    let source = Counted(std::cell::Cell::new(0));
    let medium = LinearDensityHydrostaticMedium::new(&source, 400.0, gravity(), origin()).unwrap();
    assert!(std::ptr::eq(medium.surface_source(), &source));
    close(medium.surface_transport().density_kg_m3(), 1000.0);
    let material = body(500.0, Vec3::new(2.0, 2.0, 2.0));
    let body = material.principal_body(1e-12, || false).unwrap();
    let state = body
        .state_from_kinematics(Se3::identity(), Twist::zero())
        .unwrap();
    let mesh = cube();
    let mut calls = 0;
    let expected = medium
        .load(&mesh, &body, state, || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in 1..=calls {
        let mut count = 0;
        assert!(matches!(
            medium.load(&mesh, &body, state, || {
                count += 1;
                count == stop
            }),
            Err(HydrostaticError::Cancelled)
        ));
    }
    assert_eq!(
        medium.load(&mesh, &body, state, || false).unwrap(),
        expected
    );
    assert_eq!(source.0.get(), 1);
    for gradient in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(
            LinearDensityHydrostaticMedium::new(&source, gradient, gravity(), origin()).is_err()
        );
    }
    for depth in [-1.0, f64::INFINITY, f64::NAN, f64::MAX] {
        assert!(medium.density_at_depth_kg_m3(depth).is_err());
    }
}
