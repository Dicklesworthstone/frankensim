//! G1/G3 checks of the variable-eddy-viscosity stress, independently of the
//! equilibrium wall-law fit. Frozen closure fields isolate momentum physics.

use super::*;

fn fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0,
        kinematic_viscosity_m2_s: 1.0,
    }
}

#[test]
fn eddy_viscosity_gradient_supplies_reciprocal_shear_force() {
    // For u=(U*y,0,0), nu_t=b*(1+cos(pi*x)), the COMPLETE stress has
    // div(tau)=(0,-rho*b*pi*sin(pi*x)*U,0). Balancing that physical force
    // must leave Couette flow. A scalar viscosity Laplacian misses it.
    let mut errors = Vec::new();
    for n in [6, 12, 24] {
        let domain = VoxelDomain::new(n, n, 1, 1.0 / n as f64).unwrap();
        let (speed, contrast) = (0.001, 0.1);
        let mut faces = [FvBoundary::Symmetry; 6];
        faces[0] = FvBoundary::Outlet;
        faces[1] = FvBoundary::Outlet;
        faces[2] = FvBoundary::wall();
        faces[3] = FvBoundary::Wall {
            velocity: [speed, 0.0, 0.0],
        };
        let mut config = SimpleConfig::new(faces);
        config.tolerance = 1e-11;
        config.pressure_tolerance = 1e-12;
        config.scheme = ConvectionScheme::Upwind;
        let fluid = fluid();
        let mut solver = Solver::new(&domain, &fluid, &config);
        // Prescribe the eddy field rather than applying a wall-law iteration:
        // this test verifies the stress operator, not the LVEL model itself.
        solver.eddy_viscosity = (0..domain.cell_count())
            .map(|c| {
                let x = (domain.coords(c)[0] as f64 + 0.5) / n as f64;
                contrast * (1.0 + (core::f64::consts::PI * x).cos())
            })
            .collect();
        solver.wall_ratio = solver
            .eddy_viscosity
            .iter()
            .map(|nu_t| 1.0 + nu_t)
            .collect();
        solver.wall_adjacent = vec![false; domain.cell_count()];
        solver.set_force(
            (0..domain.cell_count())
                .map(|c| {
                    let x = (domain.coords(c)[0] as f64 + 0.5) / n as f64;
                    [
                        0.0,
                        contrast
                            * core::f64::consts::PI
                            * (core::f64::consts::PI * x).sin()
                            * speed,
                        0.0,
                    ]
                })
                .collect(),
        );
        let flow = solve_steady(solver, &fluid, &config, &CancelGate::new()).unwrap();
        let error = flow
            .velocity_m_s
            .iter()
            .enumerate()
            .fold(0.0_f64, |error, (c, v)| {
                let y = (domain.coords(c)[1] as f64 + 0.5) / n as f64;
                error.max((v[0] - speed * y).abs()).max(v[1].abs())
            });
        assert!(
            error < 0.002 * speed,
            "full eddy-stress error {error} at n={n}"
        );
        assert!(flow.report.max_divergence_m3_s < 1e-12);
        errors.push(error);
    }
    assert!(
        errors.windows(2).all(|e| e[1] < 0.65 * e[0]),
        "refinement: {errors:?}"
    );
}

#[test]
fn eddy_stress_uses_one_coefficient_for_reciprocal_components() {
    let uniform = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let graded =
        VoxelDomain::graded([vec![1.0, 2.0, 4.0], vec![2.0, 3.0, 1.0], vec![1.0]]).unwrap();
    for (domain, expected) in [(uniform, 3.0), (graded, 52.0 / 15.0)] {
        let config = SimpleConfig::new([FvBoundary::wall(); 6]);
        let fluid = fluid();
        let mut molecular = Solver::new(&domain, &fluid, &config);
        let mut eddy = Solver::new(&domain, &fluid, &config);
        let total: Vec<_> = (1..=9).map(|c| c as f64).collect();
        molecular.molecular_viscosity = total.clone();
        eddy.eddy_viscosity = total.iter().map(|mu| mu - 1.0).collect();
        eddy.wall_ratio = total;
        assert!(eddy.symmetric_stress());
        let ab = eddy.edge_mu([0, 0, 0], [1, 0, 0], 1, true);
        let ba = eddy.edge_mu([0, 0, 0], [0, 1, 0], 0, true);
        assert_eq!(
            ab, ba,
            "reciprocal components must share one stress coefficient"
        );
        assert!((ab - expected).abs() < 1e-15);
        assert_eq!(
            eddy.wall_mu([0, 0, 0], [1, 0, 0]),
            molecular.wall_mu([0, 0, 0], [1, 0, 0])
        );
    }
}

#[test]
fn lvel_refuses_an_unbounded_wall_distance() {
    let domain = VoxelDomain::new(3, 3, 3, 1.0).unwrap();
    let mut config = SimpleConfig::new([FvBoundary::Symmetry; 6]);
    config.faces[0] = FvBoundary::Inlet {
        velocity: [1.0, 0.0, 0.0],
    };
    config.faces[1] = FvBoundary::Outlet;
    config.turbulence = Turbulence::Lvel;
    assert!(matches!(
        simple_flow(&domain, &fluid(), &config, &CancelGate::new()),
        Err(ChtError::InvalidInput {
            field: "simple.turbulence",
            ..
        })
    ));
    config.faces[2] = FvBoundary::wall();
    assert!(admit(&domain, &fluid(), &config).is_ok());
    config.faces[2] = FvBoundary::Symmetry;
    config.turbulence = Turbulence::Laminar;
    assert!(admit(&domain, &fluid(), &config).is_ok());
}

#[test]
fn comoving_lvel_walls_produce_no_turbulent_stress_or_heat_transport() {
    // G0/G3: an exact translation is the same physical state as rest. All
    // wall identities/axes must carry the same velocity reference as distance.
    let domain = VoxelDomain::graded([
        vec![0.125, 0.25, 0.625],
        vec![0.25, 0.5, 0.25],
        vec![0.5, 0.5],
    ])
    .unwrap();
    for axis in 0..3 {
        for speed in [-8.0, 0.0, 8.0] {
            let mut velocity = [0.0; 3];
            velocity[axis] = speed;
            let mut config = SimpleConfig::new([FvBoundary::Wall { velocity }; 6]);
            config.faces[2 * axis] = FvBoundary::Outlet;
            config.faces[2 * axis + 1] = FvBoundary::Outlet;
            config.turbulence = Turbulence::Lvel;
            let fluid = fluid();
            let mut solver = Solver::new(&domain, &fluid, &config);
            solver.vel[axis].fill(speed);
            solver.update_turbulence();
            assert!(solver.eddy_viscosity.iter().all(|&nu| nu == 0.0));
            assert!(solver.wall_ratio.iter().all(|&ratio| ratio == 1.0));
            assert!(solver.heat_eddy_viscosity().iter().all(|&nu| nu == 0.0));
            assert!(
                solver
                    .wall_boundary
                    .iter()
                    .all(|&side| { side.is_some_and(|side| usize::from(side) / 2 != axis) })
            );
            let flow = solver.finish_projected(&fluid, 0, (0.0, 0.0));
            assert!(flow.eddy_conductivity(&fluid).iter().all(|&k| k == 0.0));
        }
    }
}

#[test]
fn lvel_moving_wall_couette_is_invariant_under_galilean_boosts() {
    // G3: the dyadic geometry makes the uniform boosts exact in f64, so
    // changes in the closure cannot be excused as cancellation roundoff.
    let domain = VoxelDomain::graded([
        vec![0.25, 0.75],
        vec![0.125, 0.125, 0.25, 0.5],
        vec![0.5, 0.5],
    ])
    .unwrap();
    let fluid = FluidProperties {
        kinematic_viscosity_m2_s: 0.001,
        ..fluid()
    };
    let mut reference = None;
    for boost in [0.0, -8.0, 8.0] {
        let mut faces = [FvBoundary::Symmetry; 6];
        faces[0] = FvBoundary::Outlet;
        faces[1] = FvBoundary::Outlet;
        faces[2] = FvBoundary::Wall {
            velocity: [boost + 0.25, 0.0, 0.0],
        };
        faces[3] = FvBoundary::Wall {
            velocity: [boost + 1.25, 0.0, 0.0],
        };
        let mut config = SimpleConfig::new(faces);
        config.turbulence = Turbulence::Lvel;
        let mut solver = Solver::new(&domain, &fluid, &config);
        for face in 0..solver.vel[0].len() {
            let at = solver.comps[0].coords(face);
            let y = domain.face_coord(1, at[1]) + 0.5 * domain.width(1, at[1]);
            solver.vel[0][face] = boost + 0.25 + y;
        }
        solver.update_turbulence();
        for c in 0..domain.cell_count() {
            let at = domain.coords(c);
            let y = domain.center(at[0], at[1], at[2])[1];
            let distance = y.min(1.0 - y);
            assert_eq!(solver.wall_distance[c], distance);
            assert_eq!(solver.wall_boundary[c], Some(if y < 0.5 { 2 } else { 3 }));
            let (tangent, secant) =
                law_of_the_wall_ratios(distance * distance / fluid.kinematic_viscosity_m2_s);
            assert_eq!(
                solver.eddy_viscosity[c],
                tangent * fluid.kinematic_viscosity_m2_s
            );
            assert_eq!(solver.wall_ratio[c], secant);
        }
        let observed = (
            solver.eddy_viscosity.clone(),
            solver.wall_ratio.clone(),
            solver.heat_eddy_viscosity(),
        );
        if let Some(expected) = &reference {
            assert_eq!(&observed, expected, "Galilean boost {boost}");
        } else {
            assert!(observed.0.iter().any(|&nu| nu > 0.0));
            reference = Some(observed);
        }
    }
}

#[test]
fn lvel_uses_stationary_internal_solids_and_the_nearest_moving_domain_wall() {
    // G0: a moving outer wall must not move the frame of a nearby solid.
    let domain = VoxelDomain::from_fn(3, 5, 1, 1.0, |p| {
        if p[1] < 1.0 {
            crate::conjugate::Voxel::Solid(0)
        } else {
            crate::conjugate::Voxel::Fluid
        }
    })
    .unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    faces[3] = FvBoundary::Wall {
        velocity: [3.0, 0.0, 0.0],
    };
    let mut config = SimpleConfig::new(faces);
    config.turbulence = Turbulence::Lvel;
    let fluid = fluid();
    let mut solver = Solver::new(&domain, &fluid, &config);
    solver.vel[0].fill(3.0);
    solver.update_turbulence();
    for c in 0..domain.cell_count() {
        let y = domain.coords(c)[1];
        if y == 0 {
            assert_eq!(solver.wall_distance[c], 0.0);
            assert_eq!(solver.eddy_viscosity[c], 0.0);
        } else if y <= 2 {
            assert_eq!(solver.wall_boundary[c], None);
            let distance = y as f64 - 0.5;
            assert_eq!(solver.wall_distance[c], distance);
            assert_eq!(
                solver.eddy_viscosity[c],
                crate::conjugate::turbulence::eddy_viscosity_ratio(3.0 * distance)
            );
            assert!(solver.eddy_viscosity[c] > 0.0);
        } else {
            assert_eq!(solver.wall_boundary[c], Some(3));
            assert_eq!(solver.wall_distance[c], 4.5 - y as f64);
            assert_eq!(solver.eddy_viscosity[c], 0.0);
            assert_eq!(solver.wall_ratio[c], 1.0);
        }
    }
}

#[test]
fn graded_obstacle_surface_controls_lvel_distance_and_wall_velocity_frame() {
    // A small moving domain-wall ghost has a nearer centre than the wide
    // stationary obstacle. Its plane is farther away, however. Using the
    // nearest centre would select the comoving wall, erase nu_t, and use
    // L=2.5 instead of the actual stationary-solid distance L=1.5.
    let mut domain =
        VoxelDomain::graded([vec![1.0, 1.0, 8.0, 0.5], vec![1.0, 2.0], vec![1.0]]).unwrap();
    domain.set(2, 0, 0, crate::conjugate::Voxel::Solid(0));
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    faces[3] = FvBoundary::Wall {
        velocity: [3.0, 0.0, 0.0],
    };
    let mut config = SimpleConfig::new(faces);
    config.turbulence = Turbulence::Lvel;
    let fluid = fluid();
    let mut solver = Solver::new(&domain, &fluid, &config);
    solver.vel[0].fill(3.0);
    solver.update_turbulence();
    let cell = domain.index(0, 0, 0);
    assert_eq!(solver.wall_distance[cell], 1.5);
    assert_eq!(solver.wall_boundary[cell], None);
    let (tangent, secant) = law_of_the_wall_ratios(3.0 * 1.5);
    assert_eq!(solver.eddy_viscosity[cell], tangent);
    assert_eq!(solver.wall_ratio[cell], secant);
    assert!(solver.eddy_viscosity[cell] > 0.0);
}

#[test]
fn moving_wall_channel_recovers_plug_flow_without_artificial_heat_diffusion() {
    // G0: exercise the public solve and returned thermal coefficients, not
    // only a manufactured closure field. The belts and inflow move together.
    let domain = VoxelDomain::new(6, 4, 1, 0.25).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Inlet {
        velocity: [0.5, 0.0, 0.0],
    };
    faces[1] = FvBoundary::Outlet;
    faces[2] = FvBoundary::Wall {
        velocity: [0.5, 0.0, 0.0],
    };
    faces[3] = faces[2];
    let mut config = SimpleConfig::new(faces);
    config.turbulence = Turbulence::Lvel;
    config.tolerance = 1e-10;
    config.pressure_tolerance = 1e-12;
    let fluid = FluidProperties {
        kinematic_viscosity_m2_s: 0.01,
        ..fluid()
    };
    let flow = simple_flow(&domain, &fluid, &config, &CancelGate::new()).unwrap();
    for velocity in &flow.velocity_m_s {
        assert!((velocity[0] - 0.5).abs() < 1e-9, "{velocity:?}");
        assert!(velocity[1].abs() < 1e-9 && velocity[2].abs() < 1e-9);
    }
    assert!(flow.eddy_viscosity_m2_s.iter().all(|&nu| nu < 1e-9));
    assert!(flow.eddy_conductivity(&fluid).iter().all(|&k| k < 1e-9));
    assert!(flow.report.max_divergence_m3_s < 1e-11);
}
