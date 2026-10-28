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
