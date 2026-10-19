//! Analytic shear and Darcy balances through the actual MAC momentum owner.
use super::*;

fn fluid() -> FluidProperties {
    FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0,
        kinematic_viscosity_m2_s: 1.0,
    }
}

fn config(faces: [FvBoundary; 6]) -> SimpleConfig {
    let mut config = SimpleConfig::new(faces);
    config.tolerance = 1e-11;
    config.pressure_tolerance = 1e-12;
    config.scheme = ConvectionScheme::Upwind;
    config
}

#[test]
fn layered_couette_flow_continuously_transmits_shear_traction() {
    let domain = VoxelDomain::new(4, 16, 1, 1.0 / 16.0).unwrap();
    let speed = 0.01;
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    faces[2] = FvBoundary::wall();
    faces[3] = FvBoundary::Wall {
        velocity: [speed, 0.0, 0.0],
    };
    let config = config(faces);
    let mu: Vec<_> = (0..domain.cell_count())
        .map(|c| if domain.coords(c)[1] < 8 { 1.0 } else { 4.0 })
        .collect();
    let flow =
        simple_flow_with_viscosity(&domain, &fluid(), &mu, &config, &CancelGate::new()).unwrap();
    // tau = U / integral_0^H (1/mu) dy. The piecewise linear profile is
    // exactly represented by the half-cell resistances of the MAC stencil.
    let traction = speed / (0.5 / 1.0 + 0.5 / 4.0);
    for (c, velocity) in flow.velocity_m_s.iter().enumerate() {
        let y = (domain.coords(c)[1] as f64 + 0.5) / 16.0;
        let compliance = if y <= 0.5 { y } else { 0.5 + (y - 0.5) / 4.0 };
        let expected = traction * compliance;
        assert!(
            (velocity[0] - expected).abs() < 2e-8 * speed,
            "{c}: {velocity:?} != {expected}"
        );
        assert!(velocity[1].abs() < 2e-8 * speed);
    }
    assert!(flow.report.max_divergence_m3_s < 1e-12);
}

#[test]
fn viscosity_gradient_supplies_reciprocal_shear_force() {
    // u=(U*y,0,0), mu=1+b*cos(pi*x) has div(tau)=(0,mu'(x)*U,0). An opposing
    // prescribed body force makes it a steady Couette solution at p=0.
    // Scalar component diffusion alone omits this entire transverse load.
    let mut errors = Vec::new();
    for n in [6, 12] {
        let domain = VoxelDomain::new(n, n, 1, 1.0 / n as f64).unwrap();
        let (speed, contrast) = (0.001, 0.1);
        let mut faces = [FvBoundary::Symmetry; 6];
        faces[0] = FvBoundary::Outlet;
        faces[1] = FvBoundary::Outlet;
        faces[2] = FvBoundary::wall();
        faces[3] = FvBoundary::Wall {
            velocity: [speed, 0.0, 0.0],
        };
        let config = config(faces);
        let fluid = fluid();
        let mut solver = Solver::new(&domain, &fluid, &config);
        solver.molecular_viscosity = (0..domain.cell_count())
            .map(|c| {
                let x = (domain.coords(c)[0] as f64 + 0.5) / n as f64;
                1.0 + contrast * (core::f64::consts::PI * x).cos()
            })
            .collect();
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
        assert!(error < 0.002 * speed, "full-stress error {error} at n={n}");
        errors.push(error);
    }
    assert!(errors[1] < 0.65 * errors[0], "refinement errors {errors:?}");
}

#[test]
fn local_material_viscosity_also_updates_original_darcy_resistance() {
    let domain = VoxelDomain::new(4, 1, 1, 1.0).unwrap();
    let mut faces = [FvBoundary::Symmetry; 6];
    faces[0] = FvBoundary::Outlet;
    faces[1] = FvBoundary::Outlet;
    let mut config = config(faces);
    config.resistances.push(FlowResistance::Volume {
        lo: [0, 0, 0],
        hi: [4, 1, 1],
        permeability_m2: [1.0; 3],
        inertial_per_m: [0.0; 3],
    });
    config.internal_fans.push(InternalFan {
        patch: FacePatch {
            axis: 0,
            index: 2,
            lo: [0, 0, 0],
            hi: [4, 1, 1],
        },
        blows_positive: true,
        curve: FanCurve::new(&[(0.0, 8.0), (100.0, 8.0)]).unwrap(),
    });
    for factor in [1.0, 2.0] {
        let mu = [1.0, 2.0, 3.0, 4.0].map(|mu| factor * mu);
        let flow = simple_flow_with_viscosity(&domain, &fluid(), &mu, &config, &CancelGate::new())
            .unwrap();
        let expected = 8.0 / mu.iter().sum::<f64>();
        for face in 0..=4 {
            assert!((flow.field.flux_x(face, 0, 0) - expected).abs() < 1e-9 * expected);
        }
        assert!((flow.report.max_cell_reynolds - expected / factor).abs() < 1e-9);
    }
}

#[test]
fn shear_edge_uses_one_coefficient_for_reciprocal_stresses_and_refuses_foreign_closures() {
    let domain = VoxelDomain::new(3, 3, 1, 1.0).unwrap();
    let mut config = config([FvBoundary::wall(); 6]);
    let fluid = fluid();
    let mut solver = Solver::new(&domain, &fluid, &config);
    solver.molecular_viscosity = (1..=9).map(|value| value as f64).collect();
    assert_eq!(
        solver.edge_mu([0, 0, 0], [1, 0, 0], 1, true),
        solver.edge_mu([0, 0, 0], [0, 1, 0], 0, true)
    );
    drop(solver);
    let gate = CancelGate::new();
    assert!(simple_flow_with_viscosity(&domain, &fluid, &[1.0; 8], &config, &gate).is_err());
    config.turbulence = Turbulence::Lvel;
    assert!(simple_flow_with_viscosity(&domain, &fluid, &[1.0; 9], &config, &gate).is_err());
}
