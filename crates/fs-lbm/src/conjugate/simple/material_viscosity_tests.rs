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

#[test]
fn closed_rotating_vortex_has_the_analytic_geostrophic_pressure_and_no_extra_heat() {
    use super::viscous_heat::{ViscousHeatControls, creeping_flow_with_body_force_heat};

    // A discrete streamfunction psi=A*sin(pi*x)*sin(pi*y) on MAC vertices
    // gives exactly divergence-free free-slip velocities. The discrete
    // Laplacian eigenvalue is 8*sin(pi/(2*n))^2/h^2. Its Coriolis force is
    // exactly balanced by p=-2*rho*Omega*A*cos(pi/(2*n))^2*sin(pi*x)*sin(pi*y)
    // at cell centres. The forcing below supplies only viscous loss.
    let n = 6;
    let h = 1.0 / n as f64;
    let domain = VoxelDomain::new(n, n, 1, h).unwrap();
    let (rho, mu, amplitude) = (1.7, 0.6, 2e-4);
    let fluid = FluidProperties {
        density_kg_m3: rho,
        kinematic_viscosity_m2_s: mu / rho,
        ..fluid()
    };
    let theta = core::f64::consts::PI / (2 * n) as f64;
    let (s, c) = (theta.sin(), theta.cos());
    let velocity_scale = amplitude * 2.0 * s / h;
    let eigenvalue = 8.0 * s * s / (h * h);
    let mode = |cell| {
        let at = domain.coords(cell);
        let x = core::f64::consts::PI * (at[0] as f64 + 0.5) * h;
        let y = core::f64::consts::PI * (at[1] as f64 + 0.5) * h;
        (x, y)
    };
    let force: Vec<_> = (0..domain.cell_count())
        .map(|cell| {
            let (x, y) = mode(cell);
            let factor = mu * eigenvalue * velocity_scale / c;
            [factor * x.sin() * y.cos(), -factor * x.cos() * y.sin(), 0.0]
        })
        .collect();
    let mut unrotated_heat: Option<f64> = None;
    for omega in [0.0, 0.7, -0.7] {
        let mut config = config([FvBoundary::Symmetry; 6]);
        config.coriolis_angular_velocity_rad_s = [0.0, 0.0, omega];
        let (flow, power) = creeping_flow_with_body_force_heat(
            &domain, &fluid, &vec![mu; domain.cell_count()], &force, None, &config,
            ViscousHeatControls { maximum_reynolds: 0.1, power_tolerance_w: 1e-11 },
            &CancelGate::new(),
        ).unwrap();
        let pressure = |cell| {
            let (x, y) = mode(cell);
            -2.0 * rho * omega * amplitude * c * c * x.sin() * y.sin()
        };
        for (cell, velocity) in flow.velocity_m_s.iter().enumerate() {
            let (x, y) = mode(cell);
            let expected = [
                velocity_scale * c * x.sin() * y.cos(),
                -velocity_scale * c * x.cos() * y.sin(),
                0.0,
            ];
            for axis in 0..3 {
                assert!((velocity[axis] - expected[axis]).abs() < 2e-9 * velocity_scale,
                    "Omega={omega}, cell={cell}: {velocity:?} != {expected:?}");
            }
            let actual = flow.pressure_pa[cell] - flow.pressure_pa[0];
            assert!((actual - (pressure(cell) - pressure(0))).abs() < 2e-11,
                "Omega={omega}, pressure cell={cell}: {actual}");
        }
        assert!(flow.report.max_divergence_m3_s < 1e-13);
        assert!(flow.report.momentum_residual <= config.tolerance);
        assert!((power.supplied_body_power_w - power.heat_w).abs() < 1e-11);
        assert_eq!(power.darcy_heat_w, 0.0);
        assert_eq!(power.forchheimer_heat_w, 0.0);
        if let Some(expected) = unrotated_heat {
            assert!((power.heat_w - expected).abs() < 1e-11);
        } else {
            unrotated_heat = Some(power.heat_w);
        }
    }
}

#[test]
fn arbitrary_axis_coriolis_uses_the_graded_mac_transpose_without_work() {
    let mut domain = VoxelDomain::graded([
        vec![0.1, 0.4, 0.2], vec![0.3, 0.2, 0.6], vec![0.2, 0.7],
    ]).unwrap();
    domain.set(1, 1, 0, crate::conjugate::Voxel::Solid(0));
    let mut config = config([FvBoundary::wall(); 6]);
    config.coriolis_angular_velocity_rad_s = [0.7, -1.1, 0.3];
    let mut fluid = fluid();
    fluid.density_kg_m3 = 2.3;
    admit(&domain, &fluid, &config).unwrap();
    let mut solver = Solver::new(&domain, &fluid, &config);
    let mut x = solver.vel.clone();
    let mut y = solver.vel.clone();
    for axis in 0..3 {
        for &index in &solver.comps[axis].unknowns {
            x[axis][index] = ((index * 7 + axis * 3) % 17) as f64 / 13.0 - 0.5;
            y[axis][index] = ((index * 11 + axis * 5) % 19) as f64 / 17.0 - 0.4;
        }
    }
    let pairing = |solver: &Solver<'_>, test: &[Vec<f64>; 3]| {
        let mut sum = 0.0;
        let mut absolute = 0.0;
        for axis in 0..3 {
            for &index in &solver.comps[axis].unknowns {
                let plus = solver.comps[axis].coords(index);
                let mut minus = plus;
                minus[axis] -= 1;
                let value = test[axis][index] * solver.coriolis_face_force(axis, minus, plus);
                sum += value;
                absolute += value.abs();
            }
        }
        (sum, absolute)
    };
    solver.vel = y.clone();
    let (xy, xy_scale) = pairing(&solver, &x);
    let (yy, yy_scale) = pairing(&solver, &y);
    let y_cells: Vec<_> = (0..domain.cell_count()).map(|c| solver.cell_velocity(c)).collect();
    solver.vel = x.clone();
    let (yx, yx_scale) = pairing(&solver, &y);
    let (xx, xx_scale) = pairing(&solver, &x);
    let omega = config.coriolis_angular_velocity_rad_s;
    let cell_pairing: f64 = y_cells.iter().enumerate().map(|(cell, v)| {
        let u = solver.cell_velocity(cell);
        let cross = [
            omega[1] * v[2] - omega[2] * v[1],
            omega[2] * v[0] - omega[0] * v[2],
            omega[0] * v[1] - omega[1] * v[0],
        ];
        2.0 * fluid.density_kg_m3 * domain.volume(cell)
            * (u[0] * cross[0] + u[1] * cross[1] + u[2] * cross[2])
    }).sum();
    let tolerance = 64.0 * f64::EPSILON;
    assert!(xy_scale > 0.01);
    assert!((xy - cell_pairing).abs() <= tolerance * xy_scale);
    assert!((xy + yx).abs() <= tolerance * (xy_scale + yx_scale));
    assert!(xx.abs() <= tolerance * xx_scale);
    assert!(yy.abs() <= tolerance * yy_scale);
}

#[test]
fn coriolis_refuses_nonfinite_open_and_foreign_inertia_without_mutating_faces() {
    let domain = VoxelDomain::new(3, 3, 2, 0.25).unwrap();
    let fluid = fluid();
    let mut config = config([FvBoundary::wall(); 6]);
    for omega in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, f64::MAX] {
        config.coriolis_angular_velocity_rad_s = [0.0, omega, 0.0];
        assert!(admit(&domain, &fluid, &config).is_err());
    }
    config.faces[0] = FvBoundary::Outlet;
    config.coriolis_angular_velocity_rad_s = [0.1, 0.0, 0.0];
    assert!(admit(&domain, &fluid, &config).is_err());
    config.coriolis_angular_velocity_rad_s = [-0.0, 0.0, 0.0];
    admit(&domain, &fluid, &config).unwrap();
    config.faces[0] = FvBoundary::wall();
    config.coriolis_angular_velocity_rad_s = [0.1, -0.2, 0.3];
    let mut solver = Solver::new(&domain, &fluid, &config);
    solver.start_time(0.1, TimeScheme::BackwardEuler);
    solver.time.as_mut().unwrap().face_mass_kg =
        Some(std::array::from_fn(|a| vec![1.0; solver.vel[a].len()]));
    let original = solver.vel.clone();
    let pressure = solver.pressure.clone();
    assert!(matches!(solver.momentum(0, true, &CancelGate::new()),
        Err(ChtError::InvalidInput { field: "simple.coriolis_inertia", .. })));
    assert_eq!(solver.vel, original);
    assert_eq!(solver.pressure, pressure);
    solver.time.as_mut().unwrap().face_mass_kg = None;
    assert_eq!(solver.momentum(0, true, &CancelGate::new()).unwrap().relative, 0.0);
}

#[cfg(feature = "enthalpy-coupling")]
#[test]
fn rotating_material_momentum_retries_after_an_unconverged_step() {
    use crate::conjugate::unsteady::MaterialMomentum;
    let domain = VoxelDomain::new(3, 3, 2, 0.25).unwrap();
    let fluid = fluid();
    let mut faces = [FvBoundary::wall(); 6];
    faces[3] = FvBoundary::Wall { velocity: [0.01, 0.0, 0.0] };
    let mut config = config(faces);
    config.coriolis_angular_velocity_rad_s = [0.3, -0.4, 0.5];
    let gate = CancelGate::new();
    let mut owner = MaterialMomentum::at_rest(&domain, &fluid, &config, 0.0, &gate).unwrap();
    let initial = owner.clone();
    let viscosity = vec![1.0; domain.cell_count()];
    let advance = |owner: &mut MaterialMomentum, config: &SimpleConfig| owner.advance(
        &domain, &fluid, config, &viscosity, None, None, None, None, None,
        0.0, 0.1, 0.1, &gate,
    );
    let mut short = config.clone();
    short.max_iterations = 1;
    assert!(advance(&mut owner, &short).is_err());
    assert_eq!(owner.time_s, initial.time_s);
    assert_eq!(owner.flow, initial.flow);
    assert_eq!(owner.record, initial.record);
    advance(&mut owner, &config).unwrap();
    let mut fresh = initial;
    advance(&mut fresh, &config).unwrap();
    assert_eq!(owner.flow, fresh.flow);
    assert_eq!(owner.record, fresh.record);
    assert_eq!(owner.time_s, 0.1);
}
