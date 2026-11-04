//! Actual constitutive momentum on rigid moving charts, using original owners.
use super::*;
use fs_flux::{MaterialFluxError, MaterialFluxRheologyPoint, MaterialFluxSystem};
use fs_material::{MaterialError, fluid_transport::FluidTransport};
use std::cell::RefCell;

struct Response {
    query: MaterialFluxRheologyPoint,
    transport: FluidTransportState,
}
impl FluidTransport for Response {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        Ok(self.transport)
    }
}

fn profile(y: f64, offset: f64, slopes: [f64; 2]) -> f64 {
    offset + slopes[0] * y.min(0.5) + slopes[1] * (y - 0.5).max(0.0)
}

#[test]
fn g1_changing_rigid_motion_and_layered_rheology_recover_constant_physical_shear_stress() {
    let mesh = mesh(2);
    let system = FluxSystem::new(&mesh);
    let pivot = [0.3, -0.2];
    let old_motion = RigidAleMotion {
        linear_velocity_world_m_s: [-0.05, 0.02],
        angular_velocity_rad_s: 0.25,
    };
    let motion = RigidAleMotion {
        linear_velocity_world_m_s: [0.13, -0.07],
        angular_velocity_rad_s: 0.8,
    };
    let old_slopes = [0.15, 0.35];
    let slopes: [f64; 2] = [0.4, 0.2];
    let tau = 0.12;
    let dt = 0.04;
    let old_field = |y| {
        add(
            old_motion.linear_velocity_world_m_s,
            add(
                scale(spin(sub(y, pivot)), old_motion.angular_velocity_rad_s),
                [profile(y[1], 0.04, old_slopes), 0.0],
            ),
        )
    };
    let initial = history(&mesh, old_field);
    let bottom: Vec<_> = edges(&mesh)
        .into_iter()
        .filter(|&edge| mesh.edges[edge].normal[1] < -0.5)
        .collect();
    // The exact piecewise affine profile is continuous at y=1/2. Its actual
    // dynamic shear stress is tau on BOTH sides despite different shear rates.
    // These two laws thin/thicken with strain; neither freezes the initial mu.
    for (rho, exponent) in [(2.0, -0.25_f64), (7.0, 0.2_f64)] {
        let source = material(rho);
        let mut owner = RigidAleFlow::new(
            &system,
            &source,
            MaterialFluxNumerics::default(),
            pivot,
            old_motion,
            &initial,
        )
        .unwrap();
        let k = owner.kinematics(dt, motion, tolerance()).unwrap();
        let r = k.orientation_cos_sin;
        let omega = k.mesh_velocity.angular_velocity_rad_s;
        let to_reference = |x| add(pivot, inverse(r, sub(x, k.pivot_world_m)));
        let absolute = |y| {
            add(
                k.mesh_velocity.linear_velocity_world_m_s,
                rotate(
                    r,
                    add(
                        scale(spin(sub(y, pivot)), omega),
                        [profile(y[1], 0.1, slopes), 0.0],
                    ),
                ),
            )
        };
        let velocity = |x| absolute(to_reference(x));
        let traction = |_: [f64; 2], normal| {
            let n = inverse(r, normal);
            rotate(r, [tau * n[1], tau * n[0]])
        };
        let force = |x| {
            let y = to_reference(x);
            let temporal = scale(sub(absolute(y), old_field(y)), 1.0 / dt);
            let relative = [profile(y[1], 0.1, slopes), 0.0];
            // grad(u)*(u-w) = omega J(relative) for this shear profile.
            let convection = rotate(r, scale(spin(relative), omega));
            scale(add(temporal, convection), rho)
        };
        let boundaries: Vec<_> = edges(&mesh)
            .into_iter()
            .map(|edge| {
                if mesh.edges[edge].normal[1].abs() > 0.5 {
                    MaterialNavierStokesBoundary::Velocity {
                        edge,
                        value: &velocity,
                    }
                } else {
                    MaterialNavierStokesBoundary::OpenReservoir {
                        edge,
                        traction: &traction,
                        velocity: &velocity,
                    }
                }
            })
            .collect();
        let first = RefCell::new(vec![None; mesh.tris.len()]);
        let resolve = |point: MaterialFluxRheologyPoint, reference: &FluidTransportState| {
            first.borrow_mut()[point.cell].get_or_insert(point.shear_rate_per_s);
            let y = to_reference(point.position_m);
            let layer = usize::from(y[1] >= 0.5);
            let gamma = slopes[layer];
            let base = tau / gamma / (1.0 + gamma * gamma).powf(exponent);
            let mu = base * (1.0 + point.shear_rate_per_s.powi(2)).powf(exponent);
            Ok(Response {
                query: point,
                transport: FluidTransportState::try_new(
                    reference.density_kg_m3(),
                    1000.0,
                    0.5,
                    mu,
                )?,
            })
        };
        let step = owner
            .step_with_rheology(
                dt,
                motion,
                tolerance(),
                &force,
                &boundaries,
                NavierStokesControls {
                    maximum_iterations: 100,
                    relaxation: 0.7,
                    ..controls()
                },
                resolve,
                || false,
            )
            .unwrap();
        let exact = |y| inverse(r, absolute(y));
        let (error, div) = system.velocity_error(step.flow.solution.velocity_dofs(), &exact);
        assert!(error < 2e-7 && div < 2e-7, "{error} {div}");
        let origin = step.point_world([0.2, -0.4]).unwrap();
        let load = step.boundary_loads_world(&bottom, origin, &[None]).unwrap();
        vector(load.stress.force_n_per_m, rotate(r, [tau, 0.0]), 2e-7);
        near(load.stress.moment_n, -0.4 * tau, 2e-7);
        near(load.slip_dissipation_w_per_m, 0.0, 1e-12);
        let transport = step.boundary_transport_world(&bottom, origin).unwrap();
        near(transport.outward_mass_flux_kg_per_m_s, 0.0, 2e-8);
        vector(transport.outward_momentum.force_n_per_m, [0.0; 2], 2e-8);
        // Independent integral of the piecewise linear profile on the unit box.
        let mean_profile = 0.1 + 0.375 * slopes[0] + 0.125 * slopes[1];
        let mean = add(
            k.mesh_velocity.linear_velocity_world_m_s,
            rotate(
                r,
                add(
                    scale(spin(sub([0.5, 0.5], pivot)), omega),
                    [mean_profile, 0.0],
                ),
            ),
        );
        vector(step.world_momentum_n_s_per_m, scale(mean, rho), 4e-7);
        assert_eq!(step.rheology_samples().len(), mesh.tris.len());
        for (cell, sample) in step.rheology_samples().iter().enumerate() {
            let layer = usize::from(mesh.centroids[cell][1] >= 0.5);
            let point = sample.point();
            assert_eq!(point.cell, cell);
            vector(
                point.position_m,
                step.point_world(mesh.centroids[cell]).unwrap(),
                1e-13,
            );
            assert_eq!(sample.material().query, point);
            near(point.shear_rate_per_s, slopes[layer], 3e-7);
            near(
                sample.transport().dynamic_viscosity_pa_s(),
                tau / slopes[layer],
                3e-7,
            );
            near(first.borrow()[cell].unwrap(), old_slopes[layer], 1e-12);
            near(
                step.flow.solution.viscous_dissipation().cell_heating_w_m3()[cell],
                tau * slopes[layer],
                3e-7,
            );
        }
        near(
            step.flow
                .solution
                .viscous_dissipation()
                .total_power_w_per_m(),
            0.5 * tau * (slopes[0] + slopes[1]),
            3e-7,
        );
    }
}

#[test]
fn g3_free_world_stream_has_zero_strain_with_moving_spatial_material_queries() {
    let mesh = mesh(2);
    let system = FluxSystem::new(&mesh);
    let source = material(2.3);
    let pivot = [0.2, -0.3];
    let u = [0.4, -0.2];
    let initial = history(&mesh, |_| u);
    let incoming = |_| u;
    let boundaries = open(&mesh, &no_stress, Some(&incoming));
    let mut owner = RigidAleFlow::new(
        &system,
        &source,
        MaterialFluxNumerics::default(),
        pivot,
        RigidAleMotion::ZERO,
        &initial,
    )
    .unwrap();
    let queries = Cell::new(0);
    for (v, omega) in [([0.3, -0.2], 1.2), ([-0.1, 0.4], -0.7), ([0.2, 0.1], 0.0)] {
        let motion = RigidAleMotion {
            linear_velocity_world_m_s: v,
            angular_velocity_rad_s: omega,
        };
        let step = owner
            .step_with_rheology(
                0.08,
                motion,
                tolerance(),
                &zero,
                &boundaries,
                controls(),
                |point, reference| {
                    queries.set(queries.get() + 1);
                    let [x, y] = point.position_m;
                    let mu = 0.2
                        * (0.1 * x + 0.2 * y).exp()
                        * (1.0 + point.shear_rate_per_s.powi(2)).powf(-0.2);
                    Ok(Response {
                        query: point,
                        transport: FluidTransportState::try_new(
                            reference.density_kg_m3(),
                            1000.0,
                            0.5,
                            mu,
                        )?,
                    })
                },
                || false,
            )
            .unwrap();
        let local = inverse(step.orientation_cos_sin, u);
        let (error, div) = system.velocity_error(step.flow.solution.velocity_dofs(), &|_| local);
        near(
            step.flow
                .solution
                .viscous_dissipation()
                .total_power_w_per_m(),
            0.0,
            2e-16,
        );
        assert!(error < 4e-9 && div < 4e-9, "{error} {div}");
        vector(step.world_momentum_n_s_per_m, scale(u, 2.3), 4e-9);
        let count = queries.get();
        assert!(count > 0);
        for (cell, sample) in step.rheology_samples().iter().enumerate() {
            let position = step.point_world(mesh.centroids[cell]).unwrap();
            vector(sample.point().position_m, position, 1e-13);
            assert_eq!(sample.material().query, sample.point());
            near(sample.point().shear_rate_per_s, 0.0, 1e-8);
            near(
                sample.transport().dynamic_viscosity_pa_s(),
                0.2 * (0.1 * position[0] + 0.2 * position[1]).exp(),
                2e-9,
            );
        }
        assert_eq!(
            queries.get(),
            count,
            "retained owners are not source re-queries"
        );
        vector(
            step.boundary_loads_world(&edges(&mesh), [0.0; 2], &[None])
                .unwrap()
                .stress
                .force_n_per_m,
            [0.0; 2],
            1e-8,
        );
    }
    // Existing Newtonian continuation remains available on the same owner.
    let step = owner
        .step(
            0.02,
            RigidAleMotion::ZERO,
            tolerance(),
            &zero,
            &boundaries,
            controls(),
            || false,
        )
        .unwrap();
    assert!(step.rheology_samples().is_empty());
    vector(step.world_momentum_n_s_per_m, scale(u, 2.3), 4e-9);
}

#[test]
fn g3_stationary_rheology_matches_original_native_material_momentum_owner() {
    let mesh = mesh(2);
    let system = FluxSystem::new(&mesh);
    let source = material(2.0);
    let exact = |x: [f64; 2]| [0.1 + 0.3 * x[1], 0.07];
    let initial = history(&mesh, exact);
    let boundaries: Vec<_> = edges(&mesh)
        .into_iter()
        .map(|edge| MaterialNavierStokesBoundary::Velocity {
            edge,
            value: &exact,
        })
        .collect();
    let force = |_| [2.0 * 0.3 * 0.07, 0.0];
    let resolve = |point: MaterialFluxRheologyPoint, reference: &FluidTransportState| {
        FluidTransportState::try_new(
            reference.density_kg_m3(),
            1000.0,
            0.5,
            0.5 * (1.0 + 4.0 * point.shear_rate_per_s.powi(2)).powf(-0.2),
        )
    };
    let original = MaterialFluxSystem::new(&system, &source, MaterialFluxNumerics::default())
        .unwrap()
        .bdf1_navier_stokes_mixed_rheology_step(
            &force,
            &boundaries,
            &initial,
            0.03,
            controls(),
            resolve,
        )
        .unwrap();
    let mut owner = RigidAleFlow::new(
        &system,
        &source,
        MaterialFluxNumerics::default(),
        [0.5, 0.5],
        RigidAleMotion::ZERO,
        &initial,
    )
    .unwrap();
    let moving = owner
        .step_with_rheology(
            0.03,
            RigidAleMotion::ZERO,
            tolerance(),
            &force,
            &boundaries,
            controls(),
            resolve,
            || false,
        )
        .unwrap();
    for (&a, &b) in original
        .step
        .solution
        .velocity_dofs()
        .iter()
        .zip(moving.flow.solution.velocity_dofs())
    {
        near(a, b, 4e-9);
    }
    for (&a, &b) in original
        .step
        .solution
        .pressure_pa()
        .iter()
        .zip(moving.flow.solution.pressure_pa())
    {
        near(a, b, 3e-8);
    }
    for (a, b) in original.samples().iter().zip(moving.rheology_samples()) {
        assert_eq!(a.point().cell, b.point().cell);
        vector(a.point().position_m, b.point().position_m, 1e-14);
        near(a.point().shear_rate_per_s, b.point().shear_rate_per_s, 2e-8);
        near(
            a.transport().dynamic_viscosity_pa_s(),
            b.transport().dynamic_viscosity_pa_s(),
            2e-8,
        );
    }
}

#[test]
fn g4_source_refusal_dropped_rheology_and_terminal_cancel_preserve_joint_state() {
    let mesh = mesh(1);
    let system = FluxSystem::new(&mesh);
    let source = material(2.0);
    let initial = history(&mesh, |_| [0.1, 0.03]);
    let incoming = |_| [0.1, 0.03];
    let boundaries = open(&mesh, &no_stress, Some(&incoming));
    let mut owner = RigidAleFlow::new(
        &system,
        &source,
        MaterialFluxNumerics::default(),
        [0.4, 0.3],
        RigidAleMotion::ZERO,
        &initial,
    )
    .unwrap();
    let before = owner.state().clone();
    let motion = RigidAleMotion {
        linear_velocity_world_m_s: [0.2, -0.1],
        angular_velocity_rad_s: 0.4,
    };
    let expected_error = FluidTransportState::try_new(-1.0, 1000.0, 0.5, 0.4).unwrap_err();
    assert!(matches!(
        owner.candidate_with_rheology(
            0.05, motion, tolerance(), &zero, &boundaries, controls(), 0.0,
            |_, _| FluidTransportState::try_new(-1.0, 1000.0, 0.5, 0.4),
            || false,
        ),
        Err(RigidAleError::Flow(MaterialFluxError::Material(error))) if error == expected_error
    ));
    assert_eq!(owner.state(), &before);
    assert!(matches!(
        owner.candidate_with_rheology(
            0.05,
            motion,
            tolerance(),
            &zero,
            &boundaries,
            controls(),
            0.0,
            |_, _| FluidTransportState::try_new(3.0, 1000.0, 0.5, 0.4),
            || false,
        ),
        Err(RigidAleError::Flow(MaterialFluxError::Invalid(
            "rheological flow must retain its constant reference density",
        )))
    ));
    assert_eq!(owner.state(), &before);
    let resolve = |point: MaterialFluxRheologyPoint, reference: &FluidTransportState| {
        FluidTransportState::try_new(
            reference.density_kg_m3(),
            1000.0,
            0.5,
            0.4 / (1.0 + point.shear_rate_per_s.powi(2)).sqrt(),
        )
    };
    let mut polls = 0;
    let expected;
    {
        let candidate = owner
            .candidate_with_rheology(
                0.05,
                motion,
                tolerance(),
                &zero,
                &boundaries,
                controls(),
                0.0,
                resolve,
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        expected = candidate.state().clone();
        assert_eq!(candidate.step().rheology_samples().len(), mesh.tris.len());
    }
    assert_eq!(
        owner.state(),
        &before,
        "dropping a rheology candidate does not commit"
    );
    let mut count = 0;
    assert!(matches!(
        owner.candidate_with_rheology(
            0.05,
            motion,
            tolerance(),
            &zero,
            &boundaries,
            controls(),
            0.0,
            resolve,
            || {
                count += 1;
                count == polls
            },
        ),
        Err(RigidAleError::Cancelled)
    ));
    assert_eq!(count, polls);
    assert_eq!(owner.state(), &before);
    let accepted = owner
        .step_with_rheology(
            0.05,
            motion,
            tolerance(),
            &zero,
            &boundaries,
            controls(),
            resolve,
            || false,
        )
        .unwrap();
    assert_eq!(owner.state(), &expected);
    assert_eq!(accepted.rheology_samples().len(), mesh.tris.len());
}
