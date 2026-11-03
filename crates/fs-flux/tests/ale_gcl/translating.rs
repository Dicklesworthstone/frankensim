//! G1/G3/G4: pressure, Galilean transport and transactional translated flow.

use fs_flux::ale::translating::{TranslatingAleError, TranslatingAleFlow};
use fs_flux::bdm::edge_gauss_pub;
use fs_flux::{
    FluxParams, FluxSystem, MaterialFluxNumerics, MaterialNavierStokesBoundary,
    NavierStokesBoundary, NavierStokesControls, OpenBoundaryPolicy,
    StokesPressureReference, TriMesh,
};
use fs_material::fluid_transport::FluidTransportState;
use fs_solid::Mesh2;

fn mesh() -> TriMesh {
    TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, 3, 3))
}

fn source() -> FluidTransportState {
    FluidTransportState::try_new(2.0, 1000.0, 0.5, 0.4).unwrap()
}

fn zero(_: [f64; 2]) -> [f64; 2] {
    [0.0; 2]
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance, "{actual} versus {expected}");
}

fn controls() -> NavierStokesControls {
    NavierStokesControls {
        relative_tolerance: 1e-10,
        ..NavierStokesControls::default()
    }
}

// Independently integrate both native edge moments instead of reusing the
// production mean-moment boost conversion.
fn history(mesh: &TriMesh, field: impl Fn([f64; 2]) -> [f64; 2]) -> Vec<f64> {
    let mut values = vec![0.0; 2 * mesh.edges.len()];
    for (edge, geometry) in mesh.edges.iter().enumerate() {
        let a = mesh.verts[geometry.verts.0];
        let b = mesh.verts[geometry.verts.1];
        for (point, weight) in edge_gauss_pub(a, b) {
            let u = field(point);
            let normal = u[0] * geometry.normal[0] + u[1] * geometry.normal[1];
            let s = ((point[0] - a[0]) * (b[0] - a[0])
                + (point[1] - a[1]) * (b[1] - a[1])) / geometry.len.powi(2) - 0.5;
            values[2 * edge] += weight * normal / geometry.len;
            values[2 * edge + 1] += weight * normal * s / geometry.len;
        }
    }
    values
}

#[test]
fn accelerating_closed_container_has_physical_pressure_gradient_and_reaction() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let material = source();
    let mut owner = TranslatingAleFlow::new(
        &system, &material, MaterialFluxNumerics::default(), [0.0; 2],
        &vec![0.0; system.n_u],
    ).unwrap();
    let speed = [0.04, -0.03];
    let wall = |_: [f64; 2]| speed;
    let edges: Vec<_> = (0..mesh.edges.len()).filter(|&e| mesh.is_boundary(e)).collect();
    let boundaries: Vec<_> = edges.iter().map(|&edge| {
        MaterialNavierStokesBoundary::Velocity { edge, value: &wall }
    }).collect();
    let dt = 0.025;
    let step = owner.step(dt, speed, 1e-12, &zero, &boundaries, controls(), || false).unwrap();
    assert_eq!(step.relative_flow.solution.pressure_reference(), StokesPressureReference::MeanZero);
    let acceleration = speed.map(|v| v / dt);
    let (relative_error, divergence) = system.velocity_error(
        &step.relative_flow.solution.discrete().x, &zero,
    );
    assert!(relative_error < 2e-9 && divergence < 2e-9);
    let (world_error, _) = system.velocity_error(&step.world_velocity_dofs, &wall);
    assert!(world_error < 2e-9);
    for (pressure, point) in step.relative_flow.solution.pressure_pa().iter().zip(&mesh.centroids) {
        let expected = -2.0 * (acceleration[0] * (point[0] - 0.5)
            + acceleration[1] * (point[1] - 0.5));
        near(*pressure, expected, 3e-8);
    }
    let origin = [step.translation_m[0] + 0.2, step.translation_m[1] + 0.1];
    let load = step.boundary_loads_world(&edges, origin, &[Some(100.0)]).unwrap();
    let force = acceleration.map(|a| -2.0 * a);
    for axis in 0..2 {
        near(load.stress.force_n_per_m[axis], force[axis], 3e-8);
        near(step.translation_m[axis], dt * speed[axis], 1e-14);
        near(step.world_momentum_n_s_per_m[axis], 2.0 * speed[axis], 2e-10);
    }
    near(load.stress.moment_n, 0.3 * force[1] - 0.4 * force[0], 3e-8);
    near(step.world_kinetic_energy_j_per_m, speed[0].powi(2) + speed[1].powi(2), 2e-10);
    near(step.mass_kg_per_m, 2.0, 2e-13);
    assert_eq!(load.origin_m, origin);
    assert!(step.boundary_loads_world(&edges, origin, &[None]).is_err());
    assert!(std::ptr::eq(owner.material(), &material));
}

#[test]
fn uniform_world_flow_survives_mesh_acceleration_and_uses_relative_crossing_flux() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let material = source();
    let world = [0.12, -0.07];
    let initial_mesh = [0.03, 0.01];
    let field = |_: [f64; 2]| world;
    let mut owner = TranslatingAleFlow::new(
        &system, &material, MaterialFluxNumerics::default(), initial_mesh,
        &history(&mesh, field),
    ).unwrap();
    let boundaries: Vec<_> = (0..mesh.edges.len()).filter(|&e| mesh.is_boundary(e))
        .map(|edge| MaterialNavierStokesBoundary::Velocity { edge, value: &field }).collect();
    let right: Vec<_> = mesh.edges.iter().enumerate()
        .filter(|(e, g)| mesh.is_boundary(*e) && g.normal[0] > 0.5)
        .map(|(e, _)| e).collect();
    for speed in [[0.05, -0.01], [-0.02, 0.03], [0.01, -0.04]] {
        let step = owner.step(0.01, speed, 1e-12, &zero, &boundaries, controls(), || false).unwrap();
        let (error, divergence) = system.velocity_error(&step.world_velocity_dofs, &field);
        assert!(error < 3e-9 && divergence < 3e-9);
        for pressure in step.relative_flow.solution.pressure_pa() {
            near(*pressure, 0.0, 3e-8);
        }
        for axis in 0..2 {
            near(step.world_momentum_n_s_per_m[axis], 2.0 * world[axis], 3e-9);
        }
        near(step.world_kinetic_energy_j_per_m, world[0].powi(2) + world[1].powi(2), 2e-9);
        let load = step.boundary_loads_world(&right, [0.0; 2], &[Some(0.0)]).unwrap();
        let crossing = world[0] - step.mesh_velocity_m_s[0];
        let force = world.map(|v| 2.0 * crossing * v);
        for axis in 0..2 {
            near(load.outward_momentum.force_n_per_m[axis], force[axis], 3e-9);
        }
        let center = [1.0 + step.translation_m[0], 0.5 + step.translation_m[1]];
        near(load.outward_momentum.moment_n, center[0] * force[1] - center[1] * force[0], 3e-9);
    }
}

#[test]
fn nonlinear_mixed_flow_and_world_sampling_are_galilean_covariant() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let material = source();
    let dt = 0.02;
    let relative = |p: [f64; 2]| [0.2 * p[0] + 0.3 * p[1] + 0.4, -0.2 * p[1] + 0.1];
    let mut solutions = Vec::new();
    for boost in [[0.0; 2], [0.7, -0.4]] {
        let translation = boost.map(|v| dt * v);
        let world_velocity = |x: [f64; 2]| {
            let v = relative([x[0] - translation[0], x[1] - translation[1]]);
            [v[0] + boost[0], v[1] + boost[1]]
        };
        let force = |x: [f64; 2]| {
            let y = [x[0] - translation[0], x[1] - translation[1]];
            [2.0 * (0.04 * y[0] + 0.11) + 2.0, 2.0 * (0.04 * y[1] - 0.02) - 3.0]
        };
        let traction = |x: [f64; 2], n: [f64; 2]| {
            let y = [x[0] - translation[0], x[1] - translation[1]];
            let pressure = 7.0 + 2.0 * y[0] - 3.0 * y[1];
            [(0.4 * 0.4 - pressure) * n[0] + 0.3 * 0.4 * n[1],
             0.3 * 0.4 * n[0] + (-0.4 * 0.4 - pressure) * n[1]]
        };
        let initial = history(&mesh, |p| {
            let v = relative(p);
            [v[0] + boost[0], v[1] + boost[1]]
        });
        let mut owner = TranslatingAleFlow::new(
            &system, &material, MaterialFluxNumerics::default(), boost, &initial,
        ).unwrap();
        let boundaries: Vec<_> = mesh.edges.iter().enumerate()
            .filter(|(e, _)| mesh.is_boundary(*e)).map(|(edge, geometry)| {
                if geometry.normal[0].abs() > 0.5 {
                    MaterialNavierStokesBoundary::Open {
                        edge, traction: &traction, backflow: OpenBoundaryPolicy::NaturalTraction,
                    }
                } else {
                    MaterialNavierStokesBoundary::Velocity { edge, value: &world_velocity }
                }
            }).collect();
        let step = owner.step(dt, boost, 1e-12, &force, &boundaries, controls(), || false).unwrap();
        let (error, divergence) = system.velocity_error(&step.relative_flow.solution.discrete().x, &relative);
        assert!(error < 4e-9 && divergence < 4e-9);
        assert_eq!(step.relative_flow.solution.pressure_reference(), StokesPressureReference::TractionDatum);
        for (actual, p) in step.relative_flow.solution.pressure_pa().iter().zip(&mesh.centroids) {
            near(*actual, 7.0 + 2.0 * p[0] - 3.0 * p[1], 4e-8);
        }
        assert!(step.relative_flow.nonlinear_residual < 1e-10);
        solutions.push(step.relative_flow.solution.velocity_dofs().to_vec());
    }
    for (stationary, moving) in solutions[0].iter().zip(&solutions[1]) {
        near(*stationary, *moving, 3e-9);
    }
}

#[test]
fn mass_flux_receipt_retains_the_original_advector_instead_of_the_endpoint_trace() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let advector = history(&mesh, |_| [0.1, 0.05]);
    let previous = history(&mesh, |_| [0.4, -0.2]);
    let traction = |_: [f64; 2], n: [f64; 2]| n.map(|v| -7.0 * v);
    let boundaries: Vec<_> = (0..mesh.edges.len()).filter(|&e| mesh.is_boundary(e))
        .map(|edge| NavierStokesBoundary::Open {
            edge, traction: &traction, backflow: OpenBoundaryPolicy::NaturalTraction,
        }).collect();
    let step = system.assemble_navier_stokes_mixed_bdf1(
        FluxParams { nu: 0.2, ..FluxParams::default() }, &zero, &boundaries,
        &advector, &previous, 0.01,
    ).unwrap().solve().unwrap();
    let right: Vec<_> = mesh.edges.iter().enumerate()
        .filter(|(e, g)| mesh.is_boundary(*e) && g.normal[0] > 0.5)
        .map(|(e, _)| e).collect();
    let (flux, first) = step.reactions.volume_flux_moments(&right, [0.2, 0.3]).unwrap();
    near(flux, 0.1, 2e-12);
    near(first[0], 0.08, 2e-12);
    near(first[1], 0.02, 2e-12);
    let (error, _) = system.velocity_error(&step.field.x, &|_| [0.4, -0.2]);
    assert!(error < 2e-9);
    assert!((flux - 0.4).abs() > 0.2);
    assert!(step.reactions.volume_flux_moments(&[right[0], right[0]], [0.0; 2]).is_err());
}

#[test]
fn final_cancellation_and_unsupported_frame_laws_preserve_continuation() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let material = source();
    let make = || TranslatingAleFlow::new(
        &system, &material, MaterialFluxNumerics::default(), [0.0; 2],
        &vec![0.0; system.n_u],
    ).unwrap();
    let mut owner = make();
    let mut reference = make();
    let speed = [0.02, -0.01];
    let wall = |_: [f64; 2]| speed;
    let edges: Vec<_> = (0..mesh.edges.len()).filter(|&e| mesh.is_boundary(e)).collect();
    let boundaries: Vec<_> = edges.iter().map(|&edge| {
        MaterialNavierStokesBoundary::Velocity { edge, value: &wall }
    }).collect();
    for flow in [&mut owner, &mut reference] {
        flow.step(0.01, speed, 1e-12, &zero, &boundaries, controls(), || false).unwrap();
    }
    let before = owner.state().clone();
    let mut polls = 0;
    let expected = reference.step(0.01, speed, 1e-12, &zero, &boundaries, controls(), || {
        polls += 1;
        false
    }).unwrap();
    let mut count = 0;
    let failed = owner.step(0.01, speed, 1e-12, &zero, &boundaries, controls(), || {
        count += 1;
        count == polls
    });
    assert!(matches!(failed, Err(TranslatingAleError::Cancelled)));
    assert_eq!(owner.state(), &before);
    let actual = owner.step(0.01, speed, 1e-12, &zero, &boundaries, controls(), || false).unwrap();
    assert_eq!(owner.state(), reference.state());
    assert_eq!(actual.world_velocity_dofs, expected.world_velocity_dofs);
    let before = owner.state().clone();
    let traction = |_: [f64; 2], _: [f64; 2]| [0.0; 2];
    let unsupported = [MaterialNavierStokesBoundary::Open {
        edge: edges[0], traction: &traction, backflow: OpenBoundaryPolicy::DissipativeBackflow,
    }];
    assert!(matches!(
        owner.step(0.01, speed, 1e-12, &zero, &unsupported, controls(), || false),
        Err(TranslatingAleError::Invalid(_))
    ));
    assert!(owner.step(0.01, [f64::MAX, 0.0], 1e-12, &zero, &boundaries, controls(), || false).is_err());
    assert_eq!(owner.state(), &before);
}

#[test]
fn reservoir_inflow_and_world_momentum_are_covariant_on_translating_domains() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let material = source();
    let density = 2.0;
    let dt = 0.02;
    let relative = |p: [f64; 2]| [0.2 + 0.1 * p[1], 0.03];
    let left: Vec<_> = mesh
        .edges
        .iter()
        .enumerate()
        .filter(|(edge, e)| mesh.is_boundary(*edge) && e.normal[0] < -0.5)
        .map(|(edge, _)| edge)
        .collect();
    let mut relative_solutions = Vec::new();
    for boost in [[0.0; 2], [0.7, -0.4]] {
        let translation = boost.map(|v| dt * v);
        let queries = std::cell::Cell::new(0);
        let reservoir = |world: [f64; 2]| {
            queries.set(queries.get() + 1);
            let u = relative([world[0] - translation[0], world[1] - translation[1]]);
            [u[0] + boost[0], u[1] + boost[1]]
        };
        // grad(u)*u = [.003,0], and div(2 mu D)=0. The original material
        // density multiplies inertia/advection; the physical pressure is 7+2x-3y.
        let force = |_: [f64; 2]| [density * 0.003 + 2.0, -3.0];
        let traction = |world: [f64; 2], n: [f64; 2]| {
            let p = [world[0] - translation[0], world[1] - translation[1]];
            let pressure = 7.0 + 2.0 * p[0] - 3.0 * p[1];
            [
                -pressure * n[0] + 0.04 * n[1],
                0.04 * n[0] - pressure * n[1],
            ]
        };
        let initial = history(&mesh, |point| {
            let u = relative(point);
            [u[0] + boost[0], u[1] + boost[1]]
        });
        let mut owner = TranslatingAleFlow::new(
            &system,
            &material,
            MaterialFluxNumerics::default(),
            boost,
            &initial,
        )
        .unwrap();
        let boundaries: Vec<_> = (0..mesh.edges.len())
            .filter(|&edge| mesh.is_boundary(edge))
            .map(|edge| MaterialNavierStokesBoundary::OpenReservoir {
                edge,
                traction: &traction,
                velocity: &reservoir,
            })
            .collect();
        let step = owner
            .step(dt, boost, 1e-12, &force, &boundaries, controls(), || false)
            .unwrap();
        assert_eq!(queries.get(), 2 * boundaries.len());
        let (error, divergence) =
            system.velocity_error(&step.relative_flow.solution.discrete().x, &relative);
        assert!(error < 4e-9 && divergence < 4e-9, "{error}, {divergence}");
        let load = step
            .boundary_loads_world(&left, translation, &[None])
            .unwrap();
        // At x=0, q=.2+.1y crosses inward with speed q. World incoming
        // velocity is (q+boost_x,.03+boost_y); integrate the physical flux.
        near(
            load.outward_momentum.force_n_per_m[0],
            -density * (19.0 / 300.0 + 0.25 * boost[0]),
            4e-8,
        );
        near(
            load.outward_momentum.force_n_per_m[1],
            -density * 0.25 * (0.03 + boost[1]),
            4e-8,
        );
        near(
            load.outward_momentum.moment_n,
            density * (43.0 / 1200.0 + 2.0 * boost[0] / 15.0),
            4e-8,
        );
        near(load.stress.force_n_per_m[0], -5.5, 4e-8);
        near(load.stress.force_n_per_m[1], 0.04, 4e-8);
        assert_eq!(load.backflow.force_n_per_m, [0.0; 2]);
        assert_eq!(queries.get(), 2 * boundaries.len());
        relative_solutions.push(step.relative_flow.solution.velocity_dofs().to_vec());
        assert!(std::ptr::eq(owner.material(), &material));
    }
    for (a, b) in relative_solutions[0].iter().zip(&relative_solutions[1]) {
        near(*a, *b, 4e-9);
    }
}

#[test]
fn authored_normal_velocity_correction_admits_pressure_balance_and_retains_small_motion() {
    // A pressure-supported gradient force with zero or tiny uniform velocity
    // is an exact nonlinear endpoint. The absolute criterion acts on the
    // unrelaxed physical normal traces; it never erases that small velocity.
    for subdivisions in [2, 4] {
        let mesh = TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, subdivisions, subdivisions));
        let system = FluxSystem::new(&mesh);
        let force = |point: [f64; 2]| [0.2 + 0.3 * point[0], 0.0];
        let previous = vec![0.0; system.n_u];
        let params = FluxParams {
            nu: 0.2,
            ..FluxParams::default()
        };
        let limits = NavierStokesControls {
            maximum_iterations: 1,
            ..controls()
        };
        for velocity in [[0.0; 2], [2e-8, -1e-8]] {
            let field = |_: [f64; 2]| velocity;
            let boundaries: Vec<_> = (0..mesh.edges.len())
                .filter(|&e| mesh.is_boundary(e))
                .map(|edge| NavierStokesBoundary::Velocity {
                    edge,
                    value: &field,
                })
                .collect();
            let endpoint = system
                .bdf1_navier_stokes_mixed_step_with_velocity_tolerance(
                    params,
                    &force,
                    &boundaries,
                    &previous,
                    0.025,
                    limits,
                    1e-7,
                )
                .unwrap();
            let (error, divergence) = system.velocity_error(&endpoint.solution.field.x, &field);
            assert!(error < 1e-12 && divergence < 1e-12);
            assert_eq!(endpoint.iterations, 1);
            assert!(endpoint.nonlinear_residual <= limits.relative_tolerance);
            if velocity[0] != 0.0 {
                // With a nonzero first correction from a zero history, the
                // legacy relative-only test is one. This is not a dependence
                // on platform-specific noise in an exactly quiescent field.
                let legacy = system.bdf1_navier_stokes_mixed_step(
                    params,
                    &force,
                    &boundaries,
                    &previous,
                    0.025,
                    limits,
                );
                let zero_absolute = system.bdf1_navier_stokes_mixed_step_with_velocity_tolerance(
                    params,
                    &force,
                    &boundaries,
                    &previous,
                    0.025,
                    limits,
                    0.0,
                );
                assert!(legacy.is_err() && zero_absolute.is_err());
                assert_eq!(
                    format!("{:?}", legacy.err()),
                    format!("{:?}", zero_absolute.err())
                );
            }
        }
    }
}

#[test]
fn absolute_velocity_acceptance_keeps_nonlinear_residuals_mandatory() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let moving_lid = |_: [f64; 2]| [0.4, 0.0];
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| NavierStokesBoundary::Velocity {
            edge,
            value: if mesh.edges[edge].normal[1] > 0.5 {
                &moving_lid
            } else {
                &zero
            },
        })
        .collect();
    let result = system.bdf1_navier_stokes_mixed_step_with_velocity_tolerance(
        FluxParams {
            nu: 0.2,
            ..FluxParams::default()
        },
        &zero,
        &boundaries,
        &vec![0.0; system.n_u],
        0.2,
        NavierStokesControls {
            maximum_iterations: 1,
            ..controls()
        },
        10.0,
    );
    // The generous correction allowance cannot turn a frozen Stokes/Oseen
    // answer into an accepted nonlinear flow; the endpoint convection differs.
    assert!(matches!(result,
        Err(fs_flux::StokesBoundaryError::Rheology(
            fs_flux::FluxRheologyError::NonlinearConvergence { residual, .. }
        )) if residual > controls().relative_tolerance));
}

#[test]
fn discarded_coupling_candidates_keep_native_history_and_commit_reuses_the_solve() {
    let mesh = mesh();
    let system = FluxSystem::new(&mesh);
    let material = source();
    let make = || {
        TranslatingAleFlow::new(
            &system,
            &material,
            MaterialFluxNumerics::default(),
            [0.0; 2],
            &vec![0.0; system.n_u],
        )
        .unwrap()
    };
    let mut owner = make();
    let mut reference = make();
    let speed = [0.12, -0.03];
    let wall = |_: [f64; 2]| speed;
    let boundaries: Vec<_> = (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .map(|edge| MaterialNavierStokesBoundary::Velocity { edge, value: &wall })
        .collect();
    let before = owner.state().clone();
    let candidate = owner
        .candidate(
            0.02,
            speed,
            1e-12,
            &zero,
            &boundaries,
            controls(),
            0.0,
            || false,
        )
        .unwrap();
    assert!(candidate.state().time_s() > before.time_s());
    assert!(candidate.step().world_momentum_n_s_per_m[0] > 0.0);
    drop(candidate);
    assert_eq!(owner.state(), &before);
    let expected = reference
        .step(0.02, speed, 1e-12, &zero, &boundaries, controls(), || false)
        .unwrap();
    let candidate = owner
        .candidate(
            0.02,
            speed,
            1e-12,
            &zero,
            &boundaries,
            controls(),
            0.0,
            || false,
        )
        .unwrap();
    let proposed = candidate.state().clone();
    let actual = candidate.commit();
    assert_eq!(owner.state(), &proposed);
    assert_eq!(owner.state(), reference.state());
    assert_eq!(actual.world_velocity_dofs, expected.world_velocity_dofs);
}
