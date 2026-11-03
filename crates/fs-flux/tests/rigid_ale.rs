//! G0/G1/G3/G4 native rigid-ALE physics. Numerical execution is required:
//! source/formatter checks alone do not establish these conformance claims.

use fs_flux::ale::rigid::{RigidAleError, RigidAleFlow, RigidAleMotion, RigidAleMotionTolerance};
use fs_flux::ale::translating::TranslatingAleFlow;
use fs_flux::{
    FluxSystem, MaterialFluxNumerics, MaterialNavierStokesBoundary, NavierStokesControls,
    OpenBoundaryPolicy, TriMesh,
};
use fs_material::fluid_transport::FluidTransportState;
use fs_solid::Mesh2;
use std::cell::Cell;

fn mesh(n: usize) -> TriMesh {
    TriMesh::from_mesh2(&Mesh2::triangles(1.0, 1.0, n, n))
}
fn material(rho: f64) -> FluidTransportState {
    FluidTransportState::try_new(rho, 1000.0, 0.5, 0.4).unwrap()
}
fn controls() -> NavierStokesControls {
    NavierStokesControls {
        relative_tolerance: 1e-10,
        ..NavierStokesControls::default()
    }
}
fn tolerance() -> RigidAleMotionTolerance {
    RigidAleMotionTolerance {
        linear_velocity_m_s: 1e-12,
        angular_velocity_rad_s: 1e-12,
    }
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a:.16e} != {b:.16e}");
}
fn vector(a: [f64; 2], b: [f64; 2], tol: f64) {
    near(a[0], b[0], tol);
    near(a[1], b[1], tol);
}
fn add(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] + b[0], a[1] + b[1]]
}
fn sub(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}
fn scale(a: [f64; 2], s: f64) -> [f64; 2] {
    a.map(|v| v * s)
}
fn spin(a: [f64; 2]) -> [f64; 2] {
    [-a[1], a[0]]
}
fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[1] - a[1] * b[0]
}
fn rotate([c, s]: [f64; 2], a: [f64; 2]) -> [f64; 2] {
    [c * a[0] - s * a[1], s * a[0] + c * a[1]]
}
fn inverse([c, s]: [f64; 2], a: [f64; 2]) -> [f64; 2] {
    rotate([c, -s], a)
}
fn zero(_: [f64; 2]) -> [f64; 2] {
    [0.0; 2]
}
fn no_stress(_: [f64; 2], _: [f64; 2]) -> [f64; 2] {
    [0.0; 2]
}
fn edges(mesh: &TriMesh) -> Vec<usize> {
    (0..mesh.edges.len())
        .filter(|&e| mesh.is_boundary(e))
        .collect()
}
fn history(mesh: &TriMesh, field: impl Fn([f64; 2]) -> [f64; 2]) -> Vec<f64> {
    // Independent exact integrals of affine normal velocity along each edge.
    let mut u = vec![0.0; 2 * mesh.edges.len()];
    for (i, e) in mesh.edges.iter().enumerate() {
        let a = field(mesh.verts[e.verts.0]);
        let b = field(mesh.verts[e.verts.1]);
        for axis in 0..2 {
            u[2 * i] += 0.5 * (a[axis] + b[axis]) * e.normal[axis];
            u[2 * i + 1] += (b[axis] - a[axis]) * e.normal[axis] / 12.0;
        }
    }
    u
}
fn open<'a>(
    mesh: &TriMesh,
    traction: &'a dyn Fn([f64; 2], [f64; 2]) -> [f64; 2],
    incoming: Option<&'a dyn Fn([f64; 2]) -> [f64; 2]>,
) -> Vec<MaterialNavierStokesBoundary<'a>> {
    edges(mesh)
        .into_iter()
        .map(|edge| match incoming {
            Some(velocity) => MaterialNavierStokesBoundary::OpenReservoir {
                edge,
                traction,
                velocity,
            },
            None => MaterialNavierStokesBoundary::Open {
                edge,
                traction,
                backflow: OpenBoundaryPolicy::NaturalTraction,
            },
        })
        .collect()
}

#[test]
fn g1_g3_free_world_stream_survives_changing_rigid_frames_and_transports_real_momentum() {
    let mesh = mesh(2);
    let system = FluxSystem::new(&mesh);
    let source = material(2.3);
    let rho = 2.3;
    let u = [0.4, -0.2];
    let pivot_ref = [0.2, -0.3];
    let initial = history(&mesh, |_| u);
    let reservoir = |_| u;
    let boundaries = open(&mesh, &no_stress, Some(&reservoir));
    let all = edges(&mesh);
    let mut owner = RigidAleFlow::new(
        &system,
        &source,
        MaterialFluxNumerics::default(),
        pivot_ref,
        RigidAleMotion::ZERO,
        &initial,
    )
    .unwrap();
    for (v, omega) in [([0.3, -0.2], 1.2), ([-0.1, 0.4], -0.7), ([0.2, 0.1], 0.0)] {
        let requested = RigidAleMotion {
            linear_velocity_world_m_s: v,
            angular_velocity_rad_s: omega,
        };
        let preview = owner.kinematics(0.08, requested, tolerance()).unwrap();
        let step = owner
            .step(
                0.08,
                requested,
                tolerance(),
                &zero,
                &boundaries,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(step.pivot_world_m, preview.pivot_world_m);
        assert_eq!(step.orientation_cos_sin, preview.orientation_cos_sin);
        assert_eq!(step.mesh_velocity, preview.mesh_velocity);
        let local = inverse(step.orientation_cos_sin, u);
        let (error, div) = system.velocity_error(step.flow.solution.velocity_dofs(), &|_| local);
        assert!(error < 3e-9 && div < 3e-9, "{error} {div}");
        near(step.mass_kg_per_m, rho, 2e-12);
        vector(step.world_momentum_n_s_per_m, scale(u, rho), 3e-9);
        near(
            step.world_kinetic_energy_j_per_m,
            0.5 * rho * (u[0] * u[0] + u[1] * u[1]),
            3e-9,
        );
        let origin = [-0.4, 0.6];
        let center = add(
            step.pivot_world_m,
            rotate(step.orientation_cos_sin, sub([0.5, 0.5], pivot_ref)),
        );
        near(
            step.world_angular_momentum(origin).unwrap(),
            cross(sub(center, origin), scale(u, rho)),
            3e-9,
        );
        vector(step.point_world([0.5, 0.5]).unwrap(), center, 1e-14);
        let transport = step.boundary_transport_world(&all, origin).unwrap();
        near(transport.outward_mass_flux_kg_per_m_s, 0.0, 3e-9);
        vector(transport.outward_momentum.force_n_per_m, [0.0; 2], 3e-9);
        let mean_mesh = add(
            step.mesh_velocity.linear_velocity_world_m_s,
            scale(
                spin(sub(center, step.pivot_world_m)),
                step.mesh_velocity.angular_velocity_rad_s,
            ),
        );
        near(
            transport.outward_momentum.moment_n,
            -rho * cross(mean_mesh, u),
            3e-9,
        );
        near(transport.outward_kinetic_energy_w_per_m, 0.0, 3e-9);
        near(transport.upwind_kinetic_dissipation_w_per_m, 0.0, 3e-9);
        vector(
            step.boundary_loads_world(&all, origin, &[None])
                .unwrap()
                .stress
                .force_n_per_m,
            [0.0; 2],
            3e-9,
        );
    }
}

#[test]
fn g0_cross_mass_history_keeps_world_momentum_with_discontinuous_tangential_traces() {
    let mesh = mesh(1);
    let system = FluxSystem::new(&mesh);
    let source = material(2.0);
    assert_eq!(mesh.tris.len(), 2);
    let diagonal = mesh.edges.iter().find(|e| e.tris.1 != usize::MAX).unwrap();
    let tangent = spin(diagonal.normal);
    let fields = [scale(tangent, 0.2), scale(tangent, -0.1)];
    let mut initial = vec![0.0; system.n_u];
    for (i, e) in mesh.edges.iter().enumerate() {
        // The common interface normal is zero on both sides, while tangential
        // traces differ. A pointwise rotated old field need not itself belong
        // to the current global H(div) space; only its mass load is transferred.
        if e.tris.1 == usize::MAX {
            let v = fields[e.tris.0];
            initial[2 * i] = v[0] * e.normal[0] + v[1] * e.normal[1];
        }
    }
    let before = scale(
        add(
            scale(fields[0], mesh.areas[0]),
            scale(fields[1], mesh.areas[1]),
        ),
        2.0,
    );
    let boundaries = open(&mesh, &no_stress, None);
    let mut owner = RigidAleFlow::new(
        &system,
        &source,
        MaterialFluxNumerics::default(),
        [0.5, 0.5],
        RigidAleMotion::ZERO,
        &initial,
    )
    .unwrap();
    let dt = 0.02;
    let step = owner
        .step(
            dt,
            RigidAleMotion {
                linear_velocity_world_m_s: [0.1, -0.04],
                angular_velocity_rad_s: 7.0,
            },
            tolerance(),
            &zero,
            &boundaries,
            controls(),
            || false,
        )
        .unwrap();
    let flux = step
        .boundary_transport_world(&edges(&mesh), [0.0; 2])
        .unwrap();
    vector(
        add(
            sub(step.world_momentum_n_s_per_m, before),
            scale(flux.outward_momentum.force_n_per_m, dt),
        ),
        [0.0; 2],
        3e-9,
    );
    vector(
        step.boundary_loads_world(&edges(&mesh), [0.0; 2], &[None])
            .unwrap()
            .stress
            .force_n_per_m,
        [0.0; 2],
        3e-9,
    );
}

#[test]
fn g3_zero_rotation_matches_original_translating_nonlinear_owner() {
    let mesh = mesh(2);
    let system = FluxSystem::new(&mesh);
    let source = material(2.0);
    let exact = |x: [f64; 2]| [0.15 * x[1], 0.08];
    let initial = history(&mesh, exact);
    let boundaries: Vec<_> = edges(&mesh)
        .into_iter()
        .map(|edge| MaterialNavierStokesBoundary::Velocity {
            edge,
            value: &exact,
        })
        .collect();
    let mut translation = TranslatingAleFlow::new(
        &system,
        &source,
        MaterialFluxNumerics::default(),
        [0.0; 2],
        &initial,
    )
    .unwrap();
    let mut rigid = RigidAleFlow::new(
        &system,
        &source,
        MaterialFluxNumerics::default(),
        [0.5, 0.5],
        RigidAleMotion::ZERO,
        &initial,
    )
    .unwrap();
    let force = |_| [2.0 * 0.15 * 0.08, 0.0];
    for velocity in [[0.3, -0.2], [-0.4, 0.1]] {
        let a = translation
            .step(
                0.03,
                velocity,
                1e-12,
                &force,
                &boundaries,
                controls(),
                || false,
            )
            .unwrap();
        let b = rigid
            .step(
                0.03,
                RigidAleMotion {
                    linear_velocity_world_m_s: velocity,
                    angular_velocity_rad_s: 0.0,
                },
                tolerance(),
                &force,
                &boundaries,
                controls(),
                || false,
            )
            .unwrap();
        for (&x, &y) in a
            .world_velocity_dofs
            .iter()
            .zip(b.flow.solution.velocity_dofs())
        {
            near(x, y, 3e-9);
        }
        vector(a.world_momentum_n_s_per_m, b.world_momentum_n_s_per_m, 3e-9);
        near(
            a.world_kinetic_energy_j_per_m,
            b.world_kinetic_energy_j_per_m,
            3e-9,
        );
        for (&x, &y) in a
            .relative_flow
            .solution
            .pressure_pa()
            .iter()
            .zip(b.flow.solution.pressure_pa())
        {
            near(x, y, 3e-8);
        }
    }
}

#[test]
fn g1_rotating_couette_preserves_physical_slip_stress_torque_and_density() {
    let mesh = mesh(2);
    let system = FluxSystem::new(&mesh);
    let bottom: Vec<_> = edges(&mesh)
        .into_iter()
        .filter(|&e| mesh.edges[e].normal[1] < -0.5)
        .collect();
    let pivot_ref = [0.5, 0.5];
    let (mu, beta, u0, u1, pressure) = (0.4, 1.6, 0.1, 0.4, 2.5);
    let tau = mu * (u1 - u0) / (1.0 + 2.0 * mu / beta);
    let slope = tau / mu;
    let relative = |y: [f64; 2]| [u0 + tau / beta + slope * y[1], 0.0];
    let motion = RigidAleMotion {
        linear_velocity_world_m_s: [0.13, -0.07],
        angular_velocity_rad_s: 0.8,
    };
    let dt = 0.06;
    for rho in [2.0, 7.0] {
        let source = material(rho);
        let initial = history(&mesh, |y| {
            add(
                motion.linear_velocity_world_m_s,
                add(
                    scale(spin(sub(y, pivot_ref)), motion.angular_velocity_rad_s),
                    relative(y),
                ),
            )
        });
        let mut owner = RigidAleFlow::new(
            &system,
            &source,
            MaterialFluxNumerics::default(),
            pivot_ref,
            motion,
            &initial,
        )
        .unwrap();
        let k = owner.kinematics(dt, motion, tolerance()).unwrap();
        let r = k.orientation_cos_sin;
        let omega = k.mesh_velocity.angular_velocity_rad_s;
        let to_reference = |x| add(pivot_ref, inverse(r, sub(x, k.pivot_world_m)));
        let wall = |x: [f64; 2]| {
            let y = to_reference(x);
            let speed = if y[1] < 0.5 { u0 } else { u1 };
            add(
                k.mesh_velocity.linear_velocity_world_m_s,
                rotate(r, add(scale(spin(sub(y, pivot_ref)), omega), [speed, 0.0])),
            )
        };
        let traction = |_: [f64; 2], n: [f64; 2]| {
            let n = inverse(r, n);
            rotate(
                r,
                [-pressure * n[0] + tau * n[1], tau * n[0] - pressure * n[1]],
            )
        };
        // Independent exact backward-Euler world acceleration at paired
        // material locations, plus grad(u)*(u-w_mesh). The manufactured load
        // isolates frame/history assembly from time-truncation error.
        let force = |x: [f64; 2]| {
            let y = to_reference(x);
            let base = add(scale(spin(sub(y, pivot_ref)), omega), relative(y));
            let temporal = scale(sub(rotate(r, base), base), 1.0 / dt);
            let advective = rotate(r, scale(spin(relative(y)), omega));
            scale(add(temporal, advective), rho)
        };
        let boundaries: Vec<_> = edges(&mesh)
            .into_iter()
            .map(|edge| {
                if mesh.edges[edge].normal[1].abs() > 0.5 {
                    MaterialNavierStokesBoundary::NavierSlip {
                        edge,
                        value: &wall,
                        tangential_resistance_pa_s_per_m: beta,
                    }
                } else {
                    MaterialNavierStokesBoundary::Open {
                        edge,
                        traction: &traction,
                        backflow: OpenBoundaryPolicy::NaturalTraction,
                    }
                }
            })
            .collect();
        let step = owner
            .step(
                dt,
                motion,
                tolerance(),
                &force,
                &boundaries,
                controls(),
                || false,
            )
            .unwrap();
        let exact = |y| {
            add(
                inverse(r, k.mesh_velocity.linear_velocity_world_m_s),
                add(scale(spin(sub(y, pivot_ref)), omega), relative(y)),
            )
        };
        let (error, div) = system.velocity_error(step.flow.solution.velocity_dofs(), &exact);
        assert!(error < 3e-8 && div < 3e-8, "{error} {div}");
        let origin = step.point_world([0.2, -0.4]).unwrap();
        let load = step.boundary_loads_world(&bottom, origin, &[None]).unwrap();
        vector(load.stress.force_n_per_m, rotate(r, [tau, -pressure]), 6e-8);
        near(load.stress.moment_n, -0.3 * pressure - 0.4 * tau, 6e-8);
        near(load.slip_dissipation_w_per_m, tau * tau / beta, 3e-8);
        let flux = step.boundary_transport_world(&bottom, origin).unwrap();
        near(flux.outward_mass_flux_kg_per_m_s, 0.0, 3e-9);
        vector(flux.outward_momentum.force_n_per_m, [0.0; 2], 3e-9);
    }
}

#[test]
fn g1_angular_inventory_and_velocity_refine_for_physical_rigid_rotation() {
    let mesh = mesh(2);
    let system = FluxSystem::new(&mesh);
    let source = material(2.0);
    let pivot = [0.5, 0.5];
    let alpha = 0.2;
    let exact = |y| scale(spin(sub(y, pivot)), alpha);
    let initial = history(&mesh, exact);
    let boundaries = open(&mesh, &no_stress, None);
    // A steady physical solid rotation is sustained by its exact centripetal
    // body force, with zero strain and zero pressure. This is an independent
    // continuous solution; unlike the Couette test, no dt correction is added.
    let force = |x| scale(sub(x, pivot), -2.0 * alpha * alpha);
    let mut velocity_errors = [0.0; 3];
    let mut angular_errors = [0.0; 3];
    for (slot, count) in [4, 8, 16].into_iter().enumerate() {
        let mut owner = RigidAleFlow::new(
            &system,
            &source,
            MaterialFluxNumerics::default(),
            pivot,
            RigidAleMotion::ZERO,
            &initial,
        )
        .unwrap();
        let mut endpoint = None;
        for _ in 0..count {
            endpoint = Some(
                owner
                    .step(
                        0.2 / f64::from(count),
                        RigidAleMotion {
                            linear_velocity_world_m_s: [0.0; 2],
                            angular_velocity_rad_s: 1.2,
                        },
                        tolerance(),
                        &force,
                        &boundaries,
                        controls(),
                        || false,
                    )
                    .unwrap(),
            );
        }
        let step = endpoint.unwrap();
        let (error, div) = system.velocity_error(step.flow.solution.velocity_dofs(), &exact);
        velocity_errors[slot] = error;
        angular_errors[slot] =
            (step.world_angular_momentum(pivot).unwrap() - 2.0 * alpha / 6.0).abs();
        assert!(div < 3e-8);
    }
    for errors in [velocity_errors, angular_errors] {
        assert!(
            errors[0] / errors[1] > 1.6 && errors[1] / errors[2] > 1.6,
            "{errors:?}"
        );
        assert!(errors[2] < 2e-3, "{errors:?}");
    }
}

#[test]
fn g4_motion_preview_dropped_candidate_and_terminal_cancel_preserve_one_owner() {
    let mesh = mesh(1);
    let system = FluxSystem::new(&mesh);
    let source = material(2.0);
    let initial = history(&mesh, |_| [0.1, 0.03]);
    let calls = Cell::new(0);
    let traction = |_: [f64; 2], _: [f64; 2]| {
        calls.set(calls.get() + 1);
        [0.0; 2]
    };
    let boundaries = open(&mesh, &traction, None);
    let mut owner = RigidAleFlow::new(
        &system,
        &source,
        MaterialFluxNumerics::default(),
        [0.4, 0.3],
        RigidAleMotion::ZERO,
        &initial,
    )
    .unwrap();
    let old = owner.state().clone();
    let requested = RigidAleMotion {
        linear_velocity_world_m_s: [0.2, -0.1],
        angular_velocity_rad_s: 0.4,
    };
    for bad in [f64::NAN, core::f64::consts::PI, 4.0] {
        assert!(
            owner
                .kinematics(
                    1.0,
                    RigidAleMotion {
                        angular_velocity_rad_s: bad,
                        ..requested
                    },
                    tolerance()
                )
                .is_err()
        );
    }
    assert_eq!(calls.get(), 0);
    assert_eq!(owner.state(), &old);
    let preview = owner.kinematics(0.05, requested, tolerance()).unwrap();
    let mut polls = 0;
    let expected;
    {
        let candidate = owner
            .candidate(
                0.05,
                requested,
                tolerance(),
                &zero,
                &boundaries,
                controls(),
                0.0,
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        expected = candidate.state().clone();
        assert_eq!(candidate.step().pivot_world_m, preview.pivot_world_m);
        assert_eq!(
            candidate.step().orientation_cos_sin,
            preview.orientation_cos_sin
        );
        assert_eq!(candidate.step().mesh_velocity, preview.mesh_velocity);
        // Drop, without commit.
    }
    assert_eq!(owner.state(), &old);
    let mut observed = 0;
    assert!(matches!(
        owner.candidate(
            0.05,
            requested,
            tolerance(),
            &zero,
            &boundaries,
            controls(),
            0.0,
            || {
                observed += 1;
                observed == polls
            }
        ),
        Err(RigidAleError::Cancelled)
    ));
    assert_eq!(observed, polls);
    assert_eq!(owner.state(), &old);
    let step = owner
        .step(
            0.05,
            requested,
            tolerance(),
            &zero,
            &boundaries,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(owner.state(), &expected);
    assert_eq!(step.rotation_increment_rad, preview.rotation_increment_rad);
    assert_eq!(step.time_s, 0.05);
    let fixed = owner
        .kinematics(0.05, RigidAleMotion::ZERO, tolerance())
        .unwrap();
    assert_eq!(
        fixed.orientation_cos_sin,
        owner.state().orientation_cos_sin()
    );
    assert_eq!(fixed.mesh_velocity.angular_velocity_rad_s, 0.0);
    let unsupported: Vec<_> = edges(&mesh)
        .into_iter()
        .map(|edge| MaterialNavierStokesBoundary::Open {
            edge,
            traction: &traction,
            backflow: OpenBoundaryPolicy::DissipativeBackflow,
        })
        .collect();
    let accepted = owner.state().clone();
    let before_calls = calls.get();
    assert!(
        owner
            .step(
                0.05,
                requested,
                tolerance(),
                &zero,
                &unsupported,
                controls(),
                || false
            )
            .is_err()
    );
    assert_eq!(calls.get(), before_calls);
    assert_eq!(owner.state(), &accepted);
}
