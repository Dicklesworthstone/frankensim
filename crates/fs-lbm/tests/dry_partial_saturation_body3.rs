//! Positive dry inertia, real streamed interior feedback and late rollback.
use fs_exec::CancelGate;
use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_lbm::d3q19::{
    BodyFrameEllipsoid3, BodyFrameSampler3, BoundaryGrid3, BoundarySpec3, DryBodyCouplingConfig3,
    DryBodyExchange3, E3, FaceBoundary3, MovingPartialSaturationBody3, PartialSaturationError3, Q3,
    equilibrium3,
};
use fs_time::lie::gyrostat::InertiaTensor3;

const N: usize = 8;

fn grid() -> BoundaryGrid3 {
    BoundaryGrid3::new(N, N, N, 0.83, [0.0; 3], BoundarySpec3::periodic())
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn snapshot(grid: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    let mut result = Vec::new();
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                result.push(grid.populations(x, y, z));
            }
        }
    }
    result
}
fn body(axes: [f64; 3], density: f64) -> MovingPartialSaturationBody3 {
    let mass = density * 4.0 * core::f64::consts::PI * axes.into_iter().product::<f64>() / 3.0;
    let diagonal: [f64; 3] = core::array::from_fn(|i| {
        mass * (axes[(i + 1) % 3].powi(2) + axes[(i + 2) % 3].powi(2)) / 5.0
    });
    MovingPartialSaturationBody3 {
        center: [4.0; 3],
        orientation: So3::identity(),
        mass,
        inertia: InertiaTensor3::new(Mat3 {
            m: [
                diagonal[0],
                0.0,
                0.0,
                0.0,
                diagonal[1],
                0.0,
                0.0,
                0.0,
                diagonal[2],
            ],
        })
        .unwrap(),
        velocity: [0.0; 3],
        angular_momentum_body: Vec3::new(0.0, 0.0, 0.0),
    }
}
fn momentum(
    grid: &BoundaryGrid3,
    fractions: &[f64],
    body: MovingPartialSaturationBody3,
    origin: [f64; 3],
) -> [f64; 6] {
    // Independent physical account: actual all-fluid MINUS sampled interior,
    // PLUS original dry momentum. No analytic rho*V*v substitution.
    let mut total = [0.0; 6];
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let f = grid.populations(x, y, z);
                let exterior = 1.0 - fractions[(z * N + y) * N + x];
                let force = grid.force_density(x, y, z);
                let p: [f64; 3] = core::array::from_fn(|i| {
                    exterior
                        * ((0..Q3)
                            .map(|q| {
                                let e = [E3[q].0, E3[q].1, E3[q].2];
                                f64::from(e[i]) * f[q]
                            })
                            .sum::<f64>()
                            + 0.5 * force[i])
                });
                let angular = cross(
                    [
                        x as f64 + 0.5 - origin[0],
                        y as f64 + 0.5 - origin[1],
                        z as f64 + 0.5 - origin[2],
                    ],
                    p,
                );
                for i in 0..3 {
                    total[i] += p[i];
                    total[i + 3] += angular[i];
                }
            }
        }
    }
    let p = body.velocity.map(|x| body.mass * x);
    let spin = body.orientation.rotate(body.angular_momentum_body).unwrap();
    let orbital = cross(core::array::from_fn(|i| body.center[i] - origin[i]), p);
    let spin = [spin.x, spin.y, spin.z];
    for i in 0..3 {
        total[i] += p[i];
        total[i + 3] += spin[i] + orbital[i];
    }
    total
}

fn fluid_forcing(
    grid: &BoundaryGrid3,
    fractions: &[f64],
    origin: [f64; 3],
) -> ([f64; 6], [f64; 6]) {
    let (mut total, mut interior) = ([0.0; 6], [0.0; 6]);
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let force = grid.force_density(x, y, z);
                let moment = cross(
                    [
                        x as f64 + 0.5 - origin[0],
                        y as f64 + 0.5 - origin[1],
                        z as f64 + 0.5 - origin[2],
                    ],
                    force,
                );
                for i in 0..6 {
                    let value = if i < 3 { force[i] } else { moment[i - 3] };
                    total[i] += value;
                    interior[i] += fractions[(z * N + y) * N + x] * value;
                }
            }
        }
    }
    (total, interior)
}

fn boundary_impulse(report: &DryBodyExchange3, axis: usize) -> f64 {
    if axis < 3 {
        report.native.coupling.walls.fluid_population_impulse[axis]
    } else {
        report
            .native
            .coupling
            .walls
            .fluid_population_angular_impulse[axis - 3]
            + report.native.coupling.periodic.angular_momentum_change[axis - 3]
    }
}

#[test]
fn forced_dry_body_closes_physical_momentum_with_nonuniform_force_and_torque() {
    let axes = [1.2, 0.9, 0.7];
    let shape = BodyFrameEllipsoid3::new(Vec3::new(axes[0], axes[1], axes[2])).unwrap();
    let origin = [-0.3, 0.7, -1.1];
    let gate = CancelGate::new();
    for density in [0.2, 1.0, 3.0] {
        let boundaries = if density == 1.0 {
            BoundarySpec3::new([FaceBoundary3::stationary_wall(); 6])
        } else {
            BoundarySpec3::periodic()
        };
        let mut grid = BoundaryGrid3::new(N, N, N, 0.83, [2e-6, -1e-6, 3e-6], boundaries);
        grid.rebind_force_field(|[x, y, z]| {
            [x as f64 * 3e-6, y as f64 * -2e-6, (z as f64 - 3.0) * 1e-6]
        });
        let mut body = body(axes, density);
        body.center = [3.83, 4.12, 3.91];
        body.velocity = [0.002, -0.001, 0.0015];
        body.orientation = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.15, 0.1))).unwrap();
        body.angular_momentum_body = Vec3::new(1e-4, -2e-4, 1.5e-4);
        let initial = body;
        let mut sampler = BodyFrameSampler3::new(&grid, 4).unwrap();
        let mut midpoint_sampler = BodyFrameSampler3::new(&grid, 4).unwrap();
        for _ in 0..3 {
            sampler
                .sample(body.center, body.orientation, &shape, &gate)
                .unwrap();
            let before = momentum(&grid, sampler.fractions(), body, origin);
            let center = core::array::from_fn(|i| body.center[i] + 0.5 * body.velocity[i]);
            let (orientation, _) = body
                .inertia
                .momentum_step(
                    body.orientation,
                    body.angular_momentum_body,
                    Vec3::new(0.0, 0.0, 0.0),
                    0.5,
                )
                .unwrap();
            midpoint_sampler
                .sample(center, orientation, &shape, &gate)
                .unwrap();
            let (full_force, inside_force) =
                fluid_forcing(&grid, midpoint_sampler.fractions(), origin);
            let external = [2e-6, -3e-6, 4e-6, 5e-6, -2e-6, 3e-6];
            let external_torque = cross(
                core::array::from_fn(|i| center[i] - origin[i]),
                external[..3].try_into().unwrap(),
            );
            let report = grid
                .step_with_dry_partial_saturation_body(
                    &mut sampler,
                    &shape,
                    &mut body,
                    external,
                    origin,
                    DryBodyCouplingConfig3::default(),
                    &gate,
                )
                .unwrap();
            let after = momentum(&grid, sampler.fractions(), body, origin);
            for i in 0..6 {
                let authored = external[i] + if i < 3 { 0.0 } else { external_torque[i - 3] };
                close(report.interior_forcing_impulse[i], inside_force[i], 2e-16);
                close(
                    report.represented_fluid_forcing_impulse[i],
                    full_force[i] - inside_force[i],
                    2e-15,
                );
                close(
                    after[i] - before[i],
                    authored + full_force[i] - inside_force[i] + boundary_impulse(&report, i),
                    3e-10,
                );
                close(report.body_momentum_residual[i], 0.0, 3e-10);
                let actual_correction = report.interior[1].momentum[i]
                    - report.interior[0].momentum[i]
                    - inside_force[i];
                close(
                    report.interior_correction_impulse[i],
                    actual_correction,
                    3e-10,
                );
            }
            assert!(inside_force[..3].iter().any(|x| x.abs() > 1e-6));
            assert!(inside_force[3..].iter().any(|x| x.abs() > 1e-6));
            assert_eq!(body.mass, initial.mass);
            assert_eq!(body.inertia, initial.inertia);
            close(
                report.native.coupling.external_work,
                report.physical_external_work + report.interior_correction_work,
                2e-15,
            );
        }
        assert_ne!(body.center, initial.center);
        assert_ne!(body.angular_momentum_body, initial.angular_momentum_body);
    }
}

#[test]
fn force_on_fictitious_interior_is_not_an_external_physical_body_load() {
    let mut grid = grid();
    let force = [2e-5, -3e-5, 4e-5];
    grid.rebind_force_field(|cell| if cell == [4, 4, 4] { force } else { [0.0; 3] });
    let mut body = body([1.2; 3], 0.5);
    body.center = [4.5; 3];
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let origin = [-0.8, 0.6, 1.3];
    let gate = CancelGate::new();
    let mut sampler = BodyFrameSampler3::new(&grid, 4).unwrap();
    sampler
        .sample(body.center, body.orientation, &shape, &gate)
        .unwrap();
    assert_eq!(sampler.fractions()[(4 * N + 4) * N + 4], 1.0);
    let before = momentum(&grid, sampler.fractions(), body, origin);
    let populations = snapshot(&grid);
    let report = grid
        .step_with_dry_partial_saturation_body(
            &mut sampler,
            &shape,
            &mut body,
            [0.0; 6],
            origin,
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    let after = momentum(&grid, sampler.fractions(), body, origin);
    let torque = cross([4.5 - origin[0], 4.5 - origin[1], 4.5 - origin[2]], force);
    for i in 0..6 {
        let imposed = if i < 3 { force[i] } else { torque[i - 3] };
        close(report.interior_forcing_impulse[i], imposed, 1e-17);
        close(report.represented_fluid_forcing_impulse[i], 0.0, 1e-17);
        close(after[i] - before[i], boundary_impulse(&report, i), 2e-10);
        // Rebinding held the physical fluid at rest. Raw Guo momentum alone
        // would incorrectly report -F/2 inside this fully occupied cell.
        close(report.interior[0].momentum[i], 0.0, 2e-14);
        close(report.fluid[0].momentum[i], 0.0, 2e-12);
    }
    assert!(body.velocity.iter().any(|v| v.abs() > 1e-8));
    assert_ne!(snapshot(&grid), populations);
}

#[test]
fn forced_candidate_budget_and_cancellation_preserve_state_and_exact_retry() {
    let fixture = || {
        let mut grid = grid();
        grid.rebind_force_field(|[x, y, z]| [x as f64 * 1e-6, y as f64 * -2e-6, z as f64 * 3e-6]);
        (grid, body([1.2; 3], 0.5))
    };
    let (mut grid, mut body) = fixture();
    let (mut control_grid, mut control_body) = fixture();
    let original = body;
    let populations = snapshot(&grid);
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let mut sampler = BodyFrameSampler3::new(&grid, 4).unwrap();
    let mut control_sampler = BodyFrameSampler3::new(&control_grid, 4).unwrap();
    let impulse = [0.0, 0.0, 1e-4, 0.0, 0.0, 0.0];
    let gate = CancelGate::new();
    let error = grid
        .step_with_dry_partial_saturation_body(
            &mut sampler,
            &shape,
            &mut body,
            impulse,
            [0.0; 3],
            DryBodyCouplingConfig3 {
                maximum_evaluations: 1,
                ..DryBodyCouplingConfig3::default()
            },
            &gate,
        )
        .unwrap_err();
    assert!(
        matches!(error, PartialSaturationError3::Invalid { reason } if reason.contains("budget"))
    );
    assert_eq!(body, original);
    assert_eq!(snapshot(&grid), populations);
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_eq!(
        grid.step_with_dry_partial_saturation_body(
            &mut sampler,
            &shape,
            &mut body,
            impulse,
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &cancelled,
        ),
        Err(PartialSaturationError3::Cancelled)
    );
    assert_eq!(body, original);
    assert_eq!(snapshot(&grid), populations);
    grid.voxelize_sdf(|_| 1.0);
    let actual = grid
        .step_with_dry_partial_saturation_body(
            &mut sampler,
            &shape,
            &mut body,
            impulse,
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    let expected = control_grid
        .step_with_dry_partial_saturation_body(
            &mut control_sampler,
            &shape,
            &mut control_body,
            impulse,
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(body, control_body);
    assert_eq!(snapshot(&grid), snapshot(&control_grid));
}

#[test]
fn neutrally_buoyant_uniform_translation_preserves_original_dry_momentum() {
    let mut grid = grid();
    let velocity = [0.004, 0.0, 0.0];
    let equilibrium = equilibrium3(1.0, velocity);
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                grid.set_populations(x, y, z, &equilibrium);
            }
        }
    }
    let mut body = body([1.2; 3], 1.0);
    body.velocity = velocity;
    let original = body;
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let gate = CancelGate::new();
    let origin = [-0.6, 1.2, 0.8];
    let mut sampler = BodyFrameSampler3::new(&grid, 4).unwrap();
    sampler
        .sample(body.center, body.orientation, &shape, &gate)
        .unwrap();
    let before = momentum(&grid, sampler.fractions(), body, origin);
    let report = grid
        .step_with_dry_partial_saturation_body(
            &mut sampler,
            &shape,
            &mut body,
            [0.0; 6],
            origin,
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    assert_eq!(body.mass, original.mass);
    assert_eq!(body.inertia, original.inertia);
    assert!(body.mass * body.velocity[0] > 0.0);
    for i in 0..3 {
        close(body.velocity[i], velocity[i], 2e-13);
        close(body.center[i] - original.center[i], velocity[i], 2e-13);
    }
    sampler
        .sample(body.center, body.orientation, &shape, &gate)
        .unwrap();
    let after = momentum(&grid, sampler.fractions(), body, origin);
    for i in 0..6 {
        let boundary = if i < 3 {
            report.native.coupling.walls.fluid_population_impulse[i]
        } else {
            report
                .native
                .coupling
                .walls
                .fluid_population_angular_impulse[i - 3]
                + report.native.coupling.periodic.angular_momentum_change[i - 3]
        };
        close(after[i] - before[i], boundary, 2e-11);
        close(report.combined_momentum_residual[i], 0.0, 2e-11);
    }
    assert!(report.interior_slip_rms[1] < 2e-13);
}

#[test]
fn light_rise_heavy_sink_and_rotated_torque_close_actual_fixed_origin_momentum() {
    let axes = [1.2, 0.9, 0.7];
    let shape = BodyFrameEllipsoid3::new(Vec3::new(axes[0], axes[1], axes[2])).unwrap();
    let volume = 4.0 * core::f64::consts::PI * axes.into_iter().product::<f64>() / 3.0;
    let origin = [0.3, -0.7, 1.1];
    for density in [0.2, 1.0, 3.0] {
        let mut grid = grid();
        let mut body = body(axes, density);
        body.center = [3.83, 4.12, 3.91];
        body.orientation = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.15, 0.1))).unwrap();
        let original = body;
        let impulse = [0.0, 0.0, -(body.mass - volume) * 1e-5, 2e-6, -3e-6, 1e-6];
        let gate = CancelGate::new();
        let mut sampler = BodyFrameSampler3::new(&grid, 4).unwrap();
        sampler
            .sample(body.center, body.orientation, &shape, &gate)
            .unwrap();
        let before = momentum(&grid, sampler.fractions(), body, origin);
        let report = grid
            .step_with_dry_partial_saturation_body(
                &mut sampler,
                &shape,
                &mut body,
                impulse,
                origin,
                DryBodyCouplingConfig3::default(),
                &gate,
            )
            .unwrap();
        assert_eq!(body.mass, original.mass);
        assert_eq!(body.inertia, original.inertia);
        if density != 1.0 {
            assert_eq!(body.velocity[2].is_sign_positive(), density < 1.0);
        }
        assert!(report.candidate_evaluations > 1);
        assert!(report.candidate_evaluations <= 128);
        assert!(
            report
                .interior_correction_impulse
                .iter()
                .any(|x| x.abs() > 1e-10)
        );
        sampler
            .sample(body.center, body.orientation, &shape, &gate)
            .unwrap();
        let after = momentum(&grid, sampler.fractions(), body, origin);
        for i in 0..6 {
            let boundary = if i < 3 {
                report.native.coupling.walls.fluid_population_impulse[i]
            } else {
                report
                    .native
                    .coupling
                    .walls
                    .fluid_population_angular_impulse[i - 3]
                    + report.native.coupling.periodic.angular_momentum_change[i - 3]
            };
            close(
                after[i] - before[i],
                report.physical_external_impulse[i] + boundary,
                2e-10,
            );
            close(report.body_momentum_residual[i], 0.0, 2e-10);
        }
        let expected_kinetic = report.native.kinetic_energy_change + report.fluid[1].kinetic_energy
            - report.fluid[0].kinetic_energy
            - report.interior[1].kinetic_energy
            + report.interior[0].kinetic_energy;
        close(
            report.represented_kinetic_energy_change,
            expected_kinetic,
            2e-16,
        );
    }
}

#[test]
fn exhausted_after_real_candidate_preserves_populations_dry_body_and_topology() {
    let mut grid = grid();
    let mut body = body([1.2; 3], 0.5);
    let original = body;
    let original_populations = snapshot(&grid);
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let gate = CancelGate::new();
    let mut sampler = BodyFrameSampler3::new(&grid, 4).unwrap();
    let impulse = [0.0, 0.0, 1e-4, 0.0, 0.0, 0.0];
    let controls = DryBodyCouplingConfig3 {
        maximum_evaluations: 1,
        ..DryBodyCouplingConfig3::default()
    };
    let error = grid
        .step_with_dry_partial_saturation_body(
            &mut sampler,
            &shape,
            &mut body,
            impulse,
            [0.0; 3],
            controls,
            &gate,
        )
        .unwrap_err();
    assert!(matches!(error, PartialSaturationError3::Invalid { reason }
        if reason.contains("budget")));
    assert_eq!(body, original);
    assert_eq!(snapshot(&grid), original_populations);
    // Topology is still unlocked: rejection did not publish a hidden first tick.
    grid.voxelize_sdf(|_| 1.0);
    let report = grid
        .step_with_dry_partial_saturation_body(
            &mut sampler,
            &shape,
            &mut body,
            impulse,
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    assert!(report.candidate_evaluations > 1);
    assert!(body.velocity[2] > 0.0);
    assert_ne!(snapshot(&grid), original_populations);
}
