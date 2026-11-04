//! Simultaneous physical dry inertia, shared-fluid feedback and atomic refusal.
use fs_exec::CancelGate;
use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_lbm::d3q19::{
    BodyFrameEllipsoid3, BodyFrameSampler3, BodyFrameShape3, BoundaryGrid3, BoundarySpec3,
    DryBodiesExchange3, DryBodyCouplingConfig3, E3, FaceBoundary3, MovingPartialSaturationBody3,
    PartialSaturationError3, Q3, equilibrium3,
};
use fs_time::lie::gyrostat::InertiaTensor3;

fn grid() -> BoundaryGrid3 {
    BoundaryGrid3::new(12, 8, 8, 0.83, [0.0; 3], BoundarySpec3::periodic())
}
fn body(center: [f64; 3], density: f64) -> MovingPartialSaturationBody3 {
    let mass = density * 4.0 * core::f64::consts::PI * 1.2_f64.powi(3) / 3.0;
    MovingPartialSaturationBody3 {
        center,
        orientation: So3::identity(),
        mass,
        inertia: InertiaTensor3::new(Mat3 {
            m: [0.4, 0.02, -0.01, 0.02, 0.43, 0.015, -0.01, 0.015, 0.46]
                .map(|x| mass * 1.2_f64.powi(2) * x),
        })
        .unwrap(),
        velocity: [0.0; 3],
        angular_momentum_body: Vec3::new(0.0, 0.0, 0.0),
    }
}
fn pair() -> [MovingPartialSaturationBody3; 2] {
    [body([3.5, 4.5, 4.5], 0.3), body([6.5, 4.5, 4.5], 1.7)]
}
fn samplers<const N: usize>(grid: &BoundaryGrid3) -> [BodyFrameSampler3; N] {
    core::array::from_fn(|_| BodyFrameSampler3::new(grid, 4).unwrap())
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
    let [nx, ny, nz] = grid.dimensions();
    let mut result = Vec::new();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                result.push(grid.populations(x, y, z));
            }
        }
    }
    result
}
fn sample<const N: usize>(
    samplers: &mut [BodyFrameSampler3; N],
    bodies: [MovingPartialSaturationBody3; N],
    shape: &impl BodyFrameShape3,
    midpoint: bool,
    gate: &CancelGate,
) {
    for a in 0..N {
        let b = bodies[a];
        let (center, orientation) = if midpoint {
            let (q, _) = b
                .inertia
                .momentum_step(
                    b.orientation,
                    b.angular_momentum_body,
                    Vec3::new(0.0, 0.0, 0.0),
                    0.5,
                )
                .unwrap();
            (
                core::array::from_fn(|i| b.center[i] + 0.5 * b.velocity[i]),
                q,
            )
        } else {
            (b.center, b.orientation)
        };
        samplers[a]
            .sample(center, orientation, shape, gate)
            .unwrap();
    }
}
fn momentum<const N: usize>(
    grid: &BoundaryGrid3,
    fractions: &[BodyFrameSampler3; N],
    bodies: [MovingPartialSaturationBody3; N],
    origin: [f64; 3],
) -> [f64; 6] {
    // Independent row-major exterior PHYSICAL fluid plus original dry inertia.
    let [nx, ny, nz] = grid.dimensions();
    let mut total = [0.0; 6];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let index = (z * ny + y) * nx + x;
                let exterior = 1.0 - fractions.iter().map(|s| s.fractions()[index]).sum::<f64>();
                let f = grid.populations(x, y, z);
                let force = grid.force_density(x, y, z);
                let density: [f64; 3] = core::array::from_fn(|i| {
                    (0..Q3)
                        .map(|q| {
                            let e = [E3[q].0, E3[q].1, E3[q].2];
                            f64::from(e[i]) * f[q]
                        })
                        .sum::<f64>()
                        + 0.5 * force[i]
                });
                let p = density.map(|p| exterior * p);
                let angular = cross(
                    [
                        x as f64 + 0.5 - origin[0],
                        y as f64 + 0.5 - origin[1],
                        z as f64 + 0.5 - origin[2],
                    ],
                    p,
                );
                let first: [f64; 3] = core::array::from_fn(|i| {
                    fractions
                        .iter()
                        .map(|sampler| {
                            sampler
                                .cut_cell_first_moments()
                                .iter()
                                .find(|&&(cell, _)| cell == index)
                                .map_or(0.0, |&(_, moment)| moment[i])
                        })
                        .sum()
                });
                let interior_cut = cross(first, density);
                for i in 0..3 {
                    total[i] += p[i];
                    total[i + 3] += angular[i] - interior_cut[i];
                }
            }
        }
    }
    for b in bodies {
        let p = b.velocity.map(|v| b.mass * v);
        let spin = b.orientation.rotate(b.angular_momentum_body).unwrap();
        let orbital = cross(core::array::from_fn(|i| b.center[i] - origin[i]), p);
        let spin = [spin.x, spin.y, spin.z];
        for i in 0..3 {
            total[i] += p[i];
            total[i + 3] += orbital[i] + spin[i];
        }
    }
    total
}
fn forcing<const N: usize>(
    grid: &BoundaryGrid3,
    midpoint: &[BodyFrameSampler3; N],
    origin: [f64; 3],
) -> ([f64; 6], [[f64; 6]; N]) {
    let [nx, ny, nz] = grid.dimensions();
    let mut full = [0.0; 6];
    let mut inside = [[0.0; 6]; N];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let force = grid.force_density(x, y, z);
                let torque = cross(
                    [
                        x as f64 + 0.5 - origin[0],
                        y as f64 + 0.5 - origin[1],
                        z as f64 + 0.5 - origin[2],
                    ],
                    force,
                );
                for i in 0..6 {
                    let value = if i < 3 { force[i] } else { torque[i - 3] };
                    full[i] += value;
                    for a in 0..N {
                        inside[a][i] += midpoint[a].fractions()[(z * ny + y) * nx + x] * value;
                    }
                }
                for a in 0..N {
                    if let Some(&(_, first)) = midpoint[a]
                        .cut_cell_first_moments()
                        .iter()
                        .find(|&&(cell, _)| cell == (z * ny + y) * nx + x)
                    {
                        let cut_torque = cross(first, force);
                        for i in 0..3 {
                            inside[a][i + 3] += cut_torque[i];
                        }
                    }
                }
            }
        }
    }
    (full, inside)
}
fn boundary<const N: usize>(r: &DryBodiesExchange3<N>, i: usize) -> f64 {
    if i < 3 {
        r.native.coupling.walls.fluid_population_impulse[i]
    } else {
        r.native.coupling.walls.fluid_population_angular_impulse[i - 3]
            + r.native.coupling.periodic.angular_momentum_change[i - 3]
    }
}

#[test]
fn singleton_delegates_exactly_including_forcing_and_work() {
    let mut one = grid();
    let mut assembly = grid();
    for grid in [&mut one, &mut assembly] {
        grid.rebind_force_field(|[x, y, z]| [x as f64 * 1e-6, y as f64 * -2e-6, z as f64 * 3e-6]);
    }
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let mut b = pair()[0];
    b.velocity = [0.001, -0.0005, 0.0003];
    let mut bodies = [b];
    let mut s = BodyFrameSampler3::new(&one, 4).unwrap();
    let mut fields = samplers(&assembly);
    let origin = [-0.7, 0.4, 1.1];
    let load = [1e-5, -2e-5, 3e-5, 2e-6, 3e-6, -4e-6];
    let gate = CancelGate::new();
    for _ in 0..3 {
        let single = one
            .step_with_dry_partial_saturation_body(
                &mut s,
                &shape,
                &mut b,
                load,
                origin,
                DryBodyCouplingConfig3::default(),
                &gate,
            )
            .unwrap();
        let multi = assembly
            .step_with_dry_partial_saturation_bodies(
                &mut fields,
                [&shape],
                &mut bodies,
                [load],
                origin,
                DryBodyCouplingConfig3::default(),
                &gate,
            )
            .unwrap();
        assert_eq!(b, bodies[0]);
        assert_eq!(snapshot(&one), snapshot(&assembly));
        assert_eq!(s.fractions(), fields[0].fractions());
        assert_eq!(single.fluid, multi.fluid);
        assert_eq!(single.interior, multi.interior[0]);
        assert_eq!(
            single.interior_correction_impulse,
            multi.interior_correction_impulses[0]
        );
        assert_eq!(
            single.interior_forcing_impulse,
            multi.interior_forcing_impulses[0]
        );
        assert_eq!(
            single.represented_fluid_forcing_impulse,
            multi.represented_fluid_forcing_impulse
        );
        assert_eq!(
            single.corrected_hydrodynamic_impulse,
            multi.corrected_hydrodynamic_impulses[0]
        );
        assert_eq!(single.body_momentum_change, multi.body_momentum_changes[0]);
        assert_eq!(
            single.combined_momentum_residual,
            multi.combined_momentum_residual
        );
        assert_eq!(
            single.physical_external_work,
            multi.physical_external_work[0]
        );
        assert_eq!(
            single.interior_correction_work,
            multi.interior_correction_work[0]
        );
        assert_eq!(
            single.represented_kinetic_energy_change,
            multi.represented_kinetic_energy_change
        );
        assert_eq!(single.candidate_evaluations, multi.candidate_evaluations);
        assert_eq!(
            single.native.coupling.forcing,
            multi.native.coupling.forcing
        );
        assert_eq!(single.native.coupling.walls, multi.native.coupling.walls);
        assert_eq!(multi.native.bodies[0].coupling.forcing, Default::default());
        assert_eq!(multi.native.bodies[0].coupling.walls, Default::default());
        assert_eq!(multi.native.bodies[0].coupling.periodic, Default::default());
    }
}

#[test]
fn forced_full_tensor_dry_assembly_closes_at_translated_origin_over_repeated_steps() {
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let origin = [-0.4, 0.8, -1.3];
    let gate = CancelGate::new();
    let loads = [
        [2e-6, -3e-6, 4e-6, 1e-6, -2e-6, 3e-6],
        [-5e-6, 3e-6, 2e-6, -3e-6, 4e-6, -2e-6],
    ];
    for density in [0.2, 1.0, 3.0] {
        let boundaries = if density == 1.0 {
            BoundarySpec3::new([FaceBoundary3::stationary_wall(); 6])
        } else {
            BoundarySpec3::periodic()
        };
        let mut grid = BoundaryGrid3::new(12, 8, 8, 0.83, [2e-6, -1e-6, 3e-6], boundaries);
        grid.rebind_force_field(|[x, y, z]| {
            [x as f64 * 3e-6, y as f64 * -2e-6, (z as f64 - 3.0) * 1e-6]
        });
        let mut bodies = [
            body([3.53, 4.13, 3.91], density),
            body([6.57, 4.19, 3.89], 1.7 * density),
        ];
        for (a, b) in bodies.iter_mut().enumerate() {
            b.velocity = [0.0007, -0.0003, 0.0004 * (a + 1) as f64];
            b.orientation = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.15, 0.1))).unwrap();
            b.angular_momentum_body = Vec3::new(1e-4, -2e-4, 1.5e-4);
        }
        let initial = bodies;
        let mut fields = samplers(&grid);
        let mut midpoint = samplers(&grid);
        for _ in 0..3 {
            sample(&mut fields, bodies, &shape, false, &gate);
            let before = momentum(&grid, &fields, bodies, origin);
            sample(&mut midpoint, bodies, &shape, true, &gate);
            let (full, inside) = forcing(&grid, &midpoint, origin);
            let expected_loads: [[f64; 6]; 2] = core::array::from_fn(|a| {
                let arm = core::array::from_fn(|i| {
                    bodies[a].center[i] + 0.5 * bodies[a].velocity[i] - origin[i]
                });
                let torque = cross(arm, loads[a][..3].try_into().unwrap());
                core::array::from_fn(|i| loads[a][i] + if i < 3 { 0.0 } else { torque[i - 3] })
            });
            let report = grid
                .step_with_dry_partial_saturation_bodies(
                    &mut fields,
                    [&shape, &shape],
                    &mut bodies,
                    loads,
                    origin,
                    DryBodyCouplingConfig3::default(),
                    &gate,
                )
                .unwrap();
            let after = momentum(&grid, &fields, bodies, origin);
            for i in 0..6 {
                let exterior_force = full[i] - inside.iter().map(|f| f[i]).sum::<f64>();
                close(
                    report.represented_fluid_forcing_impulse[i],
                    exterior_force,
                    2e-14,
                );
                close(
                    after[i] - before[i],
                    expected_loads.iter().map(|f| f[i]).sum::<f64>()
                        + exterior_force
                        + boundary(&report, i),
                    5e-10,
                );
                close(report.combined_momentum_residual[i], 0.0, 5e-10);
                for a in 0..2 {
                    close(report.interior_forcing_impulses[a][i], inside[a][i], 2e-15);
                    close(
                        report.physical_external_impulses[a][i],
                        expected_loads[a][i],
                        2e-15,
                    );
                    close(report.body_momentum_residuals[a][i], 0.0, 3e-10);
                    close(
                        report.interior_correction_impulses[a][i],
                        report.interior[a][1].momentum[i]
                            - report.interior[a][0].momentum[i]
                            - inside[a][i],
                        3e-10,
                    );
                }
            }
            for a in 0..2 {
                assert_eq!(bodies[a].mass, initial[a].mass);
                assert_eq!(bodies[a].inertia, initial[a].inertia);
                close(
                    report.native.bodies[a].coupling.external_work,
                    report.physical_external_work[a] + report.interior_correction_work[a],
                    2e-15,
                );
                assert!(inside[a][3..].iter().any(|x| x.abs() > 1e-6));
            }
        }
        assert_ne!(bodies[0].center, initial[0].center);
        assert_ne!(bodies[1].orientation, initial[1].orientation);
    }
}

#[test]
fn full_cross_body_jacobian_transmits_streamed_impulse_in_one_newton_correction() {
    let mut grid = grid();
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let mut bodies = pair();
    let mut fields = samplers(&grid);
    let config = DryBodyCouplingConfig3 {
        maximum_iterations: 1,
        maximum_evaluations: 14,
        ..DryBodyCouplingConfig3::default()
    };
    let report = grid
        .step_with_dry_partial_saturation_bodies(
            &mut fields,
            [&shape, &shape],
            &mut bodies,
            [[1e-4, 2e-5, -3e-5, 1e-5, -2e-5, 3e-5], [0.0; 6]],
            [0.0; 3],
            config,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(report.iterations, 1);
    assert_eq!(report.candidate_evaluations, 14);
    // Body two has no authored load or incoming slip. Its motion comes from
    // interior fluid streamed from body one in this SAME composite candidate.
    assert!(bodies[1].velocity.iter().any(|v| v.abs() > 1e-9));
    assert!(
        report.interior_correction_impulses[1]
            .iter()
            .any(|j| j.abs() > 1e-8)
    );
    for defect in report.body_momentum_residuals.into_iter().flatten() {
        close(defect, 0.0, 1e-10);
    }
}

#[test]
fn each_fictitious_interior_source_is_removed_once_from_physical_assembly() {
    let mut grid = grid();
    grid.rebind_force_field(|cell| match cell {
        [3, 4, 4] => [2e-5, -3e-5, 4e-5],
        [6, 4, 4] => [-1e-5, 2e-5, -3e-5],
        _ => [0.0; 3],
    });
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let mut bodies = pair();
    let mut fields = samplers(&grid);
    let gate = CancelGate::new();
    let origin = [-0.8, 0.6, 1.3];
    sample(&mut fields, bodies, &shape, false, &gate);
    assert_eq!(fields[0].fractions()[(4 * 8 + 4) * 12 + 3], 1.0);
    assert_eq!(fields[1].fractions()[(4 * 8 + 4) * 12 + 6], 1.0);
    let before = momentum(&grid, &fields, bodies, origin);
    let (full, inside) = forcing(&grid, &fields, origin);
    let report = grid
        .step_with_dry_partial_saturation_bodies(
            &mut fields,
            [&shape, &shape],
            &mut bodies,
            [[0.0; 6]; 2],
            origin,
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    let after = momentum(&grid, &fields, bodies, origin);
    for i in 0..6 {
        close(full[i] - inside[0][i] - inside[1][i], 0.0, 1e-17);
        close(report.represented_fluid_forcing_impulse[i], 0.0, 1e-17);
        close(after[i] - before[i], boundary(&report, i), 2e-10);
        close(report.fluid[0].momentum[i], 0.0, 2e-12);
        for a in 0..2 {
            close(report.interior_forcing_impulses[a][i], inside[a][i], 1e-17);
            close(report.interior[a][0].momentum[i], 0.0, 2e-14);
        }
    }
}

#[test]
fn permutation_changes_only_body_order_and_population_roundoff() {
    let mut grid_a = grid();
    let mut grid_b = grid();
    for grid in [&mut grid_a, &mut grid_b] {
        grid.rebind_force_field(|[x, y, z]| [x as f64 * 1e-6, y as f64 * -2e-6, z as f64 * 3e-6]);
    }
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let mut a = pair();
    let mut b = [a[1], a[0]];
    let mut sa = samplers(&grid_a);
    let mut sb = samplers(&grid_b);
    let loads = [
        [1e-4, 2e-5, -3e-5, 1e-5, -2e-5, 3e-5],
        [-3e-5, 1e-5, 2e-5, -2e-5, 1e-5, 4e-5],
    ];
    let gate = CancelGate::new();
    for _ in 0..2 {
        let _ = grid_a
            .step_with_dry_partial_saturation_bodies(
                &mut sa,
                [&shape, &shape],
                &mut a,
                loads,
                [-0.5, 0.7, 1.2],
                DryBodyCouplingConfig3::default(),
                &gate,
            )
            .unwrap();
        let _ = grid_b
            .step_with_dry_partial_saturation_bodies(
                &mut sb,
                [&shape, &shape],
                &mut b,
                [loads[1], loads[0]],
                [-0.5, 0.7, 1.2],
                DryBodyCouplingConfig3::default(),
                &gate,
            )
            .unwrap();
        for i in 0..2 {
            let left = a[i];
            let right = b[1 - i];
            for (x, y) in left
                .center
                .into_iter()
                .chain(left.velocity)
                .chain(left.orientation.matrix().m)
                .chain([
                    left.angular_momentum_body.x,
                    left.angular_momentum_body.y,
                    left.angular_momentum_body.z,
                ])
                .zip(
                    right
                        .center
                        .into_iter()
                        .chain(right.velocity)
                        .chain(right.orientation.matrix().m)
                        .chain([
                            right.angular_momentum_body.x,
                            right.angular_momentum_body.y,
                            right.angular_momentum_body.z,
                        ]),
                )
            {
                close(x, y, 3e-11);
            }
        }
        for (left, right) in snapshot(&grid_a)
            .into_iter()
            .flatten()
            .zip(snapshot(&grid_b).into_iter().flatten())
        {
            close(left, right, 3e-12);
        }
    }
}

#[test]
fn three_translating_dry_bodies_keep_their_original_mass_and_full_inertia() {
    let mut grid = BoundaryGrid3::new(16, 8, 8, 0.83, [0.0; 3], BoundarySpec3::periodic());
    let velocity = [0.002, -0.001, 0.0005];
    let equilibrium = equilibrium3(1.0, velocity);
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..16 {
                grid.set_populations(x, y, z, &equilibrium);
            }
        }
    }
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let mut bodies = [
        body([3.5, 4.5, 4.5], 0.2),
        body([7.0, 4.5, 4.5], 1.0),
        body([10.5, 4.5, 4.5], 3.0),
    ];
    for b in &mut bodies {
        b.velocity = velocity;
    }
    let before = bodies;
    let mut fields = samplers(&grid);
    let report = grid
        .step_with_dry_partial_saturation_bodies(
            &mut fields,
            [&shape, &shape, &shape],
            &mut bodies,
            [[0.0; 6]; 3],
            [-0.7, 0.9, 1.2],
            DryBodyCouplingConfig3::default(),
            &CancelGate::new(),
        )
        .unwrap();
    for a in 0..3 {
        assert_eq!(bodies[a].mass, before[a].mass);
        assert_eq!(bodies[a].inertia, before[a].inertia);
        for i in 0..3 {
            close(bodies[a].velocity[i], velocity[i], 3e-13);
            close(
                bodies[a].center[i] - before[a].center[i],
                velocity[i],
                3e-13,
            );
        }
        close(report.interior_slip_rms[a][1], 0.0, 3e-13);
    }
    assert_eq!(report.candidate_evaluations, 1);
}

#[test]
fn failed_probes_conditioning_and_overlap_leave_all_state_retryable() {
    let mut grid = grid();
    let mut control = crate::grid();
    for grid in [&mut grid, &mut control] {
        grid.rebind_force_field(|[x, y, z]| [x as f64 * 1e-6, y as f64 * -2e-6, z as f64 * 3e-6]);
    }
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let mut bodies = pair();
    let initial = bodies;
    let mut control_bodies = bodies;
    let mut fields = samplers(&grid);
    let mut control_fields = samplers(&control);
    let populations = snapshot(&grid);
    let loads = [[1e-4, 2e-5, -3e-5, 1e-5, -2e-5, 3e-5], [0.0; 6]];
    let gate = CancelGate::new();
    for (config, reason) in [
        (
            DryBodyCouplingConfig3 {
                maximum_evaluations: 2,
                ..DryBodyCouplingConfig3::default()
            },
            "budget",
        ),
        (
            DryBodyCouplingConfig3 {
                maximum_condition_number: 1.0,
                ..DryBodyCouplingConfig3::default()
            },
            "ill-conditioned",
        ),
    ] {
        let error = grid
            .step_with_dry_partial_saturation_bodies(
                &mut fields,
                [&shape, &shape],
                &mut bodies,
                loads,
                [0.0; 3],
                config,
                &gate,
            )
            .unwrap_err();
        assert!(
            matches!(error,PartialSaturationError3::Invalid {reason:actual} if actual.contains(reason)),
            "{error:?}"
        );
        assert_eq!(bodies, initial);
        assert_eq!(snapshot(&grid), populations);
    }
    let mut overlapping = [initial[0]; 2];
    let error = grid
        .step_with_dry_partial_saturation_bodies(
            &mut fields,
            [&shape, &shape],
            &mut overlapping,
            loads,
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap_err();
    assert!(
        matches!(error,PartialSaturationError3::Invalid {reason} if reason.contains("ownership"))
    );
    assert_eq!(overlapping, [initial[0]; 2]);
    assert_eq!(snapshot(&grid), populations);
    // Failed candidates also leave topology admission open.
    grid.voxelize_sdf(|_| 1.0);
    let actual = grid
        .step_with_dry_partial_saturation_bodies(
            &mut fields,
            [&shape, &shape],
            &mut bodies,
            loads,
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    let expected = control
        .step_with_dry_partial_saturation_bodies(
            &mut control_fields,
            [&shape, &shape],
            &mut control_bodies,
            loads,
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(bodies, control_bodies);
    assert_eq!(snapshot(&grid), snapshot(&control));
}

#[test]
fn cancellation_after_streaming_publishes_no_member_and_allows_exact_retry() {
    use std::cell::Cell;
    struct CancellingShape<'a> {
        shape: BodyFrameEllipsoid3,
        gate: &'a CancelGate,
        calls: Cell<usize>,
        cancel_at: Cell<usize>,
    }
    impl BodyFrameShape3 for CancellingShape<'_> {
        fn bounding_radius(&self) -> f64 {
            self.shape.bounding_radius()
        }
        fn contains(&self, p: Vec3) -> bool {
            let calls = self.calls.get() + 1;
            self.calls.set(calls);
            if calls == self.cancel_at.get() {
                self.gate.request();
            }
            self.shape.contains(p)
        }
    }
    let mut grid = grid();
    let mut control = crate::grid();
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let gate = CancelGate::new();
    let cancelling = CancellingShape {
        shape,
        gate: &gate,
        calls: Cell::new(0),
        cancel_at: Cell::new(usize::MAX),
    };
    let mut bodies = pair();
    let initial = bodies;
    let mut control_bodies = bodies;
    let mut fields = samplers(&grid);
    let mut control_fields = samplers(&control);
    fields[0]
        .sample(bodies[0].center, bodies[0].orientation, &cancelling, &gate)
        .unwrap();
    let calls = cancelling.calls.get();
    assert!(calls > 0);
    // Initial inventory plus native start/midpoint/endpoint sampling complete
    // before the dry owner resamples the streamed endpoint inventory.
    cancelling.calls.set(0);
    cancelling.cancel_at.set(4 * calls + calls / 2);
    let populations = snapshot(&grid);
    let loads = [[1e-4, 2e-5, -3e-5, 1e-5, -2e-5, 3e-5], [0.0; 6]];
    let error = grid
        .step_with_dry_partial_saturation_bodies(
            &mut fields,
            [&cancelling, &shape],
            &mut bodies,
            loads,
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap_err();
    assert_eq!(error, PartialSaturationError3::Cancelled);
    assert!(cancelling.calls.get() > 4 * calls);
    assert_eq!(bodies, initial);
    assert_eq!(snapshot(&grid), populations);
    grid.voxelize_sdf(|_| 1.0);
    let gate = CancelGate::new();
    let actual = grid
        .step_with_dry_partial_saturation_bodies(
            &mut fields,
            [&shape, &shape],
            &mut bodies,
            loads,
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    let expected = control
        .step_with_dry_partial_saturation_bodies(
            &mut control_fields,
            [&shape, &shape],
            &mut control_bodies,
            loads,
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(bodies, control_bodies);
    assert_eq!(snapshot(&grid), snapshot(&control));
}

#[test]
fn empty_and_oversized_dry_assemblies_are_refused_before_publication() {
    let mut grid = grid();
    let populations = snapshot(&grid);
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let mut bodies = [pair()[0]; 9];
    let initial = bodies;
    let mut fields = samplers(&grid);
    assert!(
        grid.step_with_dry_partial_saturation_bodies(
            &mut fields,
            [&shape; 9],
            &mut bodies,
            [[0.0; 6]; 9],
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, initial);
    assert_eq!(snapshot(&grid), populations);
    assert!(
        grid.step_with_dry_partial_saturation_bodies::<0>(
            &mut [],
            [],
            &mut [],
            [],
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(snapshot(&grid), populations);
    grid.voxelize_sdf(|_| 1.0);
}

fn point_quadrature(
    grid: &BoundaryGrid3,
    body: MovingPartialSaturationBody3,
    shape: &impl BodyFrameShape3,
    samples: usize,
    origin: [f64; 3],
    force_only: bool,
    collapsed_to_center: bool,
) -> [f64; 6] {
    let [nx, ny, nz] = grid.dimensions();
    let weight = 1.0 / (samples * samples * samples) as f64;
    let mut total = [0.0; 6];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let force = grid.force_density(x, y, z);
                let f = grid.populations(x, y, z);
                let p: [f64; 3] = if force_only {
                    force
                } else {
                    core::array::from_fn(|i| {
                        (0..Q3)
                            .map(|q| {
                                let e = [E3[q].0, E3[q].1, E3[q].2];
                                f64::from(e[i]) * f[q]
                            })
                            .sum::<f64>()
                            + 0.5 * force[i]
                    })
                };
                for sz in 0..samples {
                    for sy in 0..samples {
                        for sx in 0..samples {
                            let cell = [x, y, z];
                            let subcell = [sx, sy, sz];
                            let point: [f64; 3] = core::array::from_fn(|i| {
                                cell[i] as f64 + (subcell[i] as f64 + 0.5) / samples as f64
                            });
                            let local = body
                                .orientation
                                .inverse()
                                .rotate(Vec3::new(
                                    point[0] - body.center[0],
                                    point[1] - body.center[1],
                                    point[2] - body.center[2],
                                ))
                                .unwrap();
                            if shape.contains(local) {
                                let arm = core::array::from_fn(|i| {
                                    (if collapsed_to_center {
                                        cell[i] as f64 + 0.5
                                    } else {
                                        point[i]
                                    }) - origin[i]
                                });
                                let angular = cross(arm, p);
                                for i in 0..3 {
                                    total[i] += weight * p[i];
                                    total[i + 3] += weight * angular[i];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    total
}

#[test]
fn cut_cell_dry_assembly_matches_point_torques_and_origin_covariance() {
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let gate = CancelGate::new();
    for forced in [false, true] {
        let mut reference: Option<([MovingPartialSaturationBody3; 2], [[f64; 6]; 2])> = None;
        for origin in [[0.0; 3], [4.3, -2.7, 1.9]] {
            let mut grid = grid();
            for z in 0..8 {
                for y in 0..8 {
                    for x in 0..12 {
                        grid.set_populations(
                            x,
                            y,
                            z,
                            &equilibrium3(
                                1.0 + 0.002 * (x + y + z) as f64,
                                [0.002 * y as f64, -0.001 * z as f64, 0.0015 * x as f64],
                            ),
                        );
                    }
                }
            }
            if forced {
                grid.rebind_force_field(|[x, y, z]| {
                    [x as f64 * 3e-6, y as f64 * -2e-6, (z as f64 - 3.0) * 1e-6]
                });
            }
            let mut bodies = [body([3.53, 4.13, 3.91], 0.3), body([6.57, 4.19, 3.89], 1.7)];
            for b in &mut bodies {
                b.velocity = [0.0007, -0.0003, 0.0004];
                b.angular_momentum_body = Vec3::new(1e-4, -2e-4, 1.5e-4);
            }
            let mut fields = samplers(&grid);
            let before =
                bodies.map(|b| point_quadrature(&grid, b, &shape, 4, origin, false, false));
            let collapsed =
                bodies.map(|b| point_quadrature(&grid, b, &shape, 4, origin, false, true));
            assert!(
                (0..2).any(|a| { (3..6).any(|i| (before[a][i] - collapsed[a][i]).abs() > 1e-5) }),
                "cut-cell angular momentum must differ from a cell-center lump"
            );
            let midpoint = bodies.map(|mut b| {
                b.center = core::array::from_fn(|i| b.center[i] + 0.5 * b.velocity[i]);
                b.orientation = b
                    .inertia
                    .momentum_step(
                        b.orientation,
                        b.angular_momentum_body,
                        Vec3::new(0.0, 0.0, 0.0),
                        0.5,
                    )
                    .unwrap()
                    .0;
                b
            });
            let source =
                midpoint.map(|b| point_quadrature(&grid, b, &shape, 4, origin, true, false));
            let report = grid
                .step_with_dry_partial_saturation_bodies(
                    &mut fields,
                    [&shape, &shape],
                    &mut bodies,
                    [[0.0; 6]; 2],
                    origin,
                    DryBodyCouplingConfig3::default(),
                    &gate,
                )
                .unwrap();
            for a in 0..2 {
                let after = point_quadrature(&grid, bodies[a], &shape, 4, origin, false, false);
                for i in 0..6 {
                    close(report.interior[a][0].momentum[i], before[a][i], 3e-12);
                    close(report.interior[a][1].momentum[i], after[i], 3e-12);
                    close(report.interior_forcing_impulses[a][i], source[a][i], 3e-16);
                    close(report.combined_momentum_residual[i], 0.0, 5e-10);
                }
            }
            if let Some((first_bodies, first_source)) = reference {
                for a in 0..2 {
                    let shift = cross(origin, first_source[a][..3].try_into().unwrap());
                    for i in 0..3 {
                        close(source[a][i + 3], first_source[a][i + 3] - shift[i], 3e-16);
                        close(bodies[a].velocity[i], first_bodies[a].velocity[i], 2e-11);
                    }
                    let delta =
                        bodies[a].angular_momentum_body - first_bodies[a].angular_momentum_body;
                    close(delta.norm(), 0.0, 2e-11);
                }
            } else {
                reference = Some((bodies, source));
            }
        }
    }
}

fn history_field([x, y, z]: [usize; 3]) -> [f64; 3] {
    [(y as f64 - 3.1) * 3e-6, x as f64 * -2e-6, (z as f64 - 2.0) * 1e-6]
}

fn force_snapshot(grid: &BoundaryGrid3) -> Vec<[f64; 3]> {
    let [nx, ny, nz] = grid.dimensions();
    (0..nx * ny * nz)
        .map(|i| grid.force_density(i % nx, (i / nx) % ny, i / (nx * ny)))
        .collect()
}

#[test]
fn endpoint_history_closes_independent_dry_assembly_momentum_and_interior_sources() {
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let origin = [-0.4, 0.8, -1.3];
    let gate = CancelGate::new();
    let uniform = [2e-6, -1e-6, 3e-6];
    let loads = [
        [2e-6, -3e-6, 4e-6, 1e-6, -2e-6, 3e-6],
        [-5e-6, 3e-6, 2e-6, -3e-6, 4e-6, -2e-6],
    ];
    for boundaries in [
        BoundarySpec3::periodic(),
        BoundarySpec3::new([FaceBoundary3::stationary_wall(); 6]),
        BoundarySpec3::velocity_pressure_duct_z([0.001, -0.002, 0.012], 1.002),
    ] {
        let mut grid = BoundaryGrid3::new(12, 8, 8, 0.83, uniform, boundaries);
        grid.rebind_force_field(history_field);
        let mut bodies = pair();
        bodies[0].center = [3.41, 3.92, 4.17];
        bodies[1].center = [6.61, 4.09, 3.91];
        for (a, b) in bodies.iter_mut().enumerate() {
            b.velocity = [0.0001, -0.0002, 0.00015];
            b.orientation = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.15, 0.1))).unwrap();
            b.angular_momentum_body = Vec3::new(1e-6, -2e-6, (a + 1) as f64 * 1.5e-6);
        }
        let initial = bodies;
        let mut fields = samplers(&grid);
        let mut independent = samplers(&grid);
        let mut midpoint = samplers(&grid);
        for scale in [2.0, 0.0, -1.0, 0.5] {
            sample(&mut independent, bodies, &shape, false, &gate);
            let before = momentum(&grid, &independent, bodies, origin);
            sample(&mut midpoint, bodies, &shape, true, &gate);
            let (old_full, old_inside) = forcing(&grid, &midpoint, origin);
            let authored: [[f64; 6]; 2] = core::array::from_fn(|a| {
                let arm = core::array::from_fn(|i| {
                    bodies[a].center[i] + 0.5 * bodies[a].velocity[i] - origin[i]
                });
                let torque = cross(arm, loads[a][..3].try_into().unwrap());
                core::array::from_fn(|i| loads[a][i] + if i < 3 { 0.0 } else { torque[i - 3] })
            });
            let report = grid
                .step_with_dry_partial_saturation_bodies_and_endpoint_force_scale(
                    &mut fields, [&shape, &shape], &mut bodies, loads, scale,
                    origin, DryBodyCouplingConfig3::default(), &gate,
                )
                .unwrap();
            sample(&mut independent, bodies, &shape, false, &gate);
            let after = momentum(&grid, &independent, bodies, origin);
            let (new_full, new_inside) = forcing(&grid, &midpoint, origin);
            for i in 0..6 {
                let full = 0.5 * old_full[i] + 0.5 * new_full[i];
                let inside: [f64; 2] =
                    core::array::from_fn(|a| 0.5 * old_inside[a][i] + 0.5 * new_inside[a][i]);
                let exterior = full - inside.iter().sum::<f64>();
                let forcing = report.native.coupling.forcing;
                let open = report.native.coupling.open;
                let source = if i < 3 { forcing.impulse[i] } else { forcing.angular_impulse[i - 3] };
                let reservoir = if i < 3 { open.impulse[i] } else { open.angular_impulse[i - 3] };
                close(source, full, 3e-15);
                close(report.represented_fluid_forcing_impulse[i], exterior, 4e-15);
                close(after[i] - before[i],
                    authored.iter().map(|j| j[i]).sum::<f64>()
                        + exterior + boundary(&report, i) + reservoir, 5e-10);
                for a in 0..2 {
                    close(report.interior_forcing_impulses[a][i], inside[a], 5e-16);
                    close(report.interior_correction_impulses[a][i],
                        report.interior[a][1].momentum[i]
                            - report.interior[a][0].momentum[i] - inside[a], 4e-10);
                    close(report.body_momentum_residuals[a][i], 0.0, 4e-10);
                }
                close(report.combined_momentum_residual[i], 0.0, 5e-10);
            }
            for (i, actual) in force_snapshot(&grid).into_iter().enumerate() {
                let base = history_field([i % 12, (i / 12) % 8, i / 96]);
                assert_eq!(actual, core::array::from_fn(|axis| uniform[axis] + scale * base[axis]));
            }
            for a in 0..2 {
                assert_eq!(bodies[a].mass, initial[a].mass);
                assert_eq!(bodies[a].inertia, initial[a].inertia);
            }
        }
        assert_ne!(bodies[0].center, initial[0].center);
        assert_ne!(bodies[1].angular_momentum_body, initial[1].angular_momentum_body);
    }
}

#[test]
fn zero_and_constant_endpoint_fields_retain_held_dry_assembly_bits() {
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let gate = CancelGate::new();
    for nonzero in [false, true] {
        let make = || {
            let mut grid = grid();
            grid.rebind_force_field(|c| if nonzero { history_field(c) } else { [0.0; 3] });
            grid
        };
        let (mut held, mut endpoint) = (make(), make());
        let (mut a, mut b) = (pair(), pair());
        let (mut sa, mut sb) = (samplers(&held), samplers(&endpoint));
        for _ in 0..2 {
            let external = [[1e-5, -2e-5, 3e-5, 2e-6, -3e-6, 4e-6]; 2];
            let expected = held.step_with_dry_partial_saturation_bodies(
                &mut sa, [&shape, &shape], &mut a, external, [0.0; 3],
                DryBodyCouplingConfig3::default(), &gate,
            ).unwrap();
            let actual = endpoint.step_with_dry_partial_saturation_bodies_and_endpoint_force_scale(
                &mut sb, [&shape, &shape], &mut b, external,
                if nonzero { 1.0 } else { 0.0 }, [0.0; 3],
                DryBodyCouplingConfig3::default(), &gate,
            ).unwrap();
            assert_eq!(actual, expected);
            assert_eq!(a, b);
            assert_eq!(snapshot(&held), snapshot(&endpoint));
            assert_eq!(force_snapshot(&held), force_snapshot(&endpoint));
        }
    }
}

#[test]
fn endpoint_dry_singleton_replays_the_original_single_owner_exactly() {
    let make = || { let mut g = grid(); g.rebind_force_field(history_field); g };
    let (mut one, mut multi) = (make(), make());
    let (mut body, mut bodies) = (pair()[0], [pair()[0]]);
    let (mut field, mut fields) = (BodyFrameSampler3::new(&one, 4).unwrap(), samplers(&multi));
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let gate = CancelGate::new();
    for scale in [0.0, -0.5, 2.0] {
        let a = one.step_with_dry_partial_saturation_body_and_endpoint_force_scale(
            &mut field, &shape, &mut body, [0.0; 6], scale, ORIGIN_FOR_HISTORY,
            DryBodyCouplingConfig3::default(), &gate,
        ).unwrap();
        let b = multi.step_with_dry_partial_saturation_bodies_and_endpoint_force_scale(
            &mut fields, [&shape], &mut bodies, [[0.0; 6]], scale, ORIGIN_FOR_HISTORY,
            DryBodyCouplingConfig3::default(), &gate,
        ).unwrap();
        assert_eq!(body, bodies[0]);
        assert_eq!(snapshot(&one), snapshot(&multi));
        assert_eq!(force_snapshot(&one), force_snapshot(&multi));
        assert_eq!(a.fluid, b.fluid);
        assert_eq!(a.interior, b.interior[0]);
        assert_eq!(a.interior_forcing_impulse, b.interior_forcing_impulses[0]);
        assert_eq!(a.interior_correction_impulse, b.interior_correction_impulses[0]);
        assert_eq!(a.combined_momentum_residual, b.combined_momentum_residual);
        assert_eq!(a.candidate_evaluations, b.candidate_evaluations);
        assert_eq!(a.represented_kinetic_energy_change, b.represented_kinetic_energy_change);
    }
}
const ORIGIN_FOR_HISTORY: [f64; 3] = [-0.7, 0.4, 1.1];

#[test]
fn endpoint_cancellation_after_streaming_preserves_source_and_all_members_for_exact_retry() {
    use std::cell::Cell;
    struct CancellingShape<'a> {
        shape: BodyFrameEllipsoid3,
        gate: &'a CancelGate,
        calls: Cell<usize>,
        cancel_at: Cell<usize>,
    }
    impl BodyFrameShape3 for CancellingShape<'_> {
        fn bounding_radius(&self) -> f64 {
            self.shape.bounding_radius()
        }
        fn contains(&self, p: Vec3) -> bool {
            let calls = self.calls.get() + 1;
            self.calls.set(calls);
            if calls == self.cancel_at.get() {
                self.gate.request();
            }
            self.shape.contains(p)
        }
    }
    let mut grid = grid();
    let mut control = crate::grid();
    grid.rebind_force_field(history_field);
    control.rebind_force_field(history_field);
    let accepted_force = force_snapshot(&grid);
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let gate = CancelGate::new();
    let cancelling = CancellingShape {
        shape,
        gate: &gate,
        calls: Cell::new(0),
        cancel_at: Cell::new(usize::MAX),
    };
    let mut bodies = pair();
    let initial = bodies;
    let mut control_bodies = bodies;
    let mut fields = samplers(&grid);
    let mut control_fields = samplers(&control);
    fields[0]
        .sample(bodies[0].center, bodies[0].orientation, &cancelling, &gate)
        .unwrap();
    let calls = cancelling.calls.get();
    assert!(calls > 0);
    // Initial inventory plus native start/midpoint/endpoint sampling complete
    // before the dry owner resamples the streamed endpoint inventory.
    cancelling.calls.set(0);
    cancelling.cancel_at.set(4 * calls + calls / 2);
    let populations = snapshot(&grid);
    let loads = [[1e-4, 2e-5, -3e-5, 1e-5, -2e-5, 3e-5], [0.0; 6]];
    let error = grid
        .step_with_dry_partial_saturation_bodies_and_endpoint_force_scale(
            &mut fields,
            [&cancelling, &shape],
            &mut bodies,
            loads,
            -0.5,
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap_err();
    assert_eq!(error, PartialSaturationError3::Cancelled);
    assert!(cancelling.calls.get() > 4 * calls);
    assert_eq!(bodies, initial);
    assert_eq!(snapshot(&grid), populations);
    assert_eq!(force_snapshot(&grid), accepted_force);
    grid.voxelize_sdf(|_| 1.0);
    let gate = CancelGate::new();
    let actual = grid
        .step_with_dry_partial_saturation_bodies_and_endpoint_force_scale(
            &mut fields,
            [&shape, &shape],
            &mut bodies,
            loads,
            -0.5,
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    let expected = control
        .step_with_dry_partial_saturation_bodies_and_endpoint_force_scale(
            &mut control_fields,
            [&shape, &shape],
            &mut control_bodies,
            loads,
            -0.5,
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(bodies, control_bodies);
    assert_eq!(snapshot(&grid), snapshot(&control));
}

#[test]
fn endpoint_invalid_scale_and_exhausted_real_candidate_preserve_all_accepted_state() {
    let mut grid = grid();
    grid.rebind_force_field(history_field);
    let mut bodies = pair();
    let original = bodies;
    let populations = snapshot(&grid);
    let force = force_snapshot(&grid);
    let shape = BodyFrameEllipsoid3::sphere(1.2).unwrap();
    let mut fields = samplers(&grid);
    for scale in [f64::NAN, -1.0] {
        let error = grid.step_with_dry_partial_saturation_bodies_and_endpoint_force_scale(
            &mut fields, [&shape, &shape], &mut bodies, [[0.0, 0.0, 1e-4, 0.0, 0.0, 0.0]; 2],
            scale, [0.0; 3],
            DryBodyCouplingConfig3 { maximum_evaluations: 1, ..DryBodyCouplingConfig3::default() },
            &CancelGate::new(),
        ).unwrap_err();
        assert!(matches!(error, PartialSaturationError3::Invalid { reason }
            if reason.contains(if scale.is_nan() { "finite" } else { "budget" })));
        assert_eq!(bodies, original);
        assert_eq!(snapshot(&grid), populations);
        assert_eq!(force_snapshot(&grid), force);
    }
}
