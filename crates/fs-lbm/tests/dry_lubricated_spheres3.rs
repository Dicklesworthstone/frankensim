//! Physical dry inertia composed with the native reciprocal film operator.
use fs_exec::CancelGate;
use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_lbm::d3q19::{
    BodyFrameEllipsoid3, BodyFrameSampler3, BoundaryGrid3, BoundarySpec3, DryBodiesExchange3,
    DryBodyCouplingConfig3, E3, LubricationPlane3, MovingPartialSaturationBody3,
    PartialSaturationError3, Q3, SphereLubrication3,
};
use fs_time::lie::gyrostat::InertiaTensor3;

const RADIUS: f64 = 1.2;
const MU: f64 = (0.83 - 0.5) / 3.0;

fn grid() -> BoundaryGrid3 {
    BoundaryGrid3::new(12, 8, 8, 0.83, [0.0; 3], BoundarySpec3::periodic())
}
fn body(center: [f64; 3], density: f64) -> MovingPartialSaturationBody3 {
    let mass = density * 4.0 * core::f64::consts::PI * RADIUS.powi(3) / 3.0;
    MovingPartialSaturationBody3 {
        center,
        orientation: So3::exp(So3Tangent::new(Vec3::new(0.2, -0.15, 0.1))).unwrap(),
        mass,
        inertia: InertiaTensor3::new(Mat3 {
            m: [0.4, 0.02, -0.01, 0.02, 0.43, 0.015, -0.01, 0.015, 0.46]
                .map(|x| mass * RADIUS.powi(2) * x),
        })
        .unwrap(),
        velocity: [0.0; 3],
        angular_momentum_body: Vec3::new(0.0, 0.0, 0.0),
    }
}
fn pair(density: f64) -> [MovingPartialSaturationBody3; 2] {
    [body([3.5, 4.0, 4.0], density), body([5.91, 4.0, 4.0], 1.7)]
}
fn samplers<const N: usize>(grid: &BoundaryGrid3) -> [BodyFrameSampler3; N] {
    core::array::from_fn(|_| BodyFrameSampler3::new(grid, 4).unwrap())
}
fn model() -> SphereLubrication3 {
    SphereLubrication3::new(MU, 0.25, 0.001, 8).unwrap()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol,
        "{a:.17e} != {b:.17e}, error={:.3e}",
        a - b
    );
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}
fn snapshot(grid: &BoundaryGrid3) -> Vec<[u64; Q3]> {
    let [nx, ny, nz] = grid.dimensions();
    let mut all = Vec::new();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                all.push(grid.populations(x, y, z).map(f64::to_bits));
            }
        }
    }
    all
}
fn sample<const N: usize>(
    fields: &mut [BodyFrameSampler3; N],
    bodies: [MovingPartialSaturationBody3; N],
    midpoint: bool,
) {
    let shape = BodyFrameEllipsoid3::sphere(RADIUS).unwrap();
    for a in 0..N {
        let b = bodies[a];
        let (center, q) = if midpoint {
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
        fields[a]
            .sample(center, q, &shape, &CancelGate::new())
            .unwrap();
    }
}
fn momentum<const N: usize>(
    grid: &BoundaryGrid3,
    fields: &[BodyFrameSampler3; N],
    bodies: [MovingPartialSaturationBody3; N],
    origin: [f64; 3],
) -> [f64; 6] {
    let [nx, ny, nz] = grid.dimensions();
    let mut total = [0.0; 6];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let index = (z * ny + y) * nx + x;
                let exterior = 1.0 - fields.iter().map(|s| s.fractions()[index]).sum::<f64>();
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
                let l = cross(
                    [
                        x as f64 + 0.5 - origin[0],
                        y as f64 + 0.5 - origin[1],
                        z as f64 + 0.5 - origin[2],
                    ],
                    p,
                );
                let first: [f64; 3] = core::array::from_fn(|i| {
                    fields
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
                    total[i + 3] += l[i] - interior_cut[i];
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
    fields: &[BodyFrameSampler3; N],
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
                        inside[a][i] += fields[a].fractions()[(z * ny + y) * nx + x] * value;
                    }
                }
                for a in 0..N {
                    if let Some(&(_, first)) = fields[a]
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
fn kinetic(b: MovingPartialSaturationBody3) -> f64 {
    let spin = b.orientation.rotate(b.angular_momentum_body).unwrap();
    0.5 * b.mass * dot(b.velocity, b.velocity)
        + 0.5
            * dot(
                [spin.x, spin.y, spin.z],
                b.angular_velocity_world().unwrap(),
            )
}

#[test]
fn inactive_films_empty_planes_and_singleton_preserve_original_dry_bits() {
    let shape = BodyFrameEllipsoid3::sphere(RADIUS).unwrap();
    let mut plain = grid();
    let mut filmed = grid();
    let mut with_planes = grid();
    let mut a = pair(0.3);
    a[1].center[0] = 7.0;
    a[0].velocity = [1e-4, -2e-4, 1e-4];
    a[1].velocity = [-2e-4, 1e-4, -1e-4];
    let mut b = a;
    let mut c = a;
    let mut sa = samplers(&plain);
    let mut sb = samplers(&filmed);
    let mut sc = samplers(&with_planes);
    let loads = [[2e-5, -3e-5, 1e-5, 3e-6, -2e-6, 1e-6]; 2];
    let origin = [-0.6, 0.8, 1.3];
    let config = DryBodyCouplingConfig3::default();
    let gate = CancelGate::new();
    for _ in 0..2 {
        let ra = plain
            .step_with_dry_partial_saturation_bodies(
                &mut sa,
                [&shape, &shape],
                &mut a,
                loads,
                origin,
                config,
                &gate,
            )
            .unwrap();
        let rb = filmed
            .step_with_dry_lubricated_spheres(
                &mut sb,
                [RADIUS; 2],
                &mut b,
                loads,
                model(),
                origin,
                config,
                &gate,
            )
            .unwrap();
        let rc = with_planes
            .step_with_dry_lubricated_spheres_and_planes(
                &mut sc,
                [RADIUS; 2],
                &mut c,
                loads,
                model(),
                &[],
                origin,
                config,
                &gate,
            )
            .unwrap();
        assert_eq!(ra, rb);
        assert_eq!(rb, rc.motion);
        assert!(rc.planes.is_empty());
        assert_eq!(a, b);
        assert_eq!(b, c);
        assert_eq!(snapshot(&plain), snapshot(&filmed));
        assert_eq!(snapshot(&filmed), snapshot(&with_planes));
    }
    let mut plain = grid();
    let mut filmed = grid();
    let mut a = pair(0.2)[0];
    let mut b = [a];
    let mut sa = BodyFrameSampler3::new(&plain, 4).unwrap();
    let mut sb = samplers(&filmed);
    let one = plain
        .step_with_dry_partial_saturation_body(
            &mut sa, &shape, &mut a, loads[0], origin, config, &gate,
        )
        .unwrap();
    let many = filmed
        .step_with_dry_lubricated_spheres(
            &mut sb,
            [RADIUS],
            &mut b,
            [loads[0]],
            model(),
            origin,
            config,
            &gate,
        )
        .unwrap();
    assert_eq!(a, b[0]);
    assert_eq!(snapshot(&plain), snapshot(&filmed));
    assert_eq!(one.fluid, many.fluid);
    assert_eq!(
        one.interior_correction_impulse,
        many.interior_correction_impulses[0]
    );
    assert_eq!(one.candidate_evaluations, many.candidate_evaluations);
}

#[test]
fn unequal_density_dry_pair_has_reciprocal_resistance_and_physical_momentum() {
    let origin = [-0.7, 1.1, -0.9];
    let config = DryBodyCouplingConfig3::default();
    let gate = CancelGate::new();
    let shape = BodyFrameEllipsoid3::sphere(RADIUS).unwrap();
    for density in [0.2, 1.0, 3.0] {
        let mut grid = grid();
        let mut uncorrected = crate::grid();
        let mut bodies = pair(density);
        bodies[0].velocity = [1e-4, 0.0, 0.0];
        bodies[1].velocity = [-6e-5, 0.0, 0.0];
        let initial = bodies;
        let mut bare = bodies;
        let mut fields = samplers(&grid);
        let mut bare_fields = samplers(&uncorrected);
        sample(&mut fields, bodies, false);
        let before = momentum(&grid, &fields, bodies, origin);
        let receipt = grid
            .step_with_dry_lubricated_spheres(
                &mut fields,
                [RADIUS; 2],
                &mut bodies,
                [[0.0; 6]; 2],
                model(),
                origin,
                config,
                &gate,
            )
            .unwrap();
        let baseline = uncorrected
            .step_with_dry_partial_saturation_bodies(
                &mut bare_fields,
                [&shape, &shape],
                &mut bare,
                [[0.0; 6]; 2],
                origin,
                config,
                &gate,
            )
            .unwrap();
        let after = momentum(&grid, &fields, bodies, origin);
        let film = receipt.native.coupling.lubrication;
        assert_eq!(film.active_pairs, 1);
        assert_eq!(film.active_resistance_modes, 1);
        assert!(film.dissipated_energy > 0.0);
        assert_eq!(film.plane_impulse, [0.0; 3]);
        let mid: [[f64; 3]; 2] = core::array::from_fn(|a| {
            core::array::from_fn(|i| initial[a].center[i] + 0.5 * initial[a].velocity[i])
        });
        let gap = mid[1][0] - mid[0][0] - 2.0 * RADIUS;
        let reduced = RADIUS / 2.0;
        let zeta = 6.0
            * core::f64::consts::PI
            * MU
            * reduced.powi(2)
            * (1.0 / gap.max(0.001) - 1.0 / 0.25);
        let slip = receipt.native.coupling.bodies[0].velocity_midpoint[0]
            - receipt.native.coupling.bodies[1].velocity_midpoint[0];
        let bare_slip = baseline.native.coupling.bodies[0].velocity_midpoint[0]
            - baseline.native.coupling.bodies[1].velocity_midpoint[0];
        close(film.body_impulses[0][0], -zeta * slip, 2e-12);
        assert!(slip.abs() < bare_slip.abs());
        for i in 0..6 {
            close(after[i] - before[i], boundary(&receipt, i), 5e-10);
            let net = if i < 3 {
                film.body_impulses[0][i] + film.body_impulses[1][i]
            } else {
                film.body_angular_impulses[0][i - 3] + film.body_angular_impulses[1][i - 3]
            };
            close(net, 0.0, 2e-13);
            for a in 0..2 {
                close(receipt.body_momentum_residuals[a][i], 0.0, 4e-10);
            }
        }
        close(
            film.body_work.iter().sum::<f64>() + film.dissipated_energy,
            0.0,
            2e-13,
        );
        for a in 0..2 {
            assert_eq!(bodies[a].mass, initial[a].mass);
            assert_eq!(bodies[a].inertia, initial[a].inertia);
            let native = receipt.native.bodies[a];
            let work = receipt.physical_external_work[a]
                + receipt.interior_correction_work[a]
                + native.coupling.hydrodynamic.body_work
                + film.body_work[a]
                + native.free_rotation_energy_change;
            close(kinetic(bodies[a]) - kinetic(initial[a]), work, 2e-13);
        }
    }
}

#[test]
fn translated_origin_dry_plane_retains_tangent_spin_couples_and_reservoir_work() {
    let mut grid = grid();
    let mut bodies = [body([4.0, 4.0, 4.0], 0.3)];
    bodies[0].velocity = [1e-4, -2e-5, -7e-5];
    bodies[0].angular_momentum_body = Vec3::new(1e-4, -2e-4, 3e-4);
    let initial = bodies;
    let plane = LubricationPlane3::new(
        [0.0, 2.79, 0.0],
        [0.0, 1.0, 0.0],
        [2e-5, 1e-5, -3e-5],
        0.25,
        0.001,
    )
    .unwrap()
    .with_tangential(0.02, 0.001)
    .unwrap();
    let origin = [-8.7, 5.1, 11.3];
    let mut fields = samplers(&grid);
    sample(&mut fields, bodies, false);
    let before = momentum(&grid, &fields, bodies, origin);
    let report = grid
        .step_with_dry_lubricated_spheres_and_planes(
            &mut fields,
            [RADIUS],
            &mut bodies,
            [[0.0; 6]],
            model(),
            &[plane],
            origin,
            DryBodyCouplingConfig3::default(),
            &CancelGate::new(),
        )
        .unwrap();
    let after = momentum(&grid, &fields, bodies, origin);
    let r = report.motion;
    let film = r.native.coupling.lubrication;
    assert_eq!(film.active_plane_pairs, 1);
    assert_eq!(film.active_tangential_plane_pairs, 1);
    assert_eq!(film.active_resistance_modes, 5);
    assert_eq!(report.planes.len(), 1);
    let received = report.planes[0];
    assert!(received.couple_impulse.iter().any(|j| j.abs() > 1e-8));
    let mid: [f64; 3] =
        core::array::from_fn(|i| initial[0].center[i] + 0.5 * initial[0].velocity[i]);
    let plane_mid = plane.translated(0.5).unwrap();
    let foot = [mid[0], plane_mid.point()[1], mid[2]];
    let orbital = cross(
        core::array::from_fn(|i| foot[i] - origin[i]),
        received.impulse,
    );
    for i in 0..3 {
        close(
            received.angular_impulse[i],
            orbital[i] + received.couple_impulse[i],
            3e-13,
        );
        close(film.body_impulses[0][i] + received.impulse[i], 0.0, 3e-13);
        close(
            film.body_angular_impulses[0][i] + received.angular_impulse[i],
            0.0,
            3e-13,
        );
    }
    for i in 0..6 {
        let plane_reaction = if i < 3 {
            received.impulse[i]
        } else {
            received.angular_impulse[i - 3]
        };
        close(
            after[i] - before[i],
            boundary(&r, i) - plane_reaction,
            8e-10,
        );
        close(r.body_momentum_residuals[0][i], 0.0, 8e-10);
        close(r.combined_momentum_residual[i], 0.0, 8e-10);
    }
    close(
        received.work,
        dot(received.impulse, plane.velocity()),
        3e-16,
    );
    close(
        film.body_work[0] + received.work + film.dissipated_energy,
        0.0,
        3e-13,
    );
    assert!(film.dissipated_energy > 0.0);
    let response = r.native.coupling.bodies[0];
    close(
        film.body_work[0],
        dot(film.body_impulses[0], response.velocity_midpoint)
            + dot(
                film.body_spin_impulses[0],
                response.angular_velocity_midpoint,
            ),
        3e-16,
    );
    close(
        kinetic(bodies[0]) - kinetic(initial[0]),
        r.physical_external_work[0]
            + r.interior_correction_work[0]
            + response.hydrodynamic.body_work
            + film.body_work[0]
            + r.native.bodies[0].free_rotation_energy_change,
        3e-13,
    );
}

#[test]
fn shared_pair_and_plane_modes_preserve_physical_forcing_and_distinct_interior_sources() {
    let mut grid = grid();
    grid.rebind_force_field(|[x, y, z]| {
        [x as f64 * 2e-7, y as f64 * -3e-7, (z as f64 - 2.0) * 1e-7]
    });
    let mut bodies = pair(0.2);
    bodies[0].velocity = [5e-5, -2e-5, -3e-5];
    bodies[1].velocity = [-3e-5, 1e-5, 4e-5];
    bodies[0].angular_momentum_body = Vec3::new(1e-4, -2e-4, 3e-4);
    bodies[1].angular_momentum_body = Vec3::new(-2e-4, 1e-4, 4e-4);
    let initial = bodies;
    let model = model().with_tangential_pairs(0.02, 0.001).unwrap();
    let mut plane = LubricationPlane3::new(
        [0.0, 2.79, 0.0],
        [0.0, 1.0, 0.0],
        [2e-5, 1e-5, -3e-5],
        0.25,
        0.001,
    )
    .unwrap()
    .with_tangential(0.02, 0.001)
    .unwrap();
    let loads = [
        [2e-6, -3e-6, 4e-6, 1e-6, -2e-6, 3e-6],
        [-3e-6, 4e-6, -2e-6, 2e-6, 3e-6, -1e-6],
    ];
    let origin = [-4.7, 2.1, 7.3];
    let mut fields = samplers(&grid);
    let mut midpoint = samplers(&grid);
    let gate = CancelGate::new();
    for _ in 0..2 {
        let incoming = bodies;
        sample(&mut fields, bodies, false);
        let before = momentum(&grid, &fields, bodies, origin);
        sample(&mut midpoint, bodies, true);
        let (full, inside) = forcing(&grid, &midpoint, origin);
        let authored: [[f64; 6]; 2] = core::array::from_fn(|a| {
            let arm = core::array::from_fn(|i| {
                bodies[a].center[i] + 0.5 * bodies[a].velocity[i] - origin[i]
            });
            let torque = cross(arm, loads[a][..3].try_into().unwrap());
            core::array::from_fn(|i| loads[a][i] + if i < 3 { 0.0 } else { torque[i - 3] })
        });
        let report = grid
            .step_with_dry_lubricated_spheres_and_planes(
                &mut fields,
                [RADIUS; 2],
                &mut bodies,
                loads,
                model,
                &[plane],
                origin,
                DryBodyCouplingConfig3::default(),
                &gate,
            )
            .unwrap();
        let after = momentum(&grid, &fields, bodies, origin);
        let r = report.motion;
        let film = r.native.coupling.lubrication;
        assert_eq!(film.active_pairs, 1);
        assert_eq!(film.active_tangential_pairs, 1);
        assert_eq!(film.active_plane_pairs, 2);
        assert_eq!(film.active_tangential_plane_pairs, 2);
        assert_eq!(film.active_resistance_modes, 15);
        assert_eq!(report.planes.len(), 1);
        let received = report.planes[0];
        for i in 0..6 {
            let external = authored.iter().map(|j| j[i]).sum::<f64>();
            let exterior = full[i] - inside.iter().map(|j| j[i]).sum::<f64>();
            let reaction = if i < 3 {
                received.impulse[i]
            } else {
                received.angular_impulse[i - 3]
            };
            close(
                after[i] - before[i],
                external + exterior + boundary(&r, i) - reaction,
                8e-10,
            );
            close(r.represented_fluid_forcing_impulse[i], exterior, 2e-14);
            close(r.combined_momentum_residual[i], 0.0, 8e-10);
            let sum_film: f64 = (0..2)
                .map(|a| {
                    if i < 3 {
                        film.body_impulses[a][i]
                    } else {
                        film.body_angular_impulses[a][i - 3]
                    }
                })
                .sum();
            close(sum_film + reaction, 0.0, 3e-13);
            for a in 0..2 {
                close(r.interior_forcing_impulses[a][i], inside[a][i], 2e-15);
                close(
                    r.interior_correction_impulses[a][i],
                    r.interior[a][1].momentum[i] - r.interior[a][0].momentum[i] - inside[a][i],
                    5e-10,
                );
                close(r.body_momentum_residuals[a][i], 0.0, 5e-10);
            }
        }
        for a in 0..2 {
            assert_eq!(bodies[a].mass, initial[a].mass);
            assert_eq!(bodies[a].inertia, initial[a].inertia);
            let native = r.native.bodies[a];
            close(
                kinetic(bodies[a]) - kinetic(incoming[a]),
                r.physical_external_work[a]
                    + r.interior_correction_work[a]
                    + native.coupling.hydrodynamic.body_work
                    + film.body_work[a]
                    + native.free_rotation_energy_change,
                3e-13,
            );
        }
        close(film.plane_work, received.work, 3e-16);
        close(
            received.work,
            dot(received.impulse, plane.velocity()),
            3e-16,
        );
        close(
            film.body_work.iter().sum::<f64>() + received.work + film.dissipated_energy,
            0.0,
            3e-13,
        );
        assert!(film.dissipated_energy > 0.0);
        plane = plane.translated(1.0).unwrap();
    }
}

#[test]
fn film_limits_failed_probes_and_late_plane_contact_preserve_atomic_retry() {
    let mut grid = grid();
    let mut control = crate::grid();
    for grid in [&mut grid, &mut control] {
        grid.rebind_force_field(|[x, y, z]| {
            [x as f64 * 2e-7, y as f64 * -3e-7, (z as f64 - 2.0) * 1e-7]
        });
    }
    let mut bodies = pair(0.3);
    let initial = bodies;
    let mut control_bodies = bodies;
    let mut fields = samplers(&grid);
    let mut control_fields = samplers(&control);
    let plane = LubricationPlane3::new([0.0, 2.79, 0.0], [0.0, 1.0, 0.0], [0.0; 3], 0.25, 0.001)
        .unwrap()
        .with_tangential(0.02, 0.001)
        .unwrap();
    let loads = [[1e-5, -2e-5, 3e-5, 1e-6, -2e-6, 3e-6], [0.0; 6]];
    let model = model().with_tangential_pairs(0.02, 0.001).unwrap();
    let populations = snapshot(&grid);
    let gate = CancelGate::new();
    for (law, config, reason) in [
        (
            model,
            DryBodyCouplingConfig3 {
                maximum_evaluations: 2,
                ..DryBodyCouplingConfig3::default()
            },
            "budget",
        ),
        (
            SphereLubrication3::new(MU, 0.25, 0.001, 1).unwrap(),
            DryBodyCouplingConfig3::default(),
            "active-pair limit",
        ),
        (
            model.with_resistance_mode_limit(2).unwrap(),
            DryBodyCouplingConfig3::default(),
            "resistance-mode limit",
        ),
    ] {
        let error = grid
            .step_with_dry_lubricated_spheres_and_planes(
                &mut fields,
                [RADIUS; 2],
                &mut bodies,
                loads,
                law,
                &[plane],
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
    // Positive initial/midpoint gaps (0.04/0.01), but the prescribed plane
    // crosses the endpoint. All films are outside the tiny cutoff. This refusal
    // occurs in native endpoint admission AFTER its fluid stream is prepared.
    let crossing = LubricationPlane3::new(
        [0.0, 2.76, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.06, 0.0],
        1e-4,
        1e-5,
    )
    .unwrap();
    let error = grid
        .step_with_dry_lubricated_spheres_and_planes(
            &mut fields,
            [RADIUS; 2],
            &mut bodies,
            [[0.0; 6]; 2],
            SphereLubrication3::new(MU, 1e-4, 1e-5, 8).unwrap(),
            &[crossing],
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap_err();
    assert!(
        matches!(error,PartialSaturationError3::Invalid {reason} if reason.contains("meets prescribed lubrication plane")),
        "{error:?}"
    );
    assert_eq!(bodies, initial);
    assert_eq!(snapshot(&grid), populations);
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_eq!(
        grid.step_with_dry_lubricated_spheres_and_planes(
            &mut fields,
            [RADIUS; 2],
            &mut bodies,
            loads,
            model,
            &[plane],
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &cancelled
        ),
        Err(PartialSaturationError3::Cancelled)
    );
    assert_eq!(bodies, initial);
    assert_eq!(snapshot(&grid), populations);
    grid.voxelize_sdf(|_| 1.0);
    let actual = grid
        .step_with_dry_lubricated_spheres_and_planes(
            &mut fields,
            [RADIUS; 2],
            &mut bodies,
            loads,
            model,
            &[plane],
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    let expected = control
        .step_with_dry_lubricated_spheres_and_planes(
            &mut control_fields,
            [RADIUS; 2],
            &mut control_bodies,
            loads,
            model,
            &[plane],
            [0.0; 3],
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    assert_eq!(actual, expected);
    assert_eq!(bodies, control_bodies);
    assert_eq!(snapshot(&grid), snapshot(&control));
}
