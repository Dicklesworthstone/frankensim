//! Shared rigid identity survives actual fluid/body motion and dry correction.
use fs_exec::CancelGate;
use fs_ga::{Mat3, So3, So3Tangent, Vec3};
use fs_lbm::d3q19::{
    BodyFrameSampler3, BoundaryGrid3, BoundarySpec3, DryBodyCouplingConfig3, E3,
    LubricationPlanePatch3, MovingCarrierSurface3, MovingFiniteLubricationPlane3,
    MovingLubricationSurface3, MovingPartialSaturationBody3, Q3, SphereLubrication3,
};
use fs_time::lie::gyrostat::InertiaTensor3;

const RADIUS: f64 = 0.7;
const ORIGIN: [f64; 3] = [-2.0, 0.7, 1.2];
const LOAD: [f64; 6] = [2e-5, -1e-5, 3e-5, 1e-6, -2e-6, 1e-6];
const CARRIER_LOAD: [f64; 6] = [-3e-5, 2e-5, -1e-5, -2e-6, 1e-6, 3e-6];

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn grid() -> BoundaryGrid3 {
    BoundaryGrid3::new(8, 8, 8, 0.83, [0.0; 3], BoundarySpec3::periodic())
}
fn body(mass: f64, velocity: [f64; 3]) -> MovingPartialSaturationBody3 {
    MovingPartialSaturationBody3 {
        center: [3.5, 4.5, 4.5],
        mass,
        velocity,
        orientation: So3::exp(So3Tangent::new(Vec3::new(0.2, -0.15, 0.1))).unwrap(),
        inertia: InertiaTensor3::new(Mat3 {
            m: [0.12, 0.009, -0.004, 0.009, 0.17, 0.006, -0.004, 0.006, 0.21].map(|x| mass * x),
        })
        .unwrap(),
        angular_momentum_body: Vec3::new(2e-6, -3e-6, 4e-6),
    }
}
fn sphere() -> MovingPartialSaturationBody3 {
    body(3.0, [-8e-5, 4e-5, -2e-5])
}
fn carrier() -> MovingPartialSaturationBody3 {
    body(5.0, [2e-5, -1e-5, 3e-5])
}
fn patch(body: MovingPartialSaturationBody3, sign: f64) -> LubricationPlanePatch3 {
    LubricationPlanePatch3::new(
        array(
            body.orientation
                .inverse()
                .rotate(Vec3::new(-sign * 0.71, 0.0, 0.0))
                .unwrap(),
        ),
        array(
            body.orientation
                .inverse()
                .rotate(Vec3::new(sign, 0.0, 0.0))
                .unwrap(),
        ),
        0.03,
        0.001,
    )
    .unwrap()
    .with_tangential(0.013, 0.001)
    .unwrap()
}
fn surfaces(body: MovingPartialSaturationBody3) -> [MovingCarrierSurface3; 2] {
    [1.0, -1.0].map(|sign| MovingCarrierSurface3::Attached {
        carrier: 0,
        patch_body: patch(body, sign),
    })
}
fn model() -> SphereLubrication3 {
    SphereLubrication3::new(0.08, 0.03, 0.001, 8)
        .unwrap()
        .with_resistance_mode_limit(32)
        .unwrap()
}
fn snapshot(g: &BoundaryGrid3) -> Vec<[f64; Q3]> {
    (0..8)
        .flat_map(|z| (0..8).flat_map(move |y| (0..8).map(move |x| g.populations(x, y, z))))
        .collect()
}
fn momentum(b: MovingPartialSaturationBody3) -> [f64; 6] {
    let p = b.velocity.map(|v| b.mass * v);
    let orbit = cross(core::array::from_fn(|i| b.center[i] - ORIGIN[i]), p);
    let spin = array(b.orientation.rotate(b.angular_momentum_body).unwrap());
    core::array::from_fn(|i| {
        if i < 3 {
            p[i]
        } else {
            orbit[i - 3] + spin[i - 3]
        }
    })
}
fn applied(center: [f64; 3], impulse: [f64; 6]) -> [f64; 6] {
    let arm = core::array::from_fn(|i| center[i] - ORIGIN[i]);
    let orbit = cross(arm, impulse[..3].try_into().unwrap());
    core::array::from_fn(|i| {
        if i < 3 {
            impulse[i]
        } else {
            impulse[i] + orbit[i - 3]
        }
    })
}
fn energy(b: MovingPartialSaturationBody3) -> f64 {
    0.5 * b.mass * dot(b.velocity, b.velocity)
        + 0.5
            * dot(
                array(b.angular_momentum_body),
                array(b.inertia.velocity(b.angular_momentum_body).unwrap()),
            )
}

#[test]
fn several_faces_share_one_actual_pose_momentum_and_energy_account() {
    let gate = CancelGate::new();
    let mut g = grid();
    let mut samplers = [BodyFrameSampler3::new(&g, 4).unwrap()];
    let mut bodies = [sphere()];
    let mut carriers = [carrier()];
    let mut patches = surfaces(carriers[0]);
    let local = patches;
    for _ in 0..3 {
        let before = snapshot(&g);
        let old_sphere = bodies[0];
        let old_carrier = carriers[0];
        let report = g
            .step_with_lubricated_spheres_and_carriers(
                &mut samplers,
                [RADIUS],
                &mut bodies,
                [LOAD],
                model(),
                &mut carriers,
                &[CARRIER_LOAD],
                &mut patches,
                ORIGIN,
                &gate,
            )
            .unwrap();
        assert_eq!(report.carriers.len(), 1);
        assert_eq!(report.surfaces.len(), 2);
        assert_eq!(patches, local);
        assert_ne!(carriers[0].center, old_carrier.center);
        assert_ne!(carriers[0].orientation, old_carrier.orientation);
        let account = report.carriers[0];
        for i in 0..3 {
            close(
                account.film.impulse[i],
                report.surfaces.iter().map(|r| r.impulse[i]).sum(),
                2e-18,
            );
            close(
                account.film.angular_impulse[i],
                report.surfaces.iter().map(|r| r.angular_impulse[i]).sum(),
                2e-18,
            );
        }
        close(
            account.film.work,
            report.surfaces.iter().map(|r| r.work).sum(),
            2e-21,
        );
        close(
            account.motion.kinetic_energy_change,
            energy(carriers[0]) - energy(old_carrier),
            3e-21,
        );
        close(account.motion.work_residual, 0.0, 4e-20);
        for &residual in &account.motion.momentum_residual {
            close(residual, 0.0, 4e-16);
        }
        let mut change = [0.0; 6];
        for (cell, (old, new)) in before.iter().zip(snapshot(&g)).enumerate() {
            let dp = core::array::from_fn(|i| {
                (0..Q3)
                    .map(|q| {
                        let e = [E3[q].0, E3[q].1, E3[q].2];
                        f64::from(e[i]) * (new[q] - old[q])
                    })
                    .sum()
            });
            let arm = [
                (cell % 8) as f64 + 0.5 - ORIGIN[0],
                ((cell / 8) % 8) as f64 + 0.5 - ORIGIN[1],
                (cell / 64) as f64 + 0.5 - ORIGIN[2],
            ];
            let dl = cross(arm, dp);
            for i in 0..3 {
                change[i] += dp[i];
                change[i + 3] += dl[i];
            }
        }
        for (old, new, load, center) in [
            (
                old_sphere,
                bodies[0],
                LOAD,
                report.motion.bodies[0].midpoint_center,
            ),
            (
                old_carrier,
                carriers[0],
                CARRIER_LOAD,
                account.motion.midpoint_center,
            ),
        ] {
            let p0 = momentum(old);
            let p1 = momentum(new);
            let source = applied(center, load);
            for i in 0..6 {
                change[i] += p1[i] - p0[i] - source[i];
            }
        }
        for i in 0..3 {
            close(change[i], 0.0, 3e-12);
            close(
                change[i + 3] - report.motion.coupling.periodic.angular_momentum_change[i],
                0.0,
                1e-11,
            );
        }
        for surface in patches {
            let MovingCarrierSurface3::Attached { patch_body, .. } = surface else {
                unreachable!()
            };
            let world = surface.geometry(&carriers).unwrap();
            let arm = array(
                carriers[0]
                    .orientation
                    .rotate(Vec3::new(
                        patch_body.point()[0],
                        patch_body.point()[1],
                        patch_body.point()[2],
                    ))
                    .unwrap(),
            );
            for i in 0..3 {
                close(world.point()[i], carriers[0].center[i] + arm[i], 1e-15);
            }
        }
        close(report.motion.coupling.lubrication.work_residual, 0.0, 4e-20);
        assert!(report.motion.coupling.lubrication.dissipated_energy >= 0.0);
    }
}

#[test]
fn one_patch_and_empty_registry_retain_legacy_dry_results() {
    for attached in [false, true] {
        let mut a = grid();
        let mut b = grid();
        let mut sa = [BodyFrameSampler3::new(&a, 4).unwrap()];
        let mut sb = [BodyFrameSampler3::new(&b, 4).unwrap()];
        let mut ba = [sphere()];
        let mut bb = ba;
        let mut c = carrier();
        c.orientation = So3::identity();
        let p = patch(c, 1.0);
        let legacy = MovingFiniteLubricationPlane3::new(
            c,
            p.point(),
            p.normal(),
            p.cutoff_gap(),
            p.roughness_gap(),
        )
        .unwrap()
        .with_tangential(0.013, 0.001)
        .unwrap();
        let mut walls = if attached {
            vec![MovingLubricationSurface3::Finite(legacy)]
        } else {
            Vec::new()
        };
        let mut carriers = if attached { vec![c] } else { Vec::new() };
        let mut patches = if attached {
            vec![MovingCarrierSurface3::Attached {
                carrier: 0,
                patch_body: p,
            }]
        } else {
            Vec::new()
        };
        let loads = if attached {
            vec![CARRIER_LOAD]
        } else {
            Vec::new()
        };
        let old = a
            .step_with_dry_lubricated_spheres_and_surfaces(
                &mut sa,
                [RADIUS],
                &mut ba,
                [LOAD],
                model(),
                &mut walls,
                &loads,
                ORIGIN,
                DryBodyCouplingConfig3::default(),
                &CancelGate::new(),
            )
            .unwrap();
        let new = b
            .step_with_dry_lubricated_spheres_and_carriers(
                &mut sb,
                [RADIUS],
                &mut bb,
                [LOAD],
                model(),
                &mut carriers,
                &loads,
                &mut patches,
                ORIGIN,
                DryBodyCouplingConfig3::default(),
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(ba, bb);
        assert_eq!(snapshot(&a), snapshot(&b));
        assert_eq!(old.motion, new.motion);
        if attached {
            let MovingLubricationSurface3::Finite(wall) = walls[0] else {
                unreachable!()
            };
            assert_eq!(wall.body(), carriers[0]);
            assert_eq!(old.surfaces[0].finite.unwrap(), new.carriers[0].motion);
            assert_eq!(
                old.surfaces[0].film.finite_body.unwrap(),
                new.carriers[0].film.motion
            );
        } else {
            assert!(new.surfaces.is_empty() && new.carriers.is_empty());
        }
    }
}

#[test]
fn shared_dry_refusal_preserves_every_owner_and_exact_retry() {
    let gate = CancelGate::new();
    let mut g = grid();
    let mut samplers = [BodyFrameSampler3::new(&g, 4).unwrap()];
    let mut bodies = [sphere()];
    let mut carriers = [carrier()];
    let mut patches = surfaces(carriers[0]);
    let initial = bodies;
    let carrier_initial = carriers;
    let patch_initial = patches;
    let fluid_initial = snapshot(&g);
    let limited = DryBodyCouplingConfig3 {
        maximum_evaluations: 1,
        ..DryBodyCouplingConfig3::default()
    };
    assert!(
        g.step_with_dry_lubricated_spheres_and_carriers(
            &mut samplers,
            [RADIUS],
            &mut bodies,
            [LOAD],
            model(),
            &mut carriers,
            &[CARRIER_LOAD],
            &mut patches,
            ORIGIN,
            limited,
            &gate,
        )
        .is_err()
    );
    assert_eq!(bodies, initial);
    assert_eq!(carriers, carrier_initial);
    assert_eq!(patches, patch_initial);
    assert_eq!(snapshot(&g), fluid_initial);
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        g.step_with_dry_lubricated_spheres_and_carriers(
            &mut samplers,
            [RADIUS],
            &mut bodies,
            [LOAD],
            model(),
            &mut carriers,
            &[CARRIER_LOAD],
            &mut patches,
            ORIGIN,
            DryBodyCouplingConfig3::default(),
            &cancelled,
        )
        .is_err()
    );
    assert_eq!(bodies, initial);
    assert_eq!(carriers, carrier_initial);
    assert_eq!(patches, patch_initial);
    assert_eq!(snapshot(&g), fluid_initial);
    let report = g
        .step_with_dry_lubricated_spheres_and_carriers(
            &mut samplers,
            [RADIUS],
            &mut bodies,
            [LOAD],
            model(),
            &mut carriers,
            &[CARRIER_LOAD],
            &mut patches,
            ORIGIN,
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    assert!(report.motion.iterations > 0);
    let mut fresh = grid();
    let mut fields = [BodyFrameSampler3::new(&fresh, 4).unwrap()];
    let mut b2 = initial;
    let mut c2 = carrier_initial;
    let mut p2 = patch_initial;
    let repeated = fresh
        .step_with_dry_lubricated_spheres_and_carriers(
            &mut fields,
            [RADIUS],
            &mut b2,
            [LOAD],
            model(),
            &mut c2,
            &[CARRIER_LOAD],
            &mut p2,
            ORIGIN,
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    assert_eq!(report, repeated);
    assert_eq!(bodies, b2);
    assert_eq!(carriers, c2);
    assert_eq!(patches, p2);
    assert_eq!(snapshot(&g), snapshot(&fresh));
}

fn endpoint_field([x, y, z]: [usize; 3]) -> [f64; 3] {
    [
        (y as f64 - 2.3) * 2e-6,
        (z as f64 + 0.7) * -1e-6,
        (x as f64 - 1.1) * 3e-6,
    ]
}
fn force_snapshot(g: &BoundaryGrid3) -> Vec<[f64; 3]> {
    (0..512)
        .map(|i| g.force_density(i % 8, (i / 8) % 8, i / 64))
        .collect()
}
fn sample_body(
    sampler: &mut BodyFrameSampler3,
    body: MovingPartialSaturationBody3,
    midpoint: bool,
) {
    let shape = fs_lbm::d3q19::BodyFrameEllipsoid3::sphere(RADIUS).unwrap();
    let (center, orientation) = if midpoint {
        let (q, _) = body
            .inertia
            .momentum_step(
                body.orientation,
                body.angular_momentum_body,
                Vec3::default(),
                0.5,
            )
            .unwrap();
        (
            core::array::from_fn(|i| body.center[i] + 0.5 * body.velocity[i]),
            q,
        )
    } else {
        (body.center, body.orientation)
    };
    sampler
        .sample(center, orientation, &shape, &CancelGate::new())
        .unwrap();
}
fn field_inventory(
    g: &BoundaryGrid3,
    sampler: &BodyFrameSampler3,
    force: bool,
) -> ([f64; 6], [f64; 6]) {
    let mut full = [0.0; 6];
    let mut inside = [0.0; 6];
    for cell in 0..512 {
        let [x, y, z] = [cell % 8, (cell / 8) % 8, cell / 64];
        let f = g.populations(x, y, z);
        let source = g.force_density(x, y, z);
        let p: [f64; 3] = if force {
            source
        } else {
            core::array::from_fn(|i| {
                (0..Q3)
                    .map(|q| {
                        let e = [E3[q].0, E3[q].1, E3[q].2];
                        f64::from(e[i]) * f[q]
                    })
                    .sum::<f64>()
                    + 0.5 * source[i]
            })
        };
        let arm = [
            x as f64 + 0.5 - ORIGIN[0],
            y as f64 + 0.5 - ORIGIN[1],
            z as f64 + 0.5 - ORIGIN[2],
        ];
        let angular = cross(arm, p);
        let epsilon = sampler.fractions()[cell];
        let first = sampler
            .cut_cell_first_moments()
            .iter()
            .find(|&&(index, _)| index == cell)
            .map_or([0.0; 3], |&(_, m)| m);
        let cut = cross(first, p);
        for i in 0..3 {
            full[i] += p[i];
            full[i + 3] += angular[i];
            inside[i] += epsilon * p[i];
            inside[i + 3] += epsilon * angular[i] + cut[i];
        }
    }
    (full, inside)
}
fn physical_inventory(
    g: &BoundaryGrid3,
    sampler: &BodyFrameSampler3,
    sphere: MovingPartialSaturationBody3,
    carriers: &[MovingPartialSaturationBody3],
) -> [f64; 6] {
    let (all, inside) = field_inventory(g, sampler, false);
    let body = momentum(sphere);
    core::array::from_fn(|i| {
        all[i] - inside[i]
            + body[i]
            + carriers
                .iter()
                .map(|&carrier| momentum(carrier)[i])
                .sum::<f64>()
    })
}

#[test]
fn endpoint_dry_history_closes_exterior_fluid_sphere_and_unique_carrier_inventories() {
    let uniform = [2e-6, -1e-6, 3e-6];
    let carrier_loads = [CARRIER_LOAD, [2e-5, -3e-5, 1e-5, 3e-6, 1e-6, -2e-6]];
    for boundaries in [
        BoundarySpec3::periodic(),
        BoundarySpec3::new([fs_lbm::d3q19::FaceBoundary3::stationary_wall(); 6]),
        BoundarySpec3::velocity_pressure_duct_z([0.001, -0.002, 0.012], 1.002),
    ] {
        let mut g = BoundaryGrid3::new(8, 8, 8, 0.83, uniform, boundaries);
        g.rebind_force_field(endpoint_field);
        let mut samplers = [BodyFrameSampler3::new(&g, 4).unwrap()];
        let mut independent = BodyFrameSampler3::new(&g, 4).unwrap();
        let mut midpoint = BodyFrameSampler3::new(&g, 4).unwrap();
        let mut bodies = [sphere()];
        let mut other = body(7.0, [-1e-5, 3e-5, -2e-5]);
        other.center = [3.3, 3.6, 4.8];
        let mut carriers = [carrier(), other];
        let bottom = LubricationPlanePatch3::new(
            array(
                other
                    .orientation
                    .inverse()
                    .rotate(Vec3::new(0.2, 0.19, -0.3))
                    .unwrap(),
            ),
            array(
                other
                    .orientation
                    .inverse()
                    .rotate(Vec3::new(0.0, 1.0, 0.0))
                    .unwrap(),
            ),
            0.03,
            0.001,
        )
        .unwrap()
        .with_tangential(0.013, 0.001)
        .unwrap();
        let prescribed = fs_lbm::d3q19::LubricationPlane3::new(
            [3.5, 4.5, 5.21],
            [0.0, 0.0, -1.0],
            [1e-5, -2e-5, 1e-5],
            0.03,
            0.001,
        )
        .unwrap()
        .with_tangential(0.013, 0.001)
        .unwrap()
        .with_axial_spin(2e-5)
        .unwrap();
        let pair = surfaces(carriers[0]);
        let mut patches = [
            pair[0],
            pair[1],
            MovingCarrierSurface3::Attached {
                carrier: 1,
                patch_body: bottom,
            },
            MovingCarrierSurface3::Prescribed(prescribed),
        ];
        let initial_carriers = carriers;
        for scale in [0.0, -0.5, 1.5] {
            let incoming = bodies[0];
            let incoming_carriers = carriers;
            let incoming_patches = patches;
            sample_body(&mut independent, incoming, false);
            let before = physical_inventory(&g, &independent, incoming, &carriers);
            sample_body(&mut midpoint, incoming, true);
            let (old_full, old_inside) = field_inventory(&g, &midpoint, true);
            let report = g
                .step_with_dry_lubricated_spheres_and_carriers_and_endpoint_force_scale(
                    &mut samplers,
                    [RADIUS],
                    &mut bodies,
                    [LOAD],
                    model(),
                    &mut carriers,
                    &carrier_loads,
                    &mut patches,
                    scale,
                    ORIGIN,
                    DryBodyCouplingConfig3::default(),
                    &CancelGate::new(),
                )
                .unwrap();
            sample_body(&mut independent, bodies[0], false);
            let after = physical_inventory(&g, &independent, bodies[0], &carriers);
            let (new_full, new_inside) = field_inventory(&g, &midpoint, true);
            assert_eq!(report.carriers.len(), 2);
            assert_eq!(report.surfaces.len(), 4);
            assert_eq!(
                report.motion.native.coupling.lubrication.active_plane_pairs,
                4
            );
            assert_eq!(
                report
                    .motion
                    .native
                    .coupling
                    .lubrication
                    .active_resistance_modes,
                20
            );
            assert!(report.motion.candidate_evaluations > 1);
            let mut authored = applied(
                core::array::from_fn(|i| incoming.center[i] + 0.5 * incoming.velocity[i]),
                LOAD,
            );
            for a in 0..2 {
                let load = applied(
                    core::array::from_fn(|i| {
                        incoming_carriers[a].center[i] + 0.5 * incoming_carriers[a].velocity[i]
                    }),
                    carrier_loads[a],
                );
                for i in 0..6 {
                    authored[i] += load[i];
                }
                let old = momentum(incoming_carriers[a]);
                let new = momentum(carriers[a]);
                let account = report.carriers[a];
                for i in 0..6 {
                    let faces: f64 = report
                        .surfaces
                        .iter()
                        .zip(incoming_patches)
                        .filter_map(|(r, patch)| match patch {
                            MovingCarrierSurface3::Attached { carrier, .. } if carrier == a => {
                                Some(if i < 3 {
                                    r.impulse[i]
                                } else {
                                    r.angular_impulse[i - 3]
                                })
                            }
                            _ => None,
                        })
                        .sum();
                    let aggregate = if i < 3 {
                        account.film.impulse[i]
                    } else {
                        account.film.angular_impulse[i - 3]
                    };
                    close(aggregate, faces, 3e-17);
                    close(new[i] - old[i], load[i] + faces, 5e-16);
                    close(account.motion.momentum_change[i], new[i] - old[i], 5e-18);
                }
                close(
                    energy(carriers[a]) - energy(incoming_carriers[a]),
                    account.film.work
                        + account.film.motion.external_work
                        + account.motion.free_rotation_energy_change,
                    2e-19,
                );
                assert_eq!(carriers[a].mass, initial_carriers[a].mass);
                assert_eq!(carriers[a].inertia, initial_carriers[a].inertia);
            }
            let native = report.motion.native.coupling;
            let prescribed = report.surfaces[3];
            for i in 0..6 {
                let full = 0.5 * old_full[i] + 0.5 * new_full[i];
                let inside = 0.5 * old_inside[i] + 0.5 * new_inside[i];
                let source = if i < 3 {
                    native.forcing.impulse[i]
                } else {
                    native.forcing.angular_impulse[i - 3]
                };
                let boundary = if i < 3 {
                    native.walls.fluid_population_impulse[i] + native.open.impulse[i]
                        - prescribed.impulse[i]
                } else {
                    native.walls.fluid_population_angular_impulse[i - 3]
                        + native.periodic.angular_momentum_change[i - 3]
                        + native.open.angular_impulse[i - 3]
                        - prescribed.angular_impulse[i - 3]
                };
                close(
                    after[i] - before[i],
                    authored[i] + full - inside + boundary,
                    5e-10,
                );
                close(source, full, 2e-15);
                close(report.motion.interior_forcing_impulses[0][i], inside, 3e-16);
                close(
                    report.motion.represented_fluid_forcing_impulse[i],
                    full - inside,
                    2e-15,
                );
                close(
                    report.motion.interior_correction_impulses[0][i],
                    report.motion.interior[0][1].momentum[i]
                        - report.motion.interior[0][0].momentum[i]
                        - inside,
                    3e-10,
                );
            }
            for (i, actual) in force_snapshot(&g).into_iter().enumerate() {
                let shape = endpoint_field([i % 8, (i / 8) % 8, i / 64]);
                assert_eq!(
                    actual,
                    core::array::from_fn(|axis| uniform[axis] + scale * shape[axis])
                );
            }
            for i in 0..3 {
                assert_eq!(patches[i], incoming_patches[i]);
            }
            let MovingCarrierSurface3::Prescribed(old) = incoming_patches[3] else {
                unreachable!()
            };
            let MovingCarrierSurface3::Prescribed(new) = patches[3] else {
                unreachable!()
            };
            assert_eq!(new, old.translated(1.0).unwrap());
        }
        assert_ne!(carriers[0].orientation, initial_carriers[0].orientation);
        assert_ne!(carriers[1].velocity, initial_carriers[1].velocity);
    }
}

#[test]
fn constant_and_zero_endpoint_fields_keep_shared_carrier_transaction_bits() {
    for dry in [false, true] {
        for nonzero in [false, true] {
            let make = || {
                let mut g = grid();
                g.rebind_force_field(|p| if nonzero { endpoint_field(p) } else { [0.0; 3] });
                g
            };
            let (mut held, mut endpoint) = (make(), make());
            let (mut a, mut b) = ([sphere()], [sphere()]);
            let (mut ca, mut cb) = ([carrier()], [carrier()]);
            let (mut pa, mut pb) = (surfaces(ca[0]), surfaces(cb[0]));
            let (mut fa, mut fb) = (
                [BodyFrameSampler3::new(&held, 4).unwrap()],
                [BodyFrameSampler3::new(&endpoint, 4).unwrap()],
            );
            let scale = if nonzero { 1.0 } else { 0.0 };
            let gate = CancelGate::new();
            if dry {
                let expected = held
                    .step_with_dry_lubricated_spheres_and_carriers(
                        &mut fa,
                        [RADIUS],
                        &mut a,
                        [LOAD],
                        model(),
                        &mut ca,
                        &[CARRIER_LOAD],
                        &mut pa,
                        ORIGIN,
                        DryBodyCouplingConfig3::default(),
                        &gate,
                    )
                    .unwrap();
                let actual = endpoint
                    .step_with_dry_lubricated_spheres_and_carriers_and_endpoint_force_scale(
                        &mut fb,
                        [RADIUS],
                        &mut b,
                        [LOAD],
                        model(),
                        &mut cb,
                        &[CARRIER_LOAD],
                        &mut pb,
                        scale,
                        ORIGIN,
                        DryBodyCouplingConfig3::default(),
                        &gate,
                    )
                    .unwrap();
                assert_eq!(actual, expected);
            } else {
                let expected = held
                    .step_with_lubricated_spheres_and_carriers(
                        &mut fa,
                        [RADIUS],
                        &mut a,
                        [LOAD],
                        model(),
                        &mut ca,
                        &[CARRIER_LOAD],
                        &mut pa,
                        ORIGIN,
                        &gate,
                    )
                    .unwrap();
                let actual = endpoint
                    .step_with_lubricated_spheres_and_carriers_and_endpoint_force_scale(
                        &mut fb,
                        [RADIUS],
                        &mut b,
                        [LOAD],
                        model(),
                        &mut cb,
                        &[CARRIER_LOAD],
                        &mut pb,
                        scale,
                        ORIGIN,
                        &gate,
                    )
                    .unwrap();
                assert_eq!(actual, expected);
            }
            assert_eq!(a, b);
            assert_eq!(ca, cb);
            assert_eq!(pa, pb);
            assert_eq!(snapshot(&held), snapshot(&endpoint));
            assert_eq!(force_snapshot(&held), force_snapshot(&endpoint));
        }
    }
}

#[test]
fn endpoint_retained_shared_faces_close_actual_full_fluid_and_one_carrier_momentum() {
    let mut g = grid();
    g.rebind_force_field(endpoint_field);
    let mut samplers = [BodyFrameSampler3::new(&g, 4).unwrap()];
    let mut bodies = [sphere()];
    let mut carriers = [carrier()];
    let mut patches = surfaces(carriers[0]);
    let gate = CancelGate::new();
    for scale in [0.0, -0.5, 1.5] {
        let old_sphere = bodies[0];
        let old_carrier = carriers[0];
        let old_fluid = field_inventory(&g, &samplers[0], false).0;
        let old_source = field_inventory(&g, &samplers[0], true).0;
        let report = g
            .step_with_lubricated_spheres_and_carriers_and_endpoint_force_scale(
                &mut samplers,
                [RADIUS],
                &mut bodies,
                [LOAD],
                model(),
                &mut carriers,
                &[CARRIER_LOAD],
                &mut patches,
                scale,
                ORIGIN,
                &gate,
            )
            .unwrap();
        let new_fluid = field_inventory(&g, &samplers[0], false).0;
        let new_source = field_inventory(&g, &samplers[0], true).0;
        let p0 = momentum(old_sphere);
        let p1 = momentum(bodies[0]);
        let c0 = momentum(old_carrier);
        let c1 = momentum(carriers[0]);
        let applied_sphere = applied(
            core::array::from_fn(|i| old_sphere.center[i] + 0.5 * old_sphere.velocity[i]),
            LOAD,
        );
        let applied_carrier = applied(
            core::array::from_fn(|i| old_carrier.center[i] + 0.5 * old_carrier.velocity[i]),
            CARRIER_LOAD,
        );
        for i in 0..6 {
            let source = 0.5 * old_source[i] + 0.5 * new_source[i];
            let flux = if i < 3 {
                0.0
            } else {
                report.motion.coupling.periodic.angular_momentum_change[i - 3]
            };
            close(
                new_fluid[i] - old_fluid[i] + p1[i] - p0[i] + c1[i] - c0[i],
                source + applied_sphere[i] + applied_carrier[i] + flux,
                1e-11,
            );
            let faces: f64 = report
                .surfaces
                .iter()
                .map(|r| {
                    if i < 3 {
                        r.impulse[i]
                    } else {
                        r.angular_impulse[i - 3]
                    }
                })
                .sum();
            close(c1[i] - c0[i], applied_carrier[i] + faces, 4e-16);
        }
        assert_eq!(report.carriers.len(), 1);
        assert_eq!(report.motion.coupling.lubrication.active_plane_pairs, 2);
        for (i, force) in force_snapshot(&g).into_iter().enumerate() {
            assert_eq!(
                force,
                endpoint_field([i % 8, (i / 8) % 8, i / 64]).map(|f| scale * f)
            );
        }
    }
}

#[test]
fn endpoint_shared_dry_refusal_preserves_source_registry_and_exact_retry() {
    let gate = CancelGate::new();
    let mut g = grid();
    g.rebind_force_field(endpoint_field);
    let initial_force = force_snapshot(&g);
    let mut samplers = [BodyFrameSampler3::new(&g, 4).unwrap()];
    let mut bodies = [sphere()];
    let mut carriers = [carrier()];
    let mut patches = surfaces(carriers[0]);
    let initial = bodies;
    let carrier_initial = carriers;
    let patch_initial = patches;
    let fluid_initial = snapshot(&g);
    let limited = DryBodyCouplingConfig3 {
        maximum_evaluations: 1,
        ..DryBodyCouplingConfig3::default()
    };
    assert!(
        g.step_with_dry_lubricated_spheres_and_carriers_and_endpoint_force_scale(
            &mut samplers,
            [RADIUS],
            &mut bodies,
            [LOAD],
            model(),
            &mut carriers,
            &[CARRIER_LOAD],
            &mut patches,
            -0.5,
            ORIGIN,
            limited,
            &gate,
        )
        .is_err()
    );
    assert_eq!(bodies, initial);
    assert_eq!(carriers, carrier_initial);
    assert_eq!(patches, patch_initial);
    assert_eq!(snapshot(&g), fluid_initial);
    assert_eq!(force_snapshot(&g), initial_force);
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        g.step_with_dry_lubricated_spheres_and_carriers_and_endpoint_force_scale(
            &mut samplers,
            [RADIUS],
            &mut bodies,
            [LOAD],
            model(),
            &mut carriers,
            &[CARRIER_LOAD],
            &mut patches,
            -0.5,
            ORIGIN,
            DryBodyCouplingConfig3::default(),
            &cancelled,
        )
        .is_err()
    );
    assert_eq!(bodies, initial);
    assert_eq!(carriers, carrier_initial);
    assert_eq!(patches, patch_initial);
    assert_eq!(snapshot(&g), fluid_initial);
    assert_eq!(force_snapshot(&g), initial_force);
    let report = g
        .step_with_dry_lubricated_spheres_and_carriers_and_endpoint_force_scale(
            &mut samplers,
            [RADIUS],
            &mut bodies,
            [LOAD],
            model(),
            &mut carriers,
            &[CARRIER_LOAD],
            &mut patches,
            -0.5,
            ORIGIN,
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    assert!(report.motion.iterations > 0);
    let mut fresh = grid();
    fresh.rebind_force_field(endpoint_field);
    let mut fields = [BodyFrameSampler3::new(&fresh, 4).unwrap()];
    let mut b2 = initial;
    let mut c2 = carrier_initial;
    let mut p2 = patch_initial;
    let repeated = fresh
        .step_with_dry_lubricated_spheres_and_carriers_and_endpoint_force_scale(
            &mut fields,
            [RADIUS],
            &mut b2,
            [LOAD],
            model(),
            &mut c2,
            &[CARRIER_LOAD],
            &mut p2,
            -0.5,
            ORIGIN,
            DryBodyCouplingConfig3::default(),
            &gate,
        )
        .unwrap();
    assert_eq!(report, repeated);
    assert_eq!(bodies, b2);
    assert_eq!(carriers, c2);
    assert_eq!(patches, p2);
    assert_eq!(snapshot(&g), snapshot(&fresh));
}
