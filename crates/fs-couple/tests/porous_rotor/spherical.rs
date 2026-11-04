//! Reuse the original population/inventory fixture for unrestricted rotation.
use super::*;
use fs_couple::porous_momentum::spherical::SphericalPorousBody;
use fs_mbd::{Gravity, RigidBodyIntegrator, Wrench};

fn radial(r: f64) -> PorousRotorResistance {
    PorousRotorResistance {
        darcy: 0.12 * (1.0 + r),
        forchheimer: 4.0 + 2.0 * r,
    }
}
fn gyro(orientation: UnitQuaternion) -> (RigidBodyState, MassProperties) {
    let principal = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -1.0, 2.0), 0.47).unwrap();
    let properties =
        MassProperties::from_principal_frame(16.0, Vec3::ZERO, Vec3::new(1.2, 1.6, 2.0), principal)
            .unwrap();
    let spin = properties
        .angular_momentum_body_checked(Vec3::new(0.015, -0.02, 0.06))
        .unwrap();
    (
        RigidBodyState::new(Pose::new(COM, orientation).unwrap(), Vec3::ZERO, spin).unwrap(),
        properties,
    )
}
fn spinning_grid() -> BoundaryGrid3 {
    let mut grid = grid(false);
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let arm = position([x, y, z]).sub(COM);
                let u = Vec3::new(0.025, -0.015, 0.01)
                    .cross(arm)
                    .add(Vec3::new(0.004, -0.003, 0.005));
                let rho = 1.0 + 0.002 * x as f64 - 0.001 * z as f64;
                grid.set_populations(
                    x,
                    y,
                    z,
                    &equilibrium3(rho, a(u.scale(units().timestep_s / units().cell_size_m))),
                );
            }
        }
    }
    grid
}
fn spherical(nonuniform: bool) -> SphericalPorousBody {
    let (body, properties) = gyro(UnitQuaternion::IDENTITY);
    SphericalPorousBody::new(
        if nonuniform {
            spinning_grid()
        } else {
            grid(false)
        },
        body,
        properties,
        units(),
        cylinder().radius_m,
        CELLS,
        &CancelGate::new(),
    )
    .unwrap()
}
fn world_spin(body: RigidBodyState) -> Vec3 {
    body.pose()
        .orientation()
        .rotate_body_to_world(body.angular_momentum_body())
}

#[test]
fn torque_free_asymmetric_body_retains_the_original_gyro_and_spatial_invariants() {
    let mut motion = spherical(false);
    let initial = motion.body();
    let angular = world_spin(initial);
    let energy = initial.kinetic_energy(motion.properties()).unwrap();
    for _ in 0..4 {
        let expected = RigidBodyIntegrator::new(Gravity::ZERO)
            .step_cayley(
                motion.body(),
                motion.properties(),
                Wrench {
                    force_world: Vec3::ZERO,
                    torque_body: Vec3::ZERO,
                },
                units().timestep_s,
            )
            .unwrap();
        let step = motion
            .step_candidate(
                |_| PorousRotorResistance::default(),
                |_| [0.0; 3],
                Vec3::ZERO,
                controls(),
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(step.mechanics, expected);
        close_v(step.guide_impulse_kg_m_s, Vec3::ZERO, 0.0);
        close(step.relative_dissipation_j, 0.0, 0.0);
        motion = step.motion;
        close_v(world_spin(motion.body()), angular, 5e-14);
        close(
            motion.body().kinetic_energy(motion.properties()).unwrap(),
            energy,
            5e-15,
        );
    }
    assert!(
        norm(
            motion
                .body()
                .angular_momentum_body()
                .sub(initial.angular_momentum_body())
        ) > 1e-4
    );
    assert_ne!(
        motion.body().pose().orientation(),
        initial.pose().orientation()
    );
    assert_eq!(motion.accepted_steps(), 4);
    close(motion.time_s(), 2.0, 0.0);
}

#[test]
fn nonlinear_off_axis_fluid_torque_reaches_all_native_spin_components_and_work() {
    let motion = spherical(true);
    let c = controls();
    let external = Vec3::new(0.004, -0.003, 0.002);
    let mut samples = 0;
    let mut sources = 0;
    let step = motion
        .step_candidate(
            |r| {
                samples += 1;
                radial(r)
            },
            |cell| {
                sources += 1;
                forcing(cell)
            },
            external,
            c,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(samples, 56);
    assert_eq!(sources, CELLS);
    assert!(step.native_trials > 2);
    let u = units();
    let speed = u.cell_size_m / u.timestep_s;
    let jscale = u.mass_per_density_cell_kg * speed;
    let lscale = jscale * u.cell_size_m;
    let mut drag = Vec3::ZERO;
    let mut drag_moment = Vec3::ZERO;
    let mut loss = 0.0;
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let arm = position([x, y, z]).sub(COM);
                let radius = norm(arm);
                if radius > motion.radius_m() {
                    continue;
                }
                let law = radial(radius);
                let imposed = step
                    .angular_velocity_world_rad_s
                    .cross(arm)
                    .scale(1.0 / speed);
                let (rho0, p0) = cell(motion.fluid(), x, y, z);
                let (rho1, p1) = cell(step.motion.fluid(), x, y, z);
                let s0 = p0.scale(1.0 / rho0).sub(imposed);
                let s1 = p1.scale(1.0 / rho1).sub(imposed);
                let j = s0
                    .scale(law.darcy + law.forchheimer * norm(s0))
                    .add(s1.scale(law.darcy + law.forchheimer * norm(s1)))
                    .scale(0.5 * jscale);
                drag = drag.add(j);
                drag_moment = drag_moment.add(arm.cross(j));
                loss += j.dot(s0.add(s1).scale(0.5 * speed));
            }
        }
    }
    close_v(v(step.exchange.skeleton_impulse).scale(jscale), drag, 3e-12);
    close_v(
        v(step.exchange.skeleton_angular_impulse).scale(lscale),
        drag_moment,
        3e-12,
    );
    close(step.relative_dissipation_j, loss, 3e-12);
    close_v(step.guide_impulse_kg_m_s, drag.scale(-1.0), 3e-12);
    let change = world_spin(step.motion.body()).sub(world_spin(motion.body()));
    close_v(change, drag_moment.add(external.scale(u.timestep_s)), 3e-10);
    assert!(norm(drag_moment.cross(world_spin(motion.body()))) > 1e-5);
    assert!(change.x.abs() > 1e-4 && change.y.abs() > 1e-4 && change.z.abs() > 1e-4);
    let kinetic = step
        .motion
        .body()
        .kinetic_energy(motion.properties())
        .unwrap()
        - motion.body().kinetic_energy(motion.properties()).unwrap();
    close(
        kinetic + step.fluid_drag_work_j + loss
            - step.angular_velocity_world_rad_s.dot(external) * u.timestep_s,
        0.0,
        3e-12,
    );
    let initial_omega = motion.angular_velocity_world_rad_s().unwrap();
    let final_omega = step.motion.angular_velocity_world_rad_s().unwrap();
    // This distinction matters for a tumbling body: averaging endpoint world
    // velocities would use the wrong work-conjugate source field.
    assert!(
        norm(
            initial_omega
                .add(final_omega)
                .scale(0.5)
                .sub(step.angular_velocity_world_rad_s)
        ) > 1e-7
    );
    let f0 = physical(motion.fluid(), c.flow.balance_origin_m);
    let f1 = physical(step.motion.fluid(), c.flow.balance_origin_m);
    let force = v(step.exchange.external_impulse).scale(jscale);
    let wall = v(step.exchange.walls.fluid_population_impulse).scale(jscale);
    let shift =
        |moment: Vec3, impulse: Vec3| moment.add(COM.sub(c.flow.balance_origin_m).cross(impulse));
    close_v(f1.0.sub(f0.0), force.add(wall).sub(drag), 3e-10);
    let supplied = shift(
        v(step.exchange.external_angular_impulse).scale(lscale),
        force,
    )
    .add(shift(
        v(step.exchange.walls.fluid_population_angular_impulse).scale(lscale),
        wall,
    ))
    .add(COM.sub(c.flow.balance_origin_m).cross(drag.scale(-1.0)))
    .add(external.scale(u.timestep_s));
    close_v(f1.1.sub(f0.1).add(change), supplied, 3e-10);
    close(
        step.other_fluid_kinetic_change_j,
        f1.2 - f0.2 - step.fluid_external_work_j - step.fluid_drag_work_j,
        3e-12,
    );
    let mut replay = motion.fluid().clone();
    let links = replay.boundary_links();
    let receipt = replay
        .step_with_darcy_forchheimer_and_wall_momentum(
            |cell| {
                let arm = position(cell).sub(COM);
                let inside = norm(arm) <= motion.radius_m();
                let law = if inside {
                    radial(norm(arm))
                } else {
                    PorousRotorResistance::default()
                };
                DarcyForchheimer3 {
                    linear: LinearResistance3 {
                        principal: [law.darcy; 3],
                        skeleton_velocity: if inside {
                            a(step
                                .angular_velocity_world_rad_s
                                .cross(arm.scale(1.0 / speed)))
                        } else {
                            [0.0; 3]
                        },
                        body_force: forcing(cell),
                        ..Default::default()
                    },
                    quadratic: law.forchheimer,
                }
            },
            &links,
            a(COM.sub(u.grid_origin_m).scale(1.0 / u.cell_size_m)),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(snapshot(&replay), snapshot(step.motion.fluid()));
    assert_eq!(receipt, step.exchange);
    let shifted = motion
        .step_candidate(
            radial,
            forcing,
            external,
            PorousRotorControls {
                flow: PorousMomentumControls {
                    balance_origin_m: Vec3::new(2.0, -1.0, 0.5),
                    ..c.flow
                },
                ..c
            },
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(
        snapshot(shifted.motion.fluid()),
        snapshot(step.motion.fluid())
    );
    assert_eq!(shifted.motion.body(), step.motion.body());
}

#[test]
fn cubic_frame_rotation_transforms_flow_full_inertia_and_external_couple_together() {
    // A proper cubic rotation is exact for D3Q19; no arbitrary continuous-grid
    // rotational invariance is assumed by this comparison.
    let rotate = |p: Vec3| Vec3::new(p.z, p.x, p.y);
    let translation = Vec3::new(2.0, -1.0, 0.5);
    let motion = spherical(true);
    let c = controls();
    let external = Vec3::new(0.004, -0.003, 0.002);
    let base = motion
        .step_candidate(radial, forcing, external, c, &CancelGate::new())
        .unwrap();
    let mut mapped = grid(false);
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let mut f = [0.0; Q3];
                for (q, e) in E3.into_iter().enumerate() {
                    let target = E3
                        .iter()
                        .position(|t| (t.0, t.1, t.2) == (e.2, e.0, e.1))
                        .unwrap();
                    f[target] = motion.fluid().populations(x, y, z)[q];
                }
                mapped.set_populations(z, x, y, &f);
            }
        }
    }
    let orientation = UnitQuaternion::new(0.5, 0.5, 0.5, 0.5).unwrap();
    let moved = RigidBodyState::new(
        Pose::new(rotate(COM).add(translation), orientation).unwrap(),
        Vec3::ZERO,
        motion.body().angular_momentum_body(),
    )
    .unwrap();
    let mapped_units = PorousLatticeUnits3 {
        grid_origin_m: rotate(units().grid_origin_m).add(translation),
        ..units()
    };
    let rotated = SphericalPorousBody::new(
        mapped,
        moved,
        motion.properties(),
        mapped_units,
        motion.radius_m(),
        CELLS,
        &CancelGate::new(),
    )
    .unwrap();
    let transformed = rotated
        .step_candidate(
            radial,
            |[x, y, z]| a(rotate(v(forcing([y, z, x])))),
            rotate(external),
            PorousRotorControls {
                flow: PorousMomentumControls {
                    balance_origin_m: rotate(c.flow.balance_origin_m).add(translation),
                    ..c.flow
                },
                ..c
            },
            &CancelGate::new(),
        )
        .unwrap();
    close_v(
        world_spin(transformed.motion.body()),
        rotate(world_spin(base.motion.body())),
        2e-10,
    );
    close_v(
        transformed.angular_velocity_world_rad_s,
        rotate(base.angular_velocity_world_rad_s),
        2e-10,
    );
    close_v(
        transformed.guide_impulse_kg_m_s,
        rotate(base.guide_impulse_kg_m_s),
        2e-10,
    );
    close(
        transformed.relative_dissipation_j,
        base.relative_dissipation_j,
        3e-12,
    );
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let (r0, p0) = cell(base.motion.fluid(), x, y, z);
                let (r1, p1) = cell(transformed.motion.fluid(), z, x, y);
                close(r0, r1, 2e-12);
                close_v(p1, rotate(p0), 2e-11);
            }
        }
    }
}

#[test]
fn spherical_budget_cancel_and_late_rim_refusal_keep_both_accepted_owners() {
    let motion = spherical(true);
    let before = snapshot(motion.fluid());
    let body = motion.body();
    let c = controls();
    let limited = PorousRotorControls {
        flow: PorousMomentumControls {
            maximum_native_trials: 2,
            ..c.flow
        },
        ..c
    };
    assert!(matches!(
        motion.step_candidate(radial, forcing, Vec3::ZERO, limited, &CancelGate::new()),
        Err(PorousMomentumError::Budget { .. })
    ));
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        motion.step_candidate(radial, forcing, Vec3::ZERO, c, &gate),
        Err(PorousMomentumError::Cancelled)
    ));
    assert!(matches!(
        motion.step_candidate(
            |_| PorousRotorResistance::default(),
            |_| [0.0; 3],
            Vec3::new(0.4, 0.0, 0.0),
            c,
            &CancelGate::new()
        ),
        Err(PorousMomentumError::Invalid(
            "accepted spherical rim speed exceeds the declared domain"
        ))
    ));
    assert_eq!(snapshot(motion.fluid()), before);
    assert_eq!(motion.body(), body);
    assert_eq!(motion.accepted_steps(), 0);
}
