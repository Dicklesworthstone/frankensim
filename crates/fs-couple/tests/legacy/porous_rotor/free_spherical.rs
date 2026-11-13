//! Free six-component recoil with independent population and rigid inventories.

use super::*;
use fs_couple::porous_momentum::free_spherical::{
    FreeSphericalPorousBody, FreeSphericalPorousStep,
};
use fs_mbd::{Gravity, RigidBodyIntegrator, Wrench};

fn full_body(position: Vec3, velocity: Vec3, omega: Vec3) -> (RigidBodyState, MassProperties) {
    let principal = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -1.0, 2.0), 0.47).unwrap();
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(2.0, 1.0, -1.0), 0.31).unwrap();
    let properties =
        MassProperties::from_principal_frame(16.0, Vec3::ZERO, Vec3::new(1.2, 1.6, 2.0), principal)
            .unwrap();
    let body = RigidBodyState::new(
        Pose::new(position, orientation).unwrap(),
        velocity.scale(properties.mass()),
        properties.angular_momentum_body_checked(omega).unwrap(),
    )
    .unwrap();
    (body, properties)
}

fn fixture(velocity: Vec3, omega: Vec3, nonuniform: bool) -> FreeSphericalPorousBody {
    let (body, properties) = full_body(COM, velocity, omega);
    FreeSphericalPorousBody::new(
        grid(nonuniform),
        body,
        properties,
        units(),
        cylinder().radius_m,
        CELLS,
        &CancelGate::new(),
    )
    .unwrap()
}

fn controls() -> PorousRotorControls {
    let c = super::controls();
    PorousRotorControls {
        flow: PorousMomentumControls {
            maximum_native_trials: 256,
            maximum_trial_cells: CELLS * 256,
            ..c.flow
        },
        ..c
    }
}

fn law() -> PorousRotorResistance {
    PorousRotorResistance {
        darcy: 0.12,
        forchheimer: 4.0,
    }
}

fn world_spin(body: RigidBodyState) -> Vec3 {
    body.pose()
        .orientation()
        .rotate_body_to_world(body.angular_momentum_body())
}

fn body_angular(body: RigidBodyState, origin: Vec3) -> Vec3 {
    body.pose()
        .position_world()
        .sub(origin)
        .cross(body.linear_momentum_world())
        .add(world_spin(body))
}

fn replay(
    before: &FreeSphericalPorousBody,
    step: &FreeSphericalPorousStep,
    resistance: PorousRotorResistance,
    centre: Vec3,
    force: impl Fn([usize; 3]) -> [f64; 3],
) -> (BoundaryGrid3, fs_lbm::d3q19::ResistanceMomentumExchange3) {
    let mut fluid = before.fluid().clone();
    let links = fluid.boundary_links();
    let speed = units().cell_size_m / units().timestep_s;
    let receipt = fluid
        .step_with_darcy_forchheimer_and_wall_momentum(
            |cell| {
                let arm = position(cell).sub(centre);
                let inside = norm(arm) <= before.radius_m();
                let law = if inside {
                    resistance
                } else {
                    PorousRotorResistance::default()
                };
                DarcyForchheimer3 {
                    linear: LinearResistance3 {
                        principal: [law.darcy; 3],
                        skeleton_velocity: if inside {
                            a(step
                                .midpoint_velocity_m_s
                                .add(step.angular_velocity_world_rad_s.cross(arm)))
                            .map(|u| u / speed)
                        } else {
                            [0.0; 3]
                        },
                        body_force: force(cell),
                        ..Default::default()
                    },
                    quadratic: law.forchheimer,
                }
            },
            &links,
            a(centre
                .sub(units().grid_origin_m)
                .scale(1.0 / units().cell_size_m)),
            &CancelGate::new(),
        )
        .unwrap();
    (fluid, receipt)
}

#[test]
fn zero_drag_keeps_free_translation_and_original_asymmetric_gyroscopic_motion() {
    let velocity = Vec3::new(0.012, -0.009, 0.007);
    let mut motion = fixture(velocity, Vec3::new(0.015, -0.02, 0.06), false);
    let initial = motion.body();
    let spin = world_spin(initial);
    let energy = initial.kinetic_energy(motion.properties()).unwrap();
    let origin = controls().flow.balance_origin_m;
    let angular = body_angular(initial, origin);
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
                PorousRotorResistance::default(),
                |_| [0.0; 3],
                Vec3::ZERO,
                Vec3::ZERO,
                controls(),
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(step.mechanics, expected);
        close_v(step.midpoint_velocity_m_s, velocity, 2e-17);
        close(step.relative_dissipation_j, 0.0, 0.0);
        close_v(step.body_linear_change_kg_m_s, Vec3::ZERO, 0.0);
        motion = step.motion;
        close_v(world_spin(motion.body()), spin, 5e-14);
        close_v(body_angular(motion.body(), origin), angular, 5e-14);
        close(
            motion.body().kinetic_energy(motion.properties()).unwrap(),
            energy,
            5e-15,
        );
    }
    close_v(
        motion.body().pose().position_world(),
        COM.add(velocity.scale(2.0)),
        5e-16,
    );
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
fn nonlinear_force_and_torque_recoil_close_about_a_nonzero_origin_without_guides() {
    let motion = fixture(
        Vec3::new(-0.006, 0.008, -0.004),
        Vec3::new(0.015, -0.02, 0.06),
        true,
    );
    let c = controls();
    let external_force = Vec3::new(0.009, -0.004, 0.006);
    let external_couple = Vec3::new(0.004, -0.003, 0.002);
    let mut sources = 0;
    let step = motion
        .step_candidate(
            law(),
            |cell| {
                sources += 1;
                forcing(cell)
            },
            external_force,
            external_couple,
            c,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(sources, CELLS);
    assert!(step.native_trials > 2);
    let u = units();
    let speed = u.cell_size_m / u.timestep_s;
    let jscale = u.mass_per_density_cell_kg * speed;
    let lscale = jscale * u.cell_size_m;
    let centre = motion
        .body()
        .pose()
        .position_world()
        .add(step.motion.body().pose().position_world())
        .scale(0.5);
    close_v(step.exchange_origin_m, centre, 0.0);
    assert!(norm(centre.sub(COM)) > 1e-4);
    let mut drag = Vec3::ZERO;
    let mut drag_moment = Vec3::ZERO;
    let mut loss = 0.0;
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let arm = position([x, y, z]).sub(centre);
                if norm(arm) > motion.radius_m() {
                    continue;
                }
                let imposed = step
                    .midpoint_velocity_m_s
                    .add(step.angular_velocity_world_rad_s.cross(arm))
                    .scale(1.0 / speed);
                let (rho0, p0) = cell(motion.fluid(), x, y, z);
                let (rho1, p1) = cell(step.motion.fluid(), x, y, z);
                let s0 = p0.scale(1.0 / rho0).sub(imposed);
                let s1 = p1.scale(1.0 / rho1).sub(imposed);
                let j = s0
                    .scale(law().darcy + law().forchheimer * norm(s0))
                    .add(s1.scale(law().darcy + law().forchheimer * norm(s1)))
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
    let linear_change = step
        .motion
        .body()
        .linear_momentum_world()
        .sub(motion.body().linear_momentum_world());
    let spin_change = world_spin(step.motion.body()).sub(world_spin(motion.body()));
    close_v(
        linear_change,
        drag.add(external_force.scale(u.timestep_s)),
        3e-10,
    );
    close_v(
        spin_change,
        drag_moment.add(external_couple.scale(u.timestep_s)),
        3e-10,
    );
    assert!(norm(linear_change) > 0.01);
    assert!(norm(spin_change) > 1e-3);
    assert!(norm(drag_moment.cross(world_spin(motion.body()))) > 1e-5);
    let body_change = body_angular(step.motion.body(), c.flow.balance_origin_m)
        .sub(body_angular(motion.body(), c.flow.balance_origin_m));
    close_v(step.body_angular_change_kg_m2_s, body_change, 3e-14);
    let f0 = physical(motion.fluid(), c.flow.balance_origin_m);
    let f1 = physical(step.motion.fluid(), c.flow.balance_origin_m);
    let force = v(step.exchange.external_impulse).scale(jscale);
    let wall = v(step.exchange.walls.fluid_population_impulse).scale(jscale);
    let shift = |moment: Vec3, impulse: Vec3| {
        moment.add(centre.sub(c.flow.balance_origin_m).cross(impulse))
    };
    close_v(
        f1.0.sub(f0.0).add(linear_change),
        force.add(wall).add(external_force.scale(u.timestep_s)),
        3e-10,
    );
    let supplied = shift(
        v(step.exchange.external_angular_impulse).scale(lscale),
        force,
    )
    .add(shift(
        v(step.exchange.walls.fluid_population_angular_impulse).scale(lscale),
        wall,
    ))
    .add(shift(
        external_couple.scale(u.timestep_s),
        external_force.scale(u.timestep_s),
    ));
    close_v(f1.1.sub(f0.1).add(body_change), supplied, 3e-10);
    let kinetic = step
        .motion
        .body()
        .kinetic_energy(motion.properties())
        .unwrap()
        - motion.body().kinetic_energy(motion.properties()).unwrap();
    let external_work = (step.midpoint_velocity_m_s.dot(external_force)
        + step.angular_velocity_world_rad_s.dot(external_couple))
        * u.timestep_s;
    close(step.body_kinetic_change_j, kinetic, 3e-14);
    close(
        kinetic + step.fluid_drag_work_j + loss - external_work,
        0.0,
        3e-12,
    );
    close(
        step.other_fluid_kinetic_change_j,
        f1.2 - f0.2 - step.fluid_external_work_j - step.fluid_drag_work_j,
        3e-12,
    );
    let (native, receipt) = replay(&motion, &step, law(), centre, forcing);
    assert_eq!(snapshot(&native), snapshot(step.motion.fluid()));
    assert_eq!(receipt, step.exchange);
    let shifted = motion
        .step_candidate(
            law(),
            forcing,
            external_force,
            external_couple,
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
fn moving_support_uses_the_solved_midpoint_and_retargets_existing_fluid_forces() {
    let mut motion = fixture(Vec3::new(0.05, 0.0, 0.0), Vec3::ZERO, false);
    let resistance = PorousRotorResistance {
        darcy: 0.02,
        forchheimer: 0.0,
    };
    let mut changed = false;
    let mut released_forced_cell = false;
    for _ in 0..2 {
        let old_centre = motion.body().pose().position_world();
        let step = motion
            .step_candidate(
                resistance,
                |_| [0.0; 3],
                Vec3::ZERO,
                Vec3::ZERO,
                controls(),
                &CancelGate::new(),
            )
            .unwrap();
        let mut moved_mask = false;
        for z in 0..N {
            for y in 0..N {
                for x in 0..N {
                    let old = norm(position([x, y, z]).sub(old_centre)) <= motion.radius_m();
                    let new =
                        norm(position([x, y, z]).sub(step.exchange_origin_m)) <= motion.radius_m();
                    moved_mask |= old != new;
                    if !new && norm(v(motion.fluid().force_density(x, y, z))) > 0.0 {
                        released_forced_cell = true;
                        assert_eq!(step.motion.fluid().force_density(x, y, z), [0.0; 3]);
                    }
                }
            }
        }
        let (native, receipt) = replay(&motion, &step, resistance, step.exchange_origin_m, |_| {
            [0.0; 3]
        });
        assert_eq!(snapshot(&native), snapshot(step.motion.fluid()));
        assert_eq!(receipt, step.exchange);
        if moved_mask {
            let (frozen, _) = replay(&motion, &step, resistance, old_centre, |_| [0.0; 3]);
            assert_ne!(snapshot(&frozen), snapshot(step.motion.fluid()));
        }
        changed |= moved_mask;
        let origin = controls().flow.balance_origin_m;
        let f0 = physical(motion.fluid(), origin);
        let f1 = physical(step.motion.fluid(), origin);
        let jscale = units().mass_per_density_cell_kg * units().cell_size_m / units().timestep_s;
        let wall = v(step.exchange.walls.fluid_population_impulse).scale(jscale);
        let wall_moment = v(step.exchange.walls.fluid_population_angular_impulse)
            .scale(jscale * units().cell_size_m)
            .add(step.exchange_origin_m.sub(origin).cross(wall));
        close_v(
            f1.0.sub(f0.0)
                .add(step.motion.body().linear_momentum_world())
                .sub(motion.body().linear_momentum_world()),
            wall,
            3e-10,
        );
        close_v(
            f1.1.sub(f0.1)
                .add(body_angular(step.motion.body(), origin))
                .sub(body_angular(motion.body(), origin)),
            wall_moment,
            3e-10,
        );
        motion = step.motion;
    }
    assert!(
        changed,
        "fixture must move lattice cells across the spherical support"
    );
    assert!(
        released_forced_cell,
        "a previously forced cell must leave the moving support"
    );
    assert_eq!(motion.accepted_steps(), 2);
}

#[test]
fn endpoint_containment_combined_speed_budget_and_cancellation_refuse_atomically() {
    let motion = fixture(
        Vec3::new(-0.006, 0.008, -0.004),
        Vec3::new(0.015, -0.02, 0.06),
        true,
    );
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
        motion.step_candidate(
            law(),
            forcing,
            Vec3::ZERO,
            Vec3::ZERO,
            limited,
            &CancelGate::new()
        ),
        Err(PorousMomentumError::Budget { .. })
    ));
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        motion.step_candidate(law(), forcing, Vec3::ZERO, Vec3::ZERO, c, &gate),
        Err(PorousMomentumError::Cancelled)
    ));
    let during_sampling = CancelGate::new();
    let mut sampled = 0;
    assert!(matches!(
        motion.step_candidate(
            law(),
            |cell| {
                sampled += 1;
                if sampled == 3 {
                    during_sampling.request();
                }
                forcing(cell)
            },
            Vec3::ZERO,
            Vec3::ZERO,
            c,
            &during_sampling,
        ),
        Err(PorousMomentumError::Cancelled)
    ));
    assert!(sampled >= 3 && sampled < CELLS);
    assert_eq!(snapshot(motion.fluid()), before);
    assert_eq!(motion.body(), body);
    assert_eq!(motion.accepted_steps(), 0);
    close(motion.time_s(), 0.0, 0.0);

    let right = units().grid_origin_m.x + N as f64 * units().cell_size_m;
    let centre = Vec3::new(right - cylinder().radius_m - 0.02, COM.y, COM.z);
    let (body, properties) = full_body(centre, Vec3::new(0.06, 0.0, 0.0), Vec3::ZERO);
    let near_wall = FreeSphericalPorousBody::new(
        grid(false),
        body,
        properties,
        units(),
        cylinder().radius_m,
        CELLS,
        &CancelGate::new(),
    )
    .unwrap();
    // The free midpoint sphere fits; its endpoint protrudes beyond the wall.
    assert!(centre.x + 0.5 * units().timestep_s * 0.06 + cylinder().radius_m < right);
    assert!(matches!(
        near_wall.step_candidate(
            PorousRotorResistance::default(),
            |_| [0.0; 3],
            Vec3::ZERO,
            Vec3::ZERO,
            c,
            &CancelGate::new(),
        ),
        Err(PorousMomentumError::Invalid(
            "complete moving porous sphere must stay inside the box"
        ))
    ));
    assert_eq!(near_wall.body(), body);
    let fast = fixture(Vec3::new(0.05, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.06), false);
    assert!(matches!(
        fast.step_candidate(
            PorousRotorResistance::default(),
            |_| [0.0; 3],
            Vec3::ZERO,
            Vec3::ZERO,
            c,
            &CancelGate::new(),
        ),
        Err(PorousMomentumError::Invalid(
            "free porous translation plus rim speed exceeds the declared domain"
        ))
    ));
}
