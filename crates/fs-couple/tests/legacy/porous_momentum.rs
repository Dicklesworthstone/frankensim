//! Actual periodic fluid populations, finite dry inertia and independent balances.
use fs_couple::porous_momentum::{
    PeriodicPorousSkeleton, PorousLatticeUnits3, PorousMomentumControls, PorousMomentumError,
};
use fs_exec::CancelGate;
use fs_lbm::d3q19::{BoundaryGrid3, BoundarySpec3, DarcyForchheimer3, E3, LinearResistance3, Q3};
use fs_lbm::equilibrium3;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const CELLS: usize = 64;
const SKELETON_MASS: f64 = 32.0;
const FLUID_V: Vec3 = Vec3::new(0.025, -0.01, 0.005);
const SKELETON_V: Vec3 = Vec3::new(-0.01, 0.0075, -0.0025);
fn a(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn v(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn close(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol,
        "{a:.17e} != {b:.17e}, tolerance {tol:e}"
    );
}
fn close_v(actual: Vec3, expected: Vec3, tol: f64) {
    for (x, y) in a(actual).into_iter().zip(a(expected)) {
        close(x, y, tol);
    }
}
fn units(dt: f64) -> PorousLatticeUnits3 {
    PorousLatticeUnits3 {
        cell_size_m: 0.25,
        timestep_s: dt,
        mass_per_density_cell_kg: 2.0,
        grid_origin_m: Vec3::new(-0.3, 0.2, 1.0),
    }
}
fn controls() -> PorousMomentumControls {
    PorousMomentumControls {
        velocity_scale_m_s: 0.05,
        maximum_speed_m_s: 0.075,
        difference_step: 1e-5,
        maximum_iterations: 24,
        maximum_native_trials: 256,
        maximum_trial_cells: CELLS * 256,
        linear_tolerance_kg_m_s: 1e-10,
        angular_tolerance_kg_m2_s: 1e-10,
        midpoint_tolerance_m_s: 1e-11,
        position_tolerance_m: 1e-12,
        work_tolerance_j: 1e-12,
        relative_mass_tolerance: 1e-12,
        balance_origin_m: Vec3::new(-0.7, 0.4, -0.3),
    }
}
fn body(velocity: Vec3) -> (RigidBodyState, MassProperties) {
    let properties =
        MassProperties::new(SKELETON_MASS, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    let pose = Pose::new(
        Vec3::new(0.1, 0.35, 1.7),
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, -0.3, 0.4), 0.4).unwrap(),
    )
    .unwrap();
    (
        RigidBodyState::new(pose, velocity.scale(SKELETON_MASS), Vec3::ZERO).unwrap(),
        properties,
    )
}
fn grid(dt: f64, nonuniform: bool) -> BoundaryGrid3 {
    let mut grid = BoundaryGrid3::new(4, 4, 4, 0.9, [0.0; 3], BoundarySpec3::periodic());
    for z in 0..4 {
        for y in 0..4 {
            for x in 0..4 {
                let (rho, u) = if nonuniform {
                    (
                        1.0 + 0.01 * x as f64 - 0.015 * z as f64,
                        Vec3::new(
                            0.02 + 0.003 * (y as f64 * std::f64::consts::FRAC_PI_2).sin(),
                            -0.01 + 0.002 * x as f64,
                            0.005 + 0.003 * (x as f64 * std::f64::consts::FRAC_PI_2).sin(),
                        ),
                    )
                } else {
                    (1.0, FLUID_V)
                };
                grid.set_populations(
                    x,
                    y,
                    z,
                    &equilibrium3(rho, a(u.scale(dt / units(dt).cell_size_m))),
                );
            }
        }
    }
    grid
}
fn fixture(dt: f64, nonuniform: bool, velocity: Vec3) -> PeriodicPorousSkeleton {
    let (body, properties) = body(velocity);
    PeriodicPorousSkeleton::new(
        grid(dt, nonuniform),
        body,
        properties,
        units(dt),
        CELLS,
        &CancelGate::new(),
    )
    .unwrap()
}
fn linear(coefficient: f64) -> DarcyForchheimer3 {
    LinearResistance3 {
        principal: [coefficient; 3],
        ..Default::default()
    }
    .into()
}
fn nonlinear() -> DarcyForchheimer3 {
    DarcyForchheimer3 {
        linear: LinearResistance3 {
            principal: [0.2, 0.8, 1.4],
            axes: [[0.8, 0.6, 0.0], [-0.6, 0.8, 0.0], [0.0, 0.0, 1.0]],
            body_force: [1e-4, -2e-4, 3e-4],
            ..Default::default()
        },
        quadratic: 4.0,
    }
}
fn snapshot(grid: &BoundaryGrid3) -> Vec<([u64; Q3], [u64; 3])> {
    (0..4)
        .flat_map(|z| {
            (0..4).flat_map(move |y| {
                (0..4).map(move |x| {
                    (
                        grid.populations(x, y, z).map(f64::to_bits),
                        grid.force_density(x, y, z).map(f64::to_bits),
                    )
                })
            })
        })
        .collect()
}

// Independent physical moments use actual populations and their held half-force,
// rather than the production coupled measurement or native velocity accessor.
fn cell(grid: &BoundaryGrid3, x: usize, y: usize, z: usize) -> (f64, Vec3) {
    let f = grid.populations(x, y, z);
    let force = grid.force_density(x, y, z);
    let rho = f.iter().sum::<f64>();
    let momentum: [f64; 3] = std::array::from_fn(|i| {
        f.iter()
            .zip(E3)
            .map(|(f, e)| f * [e.0 as f64, e.1 as f64, e.2 as f64][i])
            .sum::<f64>()
            + 0.5 * force[i]
    });
    (rho, v(momentum))
}
fn physical(grid: &BoundaryGrid3, units: PorousLatticeUnits3, origin: Vec3) -> (Vec3, Vec3, f64) {
    let velocity_scale = units.cell_size_m / units.timestep_s;
    let impulse_scale = units.mass_per_density_cell_kg * velocity_scale;
    let mut p = Vec3::ZERO;
    let mut l = Vec3::ZERO;
    let mut kinetic = 0.0;
    for z in 0..4 {
        for y in 0..4 {
            for x in 0..4 {
                let (rho, m) = cell(grid, x, y, z);
                let momentum = m.scale(impulse_scale);
                let position = units.grid_origin_m.add(
                    Vec3::new(x as f64 + 0.5, y as f64 + 0.5, z as f64 + 0.5)
                        .scale(units.cell_size_m),
                );
                p = p.add(momentum);
                l = l.add(position.sub(origin).cross(momentum));
                kinetic += 0.5
                    * units.mass_per_density_cell_kg
                    * velocity_scale
                    * velocity_scale
                    * m.dot(m)
                    / rho;
            }
        }
    }
    (p, l, kinetic)
}

#[test]
fn uniform_native_recoil_matches_two_mass_midpoint_and_actual_kinetic_loss() {
    let motion = fixture(0.5, false, SKELETON_V);
    let before = snapshot(motion.fluid());
    let original = motion.body();
    let c = controls();
    let u = motion.units();
    let f0 = physical(motion.fluid(), u, c.balance_origin_m);
    let result = motion
        .step_candidate(linear(0.8), Vec3::ZERO, c, &CancelGate::new())
        .unwrap();
    let fluid_mass = CELLS as f64 * u.mass_per_density_cell_kg;
    let gamma = CELLS as f64 * 0.8 * u.mass_per_density_cell_kg;
    let reduced = fluid_mass * SKELETON_MASS / (fluid_mass + SKELETON_MASS);
    let expected = FLUID_V
        .sub(SKELETON_V)
        .scale(gamma / (1.0 + 0.5 * gamma / reduced));
    close_v(result.momentum.linear_increment_world, expected, 2e-11);
    let next_v = result
        .motion
        .body()
        .center_of_mass_velocity_world(motion.properties())
        .unwrap();
    close_v(
        next_v,
        SKELETON_V.add(expected.scale(1.0 / SKELETON_MASS)),
        2e-12,
    );
    for z in 0..4 {
        for y in 0..4 {
            for x in 0..4 {
                let actual =
                    v(result.motion.fluid().velocity(x, y, z)).scale(u.cell_size_m / u.timestep_s);
                close_v(actual, FLUID_V.sub(expected.scale(1.0 / fluid_mass)), 2e-12);
            }
        }
    }
    let relative_after = FLUID_V.sub(expected.scale(1.0 / fluid_mass)).sub(next_v);
    assert!(relative_after.dot(FLUID_V.sub(SKELETON_V)) < 0.0); // admitted stiff midpoint reversal
    let f1 = physical(result.motion.fluid(), u, c.balance_origin_m);
    close_v(
        f1.0.add(result.motion.body().linear_momentum_world()),
        f0.0.add(original.linear_momentum_world()),
        3e-11,
    );
    let kinetic_change = f1.2 - f0.2
        + result
            .motion
            .body()
            .kinetic_energy(motion.properties())
            .unwrap()
        - original.kinetic_energy(motion.properties()).unwrap();
    close(kinetic_change + result.relative_dissipation_j, 0.0, 2e-12);
    assert!(result.relative_dissipation_j > 1e-3);
    close(result.other_fluid_kinetic_change_j, 0.0, 2e-12);
    let midpoint_position = original
        .pose()
        .position_world()
        .scale(0.5)
        .add(result.motion.body().pose().position_world().scale(0.5));
    let centroid = u.grid_origin_m.add(Vec3::new(0.5, 0.5, 0.5));
    close_v(
        result.guide_angular_impulse_kg_m2_s,
        centroid.sub(midpoint_position).cross(expected).scale(-1.0),
        2e-11,
    );
    assert!(norm(result.guide_angular_impulse_kg_m2_s) > 1e-2);
    assert_eq!(result.motion.body().angular_momentum_body(), Vec3::ZERO);
    assert_eq!(
        result.motion.body().pose().orientation(),
        original.pose().orientation()
    );
    assert_eq!(result.motion.accepted_steps(), 1);
    close(result.motion.time_s(), 0.5, 0.0);
    assert_eq!(snapshot(motion.fluid()), before);
    assert_eq!(motion.body(), original);
}

#[test]
fn nonlinear_nonuniform_flow_replays_native_populations_and_closes_guided_angular_transport() {
    let motion = fixture(0.5, true, SKELETON_V);
    let c = controls();
    let u = motion.units();
    let law = nonlinear();
    let force = Vec3::new(0.03, -0.02, 0.025);
    let result = motion
        .step_candidate(law, force, c, &CancelGate::new())
        .unwrap();
    let mut replay = motion.fluid().clone();
    let mut selected = law;
    selected.linear.skeleton_velocity = a(result
        .midpoint_velocity_m_s
        .scale(u.timestep_s / u.cell_size_m));
    let native = replay
        .step_with_darcy_forchheimer_and_wall_momentum(
            |_| selected,
            &[],
            a(result
                .exchange_origin_m
                .sub(u.grid_origin_m)
                .scale(u.cell_size_m.recip())),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(snapshot(&replay), snapshot(result.motion.fluid()));
    assert_eq!(native, result.exchange);
    assert!(result.native_trials > 3);
    let f0 = physical(motion.fluid(), u, c.balance_origin_m);
    let f1 = physical(result.motion.fluid(), u, c.balance_origin_m);
    let old = motion.body();
    let next = result.motion.body();
    let midpoint_position = old
        .pose()
        .position_world()
        .scale(0.5)
        .add(next.pose().position_world().scale(0.5));
    let jscale = u.mass_per_density_cell_kg * u.cell_size_m / u.timestep_s;
    let external_fluid = v(law.linear.body_force).scale(CELLS as f64 * jscale);
    let external_body = force.scale(u.timestep_s);
    close_v(
        f1.0.sub(f0.0).add(
            next.linear_momentum_world()
                .sub(old.linear_momentum_world()),
        ),
        external_fluid.add(external_body),
        2e-11,
    );
    // Reconstruct the distributed reaction and its physical COM moment from
    // endpoint relative velocities, independently of the aggregate receipt.
    let resistance = |velocity: Vec3| {
        let relative = velocity.sub(v(selected.linear.skeleton_velocity));
        let mut result = relative.scale(law.quadratic * norm(relative));
        for k in 0..3 {
            let axis = v(law.linear.axes[k]);
            result = result.add(axis.scale(law.linear.principal[k] * axis.dot(relative)));
        }
        result
    };
    let mut expected_guide = Vec3::ZERO;
    let mut drive_moment = Vec3::ZERO;
    for z in 0..4 {
        for y in 0..4 {
            for x in 0..4 {
                let (r0, p0) = cell(motion.fluid(), x, y, z);
                let (r1, p1) = cell(result.motion.fluid(), x, y, z);
                let j = resistance(p0.scale(r0.recip()))
                    .add(resistance(p1.scale(r1.recip())))
                    .scale(0.5 * jscale);
                let position = u.grid_origin_m.add(
                    Vec3::new(x as f64 + 0.5, y as f64 + 0.5, z as f64 + 0.5).scale(u.cell_size_m),
                );
                expected_guide = expected_guide.sub(position.sub(midpoint_position).cross(j));
                drive_moment = drive_moment.add(
                    position
                        .sub(c.balance_origin_m)
                        .cross(v(law.linear.body_force).scale(jscale)),
                );
            }
        }
    }
    close_v(result.guide_angular_impulse_kg_m2_s, expected_guide, 3e-11);
    let body_l0 = old
        .pose()
        .position_world()
        .sub(c.balance_origin_m)
        .cross(old.linear_momentum_world());
    let body_l1 = next
        .pose()
        .position_world()
        .sub(c.balance_origin_m)
        .cross(next.linear_momentum_world());
    let supplied = drive_moment
        .add(
            midpoint_position
                .sub(c.balance_origin_m)
                .cross(external_body),
        )
        .add(expected_guide)
        .add(result.periodic_angular_transport_kg_m2_s);
    close_v(f1.1.sub(f0.1).add(body_l1.sub(body_l0)), supplied, 3e-11);
    assert!(norm(result.periodic_angular_transport_kg_m2_s) > 1e-5);
    assert!(norm(expected_guide) > 1e-3);
    close(result.mechanical_work_residual_j, 0.0, 1e-12);
    close(
        result.other_fluid_kinetic_change_j,
        f1.2 - f0.2 - result.fluid_external_work_j - result.fluid_drag_work_j,
        2e-12,
    );
    let shifted = motion
        .step_candidate(
            law,
            force,
            PorousMomentumControls {
                balance_origin_m: Vec3::new(1.4, -0.8, 2.1),
                ..c
            },
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(
        snapshot(shifted.motion.fluid()),
        snapshot(result.motion.fluid())
    );
    assert_eq!(shifted.motion.body(), result.motion.body());
}

#[test]
fn quadratic_two_mass_relaxation_converges_to_the_continuous_recoil_trajectory() {
    let law = DarcyForchheimer3 {
        quadratic: 4.0,
        ..Default::default()
    };
    let fluid_mass = CELLS as f64 * units(0.5).mass_per_density_cell_kg;
    let total = fluid_mass + SKELETON_MASS;
    let center = FLUID_V
        .scale(fluid_mass / total)
        .add(SKELETON_V.scale(SKELETON_MASS / total));
    let slip = FLUID_V.sub(SKELETON_V);
    let resistance =
        CELLS as f64 * law.quadratic * units(0.5).mass_per_density_cell_kg / units(0.5).cell_size_m;
    let alpha = resistance * (1.0 / fluid_mass + 1.0 / SKELETON_MASS);
    let duration = 0.5;
    let exact =
        center.sub(slip.scale((fluid_mass / total) / (1.0 + alpha * norm(slip) * duration)));
    let run = |steps: usize| {
        let mut motion = fixture(duration / steps as f64, false, SKELETON_V);
        let initial = physical(motion.fluid(), motion.units(), controls().balance_origin_m)
            .0
            .add(motion.body().linear_momentum_world());
        for _ in 0..steps {
            motion = motion
                .step_candidate(law, Vec3::ZERO, controls(), &CancelGate::new())
                .unwrap()
                .motion;
        }
        let final_momentum = physical(motion.fluid(), motion.units(), controls().balance_origin_m)
            .0
            .add(motion.body().linear_momentum_world());
        close_v(initial, final_momentum, 2e-10);
        assert_eq!(motion.accepted_steps(), steps as u64);
        close(motion.time_s(), duration, 0.0);
        norm(
            motion
                .body()
                .center_of_mass_velocity_world(motion.properties())
                .unwrap()
                .sub(exact),
        )
    };
    let errors = [run(4), run(8), run(16)];
    for pair in errors.windows(2) {
        assert!(
            (3.0..5.0).contains(&(pair[0] / pair[1])),
            "quadratic refinement errors {errors:?}"
        );
    }
    assert!(
        errors[2] > 1e-10 && errors[2] < 2e-5,
        "fine velocity error {}",
        errors[2]
    );
}

#[test]
fn bounded_trials_cancellation_and_late_endpoint_refusal_preserve_both_original_states() {
    let motion = fixture(0.5, true, SKELETON_V);
    let before = snapshot(motion.fluid());
    let body = motion.body();
    let failure = motion
        .step_candidate(
            nonlinear(),
            Vec3::ZERO,
            PorousMomentumControls {
                maximum_native_trials: 2,
                ..controls()
            },
            &CancelGate::new(),
        )
        .err()
        .unwrap();
    assert!(matches!(failure, PorousMomentumError::Budget { .. }));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        motion.step_candidate(nonlinear(), Vec3::ZERO, controls(), &cancelled),
        Err(PorousMomentumError::Cancelled)
    ));
    assert_eq!(snapshot(motion.fluid()), before);
    assert_eq!(motion.body(), body);
    let motion = fixture(0.5, false, Vec3::new(0.07, 0.0, 0.0));
    let before = snapshot(motion.fluid());
    let body = motion.body();
    // Midpoint remains below 0.075 m/s; the converged actual endpoint exceeds it.
    let failure = motion
        .step_candidate(
            linear(0.0),
            Vec3::new(0.5, 0.0, 0.0),
            controls(),
            &CancelGate::new(),
        )
        .err()
        .unwrap();
    assert!(matches!(
        failure,
        PorousMomentumError::Invalid("accepted skeleton speed exceeds the declared domain")
    ));
    assert_eq!(snapshot(motion.fluid()), before);
    assert_eq!(motion.body(), body);
    assert_eq!(motion.accepted_steps(), 0);
}
