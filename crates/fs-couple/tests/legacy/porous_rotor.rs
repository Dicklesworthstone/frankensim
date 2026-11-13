//! Native populations, finite polar inertia, and independently reconstructed wrenches.
#[path = "porous_rotor/free_body.rs"]
mod free_body;
#[path = "porous_rotor/free_spherical.rs"]
mod free_spherical;
#[path = "porous_rotor/spherical.rs"]
mod spherical;

use fs_couple::porous_momentum::rotor::{
    AxisymmetricPorousRotor, PorousRotorControls, PorousRotorCylinder, PorousRotorResistance,
};
use fs_couple::porous_momentum::{
    PorousLatticeUnits3, PorousMomentumControls, PorousMomentumError,
};
use fs_exec::CancelGate;
use fs_lbm::d3q19::{
    BoundaryGrid3, BoundarySpec3, DarcyForchheimer3, E3, FaceBoundary3, LinearResistance3, Q3,
};
use fs_lbm::equilibrium3;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const N: usize = 8;
const CELLS: usize = N * N * N;
const INERTIA: f64 = 2.0;
const AXIS: Vec3 = Vec3::new(0.0, 0.0, 1.0);
const COM: Vec3 = Vec3::new(0.5, 1.25, 2.0);
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
fn close_v(x: Vec3, y: Vec3, tol: f64) {
    for (a, b) in a(x).into_iter().zip(a(y)) {
        close(a, b, tol);
    }
}
fn units() -> PorousLatticeUnits3 {
    PorousLatticeUnits3 {
        cell_size_m: 0.25,
        timestep_s: 0.5,
        mass_per_density_cell_kg: 2.0,
        grid_origin_m: Vec3::new(-0.5, 0.25, 1.0),
    }
}
fn controls() -> PorousRotorControls {
    PorousRotorControls {
        flow: PorousMomentumControls {
            velocity_scale_m_s: 0.04,
            maximum_speed_m_s: 0.075,
            difference_step: 1e-5,
            maximum_iterations: 24,
            maximum_native_trials: 128,
            maximum_trial_cells: CELLS * 128,
            linear_tolerance_kg_m_s: 3e-10,
            angular_tolerance_kg_m2_s: 3e-10,
            midpoint_tolerance_m_s: 1e-11,
            position_tolerance_m: 1e-12,
            work_tolerance_j: 2e-12,
            relative_mass_tolerance: 1e-12,
            balance_origin_m: Vec3::new(-0.75, 0.5, -0.25),
        },
        maximum_rotation_rad: 0.15,
    }
}
fn cylinder() -> PorousRotorCylinder {
    PorousRotorCylinder {
        axis_world: AXIS,
        radius_m: 0.55,
        half_height_m: 0.5,
    }
}
fn body(omega: f64, orientation: UnitQuaternion) -> (RigidBodyState, MassProperties) {
    let properties = MassProperties::new(16.0, Vec3::ZERO, Vec3::new(1.2, 1.6, INERTIA)).unwrap();
    (
        RigidBodyState::new(
            Pose::new(COM, orientation).unwrap(),
            Vec3::ZERO,
            AXIS.scale(INERTIA * omega),
        )
        .unwrap(),
        properties,
    )
}
fn grid(nonuniform: bool) -> BoundaryGrid3 {
    let mut result = BoundaryGrid3::new(
        N,
        N,
        N,
        0.9,
        [0.0; 3],
        BoundarySpec3::new([FaceBoundary3::stationary_wall(); 6]),
    );
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let (rho, velocity) = if nonuniform {
                    (
                        1.0 + 0.002 * x as f64 - 0.003 * z as f64,
                        Vec3::new(
                            0.018 + 0.001 * y as f64,
                            -0.014 + 0.002 * z as f64,
                            0.008 + 0.001 * x as f64,
                        ),
                    )
                } else {
                    (1.0, Vec3::ZERO)
                };
                result.set_populations(
                    x,
                    y,
                    z,
                    &equilibrium3(
                        rho,
                        a(velocity.scale(units().timestep_s / units().cell_size_m)),
                    ),
                );
            }
        }
    }
    result
}
fn fixture(omega: f64, nonuniform: bool) -> AxisymmetricPorousRotor {
    let (body, properties) = body(omega, UnitQuaternion::IDENTITY);
    AxisymmetricPorousRotor::new(
        grid(nonuniform),
        body,
        properties,
        units(),
        cylinder(),
        CELLS,
        &CancelGate::new(),
    )
    .unwrap()
}
fn coefficient(r: f64, z: f64) -> PorousRotorResistance {
    PorousRotorResistance {
        darcy: 0.08 * (1.0 + 0.8 * r + 0.6 * z),
        forchheimer: 3.0 + 0.3 * r,
    }
}
fn forcing([x, y, z]: [usize; 3]) -> [f64; 3] {
    [
        1e-5 * (1 + x) as f64,
        -2e-5 * (1 + z) as f64,
        1e-5 * y as f64,
    ]
}
fn position([x, y, z]: [usize; 3]) -> Vec3 {
    units()
        .grid_origin_m
        .add(Vec3::new(x as f64 + 0.5, y as f64 + 0.5, z as f64 + 0.5).scale(units().cell_size_m))
}
fn sample(cell: [usize; 3]) -> (PorousRotorResistance, Vec3) {
    let arm = position(cell).sub(COM);
    let tangent = AXIS.cross(arm);
    let r = norm(tangent);
    if r <= cylinder().radius_m && arm.z.abs() <= cylinder().half_height_m {
        (coefficient(r, arm.z), tangent)
    } else {
        (PorousRotorResistance::default(), Vec3::ZERO)
    }
}
fn snapshot(grid: &BoundaryGrid3) -> Vec<([u64; Q3], [u64; 3])> {
    (0..N)
        .flat_map(|z| {
            (0..N).flat_map(move |y| {
                (0..N).map(move |x| {
                    (
                        grid.populations(x, y, z).map(f64::to_bits),
                        grid.force_density(x, y, z).map(f64::to_bits),
                    )
                })
            })
        })
        .collect()
}
// Read the conserved population moment independently of the native accessor.
fn cell(grid: &BoundaryGrid3, x: usize, y: usize, z: usize) -> (f64, Vec3) {
    let f = grid.populations(x, y, z);
    let force = grid.force_density(x, y, z);
    let momentum = std::array::from_fn(|i| {
        f.iter()
            .zip(E3)
            .map(|(f, e)| f * [e.0 as f64, e.1 as f64, e.2 as f64][i])
            .sum::<f64>()
            + 0.5 * force[i]
    });
    (f.iter().sum(), v(momentum))
}
fn physical(grid: &BoundaryGrid3, origin: Vec3) -> (Vec3, Vec3, f64) {
    let u = units();
    let speed = u.cell_size_m / u.timestep_s;
    let mut result = (Vec3::ZERO, Vec3::ZERO, 0.0);
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let (rho, p) = cell(grid, x, y, z);
                let momentum = p.scale(u.mass_per_density_cell_kg * speed);
                result.0 = result.0.add(momentum);
                result.1 = result
                    .1
                    .add(position([x, y, z]).sub(origin).cross(momentum));
                result.2 += 0.5 * u.mass_per_density_cell_kg * speed * speed * p.dot(p) / rho;
            }
        }
    }
    result
}

#[test]
fn finite_rotor_slows_and_spins_the_native_fluid_with_actual_work_closure() {
    let motion = fixture(0.04, false);
    let before = snapshot(motion.fluid());
    let result = motion
        .step_candidate(
            coefficient,
            |_| [0.0; 3],
            0.0,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
    let omega = result.motion.angular_speed_rad_s().unwrap();
    assert!(omega > 0.0 && omega < 0.035, "finite recoil speed {omega}");
    assert!(result.motion.body().pose().orientation() != motion.body().pose().orientation());
    assert!(result.fluid_angular_change_kg_m2_s.z > 1e-3);
    assert!(result.relative_dissipation_j > 1e-5);
    close(
        result.rotor_kinetic_change_j,
        0.5 * INERTIA * (omega * omega - 0.04 * 0.04),
        1e-14,
    );
    close(
        result.rotor_kinetic_change_j + result.fluid_drag_work_j + result.relative_dissipation_j,
        0.0,
        2e-12,
    );
    close(result.guide_work_j, 0.0, 1e-14);
    assert_eq!(snapshot(motion.fluid()), before);
    assert_eq!(motion.accepted_steps(), 0);
    assert_eq!(result.motion.accepted_steps(), 1);
    close(result.motion.time_s(), units().timestep_s, 0.0);
}

#[test]
fn nonuniform_nonlinear_flow_uses_one_sample_and_closes_drag_bearing_and_wall_wrenches() {
    let motion = fixture(0.03, true);
    let c = controls();
    let mut coefficient_calls = 0;
    let mut force_calls = 0;
    let drive = 0.025;
    let result = motion
        .step_candidate(
            |r, z| {
                coefficient_calls += 1;
                coefficient(r, z)
            },
            |cell| {
                force_calls += 1;
                forcing(cell)
            },
            drive,
            c,
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(coefficient_calls, 64);
    assert_eq!(force_calls, CELLS);
    assert!(result.native_trials > 2);
    let u = units();
    let speed = u.cell_size_m / u.timestep_s;
    let jscale = u.mass_per_density_cell_kg * speed;
    let lscale = jscale * u.cell_size_m;
    let mut impulse = Vec3::ZERO;
    let mut moment = Vec3::ZERO;
    let mut loss = 0.0;
    for z in 0..N {
        for y in 0..N {
            for x in 0..N {
                let (law, tangent) = sample([x, y, z]);
                let imposed = tangent.scale(result.midpoint_angular_speed_rad_s / speed);
                let (r0, p0) = cell(motion.fluid(), x, y, z);
                let (r1, p1) = cell(result.motion.fluid(), x, y, z);
                let s0 = p0.scale(1.0 / r0).sub(imposed);
                let s1 = p1.scale(1.0 / r1).sub(imposed);
                let j = s0
                    .scale(law.darcy + law.forchheimer * norm(s0))
                    .add(s1.scale(law.darcy + law.forchheimer * norm(s1)))
                    .scale(0.5 * jscale);
                impulse = impulse.add(j);
                moment = moment.add(position([x, y, z]).sub(COM).cross(j));
                loss += j.dot(s0.add(s1).scale(0.5 * speed));
            }
        }
    }
    close_v(
        v(result.exchange.skeleton_impulse).scale(jscale),
        impulse,
        2e-12,
    );
    close_v(
        v(result.exchange.skeleton_angular_impulse).scale(lscale),
        moment,
        2e-12,
    );
    close(result.relative_dissipation_j, loss, 2e-12);
    close_v(result.guide_impulse_kg_m_s, impulse.scale(-1.0), 2e-12);
    close_v(
        result.guide_couple_kg_m2_s,
        AXIS.scale(moment.z).sub(moment),
        2e-12,
    );
    assert!(norm(result.guide_impulse_kg_m_s) > 1e-2);
    assert!(norm(result.guide_couple_kg_m2_s) > 1e-4);
    assert!(result.exchange.walls.measured_links > 0);
    assert!(norm(v(result.exchange.walls.fluid_population_impulse)) > 1e-2);
    assert_eq!(result.exchange.open, Default::default());
    assert_eq!(result.exchange.periodic, Default::default());
    close(result.exchange.walls.wall_work, 0.0, 0.0);
    let shift = |m: Vec3, j: Vec3| m.add(COM.sub(c.flow.balance_origin_m).cross(j));
    let p0 = physical(motion.fluid(), c.flow.balance_origin_m);
    let p1 = physical(result.motion.fluid(), c.flow.balance_origin_m);
    let ext = v(result.exchange.external_impulse).scale(jscale);
    let wall = v(result.exchange.walls.fluid_population_impulse).scale(jscale);
    close_v(
        p1.0.sub(p0.0),
        ext.add(wall).add(result.guide_impulse_kg_m_s),
        3e-10,
    );
    let world_spin = |b: RigidBodyState| {
        b.pose()
            .orientation()
            .rotate_body_to_world(b.angular_momentum_body())
    };
    let supplied = shift(
        v(result.exchange.external_angular_impulse).scale(lscale),
        ext,
    )
    .add(shift(
        v(result.exchange.walls.fluid_population_angular_impulse).scale(lscale),
        wall,
    ))
    .add(shift(
        result.guide_couple_kg_m2_s,
        result.guide_impulse_kg_m_s,
    ))
    .add(AXIS.scale(drive * u.timestep_s));
    close_v(
        p1.1.sub(p0.1)
            .add(world_spin(result.motion.body()).sub(world_spin(motion.body()))),
        supplied,
        3e-10,
    );
    close(
        result.other_fluid_kinetic_change_j,
        p1.2 - p0.2 - result.fluid_external_work_j - result.fluid_drag_work_j,
        3e-12,
    );
    // An actual native replay at the solved angular speed reproduces populations.
    let mut replay = motion.fluid().clone();
    let links = replay.boundary_links();
    let receipt = replay
        .step_with_darcy_forchheimer_and_wall_momentum(
            |cell| {
                let (law, tangent) = sample(cell);
                DarcyForchheimer3 {
                    linear: LinearResistance3 {
                        principal: [law.darcy; 3],
                        skeleton_velocity: a(
                            tangent.scale(result.midpoint_angular_speed_rad_s / speed)
                        ),
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
    assert_eq!(snapshot(&replay), snapshot(result.motion.fluid()));
    assert_eq!(receipt, result.exchange);
    let shifted = motion
        .step_candidate(
            coefficient,
            forcing,
            drive,
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
        snapshot(result.motion.fluid())
    );
    assert_eq!(shifted.motion.body(), result.motion.body());
}

#[test]
fn zero_drag_retains_native_finite_inertia_and_exact_axial_drive_phase() {
    let mut motion = fixture(0.02, false);
    let drive = 0.04;
    for step in 1..=4 {
        let result = motion
            .step_candidate(
                |_, _| PorousRotorResistance::default(),
                |_| [0.0; 3],
                drive,
                controls(),
                &CancelGate::new(),
            )
            .unwrap();
        motion = result.motion;
        let t = step as f64 * units().timestep_s;
        let expected_speed = 0.02 + drive * t / INERTIA;
        let expected_angle = 0.02 * t + 0.5 * drive * t * t / INERTIA;
        close(motion.angular_speed_rad_s().unwrap(), expected_speed, 2e-13);
        let actual = motion
            .body()
            .pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
        close_v(
            actual,
            Vec3::new(expected_angle.cos(), expected_angle.sin(), 0.0),
            2e-13,
        );
        close_v(motion.body().pose().position_world(), COM, 0.0);
        close_v(result.guide_couple_kg_m2_s, Vec3::ZERO, 0.0);
        close(result.relative_dissipation_j, 0.0, 0.0);
    }
}

#[test]
fn oblique_support_and_failed_candidates_preserve_the_native_inputs() {
    let orientation =
        UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), 0.6_f64.atan2(0.8)).unwrap();
    let (state, properties) = body(0.03, orientation);
    let axis = orientation.rotate_body_to_world(AXIS);
    let geom = PorousRotorCylinder {
        axis_world: axis,
        ..cylinder()
    };
    let admitted = AxisymmetricPorousRotor::new(
        grid(false),
        state,
        properties,
        units(),
        geom,
        CELLS,
        &CancelGate::new(),
    );
    let admitted = admitted.unwrap();
    let oblique = admitted
        .step_candidate(
            coefficient,
            |_| [0.0; 3],
            0.0,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
    assert!(oblique.motion.angular_speed_rad_s().unwrap() < 0.029);
    close_v(
        oblique
            .motion
            .body()
            .pose()
            .orientation()
            .rotate_body_to_world(AXIS),
        axis,
        2e-13,
    );
    let origin = controls().flow.balance_origin_m;
    let f0 = physical(admitted.fluid(), origin);
    let f1 = physical(oblique.motion.fluid(), origin);
    let scale = units().mass_per_density_cell_kg * units().cell_size_m / units().timestep_s;
    let drag = v(oblique.exchange.skeleton_impulse).scale(scale);
    let drag_moment =
        v(oblique.exchange.skeleton_angular_impulse).scale(scale * units().cell_size_m);
    let guide = drag.scale(-1.0);
    let guide_couple = axis.scale(axis.dot(drag_moment)).sub(drag_moment);
    let wall = v(oblique.exchange.walls.fluid_population_impulse).scale(scale);
    let wall_moment = v(oblique.exchange.walls.fluid_population_angular_impulse)
        .scale(scale * units().cell_size_m);
    close_v(f1.0.sub(f0.0), guide.add(wall), 3e-10);
    let final_body = oblique.motion.body();
    let spin = |b: RigidBodyState| {
        b.pose()
            .orientation()
            .rotate_body_to_world(b.angular_momentum_body())
    };
    close_v(
        f1.1.sub(f0.1).add(spin(final_body).sub(spin(state))),
        wall_moment
            .add(guide_couple)
            .add(COM.sub(origin).cross(guide.add(wall))),
        3e-10,
    );
    let kinetic =
        final_body.kinetic_energy(properties).unwrap() - state.kinetic_energy(properties).unwrap();
    close(
        kinetic + oblique.fluid_drag_work_j + oblique.relative_dissipation_j,
        0.0,
        2e-12,
    );
    close(oblique.guide_work_j, 0.0, 1e-14);
    let rejected = AxisymmetricPorousRotor::new(
        grid(false),
        state,
        properties,
        units(),
        PorousRotorCylinder {
            radius_m: 0.9,
            ..geom
        },
        CELLS,
        &CancelGate::new(),
    );
    assert!(matches!(rejected, Err(PorousMomentumError::Invalid(_))));
    let motion = fixture(0.04, false);
    let before = snapshot(motion.fluid());
    let old = motion.body();
    let c = controls();
    let budget = PorousRotorControls {
        flow: PorousMomentumControls {
            maximum_native_trials: 2,
            ..c.flow
        },
        ..c
    };
    assert!(matches!(
        motion.step_candidate(coefficient, |_| [0.0; 3], 0.0, budget, &CancelGate::new()),
        Err(PorousMomentumError::Budget { .. })
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        motion.step_candidate(coefficient, |_| [0.0; 3], 0.0, c, &cancelled),
        Err(PorousMomentumError::Cancelled)
    ));
    // Midpoint remains admissible, but the native driven endpoint exceeds its
    // physical rim-speed ceiling; the final acceptance must still refuse.
    assert!(matches!(
        motion.step_candidate(
            |_, _| PorousRotorResistance::default(),
            |_| [0.0; 3],
            0.5,
            c,
            &CancelGate::new()
        ),
        Err(PorousMomentumError::Invalid(_))
    ));
    assert_eq!(snapshot(motion.fluid()), before);
    assert_eq!(motion.body(), old);
    assert_eq!(motion.accepted_steps(), 0);
}
