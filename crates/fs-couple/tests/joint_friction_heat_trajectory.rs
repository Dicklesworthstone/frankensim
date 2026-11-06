//! Actual native poses and finite-cell heating along an authored physical time grid.
//! The declared coefficients are synthetic tests, not material-source claims.
use fs_contact::tangential::impact::multiple::CoupledFrictionImpulseControls;
use fs_contact::tangential::impact::multiple::relaxation::CoupledFrictionFlowControls;
use fs_couple::friction_heat::joints::trajectory::*;
use fs_couple::friction_heat::joints::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_exec::CancelGate;
use fs_mbd::rigid_rattle::distance::DistanceStepControls;
use fs_mbd::rigid_rattle::distance::graph::joint_friction::relaxation::{
    JointFrictionFlowControls, JointFrictionFlowPort,
};
use fs_mbd::rigid_rattle::distance::graph::{
    ArticulatedJoint, DistanceGraphControls, PrismaticJoint, RevoluteJoint,
};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::backward_euler::StepCoefficients;
use fs_thermal::enthalpy::*;
use std::cell::Cell;
use std::collections::BTreeMap;

const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e} within {tolerance:e}"
    );
}
struct Caloric;
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h;
        if t.is_finite() && t > 0.0 {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        Ok(Some((self.temperature_k(b)? / self.temperature_k(a)?).ln()))
    }
}
static CALORIC: Caloric = Caloric;

/// Unit-capacity, unit-length cells. Each half-cell conductance is 2 W/K.
fn thermal(temperatures: &[f64]) -> DepositionSimulation<'static> {
    let n = temperatures.len();
    let mut events: Vec<_> = temperatures
        .iter()
        .enumerate()
        .map(|(i, &temperature)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: i as u64,
                material_id: MaterialId(0),
                mass_kg: 1.0,
                specific_h_j_kg: temperature - 300.0,
                geometry: CellGeometry {
                    start_m: [i as f64 - 1.0, 0.0, 0.0],
                    end_m: [i as f64, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    events.extend((0..n - 1).map(|i| TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 10 * (i as u64 + 1),
            a: CellId(i as u64),
            b: CellId(i as u64 + 1),
            area_m2: 1.0,
            conductance_w_k: 1.0,
        }),
    }));
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), &CALORIC as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: n,
            max_links: n - 1,
            max_events: 2 * n - 1,
            max_substeps: 16,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-10,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn phase_controls() -> JointFrictionHeatControls {
    JointFrictionHeatControls {
        mechanical: JointFrictionFlowControls {
            graph: DistanceGraphControls {
                step: DistanceStepControls {
                    maximum_impulse_ns: 1000.0,
                    maximum_rotation_increment_rad: 0.25,
                    maximum_position_iterations: 30,
                    maximum_backtracks: 16,
                },
                position_tolerance_m: 1e-10,
                velocity_tolerance_mps: 1e-10,
                minimum_scaled_pivot: 1e-10,
                maximum_bodies: 3,
                maximum_constraints: 16,
                maximum_trial_evaluations: 160,
                maximum_live_bytes: 2_000_000,
            },
            flow: CoupledFrictionFlowControls {
                projection: CoupledFrictionImpulseControls {
                    maximum_iterations: 256,
                    impulse_tolerance_ns: 1e-10,
                    velocity_tolerance_m_s: 1e-10,
                    energy_tolerance_j: 1e-9,
                    relative_matrix_tolerance: 1e-12,
                    maximum_impulse_ns: 1000.0,
                },
                maximum_segments: 32,
                acceleration_tolerance_m_s2: 1e-10,
                maximum_force_n: 1000.0,
            },
            work_tolerance_j: 1e-9,
            linear_impulse_tolerance_ns: 1e-9,
            angular_impulse_tolerance_nms: 1e-9,
        },
        maximum_cells: 3,
        maximum_unknowns: 5,
        maximum_trials: 512,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-5,
        temperature_scale_j_per_k: 1.0,
        heat: FrictionHeatTolerances {
            cell_j: 2e-9,
            thermal_j: 2e-9,
            interface_j: 2e-9,
            total_j: 2e-9,
        },
        combined_energy_tolerance_j: 4e-9,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-11,
            linear_restart: 5,
            forcing_minimum: 1e-6,
            ..Default::default()
        },
        audit_entropy: true,
    }
}
fn properties(mass: f64, inertia: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap()
}
fn state(p: MassProperties, velocity: f64, spin: f64) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(Vec3::ZERO, UnitQuaternion::IDENTITY).unwrap(),
        Z.scale(velocity * p.mass()),
        p.angular_momentum_body_checked(Z.scale(spin)).unwrap(),
    )
    .unwrap()
}
fn slider(
    speed: f64,
) -> (
    [RigidBodyState; 2],
    [MassProperties; 2],
    [ArticulatedJoint; 1],
) {
    let p = [properties(2.0, 2.0); 2];
    (
        [
            state(p[0], 0.5 * speed, 0.0),
            state(p[1], -0.5 * speed, 0.0),
        ],
        p,
        [ArticulatedJoint::Prismatic(
            PrismaticJoint::new([0, 1], [Vec3::ZERO; 2], [UnitQuaternion::IDENTITY; 2], 1.0)
                .unwrap(),
        )],
    )
}
fn one_source(
    link: u64,
    joint: usize,
    capacity: f64,
    drive: f64,
    fraction: f64,
) -> JointFrictionHeatSource<()> {
    JointFrictionHeatSource {
        interfaces: vec![JointFrictionHeatInterfaceSource {
            link_id: link,
            friction: JointFrictionFlowPort::Prismatic {
                joint_index: joint,
                capacity_n: capacity,
                drive_n: drive,
            },
            conductivity_w_m_k: [1.0; 2],
            heat_fraction: fraction,
        }],
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        },
    }
}

fn controls() -> JointFrictionHeatTrajectoryControls {
    JointFrictionHeatTrajectoryControls {
        heat: phase_controls(),
        maximum_step_s: 1.0,
        maximum_steps: 64,
        maximum_thermal_trials: 4096,
        energy_residual_tolerance_j: 2e-7,
        absolute_energy_residual_tolerance_j: 2e-7,
        constraint_work_tolerance_j: 1e-8,
    }
}
fn free() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::ZERO)
}
fn relative_position(states: &[RigidBodyState]) -> f64 {
    states[0].pose().position_world().z - states[1].pose().position_world().z
}
fn relative_velocity(states: &[RigidBodyState], properties: &[MassProperties]) -> f64 {
    states[0]
        .center_of_mass_velocity_world(properties[0])
        .unwrap()
        .z
        - states[1]
            .center_of_mass_velocity_world(properties[1])
            .unwrap()
            .z
}

#[test]
fn stopping_and_reversing_motion_advance_real_poses_and_integrated_heat() {
    // The authored grids resolve the true zero-speed events. Each remaining
    // affine velocity segment therefore has exact trapezoidal native pose travel.
    for (drive, grid, endpoint, travel, loss) in [
        (0.0, vec![0.25, 0.5, 1.0, 2.0], 0.0, 0.5, 0.5),
        (
            -2.0,
            vec![1.0 / 3.0, 1.0],
            -2.0 / 3.0,
            -1.0 / 18.0,
            7.0 / 18.0,
        ),
    ] {
        let (states, p, joints) = slider(1.0);
        let old_states = states;
        let thermal = thermal(&[300.0; 2]);
        let old_cells = thermal.cells().clone();
        let route = FrictionHeatRoute::from_link(&thermal, 10, [1.0; 2]).unwrap();
        let out = advance_joint_friction_heat_trajectory(
            &states,
            &p,
            &joints,
            free(),
            &[Wrench::ZERO; 2],
            &thermal,
            &grid,
            &[JointFrictionHeatRoute {
                joint_index: 0,
                route: &route,
                temperature_guess_k: 300.0,
            }],
            controls(),
            |_, _, _, links, indices, _| {
                Ok::<_, &'static str>(one_source(links[0], indices[0], 1.0, drive, 1.0))
            },
            &CancelGate::new(),
        )
        .unwrap();
        near(relative_position(&out.states_after), travel, 2e-10);
        near(relative_velocity(&out.states_after, &p), endpoint, 2e-10);
        for i in 0..2 {
            near(
                out.thermal.temperature_k(CellId(i)).unwrap(),
                300.0 + 0.5 * loss,
                2e-8,
            );
        }
        near(out.mechanical_loss_j, loss, 2e-10);
        near(out.heat_j, loss, 2e-10);
        near(
            out.mechanical_energy_change_j,
            0.5 * (endpoint * endpoint - 1.0),
            2e-10,
        );
        near(out.energy_residual_j, 0.0, 2e-8);
        assert!(out.absolute_step_energy_residual_j < 2e-8);
        assert!(out.absolute_constraint_work_j < 1e-9);
        assert!(out.maximum_constraint_gap_m < 1e-10);
        assert_eq!(out.steps, grid.len());
        assert_eq!(out.thermal.time_s(), *grid.last().unwrap());
        assert_eq!(states, old_states);
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(thermal.time_s(), 0.0);
    }
}

#[test]
fn finite_temperature_feedback_changes_the_evolving_pose_trajectory() {
    let (states, p, joints) = slider(1.0);
    let thermal = thermal(&[300.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, 10, [1.0; 2]).unwrap();
    let grid: Vec<_> = (1..=8).map(|i| f64::from(i) / 8.0).collect();
    let out = advance_joint_friction_heat_trajectory(
        &states,
        &p,
        &joints,
        free(),
        &[Wrench::ZERO; 2],
        &thermal,
        &grid,
        &[JointFrictionHeatRoute {
            joint_index: 0,
            route: &route,
            temperature_guess_k: 300.0,
        }],
        controls(),
        |_, _, _, links, indices, temperatures| {
            Ok::<_, &'static str>(one_source(
                links[0],
                indices[0],
                1.0 - 0.5 * (temperatures[0] - 300.0),
                0.0,
                1.0,
            ))
        },
        &CancelGate::new(),
    )
    .unwrap();
    // Independent per-interval scalar solution. With current relative speed v,
    // accumulated heat Q and h=1/8: D=C*v*h-C^2*h^2/2,
    // C=1-Q/4-(1/4+1/(8h))*D. Solve its small positive quadratic root.
    let h = 0.125;
    let beta = 0.25 + 1.0 / (8.0 * h);
    let mut speed = 1.0;
    let mut accumulated_heat = 0.0;
    let mut position = 0.0;
    let mut last_capacity = 0.0;
    for _ in &grid {
        let a = 1.0 - accumulated_heat / 4.0;
        let b: f64 = 1.0 + beta * speed * h;
        let capacity = 2.0 * a / (b + (b * b - 2.0 * beta * h * h * a).sqrt());
        let next = speed - capacity * h;
        assert!(next > 0.0);
        accumulated_heat += capacity * h * (speed - 0.5 * capacity * h);
        position += 0.5 * (speed + next) * h;
        speed = next;
        last_capacity = capacity;
    }
    near(relative_position(&out.states_after), position, 2e-8);
    near(relative_velocity(&out.states_after, &p), speed, 2e-8);
    near(out.heat_j, accumulated_heat, 2e-8);
    near(out.thermal_energy_change_j, accumulated_heat, 2e-8);
    near(
        1.0 - 0.5 * (out.last_heat.interfaces[0].interface_temperature_k - 300.0),
        last_capacity,
        2e-8,
    );
    // With the initial coefficient held at one, this slider stops at t=1 and
    // travels exactly 1/2 m. The resolved thermal source gives different motion.
    assert!(speed > 0.1);
    assert!(position > 0.5);
    assert!(out.thermal_trials > grid.len());
    near(out.energy_residual_j, 0.0, 2e-8);
}

#[test]
fn authored_drive_switch_keeps_wrench_gravity_heat_and_exports_separate() {
    let (states, p, joints) = slider(1.0);
    let thermal = thermal(&[300.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, 10, [1.0; 2]).unwrap();
    let integrator = RigidBodyIntegrator::new(Gravity::new(Vec3::new(0.0, -0.5, 0.0)).unwrap());
    let wrenches = [Wrench {
        force_world: Vec3::new(2.0, 0.0, 0.0),
        torque_body: Vec3::ZERO,
    }; 2];
    let saw_first = Cell::new(false);
    let saw_last = Cell::new(false);
    let out = advance_joint_friction_heat_trajectory(
        &states,
        &p,
        &joints,
        integrator,
        &wrenches,
        &thermal,
        &[0.5, 2.0 / 3.0, 1.0],
        &[JointFrictionHeatRoute {
            joint_index: 0,
            route: &route,
            temperature_guess_k: 300.0,
        }],
        controls(),
        |time, midpoint, _, links, indices, _| {
            // The source sees the actually moved native bodies, including the
            // common transverse force and gravity, at the physical midpoint.
            near(
                midpoint[0].pose().position_world().x,
                0.5 * time * time,
                2e-10,
            );
            near(
                midpoint[0].pose().position_world().y,
                -0.25 * time * time,
                2e-10,
            );
            let drive = if time < 0.5 {
                saw_first.set(true);
                0.0
            } else {
                saw_last.set(true);
                -2.0
            };
            let mut source = one_source(links[0], indices[0], 1.0, drive, 0.5);
            source.transport.sources_w.insert(CellId(0), 0.2);
            Ok::<_, &'static str>(source)
        },
        &CancelGate::new(),
    )
    .unwrap();
    assert!(saw_first.get() && saw_last.get());
    near(relative_position(&out.states_after), 13.0 / 36.0, 2e-10);
    near(relative_velocity(&out.states_after, &p), -1.0 / 3.0, 2e-10);
    near(out.wrench_work_j, 2.0, 2e-10);
    near(out.drive_work_j, 1.0 / 36.0, 2e-10);
    near(out.mechanical_loss_j, 17.0 / 36.0, 2e-10);
    near(out.heat_j, 17.0 / 72.0, 2e-10);
    near(out.nonthermal_export_j, 17.0 / 72.0, 2e-10);
    near(out.external_heat_j, 0.2, 2e-10);
    near(out.thermal_energy_change_j, 17.0 / 72.0 + 0.2, 2e-8);
    near(out.mechanical_energy_change_j, 2.0 - 4.0 / 9.0, 2e-10);
    near(out.ports[0].mechanical_loss_j, out.mechanical_loss_j, 0.0);
    near(out.energy_residual_j, 0.0, 2e-8);
    for state in &out.states_after {
        near(state.pose().position_world().x, 0.5, 2e-10);
        near(state.pose().position_world().y, -0.25, 2e-10);
    }
}

#[test]
fn late_source_failure_or_cancellation_preserves_the_entire_accepted_trajectory() {
    let (states, p, joints) = slider(1.0);
    let old_states = states;
    let thermal = thermal(&[300.0; 2]);
    let old_cells = thermal.cells().clone();
    let route = FrictionHeatRoute::from_link(&thermal, 10, [1.0; 2]).unwrap();
    let bindings = [JointFrictionHeatRoute {
        joint_index: 0,
        route: &route,
        temperature_guess_k: 300.0,
    }];
    for cancel in [false, true] {
        let gate = CancelGate::new();
        let latest = Cell::new(0.0_f64);
        let out = advance_joint_friction_heat_trajectory(
            &states,
            &p,
            &joints,
            free(),
            &[Wrench::ZERO; 2],
            &thermal,
            &[0.25, 0.5, 0.75, 1.0],
            &bindings,
            controls(),
            |time, _, _, links, indices, _| {
                latest.set(latest.get().max(time));
                if time > 0.5 {
                    if cancel {
                        gate.request();
                    } else {
                        return Err("late source domain");
                    }
                }
                Ok(one_source(links[0], indices[0], 0.25, 0.0, 1.0))
            },
            &gate,
        );
        assert!(out.is_err());
        if !cancel {
            assert!(matches!(
                out,
                Err(JointFrictionHeatTrajectoryError::Heat(
                    JointFrictionHeatError::Source("late source domain")
                ))
            ));
        }
        assert!(latest.get() > 0.5);
        assert_eq!(states, old_states);
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(thermal.time_s(), 0.0);
    }
    let calls = Cell::new(0);
    let invalid_grid = advance_joint_friction_heat_trajectory(
        &states,
        &p,
        &joints,
        free(),
        &[Wrench::ZERO; 2],
        &thermal,
        &[0.25, 0.25],
        &bindings,
        controls(),
        |_, _, _, links, indices, _| {
            calls.set(calls.get() + 1);
            Ok::<_, &'static str>(one_source(links[0], indices[0], 1.0, 0.0, 1.0))
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        invalid_grid,
        Err(JointFrictionHeatTrajectoryError::Invalid(_))
    ));
    assert_eq!(calls.get(), 0);
}
