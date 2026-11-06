//! Native axial motion and original finite-cell heat, with analytical fixtures.
//! The declared coefficients are synthetic tests, not material-source claims.
use fs_contact::tangential::impact::multiple::CoupledFrictionImpulseControls;
use fs_contact::tangential::impact::multiple::relaxation::CoupledFrictionFlowControls;
use fs_couple::friction_heat::joints::*;
use fs_couple::friction_heat::{FrictionHeatRoute, FrictionHeatTolerances};
use fs_exec::CancelGate;
use fs_mbd::rigid_rattle::distance::DistanceStepControls;
use fs_mbd::rigid_rattle::distance::graph::joint_friction::relaxation::{
    JointFrictionFlowControls, JointFrictionFlowMotion, JointFrictionFlowPort,
};
use fs_mbd::rigid_rattle::distance::graph::{
    ArticulatedJoint, DistanceGraphControls, PrismaticJoint, RevoluteJoint,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
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
fn controls() -> JointFrictionHeatControls {
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

use fs_couple::friction_heat::joints::trajectory::{
    JointFrictionHeatTrajectoryControls, advance_joint_friction_heat_trajectory,
};
use fs_mbd::{Gravity, RigidBodyIntegrator, Wrench};

fn stiction_source(
    link: u64,
    joint: usize,
    static_capacity: f64,
    kinetic_capacity: f64,
    drive: f64,
) -> JointFrictionHeatSource<()> {
    let mut source = one_source(link, joint, kinetic_capacity, drive, 1.0);
    source.interfaces[0].friction = JointFrictionFlowPort::PrismaticStiction {
        joint_index: joint,
        static_capacity_n: static_capacity,
        kinetic_capacity_n: kinetic_capacity,
        drive_n: drive,
    };
    source
}

#[test]
fn identical_load_holds_a_resting_joint_but_drives_a_sliding_joint_and_its_heat() {
    for (speed, final_speed, loss, work) in [(0.0, 0.0, 0.0, 0.0), (1.0, 1.5, 1.25, 1.875)] {
        let (states, p, joints) = slider(speed);
        let thermal = thermal(&[300.0; 2]);
        let original = thermal.cells().clone();
        let route = FrictionHeatRoute::from_link(&thermal, 10, [1.0; 2]).unwrap();
        let out = propose_joint_friction_heat(
            &states,
            &p,
            &joints,
            &thermal,
            0.0,
            1.0,
            &[JointFrictionHeatRoute {
                joint_index: 0,
                route: &route,
                temperature_guess_k: 300.0,
            }],
            controls(),
            &[0.0; 2],
            |_, links, indices, _| {
                Ok::<_, &'static str>(stiction_source(links[0], indices[0], 2.0, 1.0, 1.5))
            },
            &CancelGate::new(),
        )
        .unwrap();
        let JointFrictionFlowMotion::Prismatic {
            velocity_after_m_s, ..
        } = out.mechanical.ports[0].motion
        else {
            panic!()
        };
        near(velocity_after_m_s, final_speed, 2e-10);
        near(out.mechanical.dissipated_j, loss, 2e-10);
        near(out.heat.dissipated_j, loss, 2e-8);
        near(out.mechanical.drive_work_j, work, 2e-10);
        assert!(out.mechanical.stiction.is_some());
        for i in 0..2 {
            near(
                out.thermal.temperature_k(CellId(i)).unwrap(),
                300.0 + 0.5 * loss,
                2e-8,
            );
        }
        near(out.combined_energy_residual_j, 0.0, 2e-8);
        assert_eq!(thermal.cells(), &original);
        assert_eq!(thermal.time_s(), 0.0);
    }
}

#[test]
fn released_joint_resolves_temperature_dependent_sliding_loss_with_finite_heat() {
    // W=1, v0=0, drive=2.5, dt=.5. Cs<drive throughout. Kinetic C=.8-.1*(Ti-300).
    // D=C*(2.5-C)/8; original equal cells/conductances give Ti-300=D.
    // Thus C^2-82.5C+64=0; choose the source-admissible small root without cancellation.
    let (states, p, joints) = slider(0.0);
    let thermal = thermal(&[300.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, 10, [1.0; 2]).unwrap();
    let out = propose_joint_friction_heat(
        &states,
        &p,
        &joints,
        &thermal,
        0.0,
        0.5,
        &[JointFrictionHeatRoute {
            joint_index: 0,
            route: &route,
            temperature_guess_k: 300.0,
        }],
        controls(),
        &[0.0; 2],
        |_, links, indices, t| {
            Ok::<_, &'static str>(stiction_source(
                links[0],
                indices[0],
                1.8 - 0.02 * (t[0] - 300.0),
                0.8 - 0.1 * (t[0] - 300.0),
                2.5,
            ))
        },
        &CancelGate::new(),
    )
    .unwrap();
    let coefficient = 128.0 / (82.5 + 6550.25_f64.sqrt());
    let loss = coefficient * (2.5 - coefficient) / 8.0;
    let JointFrictionFlowPort::PrismaticStiction {
        kinetic_capacity_n, ..
    } = out.interfaces[0].source.friction
    else {
        panic!()
    };
    let JointFrictionFlowMotion::Prismatic {
        velocity_after_m_s, ..
    } = out.mechanical.ports[0].motion
    else {
        panic!()
    };
    near(kinetic_capacity_n, coefficient, 2e-9);
    near(velocity_after_m_s, 0.5 * (2.5 - coefficient), 2e-9);
    near(out.mechanical.dissipated_j, loss, 2e-9);
    near(
        out.heat.interfaces[0].interface_temperature_k,
        300.0 + loss,
        2e-8,
    );
    near(out.combined_energy_residual_j, 0.0, 2e-8);
    assert!(out.mechanical.stiction.is_some());
}

#[test]
fn heated_motion_holds_then_breaks_away_and_continues_below_the_static_limit() {
    // First half-second: rest under1.5N<Cs. Next:2.5N starts motion, v=.75.
    // Final:1.5N<Cs but>Ck sustains sliding to v=1, total q=.625 and D=.625.
    let (states, p, joints) = slider(0.0);
    let thermal = thermal(&[300.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, 10, [1.0; 2]).unwrap();
    let out = advance_joint_friction_heat_trajectory(
        &states,
        &p,
        &joints,
        RigidBodyIntegrator::new(Gravity::ZERO),
        &[Wrench::ZERO; 2],
        &thermal,
        &[0.5, 1.0, 1.5],
        &[JointFrictionHeatRoute {
            joint_index: 0,
            route: &route,
            temperature_guess_k: 300.0,
        }],
        JointFrictionHeatTrajectoryControls {
            heat: controls(),
            maximum_step_s: 0.5,
            maximum_steps: 8,
            maximum_thermal_trials: 2000,
            energy_residual_tolerance_j: 2e-8,
            absolute_energy_residual_tolerance_j: 4e-8,
            constraint_work_tolerance_j: 1e-8,
        },
        |time, _, _, links, indices, _| {
            let drive = if time > 0.5 && time < 1.0 { 2.5 } else { 1.5 };
            Ok::<_, &'static str>(stiction_source(links[0], indices[0], 2.0, 1.0, drive))
        },
        &CancelGate::new(),
    )
    .unwrap();
    let a = out.states_after[0];
    let b = out.states_after[1];
    near(
        a.pose().position_world().z - b.pose().position_world().z,
        0.625,
        2e-10,
    );
    near(
        a.center_of_mass_velocity_world(p[0]).unwrap().z
            - b.center_of_mass_velocity_world(p[1]).unwrap().z,
        1.0,
        2e-10,
    );
    near(out.mechanical_loss_j, 0.625, 2e-10);
    near(out.heat_j, 0.625, 2e-8);
    near(out.drive_work_j, 1.125, 2e-10);
    near(out.mechanical_energy_change_j, 0.5, 2e-10);
    near(out.energy_residual_j, 0.0, 2e-8);
    assert_eq!(out.steps, 3);
    for i in 0..2 {
        near(
            out.thermal.temperature_k(CellId(i)).unwrap(),
            300.3125,
            2e-8,
        );
    }
    assert_eq!(thermal.time_s(), 0.0);
}
