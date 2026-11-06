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

#[test]
fn reversing_slider_heats_by_the_actual_continuous_path_loss() {
    let (states, p, joints) = slider(1.0);
    let old_states = states;
    let thermal = thermal(&[300.0; 2]);
    let old_cells = thermal.cells().clone();
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
            Ok::<_, &'static str>(one_source(links[0], indices[0], 1.0, -2.0, 1.0))
        },
        &CancelGate::new(),
    )
    .unwrap();
    let JointFrictionFlowMotion::Prismatic {
        velocity_before_m_s,
        velocity_after_m_s,
        friction_impulse_ns,
        ..
    } = out.mechanical.ports[0].motion
    else {
        panic!("slider motion");
    };
    near(velocity_after_m_s, -2.0 / 3.0, 2e-11);
    // The full-step endpoint pairing has the wrong sign for friction heating.
    near(
        0.5 * (velocity_before_m_s + velocity_after_m_s) * friction_impulse_ns,
        1.0 / 18.0,
        2e-11,
    );
    near(out.mechanical.dissipated_j, 7.0 / 18.0, 2e-11);
    near(out.heat.dissipated_j, 7.0 / 18.0, 2e-11);
    near(out.mechanical.drive_work_j, 1.0 / 9.0, 2e-11);
    near(out.mechanical.kinetic_energy_change_j, -5.0 / 18.0, 2e-11);
    for i in 0..2 {
        near(
            out.thermal.temperature_k(CellId(i)).unwrap(),
            300.0 + 7.0 / 36.0,
            2e-9,
        );
        near(
            out.heat.interfaces[0].body_heat_j[i as usize],
            7.0 / 36.0,
            2e-9,
        );
    }
    near(
        out.heat.interfaces[0].interface_temperature_k,
        300.0 + 7.0 / 24.0,
        2e-9,
    );
    near(out.nonthermal_export_j, 0.0, 0.0);
    near(out.combined_energy_residual_j, 0.0, 2e-9);
    let entropy = out.heat.entropy.unwrap();
    near(
        entropy.friction_generation_j_k,
        (7.0 / 18.0) / (300.0 + 7.0 / 24.0),
        2e-11,
    );
    assert!(entropy.interface_conduction_generation_j_k > 0.0);
    assert_eq!(states, old_states);
    assert_eq!(thermal.cells(), &old_cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(out.thermal.time_s(), 1.0);
}

#[test]
fn endpoint_temperature_changes_the_native_friction_trajectory() {
    // v0=1, W=1, dt=1/4, C=1-(T_interface-300)/2.
    // D=C/4-C^2/32 and T_interface-300=3D/2, giving 3C^2-152C+128=0.
    let (states, p, joints) = slider(1.0);
    let thermal = thermal(&[300.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, 10, [1.0; 2]).unwrap();
    let calls = Cell::new(0);
    let out = propose_joint_friction_heat(
        &states,
        &p,
        &joints,
        &thermal,
        0.0,
        0.25,
        &[JointFrictionHeatRoute {
            joint_index: 0,
            route: &route,
            temperature_guess_k: 300.0,
        }],
        controls(),
        &[0.0; 2],
        |_, links, indices, temperatures| {
            calls.set(calls.get() + 1);
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
    let expected_capacity = 256.0 / (152.0 + 21568.0_f64.sqrt());
    let expected_loss = expected_capacity / 4.0 - expected_capacity * expected_capacity / 32.0;
    let JointFrictionFlowPort::Prismatic { capacity_n, .. } = out.interfaces[0].source.friction
    else {
        panic!("slider source");
    };
    near(capacity_n, expected_capacity, 2e-9);
    near(out.mechanical.dissipated_j, expected_loss, 2e-9);
    near(
        out.heat.interfaces[0].interface_temperature_k,
        300.0 + 1.5 * expected_loss,
        2e-9,
    );
    let JointFrictionFlowMotion::Prismatic {
        velocity_after_m_s, ..
    } = out.mechanical.ports[0].motion
    else {
        panic!("slider motion");
    };
    near(velocity_after_m_s, 1.0 - expected_capacity / 4.0, 2e-9);
    assert!(calls.get() > 2);
    assert!(capacity_n < 0.9);
    near(out.combined_energy_residual_j, 0.0, 2e-9);
}

#[test]
fn explicit_nonthermal_exports_preserve_signed_hot_to_cold_exchange() {
    let (states, p, joints) = slider(1.0);
    for fraction in [0.0, 0.25, 1.0] {
        let thermal = thermal(&[310.0, 290.0]);
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
                Ok::<_, &'static str>(one_source(links[0], indices[0], 1.0, -2.0, fraction))
            },
            &CancelGate::new(),
        )
        .unwrap();
        let loss = 7.0 / 18.0;
        let heat = fraction * loss;
        let interface_temperature = 300.0 + 0.75 * heat;
        near(out.nonthermal_export_j, loss - heat, 2e-11);
        near(out.heat.dissipated_j, heat, 2e-11);
        near(
            out.heat.interfaces[0].interface_temperature_k,
            interface_temperature,
            2e-9,
        );
        for (i, old) in [310.0, 290.0].into_iter().enumerate() {
            near(
                out.thermal.temperature_k(CellId(i as u64)).unwrap(),
                (old + 2.0 * interface_temperature) / 3.0,
                2e-9,
            );
        }
        assert!(out.heat.interfaces[0].body_heat_j[0] < 0.0);
        assert!(out.heat.interfaces[0].body_heat_j[1] > 0.0);
        near(out.heat.interfaces[0].body_heat_j.iter().sum(), heat, 2e-9);
        near(out.loss_partition_residual_j, 0.0, 2e-11);
        near(out.combined_energy_residual_j, 0.0, 2e-9);
    }
}

#[test]
fn two_coupled_hinges_share_the_same_finite_middle_cell() {
    let p = [
        properties(1.0, 2.0),
        properties(2.0, 3.0),
        properties(3.0, 5.0),
    ];
    let states = std::array::from_fn::<_, 3, _>(|i| state(p[i], 0.0, [3.0, -1.0, 2.0][i]));
    let joints = std::array::from_fn::<_, 2, _>(|i| {
        ArticulatedJoint::Revolute(
            RevoluteJoint::new([i, i + 1], [Vec3::ZERO; 2], [Z; 2], [0.25, 2.0][i]).unwrap(),
        )
    });
    let thermal = thermal(&[300.0; 3]);
    let routes = [10, 20].map(|id| FrictionHeatRoute::from_link(&thermal, id, [1.0; 2]).unwrap());
    let bindings = [
        JointFrictionHeatRoute {
            joint_index: 1,
            route: &routes[1],
            temperature_guess_k: 300.0,
        },
        JointFrictionHeatRoute {
            joint_index: 0,
            route: &routes[0],
            temperature_guess_k: 300.0,
        },
    ];
    let source =
        |_: &DepositionSimulation<'_>, links: &[u64], indices: &[usize], temperatures: &[f64]| {
            assert_eq!(links, &[10, 20]);
            assert_eq!(indices, &[0, 1]);
            Ok::<_, &'static str>(JointFrictionHeatSource {
                interfaces: (0..2)
                    .map(|i| JointFrictionHeatInterfaceSource {
                        link_id: links[i],
                        friction: JointFrictionFlowPort::Revolute {
                            joint_index: indices[i],
                            capacity_nm: [1.0, 2.0][i],
                            drive_nm: 0.0,
                        },
                        conductivity_w_m_k: [1.0; 2],
                        heat_fraction: 1.0,
                    })
                    .collect(),
                transport: StepCoefficients {
                    sources_w: BTreeMap::new(),
                    conductances_w_k: BTreeMap::new(),
                    payload: temperatures.to_vec(),
                },
            })
        };
    let out = propose_joint_friction_heat(
        &states,
        &p,
        &joints,
        &thermal,
        0.0,
        5.0,
        &bindings,
        controls(),
        &[0.0; 3],
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let losses = [787.0 / 140.0, 45.0 / 7.0];
    let total = 241.0 / 20.0;
    let middle_rise = 11.0 * total / 32.0;
    let rises = [
        (0.5 * losses[0] + 5.0 * middle_rise) / 6.0,
        middle_rise,
        (0.5 * losses[1] + 5.0 * middle_rise) / 6.0,
    ];
    for i in 0..3 {
        near(
            out.mechanical.states_after[i]
                .point_kinematics(p[i], Vec3::ZERO)
                .unwrap()
                .angular_velocity_world
                .z,
            1.3,
            2e-11,
        );
        near(
            out.thermal.temperature_k(CellId(i as u64)).unwrap(),
            300.0 + rises[i],
            3e-9,
        );
    }
    for i in 0..2 {
        near(out.interfaces[i].mechanical_loss_j, losses[i], 2e-11);
        near(out.heat.interfaces[i].dissipated_j, losses[i], 2e-11);
        near(
            out.source[i],
            out.heat.interfaces[i].interface_temperature_k,
            0.0,
        );
    }
    near(out.heat.thermal.stored_energy_change_j, total, 3e-9);
    near(out.combined_energy_residual_j, 0.0, 3e-9);
    let reordered = propose_joint_friction_heat(
        &states,
        &p,
        &joints,
        &thermal,
        0.0,
        5.0,
        &[bindings[1], bindings[0]],
        controls(),
        &[0.0; 3],
        source,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(out.mechanical, reordered.mechanical);
    assert_eq!(out.thermal.cells(), reordered.thermal.cells());
}

#[test]
fn held_drive_has_no_friction_heat_when_the_joint_does_not_move() {
    let (states, p, joints) = slider(0.0);
    let thermal = thermal(&[300.0; 2]);
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
            Ok::<_, &'static str>(one_source(links[0], indices[0], 1.0, 0.5, 1.0))
        },
        &CancelGate::new(),
    )
    .unwrap();
    let JointFrictionFlowMotion::Prismatic {
        friction_impulse_ns,
        drive_impulse_ns,
        distance_m,
        ..
    } = out.mechanical.ports[0].motion
    else {
        panic!("slider motion");
    };
    near(friction_impulse_ns, -0.5, 2e-11);
    near(drive_impulse_ns, 0.5, 2e-11);
    near(distance_m, 0.0, 0.0);
    near(out.heat.dissipated_j, 0.0, 0.0);
    assert_eq!(out.mechanical.states_after, states);
    assert_eq!(out.thermal.cells(), thermal.cells());
}

#[test]
fn refusals_and_mid_trial_cancellation_publish_no_body_or_caloric_state() {
    let (states, p, joints) = slider(1.0);
    let accepted_states = states;
    let thermal = thermal(&[300.0; 2]);
    let accepted_cells = thermal.cells().clone();
    let route = FrictionHeatRoute::from_link(&thermal, 10, [1.0; 2]).unwrap();
    let bindings = [JointFrictionHeatRoute {
        joint_index: 0,
        route: &route,
        temperature_guess_k: 300.0,
    }];
    let source_error = propose_joint_friction_heat(
        &states,
        &p,
        &joints,
        &thermal,
        0.0,
        1.0,
        &bindings,
        controls(),
        &[0.0; 2],
        |_, _, _, _| Err::<JointFrictionHeatSource<()>, _>("source domain"),
        &CancelGate::new(),
    );
    assert!(matches!(
        source_error,
        Err(JointFrictionHeatError::Source("source domain"))
    ));
    for (link, fraction) in [(10, 1.1), (20, 1.0)] {
        let invalid = propose_joint_friction_heat(
            &states,
            &p,
            &joints,
            &thermal,
            0.0,
            1.0,
            &bindings,
            controls(),
            &[0.0; 2],
            |_, _, _, _| Ok::<_, &'static str>(one_source(link, 0, 1.0, -2.0, fraction)),
            &CancelGate::new(),
        );
        assert!(matches!(invalid, Err(JointFrictionHeatError::Invalid(_))));
    }
    let gate = CancelGate::new();
    let calls = Cell::new(0);
    let cancelled = propose_joint_friction_heat(
        &states,
        &p,
        &joints,
        &thermal,
        0.0,
        1.0,
        &bindings,
        controls(),
        &[0.0; 2],
        |_, links, indices, _| {
            calls.set(calls.get() + 1);
            gate.request();
            Ok::<_, &'static str>(one_source(links[0], indices[0], 1.0, -2.0, 1.0))
        },
        &gate,
    );
    assert!(cancelled.is_err());
    assert_eq!(calls.get(), 1);
    let mut limited = controls();
    limited.maximum_trials = 1;
    let exhausted = propose_joint_friction_heat(
        &states,
        &p,
        &joints,
        &thermal,
        0.0,
        1.0,
        &bindings,
        limited,
        &[0.0; 2],
        |_, links, indices, _| {
            Ok::<_, &'static str>(one_source(links[0], indices[0], 1.0, -2.0, 1.0))
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        exhausted,
        Err(JointFrictionHeatError::TrialBudget)
    ));
    assert_eq!(states, accepted_states);
    assert_eq!(thermal.cells(), &accepted_cells);
    assert_eq!(thermal.time_s(), 0.0);
}
