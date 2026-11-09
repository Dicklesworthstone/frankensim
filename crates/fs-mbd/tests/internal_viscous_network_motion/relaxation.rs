//! Corotated histories under finite native motion, original work and heat.

#[path = "relaxation/kinetics.rs"]
mod kinetics;
use super::*;
use fs_mbd::spatial_loads::internal_impulse::viscous::network::{
    driven::DrivenViscousNetworkLimits,
    relaxation::{RelaxationNetworkLimits, RelaxationNetworkLink, SpatialRelaxationBranch},
};
use fs_mbd::spatial_loads::internal_motion::viscous::{
    driven::HeldInternalEffort,
    network::relaxation::{
        RelaxationNetworkMotionControls, RelaxationNetworkMotionError, step_relaxation_network,
        step_relaxation_network_with_loads,
    },
};

fn memory_controls() -> RelaxationNetworkMotionControls {
    RelaxationNetworkMotionControls {
        source: controls().source,
        limits: RelaxationNetworkLimits {
            source: DrivenViscousNetworkLimits {
                maximum_bodies: 8,
                maximum_connectors: 8,
                maximum_live_bytes: 128_000_000,
                maximum_preparation_terms: 160_000_000_000,
            },
            maximum_branches: 8,
        },
    }
}
fn memory_links(viscous: [ViscousNetworkLink; 3]) -> [RelaxationNetworkLink; 3] {
    core::array::from_fn(|e| {
        let diagonal = [1.1, 0.9, 0.7, 0.3, 0.5, 0.4];
        let rows = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                if i == j {
                    diagonal[i]
                } else if j == (i + 1) % 6 {
                    0.12 * (e + 1) as f64
                } else {
                    0.0
                }
            })
        });
        RelaxationNetworkLink {
            viscous: viscous[e].clone(),
            branches: vec![
                SpatialRelaxationBranch::from_gram_rows(rows, [0.3, 0.6, 1.1][e]).unwrap(),
            ],
        }
    })
}
fn memory() -> Vec<Vec<[f64; 6]>> {
    (0..3)
        .map(|e| {
            vec![core::array::from_fn(|c| {
                0.025 * (e + 1) as f64 * (c + 1) as f64 * if c % 2 == 0 { 1.0 } else { -1.0 }
            })]
        })
        .collect()
}
fn drives() -> [[f64; 6]; 3] {
    [
        [0.2, -0.1, 0.3, 0.04, -0.03, 0.02],
        [-0.3, 0.2, 0.1, -0.02, 0.05, 0.01],
        [0.1, 0.15, -0.2, 0.02, 0.01, -0.04],
    ]
}

#[derive(Clone, Copy)]
struct HistoryReference {
    body: Reference,
    memory: [[f64; 6]; 3],
    work: [f64; 3],
}
impl HistoryReference {
    fn shifted(self, rhs: Self, h: f64) -> Self {
        Self {
            body: self.body.shifted(rhs.body, h),
            memory: core::array::from_fn(|e| {
                core::array::from_fn(|c| self.memory[e][c] + h * rhs.memory[e][c])
            }),
            work: core::array::from_fn(|e| self.work[e] + h * rhs.work[e]),
        }
    }
}
fn history_rhs(
    s: HistoryReference,
    properties: [MassProperties; 3],
    links: &[RelaxationNetworkLink; 3],
    viscous: &[ViscousNetworkLink; 3],
    efforts: &[[f64; 6]; 3],
) -> HistoryReference {
    // Parent reference already supplies the independent full-inertia native
    // Newton-Euler/Kelvin equations. Add original Maxwell efforts and histories.
    let mut out = HistoryReference {
        body: derivative(s.body, properties, viscous),
        memory: [[0.0; 6]; 3],
        work: [0.0; 3],
    };
    let q = s
        .body
        .q
        .map(|q| UnitQuaternion::new(q[0], q[1], q[2], q[3]).unwrap());
    let spin: [Vec3; 3] = core::array::from_fn(|i| {
        q[i].rotate_body_to_world(
            properties[i]
                .angular_velocity_body_checked(s.body.h[i])
                .unwrap(),
        )
    });
    let velocity: [Vec3; 3] =
        core::array::from_fn(|i| s.body.p[i].scale(1.0 / properties[i].mass()));
    for (e, source) in links.iter().enumerate() {
        let link = &source.viscous;
        let [a, b] = link.bodies;
        let rb = q[b].rotate_body_to_world(link.arms_body_m[1]);
        let arm = s.body.x[b].sub(s.body.x[a]).add(rb);
        let slip = q[a].rotate_world_to_body(
            velocity[b]
                .sub(velocity[a])
                .add(spin[b].cross(rb))
                .sub(spin[a].cross(arm)),
        );
        let angular = q[a].rotate_world_to_body(spin[b].sub(spin[a]));
        let twist = [slip.x, slip.y, slip.z, angular.x, angular.y, angular.z];
        let branch = &source.branches[0];
        let l = branch.gram_rows();
        let wrench: [f64; 6] = core::array::from_fn(|c| {
            efforts[e][c] - (0..6).map(|r| l[r][c] * s.memory[e][r]).sum::<f64>()
        });
        let force = q[a].rotate_body_to_world(Vec3::new(wrench[0], wrench[1], wrench[2]));
        let couple = q[a].rotate_body_to_world(Vec3::new(wrench[3], wrench[4], wrench[5]));
        out.body.p[a] = out.body.p[a].sub(force);
        out.body.p[b] = out.body.p[b].add(force);
        out.body.h[a] = out.body.h[a].sub(q[a].rotate_world_to_body(arm.cross(force).add(couple)));
        out.body.h[b] = out.body.h[b].add(q[b].rotate_world_to_body(rb.cross(force).add(couple)));
        for c in 0..6 {
            out.memory[e][c] = (0..6).map(|j| l[c][j] * twist[j]).sum::<f64>()
                - s.memory[e][c] / branch.relaxation_time_s();
            out.body.heat[e] += s.memory[e][c] * s.memory[e][c] / branch.relaxation_time_s();
            out.work[e] += efforts[e][c] * twist[c];
        }
    }
    out
}

#[test]
fn g1_complete_moving_memory_graph_refines_to_independent_corotated_newton_euler_ode() {
    let (initial, properties, viscous) = full();
    let links = memory_links(viscous.clone());
    let memory = memory();
    let efforts = drives();
    let duration = 0.06;
    let mut expected = HistoryReference {
        body: Reference::from_states(initial),
        memory: core::array::from_fn(|e| memory[e][0]),
        work: [0.0; 3],
    };
    let h = duration / 2048.0;
    for _ in 0..2048 {
        let k1 = history_rhs(expected, properties, &links, &viscous, &efforts);
        let k2 = history_rhs(
            expected.shifted(k1, 0.5 * h),
            properties,
            &links,
            &viscous,
            &efforts,
        );
        let k3 = history_rhs(
            expected.shifted(k2, 0.5 * h),
            properties,
            &links,
            &viscous,
            &efforts,
        );
        let k4 = history_rhs(
            expected.shifted(k3, h),
            properties,
            &links,
            &viscous,
            &efforts,
        );
        expected = expected
            .shifted(k1, h / 6.0)
            .shifted(k2, h / 3.0)
            .shifted(k3, h / 3.0)
            .shifted(k4, h / 6.0);
    }
    let mut errors = [0.0; 3];
    for (rung, steps) in [4, 8, 16].into_iter().enumerate() {
        let mut states = initial.to_vec();
        let mut history = memory.clone();
        let (mut heat, mut work) = ([0.0; 3], [0.0; 3]);
        for _ in 0..steps {
            let out = step_relaxation_network(
                &states,
                &properties,
                &links,
                &history,
                &efforts.map(HeldInternalEffort::MaterialFirst),
                duration / steps as f64,
                memory_controls(),
                || false,
            )
            .unwrap();
            near(out.work_heat_residual_j, 0.0, 3e-10);
            for e in 0..3 {
                heat[e] += out.connector_dissipation_j[e];
                work[e] += out.connector_actuator_work_j[e];
                assert!(out.connector_dissipation_j[e] >= 0.0);
            }
            (states, history) = out.into_state_and_memory();
        }
        for e in 0..3 {
            errors[rung] += norm(states[e].pose().position_world().sub(expected.body.x[e]))
                + norm(states[e].linear_momentum_world().sub(expected.body.p[e]))
                + norm(states[e].angular_momentum_body().sub(expected.body.h[e]))
                + (heat[e] - expected.body.heat[e]).abs()
                + (work[e] - expected.work[e]).abs();
            for c in 0..6 {
                errors[rung] += (history[e][0][c] - expected.memory[e][c]).abs();
            }
            for (a, b) in states[e]
                .pose()
                .orientation()
                .components()
                .into_iter()
                .zip(expected.body.q[e])
            {
                errors[rung] += (a - b).abs();
            }
        }
    }
    for pair in errors.windows(2) {
        assert!((3.8..4.2).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

fn compose(a: UnitQuaternion, b: UnitQuaternion) -> UnitQuaternion {
    let [w, x, y, z] = a.components();
    let [p, q, r, s] = b.components();
    UnitQuaternion::new(
        w * p - x * q - y * r - z * s,
        w * q + x * p + y * s - z * r,
        w * r - x * s + y * p + z * q,
        w * s + x * r - y * q + z * p,
    )
    .unwrap()
}

#[test]
fn g1_loaded_material_memory_is_objective_under_a_common_world_rotation() {
    let (states, properties, viscous) = full();
    let links = memory_links(viscous);
    let history = memory();
    let efforts = drives().map(HeldInternalEffort::MaterialFirst);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.5), 0.7).unwrap();
    let rotated: [RigidBodyState; 3] = core::array::from_fn(|i| {
        RigidBodyState::new(
            Pose::new(
                rotation.rotate_body_to_world(states[i].pose().position_world()),
                compose(rotation, states[i].pose().orientation()),
            )
            .unwrap(),
            rotation.rotate_body_to_world(states[i].linear_momentum_world()),
            states[i].angular_momentum_body(),
        )
        .unwrap()
    });
    let wrenches = [
        Wrench {
            force_world: Vec3::new(0.1, 0.2, -0.1),
            torque_body: Vec3::new(0.01, 0.0, 0.02),
        },
        Wrench {
            force_world: Vec3::new(-0.1, 0.05, 0.2),
            torque_body: Vec3::new(0.0, 0.03, -0.01),
        },
        Wrench {
            force_world: Vec3::new(0.02, -0.03, 0.04),
            torque_body: Vec3::new(-0.01, 0.02, 0.01),
        },
    ];
    let gravity = Vec3::new(0.0, 0.0, -0.2);
    let transformed = wrenches.map(|w| Wrench {
        force_world: rotation.rotate_body_to_world(w.force_world),
        torque_body: w.torque_body,
    });
    let mut control = memory_controls();
    control.source.impulse.angular_origin_world_m = Vec3::ZERO;
    let original = step_relaxation_network_with_loads(
        &states,
        &properties,
        &links,
        &history,
        &efforts,
        HeldExternalLoads {
            wrenches: &wrenches,
            gravity: Gravity::new(gravity).unwrap(),
        },
        0.002,
        control,
        || false,
    )
    .unwrap();
    let other = step_relaxation_network_with_loads(
        &rotated,
        &properties,
        &links,
        &history,
        &efforts,
        HeldExternalLoads {
            wrenches: &transformed,
            gravity: Gravity::new(rotation.rotate_body_to_world(gravity)).unwrap(),
        },
        0.002,
        control,
        || false,
    )
    .unwrap();
    for i in 0..3 {
        vector(
            other.states_after()[i].pose().position_world(),
            rotation.rotate_body_to_world(original.states_after()[i].pose().position_world()),
            4e-11,
        );
        vector(
            other.states_after()[i].linear_momentum_world(),
            rotation.rotate_body_to_world(original.states_after()[i].linear_momentum_world()),
            4e-11,
        );
        vector(
            other.states_after()[i].angular_momentum_body(),
            original.states_after()[i].angular_momentum_body(),
            4e-11,
        );
        for c in 0..6 {
            near(
                other.memory_after()[i][0][c],
                original.memory_after()[i][0][c],
                4e-11,
            );
        }
        near(
            other.connector_dissipation_j[i],
            original.connector_dissipation_j[i],
            4e-11,
        );
        near(
            other.connector_actuator_work_j[i],
            original.connector_actuator_work_j[i],
            4e-11,
        );
    }
    near(
        other.memory_energy_change_j,
        original.memory_energy_change_j,
        4e-11,
    );
    near(
        other.external_loads.applied_work_j,
        original.external_loads.applied_work_j,
        4e-11,
    );
    near(
        other.gravitational_potential_energy_change_j,
        original.gravitational_potential_energy_change_j,
        4e-11,
    );
    near(original.work_heat_residual_j, 0.0, 3e-10);
    near(other.work_heat_residual_j, 0.0, 3e-10);
}

#[test]
fn g0_whole_moving_history_transaction_checks_bounds_and_late_cancellation() {
    let (states, properties, viscous) = full();
    let links = memory_links(viscous);
    let memory = memory();
    let preserved = memory.clone();
    let efforts = drives().map(HeldInternalEffort::MaterialFirst);
    let mut count = 0;
    let expected = step_relaxation_network(
        &states,
        &properties,
        &links,
        &memory,
        &efforts,
        0.002,
        memory_controls(),
        || {
            count += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    assert_eq!(
        step_relaxation_network(
            &states,
            &properties,
            &links,
            &memory,
            &efforts,
            0.002,
            memory_controls(),
            || {
                seen += 1;
                seen == count
            }
        ),
        Err(RelaxationNetworkMotionError::Cancelled)
    );
    assert_eq!(memory, preserved);
    assert_eq!(
        step_relaxation_network(
            &states,
            &properties,
            &links,
            &memory,
            &efforts,
            0.002,
            memory_controls(),
            || false
        )
        .unwrap(),
        expected
    );
    for kind in 0..3 {
        let mut limits = memory_controls().limits;
        match kind {
            0 => limits.maximum_branches = 2,
            1 => limits.source.maximum_live_bytes = expected.planned_live_bytes - 1,
            _ => limits.source.maximum_preparation_terms = expected.planned_preparation_terms - 1,
        }
        assert!(
            step_relaxation_network(
                &states,
                &properties,
                &links,
                &memory,
                &efforts,
                0.002,
                RelaxationNetworkMotionControls {
                    limits,
                    ..memory_controls()
                },
                || false
            )
            .is_err()
        );
    }
}
