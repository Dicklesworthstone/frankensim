//! Complete Maxwell memory, exact native graph recoil and independent heat.
use super::*;
use fs_mbd::spatial_loads::internal_impulse::viscous::network::{
    driven::{DrivenViscousNetworkLimits, PreparedDrivenViscousNetwork},
    relaxation::{
        PreparedRelaxationNetwork, RelaxationNetworkError, RelaxationNetworkLimits,
        RelaxationNetworkLink, SpatialRelaxationBranch, relaxation_memory_energy_change,
    },
};

fn memory_limits() -> RelaxationNetworkLimits {
    RelaxationNetworkLimits {
        source: DrivenViscousNetworkLimits {
            maximum_bodies: 8,
            maximum_connectors: 8,
            maximum_live_bytes: 128_000_000,
            maximum_preparation_terms: 80_000_000_000,
        },
        maximum_branches: 8,
    }
}
fn branch(stiffness: f64, tau: f64) -> SpatialRelaxationBranch {
    let mut rows = [[0.0; 6]; 6];
    rows[0][0] = stiffness.sqrt();
    SpatialRelaxationBranch::from_gram_rows(rows, tau).unwrap()
}

#[test]
fn g1_original_single_maxwell_arm_matches_analytic_recoil_storage_and_heat() {
    let states = line_states([-0.5, 0.5, 0.0]);
    let properties = line_properties();
    let links = [RelaxationNetworkLink {
        viscous: axial_link([0, 1], 0.0),
        branches: vec![branch(4.0, 1.0)],
    }];
    let memory = vec![vec![[0.0; 6]]];
    let dt = 0.7_f64;
    let mut source = PreparedRelaxationNetwork::new(
        &states.map(|s| s.pose()),
        &properties,
        &links,
        dt,
        memory_limits(),
        || false,
    )
    .unwrap();
    let out = source
        .advance(&states, &memory, &[[0.0; 6]], controls(), || false)
        .unwrap();
    // Original reduced mass is one: u'=-2a, a'=2u-a.
    let omega = 3.75_f64.sqrt();
    let e = (-0.5 * dt).exp();
    let u = e * ((omega * dt).cos() + 0.5 / omega * (omega * dt).sin());
    let a = e * 2.0 / omega * (omega * dt).sin();
    let b = 2.0 * omega;
    let cosine_integral =
        (1.0 + (-dt).exp() * (-(b * dt).cos() + b * (b * dt).sin())) / (1.0 + b * b);
    let heat = 2.0 / (omega * omega) * (-(-dt).exp_m1() - cosine_integral);
    near(out.states_after[0].linear_momentum_world().x, -u, 3e-12);
    near(out.states_after[1].linear_momentum_world().x, u, 3e-12);
    near(out.memory_after[0][0][0], a, 3e-12);
    near(out.memory_energy_change_j, 0.5 * a * a, 3e-12);
    near(out.dissipated_energy_j, heat, 3e-12);
    near(out.kinetic_memory_balance_residual_j, 0.0, 3e-12);
    near(out.connector_memory_balance_residual_j, 0.0, 3e-12);
    near(out.source_balance_residual_j, 0.0, 3e-12);
    assert_eq!(out.states_after[2], states[2]);
    assert!(out.memory_energy_change_j > 0.0);
    assert!(out.dissipated_energy_j > 0.0);
    assert_eq!(out.actuator_work_j, 0.0);
}

#[test]
fn g1_unobservable_and_equal_time_histories_retain_their_original_heat() {
    let states = line_states([0.0; 3]);
    let properties = line_properties();
    let mut rows = [[0.0; 6]; 6];
    rows[0][0] = 1.0;
    rows[1][0] = 1.0;
    let arm = SpatialRelaxationBranch::from_gram_rows(rows, 0.4).unwrap();
    let links = [
        RelaxationNetworkLink {
            viscous: axial_link([0, 1], 0.0),
            branches: vec![arm.clone(), arm.clone()],
        },
        RelaxationNetworkLink {
            viscous: axial_link([1, 2], 0.0),
            branches: vec![arm],
        },
    ];
    // Opposite stress coordinates cancel mechanically. Their stored energy
    // must still relax, including the tiny third history on its own receiver.
    let memory = vec![
        vec![
            [2.0, -2.0, 0.0, 0.0, 0.0, 0.0],
            [1.0, -1.0, 0.0, 0.0, 0.0, 0.0],
        ],
        vec![[1e-12, -1e-12, 0.0, 0.0, 0.0, 0.0]],
    ];
    let t = 0.3_f64;
    let mut source = PreparedRelaxationNetwork::new(
        &states.map(|s| s.pose()),
        &properties,
        &links,
        t,
        memory_limits(),
        || false,
    )
    .unwrap();
    let out = source
        .advance(&states, &memory, &[[0.0; 6]; 2], controls(), || false)
        .unwrap();
    for (before, after) in memory
        .iter()
        .flatten()
        .zip(out.memory_after.iter().flatten())
    {
        for c in 0..6 {
            near(
                after[c],
                before[c] * (-t / 0.4).exp(),
                3e-12 * before[c].abs(),
            );
        }
    }
    let fraction = -(-2.0 * t / 0.4).exp_m1();
    near(out.connector_dissipation_j[0], 5.0 * fraction, 3e-12);
    assert!(out.connector_dissipation_j[1] > 0.0);
    near(out.connector_dissipation_j[1] / 1e-24, fraction, 3e-12);
    near(out.memory_energy_change_j, -5.0 * fraction, 3e-12);
    for (a, b) in states.iter().zip(&out.states_after) {
        vector(a.linear_momentum_world(), b.linear_momentum_world(), 3e-13);
        vector(a.angular_momentum_body(), b.angular_momentum_body(), 3e-13);
    }
    // The exact null-stress initial history also superposes unchanged on a
    // nonzero drive response; it cannot exert force or change actuator work.
    let efforts = [
        [0.2, 0.0, 0.0, 0.0, 0.0, 0.0],
        [-0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    let unstrained = vec![vec![[0.0; 6]; 2], vec![[0.0; 6]]];
    let visible = source
        .advance(&states, &unstrained, &efforts, controls(), || false)
        .unwrap();
    let driven = source
        .advance(&states, &memory, &efforts, controls(), || false)
        .unwrap();
    assert_eq!(driven.states_after, visible.states_after);
    assert_eq!(
        driven.connector_actuator_work_j,
        visible.connector_actuator_work_j
    );
    for e in 0..2 {
        near(
            driven.connector_dissipation_j[e],
            visible.connector_dissipation_j[e] + out.connector_dissipation_j[e],
            3e-12,
        );
        for branch in 0..memory[e].len() {
            for c in 0..6 {
                near(
                    driven.memory_after[e][branch][c],
                    visible.memory_after[e][branch][c] + memory[e][branch][c] * (-t / 0.4).exp(),
                    3e-12 * memory[e][branch][c].abs().max(1.0),
                );
            }
        }
    }
}

#[test]
fn g3_no_memory_reproduces_the_original_full_tensor_driven_network() {
    let states = full_states();
    let properties = full_properties();
    let viscous = cycle_links();
    let links = viscous.clone().map(|viscous| RelaxationNetworkLink {
        viscous,
        branches: Vec::new(),
    });
    let efforts = [[0.2, -0.1, 0.3, 0.04, -0.02, 0.01]; 3];
    let mut original = PreparedDrivenViscousNetwork::new(
        &states.map(|s| s.pose()),
        &properties,
        &viscous,
        0.2,
        memory_limits().source,
        || false,
    )
    .unwrap();
    let expected = original
        .advance(&states, &efforts, controls(), || false)
        .unwrap();
    let mut source = PreparedRelaxationNetwork::new(
        &states.map(|s| s.pose()),
        &properties,
        &links,
        0.2,
        memory_limits(),
        || false,
    )
    .unwrap();
    let out = source
        .advance(&states, &vec![vec![]; 3], &efforts, controls(), || false)
        .unwrap();
    for i in 0..3 {
        vector(
            out.states_after[i].linear_momentum_world(),
            expected.states_after[i].linear_momentum_world(),
            4e-11,
        );
        vector(
            out.states_after[i].angular_momentum_body(),
            expected.states_after[i].angular_momentum_body(),
            4e-11,
        );
        near(
            out.connector_dissipation_j[i],
            expected.connector_dissipation_j[i],
            4e-11,
        );
        near(
            out.connector_actuator_work_j[i],
            expected.connector_actuator_work_j[i],
            4e-11,
        );
    }
    assert_eq!(out.memory_after, vec![Vec::<[f64; 6]>::new(); 3]);
    assert_eq!(out.memory_energy_change_j, 0.0);
}

fn full_memory_links() -> [RelaxationNetworkLink; 3] {
    let source = cycle_links();
    core::array::from_fn(|e| {
        let diagonal = [1.1, 0.9, 0.7, 0.3, 0.5, 0.4];
        let rows = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                if i == j {
                    diagonal[i]
                } else if j == (i + 1) % 6 {
                    0.15 * (e + 1) as f64
                } else {
                    0.0
                }
            })
        });
        RelaxationNetworkLink {
            viscous: source[e].clone(),
            branches: vec![
                SpatialRelaxationBranch::from_gram_rows(rows, [0.3, 0.6, 1.1][e]).unwrap(),
            ],
        }
    })
}

// Independent physical momentum/torque ODE on all three native bodies.
// No production kinetic map, exact propagator, source observer or energy ledger
// is used. Coordinates are body p/h, original six-row memories, link heat/work.
fn rhs(
    y: [f64; 42],
    poses: &[Pose; 3],
    links: &[RelaxationNetworkLink; 3],
    efforts: &[[f64; 6]; 3],
) -> [f64; 42] {
    let velocity: [Vec3; 3] = core::array::from_fn(|i| {
        v([y[6 * i], y[6 * i + 1], y[6 * i + 2]]).scale(1.0 / [2.0, 3.0, 5.0][i])
    });
    let spin: [Vec3; 3] = core::array::from_fn(|i| {
        poses[i]
            .orientation()
            .rotate_body_to_world(omega(v([y[6 * i + 3], y[6 * i + 4], y[6 * i + 5]]), i))
    });
    let mut out = [0.0; 42];
    for (e, source) in links.iter().enumerate() {
        let link = &source.viscous;
        let [first, second] = link.bodies;
        let qa = poses[first].orientation();
        let rb = poses[second]
            .orientation()
            .rotate_body_to_world(link.arms_body_m[1]);
        let arm = poses[second]
            .position_world()
            .sub(poses[first].position_world())
            .add(rb);
        let linear = qa.rotate_world_to_body(
            velocity[second]
                .sub(velocity[first])
                .add(spin[second].cross(rb))
                .sub(spin[first].cross(arm)),
        );
        let angular = qa.rotate_world_to_body(spin[second].sub(spin[first]));
        let s = [
            linear.x, linear.y, linear.z, angular.x, angular.y, angular.z,
        ];
        let c = link.damping.native_coefficients();
        let viscous: [f64; 6] = core::array::from_fn(|i| (0..6).map(|j| c[i][j] * s[j]).sum());
        let branch = &source.branches[0];
        let l = branch.gram_rows();
        let passive: [f64; 6] = core::array::from_fn(|i| {
            viscous[i] + (0..6).map(|j| l[j][i] * y[18 + 6 * e + j]).sum::<f64>()
        });
        let force =
            qa.rotate_body_to_world(v(core::array::from_fn(|i| efforts[e][i] - passive[i])));
        let couple = qa.rotate_body_to_world(v(core::array::from_fn(|i| {
            efforts[e][i + 3] - passive[i + 3]
        })));
        let moments = [
            arm.cross(force).add(couple).scale(-1.0),
            rb.cross(force).add(couple),
        ];
        for (side, body) in [first, second].into_iter().enumerate() {
            let moment = poses[body]
                .orientation()
                .rotate_world_to_body(moments[side]);
            for i in 0..3 {
                out[6 * body + i] += a(force)[i] * if side == 0 { -1.0 } else { 1.0 };
                out[6 * body + 3 + i] += a(moment)[i];
            }
        }
        for i in 0..6 {
            out[18 + 6 * e + i] = (0..6).map(|j| l[i][j] * s[j]).sum::<f64>()
                - y[18 + 6 * e + i] / branch.relaxation_time_s();
        }
        out[36 + e] = (0..6)
            .map(|i| {
                viscous[i] * s[i]
                    + y[18 + 6 * e + i] * y[18 + 6 * e + i] / branch.relaxation_time_s()
            })
            .sum();
        out[39 + e] = (0..6).map(|i| efforts[e][i] * s[i]).sum();
    }
    out
}

#[test]
fn g1_full_coupled_memory_cycle_matches_independent_ode_and_continuation() {
    let states = full_states();
    let properties = full_properties();
    let poses = states.map(|s| s.pose());
    let links = full_memory_links();
    let memory: Vec<_> = (0..3)
        .map(|e| {
            vec![core::array::from_fn(|c| {
                0.03 * (e + 1) as f64 * (c + 1) as f64 * if c % 2 == 0 { 1.0 } else { -1.0 }
            })]
        })
        .collect();
    let efforts = [
        [0.2, -0.1, 0.3, 0.04, -0.03, 0.02],
        [-0.3, 0.2, 0.1, -0.02, 0.05, 0.01],
        [0.1, 0.15, -0.2, 0.02, 0.01, -0.04],
    ];
    let mut y = [0.0; 42];
    y[..18].copy_from_slice(&initial_ode(&states)[..18]);
    for e in 0..3 {
        y[18 + 6 * e..18 + 6 * (e + 1)].copy_from_slice(&memory[e][0]);
    }
    let dt = 0.2 / 1024.0;
    for _ in 0..1024 {
        let k1 = rhs(y, &poses, &links, &efforts);
        let k2 = rhs(
            core::array::from_fn(|i| y[i] + 0.5 * dt * k1[i]),
            &poses,
            &links,
            &efforts,
        );
        let k3 = rhs(
            core::array::from_fn(|i| y[i] + 0.5 * dt * k2[i]),
            &poses,
            &links,
            &efforts,
        );
        let k4 = rhs(
            core::array::from_fn(|i| y[i] + dt * k3[i]),
            &poses,
            &links,
            &efforts,
        );
        for i in 0..42 {
            y[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }
    let mut source =
        PreparedRelaxationNetwork::new(&poses, &properties, &links, 0.2, memory_limits(), || false)
            .unwrap();
    let out = source
        .advance(&states, &memory, &efforts, controls(), || false)
        .unwrap();
    for e in 0..3 {
        vector(
            out.states_after[e].linear_momentum_world(),
            v(core::array::from_fn(|j| y[6 * e + j])),
            3e-10,
        );
        vector(
            out.states_after[e].angular_momentum_body(),
            v(core::array::from_fn(|j| y[6 * e + 3 + j])),
            3e-10,
        );
        for c in 0..6 {
            near(out.memory_after[e][0][c], y[18 + 6 * e + c], 3e-10);
        }
        near(out.connector_dissipation_j[e], y[36 + e], 3e-10);
        near(out.connector_actuator_work_j[e], y[39 + e], 3e-10);
        assert!(out.connector_dissipation_j[e] > 0.0);
    }
    let mut half =
        PreparedRelaxationNetwork::new(&poses, &properties, &links, 0.1, memory_limits(), || false)
            .unwrap();
    let first = half
        .advance(&states, &memory, &efforts, controls(), || false)
        .unwrap();
    let second = half
        .advance(
            &first.states_after,
            &first.memory_after,
            &efforts,
            controls(),
            || false,
        )
        .unwrap();
    for e in 0..3 {
        vector(
            second.states_after[e].linear_momentum_world(),
            out.states_after[e].linear_momentum_world(),
            3e-11,
        );
        vector(
            second.states_after[e].angular_momentum_body(),
            out.states_after[e].angular_momentum_body(),
            3e-11,
        );
        for c in 0..6 {
            near(
                second.memory_after[e][0][c],
                out.memory_after[e][0][c],
                3e-11,
            );
        }
        near(
            first.connector_dissipation_j[e] + second.connector_dissipation_j[e],
            out.connector_dissipation_j[e],
            3e-11,
        );
        near(
            first.connector_actuator_work_j[e] + second.connector_actuator_work_j[e],
            out.connector_actuator_work_j[e],
            3e-11,
        );
    }
}

#[test]
fn g0_memory_refusals_and_complete_bounds_are_atomic_with_deterministic_retry() {
    let states = line_states([-0.5, 0.5, 0.0]);
    let properties = line_properties();
    let links = [RelaxationNetworkLink {
        viscous: axial_link([0, 1], 0.3),
        branches: vec![branch(4.0, 0.6)],
    }];
    let memory = vec![vec![[0.1; 6]]];
    let original = memory.clone();
    let efforts = [[0.2, 0.0, 0.0, 0.0, 0.0, 0.0]];
    let poses = states.map(|s| s.pose());
    let mut source =
        PreparedRelaxationNetwork::new(&poses, &properties, &links, 0.2, memory_limits(), || false)
            .unwrap();
    let mut count = 0;
    let expected = source
        .advance(&states, &memory, &efforts, controls(), || {
            count += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert_eq!(
        source.advance(&states, &memory, &efforts, controls(), || {
            seen += 1;
            seen == count
        }),
        Err(RelaxationNetworkError::Cancelled)
    );
    assert_eq!(memory, original);
    assert_eq!(
        source
            .advance(&states, &memory, &efforts, controls(), || false)
            .unwrap(),
        expected
    );
    for bad in [vec![], vec![vec![]], vec![vec![[f64::NAN; 6]]]] {
        assert!(
            source
                .advance(&states, &bad, &efforts, controls(), || false)
                .is_err()
        );
    }
    let mut limits = memory_limits();
    limits.maximum_branches = 0;
    assert!(
        PreparedRelaxationNetwork::new(&poses, &properties, &links, 0.2, limits, || false).is_err()
    );
    limits = memory_limits();
    limits.source.maximum_live_bytes = source.planned_live_bytes() - 1;
    assert!(
        PreparedRelaxationNetwork::new(&poses, &properties, &links, 0.2, limits, || false).is_err()
    );
    limits = memory_limits();
    limits.source.maximum_preparation_terms = source.planned_preparation_terms() - 1;
    assert!(
        PreparedRelaxationNetwork::new(&poses, &properties, &links, 0.2, limits, || false).is_err()
    );
    let mut zero =
        PreparedRelaxationNetwork::new(&poses, &properties, &links, 0.0, memory_limits(), || false)
            .unwrap();
    let out = zero
        .advance(&states, &memory, &[[f64::MAX; 6]], controls(), || false)
        .unwrap();
    assert_eq!(out.states_after, states);
    assert_eq!(out.memory_after, memory);
    assert_eq!(out.dissipated_energy_j, 0.0);
    assert_eq!(out.memory_energy_change_j, 0.0);
    for tau in [0.0, -1.0, f64::INFINITY, f64::NAN, f64::from_bits(1)] {
        assert!(SpatialRelaxationBranch::from_gram_rows([[0.0; 6]; 6], tau).is_err());
    }
}

#[test]
fn g1_original_memory_change_retains_sub_ulp_absolute_storage_release() {
    let before: Vec<Vec<[f64; 6]>> = vec![vec![[100_000_000.0, 0.0, 0.0, 0.0, 0.0, 0.0]]];
    let mut after = before.clone();
    after[0][0][0] = f64::from_bits(before[0][0][0].to_bits() - 1);
    let a = before[0][0][0];
    let b = after[0][0][0];
    let expected = (b - a) * (0.5 * b + 0.5 * a);
    near(
        relaxation_memory_energy_change(&before, &after).unwrap(),
        expected,
        2e-16,
    );
    assert_ne!(expected, 0.5 * b * b - 0.5 * a * a);
}
