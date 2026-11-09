//! Held temperature/kinetic coefficients preserve the original graph storage.
use super::*;
use fs_mbd::spatial_loads::internal_impulse::viscous::network::relaxation::HeldRelaxationKinetics;

fn held<'a>(viscosity: &'a [f64], times: &'a [Vec<f64>]) -> HeldRelaxationKinetics<'a> {
    HeldRelaxationKinetics {
        viscosity_scales: viscosity,
        relaxation_time_scales: times,
    }
}

#[test]
fn g3_identity_kinetics_is_the_original_full_tensor_source_bit_for_bit() {
    let states = full_states();
    let properties = full_properties();
    let links = full_memory_links();
    let memory: Vec<_> = (0..3).map(|e| vec![[0.02 * (e + 1) as f64; 6]]).collect();
    let efforts = [[0.2, -0.1, 0.3, 0.04, -0.03, 0.02]; 3];
    let poses = states.map(|s| s.pose());
    let times = vec![vec![1.0]; 3];
    let kinetics = held(&[1.0; 3], &times);
    kinetics.validate(&links).unwrap();
    let mut reference =
        PreparedRelaxationNetwork::new(&poses, &properties, &links, 0.2, memory_limits(), || false)
            .unwrap();
    let mut source = PreparedRelaxationNetwork::new_with_kinetics(
        &poses,
        &properties,
        &links,
        kinetics,
        0.2,
        memory_limits(),
        || false,
    )
    .unwrap();
    assert_eq!(
        source
            .advance(&states, &memory, &efforts, controls(), || false)
            .unwrap(),
        reference
            .advance(&states, &memory, &efforts, controls(), || false)
            .unwrap()
    );
    for e in 0..3 {
        assert_eq!(source.effective_damping(e), Some(&links[e].viscous.damping));
        assert_eq!(
            source.effective_relaxation_time_s(e, 0),
            Some(links[e].branches[0].relaxation_time_s())
        );
    }
    assert_eq!(source.effective_damping(3), None);
    assert_eq!(source.effective_relaxation_time_s(3, 0), None);
    assert_eq!(source.effective_relaxation_time_s(0, 1), None);
}

#[test]
fn g1_nonuniform_kelvin_and_maxwell_shifts_match_the_original_physical_ode() {
    let states = full_states();
    let properties = full_properties();
    let poses = states.map(|s| s.pose());
    let links = full_memory_links();
    let original = links.clone();
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
    let viscosity = [0.2, 1.7, 0.0];
    let times = vec![vec![0.35], vec![2.0], vec![5.0]];
    let kinetics = held(&viscosity, &times);
    kinetics.validate(&links).unwrap();
    // Reference re-declares the explicit physical coefficients only for the
    // independent momentum/torque/history ODE, never a production prepared map.
    let shifted = core::array::from_fn(|e| {
        let mut link = links[e].clone();
        link.viscous.damping = link.viscous.damping.scaled(viscosity[e]).unwrap();
        link.branches[0] = SpatialRelaxationBranch::from_gram_rows(
            *links[e].branches[0].gram_rows(),
            links[e].branches[0].relaxation_time_s() * times[e][0],
        )
        .unwrap();
        link
    });
    let mut y = [0.0; 42];
    y[..18].copy_from_slice(&initial_ode(&states)[..18]);
    for e in 0..3 {
        y[18 + 6 * e..18 + 6 * (e + 1)].copy_from_slice(&memory[e][0]);
    }
    let h = 0.2 / 1024.0;
    for _ in 0..1024 {
        let k1 = rhs(y, &poses, &shifted, &efforts);
        let k2 = rhs(
            core::array::from_fn(|i| y[i] + 0.5 * h * k1[i]),
            &poses,
            &shifted,
            &efforts,
        );
        let k3 = rhs(
            core::array::from_fn(|i| y[i] + 0.5 * h * k2[i]),
            &poses,
            &shifted,
            &efforts,
        );
        let k4 = rhs(
            core::array::from_fn(|i| y[i] + h * k3[i]),
            &poses,
            &shifted,
            &efforts,
        );
        for i in 0..42 {
            y[i] += h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }
    let mut source = PreparedRelaxationNetwork::new_with_kinetics(
        &poses,
        &properties,
        &links,
        kinetics,
        0.2,
        memory_limits(),
        || false,
    )
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
        assert_eq!(
            source.effective_damping(e).unwrap(),
            &shifted[e].viscous.damping
        );
        assert_eq!(
            source.effective_relaxation_time_s(e, 0),
            Some(shifted[e].branches[0].relaxation_time_s())
        );
    }
    assert_eq!(links, original);
    near(out.kinetic_memory_balance_residual_j, 0.0, 3e-10);
}

#[test]
fn g1_rate_changes_continue_each_original_history_without_storage_reset_or_double_shift() {
    let states = line_states([0.0; 3]);
    let properties = line_properties();
    let mut rows = [[0.0; 6]; 6];
    rows[0][0] = 1.0;
    rows[1][0] = 1.0;
    let arm = |tau| SpatialRelaxationBranch::from_gram_rows(rows, tau).unwrap();
    let links = [
        RelaxationNetworkLink {
            viscous: axial_link([0, 1], 0.3),
            branches: vec![arm(0.2), arm(0.8)],
        },
        RelaxationNetworkLink {
            viscous: axial_link([1, 2], 0.1),
            branches: vec![arm(2.0)],
        },
    ];
    let memory = vec![
        vec![
            [2.0, -2.0, 0.0, 0.0, 0.0, 0.0],
            [0.1, -0.1, 0.0, 0.0, 0.0, 0.0],
        ],
        vec![[1e-12, -1e-12, 0.0, 0.0, 0.0, 0.0]],
    ];
    let poses = states.map(|s| s.pose());
    let first_times = vec![vec![0.5, 2.0], vec![4.0]];
    let next_times = vec![vec![2.0, 0.25], vec![0.8]];
    let first_kinetics = held(&[0.5, 2.0], &first_times);
    let next_kinetics = held(&[3.0, 0.0], &next_times);
    let mut first = PreparedRelaxationNetwork::new_with_kinetics(
        &poses,
        &properties,
        &links,
        first_kinetics,
        0.1,
        memory_limits(),
        || false,
    )
    .unwrap();
    let initial = first
        .advance(&states, &memory, &[[0.0; 6]; 2], controls(), || false)
        .unwrap();
    let mut zero = PreparedRelaxationNetwork::new_with_kinetics(
        &poses,
        &properties,
        &links,
        next_kinetics,
        0.0,
        memory_limits(),
        || false,
    )
    .unwrap();
    let rebound = zero
        .advance(
            &initial.states_after,
            &initial.memory_after,
            &[[0.0; 6]; 2],
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(rebound.states_after, initial.states_after);
    assert_eq!(rebound.memory_after, initial.memory_after);
    assert_eq!(rebound.memory_energy_change_j, 0.0);
    assert_eq!(rebound.dissipated_energy_j, 0.0);
    assert_eq!(rebound.actuator_work_j, 0.0);
    let mut next = PreparedRelaxationNetwork::new_with_kinetics(
        &poses,
        &properties,
        &links,
        next_kinetics,
        0.3,
        memory_limits(),
        || false,
    )
    .unwrap();
    let out = next
        .advance(
            &initial.states_after,
            &initial.memory_after,
            &[[0.0; 6]; 2],
            controls(),
            || false,
        )
        .unwrap();
    for e in 0..2 {
        let mut heat = 0.0;
        for b in 0..memory[e].len() {
            let original_tau = links[e].branches[b].relaxation_time_s();
            let clock =
                0.1 / (original_tau * first_times[e][b]) + 0.3 / (original_tau * next_times[e][b]);
            assert_eq!(
                next.effective_relaxation_time_s(e, b),
                Some(original_tau * next_times[e][b])
            );
            for c in 0..6 {
                let q = memory[e][b][c];
                near(
                    out.memory_after[e][b][c],
                    q * (-clock).exp(),
                    3e-12 * q.abs(),
                );
                heat += 0.5 * q * q * -(-2.0 * clock).exp_m1();
            }
        }
        let actual = initial.connector_dissipation_j[e] + out.connector_dissipation_j[e];
        assert!(actual > 0.0);
        near(actual / heat, 1.0, 3e-12);
    }
    for i in 0..3 {
        assert_eq!(out.states_after[i], states[i]);
    }
}

#[test]
fn g0_effective_coefficient_admission_matches_preparation_and_preserves_inputs() {
    let states = line_states([-0.5, 0.5, 0.0]);
    let properties = line_properties();
    let links = [RelaxationNetworkLink {
        viscous: axial_link([0, 1], 4.0),
        branches: vec![branch(4.0, 2.0)],
    }];
    let original = links.clone();
    let poses = states.map(|s| s.pose());
    let invalid = [
        (vec![], vec![vec![1.0]]),
        (vec![1.0], vec![]),
        (vec![1.0], vec![vec![]]),
        (vec![1.0], vec![vec![1.0, 1.0]]),
        (vec![-1.0], vec![vec![1.0]]),
        (vec![f64::NAN], vec![vec![1.0]]),
        (vec![f64::MAX], vec![vec![1.0]]),
        (vec![1.0], vec![vec![0.0]]),
        (vec![1.0], vec![vec![-1.0]]),
        (vec![1.0], vec![vec![f64::NAN]]),
        (vec![1.0], vec![vec![f64::INFINITY]]),
        (vec![1.0], vec![vec![f64::MAX]]),
        (vec![1.0], vec![vec![f64::from_bits(1)]]),
    ];
    for (viscosity, times) in invalid {
        let kinetics = held(&viscosity, &times);
        assert!(kinetics.validate(&links).is_err());
        assert!(
            PreparedRelaxationNetwork::new_with_kinetics(
                &poses,
                &properties,
                &links,
                kinetics,
                0.1,
                memory_limits(),
                || false
            )
            .is_err()
        );
    }
    let times = vec![vec![0.5]];
    let kinetics = held(&[0.0], &times);
    kinetics.validate(&links).unwrap();
    assert!(matches!(
        PreparedRelaxationNetwork::new_with_kinetics(
            &poses,
            &properties,
            &links,
            kinetics,
            0.1,
            memory_limits(),
            || true
        ),
        Err(RelaxationNetworkError::Cancelled)
    ));
    let mut source = PreparedRelaxationNetwork::new_with_kinetics(
        &poses,
        &properties,
        &links,
        kinetics,
        0.1,
        memory_limits(),
        || false,
    )
    .unwrap();
    let memory = vec![vec![[0.1; 6]]];
    let efforts = [[0.0; 6]];
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
    assert_eq!(
        source
            .advance(&states, &memory, &efforts, controls(), || false)
            .unwrap(),
        expected
    );
    assert_eq!(links, original);
}
