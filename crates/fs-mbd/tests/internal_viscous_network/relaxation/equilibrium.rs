//! Exact equilibrium storage beside original Maxwell and Kelvin dissipation.
use super::*;
use fs_mbd::spatial_loads::internal_impulse::viscous::network::relaxation::{
    HeldRelaxationKinetics, SpatialMemoryEvolution,
};

fn equilibrium(stiffness: f64) -> SpatialRelaxationBranch {
    let mut rows = [[0.0; 6]; 6];
    rows[0][0] = stiffness.sqrt();
    SpatialRelaxationBranch::from_equilibrium_gram_rows(rows).unwrap()
}

#[test]
fn equilibrium_channel_recovers_exact_oscillation_and_signed_drive_with_zero_heat() {
    let states = line_states([-0.5, 0.5, 0.0]);
    let properties = line_properties();
    let links = [RelaxationNetworkLink {
        viscous: axial_link([0, 1], 0.0),
        branches: vec![equilibrium(4.0)],
    }];
    assert_eq!(
        links[0].branches[0].evolution(),
        SpatialMemoryEvolution::Equilibrium
    );
    assert_eq!(links[0].branches[0].relaxation_time_s(), None);
    assert_eq!(links[0].branches[0].relaxation_rate_per_s(), 0.0);
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
    assert_eq!(
        source.effective_evolution(0, 0),
        Some(SpatialMemoryEvolution::Equilibrium)
    );
    assert_eq!(source.effective_relaxation_time_s(0, 0), None);
    assert_eq!(source.effective_evolution(0, 1), None);
    for force in [0.0, 0.3, -0.2] {
        let mut effort = [0.0; 6];
        effort[0] = force;
        let out = source
            .advance(&states, &[vec![[0.0; 6]]], &[effort], controls(), || false)
            .unwrap();
        let u = (2.0 * dt).cos() + 0.5 * force * (2.0 * dt).sin();
        let amplitude = (2.0 * dt).sin() + 0.5 * force * (1.0 - (2.0 * dt).cos());
        near(out.states_after[0].linear_momentum_world().x, -u, 3e-12);
        near(out.states_after[1].linear_momentum_world().x, u, 3e-12);
        near(out.memory_after[0][0][0], amplitude, 3e-12);
        near(
            out.average_relative_twist[0][0] * dt,
            0.5 * amplitude,
            3e-12,
        );
        near(
            out.memory_energy_change_j,
            0.5 * amplitude * amplitude,
            3e-12,
        );
        near(out.actuator_work_j, 0.5 * force * amplitude, 3e-12);
        assert_eq!(out.dissipated_energy_j, 0.0);
        assert_eq!(out.connector_dissipation_j, [0.0]);
        assert_eq!(out.states_after[2], states[2]);
    }
}

#[test]
fn complete_equilibrium_maxwell_kelvin_law_matches_independent_forced_ode_and_continuation() {
    let states = line_states([-0.5, 0.5, 0.0]);
    let properties = line_properties();
    let links = [RelaxationNetworkLink {
        viscous: axial_link([0, 1], 0.3),
        branches: vec![equilibrium(4.0), branch(2.25, 0.4)],
    }];
    let memory = vec![vec![
        [0.2, 0.0, 0.0, 0.0, 0.0, 0.0],
        [-0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
    ]];
    let scales = [vec![7.0, 0.5]];
    let kinetics = HeldRelaxationKinetics {
        viscosity_scales: &[1.6],
        relaxation_time_scales: &scales,
    };
    kinetics.validate(&links).unwrap();
    let rhs = |y: [f64; 6]| -> [f64; 6] {
        let [u, eq, arm, _, _, _] = y;
        [
            0.7 - 0.48 * u - 2.0 * eq - 1.5 * arm,
            2.0 * u,
            1.5 * u - arm / 0.2,
            0.48 * u * u + arm * arm / 0.2,
            0.7 * u,
            u,
        ]
    };
    let mut y = [1.0, 0.2, -0.1, 0.0, 0.0, 0.0];
    let h = 0.6 / 2048.0;
    for _ in 0..2048 {
        let k1 = rhs(y);
        let k2 = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * k1[i]));
        let k3 = rhs(core::array::from_fn(|i| y[i] + 0.5 * h * k2[i]));
        let k4 = rhs(core::array::from_fn(|i| y[i] + h * k3[i]));
        y = core::array::from_fn(|i| y[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    let mut source = PreparedRelaxationNetwork::new_with_kinetics(
        &states.map(|s| s.pose()),
        &properties,
        &links,
        kinetics,
        0.6,
        memory_limits(),
        || false,
    )
    .unwrap();
    assert_eq!(
        source.effective_evolution(0, 0),
        Some(SpatialMemoryEvolution::Equilibrium)
    );
    assert_eq!(source.effective_relaxation_time_s(0, 1), Some(0.2));
    let efforts = [[0.7, 0.0, 0.0, 0.0, 0.0, 0.0]];
    let out = source
        .advance(&states, &memory, &efforts, controls(), || false)
        .unwrap();
    near(out.states_after[1].linear_momentum_world().x, y[0], 4e-11);
    near(out.memory_after[0][0][0], y[1], 4e-11);
    near(out.memory_after[0][1][0], y[2], 4e-11);
    near(out.dissipated_energy_j, y[3], 4e-11);
    near(out.actuator_work_j, y[4], 4e-11);
    near(out.average_relative_twist[0][0] * 0.6, y[5], 4e-11);
    near(out.kinetic_memory_balance_residual_j, 0.0, 4e-11);
    assert!(out.dissipated_energy_j > 0.0);
    let mut half = PreparedRelaxationNetwork::new_with_kinetics(
        &states.map(|s| s.pose()),
        &properties,
        &links,
        kinetics,
        0.3,
        memory_limits(),
        || false,
    )
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
    near(
        second.states_after[1].linear_momentum_world().x,
        y[0],
        4e-11,
    );
    for b in 0..2 {
        near(second.memory_after[0][b][0], y[b + 1], 4e-11);
    }
    near(
        first.dissipated_energy_j + second.dissipated_energy_j,
        y[3],
        4e-11,
    );
}

#[test]
fn lossless_null_history_remains_complete_beside_tiny_relaxing_heat_and_shifts() {
    let states = line_states([0.0; 3]);
    let properties = line_properties();
    let mut rows = [[0.0; 6]; 6];
    rows[0][0] = 1.0;
    rows[1][0] = 1.0;
    let links = [
        RelaxationNetworkLink {
            viscous: axial_link([0, 1], 0.0),
            branches: vec![SpatialRelaxationBranch::from_equilibrium_gram_rows(rows).unwrap()],
        },
        RelaxationNetworkLink {
            viscous: axial_link([1, 2], 0.0),
            branches: vec![SpatialRelaxationBranch::from_gram_rows(rows, 0.4).unwrap()],
        },
    ];
    let memory = vec![
        vec![[2.0, -2.0, 0.0, 0.0, 0.0, 0.0]],
        vec![[1e-12, -1e-12, 0.0, 0.0, 0.0, 0.0]],
    ];
    let scales = [vec![1e200], vec![1.0]];
    let kinetics = HeldRelaxationKinetics {
        viscosity_scales: &[1.0, 1.0],
        relaxation_time_scales: &scales,
    };
    let mut source = PreparedRelaxationNetwork::new_with_kinetics(
        &states.map(|s| s.pose()),
        &properties,
        &links,
        kinetics,
        0.3,
        memory_limits(),
        || false,
    )
    .unwrap();
    let out = source
        .advance(&states, &memory, &[[0.0; 6]; 2], controls(), || false)
        .unwrap();
    assert_eq!(out.memory_after[0], memory[0]);
    near(out.memory_after[1][0][0], 1e-12 * (-0.75_f64).exp(), 3e-24);
    assert_eq!(out.connector_dissipation_j[0], 0.0);
    near(
        out.connector_dissipation_j[1] / 1e-24,
        -(-1.5_f64).exp_m1(),
        3e-12,
    );
    let invalid = [vec![0.0], vec![1.0]];
    assert!(
        HeldRelaxationKinetics {
            relaxation_time_scales: &invalid,
            ..kinetics
        }
        .validate(&links)
        .is_err()
    );
    assert!(SpatialRelaxationBranch::from_gram_rows(rows, f64::INFINITY).is_err());
    assert!(SpatialRelaxationBranch::from_gram_rows(rows, 0.0).is_err());
    let cancelled = source.advance(&states, &memory, &[[0.0; 6]; 2], controls(), || true);
    assert_eq!(cancelled.unwrap_err(), RelaxationNetworkError::Cancelled);
    let again = source
        .advance(&states, &memory, &[[0.0; 6]; 2], controls(), || false)
        .unwrap();
    assert_eq!(again, out);
}

#[test]
fn projected_kelvin_factor_retains_original_rows_across_scaling_and_dependent_modes() {
    let mut rows = [[0.0; 6]; 6];
    rows[0] = [1.0, 2.0, 0.0, -4.0, 0.0, 0.0];
    rows[1] = [0.3, 0.6, 0.0, -1.2, 0.0, 0.0];
    rows[4] = [0.0, 0.0, 1e-12, 0.0, 0.0, 0.0];
    let scales = BushingCoordinateScales {
        translation_m: 0.002,
        rotation_rad: 0.3,
    };
    let form = PassiveBushingMatrix::from_gram_rows(rows, scales).unwrap();
    assert_eq!(form.coordinate_scales(), scales);
    assert_eq!(
        SpatialRelaxationBranch::from_stiffness(&form, 0.5)
            .unwrap()
            .gram_rows(),
        &rows
    );
    for factor in [0.3_f64, 7.3, 1e-120, 1e120] {
        let shifted = form.scaled(factor).unwrap();
        let actual = SpatialRelaxationBranch::from_stiffness(&shifted, 0.5).unwrap();
        for r in 0..6 {
            for c in 0..6 {
                assert_eq!(actual.gram_rows()[r][c], rows[r][c] * factor.sqrt());
            }
        }
    }
    assert_eq!(
        SpatialRelaxationBranch::from_stiffness(&form, 0.5)
            .unwrap()
            .elastic_memory([2.0, -1.0, 0.0, 0.0, 0.0, 0.0])
            .unwrap(),
        [0.0; 6]
    );
    let mut invalid = rows;
    invalid[0][0] = f64::INFINITY;
    assert!(PassiveBushingMatrix::from_gram_rows(invalid, scales).is_err());
}
