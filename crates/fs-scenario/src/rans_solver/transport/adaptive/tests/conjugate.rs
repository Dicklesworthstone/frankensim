use super::*;
use crate::rans_solver::{
    RansConjugateHeatLoad, RansConjugateWall, RansConjugateWorkspace, RansHeatCoupling,
    RansWallLayer,
};

fn wall() -> RansConjugateWall {
    RansConjugateWall::try_new(
        &[
            RansWallLayer {
                thickness_m: 0.15,
                volumetric_heat_capacity_j_m3_k: 4.0,
                conductivity_w_m_k: 0.1,
            },
            RansWallLayer {
                thickness_m: 0.2,
                volumetric_heat_capacity_j_m3_k: 8.0,
                conductivity_w_m_k: 0.3,
            },
        ],
        &[9.0, 3.0],
        2,
    )
    .unwrap()
}

fn thermal_values(wall: &RansConjugateWall, state: &RansFieldState) -> Vec<f64> {
    wall.temperatures_k()
        .iter()
        .chain(&state.temp_k)
        .copied()
        .collect()
}

// Independent physical-unit semi-discrete generator. Evaluate its exponential
// by a convergent Taylor action, not the production implicit/eliminated solve.
fn exact_thermal(
    s: &RansSolver,
    g: &RansChannelGrid,
    w: &RansConjugateWall,
    state: &RansFieldState,
    time: f64,
) -> Vec<f64> {
    let m = w.layers().len();
    let n = state.temp_k.len();
    let mut capacity: Vec<_> = w
        .layers()
        .iter()
        .map(|p| p.thickness_m * p.volumetric_heat_capacity_j_m3_k)
        .collect();
    let mut fluid_volume = vec![0.0; n];
    for i in 0..n - 1 {
        let half = 0.5 * (g.y_nodes[i + 1] - g.y_nodes[i]);
        fluid_volume[i] += half;
        fluid_volume[i + 1] += half;
    }
    capacity.extend(fluid_volume.iter().map(|v| s.rho * s.cp * v));
    let mut links = Vec::new();
    for (i, p) in w.layers().iter().enumerate() {
        let resistance = p.thickness_m / (2.0 * p.conductivity_w_m_k)
            + w.layers()
                .get(i + 1)
                .map_or(0.0, |q| q.thickness_m / (2.0 * q.conductivity_w_m_k));
        links.push(1.0 / resistance);
    }
    links.extend(g.y_nodes.windows(2).map(|p| s.k_thermal / (p[1] - p[0])));
    let mut a = vec![vec![0.0; m + n]; m + n];
    for (i, &link) in links.iter().enumerate() {
        a[i][i] -= link / capacity[i];
        a[i][i + 1] += link / capacity[i];
        a[i + 1][i] += link / capacity[i + 1];
        a[i + 1][i + 1] -= link / capacity[i + 1];
    }
    let rate = (0..m + n).map(|i| -a[i][i]).fold(0.0_f64, f64::max);
    let pieces = (2.0 * rate * time).ceil().max(1.0) as usize;
    let h = time / pieces as f64;
    let mut values = thermal_values(w, state);
    for _ in 0..pieces {
        let mut term = values.clone();
        let mut sum = values.clone();
        for order in 1..=32 {
            term = a
                .iter()
                .map(|row| {
                    row.iter().zip(&term).map(|(a, x)| a * x).sum::<f64>() * (h / order as f64)
                })
                .collect();
            for (s, t) in sum.iter_mut().zip(&term) {
                *s += t;
            }
        }
        values = sum;
    }
    values
}

#[test]
fn g1_joint_thermal_decay_converges_to_an_independent_matrix_exponential() {
    let g = uniform(8);
    let s = solver();
    let mut initial = laminar(&g);
    initial.u.fill(0.0);
    let initial_wall = wall();
    let duration = 0.5;
    let exact = exact_thermal(&s, &g, &initial_wall, &initial, duration);
    let mut errors = Vec::new();
    let mut steps = Vec::new();
    for tolerance in [1e-3, 1e-6] {
        let mut state = initial.clone();
        let mut wall = initial_wall.clone();
        let mut work = RansConjugateAdaptiveWorkspace::try_new(&g, &wall, 9, 2).unwrap();
        let report = s
            .advance_launder_sharma_conjugate_adaptive(
                &g,
                &mut state,
                &mut wall,
                0.0,
                duration,
                RansConjugateHeatLoad::default(),
                controls(duration, tolerance),
                &mut work,
                || false,
            )
            .unwrap();
        let error = thermal_values(&wall, &state)
            .iter()
            .zip(&exact)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        errors.push(error);
        steps.push(report.fluid.accepted_intervals);
        assert!(report.fluid.rejected_intervals > 0);
        assert!(report.fluid.maximum_accepted_error_ratio <= 1.0);
        assert!(
            state
                .temp_k
                .iter()
                .chain(wall.temperatures_k())
                .all(|&t| (3.0..=9.0).contains(&t))
        );
        let heat = report.fluid.heat.unwrap();
        close(
            heat.energy_change_j_m2 + report.wall.energy_change_j_m2,
            0.0,
            1e-10,
        );
        close(
            heat.wall_heat_j_m2 - report.wall.interface_heat_j_m2,
            -report.wall.interface_equation_defect_j_m2,
            1e-12,
        );
        assert_eq!(work.numerical_bytes(), (24 * 9 + 10 * 2) * 8);
    }
    assert!(
        errors[1] < 0.1 * errors[0] && errors[1] < 2e-3,
        "{errors:?}"
    );
    assert!(steps[1] > steps[0]);
}

#[test]
fn g1_wall_time_scale_controls_trials_even_when_fluid_change_is_small() {
    let g = uniform(4);
    let mut s = solver();
    s.cp = 4000.0;
    let mut initial = laminar(&g);
    initial.u.fill(0.0);
    initial.temp_k.fill(4.0);
    let initial_wall = wall();
    let h = 0.25;
    let mut coarse = initial.clone();
    let mut fine = initial.clone();
    let mut coarse_wall = initial_wall.clone();
    let mut fine_wall = initial_wall.clone();
    let mut kernel = RansConjugateWorkspace::try_new(&g, &initial_wall, 5, 2).unwrap();
    s.advance_launder_sharma_with_conjugate_wall(
        &g,
        &mut coarse,
        &mut coarse_wall,
        0.0,
        h,
        RansConjugateHeatLoad::default(),
        &mut kernel,
        || false,
    )
    .unwrap();
    for _ in 0..2 {
        s.advance_launder_sharma_with_conjugate_wall(
            &g,
            &mut fine,
            &mut fine_wall,
            0.0,
            h / 2.0,
            RansConjugateHeatLoad::default(),
            &mut kernel,
            || false,
        )
        .unwrap();
    }
    let discrepancy = |a: &[f64], b: &[f64]| {
        a.iter()
            .zip(b)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max)
    };
    let fluid_error = discrepancy(&coarse.temp_k, &fine.temp_k);
    let wall_error = discrepancy(coarse_wall.temperatures_k(), fine_wall.temperatures_k());
    assert!(wall_error > 100.0 * fluid_error);
    let tolerance = (fluid_error * wall_error).sqrt();
    let mut state = initial;
    let mut wall = initial_wall;
    let mut work = RansConjugateAdaptiveWorkspace::try_new(&g, &wall, 5, 2).unwrap();
    let report = s
        .advance_launder_sharma_conjugate_adaptive(
            &g,
            &mut state,
            &mut wall,
            0.0,
            h,
            RansConjugateHeatLoad::default(),
            controls(h, tolerance),
            &mut work,
            || false,
        )
        .unwrap();
    assert!(report.fluid.rejected_intervals > 0 && report.fluid.accepted_intervals > 1);
}

#[test]
fn g0_accepted_half_steps_replay_both_states_and_every_heat_input() {
    let g = uniform(8);
    let s = solver();
    let initial = active(&g);
    let initial_wall = wall();
    let h = 1.0 / 64.0;
    let load = RansConjugateHeatLoad {
        exterior_heat_flux_w_m2: 0.7,
        wall_volumetric_heating_w_m3: Some(&[0.2, -0.1]),
        fluid_volumetric_heating_w_m3: Some(&[0.125; 9]),
        coupling: RansHeatCoupling::Dissipative,
        ..RansConjugateHeatLoad::default()
    };
    let mut expected = initial.clone();
    let mut expected_wall = initial_wall.clone();
    let mut kernel = RansConjugateWorkspace::try_new(&g, &initial_wall, 9, 2).unwrap();
    let first = s
        .advance_launder_sharma_with_conjugate_wall(
            &g,
            &mut expected,
            &mut expected_wall,
            -0.125,
            h / 2.0,
            load,
            &mut kernel,
            || false,
        )
        .unwrap();
    let second = s
        .advance_launder_sharma_with_conjugate_wall(
            &g,
            &mut expected,
            &mut expected_wall,
            -0.125,
            h / 2.0,
            load,
            &mut kernel,
            || false,
        )
        .unwrap();
    let mut state = initial;
    let mut wall = initial_wall;
    let mut work = RansConjugateAdaptiveWorkspace::try_new(&g, &wall, 9, 2).unwrap();
    let report = s
        .advance_launder_sharma_conjugate_adaptive(
            &g,
            &mut state,
            &mut wall,
            -0.125,
            h,
            load,
            controls(h, 1.0),
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(report.fluid.accepted_intervals, 1);
    assert_eq!(report.fluid.step_evaluations, 3);
    assert_eq!(bits(&state), bits(&expected));
    assert_eq!(wall, expected_wall);
    let heat = report.fluid.heat.unwrap();
    assert_eq!(
        report.wall.energy_change_j_m2,
        first.wall_energy_change_j_m2 + second.wall_energy_change_j_m2
    );
    assert_eq!(
        report.wall.interface_heat_j_m2,
        (h / 2.0) * first.interface_heat_flux_w_m2 + (h / 2.0) * second.interface_heat_flux_w_m2
    );
    assert_eq!(
        heat.molecular_heat_j_m2,
        (h / 2.0) * first.fluid.molecular_heat_input_w_m2
            + (h / 2.0) * second.fluid.molecular_heat_input_w_m2
    );
    assert_eq!(
        heat.turbulent_heat_j_m2,
        (h / 2.0) * first.fluid.turbulent_heat_input_w_m2
            + (h / 2.0) * second.fluid.turbulent_heat_input_w_m2
    );
    for ((actual, a), b) in rates(report.fluid.mechanical)
        .into_iter()
        .zip(rates(first.fluid.mechanical))
        .zip(rates(second.fluid.mechanical))
    {
        assert_eq!(actual, ((h / 2.0) * a + (h / 2.0) * b) / h);
    }
}

#[test]
fn g0_rejected_trials_do_not_add_exterior_volume_or_dissipative_heat() {
    let g = uniform(8);
    let s = solver();
    let initial = active(&g);
    let mut state = initial.clone();
    let mut wall = wall();
    let mut work = RansConjugateAdaptiveWorkspace::try_new(&g, &wall, 9, 2).unwrap();
    let duration = 0.25;
    let load = RansConjugateHeatLoad {
        exterior_heat_flux_w_m2: 0.7,
        wall_volumetric_heating_w_m3: Some(&[0.2, -0.1]),
        fluid_volumetric_heating_w_m3: Some(&[0.125; 9]),
        coupling: RansHeatCoupling::Dissipative,
        ..RansConjugateHeatLoad::default()
    };
    let report = s
        .advance_launder_sharma_conjugate_adaptive(
            &g,
            &mut state,
            &mut wall,
            -0.125,
            duration,
            load,
            controls(duration, 1e-5),
            &mut work,
            || false,
        )
        .unwrap();
    assert!(report.fluid.rejected_intervals > 0);
    let h = report.fluid.heat.unwrap();
    let w = report.wall;
    close(w.exterior_heat_j_m2, 0.7 * duration, 1e-14);
    close(
        w.volume_heat_j_m2,
        (0.15 * 0.2 - 0.2 * 0.1) * duration,
        1e-14,
    );
    close(h.volume_heat_j_m2, 0.125 * duration, 1e-14);
    close(
        h.energy_change_j_m2 + w.energy_change_j_m2,
        w.exterior_heat_j_m2
            + w.volume_heat_j_m2
            + h.volume_heat_j_m2
            + h.molecular_heat_j_m2
            + h.turbulent_heat_j_m2
            + w.combined_thermal_balance_defect_j_m2,
        1e-11,
    );
    let m = report.fluid.mechanical;
    let mechanical_change = s.rho
        * (m.mean_energy_after - m.mean_energy_before + m.turbulent_energy_after
            - m.turbulent_energy_before);
    close(
        mechanical_change + h.energy_change_j_m2 + w.energy_change_j_m2,
        s.rho
            * duration
            * (m.pressure_power - m.temporal_dissipation - m.turbulent_wall_loss
                + m.mean_balance_defect
                + m.turbulent_balance_defect)
            + w.exterior_heat_j_m2
            + w.volume_heat_j_m2
            + h.volume_heat_j_m2
            + w.combined_thermal_balance_defect_j_m2,
        2e-11,
    );
}

#[test]
fn g4_every_poll_rolls_back_the_complete_window_and_retries_exactly() {
    let g = uniform(4);
    let s = solver();
    let initial = active(&g);
    let initial_wall = wall();
    let h = 1.0 / 64.0;
    let c = controls(h, 1.0);
    let load = RansConjugateHeatLoad {
        exterior_heat_flux_w_m2: 0.1,
        coupling: RansHeatCoupling::Dissipative,
        ..RansConjugateHeatLoad::default()
    };
    let mut work = RansConjugateAdaptiveWorkspace::try_new(&g, &initial_wall, 5, 2).unwrap();
    let mut expected = initial.clone();
    let mut expected_wall = initial_wall.clone();
    let mut polls = 0;
    let expected_report = s
        .advance_launder_sharma_conjugate_adaptive(
            &g,
            &mut expected,
            &mut expected_wall,
            -0.125,
            2.0 * h,
            load,
            c,
            &mut work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(expected_report.fluid.accepted_intervals, 2);
    let bytes = work.numerical_bytes();
    for stop in 1..=polls {
        let mut state = initial.clone();
        let mut wall = initial_wall.clone();
        let mut seen = 0;
        let result = s.advance_launder_sharma_conjugate_adaptive(
            &g,
            &mut state,
            &mut wall,
            -0.125,
            2.0 * h,
            load,
            c,
            &mut work,
            || {
                seen += 1;
                seen == stop
            },
        );
        assert_eq!(
            result,
            Err(RansAdaptiveError::Transport(RansTransportError::Cancelled)),
            "poll {stop}/{polls}"
        );
        assert_eq!(bits(&state), bits(&initial));
        assert_eq!(wall, initial_wall);
        let retried = s
            .advance_launder_sharma_conjugate_adaptive(
                &g,
                &mut state,
                &mut wall,
                -0.125,
                2.0 * h,
                load,
                c,
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(bits(&state), bits(&expected));
        assert_eq!(wall, expected_wall);
        assert_eq!(retried, expected_report);
        assert_eq!(work.numerical_bytes(), bytes);
    }
}

#[test]
fn g4_wall_material_work_and_refinement_limits_refuse_without_partial_state() {
    let g = uniform(4);
    let s = solver();
    let initial = active(&g);
    let initial_wall = wall();
    assert!(matches!(
        RansConjugateAdaptiveWorkspace::try_new(&g, &initial_wall, 4, 2),
        Err(RansTransportError::NodeBudget { .. })
    ));
    assert!(matches!(
        RansConjugateAdaptiveWorkspace::try_new(&g, &initial_wall, 5, 1),
        Err(RansTransportError::NodeBudget { .. })
    ));
    let mut work = RansConjugateAdaptiveWorkspace::try_new(&g, &initial_wall, 5, 2).unwrap();
    let mut state = initial.clone();
    let mut wall = initial_wall.clone();
    let mut c = controls(1.0 / 64.0, 1.0);
    c.maximum_attempts = 1;
    let e = s
        .advance_launder_sharma_conjugate_adaptive(
            &g,
            &mut state,
            &mut wall,
            0.0,
            1.0 / 32.0,
            RansConjugateHeatLoad::default(),
            c,
            &mut work,
            || false,
        )
        .unwrap_err();
    assert!(matches!(e, RansAdaptiveError::AttemptBudget { elapsed_s, .. } if elapsed_s > 0.0));
    assert_eq!(bits(&state), bits(&initial));
    assert_eq!(wall, initial_wall);
    let mut c = controls(0.25, 1e-12);
    c.minimum_interval_s = 0.25;
    assert!(matches!(
        s.advance_launder_sharma_conjugate_adaptive(
            &g,
            &mut state,
            &mut wall,
            0.0,
            0.25,
            RansConjugateHeatLoad::default(),
            c,
            &mut work,
            || false
        ),
        Err(RansAdaptiveError::MinimumInterval { .. })
    ));
    assert_eq!(bits(&state), bits(&initial));
    assert_eq!(wall, initial_wall);
    let mut changed = initial_wall.layers().to_vec();
    changed[0].conductivity_w_m_k *= 2.0;
    let mut different =
        RansConjugateWall::try_new(&changed, initial_wall.temperatures_k(), 2).unwrap();
    let before = different.clone();
    assert!(matches!(
        s.advance_launder_sharma_conjugate_adaptive(
            &g,
            &mut state,
            &mut different,
            0.0,
            0.01,
            RansConjugateHeatLoad::default(),
            controls(0.01, 1.0),
            &mut work,
            || false
        ),
        Err(RansAdaptiveError::InvalidInput(
            "adaptive wall material binding"
        ))
    ));
    assert_eq!(different, before);
    assert_eq!(bits(&state), bits(&initial));
}
