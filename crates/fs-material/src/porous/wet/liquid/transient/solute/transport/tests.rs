use super::*;

fn controls<const N: usize>() -> DiluteSoluteControls<N> {
    DiluteSoluteControls {
        absolute_equation_tolerance_mol_m2: [1e-11; N],
        relative_equation_tolerance: 1e-12,
        amount_balance_tolerance_mol_m2: [1e-11; N],
        water_balance_tolerance_m: 1e-13,
        maximum_cell_species: 131072,
        maximum_storage_bytes: 64 * 1024 * 1024,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn heterogeneous_half_resistances_match_exact_two_volume_backward_euler() {
    let volumes = [0.2, 0.6];
    let widths = [1.0, 2.0];
    let tau = [2.0, 3.0];
    let d = [[0.04, 0.07], [0.08, 0.03]];
    let initial = [[0.4, 0.0], [0.0, 1.8]];
    for dt in [0.01, 1.0, 10000.0] {
        let result = transport(
            [0.0, dt],
            dt,
            &volumes,
            &volumes,
            &widths,
            &[0.0; 3],
            &initial,
            &d,
            &tau,
            [DiluteSoluteBoundary::Outflow; 2],
            controls(),
            || false,
        )
        .unwrap();
        for k in 0..2 {
            let g = 1.0
                / (0.5 / (volumes[0] / widths[0] * d[0][k] / tau[0])
                    + 1.0 / (volumes[1] / widths[1] * d[1][k] / tau[1]));
            let mean = (initial[0][k] + initial[1][k]) / (volumes[0] + volumes[1]);
            let difference = (initial[0][k] / volumes[0] - initial[1][k] / volumes[1])
                / (1.0 + dt * g * (1.0 / volumes[0] + 1.0 / volumes[1]));
            let expected = [
                mean + volumes[1] / 0.8 * difference,
                mean - volumes[0] / 0.8 * difference,
            ];
            for i in 0..2 {
                near(result.concentration_mol_m3[i][k], expected[i], 2e-12);
                assert!(result.amount_mol_m2[i][k] >= 0.0);
            }
            near(
                result.report.diffusive_face_amount_mol_m2[1][k],
                dt * g * (expected[0] - expected[1]),
                2e-12,
            );
            near(result.report.amount_change_mol_m2[k], 0.0, 2e-12);
        }
        assert_eq!(result.report.linear_solves, 2);
    }
}

#[test]
fn constant_concentration_survives_nonuniform_changing_water_storage_and_reversal() {
    let before = [0.25, 0.6, 0.3];
    let water = [0.02, 0.01, 0.015, 0.015];
    let after = std::array::from_fn::<_, 3, _>(|i| before[i] + water[i] - water[i + 1]);
    let c = [3.0, 0.2];
    let initial = before.map(|v| c.map(|x| x * v));
    let d = [[0.02, 0.04], [0.05, 0.01], [0.03, 0.09]];
    let forward = transport(
        [2.0, 3.0],
        1.0,
        &before,
        &after,
        &[1.0; 3],
        &water,
        &initial,
        &d,
        &[1.0, 2.0, 3.0],
        [
            DiluteSoluteBoundary::InflowConcentration(c),
            DiluteSoluteBoundary::Outflow,
        ],
        controls(),
        || false,
    )
    .unwrap();
    let reverse = |values: &[f64]| values.iter().rev().copied().collect::<Vec<_>>();
    let reflected = transport(
        [2.0, 3.0],
        1.0,
        &reverse(&before),
        &reverse(&after),
        &[1.0; 3],
        &water.iter().rev().map(|x| -x).collect::<Vec<_>>(),
        &initial.iter().rev().copied().collect::<Vec<_>>(),
        &d.iter().rev().copied().collect::<Vec<_>>(),
        &[3.0, 2.0, 1.0],
        [
            DiluteSoluteBoundary::Outflow,
            DiluteSoluteBoundary::InflowConcentration(c),
        ],
        controls(),
        || false,
    )
    .unwrap();
    for i in 0..3 {
        for k in 0..2 {
            near(forward.concentration_mol_m3[i][k], c[k], 3e-15);
            near(reflected.concentration_mol_m3[2 - i][k], c[k], 3e-15);
        }
    }
    for k in 0..2 {
        near(
            forward.report.boundary_amount_into_column_mol_m2[k],
            c[k] * (water[0] - water[3]),
            2e-16,
        );
        for flux in &forward.report.diffusive_face_amount_mol_m2 {
            near(flux[k], 0.0, 2e-16);
        }
    }
}

#[test]
fn diffusion_only_reservoirs_resolve_native_half_cell_and_film_resistances() {
    let volume = 0.3;
    let dt = 5.0;
    let d = 0.02;
    let tau = 2.0;
    let left = DiluteSoluteBoundary::Reservoir {
        concentration_mol_m3: [1.0],
        film_resistance_s_m: [0.0],
    };
    let right = DiluteSoluteBoundary::Reservoir {
        concentration_mol_m3: [3.0],
        film_resistance_s_m: [17.0],
    };
    let result = transport(
        [0.0, dt],
        dt,
        &[volume],
        &[volume],
        &[1.0],
        &[0.0; 2],
        &[[0.6]],
        &[[d]],
        &[tau],
        [left, right],
        controls(),
        || false,
    )
    .unwrap();
    let half = 0.5 / (volume * d / tau);
    let g = [1.0 / half, 1.0 / (half + 17.0)];
    let c = (0.6 + dt * (g[0] + 3.0 * g[1])) / (volume + dt * (g[0] + g[1]));
    near(result.concentration_mol_m3[0][0], c, 1e-15);
    near(
        result.report.diffusive_face_amount_mol_m2[0][0],
        dt * g[0] * (1.0 - c),
        1e-16,
    );
    near(
        result.report.diffusive_face_amount_mol_m2[1][0],
        dt * g[1] * (c - 3.0),
        1e-16,
    );
    assert_eq!(result.report.advective_face_amount_mol_m2, vec![[0.0]; 2]);
    assert!(
        result
            .report
            .diffusive_face_amount_mol_m2
            .iter()
            .all(|a| a[0] < 0.0)
    );
}

#[test]
fn implicit_advection_has_no_explicit_courant_limit_and_zero_d_does_not_seal_water() {
    let volume = 0.3;
    let transported_water = 0.6;
    let initial = [[0.0], [0.6], [0.0], [0.0]];
    let result = transport(
        [0.0, 1.0],
        1.0,
        &[volume; 4],
        &[volume; 4],
        &[1.0; 4],
        &[transported_water; 5],
        &initial,
        &[[0.0]; 4],
        &[1.0; 4],
        [
            DiluteSoluteBoundary::InflowConcentration([1.0]),
            DiluteSoluteBoundary::Outflow,
        ],
        controls(),
        || false,
    )
    .unwrap();
    let mut upstream = 1.0;
    for i in 0..4 {
        upstream = (initial[i][0] + transported_water * upstream) / (volume + transported_water);
        near(result.concentration_mol_m3[i][0], upstream, 3e-16);
        assert!((0.0..=2.0).contains(&upstream));
    }
    assert_eq!(result.report.diffusive_face_amount_mol_m2, vec![[0.0]; 5]);
    near(
        result.report.boundary_amount_into_column_mol_m2[0],
        transported_water * (1.0 - upstream),
        1e-16,
    );
}

#[test]
fn original_closed_diffusion_mode_has_first_order_physical_time_refinement() {
    let n = 8;
    let h = 1.0 / n as f64;
    let theta = 0.4;
    let d = 0.1;
    let duration = 1.0;
    let mode = |i: usize| (std::f64::consts::PI * (i as f64 + 0.5) / n as f64).cos();
    let initial: Vec<_> = (0..n)
        .map(|i| [theta * h * (1.0 + 0.2 * mode(i))])
        .collect();
    let lambda = 4.0 * d / (h * h) * (std::f64::consts::PI / (2.0 * n as f64)).sin().powi(2);
    let errors = [4, 8, 16].map(|steps| {
        let dt = duration / steps as f64;
        let mut amounts = initial.clone();
        for step in 0..steps {
            let result = transport(
                [step as f64 * dt, (step + 1) as f64 * dt],
                dt,
                &vec![theta * h; n],
                &vec![theta * h; n],
                &vec![h; n],
                &vec![0.0; n + 1],
                &amounts,
                &vec![[d]; n],
                &vec![1.0; n],
                [DiluteSoluteBoundary::Outflow; 2],
                controls(),
                || false,
            )
            .unwrap();
            amounts = result.amount_mol_m2;
        }
        let discrete_factor = (1.0 + lambda * dt).powi(-steps);
        for i in 0..n {
            near(
                amounts[i][0] / (theta * h),
                1.0 + 0.2 * mode(i) * discrete_factor,
                3e-14,
            );
        }
        (amounts[0][0] / (theta * h) - (1.0 + 0.2 * mode(0) * (-lambda * duration).exp())).abs()
    });
    assert!(
        errors[0] / errors[1] > 1.8 && errors[1] / errors[2] > 1.8,
        "{errors:?}"
    );
}

#[test]
fn positive_trace_stock_and_inward_response_survive_without_species_thresholds() {
    let mut limits = controls();
    limits.absolute_equation_tolerance_mol_m2 = [1e-214];
    limits.amount_balance_tolerance_mol_m2 = [1e-214];
    let result = transport(
        [0.0, 1.0],
        1.0,
        &[0.5; 2],
        &[0.5; 2],
        &[1.0; 2],
        &[0.0; 3],
        &[[1e-200], [0.0]],
        &[[0.1]; 2],
        &[1.0; 2],
        [DiluteSoluteBoundary::Outflow; 2],
        limits,
        || false,
    )
    .unwrap();
    assert!(result.amount_mol_m2[0][0] > 0.0);
    assert!(result.amount_mol_m2[1][0] > 0.0);
    near(result.amount_mol_m2[1][0] / 1e-200, 1.0 / 12.0, 2e-16);
}

#[test]
fn missing_inflow_dry_storage_water_mismatch_and_budgets_refuse_without_mutation() {
    let initial = [[0.2], [0.1]];
    let original = initial;
    let call = |before: &[f64], after: &[f64], water: &[f64], bounds, ctl| {
        transport(
            [0.0, 1.0],
            1.0,
            before,
            after,
            &[1.0; 2],
            water,
            &initial,
            &[[0.1]; 2],
            &[1.0; 2],
            bounds,
            ctl,
            || false,
        )
    };
    assert!(matches!(
        call(
            &[0.3; 2],
            &[0.3; 2],
            &[0.1; 3],
            [DiluteSoluteBoundary::Outflow; 2],
            controls()
        ),
        Err(DiluteSoluteError::Invalid(_))
    ));
    assert!(matches!(
        call(
            &[0.3, 0.0],
            &[0.3, 0.0],
            &[0.0; 3],
            [DiluteSoluteBoundary::Outflow; 2],
            controls()
        ),
        Err(DiluteSoluteError::Invalid(_))
    ));
    assert!(matches!(
        call(
            &[0.3; 2],
            &[0.3, 0.4],
            &[0.0; 3],
            [DiluteSoluteBoundary::Outflow; 2],
            controls()
        ),
        Err(DiluteSoluteError::WaterBalance { .. })
    ));
    let mut ctl = controls();
    ctl.maximum_cell_species = 1;
    assert!(matches!(
        call(
            &[0.3; 2],
            &[0.3; 2],
            &[0.0; 3],
            [DiluteSoluteBoundary::Outflow; 2],
            ctl
        ),
        Err(DiluteSoluteError::WorkBudget { .. })
    ));
    ctl = controls();
    ctl.maximum_storage_bytes = 1;
    assert!(matches!(
        call(
            &[0.3; 2],
            &[0.3; 2],
            &[0.0; 3],
            [DiluteSoluteBoundary::Outflow; 2],
            ctl
        ),
        Err(DiluteSoluteError::Workspace { .. })
    ));
    assert_eq!(initial, original);
}

#[test]
fn cancellation_during_the_actual_solve_discards_the_whole_private_result() {
    let initial = vec![[0.2]; 256];
    let before = initial.clone();
    let mut polls = 0;
    let result = transport(
        [0.0, 1.0],
        1.0,
        &vec![0.3; 256],
        &vec![0.3; 256],
        &vec![1.0; 256],
        &vec![0.0; 257],
        &initial,
        &vec![[0.1]; 256],
        &vec![1.0; 256],
        [DiluteSoluteBoundary::Outflow; 2],
        controls(),
        || {
            polls += 1;
            polls == 25
        },
    );
    assert_eq!(result, Err(DiluteSoluteError::Cancelled));
    assert_eq!(initial, before);
    let recovered = transport(
        [0.0, 1.0],
        1.0,
        &vec![0.3; 256],
        &vec![0.3; 256],
        &vec![1.0; 256],
        &vec![0.0; 257],
        &initial,
        &vec![[0.1]; 256],
        &vec![1.0; 256],
        [DiluteSoluteBoundary::Outflow; 2],
        controls(),
        || false,
    )
    .unwrap();
    for cell in recovered.concentration_mol_m3 {
        near(cell[0], 2.0 / 3.0, 5e-16);
    }
}
