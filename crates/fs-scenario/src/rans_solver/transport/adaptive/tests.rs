use super::*;
use crate::rans_card::RansCardDraft;

fn solver() -> RansSolver {
    let card = RansCardDraft::launder_sharma_channel("adaptive-transport-regression")
        .freeze()
        .unwrap();
    let mut solver = RansSolver::new(card);
    solver.nu_laminar = 0.125;
    solver.rho = 2.0;
    solver.cp = 4.0;
    solver.k_thermal = 1.0;
    solver
}

fn uniform(cells: usize) -> RansChannelGrid {
    RansChannelGrid {
        num_cells: cells,
        half_height_m: 1.0,
        y_nodes: (0..=cells).map(|i| i as f64 / cells as f64).collect(),
    }
}

fn laminar(grid: &RansChannelGrid) -> RansFieldState {
    let n = grid.num_cells + 1;
    RansFieldState {
        u: grid
            .y_nodes
            .iter()
            .map(|y| (std::f64::consts::FRAC_PI_2 * y).sin())
            .collect(),
        k: vec![0.0; n],
        eps: vec![0.0; n],
        nu_t: vec![0.0; n],
        temp_k: grid
            .y_nodes
            .iter()
            .map(|y| 4.0 + (std::f64::consts::PI * y).cos())
            .collect(),
    }
}

fn active(grid: &RansChannelGrid) -> RansFieldState {
    let mut state = laminar(grid);
    for (i, &y) in grid.y_nodes.iter().enumerate() {
        let shape = y * (2.0 - y);
        state.k[i] = 0.02 * shape * shape;
        state.eps[i] = 0.01 * shape;
    }
    state
}

fn controls(interval: f64, absolute: f64) -> RansAdaptiveControls {
    RansAdaptiveControls {
        initial_interval_s: interval,
        minimum_interval_s: 1e-10,
        maximum_interval_s: interval,
        maximum_attempts: 10_000,
        relative_tolerance: 0.0,
        absolute: RansFieldTolerances {
            velocity_m_s: absolute,
            kinetic_energy_m2_s2: absolute,
            dissipation_m2_s3: absolute,
            eddy_viscosity_m2_s: absolute,
            temperature_k: absolute,
        },
        maximum_fractional_change: 0.5,
    }
}

fn bits(state: &RansFieldState) -> [Vec<u64>; 5] {
    [&state.u, &state.k, &state.eps, &state.nu_t, &state.temp_k]
        .map(|field| field.iter().map(|value| value.to_bits()).collect())
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}"
    );
}

// State images swap roles after acceptance, so compare all 24 allocations as
// a set, retaining lengths and capacities as well as their original addresses.
fn allocations(work: &RansAdaptiveWorkspace) -> [(usize, usize, usize); 24] {
    let k = &work.kernel;
    let arrays = [
        &k.velocity,
        &k.kinetic,
        &k.dissipation,
        &k.viscosity,
        &k.production,
        &k.lower,
        &k.diagonal,
        &k.upper,
        &k.rhs,
        &work.accepted.u,
        &work.accepted.k,
        &work.accepted.eps,
        &work.accepted.nu_t,
        &work.accepted.temp_k,
        &work.coarse.u,
        &work.coarse.k,
        &work.coarse.eps,
        &work.coarse.nu_t,
        &work.coarse.temp_k,
        &work.fine.u,
        &work.fine.k,
        &work.fine.eps,
        &work.fine.nu_t,
        &work.fine.temp_k,
    ];
    let mut result = arrays.map(|v| (v.as_ptr() as usize, v.len(), v.capacity()));
    result.sort_unstable();
    result
}

fn rates(m: RansTransportReport) -> [f64; 8] {
    [
        m.pressure_power,
        m.molecular_dissipation,
        m.turbulence_production,
        m.turbulent_destruction,
        m.turbulent_wall_loss,
        m.temporal_dissipation,
        m.mean_balance_defect,
        m.turbulent_balance_defect,
    ]
}

#[test]
fn tighter_local_tolerances_improve_independent_semidiscrete_decay() {
    let grid = uniform(8);
    let solver = solver();
    let initial = laminar(&grid);
    let duration = 1.0;
    let dy = 1.0 / grid.num_cells as f64;
    // The wall/center half-volume operator has sine eigenvectors with
    // theta=pi/(2N) for no-slip/insulated velocity. Both thermal endpoints
    // are insulated, giving cosine theta=pi/N. These eigenvalues are derived
    // from the uniform second difference, independently of the native rows.
    let lambda_u = 4.0 * solver.nu_laminar / (dy * dy)
        * (std::f64::consts::PI / (4.0 * grid.num_cells as f64))
            .sin()
            .powi(2);
    let alpha = solver.k_thermal / (solver.rho * solver.cp);
    let lambda_t = 4.0 * alpha / (dy * dy)
        * (std::f64::consts::PI / (2.0 * grid.num_cells as f64))
            .sin()
            .powi(2);
    let exact_u: Vec<_> = initial
        .u
        .iter()
        .map(|v| v * (-lambda_u * duration).exp())
        .collect();
    let exact_t: Vec<_> = initial
        .temp_k
        .iter()
        .map(|t| 4.0 + (t - 4.0) * (-lambda_t * duration).exp())
        .collect();
    let mut errors = Vec::new();
    let mut intervals = Vec::new();
    for tolerance in [1e-3, 1e-6] {
        let mut state = initial.clone();
        let mut work = RansAdaptiveWorkspace::try_new(&grid, 9).unwrap();
        let report = solver
            .advance_launder_sharma_adaptive(
                &grid,
                &mut state,
                0.0,
                duration,
                Some(RansHeatLoad::default()),
                controls(0.5, tolerance),
                &mut work,
                || false,
            )
            .unwrap();
        let error = |actual: &[f64], expected: &[f64]| {
            actual
                .iter()
                .zip(expected)
                .fold(0.0_f64, |worst, (a, b)| worst.max((a - b).abs()))
        };
        errors.push([error(&state.u, &exact_u), error(&state.temp_k, &exact_t)]);
        intervals.push(report.accepted_intervals);
        assert!(report.maximum_accepted_error_ratio <= 1.0);
        assert!(report.maximum_accepted_fractional_change <= 0.5);
        assert!(report.rejected_intervals > 0);
        assert_eq!(report.duration_s, duration);
        assert!(
            state
                .k
                .iter()
                .chain(&state.eps)
                .chain(&state.nu_t)
                .all(|v| *v == 0.0)
        );
        assert!(state.temp_k.iter().all(|t| *t >= 3.0 && *t <= 5.0));
        let heat = report.heat.unwrap();
        close(heat.energy_change_j_m2, 0.0, 1e-11);
        assert_eq!(heat.wall_heat_j_m2, 0.0);
        assert_eq!(heat.volume_heat_j_m2, 0.0);
        let m = report.mechanical;
        close(
            (m.mean_energy_after - m.mean_energy_before) / duration,
            -m.molecular_dissipation - m.temporal_dissipation,
            2e-12,
        );
    }
    for field in 0..2 {
        assert!(
            errors[1][field] < 0.1 * errors[0][field],
            "errors={errors:?}"
        );
    }
    assert!(
        errors[1][0] < 2e-4 && errors[1][1] < 4e-4,
        "errors={errors:?}"
    );
    assert!(intervals[1] > intervals[0]);
}

#[test]
fn accepted_two_half_steps_and_only_their_budgets_replay_exactly() {
    let grid = uniform(8);
    let solver = solver();
    let initial = active(&grid);
    let h = 1.0 / 64.0;
    let pressure = -0.125;
    let source = vec![0.125; 9];
    let load = RansHeatLoad {
        wall_heat_flux_w_m2: 0.25,
        volumetric_heating_w_m3: Some(&source),
        ..RansHeatLoad::default()
    };
    for heat in [
        None,
        Some(load),
        Some(RansHeatLoad {
            coupling: crate::rans_solver::RansHeatCoupling::Dissipative,
            ..load
        }),
    ] {
        let mut actual = initial.clone();
        let mut work = RansAdaptiveWorkspace::try_new(&grid, 9).unwrap();
        let addresses = allocations(&work);
        assert!(addresses.windows(2).all(|p| p[0].0 != p[1].0));
        assert_eq!(work.numerical_bytes(), 24 * 9 * core::mem::size_of::<f64>());
        let mut expected = initial.clone();
        let mut native = RansTransportWorkspace::try_new(&grid, 9).unwrap();
        let mut mechanical = Vec::new();
        let mut thermal = Vec::new();
        for _ in 0..2 {
            if let Some(load) = heat {
                let report = solver
                    .advance_launder_sharma_with_heat(
                        &grid,
                        &mut expected,
                        pressure,
                        h / 2.0,
                        load,
                        &mut native,
                        || false,
                    )
                    .unwrap();
                mechanical.push(report.mechanical);
                thermal.push(report);
            } else {
                mechanical.push(
                    solver
                        .advance_launder_sharma(
                            &grid,
                            &mut expected,
                            pressure,
                            h / 2.0,
                            &mut native,
                            || false,
                        )
                        .unwrap(),
                );
            }
        }
        let mut coarse = initial.clone();
        if let Some(load) = heat {
            solver
                .advance_launder_sharma_with_heat(
                    &grid,
                    &mut coarse,
                    pressure,
                    h,
                    load,
                    &mut native,
                    || false,
                )
                .unwrap();
        } else {
            solver
                .advance_launder_sharma(&grid, &mut coarse, pressure, h, &mut native, || false)
                .unwrap();
        }
        assert_ne!(bits(&coarse), bits(&expected));
        let report = solver
            .advance_launder_sharma_adaptive(
                &grid,
                &mut actual,
                pressure,
                h,
                heat,
                controls(h, 0.1),
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(bits(&actual), bits(&expected));
        assert_eq!(allocations(&work), addresses);
        assert_eq!(
            (
                report.attempts,
                report.accepted_intervals,
                report.rejected_intervals,
                report.step_evaluations
            ),
            (1, 1, 0, 3)
        );
        assert_eq!(report.minimum_accepted_interval_s, h);
        assert_eq!(report.maximum_accepted_interval_s, h);
        assert_eq!(report.mechanical.time_step_s, h);
        assert_eq!(
            report.mechanical.mean_energy_before,
            mechanical[0].mean_energy_before
        );
        assert_eq!(
            report.mechanical.mean_energy_after,
            mechanical[1].mean_energy_after
        );
        assert_eq!(
            report.mechanical.turbulent_energy_before,
            mechanical[0].turbulent_energy_before
        );
        assert_eq!(
            report.mechanical.turbulent_energy_after,
            mechanical[1].turbulent_energy_after
        );
        for ((actual, first), second) in rates(report.mechanical)
            .into_iter()
            .zip(rates(mechanical[0]))
            .zip(rates(mechanical[1]))
        {
            // Binary intervals make this independent two-term average exact.
            assert_eq!(actual.to_bits(), (0.5 * (first + second)).to_bits());
        }
        if let Some(heat) = report.heat {
            assert_eq!(
                heat.energy_before_j_m2,
                thermal[0].thermal_energy_before_j_m2
            );
            assert_eq!(heat.energy_after_j_m2, thermal[1].thermal_energy_after_j_m2);
            assert_eq!(
                heat.energy_change_j_m2,
                thermal[0].thermal_energy_change_j_m2 + thermal[1].thermal_energy_change_j_m2
            );
            assert_eq!(heat.wall_heat_j_m2, h * load.wall_heat_flux_w_m2);
            assert_eq!(heat.volume_heat_j_m2, h * source[0] * grid.half_height_m);
            assert_eq!(
                heat.molecular_heat_j_m2,
                (h / 2.0)
                    * (thermal[0].molecular_heat_input_w_m2 + thermal[1].molecular_heat_input_w_m2)
            );
            assert_eq!(
                heat.turbulent_heat_j_m2,
                (h / 2.0)
                    * (thermal[0].turbulent_heat_input_w_m2 + thermal[1].turbulent_heat_input_w_m2)
            );
            assert_eq!(
                heat.balance_defect_j_m2,
                (h / 2.0)
                    * (thermal[0].thermal_balance_defect_w_m2
                        + thermal[1].thermal_balance_defect_w_m2)
            );
            close(
                heat.energy_change_j_m2
                    - heat.wall_heat_j_m2
                    - heat.volume_heat_j_m2
                    - heat.molecular_heat_j_m2
                    - heat.turbulent_heat_j_m2,
                heat.balance_defect_j_m2,
                2e-14,
            );
        } else {
            assert_eq!(actual.temp_k, initial.temp_k);
        }
        let m = report.mechanical;
        close(
            (m.mean_energy_after - m.mean_energy_before) / h,
            m.pressure_power
                - m.molecular_dissipation
                - m.turbulence_production
                - m.temporal_dissipation,
            2e-13,
        );
        close(
            (m.turbulent_energy_after - m.turbulent_energy_before) / h,
            m.turbulence_production - m.turbulent_destruction - m.turbulent_wall_loss,
            2e-13,
        );
    }
    // Now force actual trial rejections under nonzero wall and volume inputs.
    // The integrated imposed heat must still cover exactly the requested time.
    let mut state = initial.clone();
    let mut work = RansAdaptiveWorkspace::try_new(&grid, 9).unwrap();
    let report = solver
        .advance_launder_sharma_adaptive(
            &grid,
            &mut state,
            pressure,
            h,
            Some(load),
            controls(h, 1e-7),
            &mut work,
            || false,
        )
        .unwrap();
    assert!(report.rejected_intervals > 0);
    let heat = report.heat.unwrap();
    close(heat.wall_heat_j_m2, h * load.wall_heat_flux_w_m2, 1e-15);
    close(
        heat.volume_heat_j_m2,
        h * source[0] * grid.half_height_m,
        1e-15,
    );
    close(
        heat.energy_change_j_m2 - heat.wall_heat_j_m2 - heat.volume_heat_j_m2,
        heat.balance_defect_j_m2,
        2e-14,
    );
    let m = report.mechanical;
    close(
        (m.mean_energy_after - m.mean_energy_before) / h,
        m.pressure_power
            - m.molecular_dissipation
            - m.turbulence_production
            - m.temporal_dissipation,
        2e-12,
    );
}

#[test]
fn active_turbulence_advances_from_an_unpopulated_eddy_cache() {
    let grid = uniform(8);
    let solver = solver();
    let initial = active(&grid);
    assert!(initial.nu_t.iter().all(|v| *v == 0.0));
    let mut state = initial.clone();
    let mut work = RansAdaptiveWorkspace::try_new(&grid, 9).unwrap();
    let mut c = controls(1.0 / 1024.0, 1e-5);
    c.relative_tolerance = 0.01;
    c.absolute.eddy_viscosity_m2_s = 1e-16;
    let report = solver
        .advance_launder_sharma_adaptive(
            &grid,
            &mut state,
            -0.125,
            1.0 / 256.0,
            None,
            c,
            &mut work,
            || false,
        )
        .unwrap();
    assert!(report.accepted_intervals >= 4);
    assert_ne!(state.k, initial.k);
    assert_ne!(state.eps, initial.eps);
    assert_eq!(state.temp_k, initial.temp_k);
    for i in 1..=grid.num_cells {
        assert!(state.k[i] > 0.0 && state.eps[i] > 0.0 && state.nu_t[i] > 0.0);
        let q = state.k[i] * (state.k[i] / state.eps[i]);
        let re_t = q / solver.nu_laminar;
        let expected = solver.coeffs.c_mu * (-3.4 / (1.0 + re_t / 50.0).powi(2)).exp() * q;
        close(state.nu_t[i], expected, 2e-14 * expected);
    }
    assert_eq!(state.u[0], 0.0);
    assert_eq!(state.k[0], 0.0);
    assert_eq!(state.eps[0], 0.0);
    assert_eq!(state.nu_t[0], 0.0);
}

#[test]
fn stiff_decay_aliasing_is_refused_even_when_the_error_estimate_passes() {
    let grid = uniform(8);
    let solver = solver();
    let mut initial = laminar(&grid);
    for (u, y) in initial.u.iter_mut().zip(&grid.y_nodes) {
        *u = 10.0 * (15.0 * std::f64::consts::FRAC_PI_2 * y).sin();
    }
    let h = 1024.0;
    let mut coarse = initial.clone();
    let mut fine = initial.clone();
    let mut native = RansTransportWorkspace::try_new(&grid, 9).unwrap();
    solver
        .advance_launder_sharma(&grid, &mut coarse, 0.0, h, &mut native, || false)
        .unwrap();
    for _ in 0..2 {
        solver
            .advance_launder_sharma(&grid, &mut fine, 0.0, h / 2.0, &mut native, || false)
            .unwrap();
    }
    // Both trial endpoints have almost erased the mode; agreement by itself
    // is misleading even at a permissive absolute allowance of 1 m/s.
    let error = coarse
        .u
        .iter()
        .zip(&fine.u)
        .fold(0.0_f64, |a, (u, v)| a.max((u - v).abs()));
    assert!(error < 0.01);
    let mut c = controls(h, 1.0);
    c.minimum_interval_s = h;
    let mut state = initial.clone();
    let mut work = RansAdaptiveWorkspace::try_new(&grid, 9).unwrap();
    let refusal = solver.advance_launder_sharma_adaptive(
        &grid,
        &mut state,
        0.0,
        h,
        None,
        c,
        &mut work,
        || false,
    );
    assert!(
        matches!(refusal, Err(RansAdaptiveError::MinimumInterval { elapsed_s: 0.0, interval_s: 1024.0, indicator }) if indicator > 1.0)
    );
    assert_eq!(bits(&state), bits(&initial));
}

#[test]
fn attempt_limits_and_minimum_interval_preserve_the_entire_window() {
    let grid = uniform(4);
    let solver = solver();
    let initial = laminar(&grid);
    let h = 1.0 / 64.0;
    let mut work = RansAdaptiveWorkspace::try_new(&grid, 5).unwrap();
    let addresses = allocations(&work);
    for maximum in [0, 1] {
        let mut state = initial.clone();
        let mut c = controls(h, 0.1);
        c.maximum_attempts = maximum;
        let result = solver.advance_launder_sharma_adaptive(
            &grid,
            &mut state,
            0.0,
            2.0 * h,
            Some(RansHeatLoad::default()),
            c,
            &mut work,
            || false,
        );
        assert_eq!(
            result,
            Err(RansAdaptiveError::AttemptBudget {
                attempts: maximum,
                elapsed_s: maximum as f64 * h
            })
        );
        assert_eq!(bits(&state), bits(&initial));
        assert_eq!(allocations(&work), addresses);
    }
    let mut state = initial.clone();
    let mut c = controls(h, 1e-16);
    c.minimum_interval_s = h;
    let result = solver.advance_launder_sharma_adaptive(
        &grid,
        &mut state,
        0.0,
        2.0 * h,
        Some(RansHeatLoad::default()),
        c,
        &mut work,
        || false,
    );
    assert!(matches!(
        result,
        Err(RansAdaptiveError::MinimumInterval { elapsed_s: 0.0, .. })
    ));
    assert_eq!(bits(&state), bits(&initial));
    // A final endpoint remainder may be shorter than the retry lower bound.
    let mut c = controls(h, 0.1);
    c.minimum_interval_s = h;
    let report = solver
        .advance_launder_sharma_adaptive(
            &grid,
            &mut state,
            0.0,
            1.5 * h,
            Some(RansHeatLoad::default()),
            c,
            &mut work,
            || false,
        )
        .unwrap();
    assert_eq!(report.duration_s, 1.5 * h);
    assert_eq!(report.accepted_intervals, 2);
    assert_eq!(report.minimum_accepted_interval_s, h / 2.0);
    assert_eq!(report.maximum_accepted_interval_s, h);
    assert_eq!(allocations(&work), addresses);
}

#[test]
fn every_observed_cancellation_poll_rolls_back_and_allows_identical_retry() {
    let grid = uniform(4);
    let solver = solver();
    let initial = laminar(&grid);
    let h = 1.0 / 64.0;
    let c = controls(h, 0.1);
    let heat = Some(RansHeatLoad {
        wall_heat_flux_w_m2: 0.125,
        volumetric_heating_w_m3: None,
        ..RansHeatLoad::default()
    });
    let mut work = RansAdaptiveWorkspace::try_new(&grid, 5).unwrap();
    let addresses = allocations(&work);
    let mut expected = initial.clone();
    let mut polls = 0;
    let expected_report = solver
        .advance_launder_sharma_adaptive(
            &grid,
            &mut expected,
            -0.125,
            2.0 * h,
            heat,
            c,
            &mut work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 10);
    assert_eq!(expected_report.accepted_intervals, 2);
    for stop in 1..=polls {
        let mut state = initial.clone();
        let mut observed = 0;
        let result = solver.advance_launder_sharma_adaptive(
            &grid,
            &mut state,
            -0.125,
            2.0 * h,
            heat,
            c,
            &mut work,
            || {
                observed += 1;
                observed == stop
            },
        );
        assert_eq!(
            result,
            Err(RansAdaptiveError::Transport(RansTransportError::Cancelled)),
            "poll {stop}"
        );
        assert_eq!(bits(&state), bits(&initial), "poll {stop}");
        let retry = solver
            .advance_launder_sharma_adaptive(
                &grid,
                &mut state,
                -0.125,
                2.0 * h,
                heat,
                c,
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(bits(&state), bits(&expected), "retry after poll {stop}");
        assert_eq!(retry, expected_report);
        assert_eq!(allocations(&work), addresses);
    }
}

#[test]
fn material_scope_shapes_node_budget_and_finite_controls_are_admitted_before_publication() {
    let grid = uniform(4);
    assert!(matches!(
        RansAdaptiveWorkspace::try_new(&grid, 4),
        Err(RansTransportError::NodeBudget {
            needed: 5,
            maximum: 4
        })
    ));
    let initial = laminar(&grid);
    let h = 1.0 / 64.0;
    let mut solver = solver();
    solver.cp = f64::NAN;
    solver.k_thermal = -1.0;
    solver.pr_turbulent = f64::NAN;
    let mut work = RansAdaptiveWorkspace::try_new(&grid, 5).unwrap();
    let mut state = initial.clone();
    let report = solver
        .advance_launder_sharma_adaptive(
            &grid,
            &mut state,
            0.0,
            h,
            None,
            controls(h, 0.1),
            &mut work,
            || false,
        )
        .unwrap();
    assert!(report.heat.is_none());
    assert_eq!(state.temp_k, initial.temp_k);
    let accepted = bits(&state);
    let result = solver.advance_launder_sharma_adaptive(
        &grid,
        &mut state,
        0.0,
        h,
        Some(RansHeatLoad::default()),
        controls(h, 0.1),
        &mut work,
        || false,
    );
    assert!(matches!(
        result,
        Err(RansAdaptiveError::Transport(
            RansTransportError::InvalidInput(_)
        ))
    ));
    assert_eq!(bits(&state), accepted);
    let solver = self::solver();
    for mode in 0..7 {
        let mut state = initial.clone();
        let mut c = controls(h, 0.1);
        match mode {
            0 => c.initial_interval_s = f64::NAN,
            1 => c.minimum_interval_s = 0.0,
            2 => c.maximum_interval_s = h / 2.0,
            3 => c.absolute.temperature_k = f64::INFINITY,
            4 => c.relative_tolerance = -1.0,
            5 => c.maximum_fractional_change = 0.75,
            _ => c.maximum_attempts = usize::MAX,
        }
        let result = solver.advance_launder_sharma_adaptive(
            &grid,
            &mut state,
            0.0,
            h,
            None,
            c,
            &mut work,
            || false,
        );
        assert!(matches!(result, Err(RansAdaptiveError::InvalidInput(_))));
        assert_eq!(bits(&state), bits(&initial));
    }
    let mut malformed = initial.clone();
    malformed.eps.pop();
    let before = bits(&malformed);
    assert!(matches!(
        solver.advance_launder_sharma_adaptive(
            &grid,
            &mut malformed,
            0.0,
            h,
            None,
            controls(h, 0.1),
            &mut work,
            || false
        ),
        Err(RansAdaptiveError::InvalidInput(_))
    ));
    assert_eq!(bits(&malformed), before);
    let other = uniform(8);
    let mut state = laminar(&other);
    let before = bits(&state);
    assert!(matches!(
        solver.advance_launder_sharma_adaptive(
            &other,
            &mut state,
            0.0,
            h,
            None,
            controls(h, 0.1),
            &mut work,
            || false
        ),
        Err(RansAdaptiveError::InvalidInput(_))
    ));
    assert_eq!(bits(&state), before);
}
