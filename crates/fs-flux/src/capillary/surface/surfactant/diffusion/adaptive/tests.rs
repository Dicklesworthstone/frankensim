use super::*;

fn law(reference: f64, langmuir: bool) -> SurfactantLaw {
    if langmuir {
        SurfactantLaw::Langmuir {
            clean_tension_n_m: 10.0,
            thermal_energy_j_mol: 0.6,
            saturation_concentration_mol_m2: reference,
        }
    } else {
        SurfactantLaw::IdealDilute {
            clean_tension_n_m: 10.0,
            thermal_energy_j_mol: 0.6,
            reference_concentration_mol_m2: reference,
        }
    }
}

fn edge(i: usize, j: usize, k: f64) -> SurfacePartitionDiffusionEdge {
    SurfacePartitionDiffusionEdge {
        faces: [i, j],
        conductance_mol_s: k,
    }
}

fn controls(span: f64, absolute: f64) -> AdaptiveSurfaceDiffusionControls {
    AdaptiveSurfaceDiffusionControls {
        initial_step_s: span,
        minimum_step_s: span * 1e-10,
        maximum_step_s: span,
        absolute_concentration_tolerance_mol_m2: absolute,
        relative_tolerance: 0.0,
        maximum_trials: 10_000,
        maximum_pair_exchanges: 1_000_000,
        maximum_faces: 3,
        maximum_edges: 3,
        amount_tolerance_mol: 1e-12,
        energy_increase_tolerance_j: 1e-12,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} vs {expected:e}; tolerance {tolerance:e}"
    );
}

#[test]
fn adaptive_noncommuting_partition_graph_converges_to_independent_spectrum() {
    let laws = [law(2.0, false), law(1.0, false), law(2.0, false)];
    let areas = [0.5, 2.0, 1.5];
    let capacities = [1.0, 2.0, 3.0];
    let initial = [0.8, 0.2, 0.9];
    // k_ij = rate*C_i*C_j/sum(C) gives dq_i/dt=rate*(mean-q_i).
    // The pair generators do not commute, despite the simple exact spectrum.
    for rate in [1.0, 4096.0] {
        let edges = [
            edge(0, 1, rate / 3.0),
            edge(0, 2, rate / 2.0),
            edge(1, 2, rate),
        ];
        let duration = 4.0 / rate;
        let mut previous_error = None;
        for absolute in [1e-3, 1e-6] {
            let mut amounts = initial;
            let report = diffuse_partitioned_surface_surfactant_adaptive(
                &laws,
                &areas,
                &mut amounts,
                &edges,
                duration,
                controls(duration, absolute),
                || false,
            )
            .unwrap();
            let error = (0..3)
                .map(|i| {
                    let exact = capacities[i]
                        * (1.9 / 6.0 + (initial[i] / capacities[i] - 1.9 / 6.0) * (-4.0_f64).exp());
                    (amounts[i] - exact).abs()
                })
                .fold(0.0_f64, f64::max);
            if let Some(coarse_error) = previous_error {
                assert!(
                    error * 20.0 < coarse_error,
                    "coarse {coarse_error:e}, fine {error:e}"
                );
                assert!(error < 2e-4, "{error:e}");
            }
            previous_error = Some(error);
            assert!(report.rejected_intervals > 0);
            assert!(report.maximum_accepted_normalized_error <= 1.0);
            assert_eq!(
                report.trials,
                report.accepted_intervals + report.rejected_intervals
            );
            assert_eq!(report.evaluated_pair_exchanges, 15 * report.trials);
            assert_eq!(report.diffusion.substeps, 2 * report.accepted_intervals);
            assert_eq!(
                report.diffusion.pair_exchanges,
                10 * report.accepted_intervals
            );
            near(report.diffusion.amount_defect_mol, 0.0, 3e-14);
            assert!(report.diffusion.free_energy_change_j < 0.0);
            near(
                report.diffusion.free_energy_change_j,
                report.diffusion.free_energy_after_j - report.diffusion.free_energy_before_j,
                3e-13,
            );
        }
    }
}

#[test]
fn stiff_langmuir_fast_slow_modes_remain_positive_conservative_and_accurate() {
    let laws = [law(1.0, true), law(4.0, true), law(2.0, true)];
    let areas = [2.0, 1.0, 1.5];
    let capacities = [2.0, 4.0, 3.0];
    let initial = [0.1, 3.2, 0.6];
    let (fast_k, slow_k) = (1000.0, 2.0);
    let edges = [edge(0, 1, fast_k), edge(1, 2, slow_k)];
    let mut amounts = initial;
    let report = diffuse_partitioned_surface_surfactant_adaptive(
        &laws,
        &areas,
        &mut amounts,
        &edges,
        1.0,
        controls(1.0, 1e-8),
        || false,
    )
    .unwrap();
    // The weighted three-node chain has one conserved mode and two rates.
    // Recover both exact amplitudes from d(0) and A*d(0), independently of
    // any edge composition or numerical reference trajectory.
    let trace: f64 = fast_k * (1.0 / 2.0 + 1.0 / 4.0) + slow_k * (1.0 / 4.0 + 1.0 / 3.0);
    let product: f64 = fast_k * slow_k * (1.0 / 8.0 + 1.0 / 6.0 + 1.0 / 12.0);
    let lambda_fast = 0.5 * (trace + (trace * trace - 4.0 * product).sqrt());
    let lambda_slow = product / lambda_fast;
    let mean = initial.iter().sum::<f64>() / capacities.iter().sum::<f64>();
    let deviation = core::array::from_fn::<_, 3, _>(|i| initial[i] - capacities[i] * mean);
    let q = core::array::from_fn::<_, 3, _>(|i| deviation[i] / capacities[i]);
    let action = [
        fast_k * (q[0] - q[1]),
        fast_k * (q[1] - q[0]) + slow_k * (q[1] - q[2]),
        slow_k * (q[2] - q[1]),
    ];
    for i in 0..3 {
        let slow = (lambda_fast * deviation[i] - action[i]) / (lambda_fast - lambda_slow);
        let fast = deviation[i] - slow;
        let exact =
            capacities[i] * mean + slow * (-lambda_slow).exp() + fast * (-lambda_fast).exp();
        near(amounts[i], exact, 2e-5);
        assert!(amounts[i] > 0.0 && amounts[i] < capacities[i]);
        let occupancy = amounts[i] / capacities[i];
        assert!((0.05..=0.8).contains(&occupancy));
        SurfactantFace {
            law: laws[i],
            amount_mol: amounts[i],
        }
        .at_area(areas[i])
        .unwrap();
    }
    near(amounts.iter().sum(), initial.iter().sum(), 2e-12);
    assert!(report.diffusion.free_energy_dissipation_j > 0.0);
    assert!(report.rejected_intervals > 0);
}

#[test]
fn equilibrium_empty_graph_and_homogeneous_recovery_keep_the_original_maps() {
    let laws = [law(1.0, false), law(4.0, false)];
    let areas = [2.0, 1.0];
    let initial = [0.2, 0.4];
    let mut amounts = initial;
    let report = diffuse_partitioned_surface_surfactant_adaptive(
        &laws,
        &areas,
        &mut amounts,
        &[edge(0, 1, 1e6)],
        1.0,
        controls(1.0, 1e-12),
        || false,
    )
    .unwrap();
    assert_eq!(amounts, initial);
    assert_eq!(report.maximum_accepted_normalized_error, 0.0);
    assert_eq!(report.diffusion.free_energy_change_j, 0.0);
    assert_eq!(report.accepted_intervals, 1);
    let mut limits = controls(1.0, 1e-12);
    limits.maximum_pair_exchanges = 0;
    let report = diffuse_partitioned_surface_surfactant_adaptive(
        &laws,
        &areas,
        &mut amounts,
        &[],
        3.0,
        limits,
        || false,
    )
    .unwrap();
    assert_eq!(amounts, initial);
    assert_eq!(report.trials, 0);
    assert_eq!(report.diffusion.duration_s, 3.0);

    let common = law(2.0, false);
    let mut uniform = [0.8, 0.2];
    let mut partitioned = uniform;
    let mut limits = controls(1.0, 1e-12);
    limits.minimum_step_s = 0.1;
    // A short final interval is not stretched to the caller's minimum step.
    let report = diffuse_surface_surfactant_adaptive(
        common,
        &areas,
        &mut uniform,
        &[SurfaceDiffusionEdge {
            faces: [0, 1],
            conductance_m2_s: 0.3,
        }],
        0.01,
        limits,
        || false,
    )
    .unwrap();
    diffuse_partitioned_surface_surfactant_adaptive(
        &[common; 2],
        &areas,
        &mut partitioned,
        &[edge(0, 1, 0.6)],
        0.01,
        limits,
        || false,
    )
    .unwrap();
    assert_eq!(uniform, partitioned);
    assert_eq!(report.diffusion.duration_s, 0.01);
    let mean = 1.0 / 3.0;
    let exact = 2.0 * (mean + (0.4 - mean) * (-0.0045_f64).exp());
    near(uniform[0], exact, 2e-15);
}

#[test]
fn temporal_work_and_late_cancellation_refusals_preserve_the_entire_inventory() {
    let laws = [law(1.0, false), law(2.0, false), law(3.0, false)];
    let areas = [1.0; 3];
    let initial = [0.8, 0.1, 0.2];
    let edges = [edge(0, 1, 1.0), edge(1, 2, 2.0)];
    let mut limits = controls(1.0, 1e-30);
    limits.minimum_step_s = 1.0;
    let mut amounts = initial;
    assert!(matches!(
        diffuse_partitioned_surface_surfactant_adaptive(
            &laws,
            &areas,
            &mut amounts,
            &edges,
            1.0,
            limits,
            || false,
        ),
        Err(AdaptiveSurfaceDiffusionError::TemporalAccuracy { .. })
    ));
    assert_eq!(amounts, initial);

    for pair_budget in [false, true] {
        let mut limits = controls(0.1, 1.0);
        if pair_budget {
            limits.maximum_pair_exchanges = 17;
        } else {
            limits.maximum_trials = 1;
        }
        let error = diffuse_partitioned_surface_surfactant_adaptive(
            &laws,
            &areas,
            &mut amounts,
            &edges,
            0.3,
            limits,
            || false,
        )
        .unwrap_err();
        if pair_budget {
            assert_eq!(
                error,
                AdaptiveSurfaceDiffusionError::PairExchangeBudget {
                    completed: 9,
                    required: 9
                }
            );
        } else {
            assert_eq!(
                error,
                AdaptiveSurfaceDiffusionError::TrialBudget {
                    completed_trials: 1
                }
            );
        }
        assert_eq!(amounts, initial);
    }
    let mut calls = 0;
    diffuse_partitioned_surface_surfactant_adaptive(
        &laws,
        &areas,
        &mut amounts,
        &edges,
        0.3,
        controls(0.1, 1e-5),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    assert_ne!(amounts, initial);
    amounts = initial;
    let mut current = 0;
    let result = diffuse_partitioned_surface_surfactant_adaptive(
        &laws,
        &areas,
        &mut amounts,
        &edges,
        0.3,
        controls(0.1, 1e-5),
        || {
            current += 1;
            current == calls
        },
    );
    assert_eq!(
        result,
        Err(AdaptiveSurfaceDiffusionError::Diffusion(
            SurfaceDiffusionError::Cancelled
        ))
    );
    assert_eq!(amounts, initial);
    let mut incompatible = laws;
    incompatible[1] = law(2.0, true);
    assert!(matches!(
        diffuse_partitioned_surface_surfactant_adaptive(
            &incompatible,
            &areas,
            &mut amounts,
            &edges,
            0.3,
            controls(0.1, 1e-5),
            || false,
        ),
        Err(AdaptiveSurfaceDiffusionError::Diffusion(
            SurfaceDiffusionError::InvalidInput(_)
        ))
    ));
    assert_eq!(amounts, initial);
}

#[test]
fn concentration_error_control_is_invariant_under_joint_area_amount_conductance_scaling() {
    let laws = [law(2.0, false), law(1.0, false), law(2.0, false)];
    let areas = [0.5, 2.0, 1.5];
    let initial = [0.8, 0.2, 0.9];
    let edges = [edge(0, 1, 1.0 / 3.0), edge(0, 2, 0.5), edge(1, 2, 1.0)];
    let mut amounts = initial;
    let mut limits = controls(2.0, 1e-7);
    limits.relative_tolerance = 1e-5;
    let report = diffuse_partitioned_surface_surfactant_adaptive(
        &laws,
        &areas,
        &mut amounts,
        &edges,
        2.0,
        limits,
        || false,
    )
    .unwrap();
    let mut scaled = initial.map(|x| 8.0 * x);
    limits.amount_tolerance_mol *= 8.0;
    limits.energy_increase_tolerance_j *= 8.0;
    let scaled_report = diffuse_partitioned_surface_surfactant_adaptive(
        &laws,
        &areas.map(|x| 8.0 * x),
        &mut scaled,
        &edges.map(|e| SurfacePartitionDiffusionEdge {
            faces: e.faces,
            conductance_mol_s: 8.0 * e.conductance_mol_s,
        }),
        2.0,
        limits,
        || false,
    )
    .unwrap();
    assert_eq!(report.trials, scaled_report.trials);
    assert_eq!(report.rejected_intervals, scaled_report.rejected_intervals);
    assert_eq!(amounts, scaled.map(|x| x / 8.0));
}
