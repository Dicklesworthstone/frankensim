//! Binary stiff exchange on the same physical owner and adaptive timeline.
use super::*;

fn select(system: OsmoticNetwork) -> OsmoticNetwork {
    system
        .with_implicit_binary_diffusion(ImplicitBinaryDiffusionControls::default(), 1 << 22)
        .unwrap()
}

#[test]
fn stiff_pair_matches_independent_unequal_volume_backward_euler() {
    let model = source(1.0, vec![0.0], vec![200.0]);
    let right_volume = 3.0 * VOLUME;
    let right = OsmoticCompartment::new(
        AqueousSolution::nonideal(model.clone(), right_volume, vec![0.0]).unwrap(),
        OsmoticPressureVolume::new(right_volume, 1e5, 1e8).unwrap(),
        0.1 * right_volume,
    )
    .unwrap();
    let mut system = select(
        OsmoticNetwork::with_solution_diffusion(
            vec!["solute-0".into()],
            vec![cell(&model, &[1e-3]), right],
            vec![SolutionDiffusionLink::new(0, 1, 0.01, 0.001, law(&model)).unwrap()],
            physical(),
            1 << 22,
        )
        .unwrap(),
    );
    let saved = bits(&system);
    let dt = 1e5;
    let conductance = 0.01 * 1e-9 / 0.001;
    let stiffness = dt * conductance * (1.0 / VOLUME + 1.0 / right_volume);
    let equilibrium = VOLUME * right_volume / (VOLUME + right_volume) * 1e-3;
    let reference = equilibrium * stiffness / (1.0 + stiffness);
    let mut solver = system.implicit_diffusion.take().unwrap();
    let mut evaluations = 0;
    let transfer = solver
        .transfer(
            &system.cells[0].solution,
            &system.cells[1].solution,
            &mut system.diffusion[0],
            dt,
            physical().constitutive,
            &mut evaluations,
            &mut || false,
        )
        .unwrap();
    // Independent dilute-limit BE reference. The retained finite-solvent
    // correction is at most VBAR*1e-3 and changes this stiff root by <1e-19 mol.
    assert!(stiffness > 1000.0);
    assert!(transfer > 0.0 && transfer < equilibrium);
    assert!((transfer - reference).abs() < 1e-19);
    let nl = 1e-3 * VOLUME;
    let (cl, cr) = ((nl - transfer) / VOLUME, transfer / right_volume);
    let mean = (cl - cr) / (cl / cr).ln();
    let integrated = dt * conductance * (cl - cr) * (1.0 + VBAR * mean);
    let tolerance = 1e-20 + 1e-10 * equilibrium;
    assert!(((nl - (nl - transfer)) - integrated).abs() <= 2.0 * tolerance);
    assert!((transfer - integrated).abs() <= 2.0 * tolerance);
    assert!(evaluations > 1 && evaluations <= 97);
    let reversed = solver
        .transfer(
            &system.cells[1].solution,
            &system.cells[0].solution,
            &mut system.diffusion[0],
            dt,
            physical().constitutive,
            &mut evaluations,
            &mut || false,
        )
        .unwrap();
    assert!((reversed + transfer).abs() <= tolerance);
    assert_eq!(bits(&system), saved);
}

#[test]
fn implicit_window_crosses_stiff_relaxation_with_original_physical_gates() {
    let model = source(1.0, vec![0.0], vec![200.0]);
    let build = || {
        OsmoticNetwork::with_solution_diffusion(
            vec!["solute-0".into()],
            vec![cell(&model, &[1e-3]), cell(&model, &[0.0])],
            vec![SolutionDiffusionLink::new(0, 1, 0.01, 0.001, law(&model)).unwrap()],
            SolutionDiffusionControls {
                maximum_species_fraction: 0.9,
                ..physical()
            },
            1 << 22,
        )
        .unwrap()
    };
    let dt = 1e6;
    let temporal = OsmoticStepControls {
        maximum_step_s: dt,
        maximum_steps: 4,
        maximum_rejections: 4,
        absolute_amount_tolerance_mol: 1e-20,
        ..numerical(1e-3)
    };
    let mut explicit = build();
    let old = bits(&explicit);
    assert_eq!(explicit.implicit_binary_diffusion_controls(), None);
    assert!(matches!(
        explicit.advance(dt, temporal, || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&explicit), old);
    let mut system = select(build());
    let report = system.advance(dt, temporal, || false).unwrap();
    let lambda = 2.0 * 0.01 * 1e-9 / (0.001 * VOLUME);
    let reference = 0.5e-3 * (1.0 + (1.0 + 0.5 * dt * lambda).powi(-2));
    assert_eq!(report.accepted_steps, 1);
    assert_eq!(report.rejected_steps, 0);
    assert!(report.residual_evaluations <= 3 * 97);
    assert!((concentration(&system, 0, 0) - reference).abs() < 2e-13);
    assert!((concentration(&system, 0, 0) - 0.5e-3).abs() < 1e-11);
    assert!(report.solute_balance_residual_mol[0].abs() <= 1e-16);
    assert_eq!(report.water_balance_residual_m3, 0.0);
    assert!(report.diffusion_free_energy_change_j.unwrap() < 0.0);

    // Implicit stability never overrides the caller's withdrawal allowance.
    let mut bounded = select(network(&model, &[&[1e-3], &[0.0]], &[(0, 1, 0.01)]));
    let old = bits(&bounded);
    assert!(matches!(
        bounded.advance(dt, temporal, || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&bounded), old);
}

#[test]
fn nonlinear_implicit_window_refines_to_original_solvent_frame_dynamics() {
    let b = 0.005;
    let model = source(1.0, vec![b], vec![200.0]);
    // Independent ODE includes BOTH nonideal chemical force and the finite
    // solvent-frame factor. This is not a backward-Euler implementation copy.
    let rate = |left: f64| {
        let right = 120.0 - left;
        let mean = (left - right) / (left / right).ln();
        -0.01 / VOLUME * 1e-9 * (left - right) * (1.0 + b * mean) * (1.0 + VBAR * mean) / 0.001
    };
    let mut reference = 100.0;
    let h = 20.0 / 4000.0;
    for _ in 0..4000 {
        let a = rate(reference);
        let b = rate(reference + 0.5 * h * a);
        let c = rate(reference + 0.5 * h * b);
        let d = rate(reference + h * c);
        reference += h * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    let mut errors = Vec::new();
    for relative in [1e-3, 1e-5] {
        let mut system = select(network(&model, &[&[100.0], &[20.0]], &[(0, 1, 0.01)]));
        let report = system.advance(20.0, numerical(relative), || false).unwrap();
        assert!(report.maximum_scaled_error <= 1.0);
        assert!(report.diffusion_free_energy_change_j.unwrap() < 0.0);
        assert!(report.solute_balance_residual_mol[0].abs() <= 1e-16);
        errors.push((concentration(&system, 0, 0) - reference).abs());
    }
    assert!(
        errors[0] < 0.5 && errors[1] < 0.25 * errors[0],
        "{errors:?}"
    );
}

#[test]
fn whole_window_work_and_late_cancellation_are_atomic_and_retryable() {
    let model = source(1.0, vec![0.0], vec![200.0]);
    let empty = || network(&model, &[&[0.0], &[0.0]], &[(0, 1, 0.01)]);
    let mut system = empty()
        .with_implicit_binary_diffusion(
            ImplicitBinaryDiffusionControls {
                maximum_evaluations: 2,
                ..ImplicitBinaryDiffusionControls::default()
            },
            1 << 22,
        )
        .unwrap();
    let old = bits(&system);
    // Each empty pair needs one genuine constitutive sample; the third sweep
    // must consume the SAME complete-window budget, not a fresh per-root cap.
    assert!(matches!(
        system.advance(1.0, numerical(1e-3), || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&system), old);
    system = select(system);
    let report = system.advance(1.0, numerical(1e-3), || false).unwrap();
    assert_eq!(report.residual_evaluations, 3);
    assert_eq!(report.diffusion_free_energy_change_j, Some(0.0));
    assert!(
        system
            .compartments()
            .iter()
            .all(|c| c.solution().solute_moles() == [0.0])
    );

    let build = || select(network(&model, &[&[100.0], &[20.0]], &[(0, 1, 0.01)]));
    let mut completed = build();
    let mut polls = 0;
    completed
        .advance(1.0, numerical(1e-3), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut cancelled = build();
    let old = bits(&cancelled);
    let mut current = 0;
    assert!(matches!(
        cancelled.advance(1.0, numerical(1e-3), || {
            current += 1;
            current == polls
        }),
        Err(OsmoticNetworkError::Cancelled)
    ));
    assert_eq!(bits(&cancelled), old);
    cancelled.advance(1.0, numerical(1e-3), || false).unwrap();
    assert_eq!(bits(&cancelled), bits(&completed));
}

#[test]
fn nonlinear_and_extreme_exposure_refusals_never_publish_an_invented_root() {
    let model = source(1.0, vec![0.005], vec![200.0]);
    let mut system = network(&model, &[&[100.0], &[20.0]], &[(0, 1, 0.01)])
        .with_implicit_binary_diffusion(
            ImplicitBinaryDiffusionControls {
                absolute_residual_tolerance_mol: 1e-30,
                relative_residual_tolerance: 0.0,
                maximum_iterations: 1,
                ..ImplicitBinaryDiffusionControls::default()
            },
            1 << 22,
        )
        .unwrap();
    let old = bits(&system);
    assert!(matches!(
        system.advance(10.0, numerical(1e-3), || false),
        Err(OsmoticNetworkError::ImplicitDiffusionConvergence { .. })
    ));
    assert_eq!(bits(&system), old);
    let model = source(1.0, vec![0.0], vec![200.0]);
    let mut extreme = select(network(&model, &[&[100.0], &[20.0]], &[(0, 1, 1e100)]));
    let old = bits(&extreme);
    assert!(matches!(
        extreme.advance(
            1e300,
            OsmoticStepControls {
                maximum_step_s: 1e300,
                ..numerical(1e-3)
            },
            || false,
        ),
        Err(OsmoticNetworkError::ImplicitDiffusionConvergence { .. })
    ));
    assert_eq!(bits(&extreme), old);
}

#[test]
fn binary_selection_is_explicit_and_counts_reusable_original_solution_storage() {
    let model = source(1.0, vec![0.0], vec![200.0]);
    let system = network(&model, &[&[100.0], &[20.0]], &[(0, 1, 0.01)]);
    let explicit_bytes = system.workspace_bytes();
    assert!(matches!(
        system.clone().with_implicit_binary_diffusion(
            ImplicitBinaryDiffusionControls::default(),
            explicit_bytes
        ),
        Err(OsmoticNetworkError::Workspace { .. })
    ));
    let implicit = select(system);
    let bytes = implicit.workspace_bytes();
    assert!(bytes > explicit_bytes);
    assert_eq!(select(implicit).workspace_bytes(), bytes);
    let two = source(1.0, vec![0.0; 4], vec![200.0; 2]);
    assert!(matches!(
        network(&two, &[&[100.0, 1.0], &[20.0, 1.0]], &[(0, 1, 0.01)])
            .with_implicit_binary_diffusion(ImplicitBinaryDiffusionControls::default(), 1 << 22),
        Err(OsmoticNetworkError::Invalid(_))
    ));
}
