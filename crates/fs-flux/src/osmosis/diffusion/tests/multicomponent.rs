//! Stiff mixture exchange, counterdiffusion and original cross-force dynamics.
use super::*;

fn select(system: OsmoticNetwork) -> OsmoticNetwork {
    system
        .with_implicit_multicomponent_diffusion(
            ImplicitMulticomponentDiffusionControls::default(),
            1 << 22,
        )
        .unwrap()
}

fn counterdiffusion() -> OsmoticNetwork {
    let model = source(1.0, vec![0.0; 9], vec![200.0; 3]);
    let ids = (0..3)
        .map(|i| {
            (
                SpeciesId::new(format!("solute-{i}")).unwrap(),
                MolarMass::new(0.04),
            )
        })
        .collect();
    let mut coefficients = vec![1e-9; 16];
    for i in 0..4 {
        coefficients[i * 4 + i] = 0.0;
    }
    let law = SolutionMaxwellStefan::new(
        model.clone(),
        (SpeciesId::new("water").unwrap(), MolarMass::new(0.018)),
        ids,
        coefficients,
    )
    .unwrap();
    let right_volume = 3.0 * VOLUME;
    let right = OsmoticCompartment::new(
        AqueousSolution::nonideal(
            model.clone(),
            right_volume,
            vec![0.0, 100.0 * right_volume, 0.0],
        )
        .unwrap(),
        OsmoticPressureVolume::new(right_volume, 1e5, 1e8).unwrap(),
        0.1 * right_volume,
    )
    .unwrap();
    OsmoticNetwork::with_solution_diffusion(
        vec!["solute-0".into(), "solute-1".into(), "solute-2".into()],
        vec![cell(&model, &[100.0, 0.0, 0.0]), right],
        vec![SolutionDiffusionLink::new(0, 1, 0.01, 0.001, law).unwrap()],
        SolutionDiffusionControls {
            maximum_species_fraction: 0.9,
            ..physical()
        },
        1 << 22,
    )
    .unwrap()
}

#[test]
fn stiff_multisolute_pair_matches_exact_counterdiffusion_and_keeps_absent_species() {
    let mut system = select(counterdiffusion());
    let saved = bits(&system);
    let mut implicit = system.implicit_multicomponent.take().unwrap();
    let mut evaluations = 0;
    let mut transfer = [f64::NAN; 3];
    let dt = 1e5;
    implicit
        .transfers_into(
            &system.cells[0].solution,
            &system.cells[1].solution,
            &mut system.diffusion[0],
            dt,
            physical().constitutive,
            &mut evaluations,
            &mut transfer,
            &mut || false,
        )
        .unwrap();
    // Equal binary coefficients and opposite gradients give zero solvent
    // molar flux. Counterdiffusion is EXACTLY Fickian at finite concentration,
    // not just a dilute approximation. The third declared solute is absent.
    let lambda = 0.01 * 1e-9 / 0.001 * (1.0 / VOLUME + 1.0 / (3.0 * VOLUME));
    let equilibrium = 75.0 * VOLUME;
    let reference = equilibrium * (dt * lambda) / (1.0 + dt * lambda);
    assert!(dt * lambda > 1000.0);
    assert!((transfer[0] - reference).abs() < 2e-14, "{transfer:?}");
    assert!((transfer[1] + reference).abs() < 2e-14);
    assert_eq!(transfer[2], 0.0);
    let before = 100.0 * VOLUME;
    let left = before - transfer[0];
    let right = transfer[0];
    let integrated = dt * 0.01 * 1e-9 / 0.001 * (left / VOLUME - right / (3.0 * VOLUME));
    assert!(((before - left) - integrated).abs() < 2e-14);
    assert!((right - integrated).abs() < 2e-14);
    assert!(evaluations > 4 && evaluations < 2000);
    assert_eq!(bits(&system), saved);
}

#[test]
fn stiff_multisolute_window_uses_the_original_adaptive_and_physical_gates() {
    let dt = 1e5;
    let numerical = OsmoticStepControls {
        maximum_step_s: dt,
        maximum_steps: 4,
        maximum_rejections: 4,
        ..numerical(1e-2)
    };
    let mut explicit = counterdiffusion();
    let saved = bits(&explicit);
    assert!(matches!(
        explicit.advance(dt, numerical, || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&explicit), saved);
    let mut system = select(counterdiffusion());
    let report = system.advance(dt, numerical, || false).unwrap();
    let lambda = 0.01 * 1e-9 / 0.001 * (1.0 / VOLUME + 1.0 / (3.0 * VOLUME));
    let reference = 25.0 + 75.0 * (1.0 + 0.5 * lambda * dt).powi(-2);
    assert_eq!(report.accepted_steps, 1);
    assert_eq!(report.rejected_steps, 0);
    assert!((concentration(&system, 0, 0) - reference).abs() < 1e-7);
    assert!((concentration(&system, 0, 0) + concentration(&system, 0, 1) - 100.0).abs() < 1e-7);
    assert!(concentration(&system, 0, 1) > 0.0 && concentration(&system, 1, 0) > 0.0);
    assert_eq!(concentration(&system, 0, 2), 0.0);
    assert_eq!(concentration(&system, 1, 2), 0.0);
    assert!(report.diffusion_free_energy_change_j.unwrap() < 0.0);
    assert!(
        report.solute_balance_residual_mol[..3]
            .iter()
            .all(|r| r.abs() <= 1e-16)
    );
    assert_eq!(report.water_balance_residual_m3, 0.0);
}

fn independent_rate(left: [f64; 2]) -> [f64; 2] {
    let right = [41.0 - left[0], 101.0 - left[1]];
    let delta = [left[0] - right[0], left[1] - right[1]];
    let mean = std::array::from_fn::<_, 2, _>(|i| {
        if left[i] == right[i] {
            left[i]
        } else if left[i] == 0.0 || right[i] == 0.0 {
            0.0
        } else {
            delta[i] / (left[i] / right[i]).ln()
        }
    });
    let g = [
        delta[0] + mean[0] * (0.001 * delta[0] + 0.008 * delta[1]),
        delta[1] + mean[1] * (0.008 * delta[0] + 0.1 * delta[1]),
    ];
    // Independent 2x2 solvent-fixed MS equations, including cross friction.
    // A_ii=c_water/D_iw + sum(c_j/D_ij), A_ij=-c_i/D_ij.
    let water = 1.0 / VBAR;
    let total = water + mean[0] + mean[1];
    let a = water / 1e-9 + mean[1] / 5e-10;
    let b = -mean[0] / 5e-10;
    let c = -mean[1] / 5e-10;
    let d = water / 3e-9 + mean[0] / 5e-10;
    let rhs = [total * g[0] / 0.001, total * g[1] / 0.001];
    let determinant = a * d - b * c;
    [
        -0.001 / VOLUME * (d * rhs[0] - b * rhs[1]) / determinant,
        -0.001 / VOLUME * (a * rhs[1] - c * rhs[0]) / determinant,
    ]
}

#[test]
fn native_multicomponent_steps_refine_to_nonideal_uphill_cross_diffusion() {
    let model = source(1.0, vec![0.001, 0.008, 0.008, 0.1], vec![200.0; 2]);
    let mut exact = [20.0, 100.0];
    let h = 20.0 / 4000.0;
    for _ in 0..4000 {
        let a = independent_rate(exact);
        let b = independent_rate(std::array::from_fn(|i| exact[i] + 0.5 * h * a[i]));
        let c = independent_rate(std::array::from_fn(|i| exact[i] + 0.5 * h * b[i]));
        let d = independent_rate(std::array::from_fn(|i| exact[i] + h * c[i]));
        for i in 0..2 {
            exact[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    let mut errors = Vec::new();
    for relative in [1e-3, 1e-5] {
        let mut system = select(network(
            &model,
            &[&[20.0, 100.0], &[21.0, 1.0]],
            &[(0, 1, 0.001)],
        ));
        let report = system.advance(20.0, numerical(relative), || false).unwrap();
        assert!(concentration(&system, 0, 0) < 20.0);
        assert!(concentration(&system, 1, 0) > 21.0);
        assert!(report.diffusion_free_energy_change_j.unwrap() < 0.0);
        assert!(
            report.solute_balance_residual_mol[..2]
                .iter()
                .all(|r| r.abs() <= 1e-16)
        );
        errors.push(
            (0..2)
                .map(|i| (concentration(&system, 0, i) - exact[i]).abs())
                .fold(0.0_f64, f64::max),
        );
    }
    assert!(
        errors[0] < 1.0 && errors[1] < 0.35 * errors[0],
        "{errors:?}"
    );
}

#[test]
fn multicomponent_work_cancellation_and_extreme_exposure_leave_no_partial_state() {
    let model = source(1.0, vec![0.0; 4], vec![200.0; 2]);
    let empty = || network(&model, &[&[0.0, 0.0], &[0.0, 0.0]], &[(0, 1, 0.01)]);
    let mut system = empty()
        .with_implicit_multicomponent_diffusion(
            ImplicitMulticomponentDiffusionControls {
                maximum_evaluations: 2,
                ..ImplicitMulticomponentDiffusionControls::default()
            },
            1 << 22,
        )
        .unwrap();
    let saved = bits(&system);
    assert!(matches!(
        system.advance(1.0, numerical(1e-3), || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&system), saved);
    system = select(system);
    let report = system.advance(1.0, numerical(1e-3), || false).unwrap();
    assert_eq!(report.residual_evaluations, 3);
    assert_eq!(report.diffusion_free_energy_change_j, Some(0.0));

    let build = || {
        select(network(
            &model,
            &[&[100.0, 20.0], &[20.0, 100.0]],
            &[(0, 1, 0.01)],
        ))
    };
    let mut completed = build();
    let mut polls = 0;
    completed
        .advance(1.0, numerical(1e-3), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut cancelled = build();
    let saved = bits(&cancelled);
    let mut current = 0;
    assert!(matches!(
        cancelled.advance(1.0, numerical(1e-3), || {
            current += 1;
            current == polls
        }),
        Err(OsmoticNetworkError::Cancelled)
    ));
    assert_eq!(bits(&cancelled), saved);
    cancelled.advance(1.0, numerical(1e-3), || false).unwrap();
    assert_eq!(bits(&cancelled), bits(&completed));

    let mut extreme = select(network(
        &model,
        &[&[100.0, 20.0], &[20.0, 100.0]],
        &[(0, 1, 1e100)],
    ));
    let saved = bits(&extreme);
    assert!(
        extreme
            .advance(
                1e300,
                OsmoticStepControls {
                    maximum_step_s: 1e300,
                    ..numerical(1e-3)
                },
                || false
            )
            .is_err()
    );
    assert_eq!(bits(&extreme), saved);
}

#[test]
fn multicomponent_native_workspace_and_policy_are_explicit() {
    let system = counterdiffusion();
    let bytes = system.workspace_bytes();
    assert!(matches!(
        system.clone().with_implicit_multicomponent_diffusion(
            ImplicitMulticomponentDiffusionControls::default(),
            bytes,
        ),
        Err(OsmoticNetworkError::Workspace { .. })
    ));
    let configured = select(system);
    assert!(
        configured.workspace_bytes()
            >= bytes + ImplicitMulticomponentDiffusionControls::solver_workspace_bytes()
    );
    let new_bytes = configured.workspace_bytes();
    assert_eq!(select(configured).workspace_bytes(), new_bytes);
    let mut bad = ImplicitMulticomponentDiffusionControls::default();
    bad.newton.relative_tolerance = 0.1;
    assert!(matches!(
        counterdiffusion().with_implicit_multicomponent_diffusion(bad, 1 << 22),
        Err(OsmoticNetworkError::Invalid(_))
    ));
}
