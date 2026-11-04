//! Stiff globally coupled ionic evolution through original physical owners.
use super::*;

fn nonlinear() -> ImplicitElectroneutralDiffusionControls {
    ImplicitElectroneutralDiffusionControls {
        nonlinear: ImplicitMulticomponentDiffusionControls {
            absolute_residual_tolerance_mol: 1e-16,
            ..ImplicitMulticomponentDiffusionControls::default()
        },
        ..ImplicitElectroneutralDiffusionControls::default()
    }
}

fn stiff_triangle() -> OsmoticNetwork {
    let model = source(1.0, vec![0.0; 25], vec![200.0; 5]);
    let ids = (0..5)
        .map(|i| {
            (
                SpeciesId::new(format!("solute-{i}")).unwrap(),
                MolarMass::new(0.04),
            )
        })
        .collect();
    let mut coefficients = vec![1e-9; 36];
    for i in 0..6 {
        coefficients[i * 6 + i] = 0.0;
    }
    let law = SolutionMaxwellStefan::new(
        model.clone(),
        (SpeciesId::new("water").unwrap(), MolarMass::new(0.018)),
        ids,
        coefficients,
    )
    .unwrap();
    let rows = [
        [60.0, 30.0, 90.0, 20.0, 0.0],
        [40.0, 20.0, 60.0, 80.0, 0.0],
        [20.0, 10.0, 30.0, 140.0, 0.0],
    ];
    OsmoticNetwork::with_solution_diffusion(
        (0..5).map(|i| format!("solute-{i}")).collect(),
        rows.iter().map(|row| cell(&model, row)).collect(),
        [(0, 1), (0, 2), (1, 2)]
            .into_iter()
            .map(|(a, b)| SolutionDiffusionLink::new(a, b, 0.01, 0.001, law.clone()).unwrap())
            .collect(),
        SolutionDiffusionControls {
            maximum_species_fraction: 0.95,
            ..physical()
        },
        1 << 22,
    )
    .unwrap()
    .with_electroneutral_diffusion(
        ZeroCurrentDiffusion::new(vec![1, 1, -1, 0, 0], 1e-13, 1e-16).unwrap(),
        ElectroneutralDiffusionControls::default(),
        1 << 22,
    )
    .unwrap()
}

#[test]
fn stiff_global_backward_euler_matches_finite_concentration_triangle_eigenmodes() {
    let mut system = stiff_triangle()
        .with_implicit_electroneutral_diffusion(nonlinear(), 1 << 22)
        .unwrap();
    let before: Vec<_> = system
        .cells
        .iter()
        .map(|cell| cell.solution.solute_moles().to_vec())
        .collect();
    let controls = system.diffusion_controls.unwrap();
    let mut implicit = system.implicit_electroneutral.take().unwrap();
    let mut electrical = system.electroneutral_diffusion.take().unwrap();
    let mut evaluations = 0;
    let dt = 1e4;
    implicit
        .sweep(
            &mut electrical,
            &mut system.cells,
            &mut system.diffusion,
            dt,
            controls,
            &mut evaluations,
            &mut || false,
        )
        .unwrap();
    // Proportional charged profiles retain EACH logarithmic chain rule, and
    // neutral counterdiffusion keeps total concentration constant at 200.
    // Equal original binary coefficients therefore give exactly Fickian flux
    // with zero solvent motion and zero electric field at finite concentration.
    // The triangle Laplacian has both nonconstant eigenvalues equal to 3.
    let k = dt * 0.01 * 1e-9 / (0.001 * VOLUME);
    assert!(k > 10.0);
    for species in 0..5 {
        let mean = before.iter().map(|row| row[species] / VOLUME).sum::<f64>() / 3.0;
        let total = before.iter().map(|row| row[species]).sum::<f64>();
        for node in 0..3 {
            let original = before[node][species] / VOLUME;
            let expected = mean + (original - mean) / (1.0 + 3.0 * k);
            let actual = concentration(&system, node, species);
            assert!(
                (actual - expected).abs() < 2e-7,
                "node={node} species={species}: {actual:e} vs {expected:e}"
            );
            let sum_other = (0..3)
                .filter(|&other| other != node)
                .map(|other| concentration(&system, other, species) - actual)
                .sum::<f64>();
            let residual = VOLUME * (actual - original - k * sum_other);
            assert!(residual.abs() < 2e-14, "original ion BE row {residual:e}");
        }
        let after = system
            .cells
            .iter()
            .map(|cell| cell.solution.solute_moles()[species])
            .sum::<f64>();
        assert!((after - total).abs() <= 1e-16);
    }
    for cell in &system.cells {
        assert_eq!(cell.solution.solute_moles()[4], 0.0);
        let n = cell.solution.solute_moles();
        assert!((n[0] + n[1] - n[2]).abs() <= 1e-16);
    }
    assert!(electrical.potentials.iter().all(|value| value.abs() < 1e-7));
    assert!(evaluations > 6 && evaluations <= nonlinear().nonlinear.maximum_evaluations);
}

#[test]
fn stiff_global_window_completes_with_native_temporal_and_energy_gates() {
    let dt = 1e4;
    let time = OsmoticStepControls {
        maximum_step_s: dt,
        maximum_steps: 4,
        maximum_rejections: 4,
        ..numerical(1e-2)
    };
    let mut explicit = stiff_triangle();
    let saved = bits(&explicit);
    assert!(matches!(
        explicit.advance(dt, time, || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&explicit), saved);
    let mut system = stiff_triangle()
        .with_implicit_electroneutral_diffusion(nonlinear(), 1 << 22)
        .unwrap();
    let report = system.advance(dt, time, || false).unwrap();
    let k = dt * 0.01 * 1e-9 / (0.001 * VOLUME);
    let expected = 40.0 + 20.0 * (1.0 + 1.5 * k).powi(-2);
    assert_eq!(report.accepted_steps, 1);
    assert_eq!(report.rejected_steps, 0);
    assert_eq!(system.time_s(), dt);
    assert!((concentration(&system, 0, 0) - expected).abs() < 2e-7);
    assert!(report.diffusion_free_energy_change_j.unwrap() < 0.0);
    assert!(
        report.solute_balance_residual_mol[..5]
            .iter()
            .all(|r| r.abs() <= 1e-16)
    );
    assert_eq!(report.water_balance_residual_m3, 0.0);
    assert!(report.electroneutral_diffusion.unwrap().potential_solves > 3);
}

fn circulating_triangle() -> OsmoticNetwork {
    let model = source(1.0, vec![0.0; 4], vec![200.0; 2]);
    let links = [(0, 1, 1e-9, 3e-9), (0, 2, 3e-9, 1e-9), (1, 2, 2e-9, 1e-9)]
        .into_iter()
        .map(|(a, b, plus, minus)| {
            let law = SolutionMaxwellStefan::new(
                model.clone(),
                (SpeciesId::new("water").unwrap(), MolarMass::new(0.018)),
                vec![
                    (SpeciesId::new("solute-0").unwrap(), MolarMass::new(0.04)),
                    (SpeciesId::new("solute-1").unwrap(), MolarMass::new(0.06)),
                ],
                vec![0.0, plus, minus, plus, 0.0, 5e-10, minus, 5e-10, 0.0],
            )
            .unwrap();
            SolutionDiffusionLink::new(a, b, 0.01, 0.001, law).unwrap()
        })
        .collect();
    OsmoticNetwork::with_solution_diffusion(
        vec!["solute-0".into(), "solute-1".into()],
        [100.0, 40.0, 15.0].map(|c| cell(&model, &[c, c])).into(),
        links,
        SolutionDiffusionControls {
            maximum_species_fraction: 0.95,
            ..physical()
        },
        1 << 22,
    )
    .unwrap()
    .with_electroneutral_diffusion(
        ZeroCurrentDiffusion::new(vec![1, -1], 1e-13, 1e-16).unwrap(),
        ElectroneutralDiffusionControls::default(),
        1 << 22,
    )
    .unwrap()
    .with_implicit_electroneutral_diffusion(nonlinear(), 1 << 22)
    .unwrap()
}

#[test]
fn global_backward_euler_retains_circulating_current_and_every_original_ion_equation() {
    let mut system = circulating_triangle();
    let original: Vec<_> = system
        .cells
        .iter()
        .flat_map(|cell| cell.solution.solute_moles().iter().copied())
        .collect();
    let mut implicit = system.implicit_electroneutral.take().unwrap();
    let mut electrical = system.electroneutral_diffusion.take().unwrap();
    let physical = system.diffusion_controls.unwrap();
    let mut evaluations = 0;
    let dt = 1000.0;
    implicit
        .sweep(
            &mut electrical,
            &mut system.cells,
            &mut system.diffusion,
            dt,
            physical,
            &mut evaluations,
            &mut || false,
        )
        .unwrap();
    electrical
        .evaluate_transfers(
            &system.cells,
            &mut system.diffusion,
            dt,
            physical,
            &mut evaluations,
            &mut || false,
        )
        .unwrap();
    for (index, &integrated) in electrical.transfers_mol().iter().enumerate() {
        let amount = system.cells[index / 2].solution.solute_moles()[index % 2];
        let change = amount - original[index];
        let tolerance = nonlinear().nonlinear.absolute_residual_tolerance_mol
            + nonlinear().nonlinear.relative_residual_tolerance * change.abs();
        assert!((change - integrated).abs() <= tolerance);
    }
    assert!(electrical.report().maximum_edge_charge_rate_mol_s > 1e-10);
    assert!(electrical.report().maximum_node_charge_rate_mol_s <= 1e-16);
    assert!(
        electrical
            .potentials
            .iter()
            .any(|potential| potential.abs() > 1e-5)
    );
}

#[test]
fn global_implicit_work_refusals_and_late_cancellation_publish_no_state() {
    let mut controls = nonlinear();
    controls.nonlinear.maximum_evaluations = 11;
    let mut bounded = stiff_triangle()
        .with_implicit_electroneutral_diffusion(controls, 1 << 22)
        .unwrap();
    let saved = bits(&bounded);
    assert!(matches!(
        bounded.advance(0.1, numerical(1e-4), || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&bounded), saved);
    let mut reference = stiff_triangle()
        .with_implicit_electroneutral_diffusion(nonlinear(), 1 << 22)
        .unwrap();
    let mut polls = 0;
    reference
        .advance(0.1, numerical(1e-4), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut cancelled = stiff_triangle()
        .with_implicit_electroneutral_diffusion(nonlinear(), 1 << 22)
        .unwrap();
    let saved = bits(&cancelled);
    let mut seen = 0;
    assert!(matches!(
        cancelled.advance(0.1, numerical(1e-4), || {
            seen += 1;
            seen == polls - 1
        }),
        Err(OsmoticNetworkError::Cancelled)
    ));
    assert_eq!(bits(&cancelled), saved);
    let system = stiff_triangle();
    let budget = system.workspace_bytes();
    assert!(matches!(
        system.with_implicit_electroneutral_diffusion(nonlinear(), budget),
        Err(OsmoticNetworkError::Workspace { .. })
    ));
    let mut too_few = nonlinear();
    too_few.maximum_unknowns = 1;
    assert!(
        stiff_triangle()
            .with_implicit_electroneutral_diffusion(too_few, 1 << 22)
            .is_err()
    );
}
