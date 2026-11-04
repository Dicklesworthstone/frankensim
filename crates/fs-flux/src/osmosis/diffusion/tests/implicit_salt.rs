//! Independent salt equations, original source refresh and whole-window gates.
use super::*;
use crate::osmosis::coupled::CoupledSolutionControls;

fn policy(valences: [i32; 2]) -> ZeroCurrentDiffusion {
    ZeroCurrentDiffusion::new(valences.to_vec(), 1e-15, 1e-16).unwrap()
}

fn select(system: OsmoticNetwork) -> OsmoticNetwork {
    system
        .with_implicit_salt_diffusion(ImplicitSaltDiffusionControls::default(), 1 << 22)
        .unwrap()
}

fn binary(model: &Arc<VirialSolutionModel>, left: [f64; 2], right: [f64; 2]) -> OsmoticNetwork {
    network(model, &[&left, &right], &[(0, 1, 0.01)])
        .with_zero_current_diffusion(policy([1, -1]), 1 << 22)
        .unwrap()
}

#[test]
fn stiff_salt_matches_independent_unequal_volume_ambipolar_backward_euler() {
    let model = source(1.0, vec![0.0; 4], vec![200.0; 2]);
    let right_volume = 3.0 * VOLUME;
    let right = OsmoticCompartment::new(
        AqueousSolution::nonideal(model.clone(), right_volume, vec![2e-4 * right_volume; 2])
            .unwrap(),
        OsmoticPressureVolume::new(right_volume, 1e5, 1e8).unwrap(),
        0.1 * right_volume,
    )
    .unwrap();
    let mut system = select(
        OsmoticNetwork::with_solution_diffusion(
            vec!["solute-0".into(), "solute-1".into()],
            vec![cell(&model, &[1e-3; 2]), right],
            vec![SolutionDiffusionLink::new(0, 1, 0.01, 0.001, law(&model)).unwrap()],
            physical(),
            1 << 22,
        )
        .unwrap()
        .with_zero_current_diffusion(policy([1, -1]), 1 << 22)
        .unwrap(),
    );
    let saved = bits(&system);
    let mut implicit = system.implicit_salt.take().unwrap();
    let mut transfers = [f64::NAN; 2];
    let mut evaluations = 0;
    let dt = 1e5;
    implicit
        .transfers_into(
            &system.cells[0].solution,
            &system.cells[1].solution,
            &mut system.diffusion[0],
            dt,
            physical().constitutive,
            &mut evaluations,
            &mut || false,
            &mut transfers,
        )
        .unwrap();
    // Independent dilute Nernst--Hartley D=2*D+*D-/(D++D-), followed by
    // unequal-volume BE. The finite-solvent correction below is retained by
    // production and changes this stiff dilute root by less than 1e-19 mol.
    let conductance = 0.01 * 1.5e-9 / 0.001;
    let stiffness = dt * conductance * (1.0 / VOLUME + 1.0 / right_volume);
    let equilibrium = VOLUME * right_volume / (VOLUME + right_volume) * 8e-4;
    let reference = equilibrium * stiffness / (1.0 + stiffness);
    assert!(stiffness > 1000.0);
    assert_eq!(transfers[0], transfers[1]);
    assert!(transfers[0] > 0.0 && transfers[0] < equilibrium);
    assert!((transfers[0] - reference).abs() < 1e-19);
    let (nl, nr) = (1e-3 * VOLUME, 2e-4 * right_volume);
    let (next_left, next_right) = (nl - transfers[0], nr + transfers[0]);
    let (a, b) = (next_left / VOLUME, next_right / right_volume);
    let mean = (a - b) / (a / b).ln();
    let integrated = dt * conductance * (a - b) * (1.0 + 2.0 * VBAR * mean);
    let tolerance = 1e-20 + 1e-10 * transfers[0];
    assert!(((nl - next_left) - integrated).abs() <= 2.0 * tolerance);
    assert!(((next_right - nr) - integrated).abs() <= 2.0 * tolerance);
    assert!(evaluations > 1 && evaluations <= 97);

    let mut reverse = [f64::NAN; 2];
    implicit
        .transfers_into(
            &system.cells[1].solution,
            &system.cells[0].solution,
            &mut system.diffusion[0],
            dt,
            physical().constitutive,
            &mut evaluations,
            &mut || false,
            &mut reverse,
        )
        .unwrap();
    for i in 0..2 {
        assert!((reverse[i] + transfers[i]).abs() <= 2.0 * tolerance);
    }
    assert_eq!(bits(&system), saved);
}

#[test]
fn nonideal_mixed_valences_follow_one_extent_and_both_original_amount_equations() {
    let interactions = [0.004, 0.001, 0.001, 0.003];
    let model = source(1.0, interactions.to_vec(), vec![200.0; 2]);
    // Opposite species ordering and reducible valences exercise the declared
    // identities. No salt name, dissociation rule or neutral species is inferred.
    for (valences, nu) in [
        ([2, -1], [1.0, 2.0]),
        ([3, -2], [2.0, 3.0]),
        ([-1, 2], [2.0, 1.0]),
        ([4, -2], [1.0, 2.0]),
    ] {
        let left = nu.map(|n| 20.0 * n);
        let right = nu.map(|n| 4.0 * n);
        let electrical = policy(valences);
        let mut system = select(
            network(&model, &[&left, &right], &[(0, 1, 0.01)])
                .with_zero_current_diffusion(electrical.clone(), 1 << 22)
                .unwrap(),
        );
        let mut implicit = system.implicit_salt.take().unwrap();
        let mut transfers = [f64::NAN; 2];
        let mut evaluations = 0;
        let dt = 100.0;
        implicit
            .transfers_into(
                &system.cells[0].solution,
                &system.cells[1].solution,
                &mut system.diffusion[0],
                dt,
                physical().constitutive,
                &mut evaluations,
                &mut || false,
                &mut transfers,
            )
            .unwrap();
        let extent = transfers[0] / nu[0];
        assert!((transfers[1] / nu[1] - extent).abs() <= 2e-20);
        let next_left = std::array::from_fn::<_, 2, _>(|i| left[i] * VOLUME - transfers[i]);
        let next_right = std::array::from_fn::<_, 2, _>(|i| right[i] * VOLUME + transfers[i]);
        electrical.admit_inventories(&next_left).unwrap();
        electrical.admit_inventories(&next_right).unwrap();

        // Independently eliminate charge and electric force in the two-ion
        // friction equations. Ion/ion friction cancels on a pure salt extent;
        // the original finite-solvent and virial terms remain explicitly.
        let (a, b) = (20.0 - extent / VOLUME, 4.0 + extent / VOLUME);
        let mean = (a - b) / (a / b).ln();
        let sum = nu[0] + nu[1];
        let ambipolar = sum / (nu[0] / 1e-9 + nu[1] / 3e-9);
        let virial = nu[0] * nu[0] * interactions[0]
            + 2.0 * nu[0] * nu[1] * interactions[1]
            + nu[1] * nu[1] * interactions[3];
        let predicted = dt * 0.01 * ambipolar * (a - b) / 0.001
            * (1.0 + sum * VBAR * mean)
            * (1.0 + mean * virial / sum);
        for i in 0..2 {
            let tolerance = 1e-20 + 1e-10 * transfers[i].abs();
            assert!(
                ((left[i] * VOLUME - next_left[i]) - nu[i] * predicted).abs() <= 2.0 * tolerance
            );
            assert!(
                ((next_right[i] - right[i] * VOLUME) - nu[i] * predicted).abs() <= 2.0 * tolerance
            );
        }
    }
}

#[test]
fn stiff_salt_uses_original_adaptive_timeline_and_physical_withdrawal_gates() {
    let model = source(1.0, vec![0.0; 4], vec![200.0; 2]);
    let build = || {
        OsmoticNetwork::with_solution_diffusion(
            vec!["solute-0".into(), "solute-1".into()],
            vec![cell(&model, &[1e-3; 2]), cell(&model, &[2e-4; 2])],
            vec![SolutionDiffusionLink::new(0, 1, 0.01, 0.001, law(&model)).unwrap()],
            SolutionDiffusionControls {
                maximum_species_fraction: 0.9,
                ..physical()
            },
            1 << 22,
        )
        .unwrap()
        .with_zero_current_diffusion(policy([1, -1]), 1 << 22)
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
    assert!(matches!(
        explicit.advance(dt, temporal, || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&explicit), old);

    let mut system = select(build());
    let report = system.advance(dt, temporal, || false).unwrap();
    let lambda = 2.0 * 0.01 * 1.5e-9 / (0.001 * VOLUME);
    let reference = 6e-4 + 4e-4 * (1.0 + 0.5 * dt * lambda).powi(-2);
    assert_eq!(report.accepted_steps, 1);
    assert_eq!(report.rejected_steps, 0);
    assert!(report.residual_evaluations <= 3 * 97);
    assert_eq!(system.time_s(), dt);
    for i in 0..2 {
        assert!((concentration(&system, 0, i) - reference).abs() < 2e-13);
        assert!(report.solute_balance_residual_mol[i].abs() <= 1e-16);
    }
    assert_eq!(report.water_balance_residual_m3, 0.0);
    assert!(report.diffusion_free_energy_change_j.unwrap() < 0.0);
    for cell in system.compartments() {
        policy([1, -1])
            .admit_inventories(cell.solution().solute_moles())
            .unwrap();
    }

    let mut bounded = select(binary(&model, [1e-3; 2], [2e-4; 2]));
    let old = bits(&bounded);
    assert!(matches!(
        bounded.advance(
            dt,
            OsmoticStepControls {
                maximum_rejections: 0,
                ..temporal
            },
            || false
        ),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&bounded), old);
}

#[test]
fn implicit_salt_requeries_original_nonlinear_diffusivity_at_actual_candidate_states() {
    let model = source(1.0, vec![0.0; 4], vec![300.0; 2]);
    let material = super::source::multicomponent_material(&model);
    let mut system = select(
        OsmoticNetwork::with_zero_current_solution_transport(
            vec!["solute-0".into(), "solute-1".into()],
            vec![cell(&model, &[100.0; 2]), cell(&model, &[20.0; 2])],
            Vec::new(),
            vec![
                SolutionDiffusionLink::from_material_sources(0, 1, 0.01, 0.001, material).unwrap(),
            ],
            CoupledSolutionControls {
                diffusion: physical(),
                water_balance_tolerance_m3: 1e-18,
            },
            policy([1, -1]),
            1 << 22,
        )
        .unwrap(),
    );
    let saved = bits(&system);
    let mut implicit = system.implicit_salt.take().unwrap();
    let mut transfers = [f64::NAN; 2];
    let mut evaluations = 0;
    let dt = 50.0;
    implicit
        .transfers_into(
            &system.cells[0].solution,
            &system.cells[1].solution,
            &mut system.diffusion[0],
            dt,
            physical().constitutive,
            &mut evaluations,
            &mut || false,
            &mut transfers,
        )
        .unwrap();
    let (a, b) = (100.0 - transfers[0] / VOLUME, 20.0 + transfers[0] / VOLUME);
    let mean = (a - b) / (a / b).ln();
    let mobility = 1e-9 * (1.0 + mean / 100.0);
    let original_mobility = 1e-9 * (1.0 + (80.0 / 5.0_f64.ln()) / 100.0);
    let expected = dt * 0.01 * mobility * (a - b) * (1.0 + 2.0 * VBAR * mean) / 0.001;
    assert!((mobility - original_mobility).abs() > 1e-11);
    for transfer in transfers {
        assert!((transfer - expected).abs() <= 2.0 * (1e-20 + 1e-10 * transfer.abs()));
    }
    // The thermodynamic domain admits this input, but its original mobility
    // card does not. A source refusal cannot publish or overwrite output.
    let unsupported = cell(&model, &[201.0; 2]);
    let mut unchanged = [7.0, 9.0];
    assert!(matches!(
        implicit.transfers_into(
            unsupported.solution(),
            &system.cells[1].solution,
            &mut system.diffusion[0],
            dt,
            physical().constitutive,
            &mut evaluations,
            &mut || false,
            &mut unchanged,
        ),
        Err(OsmoticNetworkError::SolutionDiffusion(_))
    ));
    assert_eq!(unchanged, [7.0, 9.0]);
    assert_eq!(bits(&system), saved);
}

#[test]
fn salt_whole_window_work_late_cancellation_and_extreme_exposure_are_transactional() {
    let model = source(1.0, vec![0.0; 4], vec![200.0; 2]);
    let mut limited = binary(&model, [100.0; 2], [20.0; 2])
        .with_implicit_salt_diffusion(
            ImplicitSaltDiffusionControls {
                maximum_evaluations: 1,
                ..ImplicitSaltDiffusionControls::default()
            },
            1 << 22,
        )
        .unwrap();
    let old = bits(&limited);
    assert!(matches!(
        limited.advance(1.0, numerical(1e-3), || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&limited), old);

    let build = || select(binary(&model, [100.0; 2], [20.0; 2]));
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

    let mut extreme = select(
        network(&model, &[&[100.0; 2], &[20.0; 2]], &[(0, 1, 1e100)])
            .with_zero_current_diffusion(policy([1, -1]), 1 << 22)
            .unwrap(),
    );
    let old = bits(&extreme);
    assert!(matches!(
        extreme.advance(
            1e300,
            OsmoticStepControls {
                maximum_step_s: 1e300,
                ..numerical(1e-3)
            },
            || false
        ),
        Err(OsmoticNetworkError::ImplicitDiffusionConvergence { .. })
    ));
    assert_eq!(bits(&extreme), old);
}

#[test]
fn exact_empty_salt_stays_empty_and_singular_one_sided_fronts_are_not_floored() {
    let model = source(1.0, vec![0.0; 4], vec![200.0; 2]);
    let mut empty = binary(&model, [0.0; 2], [0.0; 2])
        .with_implicit_salt_diffusion(
            ImplicitSaltDiffusionControls {
                maximum_evaluations: 2,
                ..ImplicitSaltDiffusionControls::default()
            },
            1 << 22,
        )
        .unwrap();
    let old = bits(&empty);
    // Each of the full/two-half sweeps owes its own original constitutive
    // sample, even at equilibrium, against ONE complete-window work allowance.
    assert!(matches!(
        empty.advance(1.0, numerical(1e-3), || false),
        Err(OsmoticNetworkError::WorkBudget)
    ));
    assert_eq!(bits(&empty), old);
    empty = select(empty);
    let report = empty.advance(1.0, numerical(1e-3), || false).unwrap();
    assert_eq!(report.residual_evaluations, 3);
    assert_eq!(report.diffusion_free_energy_change_j, Some(0.0));
    for cell in empty.compartments() {
        assert_eq!(cell.solution().solute_moles(), &[0.0; 2]);
    }

    let mut singular = select(binary(&model, [0.0; 2], [20.0; 2]));
    let old = bits(&singular);
    assert!(matches!(
        singular.advance(1.0, numerical(1e-3), || false),
        Err(OsmoticNetworkError::SolutionDiffusion(
            SolutionDiffusionError::SingularDiffusionPotential
        ))
    ));
    assert_eq!(bits(&singular), old);
}

#[test]
fn salt_selection_is_explicit_bounded_and_refuses_independent_compositions() {
    let model = source(1.0, vec![0.0; 4], vec![200.0; 2]);
    let explicit = binary(&model, [100.0; 2], [20.0; 2]);
    let bytes = explicit.workspace_bytes();
    assert!(matches!(
        explicit
            .clone()
            .with_implicit_salt_diffusion(ImplicitSaltDiffusionControls::default(), bytes),
        Err(OsmoticNetworkError::Workspace { .. })
    ));
    let selected = select(explicit);
    let bytes = selected.workspace_bytes();
    assert_eq!(select(selected).workspace_bytes(), bytes);
    assert!(matches!(
        network(&model, &[&[100.0; 2], &[20.0; 2]], &[(0, 1, 0.01)])
            .with_implicit_salt_diffusion(ImplicitSaltDiffusionControls::default(), 1 << 22),
        Err(OsmoticNetworkError::Invalid(_))
    ));
    for controls in [
        ImplicitSaltDiffusionControls {
            relative_residual_tolerance: 0.01,
            ..ImplicitSaltDiffusionControls::default()
        },
        ImplicitSaltDiffusionControls {
            maximum_iterations: 0,
            ..ImplicitSaltDiffusionControls::default()
        },
        ImplicitSaltDiffusionControls {
            maximum_evaluations: 0,
            ..ImplicitSaltDiffusionControls::default()
        },
    ] {
        assert!(matches!(
            binary(&model, [100.0; 2], [20.0; 2]).with_implicit_salt_diffusion(controls, 1 << 22),
            Err(OsmoticNetworkError::Invalid(_))
        ));
    }
    let model = source(1.0, vec![0.0; 9], vec![200.0; 3]);
    let material = SolutionMaxwellStefan::new(
        model.clone(),
        (SpeciesId::new("water").unwrap(), MolarMass::new(0.018)),
        (0..3)
            .map(|i| {
                (
                    SpeciesId::new(format!("solute-{i}")).unwrap(),
                    MolarMass::new(0.04),
                )
            })
            .collect(),
        (0..16)
            .map(|i| if i / 4 == i % 4 { 0.0 } else { 1e-9 })
            .collect(),
    )
    .unwrap();
    let extra_species = OsmoticNetwork::with_zero_current_solution_transport(
        (0..3).map(|i| format!("solute-{i}")).collect(),
        vec![
            cell(&model, &[100.0, 100.0, 0.0]),
            cell(&model, &[20.0, 20.0, 0.0]),
        ],
        Vec::new(),
        vec![SolutionDiffusionLink::new(0, 1, 0.01, 0.001, material).unwrap()],
        CoupledSolutionControls {
            diffusion: physical(),
            water_balance_tolerance_m3: 1e-18,
        },
        ZeroCurrentDiffusion::new(vec![1, -1, 0], 1e-15, 1e-16).unwrap(),
        1 << 22,
    )
    .unwrap();
    assert!(matches!(
        extra_species
            .with_implicit_salt_diffusion(ImplicitSaltDiffusionControls::default(), 1 << 22),
        Err(OsmoticNetworkError::Invalid(_))
    ));
}
