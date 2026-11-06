use super::*;
use crate::SpeciesId;
use crate::solution::{BinaryCoexistenceOptionsV1, BinaryRegularSolutionV1};

fn model(n: usize, pairs: &[(usize, usize, f64)]) -> MulticomponentRegularSolutionV1 {
    let mut chi = vec![0.0; n * n];
    for &(i, j, value) in pairs {
        chi[i * n + j] = value;
        chi[j * n + i] = value;
    }
    MulticomponentRegularSolutionV1::new(
        (0..n)
            .map(|i| SpeciesId::new(format!("S{i:03}")).unwrap())
            .collect(),
        chi,
        vec![0.0; n * n],
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap()
}

fn symmetric_ternary(chi: f64) -> MulticomponentRegularSolutionV1 {
    model(3, &[(0, 1, chi), (0, 2, chi), (1, 2, chi)])
}

fn ternary_seed(rich: f64, weights: &[f64]) -> MulticomponentEquilibriumSeedV1 {
    let dilute = (1.0 - rich) / 2.0;
    MulticomponentEquilibriumSeedV1 {
        phase_mole_fractions: vec![
            vec![rich, dilute, dilute],
            vec![dilute, rich, dilute],
            vec![dilute, dilute, rich],
        ],
        phase_weights: weights.to_vec(),
    }
}

fn binary_seed() -> MulticomponentEquilibriumSeedV1 {
    MulticomponentEquilibriumSeedV1 {
        phase_mole_fractions: vec![vec![0.9, 0.1], vec![0.1, 0.9]],
        phase_weights: vec![1.0, 1.0],
    }
}

fn solve(
    law: &MulticomponentRegularSolutionV1,
    amounts: &[f64],
    seed: &MulticomponentEquilibriumSeedV1,
) -> MulticomponentEquilibriumCandidateV1 {
    law.equilibrium_candidate(
        Temperature::new(400.0),
        amounts,
        seed,
        MulticomponentEquilibriumControlV1::default(),
        || false,
    )
    .unwrap()
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "actual={a:e}, expected={b:e}, tolerance={tolerance:e}"
    );
}

// Independent scalar common-tangent oracle. The outer, stable root is above
// 2/3 for both tested chi values; the central 1/3 stationary root is excluded.
fn symmetric_three_phase_composition(chi: f64) -> (f64, f64) {
    let (mut low, mut high): (f64, f64) = (2.0 / 3.0, 1.0 - 1e-14);
    for _ in 0..90 {
        let rich = 0.5 * (low + high);
        let dilute = (1.0 - rich) / 2.0;
        if (rich / dilute).ln() - chi * (rich - dilute) > 0.0 {
            high = rich;
        } else {
            low = rich;
        }
    }
    let rich = 0.5 * (low + high);
    (rich, (1.0 - rich) / 2.0)
}

// Explicit free energy and chemical potentials, independent of the production
// state evaluator, the Newton gradient and the Newton Hessian.
fn oracle(law: &MulticomponentRegularSolutionV1, t: f64, x: &[f64]) -> (f64, Vec<Option<f64>>) {
    let n = x.len();
    let mut excess = 0.0;
    let mut rows = vec![0.0; n];
    for i in 0..n {
        for j in i + 1..n {
            let chi =
                law.chi_constants()[i * n + j] + law.chi_inverse_temperature_k()[i * n + j] / t;
            excess += chi * x[i] * x[j];
            rows[i] += chi * x[j];
            rows[j] += chi * x[i];
        }
    }
    let ideal: f64 = x.iter().filter(|v| **v > 0.0).map(|v| v * v.ln()).sum();
    (
        ideal + excess,
        x.iter()
            .enumerate()
            .map(|(i, &v)| (v > 0.0).then(|| v.ln() + rows[i] - excess))
            .collect(),
    )
}

fn audit(candidate: &MulticomponentEquilibriumCandidateV1) {
    let n = candidate.component_moles().len();
    let t = candidate.feed().temperature().value();
    let law = candidate.feed().model();
    let mut chemical = Vec::new();
    let mut energy = 0.0;
    for phase in candidate.phases() {
        assert_eq!(phase.thermodynamics().model(), law);
        assert_eq!(
            phase.thermodynamics().temperature(),
            candidate.feed().temperature()
        );
        assert!(phase.total_moles() > 0.0);
        near(
            phase.component_moles().iter().sum(),
            phase.total_moles(),
            2e-14 * phase.total_moles(),
        );
        for i in 0..n {
            if candidate.component_moles()[i] == 0.0 {
                assert_eq!(phase.component_moles()[i], 0.0);
                assert_eq!(phase.thermodynamics().mole_fractions()[i], 0.0);
                assert_eq!(
                    phase.thermodynamics().reduced_chemical_potentials()[i],
                    None
                );
            } else {
                assert!(phase.component_moles()[i] > 0.0);
                let reconstructed = phase.component_moles()[i] / phase.total_moles();
                near(
                    phase.thermodynamics().mole_fractions()[i] / reconstructed,
                    1.0,
                    2e-14,
                );
            }
        }
        assert!(matches!(
            phase.local_stability().kind(),
            MulticomponentLocalStabilityKindV1::LocallyStable
                | MulticomponentLocalStabilityKindV1::NoCompositionDirections
        ));
        let (g, mu) = oracle(law, t, phase.thermodynamics().mole_fractions());
        energy += phase.total_moles() * g * R * t;
        chemical.push(mu);
    }
    for i in 0..n {
        let recovered: f64 = candidate
            .phases()
            .iter()
            .map(|p| p.component_moles()[i])
            .sum();
        if candidate.component_moles()[i] > 0.0 {
            near(recovered / candidate.component_moles()[i], 1.0, 2e-14);
            for phase in 1..candidate.phases().len() {
                near(chemical[phase][i].unwrap(), chemical[0][i].unwrap(), 3e-9);
            }
        } else {
            assert_eq!(recovered, 0.0);
        }
    }
    near(
        energy,
        candidate.report().equilibrium_mixing_gibbs_j,
        2e-11 * candidate.total_moles(),
    );
    assert!(
        candidate
            .report()
            .maximum_reduced_chemical_potential_residual
            <= candidate.controls().chemical_potential_tolerance
    );
    assert!(
        candidate.report().maximum_relative_species_balance_residual
            <= candidate.controls().relative_species_balance_tolerance
    );
    assert!(
        candidate.report().equilibrium_mixing_gibbs_j
            <= candidate.report().feed_mixing_gibbs_j + 1e-8 * candidate.total_moles()
    );
    assert!(
        candidate.report().equilibrium_mixing_gibbs_j
            <= candidate.report().seed_mixing_gibbs_j + 1e-8 * candidate.total_moles()
    );
}

#[test]
fn g1_binary_reduces_to_independent_common_tangent_and_molar_lever_rule() {
    let species = model(2, &[]).species().to_vec();
    let law = MulticomponentRegularSolutionV1::new(
        species.clone(),
        vec![0.0, 1.0, 1.0, 0.0],
        vec![0.0, 800.0, 800.0, 0.0],
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap();
    let binary = BinaryRegularSolutionV1::new(
        [species[0].clone(), species[1].clone()],
        1.0,
        800.0,
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap();
    let exact = binary
        .coexistence(
            Temperature::new(400.0),
            BinaryCoexistenceOptionsV1::default(),
            || false,
        )
        .unwrap()
        .unwrap();
    let high = exact.upper_a().mole_fractions()[0];
    let low = exact.lower_a().mole_fractions()[0];
    let candidate = solve(&law, &[3.5, 6.5], &binary_seed());
    assert_eq!(candidate.component_moles(), &[3.5, 6.5]);
    assert_eq!(candidate.total_moles(), 10.0);
    assert_eq!(candidate.phases().len(), 2);
    near(
        candidate.phases()[0].thermodynamics().mole_fractions()[0],
        high,
        2e-10,
    );
    near(
        candidate.phases()[1].thermodynamics().mole_fractions()[0],
        low,
        2e-10,
    );
    near(
        candidate.phases()[0].mole_fraction_of_feed(),
        (0.35 - low) / (high - low),
        2e-10,
    );
    assert!(candidate.report().completed_iterations > 0);
    audit(&candidate);
}

#[test]
fn g1_three_phase_ternary_matches_outer_root_and_asymmetric_feed_amounts() {
    let law = symmetric_ternary(4.0);
    let seed = ternary_seed(0.9, &[2.0, 3.0, 5.0]);
    let candidate = solve(&law, &[2.0, 3.0, 5.0], &seed);
    let (rich, dilute) = symmetric_three_phase_composition(4.0);
    near(rich, 0.9539302688773259, 2e-15);
    assert_eq!(candidate.phases().len(), 3);
    for (phase, &feed) in [0.2, 0.3, 0.5].iter().enumerate() {
        for species in 0..3 {
            near(
                candidate.phases()[phase].thermodynamics().mole_fractions()[species],
                if phase == species { rich } else { dilute },
                3e-10,
            );
        }
        near(
            candidate.phases()[phase].mole_fraction_of_feed(),
            (feed - dilute) / (rich - dilute),
            3e-10,
        );
    }
    audit(&candidate);
    assert_eq!(candidate, solve(&law, &[2.0, 3.0, 5.0], &seed));
}

#[test]
fn g1_metastable_equimolar_feed_reaches_lower_energy_three_phase_coexistence() {
    let law = symmetric_ternary(2.9);
    let candidate = solve(&law, &[1.0; 3], &ternary_seed(0.8, &[1.0; 3]));
    let assessment = candidate
        .feed()
        .local_stability(MulticomponentStabilityControlV1::default(), || false)
        .unwrap();
    assert_eq!(
        assessment.kind(),
        MulticomponentLocalStabilityKindV1::LocallyStable
    );
    let (rich, dilute) = symmetric_three_phase_composition(2.9);
    near(rich, 0.7684741317715524, 4e-15);
    for (p, phase) in candidate.phases().iter().enumerate() {
        near(phase.mole_fraction_of_feed(), 1.0 / 3.0, 2e-10);
        for i in 0..3 {
            near(
                phase.thermodynamics().mole_fractions()[i],
                if i == p { rich } else { dilute },
                3e-10,
            );
        }
        assert!(phase.local_stability().smallest_eigenvalue().unwrap() > 0.0);
    }
    let delta = (candidate.report().equilibrium_mixing_gibbs_j
        - candidate.report().feed_mixing_gibbs_j)
        / (3.0 * R * 400.0);
    near(delta, -0.014813560085012667, 3e-13);
    audit(&candidate);
}

#[test]
fn g1_metastable_three_phase_tie_triangle_gives_nontrivial_lever_amounts() {
    let law = symmetric_ternary(2.9);
    let candidate = solve(&law, &[0.5, 0.3, 0.2], &ternary_seed(0.8, &[0.5, 0.3, 0.2]));
    let expected = [0.5886785262224036, 0.2822642947555193, 0.1290571790220771];
    for (phase, &fraction) in candidate.phases().iter().zip(&expected) {
        near(phase.mole_fraction_of_feed(), fraction, 4e-10);
    }
    audit(&candidate);
}

#[test]
fn g1_asymmetric_ternary_couples_solvent_partition_to_both_phase_compositions() {
    let law = model(3, &[(0, 1, 4.0), (0, 2, 0.5), (1, 2, 1.2)]);
    let seed = MulticomponentEquilibriumSeedV1 {
        phase_mole_fractions: vec![vec![0.75, 0.05, 0.2], vec![0.05, 0.75, 0.2]],
        phase_weights: vec![0.55, 0.45],
    };
    let candidate = solve(&law, &[0.45, 0.35, 0.2], &seed);
    let a = candidate.phases()[0].thermodynamics().mole_fractions();
    let b = candidate.phases()[1].thermodynamics().mole_fractions();
    assert!(a[0] > 0.7 && b[1] > 0.7, "{a:?}, {b:?}");
    assert!(
        a[2] > b[2],
        "solvent favors its less repulsive pair: {a:?}, {b:?}"
    );
    assert!((a[0] - b[1]).abs() > 0.01);
    assert!(
        candidate.report().equilibrium_mixing_gibbs_j
            < candidate.report().feed_mixing_gibbs_j - 10.0
    );
    audit(&candidate);
}

#[test]
fn g0_exact_absences_stay_absent_in_each_of_three_phases() {
    let law = model(
        4,
        &[
            (0, 1, 4.0),
            (0, 2, 4.0),
            (1, 2, 4.0),
            (0, 3, 9.0),
            (1, 3, -2.0),
        ],
    );
    let mut seed = ternary_seed(0.9, &[2.0, 3.0, 5.0]);
    for x in &mut seed.phase_mole_fractions {
        x.push(0.0);
    }
    let candidate = solve(&law, &[2.0, 3.0, 5.0, 0.0], &seed);
    assert_eq!(candidate.phases().len(), 3);
    audit(&candidate);
}

#[test]
fn g0_finite_trace_partitions_without_a_floor_or_absolute_inventory_allowance() {
    let law = model(3, &[(0, 1, 4.0), (0, 2, 0.6), (1, 2, 1.1)]);
    let seed = MulticomponentEquilibriumSeedV1 {
        phase_mole_fractions: vec![vec![0.95, 0.05, 1e-120], vec![0.05, 0.95, 1e-120]],
        phase_weights: vec![1.0, 1.0],
    };
    let candidate = solve(&law, &[0.5, 0.5, 1e-120], &seed);
    let a = candidate.phases()[0].thermodynamics().mole_fractions();
    let b = candidate.phases()[1].thermodynamics().mole_fractions();
    assert!(a[2] > b[2] && b[2] > 0.0);
    assert!(a[2] < 1e-119);
    audit(&candidate);
}

#[test]
fn g0_pure_feed_preserves_moles_with_no_nonlinear_iteration() {
    let law = symmetric_ternary(4.0);
    let candidate = law
        .equilibrium_candidate(
            Temperature::new(400.0),
            &[0.0, 4.5, 0.0],
            &MulticomponentEquilibriumSeedV1 {
                phase_mole_fractions: vec![vec![0.0, 1.0, 0.0]],
                phase_weights: vec![2.0],
            },
            MulticomponentEquilibriumControlV1 {
                max_iterations: 0,
                ..Default::default()
            },
            || false,
        )
        .unwrap();
    assert_eq!(candidate.phases()[0].component_moles(), &[0.0, 4.5, 0.0]);
    assert_eq!(candidate.report().completed_iterations, 0);
    assert_eq!(candidate.report().equilibrium_mixing_gibbs_j, 0.0);
    audit(&candidate);
}

#[test]
fn g0_stationary_spinodal_and_coincident_phases_are_not_accepted() {
    let seed = MulticomponentEquilibriumSeedV1 {
        phase_mole_fractions: vec![vec![0.5, 0.5]; 2],
        phase_weights: vec![1.0; 2],
    };
    let unstable = model(2, &[(0, 1, 4.0)]).equilibrium_candidate(
        Temperature::new(400.0),
        &[1.0, 1.0],
        &seed,
        Default::default(),
        || false,
    );
    assert!(matches!(
        unstable,
        Err(MulticomponentEquilibriumErrorV1::LocallyUnstablePhase { .. })
    ));
    let coincident = model(2, &[]).equilibrium_candidate(
        Temperature::new(400.0),
        &[1.0, 1.0],
        &seed,
        Default::default(),
        || false,
    );
    assert!(matches!(
        coincident,
        Err(MulticomponentEquilibriumErrorV1::CoincidentPhases { .. })
    ));
}

#[test]
fn g0_work_scratch_phase_count_and_disappearance_limits_refuse_atomically() {
    let law = model(2, &[(0, 1, 3.0)]);
    let feed = [0.35, 0.65];
    let seed = binary_seed();
    let run = |controls| {
        law.equilibrium_candidate(Temperature::new(400.0), &feed, &seed, controls, || false)
    };
    assert!(matches!(
        run(MulticomponentEquilibriumControlV1 {
            max_workspace_bytes: 0,
            ..Default::default()
        }),
        Err(MulticomponentEquilibriumErrorV1::WorkspaceBudget { .. })
    ));
    assert!(matches!(
        run(MulticomponentEquilibriumControlV1 {
            max_iterations: 0,
            ..Default::default()
        }),
        Err(MulticomponentEquilibriumErrorV1::IterationLimit {
            completed_iterations: 0
        })
    ));
    assert!(matches!(
        run(MulticomponentEquilibriumControlV1 {
            max_phases: 1,
            ..Default::default()
        }),
        Err(MulticomponentEquilibriumErrorV1::Invalid("phase count"))
    ));
    let disappearance = law.equilibrium_candidate(
        Temperature::new(400.0),
        &feed,
        &MulticomponentEquilibriumSeedV1 {
            phase_weights: vec![100.0, 1.0],
            ..seed.clone()
        },
        MulticomponentEquilibriumControlV1 {
            minimum_phase_fraction: 0.2,
            ..Default::default()
        },
        || false,
    );
    assert!(matches!(
        disappearance,
        Err(MulticomponentEquilibriumErrorV1::PhaseDisappeared { .. })
    ));
    assert_eq!(feed, [0.35, 0.65]);
    assert_eq!(seed, binary_seed());
}

#[test]
fn g0_seed_support_and_positive_representable_inventory_are_required() {
    let law = model(2, &[(0, 1, 3.0)]);
    let endpoint_seed = MulticomponentEquilibriumSeedV1 {
        phase_mole_fractions: vec![vec![1.0, 0.0], vec![0.0, 1.0]],
        phase_weights: vec![1.0; 2],
    };
    assert!(matches!(
        law.equilibrium_candidate(
            Temperature::new(400.0),
            &[1.0; 2],
            &endpoint_seed,
            Default::default(),
            || false,
        ),
        Err(MulticomponentEquilibriumErrorV1::Invalid(_))
    ));
    assert!(matches!(
        law.equilibrium_candidate(
            Temperature::new(400.0),
            &[1e300, 1e-300],
            &binary_seed(),
            Default::default(),
            || false,
        ),
        Err(MulticomponentEquilibriumErrorV1::Unrepresentable(_))
    ));
    let excessive_count = MulticomponentEquilibriumSeedV1 {
        phase_mole_fractions: vec![vec![0.5, 0.5]; 3],
        phase_weights: vec![1.0; 3],
    };
    assert!(matches!(
        law.equilibrium_candidate(
            Temperature::new(400.0),
            &[1.0; 2],
            &excessive_count,
            Default::default(),
            || false,
        ),
        Err(MulticomponentEquilibriumErrorV1::Invalid("phase count"))
    ));
}

#[test]
fn g0_cancellation_before_during_and_after_optimization_publishes_no_candidate() {
    let law = symmetric_ternary(4.0);
    let seed = ternary_seed(0.9, &[2.0, 3.0, 5.0]);
    let mut polls = 0;
    law.equilibrium_candidate(
        Temperature::new(400.0),
        &[2.0, 3.0, 5.0],
        &seed,
        Default::default(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    assert!(polls > 100);
    for allowed_polls in [0, 12, polls / 2, polls - 1] {
        let mut count = 0;
        let result = law.equilibrium_candidate(
            Temperature::new(400.0),
            &[2.0, 3.0, 5.0],
            &seed,
            Default::default(),
            || {
                count += 1;
                count > allowed_polls
            },
        );
        assert_eq!(result, Err(MulticomponentEquilibriumErrorV1::Cancelled));
    }
    assert_eq!(seed, ternary_seed(0.9, &[2.0, 3.0, 5.0]));
}
