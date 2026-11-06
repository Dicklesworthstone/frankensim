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

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "actual={a:e}, expected={b:e}, tolerance={tolerance:e}"
    );
}

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

#[test]
fn g1_binary_spinodal_discovers_the_independent_two_phase_common_tangent() {
    let species = model(2, &[]).species().to_vec();
    let law = MulticomponentRegularSolutionV1::new(
        species.clone(),
        vec![0.0, 1.0, 1.0, 0.0],
        vec![0.0, 800.0, 800.0, 0.0],
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap();
    let oracle = BinaryRegularSolutionV1::new(
        [species[0].clone(), species[1].clone()],
        1.0,
        800.0,
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap()
    .coexistence(
        Temperature::new(400.0),
        BinaryCoexistenceOptionsV1::default(),
        || false,
    )
    .unwrap()
    .unwrap();

    let discovered = law
        .discover_equilibrium_candidate(
            Temperature::new(400.0),
            &[3.5, 6.5],
            Default::default(),
            || false,
        )
        .unwrap();
    assert_eq!(discovered.candidate().phases().len(), 2);
    assert_eq!(
        discovered.report().feed_local_stability().kind(),
        MulticomponentLocalStabilityKindV1::Spinodal
    );
    let mut resolved: Vec<_> = discovered
        .candidate()
        .phases()
        .iter()
        .map(|phase| phase.thermodynamics().mole_fractions()[0])
        .collect();
    resolved.sort_by(f64::total_cmp);
    near(resolved[0], oracle.lower_a().mole_fractions()[0], 3e-10);
    near(resolved[1], oracle.upper_a().mole_fractions()[0], 3e-10);
    assert!(discovered.report().selected_reduced_gibbs_improvement() > 0.0);
    assert!(discovered.report().attempted_candidates() >= 2);
    assert!(discovered.report().converged_candidates() >= 1);
}

#[test]
fn g1_locally_stable_metastable_ternary_discovers_three_phase_lower_gibbs_state() {
    let law = symmetric_ternary(2.9);
    let discovered = law
        .discover_equilibrium_candidate(
            Temperature::new(400.0),
            &[1.0, 1.0, 1.0],
            Default::default(),
            || false,
        )
        .unwrap();
    assert_eq!(
        discovered.report().feed_local_stability().kind(),
        MulticomponentLocalStabilityKindV1::LocallyStable
    );
    assert_eq!(discovered.candidate().phases().len(), 3);
    assert_eq!(discovered.report().selected_phase_count(), 3);
    assert!(
        discovered
            .report()
            .most_negative_seed_tangent_plane_distance()
            .is_some_and(|value| value < 0.0)
    );
    let (rich, dilute) = symmetric_three_phase_composition(2.9);
    let mut rich_values: Vec<_> = discovered
        .candidate()
        .phases()
        .iter()
        .map(|phase| {
            phase
                .thermodynamics()
                .mole_fractions()
                .iter()
                .copied()
                .fold(0.0, f64::max)
        })
        .collect();
    rich_values.sort_by(f64::total_cmp);
    for value in rich_values {
        near(value, rich, 4e-10);
    }
    for phase in discovered.candidate().phases() {
        let minimum = phase
            .thermodynamics()
            .mole_fractions()
            .iter()
            .copied()
            .fold(1.0, f64::min);
        near(minimum, dilute, 4e-10);
        near(phase.mole_fraction_of_feed(), 1.0 / 3.0, 4e-10);
    }
    near(
        discovered.report().selected_reduced_gibbs_improvement(),
        0.014813560085012667,
        4e-13,
    );
}

#[test]
fn g1_asymmetric_ternary_discovers_a_nontrivial_two_phase_partition() {
    let law = model(3, &[(0, 1, 4.0), (0, 2, 0.5), (1, 2, 1.2)]);
    let discovered = law
        .discover_equilibrium_candidate(
            Temperature::new(400.0),
            &[0.45, 0.35, 0.20],
            Default::default(),
            || false,
        )
        .unwrap();
    assert_eq!(discovered.candidate().phases().len(), 2);
    let first = discovered.candidate().phases()[0]
        .thermodynamics()
        .mole_fractions();
    let second = discovered.candidate().phases()[1]
        .thermodynamics()
        .mole_fractions();
    assert!(
        (first[0] > first[1] && second[1] > second[0])
            || (second[0] > second[1] && first[1] > first[0])
    );
    assert!(discovered.report().selected_reduced_gibbs_improvement() > 0.0);
    for i in 0..3 {
        let recovered: f64 = discovered
            .candidate()
            .phases()
            .iter()
            .map(|phase| phase.component_moles()[i])
            .sum();
        near(recovered, [0.45, 0.35, 0.20][i], 2e-14);
    }
}

#[test]
fn g0_pure_and_ideal_feeds_select_one_phase_without_inventing_separation() {
    let pure = symmetric_ternary(4.0)
        .discover_equilibrium_candidate(
            Temperature::new(400.0),
            &[0.0, 4.5, 0.0],
            Default::default(),
            || false,
        )
        .unwrap();
    assert_eq!(pure.candidate().phases().len(), 1);
    assert_eq!(pure.candidate().phases()[0].component_moles(), [0.0, 4.5, 0.0]);

    let ideal = symmetric_ternary(0.0)
        .discover_equilibrium_candidate(
            Temperature::new(400.0),
            &[2.0, 3.0, 5.0],
            Default::default(),
            || false,
        )
        .unwrap();
    assert_eq!(ideal.candidate().phases().len(), 1);
    assert_eq!(
        ideal.report().feed_local_stability().kind(),
        MulticomponentLocalStabilityKindV1::LocallyStable
    );
    near(ideal.report().selected_reduced_gibbs_improvement(), 0.0, 2e-14);
}

#[test]
fn g3_search_is_deterministic_and_budget_exhaustion_is_explicit() {
    let law = symmetric_ternary(2.9);
    let a = law
        .discover_equilibrium_candidate(
            Temperature::new(400.0),
            &[0.5, 0.3, 0.2],
            Default::default(),
            || false,
        )
        .unwrap();
    let b = law
        .discover_equilibrium_candidate(
            Temperature::new(400.0),
            &[0.5, 0.3, 0.2],
            Default::default(),
            || false,
        )
        .unwrap();
    assert_eq!(a, b);

    let limited = law
        .discover_equilibrium_candidate(
            Temperature::new(400.0),
            &[0.5, 0.3, 0.2],
            MulticomponentPhaseDiscoveryControlV1 {
                maximum_seed_trials: 1,
                ..Default::default()
            },
            || false,
        )
        .unwrap();
    assert_eq!(limited.candidate().phases().len(), 1);
    assert!(limited.report().search_exhausted());
    assert_eq!(limited.report().attempted_candidates(), 1);
}

#[test]
fn g4_cancellation_at_every_observed_boundary_publishes_no_partial_discovery() {
    let law = symmetric_ternary(2.9);
    let controls = MulticomponentPhaseDiscoveryControlV1 {
        maximum_seed_trials: 12,
        ..Default::default()
    };
    let mut polls = 0;
    let complete = law
        .discover_equilibrium_candidate(
            Temperature::new(400.0),
            &[1.0, 1.0, 1.0],
            controls,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(polls > 20);
    for stop in [0, 1, polls / 3, polls / 2, polls - 1] {
        let mut current = 0;
        let result = law.discover_equilibrium_candidate(
            Temperature::new(400.0),
            &[1.0, 1.0, 1.0],
            controls,
            || {
                let cancel = current == stop;
                current += 1;
                cancel
            },
        );
        assert_eq!(result, Err(MulticomponentPhaseDiscoveryErrorV1::Cancelled));
    }
    assert_eq!(
        law.discover_equilibrium_candidate(
            Temperature::new(400.0),
            &[1.0, 1.0, 1.0],
            controls,
            || false
        )
        .unwrap(),
        complete
    );
}
