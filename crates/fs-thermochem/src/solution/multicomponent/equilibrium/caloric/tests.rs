use super::*;
use crate::solution::caloric::{
    BinaryCaloricHeatCapacityQualificationV1, BinaryLiquidCaloricModelV1,
};
use crate::solution::multicomponent::MulticomponentRegularSolutionV1;
use crate::solution::multicomponent::equilibrium::MulticomponentEquilibriumSeedV1;
use crate::solution::{BinaryCoexistenceOptionsV1, BinaryRegularSolutionV1};
use crate::{MolarThermalQuantityV1, SpeciesId};

fn law(n: usize, a: f64, b: f64) -> MulticomponentRegularSolutionV1 {
    let matrix = |value| {
        (0..n)
            .flat_map(|i| (0..n).map(move |j| if i == j { 0.0 } else { value }))
            .collect()
    };
    MulticomponentRegularSolutionV1::new(
        (0..n)
            .map(|i| SpeciesId::new(format!("S{i:03}")).unwrap())
            .collect(),
        matrix(a),
        matrix(b),
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap()
}

fn curves(law: &MulticomponentRegularSolutionV1) -> Vec<PureLiquidHeatCapacityV1> {
    law.species()
        .iter()
        .enumerate()
        .map(|(i, species)| {
            let scale = (i + 1) as f64;
            PureLiquidHeatCapacityV1::new(
                species.clone(),
                Temperature::new(300.0),
                [
                    (100.0, 50.0 * scale),
                    (500.0, 90.0 * scale),
                    (1000.0, 190.0 * scale),
                ]
                .into_iter()
                .map(|(t, cp)| (Temperature::new(t), MolarThermalQuantityV1::new(cp)))
                .collect(),
            )
            .unwrap()
        })
        .collect()
}

fn seed(n: usize, rich: f64) -> MulticomponentEquilibriumSeedV1 {
    let dilute = (1.0 - rich) / (n - 1) as f64;
    MulticomponentEquilibriumSeedV1 {
        phase_mole_fractions: (0..n)
            .map(|i| (0..n).map(|j| if i == j { rich } else { dilute }).collect())
            .collect(),
        phase_weights: vec![1.0; n],
    }
}

fn candidate(
    law: &MulticomponentRegularSolutionV1,
    t: f64,
    amounts: &[f64],
    seed: &MulticomponentEquilibriumSeedV1,
) -> MulticomponentEquilibriumCandidateV1 {
    law.equilibrium_candidate(
        Temperature::new(t),
        amounts,
        seed,
        MulticomponentEquilibriumControlV1 {
            chemical_potential_tolerance: 2e-12,
            ..Default::default()
        },
        || false,
    )
    .unwrap()
}

fn response<'a>(
    candidate: &'a MulticomponentEquilibriumCandidateV1,
    curves: &'a [PureLiquidHeatCapacityV1],
) -> MulticomponentCaloricResponseV1<'a> {
    candidate
        .isobaric_caloric_response(
            curves,
            MulticomponentCaloricResponseControlV1::default(),
            || false,
        )
        .unwrap()
}

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn total_enthalpy(
    candidate: &MulticomponentEquilibriumCandidateV1,
    curves: &[PureLiquidHeatCapacityV1],
) -> f64 {
    let t = candidate.feed().temperature();
    candidate
        .component_moles()
        .iter()
        .zip(curves)
        .map(|(n, cp)| {
            n * cp
                .enthalpy_change_j_per_mol(cp.reference_temperature(), t)
                .unwrap()
        })
        .sum::<f64>()
        + candidate
            .phases()
            .iter()
            .map(|phase| phase.total_moles() * phase.thermodynamics().mixing_enthalpy_j_per_mol())
            .sum::<f64>()
}

fn total_entropy(
    candidate: &MulticomponentEquilibriumCandidateV1,
    curves: &[PureLiquidHeatCapacityV1],
) -> f64 {
    let t = candidate.feed().temperature();
    candidate
        .component_moles()
        .iter()
        .zip(curves)
        .map(|(n, cp)| {
            n * cp
                .entropy_change_j_per_mol_k(cp.reference_temperature(), t)
                .unwrap()
        })
        .sum::<f64>()
        + candidate
            .phases()
            .iter()
            .map(|phase| phase.total_moles() * phase.thermodynamics().mixing_entropy_j_per_mol_k())
            .sum::<f64>()
}

#[test]
fn g1_binary_branch_response_matches_original_closed_form_caloric_owner() {
    for (a, b, t) in [(1.0, 800.0, 400.0), (5.0, -900.0, 450.0)] {
        let law = law(2, a, b);
        let cp = curves(&law);
        let split = candidate(&law, t, &[3.5, 6.5], &seed(2, 0.95));
        let result = response(&split, &cp);
        assert!(core::ptr::eq(result.candidate(), &split));
        assert!(core::ptr::eq(
            result.pure_liquid_heat_capacities(),
            cp.as_slice()
        ));
        assert_eq!(result.temperature_domain(), law.temperature_domain());
        let binary = BinaryLiquidCaloricModelV1::new(
            BinaryRegularSolutionV1::new(
                [law.species()[0].clone(), law.species()[1].clone()],
                a,
                b,
                Temperature::new(100.0),
                Temperature::new(1000.0),
            )
            .unwrap(),
            [cp[0].clone(), cp[1].clone()],
        )
        .unwrap();
        let exact = binary
            .equilibrate(
                Temperature::new(t),
                [3.5, 6.5],
                BinaryCoexistenceOptionsV1::default(),
                || false,
            )
            .unwrap();
        assert_eq!(
            exact
                .isobaric_heat_capacity_response()
                .unwrap()
                .qualification(),
            BinaryCaloricHeatCapacityQualificationV1::Coexisting
        );
        near(
            result.isobaric_heat_capacity_j_per_k(),
            exact.isobaric_heat_capacity_j_per_k().unwrap(),
            2e-8,
            3e-11,
        );
        assert!(result.redistribution_heat_capacity_j_per_k() > 0.0);
        assert!(result.maximum_scaled_residual() < 1e-14);
        assert!(
            result
                .phase_component_moles_temperature_derivatives_mol_per_k(2)
                .is_none()
        );
    }
}

#[test]
fn g1_ternary_response_matches_scalar_common_tangent_and_original_h_s_derivatives() {
    let law = law(3, 0.5, 1000.0);
    let cp = curves(&law);
    let t = 400.0;
    let amounts = [2.0, 3.0, 5.0];
    let guess = seed(3, 0.9);
    let center = candidate(&law, t, &amounts, &guess);
    let result = response(&center, &cp);
    let chi = 0.5 + 1000.0 / t;
    // Independent scalar three-phase common tangent: each phase is a
    // permutation of (r,d,d), d=(1-r)/2, ln(r/d)=chi*(r-d).
    let (mut lo, mut hi) = (2.0 / 3.0, 1.0 - 1e-12);
    for _ in 0..90 {
        let rich: f64 = 0.5 * (lo + hi);
        if (2.0 * rich / (1.0 - rich)).ln() > chi * (1.5 * rich - 0.5) {
            hi = rich;
        } else {
            lo = rich;
        }
    }
    let rich = 0.5 * (lo + hi);
    let curvature = 1.0 / rich + 1.0 / (1.0 - rich) - 1.5 * chi;
    let oracle = 10.0 * R * (1000.0 / t).powi(2) * (3.0 * rich - 1.0).powi(2) / (4.0 * curvature);
    near(
        result.redistribution_heat_capacity_j_per_k(),
        oracle,
        3e-8,
        3e-11,
    );
    let dt = 0.002;
    let cold = candidate(&law, t - dt, &amounts, &guess);
    let hot = candidate(&law, t + dt, &amounts, &guess);
    near(
        result.isobaric_heat_capacity_j_per_k(),
        (total_enthalpy(&hot, &cp) - total_enthalpy(&cold, &cp)) / (2.0 * dt),
        2e-5,
        0.0,
    );
    near(
        result.isobaric_heat_capacity_j_per_k(),
        t * (total_entropy(&hot, &cp) - total_entropy(&cold, &cp)) / (2.0 * dt),
        2e-5,
        0.0,
    );
    for phase in 0..3 {
        let rates = result
            .phase_component_moles_temperature_derivatives_mol_per_k(phase)
            .unwrap();
        for (i, &actual) in rates.iter().enumerate() {
            let expected = (hot.phases()[phase].component_moles()[i]
                - cold.phases()[phase].component_moles()[i])
                / (2.0 * dt);
            near(actual, expected, 5e-10, 2e-6);
        }
        near(
            result.phase_total_moles_temperature_derivatives_mol_per_k()[phase],
            rates.iter().sum(),
            1e-17,
            1e-14,
        );
    }
    for species in 0..3 {
        let sum: f64 = (0..3)
            .map(|phase| {
                result
                    .phase_component_moles_temperature_derivatives_mol_per_k(phase)
                    .unwrap()[species]
            })
            .sum();
        near(sum, 0.0, 2e-18, 0.0);
    }
}

#[test]
fn g1_response_preserves_inventory_scaling_and_species_permutation() {
    let original = law(3, 0.5, 1000.0);
    let cp = curves(&original);
    let guess = seed(3, 0.9);
    let amounts = [2.0, 3.0, 5.0];
    let center = candidate(&original, 400.0, &amounts, &guess);
    let actual = response(&center, &cp);
    for scale in [1e-80, 1e80] {
        let scaled = candidate(&original, 400.0, &amounts.map(|n| n * scale), &guess);
        let scaled_response = response(&scaled, &cp);
        near(
            scaled_response.isobaric_heat_capacity_j_per_k() / scale,
            actual.isobaric_heat_capacity_j_per_k(),
            2e-9,
            1e-12,
        );
        for phase in 0..3 {
            let rates = scaled_response
                .phase_component_moles_temperature_derivatives_mol_per_k(phase)
                .unwrap();
            for (i, &value) in rates.iter().enumerate() {
                near(
                    value / scale,
                    actual
                        .phase_component_moles_temperature_derivatives_mol_per_k(phase)
                        .unwrap()[i],
                    3e-14,
                    1e-10,
                );
            }
        }
    }
    let permuted = candidate(&original, 400.0, &[5.0, 2.0, 3.0], &guess);
    let permuted_cp = [
        PureLiquidHeatCapacityV1::new(
            original.species()[0].clone(),
            Temperature::new(300.0),
            vec![
                (Temperature::new(100.0), MolarThermalQuantityV1::new(150.0)),
                (Temperature::new(500.0), MolarThermalQuantityV1::new(270.0)),
                (Temperature::new(1000.0), MolarThermalQuantityV1::new(570.0)),
            ],
        )
        .unwrap(),
        PureLiquidHeatCapacityV1::new(
            original.species()[1].clone(),
            Temperature::new(300.0),
            vec![
                (Temperature::new(100.0), MolarThermalQuantityV1::new(50.0)),
                (Temperature::new(500.0), MolarThermalQuantityV1::new(90.0)),
                (Temperature::new(1000.0), MolarThermalQuantityV1::new(190.0)),
            ],
        )
        .unwrap(),
        PureLiquidHeatCapacityV1::new(
            original.species()[2].clone(),
            Temperature::new(300.0),
            vec![
                (Temperature::new(100.0), MolarThermalQuantityV1::new(100.0)),
                (Temperature::new(500.0), MolarThermalQuantityV1::new(180.0)),
                (Temperature::new(1000.0), MolarThermalQuantityV1::new(380.0)),
            ],
        )
        .unwrap(),
    ];
    near(
        response(&permuted, &permuted_cp).isobaric_heat_capacity_j_per_k(),
        actual.isobaric_heat_capacity_j_per_k(),
        3e-9,
        1e-12,
    );
}

#[test]
fn g0_single_phase_and_temperature_independent_splits_have_zero_redistribution() {
    let ideal = law(3, 0.0, 0.0);
    let cp = curves(&ideal);
    let one_phase = MulticomponentEquilibriumSeedV1 {
        phase_mole_fractions: vec![vec![0.2, 0.3, 0.5]],
        phase_weights: vec![1.0],
    };
    let homogeneous = candidate(&ideal, 400.0, &[2.0, 3.0, 5.0], &one_phase);
    let result = response(&homogeneous, &cp);
    assert_eq!(result.isobaric_heat_capacity_j_per_k(), 1840.0);
    assert_eq!(result.redistribution_heat_capacity_j_per_k(), 0.0);
    assert_eq!(
        result.phase_component_moles_temperature_derivatives_mol_per_k(0),
        Some([0.0; 3].as_slice())
    );
    assert_eq!(
        result.phase_total_moles_temperature_derivatives_mol_per_k(),
        &[0.0]
    );
    let split_law = law(3, 3.0, 0.0);
    let split = candidate(&split_law, 400.0, &[2.0, 3.0, 5.0], &seed(3, 0.9));
    let result = response(&split, &cp);
    assert_eq!(result.isobaric_heat_capacity_j_per_k(), 1840.0);
    assert_eq!(result.redistribution_heat_capacity_j_per_k(), 0.0);
    for phase in 0..3 {
        assert!(
            result
                .phase_component_moles_temperature_derivatives_mol_per_k(phase)
                .unwrap()
                .iter()
                .all(|&rate| rate == 0.0)
        );
    }
}

#[test]
fn g1_phase_derivatives_preserve_absent_species_and_a_real_trace_inventory() {
    let law = MulticomponentRegularSolutionV1::new(
        (0..3)
            .map(|i| SpeciesId::new(format!("S{i:03}")).unwrap())
            .collect(),
        vec![0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        vec![0.0, 800.0, 0.0, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap();
    let cp = curves(&law);
    for trace in [0.0, 1e-24] {
        let guess = MulticomponentEquilibriumSeedV1 {
            phase_mole_fractions: vec![
                vec![0.95, 0.05, trace / 10.0],
                vec![0.05, 0.95, trace / 10.0],
            ],
            phase_weights: vec![1.0; 2],
        };
        let split = candidate(&law, 400.0, &[3.5, 6.5, trace], &guess);
        let result = response(&split, &cp);
        let mut column = 0.0;
        for phase in 0..2 {
            let derivative = result
                .phase_component_moles_temperature_derivatives_mol_per_k(phase)
                .unwrap()[2];
            column += derivative;
            if trace == 0.0 {
                assert_eq!(derivative, 0.0);
            } else {
                assert!(derivative != 0.0);
                // This symmetric inert tracer has the same concentration in
                // both phases. Its REAL tiny inventory follows their changing
                // amounts, without a trace floor or carrier subtraction.
                near(
                    derivative / (trace / split.total_moles()),
                    result.phase_total_moles_temperature_derivatives_mol_per_k()[phase],
                    2e-14,
                    1e-10,
                );
            }
        }
        assert_eq!(column, 0.0);
    }
}

#[test]
fn g3_source_identity_datum_domain_budget_and_cancellation_refuse_atomically() {
    let law = law(2, 1.0, 800.0);
    let cp = curves(&law);
    let split = candidate(&law, 400.0, &[3.5, 6.5], &seed(2, 0.95));
    let original = split.clone();
    let control = MulticomponentCaloricResponseControlV1::default();
    assert!(
        split
            .isobaric_caloric_response(&cp[..1], control, || false)
            .is_err()
    );
    assert!(
        split
            .isobaric_caloric_response(&[cp[1].clone(), cp[0].clone()], control, || false)
            .is_err()
    );
    for (low, high, datum) in [(100.0, 1000.0, 400.0), (401.0, 1000.0, 500.0)] {
        let other = PureLiquidHeatCapacityV1::new(
            law.species()[1].clone(),
            Temperature::new(datum),
            vec![
                (Temperature::new(low), MolarThermalQuantityV1::new(100.0)),
                (Temperature::new(high), MolarThermalQuantityV1::new(100.0)),
            ],
        )
        .unwrap();
        assert!(
            split
                .isobaric_caloric_response(&[cp[0].clone(), other], control, || false)
                .is_err()
        );
    }
    assert!(matches!(
        split.isobaric_caloric_response(
            &cp,
            MulticomponentCaloricResponseControlV1 {
                max_workspace_bytes: 1,
                ..control
            },
            || false
        ),
        Err(MulticomponentCaloricResponseErrorV1::WorkspaceBudget { .. })
    ));
    assert!(matches!(
        split.isobaric_caloric_response(&cp, control, || true),
        Err(MulticomponentCaloricResponseErrorV1::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        split.isobaric_caloric_response(&cp, control, || {
            polls += 1;
            polls == 9
        }),
        Err(MulticomponentCaloricResponseErrorV1::Cancelled)
    ));
    for tolerance in [0.0, f64::NAN, 1.0] {
        assert!(matches!(
            split.isobaric_caloric_response(
                &cp,
                MulticomponentCaloricResponseControlV1 {
                    relative_residual_tolerance: tolerance,
                    ..control
                },
                || false
            ),
            Err(MulticomponentCaloricResponseErrorV1::Invalid(_))
        ));
    }
    assert_eq!(split, original);
}

#[test]
fn g3_near_critical_joint_hessian_refuses_instead_of_regularizing_cp() {
    let chi = 2.0 + 1e-8;
    let law = law(2, 0.0, 400.0 * chi);
    let cp = curves(&law);
    let binary = BinaryRegularSolutionV1::new(
        [law.species()[0].clone(), law.species()[1].clone()],
        0.0,
        400.0 * chi,
        Temperature::new(100.0),
        Temperature::new(1000.0),
    )
    .unwrap();
    let coexistence = binary
        .coexistence(
            Temperature::new(400.0),
            BinaryCoexistenceOptionsV1::default(),
            || false,
        )
        .unwrap()
        .unwrap();
    let seed = MulticomponentEquilibriumSeedV1 {
        phase_mole_fractions: vec![
            coexistence.lower_a().mole_fractions().to_vec(),
            coexistence.upper_a().mole_fractions().to_vec(),
        ],
        phase_weights: vec![1.0; 2],
    };
    let split = candidate(&law, 400.0, &[1.0, 1.0], &seed);
    assert!(split.phases().iter().all(|phase| phase.local_stability().kind() == super::super::super::stability::MulticomponentLocalStabilityKindV1::LocallyStable));
    assert!(matches!(
        split.isobaric_caloric_response(
            &cp,
            MulticomponentCaloricResponseControlV1::default(),
            || false
        ),
        Err(MulticomponentCaloricResponseErrorV1::SingularResponse)
    ));
}
