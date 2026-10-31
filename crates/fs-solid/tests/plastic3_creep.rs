//! G1/G3/G4: physical-time creep changes a loaded bar's shape and held tension.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::creep::{CreepTolerances, NortonCreep, NortonCreepParameters};
use fs_material::{IsotropicElastic, SmallStrainLaw};
use fs_solid::plastic3::creep::{CreepTetError, CreepTetTolerances};
use fs_solid::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem,
    TetAssemblyBudget,
};

const YOUNGS: f64 = 1000.0;
const POISSON: f64 = 0.25;
const STRESS: f64 = 2.0;
const RATE: f64 = 1e-4;
const LENGTH: f64 = 2.0;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 37,
                kernel_id: 79,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn mesh() -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let nodes = vec![
        [0.0, 0.0, 0.0],
        [LENGTH, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [LENGTH, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [LENGTH, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [LENGTH, 1.0, 1.0],
    ];
    let tets = vec![
        [0, 1, 3, 7],
        [0, 1, 7, 5],
        [0, 2, 7, 3],
        [0, 2, 6, 7],
        [0, 4, 5, 7],
        [0, 4, 7, 6],
    ];
    (nodes, tets)
}

fn law(exponent: f64, temperature: f64, dt: f64) -> NortonCreep {
    NortonCreep::new(
        IsotropicElastic::new(YOUNGS, POISSON, 0.05).unwrap(),
        NortonCreepParameters {
            reference_stress_pa: STRESS,
            reference_rate_s: RATE,
            stress_exponent: exponent,
            activation_energy_j_mol: 50_000.0,
            reference_temperature_k: 300.0,
            temperature_range_k: [290.0, 350.0],
        },
        temperature,
        dt,
    )
    .unwrap()
}

fn tolerances() -> CreepTetTolerances {
    CreepTetTolerances {
        material: CreepTolerances {
            stress_pa: 1e-9,
            energy_j_m3: 1e-10,
        },
        energy_j: 1e-9,
    }
}

fn rollers(nodes: &[[f64; 3]]) -> Vec<(usize, f64)> {
    nodes
        .iter()
        .enumerate()
        .flat_map(|(i, x)| {
            (0..3).filter_map(move |d| {
                if x[d] == 0.0 {
                    Some((3 * i + d, 0.0))
                } else {
                    None
                }
            })
        })
        .collect()
}

fn elastic_displacement(nodes: &[[f64; 3]]) -> Vec<f64> {
    nodes
        .iter()
        .flat_map(|x| {
            [
                STRESS / YOUNGS * x[0],
                -POISSON * STRESS / YOUNGS * x[1],
                -POISSON * STRESS / YOUNGS * x[2],
            ]
        })
        .collect()
}

// Constant traction on the two triangles of the x=L end face. Each triangle
// contributes area/3 to its vertices; their shared diagonal is 1--7.
fn end_force() -> Vec<(usize, f64)> {
    vec![
        (3, STRESS / 3.0),
        (9, STRESS / 6.0),
        (15, STRESS / 6.0),
        (21, STRESS / 3.0),
    ]
}

#[test]
fn sustained_load_changes_bar_shape_and_unloading_retains_permanent_creep() {
    let (nodes, tets) = mesh();
    let supports = rollers(&nodes);
    let loads = end_force();
    let material = law(3.0, 300.0, 1.0);
    let problem = SmallStrainTetProblem {
        nodes_m: &nodes,
        tetrahedra: &tets,
        laws: SmallStrainLawField::Uniform(&material),
        prescribed_m: &supports,
        nodal_forces_n: &loads,
        load_path: &[1.0],
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    };
    with_cx(|cx| {
        let mut displacement = elastic_displacement(&nodes);
        let mut states = vec![material.initial_state(); tets.len()];
        for step in 1..=5 {
            let result = problem
                .solve_creep_from(&displacement, &states, tolerances(), cx)
                .unwrap();
            let elapsed = step as f64;
            let axial = STRESS / YOUNGS + RATE * elapsed;
            let lateral = -POISSON * STRESS / YOUNGS - RATE * elapsed / 2.0;
            for (x, u) in nodes.iter().zip(&result.solution.displacement_m) {
                assert!((u[0] - axial * x[0]).abs() < 1e-11);
                assert!((u[1] - lateral * x[1]).abs() < 1e-11);
                assert!((u[2] - lateral * x[2]).abs() < 1e-11);
            }
            for (stress, state) in result
                .solution
                .element_stress
                .iter()
                .zip(&result.solution.element_states)
            {
                assert!((stress[0] - STRESS).abs() < 1e-8);
                assert!(stress[1..].iter().all(|v| v.abs() < 1e-8));
                assert!((state.equivalent_creep_strain - RATE * elapsed).abs() < 1e-11);
                assert!((state.creep_strain[0] - RATE * elapsed).abs() < 1e-11);
            }
            let heat = LENGTH * STRESS * RATE;
            assert!((result.energy.physical_dissipation_j - heat).abs() < 1e-10);
            assert!(result.energy.numerical_dissipation_j < 1e-16);
            assert!((result.energy.applied_work_j - heat).abs() < 1e-10);
            assert!(result.energy.balance_residual_j.abs() < 1e-10);
            displacement = result
                .solution
                .displacement_m
                .iter()
                .flatten()
                .copied()
                .collect();
            states = result.solution.element_states;
        }
        let unload = SmallStrainTetProblem {
            nodal_forces_n: &[],
            ..problem
        };
        let result = unload
            .solve_creep_from(&displacement, &states, tolerances(), cx)
            .unwrap();
        for (x, u) in nodes.iter().zip(&result.solution.displacement_m) {
            assert!((u[0] - 5.0 * RATE * x[0]).abs() < 1e-11);
            assert!((u[1] + 2.5 * RATE * x[1]).abs() < 1e-11);
        }
        for (before, after) in states.iter().zip(&result.solution.element_states) {
            assert!((after.equivalent_creep_strain - before.equivalent_creep_strain).abs() < 1e-12);
        }
        assert!(
            result
                .solution
                .element_stress
                .iter()
                .flatten()
                .all(|s| s.abs() < 1e-8)
        );
        assert!(result.energy.physical_dissipation_j < 1e-15);
    });
}

#[test]
fn held_length_relaxes_tension_with_first_order_physical_time_refinement() {
    let (nodes, tets) = mesh();
    let mut supports = rollers(&nodes);
    supports.extend(
        nodes
            .iter()
            .enumerate()
            .filter(|(_, x)| x[0] == LENGTH)
            .map(|(i, _)| (3 * i, LENGTH * STRESS / YOUNGS)),
    );
    let mut errors = Vec::new();
    let final_time = 10.0;
    let continuous = STRESS * (-YOUNGS * RATE / STRESS * final_time).exp();
    for steps in [10, 20, 40] {
        let dt = final_time / steps as f64;
        let material = law(1.0, 300.0, dt);
        let problem = SmallStrainTetProblem {
            nodes_m: &nodes,
            tetrahedra: &tets,
            laws: SmallStrainLawField::Uniform(&material),
            prescribed_m: &supports,
            nodal_forces_n: &[],
            load_path: &[1.0],
            budget: TetAssemblyBudget::standard(),
            settings: IncrementSettings::default(),
        };
        with_cx(|cx| {
            let mut displacement = elastic_displacement(&nodes);
            let mut states = vec![material.initial_state(); tets.len()];
            let mut final_stress = 0.0;
            for step in 1..=steps {
                let result = problem
                    .solve_creep_from(&displacement, &states, tolerances(), cx)
                    .unwrap();
                let discrete = STRESS / (1.0 + YOUNGS * RATE / STRESS * dt).powi(step);
                final_stress = result.solution.element_stress[0][0];
                assert!((final_stress - discrete).abs() < 1e-8);
                assert!(result.energy.applied_work_j.abs() < 1e-18);
                assert!(result.energy.support_work_j.abs() < 1e-16);
                assert!(result.energy.physical_dissipation_j > 0.0);
                assert!(result.energy.numerical_dissipation_j > 0.0);
                assert!(result.energy.final_stored_j < result.energy.initial_stored_j);
                let force: f64 = supports
                    .iter()
                    .zip(&result.solution.increments[0].reactions_n)
                    .filter(|((d, u), _)| d % 3 == 0 && *u > 0.0)
                    .map(|(_, r)| r)
                    .sum();
                assert!((force - final_stress).abs() < 1e-8);
                displacement = result
                    .solution
                    .displacement_m
                    .iter()
                    .flatten()
                    .copied()
                    .collect();
                states = result.solution.element_states;
            }
            errors.push((final_stress - continuous).abs());
        });
    }
    assert!(
        errors[0] / errors[1] > 1.9 && errors[1] / errors[2] > 1.9,
        "{errors:?}"
    );
}

#[test]
fn per_element_temperature_preserves_distinct_history_and_heat() {
    let (nodes, tets) = mesh();
    let displacement = elastic_displacement(&nodes);
    let supports: Vec<_> = displacement.iter().copied().enumerate().collect();
    let materials: Vec<_> = (0..tets.len())
        .map(|i| law(3.0, if i % 2 == 0 { 300.0 } else { 330.0 }, 1.0))
        .collect();
    let states = vec![materials[0].initial_state(); tets.len()];
    let original = states.clone();
    let problem = SmallStrainTetProblem {
        nodes_m: &nodes,
        tetrahedra: &tets,
        laws: SmallStrainLawField::PerElement(&materials),
        prescribed_m: &supports,
        nodal_forces_n: &[],
        load_path: &[1.0],
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    };
    with_cx(|cx| {
        let result = problem
            .solve_creep_from(&displacement, &states, tolerances(), cx)
            .unwrap();
        assert_eq!(states, original);
        assert_eq!(
            result
                .solution
                .displacement_m
                .iter()
                .flatten()
                .copied()
                .collect::<Vec<_>>(),
            displacement
        );
        assert!(
            result.solution.element_states[1].equivalent_creep_strain
                > result.solution.element_states[0].equivalent_creep_strain
        );
        assert!(
            result.element_energy[1].physical_dissipation_j_m3
                > result.element_energy[0].physical_dissipation_j_m3
        );
        for i in 0..tets.len() {
            assert_eq!(
                result.solution.element_states[i],
                result.solution.element_states[i % 2]
            );
        }
        assert_eq!(result.energy.applied_work_j, 0.0);
        assert_eq!(result.energy.support_work_j, 0.0);
        assert!(result.energy.balance_residual_j.abs() < 1e-10);
    });
}

#[test]
fn source_material_and_temperature_changes_drive_the_same_loaded_bar() {
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
        PropertyValue, Provenance, QueryPoint, UncertaintyModel,
    };
    use fs_material::creep::source::{NortonCreepSource, NortonSourceQuery, norton_requirements};
    use fs_material::state_point::MaterialPropertySelection;

    let (nodes, tets) = mesh();
    let supports = rollers(&nodes);
    let loads = end_force();
    let point = |temperature| QueryPoint::new().with("T", temperature).unwrap();
    let low = point(300.0);
    let high = point(350.0);
    let thermal_factor = (50_000.0_f64 / 8.314_462_618_153_24 * (1.0 / 300.0 - 1.0 / 330.0)).exp();
    let mut final_extensions = Vec::new();
    for rate_multiplier in [1.0, 2.0] {
        let values = [
            YOUNGS,
            POISSON,
            STRESS,
            RATE * rate_multiplier,
            3.0,
            50_000.0,
            300.0,
        ];
        let mut claims = ClaimSet::new();
        for (requirement, value) in norton_requirements().unwrap().iter().zip(values) {
            let dims = requirement.dims();
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(requirement.name(), dims),
                    value: PropertyValue::Scalar { value, dims },
                    validity: ValidityDomain::unconstrained().with("T", 300.0, 350.0),
                    interpolation: InterpolationPolicy::ConstantWithinValidity,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "manufactured source-to-solid creep fixture".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: format!("synthetic creep solid {rate_multiplier}"),
                phase: "solid".into(),
                process: "manufactured implementation reference".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let source = NortonCreepSource::new(
            &card,
            NortonSourceQuery {
                lower: &low,
                upper: &high,
                temperature_axis: "T",
                selection: MaterialPropertySelection::SingleClaimOnly,
                strain_limit: 0.05,
            },
            || false,
        )
        .unwrap();
        let mut displacement = elastic_displacement(&nodes);
        let mut states =
            vec![source.at_temperature(300.0, 1.0).unwrap().initial_state(); tets.len()];
        with_cx(|cx| {
            for (temperature, integrated_rate) in [(300.0, 1.0), (330.0, 1.0 + thermal_factor)] {
                let material = source.at_temperature(temperature, 1.0).unwrap();
                let problem = SmallStrainTetProblem {
                    nodes_m: &nodes,
                    tetrahedra: &tets,
                    laws: SmallStrainLawField::Uniform(&material),
                    prescribed_m: &supports,
                    nodal_forces_n: &loads,
                    load_path: &[1.0],
                    budget: TetAssemblyBudget::standard(),
                    settings: IncrementSettings::default(),
                };
                let next = problem
                    .solve_creep_from(&displacement, &states, tolerances(), cx)
                    .unwrap();
                let expected_creep = RATE * rate_multiplier * integrated_rate;
                assert!(
                    (next.solution.displacement_m[7][0]
                        - LENGTH * (STRESS / YOUNGS + expected_creep))
                        .abs()
                        < 1e-10
                );
                for state in &next.solution.element_states {
                    assert!((state.equivalent_creep_strain - expected_creep).abs() < 1e-10);
                }
                assert_eq!(material.material().card_identity(), card.content_hash());
                assert_eq!(material.material().state_coordinate("T"), Some(temperature));
                for property in material.material().properties() {
                    let original = source
                        .reference()
                        .property(property.requirement().name())
                        .unwrap();
                    assert_eq!(
                        property.answer().receipt.selected,
                        original.answer().receipt.selected
                    );
                }
                displacement = next
                    .solution
                    .displacement_m
                    .iter()
                    .flatten()
                    .copied()
                    .collect();
                states = next.solution.element_states;
            }
        });
        final_extensions.push(displacement[21] - LENGTH * STRESS / YOUNGS);
    }
    assert!((final_extensions[1] / final_extensions[0] - 2.0).abs() < 1e-8);
}

#[test]
fn invalid_step_or_failed_equilibrium_does_not_consume_history() {
    let (nodes, tets) = mesh();
    let supports = rollers(&nodes);
    let loads = end_force();
    let material = law(3.0, 300.0, 1.0);
    let displacement = elastic_displacement(&nodes);
    let states = vec![material.initial_state(); tets.len()];
    let original_u = displacement.clone();
    let original_states = states.clone();
    let problem = SmallStrainTetProblem {
        nodes_m: &nodes,
        tetrahedra: &tets,
        laws: SmallStrainLawField::Uniform(&material),
        prescribed_m: &supports,
        nodal_forces_n: &loads,
        load_path: &[1.0],
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    };
    with_cx(|cx| {
        let multiple = SmallStrainTetProblem {
            load_path: &[1.0, 1.0],
            ..problem.clone()
        };
        assert!(matches!(
            multiple.solve_creep_from(&displacement, &states, tolerances(), cx),
            Err(CreepTetError::Invalid(_))
        ));
        let mut materials = vec![material.clone(); tets.len()];
        materials[1] = law(3.0, 300.0, 2.0);
        let unequal = SmallStrainTetProblem {
            laws: SmallStrainLawField::PerElement(&materials),
            ..problem.clone()
        };
        assert!(matches!(
            unequal.solve_creep_from(&displacement, &states, tolerances(), cx),
            Err(CreepTetError::Invalid(_))
        ));
        let failed = SmallStrainTetProblem {
            settings: IncrementSettings {
                max_newton_iterations: 1,
                ..IncrementSettings::default()
            },
            ..problem.clone()
        };
        assert!(matches!(
            failed.solve_creep_from(&displacement, &states, tolerances(), cx),
            Err(CreepTetError::Mechanical(_))
        ));
        let mut invalid_history = states.clone();
        invalid_history[0].creep_strain[0] = f64::NAN;
        assert!(matches!(
            problem.solve_creep_from(&displacement, &invalid_history, tolerances(), cx),
            Err(CreepTetError::Material { .. })
        ));
    });
    let gate = CancelGate::new();
    gate.request();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 37,
                kernel_id: 79,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        assert!(matches!(
            problem.solve_creep_from(&displacement, &states, tolerances(), &cx),
            Err(CreepTetError::Mechanical(SmallStrainTetError::Cancelled))
        ));
    });
    assert_eq!(states, original_states);
    assert_eq!(displacement, original_u);
}
