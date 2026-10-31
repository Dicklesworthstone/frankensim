//! Native equilibrium with simultaneous permanent creep and isotropic yielding.
//! Manufactured constitutive references, not measured material qualification.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::creep::plastic::NortonJ2;
use fs_material::creep::{CreepTolerances, NortonCreep, NortonCreepParameters};
use fs_material::plastic::J2Plasticity;
use fs_material::{IsotropicElastic, SmallStrainLaw};
use fs_solid::plastic3::creep::{CreepTetError, CreepTetTolerances};
use fs_solid::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem,
    TetAssemblyBudget,
};

const E: f64 = 1000.0;
const NU: f64 = 0.25;
const STRESS: f64 = 2.0;
const YIELD: f64 = 1.4;
const HARDENING: f64 = 200.0;
const RATE: f64 = 1e-4;
const VOLUME: f64 = 1.0 / 6.0;
const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const ROLLERS: [(usize, f64); 9] = [
    (0, 0.0),
    (1, 0.0),
    (2, 0.0),
    (4, 0.0),
    (5, 0.0),
    (6, 0.0),
    (8, 0.0),
    (9, 0.0),
    (10, 0.0),
];

fn law(dt: f64) -> NortonJ2 {
    let elastic = IsotropicElastic::new(E, NU, 0.05).unwrap();
    let creep = NortonCreep::new(
        elastic.clone(),
        NortonCreepParameters {
            reference_stress_pa: STRESS,
            reference_rate_s: RATE,
            stress_exponent: 3.0,
            activation_energy_j_mol: 0.0,
            reference_temperature_k: 300.0,
            temperature_range_k: [280.0, 380.0],
        },
        300.0,
        dt,
    )
    .unwrap();
    NortonJ2::new(creep, J2Plasticity::new(elastic, YIELD, HARDENING).unwrap()).unwrap()
}

fn with_gate<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 29,
                kernel_id: 0x4352504c,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
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

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

#[test]
fn loaded_tetrahedron_yields_then_keeps_creeping_with_both_histories() {
    with_gate(&CancelGate::new(), |cx| {
        let material = law(1.0);
        let mut old_u = vec![0.0; 12];
        let mut states = vec![material.initial_state()];
        let mut accumulated_creep = 0.0;
        let plastic = (STRESS - YIELD) / HARDENING;
        for (step, stress) in [STRESS, STRESS, 0.5].into_iter().enumerate() {
            let force = [(3, stress * VOLUME)];
            let problem = SmallStrainTetProblem {
                nodes_m: &NODES,
                tetrahedra: &TETS,
                laws: SmallStrainLawField::Uniform(&material),
                prescribed_m: &ROLLERS,
                nodal_forces_n: &force,
                load_path: &[1.0],
                budget: TetAssemblyBudget::standard(),
                settings: IncrementSettings::default(),
            };
            let original_u = old_u.clone();
            let original_states = states.clone();
            let next = problem
                .solve_creep_from(&old_u, &states, tolerances(), cx)
                .unwrap();
            assert_eq!(old_u, original_u);
            assert_eq!(states, original_states);
            let delta_creep = RATE * (stress / STRESS).powi(3);
            accumulated_creep += delta_creep;
            let delta_plastic = if step == 0 { plastic } else { 0.0 };
            let state = &next.solution.element_states[0];
            close(
                state.creep.equivalent_creep_strain,
                accumulated_creep,
                1e-10,
            );
            close(state.plastic.alpha, plastic, 1e-10);
            close(state.creep.creep_strain[0], accumulated_creep, 1e-10);
            close(state.plastic.plastic_strain[0], plastic, 1e-10);
            close(
                next.solution.displacement_m[1][0],
                stress / E + accumulated_creep + plastic,
                1e-10,
            );
            close(
                next.solution.displacement_m[2][1],
                -NU * stress / E - 0.5 * (accumulated_creep + plastic),
                1e-10,
            );
            close(next.solution.element_stress[0][0], stress, 1e-8);
            assert!(
                next.solution.element_stress[0][1..]
                    .iter()
                    .all(|s| s.abs() < 1e-8)
            );
            let energy = next.element_energy[0];
            close(energy.creep_dissipation_j_m3, stress * delta_creep, 1e-9);
            close(energy.plastic_dissipation_j_m3, YIELD * delta_plastic, 1e-9);
            close(
                energy.physical_dissipation_j_m3,
                stress * delta_creep + YIELD * delta_plastic,
                1e-9,
            );
            close(
                energy.final_stored_j_m3,
                stress * stress / (2.0 * E) + 0.5 * HARDENING * plastic * plastic,
                1e-9,
            );
            close(next.energy.balance_residual_j, 0.0, 1e-9);
            close(next.energy.equilibrium_work_residual_j, 0.0, 1e-9);
            old_u = next.solution.displacement_m.as_flattened().to_vec();
            states = next.solution.element_states;
        }
    });
}

#[test]
fn held_shape_relaxes_below_hardened_yield_and_refusal_preserves_both_memories() {
    let material = law(1.0);
    let force = [(3, STRESS * VOLUME)];
    let initial = SmallStrainTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        laws: SmallStrainLawField::Uniform(&material),
        prescribed_m: &ROLLERS,
        nodal_forces_n: &force,
        load_path: &[1.0],
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    };
    let loaded = with_gate(&CancelGate::new(), |cx| {
        initial
            .solve_creep_from(&[0.0; 12], &[material.initial_state()], tolerances(), cx)
            .unwrap()
    });
    let old_u = loaded.solution.displacement_m.as_flattened().to_vec();
    let held: Vec<_> = old_u.iter().copied().enumerate().collect();
    let old_states = loaded.solution.element_states.clone();
    let problem = SmallStrainTetProblem {
        prescribed_m: &held,
        nodal_forces_n: &[],
        ..initial
    };
    with_gate(&CancelGate::new(), |cx| {
        let next = problem
            .solve_creep_from(&old_u, &old_states, tolerances(), cx)
            .unwrap();
        // At fixed strain the Norton branch solves this scalar equation while
        // plasticity remains inactive below the previously hardened yield.
        let shear = E / (2.0 * (1.0 + NU));
        let (mut lo, mut hi) = (0.0, STRESS);
        for _ in 0..80 {
            let q = 0.5 * (lo + hi);
            if q + 3.0 * shear * RATE * (q / STRESS).powi(3) > STRESS {
                hi = q;
            } else {
                lo = q;
            }
        }
        let q = 0.5 * (lo + hi);
        let sigma = next.solution.element_stress[0];
        close(sigma[0] - sigma[1], q, 1e-8);
        close(
            next.solution.element_states[0].plastic.alpha,
            old_states[0].plastic.alpha,
            1e-10,
        );
        assert!(
            next.solution.element_states[0]
                .creep
                .equivalent_creep_strain
                > old_states[0].creep.equivalent_creep_strain
        );
        assert_eq!(next.solution.displacement_m.as_flattened(), old_u);
        close(
            next.energy.applied_work_j + next.energy.support_work_j,
            0.0,
            1e-14,
        );
        assert!(next.energy.final_stored_j < next.energy.initial_stored_j);
        close(next.element_energy[0].plastic_dissipation_j_m3, 0.0, 1e-9);
        let invalid = SmallStrainTetProblem {
            load_path: &[1.0, 1.0],
            ..problem.clone()
        };
        assert!(matches!(
            invalid.solve_creep_from(&old_u, &old_states, tolerances(), cx),
            Err(CreepTetError::Invalid(_))
        ));
    });
    let gate = CancelGate::new();
    gate.request();
    with_gate(&gate, |cx| {
        assert!(matches!(
            problem.solve_creep_from(&old_u, &old_states, tolerances(), cx),
            Err(CreepTetError::Mechanical(SmallStrainTetError::Cancelled))
        ))
    });
    assert_eq!(loaded.solution.element_states, old_states);
    assert_eq!(loaded.solution.displacement_m.as_flattened(), old_u);
}
