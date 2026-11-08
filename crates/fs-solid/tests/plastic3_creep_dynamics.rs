//! G1/G3/G4: original tetrahedral inertia and physical creep share one clock.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::creep::{CreepTolerances, NortonCreep, NortonCreepParameters};
use fs_material::{IsotropicElastic, SmallStrainLaw};
use fs_solid::plastic3::creep::{CreepTetError, CreepTetTolerances};
use fs_solid::{
    IncrementSettings, SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem,
    TetAssemblyBudget,
};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const YOUNG: f64 = 1000.0;
const POISSON: f64 = 0.25;
const RHO: f64 = 60.0;
const DT: f64 = 0.1;
const RATE: f64 = 0.02;
const REFERENCE_STRESS: f64 = 2.0;
const VOLUME: f64 = 1.0 / 6.0;

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 37,
                kernel_id: 179,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx);
    });
}

fn law() -> NortonCreep {
    NortonCreep::new(
        IsotropicElastic::new(YOUNG, POISSON, 0.05).unwrap(),
        NortonCreepParameters {
            reference_stress_pa: REFERENCE_STRESS,
            reference_rate_s: RATE,
            stress_exponent: 1.0,
            activation_energy_j_mol: 0.0,
            reference_temperature_k: 300.0,
            temperature_range_k: [290.0, 350.0],
        },
        300.0,
        DT,
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

fn problem<'a>(
    material: &'a NortonCreep,
    supports: &'a [(usize, f64)],
    load: &'a [(usize, f64)],
) -> SmallStrainTetProblem<'a, NortonCreep> {
    SmallStrainTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        laws: SmallStrainLawField::Uniform(material),
        prescribed_m: supports,
        nodal_forces_n: load,
        load_path: &[1.0],
        budget: TetAssemblyBudget::standard(),
        settings: IncrementSettings::default(),
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn native_transient_creep_matches_independent_mass_and_deviatoric_flow_recurrence() {
    let material = law();
    let supports: Vec<_> = (0..12).filter(|&d| d != 3).map(|d| (d, 0.0)).collect();
    let load = [(3, 0.15)];
    let problem = problem(&material, &supports, &load);
    let mut u = vec![0.0; 12];
    u[3] = 0.002;
    let mut velocity = vec![[0.0; 3]; 4];
    velocity[1][0] = 0.002;
    let mut states = vec![material.initial_state()];
    // Independent one-coordinate reduction of this physical tetrahedron:
    // M_33 = rho V/10; p = K e; q_next = 2 mu(e_next-3 c_old/2)/r,
    // r = 1+3 mu dt A/sigma_ref, and sigma_x = p+2 q_next/3.
    // No production constitutive/tangent evaluator is used for the oracle.
    let mu = YOUNG / (2.0 * (1.0 + POISSON));
    let bulk = YOUNG / (3.0 * (1.0 - 2.0 * POISSON));
    let mass = RHO * VOLUME / 10.0;
    let relaxation = 1.0 + 3.0 * mu * DT * RATE / REFERENCE_STRESS;
    let stiffness = VOLUME * (bulk + 4.0 * mu / (3.0 * relaxation));
    with_cx(|_, cx| {
        for _ in 0..5 {
            let old_u = u[3];
            let old_v = velocity[1][0];
            let old_c = states[0].creep_strain[0];
            let expected_u = (mass * (old_u + DT * old_v)
                + DT * DT * (load[0].1 + VOLUME * 2.0 * mu * old_c / relaxation))
                / (mass + DT * DT * stiffness);
            let expected_v = (expected_u - old_u) / DT;
            let q = 2.0 * mu * (expected_u - 1.5 * old_c) / relaxation;
            assert!(q > 0.0);
            let dc = DT * RATE * q / REFERENCE_STRESS;
            let step = problem
                .advance_creep_backward_euler(
                    &u,
                    &velocity,
                    &states,
                    &[RHO],
                    DT,
                    tolerances(),
                    cx,
                )
                .unwrap();
            close(step.creep.solution.displacement_m[1][0], expected_u, 1e-11);
            close(step.velocity_mps[1][0], expected_v, 1e-10);
            close(step.creep.solution.element_states[0].creep_strain[0], old_c + dc, 1e-11);
            close(step.creep.solution.element_states[0].creep_strain[1], -(old_c + dc) / 2.0, 1e-11);
            close(step.creep.energy.physical_dissipation_j, VOLUME * q * dc, 1e-11);
            close(step.inertia.kinetic_before_j, 0.5 * mass * old_v * old_v, 1e-12);
            close(step.inertia.kinetic_after_j, 0.5 * mass * expected_v * expected_v, 1e-12);
            close(step.inertia.numerical_loss_j, 0.5 * mass * (expected_v - old_v).powi(2), 1e-12);
            close(step.creep.energy.applied_work_j, load[0].1 * (expected_u - old_u), 1e-11);
            close(step.creep.energy.equilibrium_work_residual_j, 0.0, 1e-10);
            close(step.creep.energy.balance_residual_j, 0.0, 1e-10);
            // A repeated Newton solve must not age the borrowed history twice.
            let repeated = problem
                .advance_creep_backward_euler(&u, &velocity, &states, &[RHO], DT, tolerances(), cx)
                .unwrap();
            assert_eq!(step, repeated);
            u = step.creep.solution.displacement_m.into_iter().flatten().collect();
            velocity = step.velocity_mps;
            states = step.creep.solution.element_states;
        }
    });
}

#[test]
fn moving_support_work_includes_mass_and_keeps_kinetic_loss_out_of_heat() {
    let material = law();
    let mut u = vec![0.0; 12];
    u[3] = 0.002;
    let mut velocity = vec![[0.0; 3]; 4];
    velocity[1][0] = 0.002;
    let supports: Vec<_> = (0..12).map(|d| (d, if d == 3 { 0.001 } else { 0.0 })).collect();
    let problem = problem(&material, &supports, &[]);
    with_cx(|_, cx| {
        let step = problem
            .advance_creep_backward_euler(&u, &velocity, &[material.initial_state()], &[RHO], DT, tolerances(), cx)
            .unwrap();
        let mu = YOUNG / (2.0 * (1.0 + POISSON));
        let bulk = YOUNG / (3.0 * (1.0 - 2.0 * POISSON));
        let q = 2.0 * mu * 0.001 / (1.0 + 3.0 * mu * DT * RATE / REFERENCE_STRESS);
        let v = (0.001 - u[3]) / DT;
        let inertial_force = RHO * VOLUME / 10.0 * (v - velocity[1][0]) / DT;
        let reaction = VOLUME * (bulk * 0.001 + 2.0 * q / 3.0) + inertial_force;
        close(step.creep.energy.support_work_j, reaction * (0.001 - u[3]), 1e-11);
        close(step.creep.element_energy[0].caloric_heat_j_m3(), q * DT * RATE * q / REFERENCE_STRESS, 1e-11);
        assert!(step.inertia.numerical_loss_j > 0.0);
        close(step.creep.energy.balance_residual_j, 0.0, 1e-10);
    });
}

#[test]
fn free_rigid_translation_does_not_creep_or_dissipate() {
    let material = law();
    let problem = problem(&material, &[], &[]);
    let velocity = [[0.001, -0.002, 0.003]; 4];
    with_cx(|_, cx| {
        let step = problem
            .advance_creep_backward_euler(&[0.0; 12], &velocity, &[material.initial_state()], &[RHO], DT, tolerances(), cx)
            .unwrap();
        for (u, v) in step.creep.solution.displacement_m.iter().zip(&step.velocity_mps) {
            for d in 0..3 {
                close(u[d], DT * velocity[0][d], 1e-12);
                close(v[d], velocity[0][d], 1e-11);
            }
        }
        close(step.creep.solution.element_states[0].equivalent_creep_strain, 0.0, 1e-12);
        close(step.creep.energy.physical_dissipation_j, 0.0, 1e-12);
        close(step.inertia.numerical_loss_j, 0.0, 1e-12);
        close(step.inertia.kinetic_after_j, RHO * VOLUME * 14.0e-6 / 2.0, 1e-12);
    });
}

#[test]
fn mismatched_clock_and_cancellation_leave_accepted_state_unchanged() {
    let material = law();
    let supports: Vec<_> = (0..12).map(|d| (d, 0.0)).collect();
    let problem = problem(&material, &supports, &[]);
    let u = [0.0; 12];
    let velocity = [[0.0; 3]; 4];
    let state = [material.initial_state()];
    let before = state.clone();
    with_cx(|gate, cx| {
        for dt in [DT / 2.0, f64::NAN, 0.0] {
            assert!(matches!(
                problem.advance_creep_backward_euler(&u, &velocity, &state, &[RHO], dt, tolerances(), cx),
                Err(CreepTetError::Invalid(_))
            ));
        }
        gate.request();
        assert!(matches!(
            problem.advance_creep_backward_euler(&u, &velocity, &state, &[RHO], DT, tolerances(), cx),
            Err(CreepTetError::Mechanical(SmallStrainTetError::Cancelled))
        ));
        assert_eq!(state, before);
    });
}
