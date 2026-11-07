//! Independent physical equations and accounts for both original wall families.

use super::super::{AcousticWallRelaxation, MechanicalTimeLimits, WallPin};
use super::*;
use crate::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_rep_mesh::TetComplex;
use std::cell::Cell;

const MECHANICAL: [u32; 3] = [0, 1, 3];
const ADMITTANCE: [u32; 3] = [0, 2, 3];
const RELEASE: [u32; 3] = [0, 1, 2];
const POSITIONS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const MEDIUM: Medium = Medium {
    density_kg_m3: 2.0,
    sound_speed_m_s: 3.0,
};

fn material<'f, 'c>(geometry: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        geometry,
        4,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap()
}

fn field<'m, 'f, 'c>(
    model: &'m HeterogeneousHelmholtzFe<'f, 'c>,
    law: impl Fn([f64; 3]) -> f64,
) -> SpatialResistiveBoundary<'m, 'f, 'c> {
    model
        .spatial_resistive_boundary(
            &[ADMITTANCE],
            4,
            &mut |_, sample| Ok::<_, AcousticFeError>(law(sample.point.position_m)),
            &mut || false,
        )
        .unwrap()
}

fn walls() -> [MechanicalBoundary; 1] {
    [MechanicalBoundary {
        faces: vec![MECHANICAL],
        wall: WallPin {
            surface_density: 0.6,
            stiffness_per_area: 3.0,
            resistance: 0.4,
        },
        relaxation: vec![AcousticWallRelaxation {
            stiffness_per_area: 2.5,
            relaxation_time_s: 0.06,
        }],
    }]
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn composed_exact_fields_and_seven_accounts_match_independent_physical_ode() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let c0 = field(&model, |_| 0.12);
    let residue = field(&model, |x| 1.0 + x[1]);
    let law = model
        .spatial_relaxing_boundary(Some(&c0), &[(&residue, 0.04)], &mut || false)
        .unwrap();
    let walls = walls();
    let ports = [ExactVelocityPort {
        profile: vec![(ADMITTANCE, [1.0; 3])],
    }];
    let duration = 0.025;
    let owner = model
        .prepare_exact_time_with_walls_and_relaxing_boundaries(
            duration,
            &[],
            &[],
            &walls,
            &[&law],
            &ports,
            &[],
            &[],
            &[RELEASE],
            ExactMechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    // Independent unit-tetra/triangle moments: M33=1/1080, K33=W33=1/12,
    // C0,33=1/100, Cj,33=1/10, and the held inward load b3=0.06/6.
    // The remaining wall/memory trace coordinates are live but start at zero.
    let coupling = det::sqrt(2.5);
    let derivative = |x: [f64; 13]| {
        let [psi, p, w, v, z, q, _, _, _, _, _, _, _] = x;
        let e = coupling * w - z;
        [
            p,
            1080.0 * (0.01 - psi / 12.0 - p / 100.0 - v / 12.0 - q / 10.0),
            v,
            (p - 0.4 * v - 3.0 * w - coupling * e) / 0.6,
            e / 0.06,
            (p - q) / 0.04,
            p * p / 100.0,
            0.4 * v * v / 12.0,
            e * e / (12.0 * 0.06),
            q * q / 10.0,
            p * v / 12.0,
            p * q / 10.0,
            0.01 * p,
        ]
    };
    let mut reference = [
        0.0, 1.0, 0.01, -0.02, 0.03, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ];
    let dt = duration / 4096.0;
    for _ in 0..4096 {
        let a = derivative(reference);
        let b = derivative(core::array::from_fn(|i| reference[i] + 0.5 * dt * a[i]));
        let c = derivative(core::array::from_fn(|i| reference[i] + 0.5 * dt * b[i]));
        let d = derivative(core::array::from_fn(|i| reference[i] + dt * c[i]));
        reference = core::array::from_fn(|i| {
            reference[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
        });
    }
    let initial = owner
        .state_with_all_memories(
            &[0.0; 4],
            &[0.0, 0.0, 0.0, 1.0],
            &[vec![0.0, 0.0, 0.01]],
            &[vec![0.0, 0.0, -0.02]],
            &[vec![vec![0.0, 0.0, 0.03]]],
            &[vec![0.0, 0.0, 0.2]],
        )
        .unwrap();
    let result = owner.step(&initial, &[0.06], &mut || false).unwrap();
    let endpoint = [
        result.state.pressure_potential_pa_s()[3],
        result.state.pressure_pa()[3],
        result.state.wall_displacement_m(0).unwrap()[2],
        result.state.wall_normal_velocity_m_s(0).unwrap()[2],
        result.state.wall_relaxation_memory(0, 0).unwrap()[2],
        result.state.admittance_memory_pressure_pa(0).unwrap()[2],
    ];
    for (&actual, &expected) in endpoint.iter().zip(&reference[..6]) {
        close(actual, expected, 3e-10);
    }
    let accounts = [
        result.work.boundary_dissipation_j,
        result.work.wall_resistance_dissipation_j,
        result.work.wall_relaxation_dissipation_j,
        result.work.admittance_dissipation_j,
        result.work.wall_work_j,
        result.work.admittance_work_j,
        result.work.input_work_j,
    ];
    for (&actual, &expected) in accounts.iter().zip(&reference[6..]) {
        close(actual, expected, 3e-12);
    }
    for residual in [
        result.work.balance_residual_j(),
        result.work.acoustic_balance_residual_j(),
        result.work.wall_balance_residual_j(),
        result.work.admittance_balance_residual_j(),
    ] {
        close(residual, 0.0, 3e-12);
    }
    // Exact interval partitioning must preserve ALL physical states and accounts.
    let half = model
        .prepare_exact_time_with_walls_and_relaxing_boundaries(
            duration / 2.0,
            &[],
            &[],
            &walls,
            &[&law],
            &ports,
            &[],
            &[],
            &[RELEASE],
            ExactMechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let s0 = half
        .state_with_all_memories(
            initial.pressure_potential_pa_s(),
            initial.pressure_pa(),
            &[initial.wall_displacement_m(0).unwrap().to_vec()],
            &[initial.wall_normal_velocity_m_s(0).unwrap().to_vec()],
            &[vec![initial.wall_relaxation_memory(0, 0).unwrap().to_vec()]],
            &[initial.admittance_memory_pressure_pa(0).unwrap().to_vec()],
        )
        .unwrap();
    let first = half.step(&s0, &[0.06], &mut || false).unwrap();
    let second = half.step(&first.state, &[0.06], &mut || false).unwrap();
    for (partitioned, whole) in [
        (
            second.state.pressure_potential_pa_s(),
            result.state.pressure_potential_pa_s(),
        ),
        (second.state.pressure_pa(), result.state.pressure_pa()),
        (
            second.state.wall_displacement_m(0).unwrap(),
            result.state.wall_displacement_m(0).unwrap(),
        ),
        (
            second.state.wall_normal_velocity_m_s(0).unwrap(),
            result.state.wall_normal_velocity_m_s(0).unwrap(),
        ),
        (
            second.state.wall_relaxation_memory(0, 0).unwrap(),
            result.state.wall_relaxation_memory(0, 0).unwrap(),
        ),
        (
            second.state.admittance_memory_pressure_pa(0).unwrap(),
            result.state.admittance_memory_pressure_pa(0).unwrap(),
        ),
    ] {
        for (&actual, &expected) in partitioned.iter().zip(whole) {
            close(actual, expected, 5e-11);
        }
    }
    let totals = [
        first.work.boundary_dissipation_j + second.work.boundary_dissipation_j,
        first.work.wall_resistance_dissipation_j + second.work.wall_resistance_dissipation_j,
        first.work.wall_relaxation_dissipation_j + second.work.wall_relaxation_dissipation_j,
        first.work.admittance_dissipation_j + second.work.admittance_dissipation_j,
        first.work.wall_work_j + second.work.wall_work_j,
        first.work.admittance_work_j + second.work.admittance_work_j,
        first.work.input_work_j + second.work.input_work_j,
    ];
    for (&actual, &expected) in totals.iter().zip(&accounts) {
        close(actual, expected, 3e-12);
    }
}

#[test]
fn both_complete_wall_families_refine_together_under_sparse_midpoint() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let c0 = field(&model, |_| 0.08);
    let residue = field(&model, |x| 0.6 * (1.0 + x[1]));
    let law = model
        .spatial_relaxing_boundary(Some(&c0), &[(&residue, 0.04)], &mut || false)
        .unwrap();
    let walls = walls();
    let duration = 0.03;
    let exact = model
        .prepare_exact_time_with_walls_and_relaxing_boundaries(
            duration,
            &[],
            &[],
            &walls,
            &[&law],
            &[],
            &[],
            &[],
            &[],
            ExactMechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let pressure = [1.0, -0.4, 2.0, 0.3];
    let w = [vec![0.02, -0.01, 0.03]];
    let v = [vec![-0.01, 0.03, 0.01]];
    let z = [vec![vec![0.005, -0.01, 0.02]]];
    let q = [vec![0.2, -0.1, 0.4]];
    let initial = exact
        .state_with_all_memories(&[0.0; 4], &pressure, &w, &v, &z, &q)
        .unwrap();
    let expected = exact.step(&initial, &[], &mut || false).unwrap();
    let mut errors = Vec::new();
    for count in [16, 32, 64] {
        let midpoint = model
            .prepare_time_with_walls_and_relaxing_boundaries(
                duration / count as f64,
                &[],
                &[],
                &walls,
                &[&law],
                &[],
                MechanicalTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let mut state = midpoint
            .state_with_all_memories(&[0.0; 4], &pressure, &w, &v, &z, &q)
            .unwrap();
        let mut heat = [0.0; 3];
        for _ in 0..count {
            let result = midpoint.step(&state, &[], &mut || false).unwrap();
            heat[0] += result.work.wall_dissipation_j - result.work.wall_relaxation_dissipation_j;
            heat[1] += result.work.wall_relaxation_dissipation_j;
            heat[2] += result.work.admittance_dissipation_j;
            state = result.state;
        }
        let mut error: f64 = 0.0;
        for (actual, expected) in [
            (
                state.pressure_potential_pa_s(),
                expected.state.pressure_potential_pa_s(),
            ),
            (state.pressure_pa(), expected.state.pressure_pa()),
            (
                state.wall_displacement_m(0).unwrap(),
                expected.state.wall_displacement_m(0).unwrap(),
            ),
            (
                state.wall_normal_velocity_m_s(0).unwrap(),
                expected.state.wall_normal_velocity_m_s(0).unwrap(),
            ),
            (
                state.wall_relaxation_memory(0, 0).unwrap(),
                expected.state.wall_relaxation_memory(0, 0).unwrap(),
            ),
            (
                state.admittance_memory_pressure_pa(0).unwrap(),
                expected.state.admittance_memory_pressure_pa(0).unwrap(),
            ),
        ] {
            error = actual
                .iter()
                .zip(expected)
                .map(|(a, b)| (a - b).abs())
                .fold(error, f64::max);
        }
        errors.push(error);
        for (actual, target) in heat.into_iter().zip([
            expected.work.wall_resistance_dissipation_j,
            expected.work.wall_relaxation_dissipation_j,
            expected.work.admittance_dissipation_j,
        ]) {
            close(
                actual,
                target,
                3.0 * exact.energy(&initial).unwrap().total_j() / (count * count) as f64,
            );
        }
    }
    assert!(
        errors[0] > 1e-5,
        "nonvacuous complete-field refinement: {errors:?}"
    );
    for pair in errors.windows(2) {
        assert!(
            pair[1] < 0.28 * pair[0],
            "second-order complete-field refinement: {errors:?}"
        );
    }
}

#[test]
fn coupled_zero_residue_release_memories_decay_without_energy_or_state_removal() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let zero = field(&model, |_| 0.0);
    let law = model
        .spatial_relaxing_boundary(None, &[(&zero, 0.0002)], &mut || false)
        .unwrap();
    let walls = walls();
    let prepare = |maximum| {
        model.prepare_exact_time_with_walls_and_relaxing_boundaries(
            0.01,
            &[],
            &[],
            &walls,
            &[&law],
            &[],
            &[],
            &[],
            &[RELEASE],
            ExactMechanicalTimeLimits {
                max_states: maximum,
                ..ExactMechanicalTimeLimits::default()
            },
            &mut || false,
        )
    };
    assert!(matches!(prepare(13), Err(HeterogeneousError::Limit { .. })));
    let owner = prepare(14).unwrap();
    assert_eq!(owner.state_count(), 14);
    assert_eq!(owner.admittance_terms(), 1);
    assert_eq!(owner.admittance_memory_dofs(0).unwrap(), &[0, 2, 3]);
    assert_eq!(owner.admittance_memory_patch_and_term(0), Some((0, 0)));
    assert_eq!(owner.admittance_relaxation_time_s(0), Some(0.0002));
    let w = [vec![0.0; 3]];
    let z = [vec![vec![0.0; 3]]];
    assert!(
        owner
            .state_with_wall_memory(&[0.0; 4], &[0.0; 4], &w, &w, &z)
            .is_err()
    );
    assert!(
        owner
            .state_with_all_memories(&[0.0; 4], &[0.0; 4], &w, &w, &z, &[vec![0.0; 2]])
            .is_err()
    );
    let initial = owner
        .state_with_all_memories(&[0.0; 4], &[0.0; 4], &w, &w, &z, &[vec![1.0, -2.0, 3.0]])
        .unwrap();
    let result = owner.step(&initial, &[], &mut || false).unwrap();
    let decay = det::exp(-50.0);
    for (&actual, original) in result
        .state
        .admittance_memory_pressure_pa(0)
        .unwrap()
        .iter()
        .zip([1.0, -2.0, 3.0])
    {
        close(actual, original * decay, 3e-10 * decay);
    }
    assert_eq!(owner.energy(&result.state).unwrap().total_j(), 0.0);
    assert_eq!(result.work.admittance_dissipation_j, 0.0);
    assert_eq!(result.work.admittance_work_j, 0.0);
    let other = prepare(14).unwrap();
    assert!(other.step(&initial, &[], &mut || false).is_err());
    let mut polls = 0;
    assert!(matches!(
        owner.step(&initial, &[], &mut || {
            polls += 1;
            polls == 3
        }),
        Err(HeterogeneousError::Cancelled)
    ));
    assert_eq!(
        initial.admittance_memory_pressure_pa(0).unwrap(),
        &[1.0, -2.0, 3.0]
    );
}

#[test]
fn all_held_sources_drive_both_families_with_correct_observation_partitions() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let calls = Cell::new(0);
    let residue = field(&model, |x| {
        calls.set(calls.get() + 1);
        1.0 + x[1]
    });
    let law = model
        .spatial_relaxing_boundary(None, &[(&residue, 0.04)], &mut || false)
        .unwrap();
    let walls = walls();
    let ports = [ExactVelocityPort {
        profile: vec![(MECHANICAL, [1.0; 3])],
    }];
    let boundary = model
        .quadratic_velocity_source(&[(MECHANICAL, [1.0; 6])], &mut || false)
        .unwrap();
    let volume = model
        .volume_velocity_source(4, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || {
            false
        })
        .unwrap();
    let owner = model
        .prepare_exact_time_with_walls_and_relaxing_boundaries(
            0.003,
            &[],
            &[],
            &walls,
            &[&law],
            &ports,
            &[&boundary],
            &[&volume],
            &[],
            ExactMechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let initial = owner.zero_state();
    let result = owner
        .step_with_volume_sources(&initial, &[0.2, -0.05], &[0.03], &mut || false)
        .unwrap();
    let face = &result.integrated_port_pressure_pa_m2_s;
    let volume = &result.integrated_volume_source_pressure_pa_s;
    assert_eq!(face.len(), 2);
    assert_eq!(volume.len(), 1);
    close(face[0], face[1], 2e-12);
    assert!(result.work.boundary_source_work_j < 0.0);
    assert!(result.work.volume_source_work_j > 0.0);
    assert!(result.work.wall_resistance_dissipation_j > 0.0);
    assert!(result.work.admittance_dissipation_j > 0.0);
    close(result.work.boundary_source_work_j, -0.05 * face[1], 2e-13);
    close(result.work.volume_source_work_j, 0.03 * volume[0], 2e-13);
    close(
        result.work.input_work_j,
        0.2 * face[0] - 0.05 * face[1] + 0.03 * volume[0],
        2e-13,
    );
    close(result.work.balance_residual_j(), 0.0, 3e-12);
    owner
        .step_with_volume_sources(&result.state, &[-0.1, 0.02], &[-0.01], &mut || false)
        .unwrap();
    assert_eq!(calls.get(), 16);
}
