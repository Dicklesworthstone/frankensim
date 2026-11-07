//! Coupled original wall laws, independently integrated physical trace moments.

use super::*;
use crate::acoustics::helmholtz::heterogeneous::time::relaxing::RelaxingTimeLimits;
use crate::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_rep_mesh::TetComplex;

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
const WALL: [u32; 3] = [0, 1, 2];
const ABSORBER: [u32; 3] = [0, 1, 3];
const DRIVE: [u32; 3] = [1, 2, 3];

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
    coefficients: [f64; 3],
) -> SpatialResistiveBoundary<'m, 'f, 'c> {
    model
        .spatial_resistive_boundary(
            &[ABSORBER],
            4,
            &mut |_, sample| {
                let [a, b, c] = coefficients;
                let p = sample.point.position_m;
                Ok::<_, AcousticFeError>(a + b * p[0] + c * p[2])
            },
            &mut || false,
        )
        .unwrap()
}

fn wall() -> MechanicalBoundary {
    MechanicalBoundary {
        faces: vec![WALL],
        wall: WallPin {
            surface_density: 0.3,
            stiffness_per_area: 7.0,
            resistance: 0.2,
        },
        relaxation: vec![AcousticWallRelaxation {
            stiffness_per_area: 2.0,
            relaxation_time_s: 0.04,
        }],
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

// Exact monomial moments on either unit coordinate triangle, area 1/2:
// integral lambda^alpha = product(alpha!)/(sum(alpha)+2)!.
fn moment(indices: &[usize]) -> f64 {
    let mut counts = [0; 3];
    for &i in indices {
        counts[i] += 1;
    }
    let factorial = |n: usize| (1..=n).product::<usize>() as f64;
    counts.iter().map(|&n| factorial(n)).product::<f64>() / factorial(indices.len() + 2)
}

fn paired(left: &[f64], right: &[f64], coefficient: [f64; 3]) -> f64 {
    let mut result = 0.0;
    for i in 0..3 {
        for j in 0..3 {
            let c = coefficient[0] * moment(&[i, j])
                + coefficient[1] * moment(&[i, j, 1])
                + coefficient[2] * moment(&[i, j, 2]);
            result += left[i] * right[j] * c;
        }
    }
    result
}

fn middle(before: &[f64], after: &[f64]) -> Vec<f64> {
    before
        .iter()
        .zip(after)
        .map(|(&a, &b)| f64::midpoint(a, b))
        .collect()
}

#[test]
fn simultaneous_wall_laws_retain_independent_physical_storage_heat_and_exchange() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let c0 = field(&model, [0.05, 0.05, 0.05]);
    let c1 = field(&model, [0.7, 0.2, 0.0]);
    let c2 = field(&model, [0.3, 0.0, 0.1]);
    let law = model
        .spatial_relaxing_boundary(Some(&c0), &[(&c1, 0.03), (&c2, 0.09)], &mut || false)
        .unwrap();
    let h = 0.001;
    let owner = model
        .prepare_time_with_walls_and_relaxing_boundaries(
            h,
            &[],
            &[],
            &[wall()],
            &[&law],
            &[],
            MechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(owner.wall_dofs(0), Some([0, 1, 2].as_slice()));
    assert_eq!(owner.admittance_memory_dofs(0), Some([0, 1, 3].as_slice()));
    assert_eq!(owner.admittance_terms(), 2);
    let mut state = owner
        .state_with_all_memories(
            &[0.0; 4],
            &[1.2, -0.2, 0.8, 2.0],
            &[vec![0.01, -0.02, 0.03]],
            &[vec![0.02, -0.01, 0.01]],
            &[vec![vec![0.01, 0.02, -0.01]]],
            &[vec![0.5, -0.1, 0.3], vec![-0.2, 0.1, 0.4]],
        )
        .unwrap();
    let unit = [1.0, 0.0, 0.0];
    let k = fs_math::det::sqrt(2.0);
    for _ in 0..20 {
        let initial_energy = owner.energy(&state).unwrap();
        let next = owner
            .step(&state, &[(DRIVE, [0.02; 3])], &mut || false)
            .unwrap();
        let p = middle(state.pressure_pa(), next.state.pressure_pa());
        let p_wall = [p[0], p[1], p[2]];
        let p_absorber = [p[0], p[1], p[3]];
        let v = middle(
            state.wall_normal_velocity_m_s(0).unwrap(),
            next.state.wall_normal_velocity_m_s(0).unwrap(),
        );
        let w = middle(
            state.wall_displacement_m(0).unwrap(),
            next.state.wall_displacement_m(0).unwrap(),
        );
        let z = middle(
            state.wall_relaxation_memory(0, 0).unwrap(),
            next.state.wall_relaxation_memory(0, 0).unwrap(),
        );
        let extension = w
            .iter()
            .zip(&z)
            .map(|(&w, &z)| k * w - z)
            .collect::<Vec<_>>();
        let kelvin_heat = h * 0.2 * paired(&v, &v, unit);
        let maxwell_heat = h / 0.04 * paired(&extension, &extension, unit);
        close(
            next.work.wall_dissipation_j,
            kelvin_heat + maxwell_heat,
            3e-13,
        );
        close(next.work.wall_relaxation_dissipation_j, maxwell_heat, 3e-13);
        close(next.work.wall_work_j, h * paired(&p_wall, &v, unit), 3e-13);
        close(
            next.work.boundary_dissipation_j,
            h * paired(&p_absorber, &p_absorber, [0.05, 0.05, 0.05]),
            3e-13,
        );
        let mut storage = 0.0;
        let mut heat = 0.0;
        let mut exchange = 0.0;
        let mut change = 0.0;
        for (arm, (coefficient, tau)) in [([0.7, 0.2, 0.0], 0.03), ([0.3, 0.0, 0.1], 0.09)]
            .into_iter()
            .enumerate()
        {
            let before = state.admittance_memory_pressure_pa(arm).unwrap();
            let after = next.state.admittance_memory_pressure_pa(arm).unwrap();
            let q = middle(before, after);
            storage += 0.5 * tau * paired(before, before, coefficient);
            heat += h * paired(&q, &q, coefficient);
            exchange += h * paired(&p_absorber, &q, coefficient);
            change += 0.5
                * tau
                * (paired(after, after, coefficient) - paired(before, before, coefficient));
        }
        close(initial_energy.admittance_memory_j, storage, 3e-13);
        close(next.work.admittance_dissipation_j, heat, 3e-13);
        close(next.work.admittance_work_j, exchange, 3e-13);
        close(next.work.admittance_delta_energy_j, change, 3e-13);
        close(
            next.work.input_work_j,
            h * 0.02 * fs_math::det::sqrt(3.0) / 6.0 * (p[1] + p[2] + p[3]),
            3e-13,
        );
        for residual in [
            next.work.balance_residual_j(),
            next.work.acoustic_balance_residual_j(),
            next.work.wall_balance_residual_j(),
            next.work.admittance_balance_residual_j(),
        ] {
            close(residual, 0.0, 5e-13);
        }
        assert!(next.work.backward_error < 1e-10);
        state = next.state;
    }
}

#[test]
fn release_shared_admittance_history_can_drive_both_fluid_and_mechanical_walls() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let residue = field(&model, [1.0, 0.2, 0.0]);
    let tau = 0.04;
    let law = model
        .spatial_relaxing_boundary(None, &[(&residue, tau)], &mut || false)
        .unwrap();
    let h = 0.002;
    let owner = model
        .prepare_time_with_walls_and_relaxing_boundaries(
            h,
            &[],
            &[],
            &[wall()],
            &[&law],
            &[DRIVE],
            MechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let initial = owner
        .state_with_all_memories(
            &[0.0; 4],
            &[0.0; 4],
            &[vec![0.0; 3]],
            &[vec![0.0; 3]],
            &[vec![vec![0.0; 3]]],
            &[vec![0.0, 1.0, 0.0]],
        )
        .unwrap();
    let result = owner.step(&initial, &[], &mut || false).unwrap();
    assert!(result.work.admittance_work_j < 0.0);
    assert!(result.work.acoustic_delta_energy_j > 0.0);
    assert!(result.work.wall_delta_energy_j > 0.0);
    assert!(result.work.admittance_dissipation_j > 0.0);
    close(
        result.state.admittance_memory_pressure_pa(0).unwrap()[1],
        (tau - h / 2.0) / (tau + h / 2.0),
        2e-14,
    );
    for &i in &[1, 2, 3] {
        assert_eq!(result.state.pressure_pa()[i], 0.0);
        assert_eq!(result.state.pressure_potential_pa_s()[i], 0.0);
    }
    close(result.work.input_work_j, 0.0, 0.0);
    close(result.work.balance_residual_j(), 0.0, 5e-13);
}

#[test]
fn zero_residue_memory_is_retained_and_both_memory_families_count_toward_the_cap() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let zero = field(&model, [0.0; 3]);
    let law = model
        .spatial_relaxing_boundary(None, &[(&zero, 0.02)], &mut || false)
        .unwrap();
    let limits = MechanicalTimeLimits {
        max_wall_relaxation_dofs: 5,
        ..MechanicalTimeLimits::default()
    };
    assert!(matches!(
        model.prepare_time_with_walls_and_relaxing_boundaries(
            0.001,
            &[],
            &[],
            &[wall()],
            &[&law],
            &[DRIVE],
            limits,
            &mut || false
        ),
        Err(HeterogeneousError::Limit { .. })
    ));
    let limits = MechanicalTimeLimits {
        max_wall_relaxation_dofs: 6,
        ..limits
    };
    let owner = model
        .prepare_time_with_walls_and_relaxing_boundaries(
            0.001,
            &[],
            &[],
            &[wall()],
            &[&law],
            &[DRIVE],
            limits,
            &mut || false,
        )
        .unwrap();
    let initial = owner
        .state_with_all_memories(
            &[0.0; 4],
            &[0.0; 4],
            &[vec![0.0; 3]],
            &[vec![0.0; 3]],
            &[vec![vec![0.0; 3]]],
            &[vec![1.0, -2.0, 3.0]],
        )
        .unwrap();
    let result = owner.step(&initial, &[], &mut || false).unwrap();
    for (&q, original) in result
        .state
        .admittance_memory_pressure_pa(0)
        .unwrap()
        .iter()
        .zip([1.0, -2.0, 3.0])
    {
        close(q, original * 0.0195 / 0.0205, 2e-14);
    }
    assert_eq!(owner.energy(&result.state).unwrap().total_j(), 0.0);
    assert_eq!(result.work.admittance_dissipation_j, 0.0);
    assert_eq!(result.work.admittance_work_j, 0.0);
    assert!(
        owner
            .state_with_wall_memory(
                &[0.0; 4],
                &[0.0; 4],
                &[vec![0.0; 3]],
                &[vec![0.0; 3]],
                &[vec![vec![0.0; 3]]]
            )
            .is_err()
    );
    assert!(
        owner
            .state_with_all_memories(
                &[0.0; 4],
                &[0.0; 4],
                &[vec![0.0; 3]],
                &[vec![0.0; 3]],
                &[vec![vec![0.0; 3]]],
                &[vec![0.0; 2]]
            )
            .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        owner.step(&initial, &[], &mut || {
            polls += 1;
            polls == 4
        }),
        Err(HeterogeneousError::Cancelled)
    ));
    assert_eq!(
        initial.admittance_memory_pressure_pa(0).unwrap(),
        &[1.0, -2.0, 3.0]
    );
    assert!(
        model
            .prepare_time_with_walls_and_relaxing_boundaries(
                0.001,
                &[],
                &[],
                &[MechanicalBoundary {
                    faces: vec![ABSORBER],
                    ..wall()
                }],
                &[&law],
                &[],
                limits,
                &mut || false
            )
            .is_err()
    );
}

#[test]
fn removing_mechanical_patches_recovers_the_original_spatial_memory_source_response() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let c0 = field(&model, [0.05, 0.01, 0.0]);
    let residue = field(&model, [0.4, 0.0, 0.2]);
    let law = model
        .spatial_relaxing_boundary(Some(&c0), &[(&residue, 0.03)], &mut || false)
        .unwrap();
    let h = 0.001;
    let combined = model
        .prepare_time_with_walls_and_relaxing_boundaries(
            h,
            &[],
            &[],
            &[],
            &[&law],
            &[],
            MechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let original = model
        .prepare_time_with_relaxing_boundaries(
            h,
            &[],
            &[],
            &[&law],
            &[],
            RelaxingTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let source = model
        .quadratic_velocity_source(&[(WALL, [1.0; 6])], &mut || false)
        .unwrap();
    let volume = model
        .volume_velocity_source(4, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || {
            false
        })
        .unwrap();
    let pressure = [1.0, -0.2, 0.3, 0.7];
    let q = vec![0.2, -0.1, 0.4];
    let mut a = combined
        .state_with_all_memories(&[0.0; 4], &pressure, &[], &[], &[], &[q.clone()])
        .unwrap();
    let mut b = original
        .state_with_memories(&[0.0; 4], &pressure, &[&q])
        .unwrap();
    for _ in 0..12 {
        let next_a = combined
            .step_with_boundary_sources(
                &a,
                &[(DRIVE, [0.01; 3])],
                &[(&source, -0.02)],
                &[(&volume, 0.005)],
                &mut || false,
            )
            .unwrap();
        let next_b = original
            .step_with_boundary_sources(
                &b,
                &[(DRIVE, [0.01; 3])],
                &[(&source, -0.02)],
                &[(&volume, 0.005)],
                &mut || false,
            )
            .unwrap();
        for (left, right) in [
            (next_a.state.pressure_pa(), next_b.state.pressure_pa()),
            (
                next_a.state.pressure_potential_pa_s(),
                next_b.state.pressure_potential_pa_s(),
            ),
            (
                next_a.state.admittance_memory_pressure_pa(0).unwrap(),
                next_b.state.memory_pressure_pa(0).unwrap(),
            ),
        ] {
            for (&x, &y) in left.iter().zip(right) {
                close(x, y, 2e-13);
            }
        }
        for (x, y) in [
            (next_a.work.admittance_work_j, next_b.work.memory_work_j),
            (
                next_a.work.admittance_dissipation_j,
                next_b.work.relaxation_dissipation_j,
            ),
            (next_a.work.input_work_j, next_b.work.input_work_j),
            (
                next_a.work.boundary_source_work_j,
                next_b.work.boundary_source_work_j,
            ),
            (
                next_a.work.volume_source_work_j,
                next_b.work.volume_source_work_j,
            ),
        ] {
            close(x, y, 2e-13);
        }
        a = next_a.state;
        b = next_b.state;
    }
}
