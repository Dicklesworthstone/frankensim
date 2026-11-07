//! G0/G1/G3: original traction work, independent equations and full-state refinement.

use super::super::{AcousticWallRelaxation, MechanicalTimeLimits, WallPin};
use super::*;
use crate::acoustics::helmholtz::{AcousticFeError, HelmholtzFe, Medium};
use fs_rep_mesh::TetComplex;

const WALL: [u32; 3] = [0, 1, 3];
const ABSORBING: [u32; 3] = [0, 2, 3];
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
    value: impl Fn([f64; 3]) -> f64,
) -> SpatialResistiveBoundary<'m, 'f, 'c> {
    model
        .spatial_resistive_boundary(
            &[ABSORBING],
            4,
            &mut |_, sample| Ok::<_, AcousticFeError>(value(sample.point.position_m)),
            &mut || false,
        )
        .unwrap()
}

fn walls() -> [MechanicalBoundary; 1] {
    [MechanicalBoundary {
        faces: vec![WALL],
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

fn original_fields(state: &ExactMechanicalTimeState<'_, '_, '_, '_>) -> Vec<f64> {
    state
        .pressure_potential_pa_s()
        .iter()
        .chain(state.pressure_pa())
        .chain(state.wall_displacement_m(0).unwrap())
        .chain(state.wall_normal_velocity_m_s(0).unwrap())
        .chain(state.wall_relaxation_memory(0, 0).unwrap())
        .chain(state.admittance_memory_pressure_pa(0).unwrap_or(&[]))
        .copied()
        .collect()
}

#[test]
fn applied_traction_and_all_wall_memories_match_independent_ode_and_interval_partition() {
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
        profile: vec![(ABSORBING, [1.0; 3])],
    }];
    let tractions = [WallTractionProfile {
        wall: 0,
        coefficients: vec![(3, 1.3)],
    }];
    let prepare = |duration| {
        model
            .prepare_exact_time_with_wall_tractions(
                duration,
                &[],
                &[],
                &walls,
                &[&law],
                &ports,
                &[],
                &[],
                &tractions,
                &[RELEASE],
                ExactMechanicalTimeLimits::default(),
                &mut || false,
            )
            .unwrap()
    };
    let duration = 0.025;
    let owner = prepare(duration);
    assert_eq!(owner.wall_traction_count(), 1);
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
    let result = owner
        .step_with_wall_tractions(&initial, &[0.06], &[], &[0.7], &mut || false)
        .unwrap();
    // Unit-simplex moments, specified independently of assembled operators:
    // M33=1/1080, K33=W33=1/12, C0,33=1/100, Cj,33=1/10.
    // The applied traction is 1.3*0.7*phi3 Pa; the held fluid load is 0.06/6.
    let coupling = 2.5_f64.sqrt();
    let derivative = |x: [f64; 15]| {
        let [psi, p, w, v, z, q, _, _, _, _, _, _, _, _, _] = x;
        let e = coupling * w - z;
        [
            p,
            1080.0 * (0.01 - psi / 12.0 - p / 100.0 - v / 12.0 - q / 10.0),
            v,
            (p + 1.3 * 0.7 - 0.4 * v - 3.0 * w - coupling * e) / 0.6,
            e / 0.06,
            (p - q) / 0.04,
            p * p / 100.0,
            0.4 * v * v / 12.0,
            e * e / (12.0 * 0.06),
            q * q / 10.0,
            p * v / 12.0,
            p * q / 10.0,
            0.01 * p,
            0.7 * 1.3 * v / 12.0,
            1.3 * v / 12.0,
        ]
    };
    let mut reference = [0.0; 15];
    reference[..6].copy_from_slice(&[0.0, 1.0, 0.01, -0.02, 0.03, 0.2]);
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
    let accounts = |step: &ExactMechanicalTimeStep<'_, '_, '_, '_>| {
        [
            step.work.boundary_dissipation_j,
            step.work.wall_resistance_dissipation_j,
            step.work.wall_relaxation_dissipation_j,
            step.work.admittance_dissipation_j,
            step.work.wall_work_j,
            step.work.admittance_work_j,
            step.work.acoustic_input_work_j,
            step.work.wall_traction_work_j,
            step.integrated_wall_traction_displacement_m3[0],
        ]
    };
    for (&actual, &expected) in accounts(&result).iter().zip(&reference[6..]) {
        close(actual, expected, 3e-12);
    }
    close(
        result.work.input_work_j,
        reference[12] + reference[13],
        3e-12,
    );
    close(
        result.integrated_wall_traction_displacement_m3[0],
        1.3 * (result.state.wall_displacement_m(0).unwrap()[2] - 0.01) / 12.0,
        3e-13,
    );
    for residual in [
        result.work.balance_residual_j(),
        result.work.acoustic_balance_residual_j(),
        result.work.wall_balance_residual_j(),
        result.work.admittance_balance_residual_j(),
    ] {
        close(residual, 0.0, 3e-12);
    }
    let half = prepare(duration / 2.0);
    let half_initial = half
        .state_with_all_memories(
            initial.pressure_potential_pa_s(),
            initial.pressure_pa(),
            &[initial.wall_displacement_m(0).unwrap().to_vec()],
            &[initial.wall_normal_velocity_m_s(0).unwrap().to_vec()],
            &[vec![initial.wall_relaxation_memory(0, 0).unwrap().to_vec()]],
            &[initial.admittance_memory_pressure_pa(0).unwrap().to_vec()],
        )
        .unwrap();
    let first = half
        .step_with_wall_tractions(&half_initial, &[0.06], &[], &[0.7], &mut || false)
        .unwrap();
    let second = half
        .step_with_wall_tractions(&first.state, &[0.06], &[], &[0.7], &mut || false)
        .unwrap();
    for (actual, expected) in original_fields(&second.state)
        .into_iter()
        .zip(original_fields(&result.state))
    {
        close(actual, expected, 5e-11);
    }
    for ((a, b), total) in accounts(&first)
        .into_iter()
        .zip(accounts(&second))
        .zip(accounts(&result))
    {
        close(a + b, total, 3e-12);
    }
}

#[test]
fn distributed_traction_refines_with_the_complete_midpoint_wall_and_admittance_state() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let c0 = field(&model, |_| 0.08);
    let residue = field(&model, |x| 0.6 * (1.0 + x[1]));
    let law = model
        .spatial_relaxing_boundary(Some(&c0), &[(&residue, 0.04)], &mut || false)
        .unwrap();
    let walls = walls();
    let profile = WallTractionProfile {
        wall: 0,
        coefficients: vec![(0, 0.7), (1, -0.2), (3, 1.1)],
    };
    let duration = 0.03;
    let exact = model
        .prepare_exact_time_with_wall_tractions(
            duration,
            &[],
            &[],
            &walls,
            &[&law],
            &[],
            &[],
            &[],
            &[profile.clone()],
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
    let expected = exact
        .step_with_wall_tractions(&initial, &[], &[], &[0.8], &mut || false)
        .unwrap();
    let target = original_fields(&expected.state);
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
        let traction = midpoint
            .prepare_wall_traction(&profile, &mut || false)
            .unwrap();
        let mut state = midpoint
            .state_with_all_memories(&[0.0; 4], &pressure, &w, &v, &z, &q)
            .unwrap();
        let mut traction_work = 0.0;
        for _ in 0..count {
            let step = midpoint
                .step_with_wall_tractions(&state, &[], &[], &[], &[(&traction, 0.8)], &mut || false)
                .unwrap();
            traction_work += step.work.wall_traction_work_j;
            state = step.state;
        }
        let actual = state
            .pressure_potential_pa_s()
            .iter()
            .chain(state.pressure_pa())
            .chain(state.wall_displacement_m(0).unwrap())
            .chain(state.wall_normal_velocity_m_s(0).unwrap())
            .chain(state.wall_relaxation_memory(0, 0).unwrap())
            .chain(state.admittance_memory_pressure_pa(0).unwrap());
        let error = actual
            .zip(&target)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        errors.push(error);
        close(
            traction_work,
            expected.work.wall_traction_work_j,
            3.0 * exact.energy(&initial).unwrap().total_j() / (count * count) as f64,
        );
    }
    assert!(
        errors[0] > 1e-5,
        "nonvacuous complete-state error: {errors:?}"
    );
    for pair in errors.windows(2) {
        assert!(pair[1] < 0.28 * pair[0], "midpoint refinement: {errors:?}");
    }
}

#[test]
fn coherent_tractions_keep_original_source_order_units_and_independent_acoustic_work() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let walls = walls();
    let ports = [ExactVelocityPort {
        profile: vec![(WALL, [1.0; 3])],
    }];
    let boundary = model
        .quadratic_velocity_source(&[(WALL, [1.0; 6])], &mut || false)
        .unwrap();
    let volume = model
        .volume_velocity_source(4, &mut |_, _| Ok::<_, HeterogeneousError>(1.0), &mut || {
            false
        })
        .unwrap();
    let profile = WallTractionProfile {
        wall: 0,
        coefficients: vec![(0, 0.5), (1, -0.4), (3, 1.2)],
    };
    let opposite = WallTractionProfile {
        wall: 0,
        coefficients: profile
            .coefficients
            .iter()
            .map(|&(i, x)| (i, -2.0 * x))
            .collect(),
    };
    let prepare = |profiles: &[WallTractionProfile]| {
        model
            .prepare_exact_time_with_wall_tractions(
                0.003,
                &[],
                &[],
                &walls,
                &[],
                &ports,
                &[&boundary],
                &[&volume],
                profiles,
                &[],
                ExactMechanicalTimeLimits::default(),
                &mut || false,
            )
            .unwrap()
    };
    let owner = prepare(&[profile.clone(), opposite]);
    let combined = prepare(&[profile]);
    let result = owner
        .step_with_wall_tractions(
            &owner.zero_state(),
            &[0.02, -0.005],
            &[0.003],
            &[2.0, -0.5],
            &mut || false,
        )
        .unwrap();
    let reference = combined
        .step_with_wall_tractions(
            &combined.zero_state(),
            &[0.02, -0.005],
            &[0.003],
            &[3.0],
            &mut || false,
        )
        .unwrap();
    for (actual, expected) in original_fields(&result.state)
        .into_iter()
        .zip(original_fields(&reference.state))
    {
        close(actual, expected, 3e-11);
    }
    let face = &result.integrated_port_pressure_pa_m2_s;
    let volume_integral = &result.integrated_volume_source_pressure_pa_s;
    let wall_integral = &result.integrated_wall_traction_displacement_m3;
    assert_eq!(
        (face.len(), volume_integral.len(), wall_integral.len()),
        (2, 1, 2)
    );
    close(face[0], face[1], 2e-12);
    close(wall_integral[1], -2.0 * wall_integral[0], 2e-12);
    close(result.work.boundary_source_work_j, -0.005 * face[1], 2e-13);
    close(
        result.work.volume_source_work_j,
        0.003 * volume_integral[0],
        2e-13,
    );
    close(
        result.work.wall_traction_work_j,
        2.0 * wall_integral[0] - 0.5 * wall_integral[1],
        2e-13,
    );
    close(
        result.work.acoustic_input_work_j,
        0.02 * face[0] - 0.005 * face[1] + 0.003 * volume_integral[0],
        2e-13,
    );
    close(
        result.work.input_work_j,
        result.work.acoustic_input_work_j + result.work.wall_traction_work_j,
        2e-13,
    );
    // Independently integrate the original triangle's shape times displacement.
    // W_ii=1/12 and W_ij=1/24; no absolute areas or normalized source gains.
    let coefficients = [0.5, -0.4, 1.2];
    let displacement = result.state.wall_displacement_m(0).unwrap();
    let mut physical_displacement = 0.0;
    for (i, &coefficient) in coefficients.iter().enumerate() {
        for (j, &w) in displacement.iter().enumerate() {
            physical_displacement += coefficient * w / if i == j { 12.0 } else { 24.0 };
        }
    }
    close(wall_integral[0], physical_displacement, 3e-13);
    // Weak fluid work survives a much larger mechanical drive. Recovering it
    // by subtracting total minus mechanical input would round it away.
    let weak = combined
        .step_with_wall_tractions(
            &combined.zero_state(),
            &[1e-21, 0.0],
            &[0.0],
            &[1.0],
            &mut || false,
        )
        .unwrap();
    let expected_weak = 1e-21 * weak.integrated_port_pressure_pa_m2_s[0];
    assert!(expected_weak != 0.0);
    assert!(expected_weak.abs() < 0.25 * f64::EPSILON * weak.work.wall_traction_work_j.abs());
    close(
        weak.work.acoustic_input_work_j,
        expected_weak,
        4.0 * f64::EPSILON * expected_weak.abs(),
    );
}

#[test]
fn release_shared_wall_forcing_and_source_admission_keep_complete_original_ownership() {
    let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&mesh, &POSITIONS, 1, MEDIUM).unwrap();
    let model = material(&geometry);
    let walls = walls();
    let profile = WallTractionProfile {
        wall: 0,
        coefficients: vec![(0, 1.0)],
    };
    let prepare = |profiles: &[WallTractionProfile], limits| {
        model.prepare_exact_time_with_wall_tractions(
            0.01,
            &[],
            &[],
            &walls,
            &[],
            &[],
            &[],
            &[],
            profiles,
            &[RELEASE],
            limits,
            &mut || false,
        )
    };
    let owner = prepare(&[profile.clone()], ExactMechanicalTimeLimits::default()).unwrap();
    let initial = owner.zero_state();
    let initial_fields = original_fields(&initial);
    let result = owner
        .step_with_wall_tractions(&initial, &[], &[], &[1.0], &mut || false)
        .unwrap();
    assert_eq!(result.state.pressure_pa()[0], 0.0);
    assert!(result.state.wall_normal_velocity_m_s(0).unwrap()[0] > 0.0);
    assert!(result.state.pressure_pa()[3].abs() > 1e-6);
    assert!(result.work.wall_traction_work_j > 0.0);
    assert_eq!(result.work.acoustic_input_work_j, 0.0);
    assert!(owner.step(&initial, &[], &mut || false).is_err());
    assert!(
        owner
            .step_with_volume_sources(&initial, &[], &[], &mut || false)
            .is_err()
    );
    assert!(
        owner
            .step_with_wall_tractions(&initial, &[], &[], &[f64::NAN], &mut || false)
            .is_err()
    );
    let other = prepare(&[profile.clone()], ExactMechanicalTimeLimits::default()).unwrap();
    assert!(
        other
            .step_with_wall_tractions(&initial, &[], &[], &[1.0], &mut || false)
            .is_err()
    );
    assert!(matches!(
        owner.step_with_wall_tractions(&initial, &[], &[], &[1.0], &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    assert_eq!(original_fields(&initial), initial_fields);
    for invalid_profile in [
        WallTractionProfile {
            wall: 1,
            coefficients: vec![],
        },
        WallTractionProfile {
            wall: 0,
            coefficients: vec![(0, 1.0), (0, -1.0)],
        },
        WallTractionProfile {
            wall: 0,
            coefficients: vec![(2, 1.0)],
        },
        WallTractionProfile {
            wall: 0,
            coefficients: vec![(0, f64::NAN)],
        },
    ] {
        assert!(prepare(&[invalid_profile], ExactMechanicalTimeLimits::default()).is_err());
    }
    let mut limited = ExactMechanicalTimeLimits::default();
    limited.preparation.max_dimension = owner.state_count() + 1;
    assert!(matches!(
        prepare(&[profile], limited),
        Err(HeterogeneousError::Limit { .. })
    ));
}
