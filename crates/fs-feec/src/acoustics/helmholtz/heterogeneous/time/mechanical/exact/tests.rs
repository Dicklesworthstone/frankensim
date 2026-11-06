use super::super::{AcousticWallRelaxation, MechanicalTimeLimits, WallPin};
use super::*;
use crate::acoustics::helmholtz::HelmholtzFe;
use crate::acoustics::helmholtz::heterogeneous::time::tests::{MEDIUM, brick, p2, plane, uniform};
use crate::fixtures::single_tet;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn equivalent(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    let scale = expected
        .iter()
        .map(|v| v.abs())
        .fold(0.0, f64::max)
        .max(1e-13);
    for (&actual, &expected) in actual.iter().zip(expected) {
        close(actual, expected, tolerance * scale);
    }
}

fn wall(
    pin: WallPin,
    faces: Vec<[u32; 3]>,
    relaxation: Vec<AcousticWallRelaxation>,
) -> MechanicalBoundary {
    MechanicalBoundary {
        faces,
        wall: pin,
        relaxation,
    }
}

#[test]
fn analytic_coupled_tetrahedron_mode_has_no_midpoint_frequency_warp() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let sigma = 0.6;
    let stiffness = 3.0;
    let boundaries = [wall(
        WallPin {
            surface_density: sigma,
            stiffness_per_area: stiffness,
            resistance: 0.0,
        },
        vec![[0, 1, 3]],
        Vec::new(),
    )];
    // The release face leaves only phi3=z in the fluid. Its original mass is
    // V/(10 rho c²), stiffness V/rho, V=1/6. The wall face area is 1/2 and
    // W33=A/6=1/12. Other wall coefficients remain live, initially zero.
    // Thus psi'' + a psi + b w' = 0, sigma w'' + K w = psi'.
    // Squared frequencies solve s²-(a+K/sigma+b/sigma)s+a K/sigma=0.
    let a = 10.0 * MEDIUM.sound_speed_m_s.powi(2);
    let b = 5.0 * MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2);
    let k = stiffness / sigma;
    let sum = a + k + b / sigma;
    let root = det::sqrt(sum * sum - 4.0 * a * k);
    let duration = 0.71;
    for squared_frequency in [0.5 * (sum + root), 2.0 * a * k / (sum + root)] {
        let omega = det::sqrt(squared_frequency);
        for count in [1, 3] {
            let owner = model
                .prepare_exact_time_with_walls(
                    duration / f64::from(count),
                    &[],
                    &boundaries,
                    &[],
                    &[],
                    &[],
                    &[[0, 1, 2]],
                    ExactMechanicalTimeLimits {
                        max_states: 8,
                        ..ExactMechanicalTimeLimits::default()
                    },
                    &mut || false,
                )
                .unwrap();
            assert_eq!(owner.state_count(), 8);
            assert_eq!(owner.wall_dofs(0).unwrap(), [0, 1, 3]);
            let w0 = 0.02;
            let p0 = sigma * (k - squared_frequency) * w0;
            let mut state = owner
                .state_with_walls(
                    &[0.0; 4],
                    &[0.0, 0.0, 0.0, p0],
                    &[vec![0.0, 0.0, w0]],
                    &[vec![0.0; 3]],
                )
                .unwrap();
            let initial_energy = 0.5 * p0 * p0
                / (60.0 * MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2))
                + stiffness * w0 * w0 / 24.0;
            close(
                owner.energy(&state).unwrap().total_j(),
                initial_energy,
                3e-12 * initial_energy,
            );
            let mut exchange = 0.0;
            for _ in 0..count {
                let result = owner.step(&state, &[], &mut || false).unwrap();
                close(
                    result.work.wall_dissipation_j(),
                    0.0,
                    1e-13 * initial_energy,
                );
                close(
                    result.work.balance_residual_j(),
                    0.0,
                    3e-10 * initial_energy,
                );
                close(
                    result.work.acoustic_balance_residual_j(),
                    0.0,
                    3e-10 * initial_energy,
                );
                close(
                    result.work.wall_balance_residual_j(),
                    0.0,
                    3e-10 * initial_energy,
                );
                exchange += result.work.wall_work_j;
                state = result.state;
            }
            let phase = omega * duration;
            let pressure = p0 * det::cos(phase);
            let potential = p0 * det::sin(phase) / omega;
            let displacement = w0 * det::cos(phase);
            let velocity = -omega * w0 * det::sin(phase);
            equivalent(state.pressure_pa(), &[0.0, 0.0, 0.0, pressure], 3e-9);
            equivalent(
                state.pressure_potential_pa_s(),
                &[0.0, 0.0, 0.0, potential],
                3e-9,
            );
            equivalent(
                state.wall_displacement_m(0).unwrap(),
                &[0.0, 0.0, displacement],
                3e-9,
            );
            equivalent(
                state.wall_normal_velocity_m_s(0).unwrap(),
                &[0.0, 0.0, velocity],
                3e-9,
            );
            // Independent physical wall energy difference, not a total deficit.
            let wall_change = (sigma * velocity * velocity
                + stiffness * (displacement * displacement - w0 * w0))
                / 24.0;
            close(exchange, wall_change, 3e-9 * initial_energy);
        }
    }
}

#[test]
fn fast_self_equilibrated_maxwell_history_decays_exponentially_and_makes_exact_heat() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let tau = 0.01;
    let duration = 8.0 * tau;
    let boundaries = [wall(
        WallPin {
            surface_density: 0.6,
            stiffness_per_area: 3.0,
            resistance: 0.2,
        },
        vec![[0, 1, 3]],
        vec![
            AcousticWallRelaxation {
                stiffness_per_area: 2.0,
                relaxation_time_s: tau,
            },
            AcousticWallRelaxation {
                stiffness_per_area: 8.0,
                relaxation_time_s: tau,
            },
        ],
    )];
    let owner = model
        .prepare_exact_time_with_walls(
            duration,
            &[],
            &boundaries,
            &[],
            &[],
            &[],
            &[[0, 1, 2]],
            ExactMechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(owner.state_count(), 14);
    let first = vec![0.3, -0.2, 0.6];
    let second = first.iter().map(|z| -0.5 * z).collect::<Vec<_>>();
    let state = owner
        .state_with_wall_memory(
            &[0.0; 4],
            &[0.0; 4],
            &[vec![0.0; 3]],
            &[vec![0.0; 3]],
            &[vec![first.clone(), second.clone()]],
        )
        .unwrap();
    // sqrt(2)e1 + sqrt(8)e2 = 0 everywhere, so fluid and wall motion remain
    // exactly zero while BOTH original histories relax. With A=1/2:
    // z^T W z = ((sum z)^2 + sum z²)/24, not an area-averaged piston.
    let moment = |values: &[f64]| {
        (values.iter().sum::<f64>().powi(2) + values.iter().map(|v| v * v).sum::<f64>()) / 24.0
    };
    let energy = 0.5 * (moment(&first) + moment(&second));
    close(
        owner.energy(&state).unwrap().wall_relaxation_j,
        energy,
        3e-12 * energy,
    );
    let result = owner.step(&state, &[], &mut || false).unwrap();
    let decay = det::exp(-duration / tau);
    for (arm, original) in [first, second].iter().enumerate() {
        equivalent(
            result.state.wall_relaxation_memory(0, arm).unwrap(),
            &original.iter().map(|z| z * decay).collect::<Vec<_>>(),
            3e-8,
        );
    }
    for &v in result
        .state
        .pressure_pa()
        .iter()
        .chain(result.state.pressure_potential_pa_s())
        .chain(result.state.wall_displacement_m(0).unwrap())
        .chain(result.state.wall_normal_velocity_m_s(0).unwrap())
    {
        close(v, 0.0, 3e-12);
    }
    close(
        result.work.wall_relaxation_dissipation_j,
        energy * (1.0 - decay * decay),
        3e-9 * energy,
    );
    close(
        owner.energy(&result.state).unwrap().total_j(),
        energy * decay * decay,
        3e-9 * energy,
    );
    close(result.work.wall_work_j, 0.0, 3e-11 * energy);
    close(
        result.work.wall_resistance_dissipation_j,
        0.0,
        3e-11 * energy,
    );
    close(result.work.balance_residual_j(), 0.0, 3e-9 * energy);
    assert!(
        owner
            .state_with_walls(&[0.0; 4], &[0.0; 4], &[vec![0.0; 3]], &[vec![0.0; 3]])
            .is_err()
    );
}

#[test]
fn continuum_constant_wall_velocity_integrates_original_maxwell_and_source_work() {
    let pin = WallPin {
        surface_density: 0.4,
        stiffness_per_area: 1.7,
        resistance: 0.3,
    };
    let relaxation = vec![
        AcousticWallRelaxation {
            stiffness_per_area: 2.5,
            relaxation_time_s: 0.03,
        },
        AcousticWallRelaxation {
            stiffness_per_area: 1.2,
            relaxation_time_s: 0.2,
        },
    ];
    for height in [0.1, 0.3] {
        let lengths = [1.0, height, 0.2];
        let area = height * lengths[2];
        let (complex, positions) = brick(1, lengths);
        let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
        let model = uniform(&geometry);
        let boundaries = [wall(pin, plane(&geometry, 1.0), relaxation.clone())];
        let duration = 0.17;
        let ports = [ExactVelocityPort {
            profile: plane(&geometry, 0.0)
                .into_iter()
                .map(|face| (face, [1.0; 3]))
                .collect(),
        }];
        let owner = model
            .prepare_exact_time_with_walls(
                duration,
                &[],
                &boundaries,
                &ports,
                &[],
                &[],
                &[],
                ExactMechanicalTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        for amplitude in [1.0, 1e-8] {
            let velocity = 0.02 * amplitude;
            let w0 = 0.03 * amplitude;
            let maxwell_resistance = relaxation
                .iter()
                .map(|arm| arm.stiffness_per_area * arm.relaxation_time_s)
                .sum::<f64>();
            let p0 = (pin.resistance + maxwell_resistance) * velocity + pin.stiffness_per_area * w0;
            let a = pin.stiffness_per_area * velocity / MEDIUM.sound_speed_m_s.powi(2);
            let b = -MEDIUM.density_kg_m3 * velocity - a;
            let potential = p2(&geometry, |x| 0.5 * a * x[0] * x[0] + b * x[0]);
            let pressure = p2(&geometry, |_| p0);
            let restrict = |field: Vec<f64>| {
                owner
                    .wall_dofs(0)
                    .unwrap()
                    .iter()
                    .map(|&dof| field[dof])
                    .collect::<Vec<_>>()
            };
            let displacement = restrict(p2(&geometry, |_| w0));
            let wall_velocity = restrict(p2(&geometry, |_| velocity));
            let memory = relaxation
                .iter()
                .map(|arm| {
                    restrict(p2(&geometry, |_| {
                        det::sqrt(arm.stiffness_per_area) * (w0 - arm.relaxation_time_s * velocity)
                    }))
                })
                .collect::<Vec<_>>();
            let initial = owner
                .state_with_wall_memory(
                    &potential,
                    &pressure,
                    &[displacement],
                    &[wall_velocity],
                    &[memory],
                )
                .unwrap();
            let inward = -b / MEDIUM.density_kg_m3;
            let result = owner.step(&initial, &[inward], &mut || false).unwrap();
            let w = w0 + duration * velocity;
            let p = p0 + pin.stiffness_per_area * velocity * duration;
            let integral_p =
                p0 * duration + 0.5 * pin.stiffness_per_area * velocity * duration * duration;
            equivalent(result.state.pressure_pa(), &p2(&geometry, |_| p), 3e-9);
            equivalent(
                result.state.pressure_potential_pa_s(),
                &p2(&geometry, |x| 0.5 * a * x[0] * x[0] + b * x[0] + integral_p),
                3e-9,
            );
            equivalent(
                result.state.wall_displacement_m(0).unwrap(),
                &restrict(p2(&geometry, |_| w)),
                3e-9,
            );
            equivalent(
                result.state.wall_normal_velocity_m_s(0).unwrap(),
                &restrict(p2(&geometry, |_| velocity)),
                3e-9,
            );
            for (index, arm) in relaxation.iter().enumerate() {
                equivalent(
                    result.state.wall_relaxation_memory(0, index).unwrap(),
                    &restrict(p2(&geometry, |_| {
                        det::sqrt(arm.stiffness_per_area) * (w - arm.relaxation_time_s * velocity)
                    })),
                    3e-9,
                );
            }
            let heat_scale = area * duration * velocity * velocity;
            close(
                result.work.wall_resistance_dissipation_j,
                pin.resistance * heat_scale,
                3e-9 * heat_scale,
            );
            close(
                result.work.wall_relaxation_dissipation_j,
                maxwell_resistance * heat_scale,
                3e-9 * heat_scale,
            );
            close(
                result.work.wall_work_j,
                area * velocity * integral_p,
                3e-9 * area * velocity * integral_p,
            );
            close(
                result.integrated_port_pressure_pa_m2_s[0],
                area * integral_p,
                3e-9 * area * integral_p,
            );
            close(
                result.work.input_work_j,
                inward * area * integral_p,
                3e-9 * inward * area * integral_p,
            );
            close(
                result.work.balance_residual_j(),
                0.0,
                3e-9 * result.work.initial_energy_j,
            );
        }
    }
}

#[test]
fn full_quadratic_trace_obeys_semigroup_and_midpoint_converges_to_it() {
    let lengths = [1.0, 0.2, 0.3];
    let area = lengths[1] * lengths[2];
    let (complex, positions) = brick(1, lengths);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let pin = WallPin {
        surface_density: 0.6,
        stiffness_per_area: 3.0,
        resistance: 0.2,
    };
    let arm = AcousticWallRelaxation {
        stiffness_per_area: 2.5,
        relaxation_time_s: 0.04,
    };
    let boundaries = [wall(pin, plane(&geometry, 1.0), vec![arm])];
    let resistive = [ResistiveBoundary {
        faces: plane(&geometry, 0.0),
        admittance_m_pa_s: 0.03,
    }];
    let duration = 0.08;
    let owner = model
        .prepare_exact_time_with_walls(
            duration,
            &resistive,
            &boundaries,
            &[],
            &[],
            &[],
            &[],
            ExactMechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let dofs = owner.wall_dofs(0).unwrap();
    assert!(dofs.iter().any(|&dof| dof >= positions.len()));
    assert_eq!(owner.state_count(), 2 * model.ndof() + 3 * dofs.len());
    let amplitude = 0.01;
    let field = p2(&geometry, |x| amplitude * (1.0 + x[1] / lengths[1]).powi(2));
    let displacement = dofs.iter().map(|&dof| field[dof]).collect::<Vec<_>>();
    let velocity = vec![0.0; dofs.len()];
    let memory = displacement
        .iter()
        .map(|&w| 0.4 * det::sqrt(arm.stiffness_per_area) * w)
        .collect::<Vec<_>>();
    let zeros = vec![0.0; model.ndof()];
    let state = owner
        .state_with_wall_memory(
            &zeros,
            &zeros,
            &[displacement.clone()],
            &[velocity.clone()],
            &[vec![memory.clone()]],
        )
        .unwrap();
    // Original fourth-order surface moment int_0^1 (1+y)^4 dy = 31/5.
    // Dropping edge bubbles or averaging this patch cannot pass this oracle.
    let expected = 0.5
        * (pin.stiffness_per_area + 0.36 * arm.stiffness_per_area)
        * area
        * amplitude
        * amplitude
        * 31.0
        / 5.0;
    close(
        owner.energy(&state).unwrap().total_j(),
        expected,
        3e-12 * expected,
    );
    let reference = owner.step(&state, &[], &mut || false).unwrap();
    assert!(reference.work.wall_work_j < 0.0);
    assert!(reference.work.wall_relaxation_dissipation_j > 0.0);
    assert!(reference.work.boundary_dissipation_j > 0.0);
    let half = model
        .prepare_exact_time_with_walls(
            0.5 * duration,
            &resistive,
            &boundaries,
            &[],
            &[],
            &[],
            &[],
            ExactMechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let start = half
        .state_with_wall_memory(
            &zeros,
            &zeros,
            &[displacement.clone()],
            &[velocity.clone()],
            &[vec![memory.clone()]],
        )
        .unwrap();
    let first = half.step(&start, &[], &mut || false).unwrap();
    let second = half.step(&first.state, &[], &mut || false).unwrap();
    equivalent(
        second.state.pressure_pa(),
        reference.state.pressure_pa(),
        3e-9,
    );
    equivalent(
        second.state.pressure_potential_pa_s(),
        reference.state.pressure_potential_pa_s(),
        3e-9,
    );
    equivalent(
        second.state.wall_displacement_m(0).unwrap(),
        reference.state.wall_displacement_m(0).unwrap(),
        3e-9,
    );
    equivalent(
        second.state.wall_normal_velocity_m_s(0).unwrap(),
        reference.state.wall_normal_velocity_m_s(0).unwrap(),
        3e-9,
    );
    equivalent(
        second.state.wall_relaxation_memory(0, 0).unwrap(),
        reference.state.wall_relaxation_memory(0, 0).unwrap(),
        3e-9,
    );
    for (full, part1, part2) in [
        (
            reference.work.boundary_dissipation_j,
            first.work.boundary_dissipation_j,
            second.work.boundary_dissipation_j,
        ),
        (
            reference.work.wall_resistance_dissipation_j,
            first.work.wall_resistance_dissipation_j,
            second.work.wall_resistance_dissipation_j,
        ),
        (
            reference.work.wall_relaxation_dissipation_j,
            first.work.wall_relaxation_dissipation_j,
            second.work.wall_relaxation_dissipation_j,
        ),
        (
            reference.work.wall_work_j,
            first.work.wall_work_j,
            second.work.wall_work_j,
        ),
    ] {
        close(part1 + part2, full, 3e-9 * expected);
    }
    let mut errors = Vec::new();
    for count in [32, 64] {
        let midpoint = model
            .prepare_time_with_walls(
                duration / f64::from(count),
                &resistive,
                &boundaries,
                MechanicalTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        assert_eq!(midpoint.wall_dofs(0).unwrap(), dofs);
        let mut state = midpoint
            .state_with_wall_memory(
                &zeros,
                &zeros,
                &[displacement.clone()],
                &[velocity.clone()],
                &[vec![memory.clone()]],
            )
            .unwrap();
        for _ in 0..count {
            state = midpoint.step(&state, &[], &mut || false).unwrap().state;
        }
        // The norm only compares convergence in the SAME physical coordinates;
        // physical storage and heat are checked separately above.
        let squared_error = [
            (state.pressure_pa(), reference.state.pressure_pa()),
            (
                state.pressure_potential_pa_s(),
                reference.state.pressure_potential_pa_s(),
            ),
            (
                state.wall_displacement_m(0).unwrap(),
                reference.state.wall_displacement_m(0).unwrap(),
            ),
            (
                state.wall_normal_velocity_m_s(0).unwrap(),
                reference.state.wall_normal_velocity_m_s(0).unwrap(),
            ),
            (
                state.wall_relaxation_memory(0, 0).unwrap(),
                reference.state.wall_relaxation_memory(0, 0).unwrap(),
            ),
        ]
        .iter()
        .flat_map(|(actual, expected)| actual.iter().zip(*expected))
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>();
        errors.push(det::sqrt(squared_error));
    }
    assert!(errors[0] > 1e-8);
    assert!(errors[1] < 0.3 * errors[0], "midpoint errors {errors:?}");
}

#[test]
fn prepared_sources_keep_signed_conjugate_work_and_refusals_leave_state_usable() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let other_model = uniform(&geometry);
    let profile = [([0, 1, 2], [0.0, 0.0, 0.0, 1.0, 0.0, 0.0])];
    let source = model
        .quadratic_velocity_source(&profile, &mut || false)
        .unwrap();
    let foreign_source = other_model
        .quadratic_velocity_source(&profile, &mut || false)
        .unwrap();
    let mut calls = 0;
    let volume = model
        .volume_velocity_source(
            4,
            &mut |_, _| {
                calls += 1;
                Ok::<_, HeterogeneousError>(6.0)
            },
            &mut || false,
        )
        .unwrap();
    let sampled = calls;
    let boundaries = [wall(
        WallPin {
            surface_density: 0.6,
            stiffness_per_area: 3.0,
            resistance: 0.2,
        },
        vec![[0, 1, 2]],
        vec![AcousticWallRelaxation {
            stiffness_per_area: 2.5,
            relaxation_time_s: 0.04,
        }],
    )];
    let resistive = [ResistiveBoundary {
        faces: vec![[1, 2, 3]],
        admittance_m_pa_s: 0.02,
    }];
    let owner = model
        .prepare_exact_time_with_walls(
            0.03,
            &resistive,
            &boundaries,
            &[],
            &[&source],
            &[&volume],
            &[],
            ExactMechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let initial = owner.zero_state();
    assert_eq!(owner.port_count(), 1);
    assert_eq!(owner.volume_source_count(), 1);
    assert!(owner.step(&initial, &[0.02], &mut || false).is_err());
    assert!(
        owner
            .step_with_volume_sources(&initial, &[f64::NAN], &[-0.003], &mut || false)
            .is_err()
    );
    assert!(matches!(
        owner.step_with_volume_sources(&initial, &[0.02], &[-0.003], &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    let result = owner
        .step_with_volume_sources(&initial, &[0.02], &[-0.003], &mut || false)
        .unwrap();
    assert!(result.state.pressure_pa().iter().any(|p| *p != 0.0));
    assert!(
        result
            .state
            .wall_relaxation_memory(0, 0)
            .unwrap()
            .iter()
            .any(|z| *z != 0.0)
    );
    let work_scale =
        result.work.boundary_source_work_j.abs() + result.work.volume_source_work_j.abs();
    close(
        result.work.boundary_source_work_j,
        0.02 * result.integrated_port_pressure_pa_m2_s[0],
        3e-12 * work_scale,
    );
    close(
        result.work.volume_source_work_j,
        -0.003 * result.integrated_volume_source_pressure_pa_s[0],
        3e-12 * work_scale,
    );
    close(
        result.work.input_work_j,
        result.work.boundary_source_work_j + result.work.volume_source_work_j,
        3e-12 * work_scale,
    );
    close(result.work.balance_residual_j(), 0.0, 3e-9 * work_scale);
    assert_eq!(
        calls, sampled,
        "warm wall steps must not resample material or source laws"
    );
    assert!(initial.pressure_pa().iter().all(|p| *p == 0.0));
    assert!(
        model
            .prepare_exact_time_with_walls(
                0.03,
                &resistive,
                &boundaries,
                &[],
                &[&foreign_source],
                &[],
                &[],
                ExactMechanicalTimeLimits::default(),
                &mut || false,
            )
            .is_err()
    );
    assert!(matches!(
        model.prepare_exact_time_with_walls(
            0.03,
            &resistive,
            &boundaries,
            &[],
            &[&source],
            &[&volume],
            &[],
            ExactMechanicalTimeLimits {
                max_states: owner.state_count() - 1,
                ..ExactMechanicalTimeLimits::default()
            },
            &mut || false,
        ),
        Err(HeterogeneousError::Limit { .. })
    ));
    assert!(
        model
            .prepare_exact_time_with_walls(
                0.03,
                &[],
                &boundaries,
                &[],
                &[],
                &[],
                &[[0, 1, 2]],
                ExactMechanicalTimeLimits::default(),
                &mut || false,
            )
            .is_err()
    );
    let mut duplicate = boundaries.to_vec();
    duplicate[0].faces.push([2, 1, 0]);
    assert!(
        model
            .prepare_exact_time_with_walls(
                0.03,
                &[],
                &duplicate,
                &[],
                &[],
                &[],
                &[],
                ExactMechanicalTimeLimits::default(),
                &mut || false,
            )
            .is_err()
    );
    let foreign_owner = model
        .prepare_exact_time_with_walls(
            0.015,
            &resistive,
            &boundaries,
            &[],
            &[&source],
            &[&volume],
            &[],
            ExactMechanicalTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    assert!(
        foreign_owner
            .step_with_volume_sources(&initial, &[0.02], &[-0.003], &mut || false)
            .is_err()
    );
}
