use super::super::tests::{MEDIUM, brick, p2, plane};
use super::super::{AcousticFeError, HeterogeneousTimeLimits};
use super::*;
use crate::acoustics::helmholtz::pressure_release::EmptyFreePressure;
use crate::acoustics::helmholtz::{HelmholtzFe, Medium};
use crate::fixtures::single_tet;
use fs_modal::SliceOptions;

fn constant<'f, 'c>(
    geometry: &'f HelmholtzFe<'c>,
    medium: Medium,
) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        geometry,
        geometry.order() + 2,
        &mut |_, _| Ok::<_, AcousticFeError>(medium),
        &mut || false,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn same_bits(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (&a, &b) in actual.iter().zip(expected) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
}

#[test]
fn pressure_release_exact_interface_mode_preserves_phase_and_original_energy() {
    let (complex, positions) = brick(2, [1.0, 0.2, 0.2]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = HeterogeneousHelmholtzFe::new(
        &geometry,
        4,
        &mut |_, point| {
            Ok::<_, AcousticFeError>(Medium {
                density_kg_m3: if point[0] < 0.5 { 1.3 } else { 2.6 },
                sound_speed_m_s: 2.0,
            })
        },
        &mut || false,
    )
    .unwrap();
    let release = plane(&geometry, 1.0);
    let modes = model
        .modes_with_pressure_release(
            (0.1, 1.0),
            &release,
            &SliceOptions {
                mass_diagonal_equilibration: true,
                ..SliceOptions::default()
            },
            &mut || false,
        )
        .unwrap()
        .eigenpairs;
    assert_eq!(modes.shapes.len(), 1);
    let shape = &modes.shapes[0];
    let maximum = shape.iter().map(|value| value.abs()).fold(0.0, f64::max);
    let omega = det::sqrt(modes.omega_squared_rad2_s2[0]);
    let duration = 0.35;
    // G1: this is the original heterogeneous semidiscrete eigenpair, with
    // every free coordinate retained. Spatial frequency error is separate.
    // Exact propagation must not acquire midpoint's temporal frequency warp.
    for count in [1, 2, 4] {
        let owner = model
            .prepare_exact_time_with_pressure_release(
                duration / f64::from(count),
                &[],
                &[],
                &[],
                &release,
                ExactTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let mut state = owner.state(&vec![0.0; model.ndof()], shape).unwrap();
        let initial = owner.energy(&state).unwrap().total_j();
        for _ in 0..count {
            let next = owner.step(&state, &[], &mut || false).unwrap();
            close(next.work.input_work_j, 0.0, 1e-14 * initial);
            close(next.work.boundary_dissipation_j, 0.0, 1e-14 * initial);
            state = next.state;
        }
        for (dof, &mode) in shape.iter().enumerate() {
            close(
                state.pressure_pa()[dof],
                mode * det::cos(omega * duration),
                2e-8 * maximum,
            );
            close(
                state.pressure_potential_pa_s()[dof],
                mode * det::sin(omega * duration) / omega,
                2e-8 * maximum / omega,
            );
            if owner.pressure_release.constrained[dof] {
                assert_eq!(state.pressure_pa()[dof], 0.0);
                assert_eq!(state.pressure_potential_pa_s()[dof], 0.0);
            }
        }
        close(
            owner.energy(&state).unwrap().total_j(),
            initial,
            2e-8 * initial,
        );
    }
}

#[test]
fn pressure_release_exact_held_sources_and_robin_heat_match_analytic_tetrahedron() {
    let (complex, positions) = single_tet();
    for medium in [
        MEDIUM,
        Medium {
            density_kg_m3: 1.2,
            sound_speed_m_s: 343.0,
        },
    ] {
        let geometry = HelmholtzFe::new(&complex, &positions, 1, medium).unwrap();
        let model = constant(&geometry, medium);
        let mut source_calls = 0;
        let volume_source = model
            .volume_velocity_source(
                3,
                &mut |_, _| {
                    source_calls += 1;
                    Ok::<_, HeterogeneousError>(6.0)
                },
                &mut || false,
            )
            .unwrap();
        let prepared_calls = source_calls;
        let rho = medium.density_kg_m3;
        let speed = medium.sound_speed_m_s;
        let admittance = 0.2 / (rho * speed);
        let duration = 0.6 / speed;
        let owner = model
            .prepare_exact_time_with_pressure_release(
                duration,
                &[ResistiveBoundary {
                    faces: vec![[0, 1, 2]],
                    admittance_m_pa_s: admittance,
                }],
                &[ExactVelocityPort {
                    profile: vec![([0, 1, 2], [1.0; 3])],
                }],
                &[&volume_source],
                &[[3, 2, 1]],
                ExactTimeLimits {
                    max_dofs: 1,
                    ..ExactTimeLimits::default()
                },
                &mut || false,
            )
            .unwrap();
        assert_eq!(owner.port_count(), 1);
        assert_eq!(owner.volume_source_count(), 1);
        assert_eq!(owner.pressure_release.free_dofs, [0]);
        // G1: the complete P1 tetrahedron after the declared Dirichlet face
        // has one free basis function phi=1-x-y-z. Independently integrated:
        // V=1/6, int(phi²)=V/10, int(|grad phi|²)=3V,
        // inlet area A=1/2, int_face(phi)=A/3, int_face(phi²)=A/6,
        // int_volume(phi/V)=1/4. Inverting full M then slicing would give
        // 16/(beta V), whereas the correct inverse of M_ff is 10/(beta V).
        let mass = 1.0 / (60.0 * rho * speed * speed);
        let stiffness = 1.0 / (2.0 * rho);
        let resistance = admittance / 12.0;
        let boundary_load = 1.0 / 6.0;
        let volume_load = 1.0 / 4.0;
        let alpha = resistance / (2.0 * mass);
        let damped_omega = det::sqrt(stiffness / mass - alpha * alpha);
        let phase = damped_omega * duration;
        let decay = det::exp(-alpha * duration);
        // Closed-form integral of exp(-2 alpha t) sin²(omega_d t).
        // This observes physical Robin heat independently of energy balance.
        let square_integral = (1.0 - decay * decay) / (4.0 * alpha)
            - (alpha
                + decay
                    * decay
                    * (-alpha * det::cos(2.0 * phase) + damped_omega * det::sin(2.0 * phase)))
                / (4.0 * (alpha * alpha + damped_omega * damped_omega));
        for scale in [1.0, 1e-8] {
            let inward_velocity = 0.02 * scale;
            let volume_flow = -0.008 * scale;
            let forcing = boundary_load * inward_velocity + volume_load * volume_flow;
            let potential = forcing / stiffness
                * (1.0 - decay * (det::cos(phase) + alpha / damped_omega * det::sin(phase)));
            let pressure_amplitude = forcing / (mass * damped_omega);
            let pressure = pressure_amplitude * decay * det::sin(phase);
            let heat = resistance * pressure_amplitude * pressure_amplitude * square_integral;
            let energy = 0.5 * (mass * pressure * pressure + stiffness * potential * potential);
            let input = forcing * potential;
            let result = owner
                .step_with_volume_sources(
                    &owner.zero_state(),
                    &[inward_velocity],
                    &[volume_flow],
                    &mut || false,
                )
                .unwrap();
            close(
                result.state.pressure_pa()[0],
                pressure,
                3e-9 * pressure_amplitude.abs(),
            );
            close(
                result.state.pressure_potential_pa_s()[0],
                potential,
                3e-9 * potential.abs(),
            );
            for dof in 1..4 {
                assert_eq!(result.state.pressure_pa()[dof], 0.0);
                assert_eq!(result.state.pressure_potential_pa_s()[dof], 0.0);
            }
            close(result.work.final_energy_j, energy, 3e-9 * energy);
            close(result.work.input_work_j, input, 3e-9 * input);
            close(result.work.boundary_dissipation_j, heat, 3e-9 * heat);
            close(result.work.balance_residual_j(), 0.0, 3e-9 * input);
            close(
                result.integrated_port_pressure_pa_m2_s[0],
                boundary_load * potential,
                3e-9 * (boundary_load * potential).abs(),
            );
            close(
                result.integrated_volume_source_pressure_pa_s[0],
                volume_load * potential,
                3e-9 * (volume_load * potential).abs(),
            );
            let extraction = volume_flow * volume_load * potential;
            close(
                result.work.volume_source_work_j,
                extraction,
                3e-9 * extraction.abs(),
            );
            assert!(result.work.volume_source_work_j < 0.0);
            assert!(result.work.boundary_dissipation_j > 0.0);
            assert_eq!(source_calls, prepared_calls);
        }
    }
}

#[test]
fn pressure_release_exact_caps_free_states_and_refuses_conflicting_laws_without_mutation() {
    let (complex, positions) = brick(1, [1.0; 3]);
    let geometry = HelmholtzFe::new(&complex, &positions, 3, MEDIUM).unwrap();
    let model = constant(&geometry, MEDIUM);
    let release = geometry.boundary_faces();
    assert!(model.ndof() > ExactTimeLimits::default().max_dofs);
    let constraints = geometry
        .pressure_release_constraints(&release, EmptyFreePressure::Refuse, &mut || false)
        .unwrap();
    let free = constraints.free_dofs.len();
    let limits = ExactTimeLimits {
        max_dofs: free,
        ..ExactTimeLimits::default()
    };
    let owner = model
        .prepare_exact_time_with_pressure_release(
            0.01,
            &[],
            &[],
            &[],
            &release,
            limits,
            &mut || false,
        )
        .unwrap();
    assert!(free < model.ndof());
    assert!(
        model
            .prepare_exact_time_with_pressure_release(
                0.01,
                &[],
                &[],
                &[],
                &release,
                ExactTimeLimits {
                    max_dofs: free - 1,
                    ..limits
                },
                &mut || false,
            )
            .is_err()
    );
    // G0: empty faces retain the original cheap augmented-state refusal,
    // even when a caller raises its own pressure-DOF allowance.
    let mut admission_polls = 0;
    assert!(matches!(
        model.prepare_exact_time(
            0.01,
            &[],
            &[],
            ExactTimeLimits {
                max_dofs: usize::MAX,
                ..limits
            },
            &mut || {
                admission_polls += 1;
                false
            },
        ),
        Err(HeterogeneousError::Limit {
            what: "exact transient augmented states",
            ..
        })
    ));
    assert_eq!(admission_polls, 1);
    let face = release[0];
    let reversed = [face[2], face[1], face[0]];
    assert!(
        model
            .prepare_exact_time_with_pressure_release(
                0.01,
                &[ResistiveBoundary {
                    faces: vec![reversed],
                    admittance_m_pa_s: 0.0
                }],
                &[],
                &[],
                &release,
                limits,
                &mut || false,
            )
            .is_err()
    );
    assert!(
        model
            .prepare_exact_time_with_pressure_release(
                0.01,
                &[],
                &[ExactVelocityPort {
                    profile: vec![(reversed, [0.0; 3])]
                }],
                &[],
                &release,
                limits,
                &mut || false,
            )
            .is_err()
    );
    assert!(matches!(
        model.prepare_exact_time_with_pressure_release(
            0.01,
            &[],
            &[],
            &[],
            &release,
            limits,
            &mut || true,
        ),
        Err(HeterogeneousError::Cancelled)
    ));
    let zero = vec![0.0; model.ndof()];
    let face_dof = constraints
        .constrained
        .iter()
        .enumerate()
        .find(|(dof, fixed)| *dof >= geometry.space.face_off && **fixed)
        .unwrap()
        .0;
    let mut incompatible = zero.clone();
    incompatible[face_dof] = f64::MIN_POSITIVE;
    assert!(owner.state(&incompatible, &zero).is_err());
    assert!(owner.state(&zero, &incompatible).is_err());
    let rigid = model
        .prepare_time(0.01, &[], HeterogeneousTimeLimits::default(), &mut || false)
        .unwrap();
    let imported = rigid.state(&incompatible, &zero).unwrap();
    assert!(owner.energy(&imported).is_err());
    assert!(owner.step(&imported, &[], &mut || false).is_err());
    assert_eq!(
        imported.pressure_potential_pa_s()[face_dof],
        f64::MIN_POSITIVE
    );
    let mut pressure = zero.clone();
    for (index, &dof) in constraints.free_dofs.iter().enumerate() {
        pressure[dof] = 0.01 * (index + 1) as f64;
    }
    let initial = owner.state(&zero, &pressure).unwrap();
    let mut polls = 0;
    assert!(matches!(
        owner.step(&initial, &[], &mut || {
            polls += 1;
            polls == 4
        }),
        Err(HeterogeneousError::Cancelled)
    ));
    same_bits(initial.pressure_pa(), &pressure);
    same_bits(initial.pressure_potential_pa_s(), &zero);
    let next = owner.step(&initial, &[], &mut || false).unwrap();
    assert_eq!(next.state.pressure_pa().len(), model.ndof());
    for (dof, &fixed) in constraints.constrained.iter().enumerate() {
        if fixed {
            assert_eq!(next.state.pressure_pa()[dof], 0.0);
            assert_eq!(next.state.pressure_potential_pa_s()[dof], 0.0);
        }
    }
}

#[test]
fn pressure_release_exact_empty_faces_preserve_full_fields_and_signed_source_work_bits() {
    let (complex, positions) = single_tet();
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = constant(&geometry, MEDIUM);
    let source = model
        .volume_velocity_source(4, &mut |_, _| Ok::<_, HeterogeneousError>(6.0), &mut || {
            false
        })
        .unwrap();
    let boundaries = [ResistiveBoundary {
        faces: vec![[0, 1, 2]],
        admittance_m_pa_s: 0.1,
    }];
    let ports = [ExactVelocityPort {
        profile: vec![([0, 1, 2], [0.5, 1.0, 1.5])],
    }];
    let legacy = model
        .prepare_exact_time_with_volume_sources(
            0.03,
            &boundaries,
            &ports,
            &[&source],
            ExactTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let explicit = model
        .prepare_exact_time_with_pressure_release(
            0.03,
            &boundaries,
            &ports,
            &[&source],
            &[],
            ExactTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let initial = legacy
        .state(
            &p2(&geometry, |point| {
                0.3 + 0.02 * point[0] * point[1] + 0.03 * point[2].powi(2)
            }),
            &p2(&geometry, |point| 0.1 + 0.02 * point[0] - 0.03 * point[1]),
        )
        .unwrap();
    let old = legacy
        .step_with_volume_sources(&initial, &[0.02], &[-0.003], &mut || false)
        .unwrap();
    let new = explicit
        .step_with_volume_sources(&initial, &[0.02], &[-0.003], &mut || false)
        .unwrap();
    same_bits(new.state.pressure_pa(), old.state.pressure_pa());
    same_bits(
        new.state.pressure_potential_pa_s(),
        old.state.pressure_potential_pa_s(),
    );
    same_bits(
        &new.integrated_port_pressure_pa_m2_s,
        &old.integrated_port_pressure_pa_m2_s,
    );
    same_bits(
        &new.integrated_volume_source_pressure_pa_s,
        &old.integrated_volume_source_pressure_pa_s,
    );
    let work = |w: ExactTimeWork| {
        [
            w.initial_energy_j,
            w.final_energy_j,
            w.delta_energy_j,
            w.input_work_j,
            w.volume_source_work_j,
            w.boundary_dissipation_j,
        ]
    };
    same_bits(&work(new.work), &work(old.work));
}
