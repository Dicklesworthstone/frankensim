use super::super::AcousticFeError;
use super::super::tests::{MEDIUM, brick, p2, plane, uniform};
use super::*;
use crate::acoustics::helmholtz::{HelmholtzFe, Medium};
use fs_modal::SliceOptions;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn projection(model: &HeterogeneousHelmholtzFe<'_, '_>, mode: &[f64], field: &[f64]) -> f64 {
    let mut product = vec![0.0; field.len()];
    model.mass().spmv(field, &mut product);
    mode.iter().zip(product).map(|(a, b)| a * b).sum()
}

#[test]
fn polynomial_compression_integrates_actual_port_pressure_even_for_weak_sources() {
    let lengths = [1.0, 0.4, 0.2];
    let volume = lengths.iter().product::<f64>();
    let (complex, positions) = brick(1, lengths);
    for medium in [
        MEDIUM,
        Medium {
            density_kg_m3: 1.2,
            sound_speed_m_s: 343.0,
        },
        Medium {
            density_kg_m3: 1000.0,
            sound_speed_m_s: 1480.0,
        },
    ] {
        let geometry = HelmholtzFe::new(&complex, &positions, 2, medium).unwrap();
        let model = HeterogeneousHelmholtzFe::new(
            &geometry,
            4,
            &mut |_, _| Ok::<_, AcousticFeError>(medium),
            &mut || false,
        )
        .unwrap();
        let port = ExactVelocityPort {
            profile: plane(&geometry, lengths[0])
                .into_iter()
                .map(|face| (face, [1.0; 3]))
                .collect(),
        };
        let c = medium.sound_speed_m_s;
        // Identical c*dt/L tests SI air/water conditioning at the same phase.
        let dt = 0.05 * lengths[0] / c;
        let owner = model
            .prepare_exact_time(dt, &[], &[port], ExactTimeLimits::default(), &mut || false)
            .unwrap();
        for a in [0.02, 2e-10] {
            // psi=a*x²/2+a*c²*t²/2, p=a*c²*t is an exact P2 solution.
            // v_in=a*L/rho on the right; all other walls are rigid.
            let mut state = owner
                .state(
                    &p2(&geometry, |x| 0.5 * a * x[0] * x[0]),
                    &vec![0.0; model.ndof()],
                )
                .unwrap();
            let kinetic = a * a * volume * lengths[0] * lengths[0] / (6.0 * medium.density_kg_m3);
            let amplitude = a * lengths[0] / medium.density_kg_m3;
            let pressure_scale = a.abs() * c * lengths[0];
            let potential_scale = a.abs() * lengths[0] * lengths[0];
            let mut accumulated_work = 0.0;
            for step in 1..=4 {
                let result = owner.step(&state, &[amplitude], &mut || false).unwrap();
                let t = f64::from(step) * dt;
                let p = a * c * c * t;
                let expected_p = p2(&geometry, |_| p);
                let expected_psi = p2(&geometry, |x| 0.5 * a * (x[0] * x[0] + c * c * t * t));
                for (actual, expected) in result.state.pressure_pa().iter().zip(expected_p) {
                    close(*actual, expected, 2e-10 * pressure_scale);
                }
                for (actual, expected) in result
                    .state
                    .pressure_potential_pa_s()
                    .iter()
                    .zip(expected_psi)
                {
                    close(*actual, expected, 2e-10 * potential_scale);
                }
                let integral = lengths[1] * lengths[2] * a * c * c * dt * (t - 0.5 * dt);
                close(
                    result.integrated_port_pressure_pa_m2_s[0],
                    integral,
                    2e-9 * integral.abs(),
                );
                close(
                    result.work.input_work_j,
                    amplitude * integral,
                    2e-9 * (amplitude * integral).abs(),
                );
                assert!(result.work.input_work_j > 0.0);
                close(result.work.boundary_dissipation_j, 0.0, 1e-12 * kinetic);
                close(result.work.balance_residual_j(), 0.0, 2e-9 * kinetic);
                accumulated_work += result.work.input_work_j;
                state = result.state;
                let energy = owner.energy(&state).unwrap();
                let compression = volume * p * p / (2.0 * medium.density_kg_m3 * c * c);
                close(energy.kinetic_j, kinetic, 2e-9 * kinetic);
                close(energy.compression_j, compression, 2e-9 * compression);
                close(accumulated_work, compression, 2e-9 * compression);
            }
        }
    }
}

#[test]
fn complete_cavity_mode_has_no_clock_warping_across_a_material_interface() {
    let (complex, positions) = brick(2, [1.0, 0.2, 0.2]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = HeterogeneousHelmholtzFe::new(
        &geometry,
        4,
        &mut |_, x| {
            Ok::<_, AcousticFeError>(if x[0] < 0.5 {
                MEDIUM
            } else {
                Medium {
                    density_kg_m3: 2.6,
                    sound_speed_m_s: 3.0,
                }
            })
        },
        &mut || false,
    )
    .unwrap();
    let modes = model
        .modes(
            (0.5, 1.7),
            &SliceOptions {
                mass_diagonal_equilibration: true,
                ..SliceOptions::default()
            },
            &mut || false,
        )
        .unwrap();
    assert_eq!(modes.shapes.len(), 1);
    let shape = &modes.shapes[0];
    let omega = det::sqrt(modes.omega_squared_rad2_s2[0]);
    let duration = 0.35;
    for count in [1, 2, 4] {
        let owner = model
            .prepare_exact_time(
                duration / f64::from(count),
                &[],
                &[],
                ExactTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let mut state = owner.state(&vec![0.0; model.ndof()], shape).unwrap();
        let initial = owner.energy(&state).unwrap().total_j();
        for _ in 0..count {
            state = owner.step(&state, &[], &mut || false).unwrap().state;
        }
        // The eigenpair is a reference for the ORIGINAL semidiscrete operator,
        // not a claim that its coarse spatial eigenfrequency is continuum-exact.
        close(
            projection(&model, shape, state.pressure_pa()),
            det::cos(omega * duration),
            2e-8,
        );
        close(
            projection(&model, shape, state.pressure_potential_pa_s()),
            det::sin(omega * duration) / omega,
            2e-8,
        );
        close(
            owner.energy(&state).unwrap().total_j(),
            initial,
            2e-8 * initial,
        );
    }
}

#[test]
fn spatial_refinement_converges_to_the_analytic_rigid_cavity_phase() {
    let duration = 0.21;
    // Continuum axial Neumann mode: omega=pi*c/L, independent of FE matrices.
    let continuum = det::cos(core::f64::consts::PI * MEDIUM.sound_speed_m_s * duration);
    let mut errors = Vec::new();
    for nx in [2, 4, 8] {
        let (complex, positions) = brick(nx, [1.0, 0.05, 0.05]);
        let geometry = HelmholtzFe::new(&complex, &positions, 1, MEDIUM).unwrap();
        let model = uniform(&geometry);
        let modes = model
            .modes(
                (0.7, 1.4),
                &SliceOptions {
                    mass_diagonal_equilibration: true,
                    ..SliceOptions::default()
                },
                &mut || false,
            )
            .unwrap();
        assert_eq!(modes.shapes.len(), 1);
        let shape = &modes.shapes[0];
        let owner = model
            .prepare_exact_time(duration, &[], &[], ExactTimeLimits::default(), &mut || {
                false
            })
            .unwrap();
        let initial = owner.state(&vec![0.0; model.ndof()], shape).unwrap();
        let result = owner.step(&initial, &[], &mut || false).unwrap();
        errors.push((projection(&model, shape, result.state.pressure_pa()) - continuum).abs());
    }
    assert!(errors[1] < 0.7 * errors[0], "spatial errors: {errors:?}");
    assert!(errors[2] < 0.7 * errors[1], "spatial errors: {errors:?}");
    assert!(errors[2] < 0.03, "spatial errors: {errors:?}");
}

// Independent exact P1 triangle integral: area/12 * (sum(p_i²)+(sum p_i)²).
// This deliberately uses original geometry/pressure, never the loss matrix or
// a subtraction of endpoint energies returned by the production integrator.
fn wall_power(geometry: &HelmholtzFe<'_>, faces: &[[u32; 3]], p: &[f64], admittance: f64) -> f64 {
    faces
        .iter()
        .map(|face| {
            let a = geometry.positions()[face[0] as usize];
            let b = geometry.positions()[face[1] as usize];
            let c = geometry.positions()[face[2] as usize];
            let u: [f64; 3] = core::array::from_fn(|i| b[i] - a[i]);
            let v: [f64; 3] = core::array::from_fn(|i| c[i] - a[i]);
            let cross = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            let area = 0.5 * det::sqrt(cross.iter().map(|x| x * x).sum());
            let pressure = face.map(|i| p[i as usize]);
            let sum = pressure.iter().sum::<f64>();
            admittance * area * (pressure.iter().map(|x| x * x).sum::<f64>() + sum * sum) / 12.0
        })
        .sum()
}

#[test]
fn boundary_heat_matches_independent_time_quadrature_of_actual_face_pressure() {
    let (complex, positions) = brick(1, [1.0, 0.4, 0.2]);
    let geometry = HelmholtzFe::new(&complex, &positions, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let faces = plane(&geometry, 1.0);
    let admittance = 0.25;
    let boundary = ResistiveBoundary {
        faces: faces.clone(),
        admittance_m_pa_s: admittance,
    };
    let duration = 0.08;
    let coarse = model
        .prepare_exact_time(
            duration,
            &[boundary.clone()],
            &[],
            ExactTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let fine = model
        .prepare_exact_time(
            duration / 32.0,
            &[boundary],
            &[],
            ExactTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let initial_pressure: Vec<f64> = positions
        .iter()
        .map(|x| 1.0 + 0.4 * x[0] - 0.2 * x[1])
        .collect();
    let initial = coarse
        .state(&vec![0.0; model.ndof()], &initial_pressure)
        .unwrap();
    let full = coarse.step(&initial, &[], &mut || false).unwrap();
    let mut state = initial.clone();
    let mut simpson = wall_power(&geometry, &faces, state.pressure_pa(), admittance);
    for step in 1..=32 {
        state = fine.step(&state, &[], &mut || false).unwrap().state;
        let weight = if step == 32 {
            1.0
        } else if step % 2 == 0 {
            2.0
        } else {
            4.0
        };
        simpson += weight * wall_power(&geometry, &faces, state.pressure_pa(), admittance);
    }
    simpson *= duration / (3.0 * 32.0);
    assert!(full.work.boundary_dissipation_j > 0.0);
    close(full.work.boundary_dissipation_j, simpson, 3e-6 * simpson);
    close(full.work.input_work_j, 0.0, 1e-15);
    for (actual, expected) in full.state.pressure_pa().iter().zip(state.pressure_pa()) {
        close(*actual, *expected, 1e-9);
    }
    close(
        full.work.balance_residual_j(),
        0.0,
        1e-9 * full.work.initial_energy_j,
    );
}

#[test]
fn admission_state_binding_and_cancellation_do_not_publish_partial_fields() {
    let (complex, positions) = brick(1, [1.0, 0.4, 0.2]);
    let geometry = HelmholtzFe::new(&complex, &positions, 1, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let another_model = uniform(&geometry);
    let limits = ExactTimeLimits::default();
    assert!(
        model
            .prepare_exact_time(
                0.01,
                &[],
                &[],
                ExactTimeLimits {
                    max_dofs: model.ndof() - 1,
                    ..limits
                },
                &mut || false
            )
            .is_err()
    );
    let mut memory = limits;
    memory.preparation.max_dense_bytes = 1;
    assert!(
        model
            .prepare_exact_time(0.01, &[], &[], memory, &mut || false)
            .is_err()
    );
    assert!(matches!(
        model.prepare_exact_time(0.01, &[], &[], limits, &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    assert!(
        model
            .prepare_exact_time(f64::NAN, &[], &[], limits, &mut || false)
            .is_err()
    );
    assert!(
        model
            .prepare_exact_time(
                0.01,
                &[],
                &[ExactVelocityPort {
                    profile: Vec::new()
                }],
                limits,
                &mut || false
            )
            .is_err()
    );
    let owner = model
        .prepare_exact_time(0.01, &[], &[], limits, &mut || false)
        .unwrap();
    let foreign = another_model
        .prepare_exact_time(0.01, &[], &[], limits, &mut || false)
        .unwrap();
    assert!(
        owner
            .step(&foreign.zero_state(), &[], &mut || false)
            .is_err()
    );
    assert!(owner.state(&[0.0], &[0.0]).is_err());
    assert!(
        owner
            .state(&vec![0.0; model.ndof()], &vec![f64::NAN; model.ndof()])
            .is_err()
    );
    let state = owner
        .state(&vec![0.0; model.ndof()], &vec![1.0; model.ndof()])
        .unwrap();
    assert!(owner.step(&state, &[1.0], &mut || false).is_err());
    let mut polls = 0;
    let result = owner.step(&state, &[], &mut || {
        polls += 1;
        polls == 4
    });
    assert!(matches!(result, Err(HeterogeneousError::Cancelled)));
    assert!(state.pressure_pa().iter().all(|&p| p == 1.0));
    assert!(state.pressure_potential_pa_s().iter().all(|&p| p == 0.0));
}
