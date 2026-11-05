//! Native acoustic state, phase and work checks on exact pressure-release traces.
use super::tests::{MEDIUM, brick, p2, plane, uniform};
use super::*;
use crate::acoustics::helmholtz::HelmholtzFe;
use fs_modal::SliceOptions;

fn velocity(faces: &[[u32; 3]], value: f64) -> Vec<RealVelocityTrace> {
    faces.iter().map(|&face| (face, [value; 3])).collect()
}

#[test]
fn pressure_release_quarter_wave_phase_converges_while_energy_is_conserved() {
    let (complex, positions) = brick(8, [1.0, 0.2, 0.2]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let release = plane(&geometry, 1.0);
    let modal = model
        .modes_with_pressure_release((0.1, 0.8), &release, &SliceOptions::default(), &mut || {
            false
        })
        .unwrap();
    assert_eq!(modal.eigenpairs.frequency_hz.len(), 1);
    // The native spatial mode converges to f=c/(4L)=0.5 Hz. Using its
    // certified frequency separates temporal error from the spatial error.
    assert!((modal.eigenpairs.frequency_hz[0] / 0.5 - 1.0).abs() < 2e-4);
    let omega = modal.eigenpairs.omega_squared_rad2_s2[0].sqrt();
    let mode = &modal.eigenpairs.shapes[0];
    let amplitude = mode[0];
    assert!(amplitude.abs() > 0.0);
    let initial_pressure: Vec<_> = mode.iter().map(|value| value / amplitude).collect();
    let duration = 0.7;
    let mut errors = Vec::new();
    for steps in [16, 32, 64] {
        let prepared = model
            .prepare_time_with_pressure_release(
                duration / steps as f64,
                &[],
                &release,
                HeterogeneousTimeLimits::default(),
                &mut || false,
            )
            .unwrap();
        let mut state = prepared
            .state(&vec![0.0; geometry.ndof()], &initial_pressure)
            .unwrap();
        let energy = prepared.energy(&state).unwrap().total_j();
        assert!(energy > 0.0);
        for _ in 0..steps {
            let step = prepared.step(&state, &[], &mut || false).unwrap();
            assert_eq!(step.work.input_work_j, 0.0);
            assert_eq!(step.work.boundary_dissipation_j, 0.0);
            assert!((step.work.final_energy_j / energy - 1.0).abs() < 2e-9);
            assert!(step.work.backward_error <= 5e-12);
            for (dof, &fixed) in prepared.pressure_release.constrained.iter().enumerate() {
                if fixed {
                    assert_eq!(step.state.pressure_pa()[dof], 0.0);
                    assert_eq!(step.state.pressure_potential_pa_s()[dof], 0.0);
                }
            }
            state = step.state;
        }
        let (sin, cos) = (omega * duration).sin_cos();
        let mut error = 0.0_f64;
        for (dof, &shape) in initial_pressure.iter().enumerate() {
            error = error.max((state.pressure_pa()[dof] - shape * cos).abs());
            error = error.max((state.pressure_potential_pa_s()[dof] * omega - shape * sin).abs());
        }
        errors.push(error);
    }
    assert!(
        errors[1] < errors[0] / 3.5 && errors[2] < errors[1] / 3.5,
        "{errors:?}"
    );
    assert!(errors[2] < 3e-4, "{errors:?}");
}

#[test]
fn pressure_release_volume_injection_and_resistive_work_match_exact_p2_flow() {
    let (complex, positions) = brick(2, [1.0, 0.2, 0.2]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let release = plane(&geometry, 1.0);
    let inlet = plane(&geometry, 0.0);
    let h = 0.03;
    let admittance = 0.03;
    let prepared = model
        .prepare_time_with_pressure_release(
            h,
            &[ResistiveBoundary {
                faces: inlet.clone(),
                admittance_m_pa_s: admittance,
            }],
            &release,
            HeterogeneousTimeLimits::default(),
            &mut || false,
        )
        .unwrap();
    let volume = 0.04;
    let mut source_queries = 0;
    let source = model
        .volume_velocity_source(
            4,
            &mut |_, _| {
                source_queries += 1;
                Ok::<_, HeterogeneousError>(1.0 / volume)
            },
            &mut || false,
        )
        .unwrap();
    let original_queries = source_queries;
    // psi=t*(1-x²), p=1-x², v_x=2*t*x/rho solve the original PDE
    // exactly with uniform volume injection s=2*t/rho. At x=0, the
    // imposed inward velocity Y*p cancels the explicit passive admittance.
    // At x=1, p=psi=0; its nonzero outgoing flow performs zero work.
    let shape = p2(&geometry, |point| 1.0 - point[0] * point[0]);
    let mut state = prepared.state(&vec![0.0; geometry.ndof()], &shape).unwrap();
    let compression = 4.0 * volume / (15.0 * MEDIUM.density_kg_m3 * MEDIUM.sound_speed_m_s.powi(2));
    for index in 0..8 {
        let middle = (index as f64 + 0.5) * h;
        let flow = 2.0 * middle * volume / MEDIUM.density_kg_m3;
        let step = prepared
            .step_with_volume_sources(
                &state,
                &velocity(&inlet, admittance),
                &[(&source, flow)],
                &mut || false,
            )
            .unwrap();
        let time = (index + 1) as f64 * h;
        for (dof, &expected) in shape.iter().enumerate() {
            assert!((step.state.pressure_pa()[dof] - expected).abs() < 2e-10);
            assert!((step.state.pressure_potential_pa_s()[dof] - time * expected).abs() < 2e-10);
        }
        let kinetic = 2.0 * volume * time * time / (3.0 * MEDIUM.density_kg_m3);
        let energy = prepared.energy(&step.state).unwrap();
        assert!((energy.compression_j / compression - 1.0).abs() < 2e-10);
        assert!((energy.kinetic_j / kinetic - 1.0).abs() < 2e-10);
        let boundary_heat = h * admittance * 0.04;
        let volume_work = h * flow * (2.0 / 3.0);
        assert!((step.work.boundary_dissipation_j / boundary_heat - 1.0).abs() < 2e-10);
        assert!((step.work.volume_source_work_j / volume_work - 1.0).abs() < 2e-10);
        assert!((step.work.input_work_j / (boundary_heat + volume_work) - 1.0).abs() < 2e-10);
        assert!(step.work.balance_residual_j().abs() < 2e-10 * step.work.input_work_j);
        state = step.state;
    }
    assert_eq!(
        source_queries, original_queries,
        "no material/source resampling during evolution"
    );
}

#[test]
fn pressure_release_time_refuses_incompatible_states_sources_and_limits_without_mutation() {
    let (complex, positions) = brick(2, [1.0, 0.2, 0.2]);
    let geometry = HelmholtzFe::new(&complex, &positions, 2, MEDIUM).unwrap();
    let model = uniform(&geometry);
    let release = plane(&geometry, 1.0);
    let limits = HeterogeneousTimeLimits::default();
    let prepared = model
        .prepare_time_with_pressure_release(0.01, &[], &release, limits, &mut || false)
        .unwrap();
    let rigid = model
        .prepare_time(0.01, &[], limits, &mut || false)
        .unwrap();
    let empty = model
        .prepare_time_with_pressure_release(0.01, &[], &[], limits, &mut || false)
        .unwrap();
    let zero = prepared.zero_state();
    let free = prepared.pressure_release.free_dofs.len();
    assert!(free > 1 && free < geometry.ndof());
    model
        .prepare_time_with_pressure_release(
            0.01,
            &[],
            &release,
            HeterogeneousTimeLimits {
                max_dofs: free,
                ..limits
            },
            &mut || false,
        )
        .unwrap();
    assert!(
        model
            .prepare_time_with_pressure_release(
                0.01,
                &[],
                &release,
                HeterogeneousTimeLimits {
                    max_dofs: free - 1,
                    ..limits
                },
                &mut || false,
            )
            .is_err()
    );
    let edge = prepared
        .pressure_release
        .constrained
        .iter()
        .enumerate()
        .find_map(|(dof, &fixed)| (fixed && dof >= positions.len()).then_some(dof))
        .unwrap();
    let mut nonzero_trace = vec![0.0; geometry.ndof()];
    nonzero_trace[edge] = f64::MIN_POSITIVE;
    let zeros = vec![0.0; geometry.ndof()];
    assert!(prepared.state(&nonzero_trace, &zeros).is_err());
    assert!(prepared.state(&zeros, &nonzero_trace).is_err());
    let incompatible = rigid.state(&nonzero_trace, &zeros).unwrap();
    assert!(prepared.energy(&incompatible).is_err());
    assert!(prepared.step(&incompatible, &[], &mut || false).is_err());
    assert_eq!(incompatible.pressure_potential_pa_s(), nonzero_trace);
    // A constant potential gauge is not silently removed on import.
    let gauge = p2(&geometry, |_| 1.0);
    assert!(prepared.state(&gauge, &zeros).is_err());
    let conflict = [ResistiveBoundary {
        faces: release.clone(),
        admittance_m_pa_s: 0.0,
    }];
    assert!(
        model
            .prepare_time_with_pressure_release(0.01, &conflict, &release, limits, &mut || false)
            .is_err()
    );
    let source = velocity(&release, 0.0);
    assert!(prepared.step(&zero, &source, &mut || false).is_err());
    assert!(
        prepared
            .step_with_volume_sources(&zero, &source, &[], &mut || false)
            .is_err()
    );
    let drive = velocity(&plane(&geometry, 0.0), 0.02);
    let expected = prepared.step(&zero, &drive, &mut || false).unwrap();
    let mut calls = 0;
    assert!(matches!(
        prepared.step(&zero, &drive, &mut || {
            calls += 1;
            calls == 20
        }),
        Err(HeterogeneousError::Cancelled)
    ));
    assert!(zero.pressure_pa().iter().all(|&p| p == 0.0));
    assert!(zero.pressure_potential_pa_s().iter().all(|&psi| psi == 0.0));
    let repeated = prepared.step(&zero, &drive, &mut || false).unwrap();
    assert_eq!(repeated.state.pressure_pa(), expected.state.pressure_pa());
    assert_eq!(
        repeated.state.pressure_potential_pa_s(),
        expected.state.pressure_potential_pa_s()
    );
    assert!(matches!(
        model.prepare_time_with_pressure_release(0.01, &[], &release, limits, &mut || true),
        Err(HeterogeneousError::Cancelled)
    ));
    let initial = rigid
        .state(&p2(&geometry, |p| p[0]), &p2(&geometry, |p| 0.5 + p[1]))
        .unwrap();
    let a = rigid.step(&initial, &drive, &mut || false).unwrap();
    let b = empty.step(&initial, &drive, &mut || false).unwrap();
    assert_eq!(a.state.pressure_pa(), b.state.pressure_pa());
    assert_eq!(
        a.state.pressure_potential_pa_s(),
        b.state.pressure_potential_pa_s()
    );
    assert_eq!(
        a.work.delta_energy_j.to_bits(),
        b.work.delta_energy_j.to_bits()
    );
    assert_eq!(
        a.work.backward_error.to_bits(),
        b.work.backward_error.to_bits()
    );
}
