//! Native trajectory, work/refinement and transactional source regressions.
use super::*;
use fs_solid::membrane::vibration::VibrationSettings;
use fs_time::galpha::{ImplicitSolveConfig, SecondOrderState};

pub(in crate::cavity_phs::membrane) fn motion_controls() -> GasMembraneMotionControls {
    GasMembraneMotionControls {
        gas: gas_controls(),
        time: ImplicitSolveConfig {
            newton: controls().membrane.newton,
            max_newton_iterations: 40,
        },
        maximum_evaluations: 1000,
        maximum_steps: 512,
        maximum_krylov_restart: 8,
        energy_tolerance_j: 1e-7,
    }
}

#[test]
fn small_motion_and_changed_held_load_follow_the_original_gas_loaded_mode() {
    let (patch, reference, fixed) = pyramid();
    let model = GasMembrane::new(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        P0,
        [0.0007, 0.0015],
        || false,
    )
    .unwrap();
    let modal = model
        .solve_vibration(
            &[[0.0; 3]; 5],
            &reference,
            controls(),
            VibrationSettings {
                max_free_dofs: 1,
                relative_residual_tolerance: 1e-8,
            },
            || false,
        )
        .unwrap();
    let omega = modal.basis().modes()[0].lambda.sqrt();
    let stiffness = omega * omega * patch.nodal_mass_kg()[0];
    let dt = 0.025 / omega;
    // Exact linear midpoint phase, independent of the nonlinear integrator.
    let theta = 2.0 * (omega * dt / 2.0).atan();
    let (sine, cosine) = theta.sin_cos();
    let amplitude = 1e-6;
    let mut initial = reference.clone();
    initial[0][2] = amplitude;
    let mut loads = [[0.0; 3]; 5];
    let mut motion = GasMembraneMotion::new(
        &model,
        &initial,
        &[[0.0; 3]; 5],
        &loads,
        motion_controls(),
        || false,
    )
    .unwrap();
    let (mut q, mut v) = (amplitude, 0.0);
    for step in 0..64 {
        if step == 32 {
            loads[0][2] = 0.002;
        }
        let equilibrium = loads[0][2] / stiffness;
        let next_q = equilibrium + (q - equilibrium) * cosine + v / omega * sine;
        let next_v = -omega * (q - equilibrium) * sine + v * cosine;
        (q, v) = (next_q, next_v);
        let frame = motion
            .step(dt, &loads, motion_controls(), || false)
            .unwrap();
        near(motion.state().q[0], q, 2e-4 * amplitude, 0.0);
        near(motion.state().v[0], v, 2e-4 * amplitude * omega, 0.0);
        assert!(
            frame.gas_trials < frame.evaluations,
            "native tangent callbacks reuse a gas trial"
        );
        assert!(frame.integration_defect_j.abs() <= motion_controls().energy_tolerance_j);
    }
    assert_eq!(motion.state().steps, 64);
    assert_eq!(motion.state().history.len(), 64);
    let velocity = motion.velocity_m_s();
    for (i, prescribed) in fixed.iter().enumerate() {
        if let Some(value) = prescribed {
            assert_eq!(motion.accepted().displacement_m[i / 3][i % 3], *value);
            assert_eq!(velocity[i / 3][i % 3], 0.0);
        }
    }
    assert_eq!(
        motion.accepted().gas.mass_kg(),
        model.reference_gas().mass_kg()
    );
}

#[test]
fn finite_amplitude_motion_refines_and_closes_measured_interval_energy_accounts() {
    let (patch, reference, fixed) = pyramid();
    let model = GasMembrane::new(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        P0,
        [0.0007, 0.0015],
        || false,
    )
    .unwrap();
    let mut initial = reference.clone();
    initial[0][2] = 0.008;
    let mut velocity = [[0.0; 3]; 5];
    velocity[0][2] = 0.15;
    let loads = [[0.0; 3]; 5];
    let controls = GasMembraneMotionControls {
        energy_tolerance_j: 1e-3,
        ..motion_controls()
    };
    let run = |steps: usize| {
        let mut motion =
            GasMembraneMotion::new(&model, &initial, &velocity, &loads, controls, || false)
                .unwrap();
        let before = motion.accepted().clone();
        let kinetic_before = motion.kinetic_energy_j().unwrap();
        let (mut defects, mut worst) = (0.0, 0.0_f64);
        for _ in 0..steps {
            let old = motion.accepted().clone();
            let frame = motion
                .step(0.0004 / steps as f64, &loads, controls, || false)
                .unwrap();
            let new = motion.accepted();
            near(
                frame.gas_energy_change_j,
                new.gas_volume.stored_energy_change_j - old.gas_volume.stored_energy_change_j,
                1e-14,
                1e-13,
            );
            near(frame.work_by_gas_j, -frame.gas_energy_change_j, 0.0, 0.0);
            let swept = patch
                .pressure_work_j(&old.displacement_m, &new.displacement_m, 1.0, || false)
                .unwrap();
            near(frame.swept_volume_change_m3, swept, 0.0, 0.0);
            near(frame.exterior_pressure_work_j, P0 * swept, 0.0, 0.0);
            near(
                frame.integration_defect_j,
                frame.kinetic_energy_change_j
                    + frame.membrane_energy_change_j
                    + frame.gas_energy_change_j
                    + frame.exterior_pressure_work_j
                    - frame.dead_load_work_j,
                1e-13,
                1e-10,
            );
            assert_eq!(frame.dead_load_work_j, 0.0);
            defects += frame.integration_defect_j;
            worst = worst.max(frame.integration_defect_j.abs());
        }
        let after = motion.accepted();
        let balance = motion.kinetic_energy_j().unwrap() - kinetic_before + after.strain_energy_j
            - before.strain_energy_j
            + after.gas_volume.stored_energy_change_j
            - before.gas_volume.stored_energy_change_j
            + P0 * patch
                .pressure_work_j(&before.displacement_m, &after.displacement_m, 1.0, || false)
                .unwrap();
        near(defects, balance, 5e-10, 2e-7);
        assert!((after.gas.volume_m3() - before.gas.volume_m3()).abs() > 1e-6);
        (motion.state().q[0], motion.state().v[0], worst)
    };
    let (coarse, fine, reference_run) = (run(16), run(32), run(64));
    let rate = ((2.0 * STIFFNESS * (STRETCH.powi(2) - 1.0)
        + GAMMA * P0 / V0 * volume_gradient().powi(2))
        / patch.nodal_mass_kg()[0])
        .sqrt();
    let error = |state: (f64, f64, f64)| {
        (state.0 - reference_run.0).hypot((state.1 - reference_run.1) / rate)
    };
    assert!(error(coarse) > 2.5 * error(fine));
    assert!(
        coarse.2 > 2.0 * fine.2 && fine.2 > 2.0 * reference_run.2,
        "the actual per-step nonlinear energy defect must refine"
    );
}

#[test]
fn curved_full_coordinate_motion_balances_native_pressure_at_the_midpoint() {
    let (patch, reference, fixed) =
        chart(Mesh2::triangles(0.3, 0.2, 3, 2), &[15, 16, 17, 18, 19, 20]);
    let model = GasMembrane::new(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        80000.0,
        [0.0007, 0.0015],
        || false,
    )
    .unwrap();
    let mut initial = reference.clone();
    initial[5] = [initial[5][0] + 0.0004, initial[5][1] - 0.0003, 0.006];
    initial[6] = [initial[6][0] - 0.0006, initial[6][1] + 0.0002, 0.009];
    let mut loads = model
        .evaluate(
            &initial,
            &vec![[0.0; 3]; initial.len()],
            gas_controls(),
            || false,
        )
        .unwrap()
        .residual_n;
    loads[5][2] += 0.3;
    loads[6][0] -= 0.2;
    let mut velocity = vec![[0.0; 3]; initial.len()];
    velocity[5] = [0.01, -0.02, 0.03];
    velocity[6] = [-0.02, 0.01, -0.01];
    let mut motion = GasMembraneMotion::new(
        &model,
        &initial,
        &velocity,
        &loads,
        motion_controls(),
        || false,
    )
    .unwrap();
    let before = motion.state().clone();
    let dt = 1e-5;
    motion
        .step(dt, &loads, motion_controls(), || false)
        .unwrap();
    let midpoint: Vec<_> = initial
        .iter()
        .zip(&motion.accepted().displacement_m)
        .map(|(a, b)| std::array::from_fn(|j| 0.5 * a[j] + 0.5 * b[j]))
        .collect();
    let stage = model
        .evaluate(&midpoint, &loads, gas_controls(), || false)
        .unwrap();
    for (i, &dof) in model.free_dofs().iter().enumerate() {
        near(
            patch.nodal_mass_kg()[dof / 3] * (motion.state().v[i] - before.v[i])
                + dt * stage.residual_n[dof / 3][dof % 3],
            0.0,
            2e-12,
            0.0,
        );
        near(
            motion.state().q[i] - before.q[i],
            dt * (0.5 * motion.state().v[i] + 0.5 * before.v[i]),
            1e-15,
            1e-10,
        );
    }
    let geometric = patch
        .pressure_tangent_apply(
            &midpoint,
            &velocity,
            stage.pressure.pressure_pa - 80000.0,
            || false,
        )
        .unwrap();
    assert!(geometric.iter().flatten().any(|value| value.abs() > 0.01));
}

pub(in crate::cavity_phs::membrane) fn real_gas() -> ClosedRealGasParcel {
    let ideal = source_gas(0.75, 0.001);
    let source = ThermochemicalPengRobinsonModel::new(
        ideal.model().clone(),
        vec![
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("A").unwrap(),
                critical_temperature: Temperature::new(190.564),
                critical_pressure: Pressure::new(4.5992e6),
                acentric_factor: 0.01142,
            },
            PengRobinsonSpeciesBinding {
                species: SpeciesId::new("B").unwrap(),
                critical_temperature: Temperature::new(150.687),
                critical_pressure: Pressure::new(4.863e6),
                acentric_factor: -0.00219,
            },
        ],
        vec![0.0, 0.03, 0.03, 0.0],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            300.0,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 5e6)
        .unwrap()
        .with("specimen-lot", 42.0)
        .unwrap();
    ClosedRealGasParcel::new(source, point, V0, PengRobinsonPhaseV1::Vapor).unwrap()
}

#[test]
fn accepted_motion_preserves_variable_caloric_and_real_gas_owners() {
    let (patch, reference, fixed) = pyramid();
    let mut velocity = [[0.0; 3]; 5];
    velocity[0][2] = 0.05;
    for source in [
        MembraneGas::Ideal(source_gas(0.75, 0.002)),
        MembraneGas::PengRobinson(real_gas()),
    ] {
        let controls = GasMembraneMotionControls {
            gas: GasEntropySolve {
                minimum_temperature_k: 280.0,
                maximum_temperature_k: 450.0,
                ..gas_controls()
            },
            ..motion_controls()
        };
        let pressure = source.pressure_response().unwrap().pressure_pa;
        let model = GasMembrane::new(
            &patch,
            source.clone(),
            &reference,
            &fixed,
            pressure,
            [0.0007, 0.0015],
            || false,
        )
        .unwrap();
        let mut motion = GasMembraneMotion::new(
            &model,
            &reference,
            &velocity,
            &[[0.0; 3]; 5],
            controls,
            || false,
        )
        .unwrap();
        for _ in 0..5 {
            motion
                .step(1e-6, &[[0.0; 3]; 5], controls, || false)
                .unwrap();
        }
        let accepted = motion.accepted();
        assert!(accepted.displacement_m[0][2] > 0.0);
        assert!(motion.state().v[0] < velocity[0][2]);
        let (direct, work) = match source.clone() {
            MembraneGas::Ideal(mut gas) => {
                let work = gas
                    .set_volume_isentropic(accepted.gas.volume_m3(), controls.gas, || false)
                    .unwrap();
                assert_eq!(gas.point().axes().get("specimen-lot"), Some(&42.0));
                (MembraneGas::Ideal(gas), work)
            }
            MembraneGas::PengRobinson(mut gas) => {
                let work = gas
                    .set_volume_isentropic(accepted.gas.volume_m3(), controls.gas, || false)
                    .unwrap();
                assert_eq!(gas.point().axes().get("specimen-lot"), Some(&42.0));
                (MembraneGas::PengRobinson(gas), work)
            }
        };
        assert_eq!(accepted.gas, direct);
        assert_eq!(accepted.gas_volume, work);
        assert_eq!(accepted.pressure, direct.pressure_response().unwrap());
        assert_eq!(accepted.gas.mass_kg(), source.mass_kg());
        assert_eq!(model.reference_gas(), &source);
    }
}

fn unchanged(
    motion: &GasMembraneMotion<'_, '_>,
    state: &SecondOrderState,
    accepted: &GasMembraneEvaluation,
) {
    assert_eq!(motion.state(), state);
    assert_eq!(motion.accepted().gas, accepted.gas);
    assert_eq!(motion.accepted().gas_volume, accepted.gas_volume);
    assert_eq!(motion.accepted().displacement_m, accepted.displacement_m);
    assert_eq!(motion.accepted().residual_n, accepted.residual_n);
    assert_eq!(
        motion.accepted().potential_energy_change_j,
        accepted.potential_energy_change_j
    );
}

#[test]
fn late_energy_cancel_and_budget_refusals_preserve_all_accepted_state() {
    let (patch, reference, fixed) = pyramid();
    let model = GasMembrane::new(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        P0,
        [0.0007, 0.0015],
        || false,
    )
    .unwrap();
    let mut initial = reference.clone();
    initial[0][2] = 0.006;
    let mut velocity = [[0.0; 3]; 5];
    velocity[0][2] = 0.2;
    let loads = [[0.0; 3]; 5];
    let controls = GasMembraneMotionControls {
        energy_tolerance_j: 1e-3,
        ..motion_controls()
    };
    let mut motion =
        GasMembraneMotion::new(&model, &initial, &velocity, &loads, controls, || false).unwrap();
    motion.step(2e-5, &loads, controls, || false).unwrap();
    let state = motion.state().clone();
    let accepted = motion.accepted().clone();
    assert!(matches!(
        motion.step(
            2e-5,
            &loads,
            GasMembraneMotionControls {
                maximum_evaluations: 1,
                ..controls
            },
            || false
        ),
        Err(GasMembraneMotionError::EvaluationBudget)
    ));
    unchanged(&motion, &state, &accepted);
    assert!(matches!(
        motion.step(
            2e-5,
            &loads,
            GasMembraneMotionControls {
                maximum_steps: state.steps,
                ..controls
            },
            || false
        ),
        Err(GasMembraneMotionError::StepBudget)
    ));
    unchanged(&motion, &state, &accepted);
    assert!(matches!(
        motion.step(
            2e-5,
            &loads,
            GasMembraneMotionControls {
                energy_tolerance_j: 1e-25,
                ..controls
            },
            || false
        ),
        Err(GasMembraneMotionError::Energy { .. })
    ));
    unchanged(&motion, &state, &accepted);
    // Cancel at the final publication checkpoint, after a matching native
    // proposal has completed its source and energy checks in a dry candidate.
    let mut probe = motion.clone();
    let mut total_polls = 0;
    probe
        .step(2e-5, &loads, controls, || {
            total_polls += 1;
            false
        })
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        motion.step(2e-5, &loads, controls, || {
            polls += 1;
            polls >= total_polls
        }),
        Err(GasMembraneMotionError::Cancelled)
    ));
    unchanged(&motion, &state, &accepted);
    motion.step(2e-5, &loads, controls, || false).unwrap();
    assert_eq!(motion.state(), probe.state());
    assert_eq!(motion.accepted().gas, probe.accepted().gas);

    let narrow = GasMembrane::new(
        &patch,
        source_gas(0.75, 0.002),
        &reference,
        &fixed,
        P0,
        [V0 - 1e-9, V0 + 1e-9],
        || false,
    )
    .unwrap();
    velocity[0][2] = 5.0;
    let mut blocked =
        GasMembraneMotion::new(&narrow, &reference, &velocity, &loads, controls, || false).unwrap();
    let state = blocked.state().clone();
    let accepted = blocked.accepted().clone();
    assert!(blocked.step(1e-4, &loads, controls, || false).is_err());
    unchanged(&blocked, &state, &accepted);
}
