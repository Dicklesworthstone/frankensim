//! Original acoustic energy/adjoint + original fallible L-BFGS, with fixed disturbances.
use fs_ascent::harmonic_control::{HarmonicActuator, HarmonicControlOptions, HarmonicControlStudy};
use fs_ascent::{LbfgsError, StopReason, StopRule};
use fs_bem::helmholtz::Medium;
use fs_bem::observation_transfers::{
    HarmonicControlLimits, HarmonicField, HarmonicVolumeQuadrature, PreparedHarmonicEnergyControl,
    PreparedHarmonicObservation, TransferError,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;

const MEDIUM: Medium = Medium {
    density: 2.0,
    sound_speed: 3.0,
};
const VOLUMES: [HarmonicVolumeQuadrature; 2] = [
    HarmonicVolumeQuadrature {
        receiver: 0,
        volume_m3: 0.4,
    },
    HarmonicVolumeQuadrature {
        receiver: 1,
        volume_m3: 0.7,
    },
];

// Original public BEM preparation: two identical unit normal-velocity columns,
// plus an independently normalized aperture-flow coordinate.
fn field() -> PreparedHarmonicObservation {
    use fs_bem::observation_transfers::{BemObservation, ReceiverSpec, SourceId};
    use fs_bem::panel3d::SpherePanels;
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let cap_area: f64 = surface
        .areas()
        .iter()
        .zip(surface.normals())
        .filter(|(_, n)| n[2] > 0.5)
        .map(|(a, _)| a)
        .sum();
    let aperture = surface
        .normals()
        .iter()
        .map(|n| if n[2] > 0.5 { 1.0 / cap_area } else { 0.0 })
        .collect();
    let traces = vec![
        vec![1.0; surface.areas().len()],
        vec![1.0; surface.areas().len()],
        aperture,
    ];
    let sources = [
        SourceId::StructuralMode("uncontrollable disturbance".into()),
        SourceId::StructuralMode("velocity actuator".into()),
        SourceId::ApertureVolumeFlow("flow actuator".into()),
    ];
    let receivers = [
        ReceiverSpec {
            id: "left".into(),
            position_m: [0.4, 0.2, 0.3],
        },
        ReceiverSpec {
            id: "right".into(),
            position_m: [-0.3, 0.1, 0.5],
        },
    ];
    BemObservation::new(&surface, MEDIUM, &sources, &traces, &receivers, 0.001)
        .unwrap()
        .prepare_harmonic(7.0, || false)
        .unwrap()
}

fn objective<'a>(
    field: &'a PreparedHarmonicObservation,
    targets: Option<&'a [HarmonicField]>,
) -> PreparedHarmonicEnergyControl<'a> {
    PreparedHarmonicEnergyControl::new(
        field,
        &VOLUMES,
        targets,
        HarmonicControlLimits {
            maximum_quadrature_points: 2,
            maximum_products: 1024,
            maximum_workspace_bytes: 16384,
        },
        || false,
    )
    .unwrap()
}
fn observe(field: &PreparedHarmonicObservation, coefficients: &[C64]) -> [HarmonicField; 2] {
    let mut values = [HarmonicField::default(); 2];
    field
        .evaluate_coherent(coefficients, &mut values, || false)
        .unwrap();
    values
}
fn inner(a: &[HarmonicField], b: &[HarmonicField]) -> C64 {
    let mut sum = C64::ZERO;
    for (volume, (a, b)) in VOLUMES.iter().zip(a.iter().zip(b)) {
        let pressure = (a.pressure_pa.conj() * b.pressure_pa)
            .scale(1.0 / (MEDIUM.density * MEDIUM.sound_speed * MEDIUM.sound_speed));
        let velocity = a
            .particle_velocity_m_s
            .iter()
            .zip(b.particle_velocity_m_s)
            .fold(C64::ZERO, |s, (&a, b)| s + a.conj() * b)
            .scale(MEDIUM.density);
        sum = sum + (pressure + velocity).scale(volume.volume_m3 / 4.0);
    }
    sum
}
fn options(field: &PreparedHarmonicObservation) -> HarmonicControlOptions {
    let unit = observe(field, &[C64::ONE, C64::ZERO, C64::ZERO]);
    HarmonicControlOptions {
        energy_scale_j: inner(&unit, &unit).re,
        memory: 7,
        maximum_iterations: 200,
        maximum_evaluations: 1000,
        maximum_workspace_bytes: 1_000_000,
    }
}
fn actuator(source: usize, scale: f64, effort: f64, reference: C64) -> HarmonicActuator {
    HarmonicActuator {
        source,
        scale_si: scale,
        effort_weight: effort,
        reference,
    }
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-10 * b.abs().max(1.0), "{a:e} != {b:e}");
}
fn complex_close(a: C64, b: C64) {
    close(a.re, b.re);
    close(a.im, b.im);
}

#[test]
fn quiet_zone_cancels_a_fixed_complex_disturbance_without_turning_it_off() {
    let field = field();
    let baseline = [C64::new(0.7, -0.4), C64::new(0.1, 0.3), C64::new(-0.0, 0.0)];
    let gate = CancelGate::new_clock_free();
    let controls = [actuator(1, 0.3, 0.0, C64::ZERO)];
    let mut study = HarmonicControlStudy::new(
        objective(&field, None),
        &baseline,
        &controls,
        &[0.0; 2],
        options(&field),
        &gate,
    )
    .unwrap();
    let initial = study.mean_residual_energy_j();
    let report = study.run(&StopRule::GradNorm(1e-12), 200, &gate).unwrap();
    assert_eq!(report.optimizer.reason, StopReason::GradNorm);
    complex_close(study.coefficients()[1], -baseline[0]);
    assert!(report.mean_residual_energy_j < initial * 1e-24);
    assert_eq!(study.coefficients()[0], baseline[0]);
    assert_eq!(
        study.coefficients()[2].re.to_bits(),
        baseline[2].re.to_bits()
    );
    assert_eq!(study.effort_penalty(), 0.0);
    assert_eq!(study.baseline(), baseline);
    assert_eq!(study.objective().sources().len(), 3);
}

#[test]
fn mixed_source_scales_and_target_phase_use_the_cartesian_hermitian_gradient_once() {
    let field = field();
    let baseline = [C64::new(0.8, 0.2), C64::new(0.1, -0.1), C64::ZERO];
    let controls = [
        actuator(2, 0.3, 0.07, C64::new(0.2, 0.1)),
        actuator(1, 2.0, 0.11, C64::new(-0.3, 0.2)),
    ];
    let x = [0.15, -0.22, 0.3, -0.45];
    let target = observe(
        &field,
        &[
            C64::new(0.2, 0.1),
            C64::new(-0.1, 0.2),
            C64::new(0.03, -0.05),
        ],
    );
    let config = options(&field);
    let study = HarmonicControlStudy::new(
        objective(&field, Some(&target)),
        &baseline,
        &controls,
        &x,
        config,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let physical = observe(&field, study.coefficients());
    let residual: Vec<_> = physical
        .iter()
        .zip(&target)
        .map(|(a, b)| HarmonicField {
            pressure_pa: a.pressure_pa - b.pressure_pa,
            particle_velocity_m_s: core::array::from_fn(|c| {
                a.particle_velocity_m_s[c] - b.particle_velocity_m_s[c]
            }),
        })
        .collect();
    let mut effort = 0.0;
    let mut expected = Vec::new();
    for control in controls {
        let mut source = [C64::ZERO; 3];
        source[control.source] = C64::ONE;
        let column = observe(&field, &source);
        let physical_gradient =
            inner(&column, &residual).scale(2.0 * control.scale_si / config.energy_scale_j);
        let offset = (study.coefficients()[control.source] - control.reference)
            .scale(1.0 / control.scale_si);
        let gradient = physical_gradient + offset.scale(control.effort_weight);
        expected.extend([gradient.re, gradient.im]);
        effort += 0.5 * control.effort_weight * offset.norm_sq();
    }
    for (&a, &b) in study.state().g.iter().zip(&expected) {
        close(a, b);
    }
    close(
        study.mean_residual_energy_j() / config.energy_scale_j,
        inner(&residual, &residual).re / config.energy_scale_j,
    );
    close(study.effort_penalty(), effort);
    close(
        study.state().f,
        inner(&residual, &residual).re / config.energy_scale_j + effort,
    );
    assert!(expected.iter().skip(1).step_by(2).any(|g| g.abs() > 1e-3));
}

#[test]
fn regularized_complex_actuation_matches_the_independent_quadratic_minimizer() {
    let field = field();
    let baseline = [C64::new(0.6, -0.5), C64::ZERO, C64::ZERO];
    let control = actuator(2, 2.3, 0.12, C64::new(-0.2, 0.35));
    let config = options(&field);
    let background = observe(&field, &baseline);
    let column = observe(&field, &[C64::ZERO, C64::ZERO, C64::ONE]);
    let a = inner(&column, &column).re / config.energy_scale_j;
    let b = inner(&column, &background).scale(1.0 / config.energy_scale_j);
    let lambda = control.effort_weight / (control.scale_si * control.scale_si);
    let expected = (control.reference.scale(lambda) - b.scale(2.0)).scale(1.0 / (2.0 * a + lambda));
    let mut study = HarmonicControlStudy::new(
        objective(&field, None),
        &baseline,
        &[control],
        &[0.0; 2],
        config,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let result = study
        .run(
            &StopRule::GradNorm(1e-11),
            200,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(result.optimizer.reason, StopReason::GradNorm);
    complex_close(study.coefficients()[2], expected);
    assert_eq!(study.coefficients()[0], baseline[0]);
    assert!(study.effort_penalty() > 0.0);
    let actual = observe(&field, study.coefficients());
    close(
        result.mean_residual_energy_j / config.energy_scale_j,
        inner(&actual, &actual).re / config.energy_scale_j,
    );
}

#[test]
fn accepted_iteration_splits_preserve_original_curvature_history_and_physical_checkpoint() {
    let field = field();
    let baseline = [C64::new(0.7, -0.4), C64::ZERO, C64::ZERO];
    let controls = [
        actuator(1, 0.3, 0.07, C64::ZERO),
        actuator(2, 0.2, 0.04, C64::new(0.1, -0.2)),
    ];
    let gate = CancelGate::new_clock_free();
    let mut whole = HarmonicControlStudy::new(
        objective(&field, None),
        &baseline,
        &controls,
        &[0.0; 4],
        options(&field),
        &gate,
    )
    .unwrap();
    let mut split = HarmonicControlStudy::new(
        objective(&field, None),
        &baseline,
        &controls,
        &[0.0; 4],
        options(&field),
        &gate,
    )
    .unwrap();
    let rule = StopRule::GradNorm(0.0);
    whole.run(&rule, 8, &gate).unwrap();
    for n in [2, 1, 5] {
        if split.run(&rule, n, &gate).unwrap().optimizer.reason != StopReason::IterationCap {
            break;
        }
    }
    assert_eq!(
        format!("{:?}", whole.state()),
        format!("{:?}", split.state())
    );
    assert_eq!(whole.coefficients(), split.coefficients());
    assert_eq!(
        whole.mean_residual_energy_j(),
        split.mean_residual_energy_j()
    );
    assert_eq!(whole.effort_penalty(), split.effort_penalty());
}

#[test]
fn cancellation_inside_the_original_physics_preserves_checkpoint_and_charges_the_attempt() {
    let field = field();
    let baseline = [C64::new(0.7, -0.4), C64::ZERO, C64::ZERO];
    let controls = [
        actuator(1, 0.3, 0.07, C64::ZERO),
        actuator(2, 0.2, 0.04, C64::ZERO),
    ];
    let gate = CancelGate::new_clock_free();
    let mut study = HarmonicControlStudy::new(
        objective(&field, None),
        &baseline,
        &controls,
        &[0.0; 4],
        options(&field),
        &gate,
    )
    .unwrap();
    let checkpoint = study.state().clone();
    let coefficients = study.coefficients().to_vec();
    let energy = study.mean_residual_energy_j();
    let mut polls = 0;
    let error = study
        .run_cancellable(&StopRule::GradNorm(1e-12), 50, || {
            polls += 1;
            polls == 12
        })
        .unwrap_err();
    assert!(matches!(
        error,
        LbfgsError::Evaluation(TransferError::Cancelled)
    ));
    assert_eq!(polls, 12);
    assert_eq!(study.state().evals, checkpoint.evals + 1);
    let mut observed = study.state().clone();
    observed.evals = checkpoint.evals;
    assert_eq!(format!("{observed:?}"), format!("{checkpoint:?}"));
    assert_eq!(study.coefficients(), coefficients);
    assert_eq!(study.mean_residual_energy_j(), energy);
    let resumed = study.run(&StopRule::GradNorm(1e-10), 200, &gate).unwrap();
    assert_eq!(resumed.optimizer.reason, StopReason::GradNorm);
}

#[test]
fn hard_limits_admit_before_physics_and_never_misreport_exhaustion_as_convergence() {
    let field = field();
    let baseline = [C64::new(0.7, -0.4), C64::ZERO, C64::ZERO];
    let controls = [
        actuator(1, 0.3, 0.07, C64::ZERO),
        actuator(2, 0.2, 0.04, C64::ZERO),
    ];
    let gate = CancelGate::new_clock_free();
    let mut config = options(&field);
    config.maximum_evaluations = 1;
    let mut study = HarmonicControlStudy::new(
        objective(&field, None),
        &baseline,
        &controls,
        &[0.0; 4],
        config,
        &gate,
    )
    .unwrap();
    let saved = format!("{:?}", study.state());
    assert_eq!(
        study
            .run(&StopRule::GradNorm(0.0), 200, &gate)
            .unwrap()
            .optimizer
            .reason,
        StopReason::Budget
    );
    assert_eq!(format!("{:?}", study.state()), saved);
    let bytes = study.workspace_bytes();
    config.maximum_workspace_bytes = bytes - 1;
    assert!(
        HarmonicControlStudy::new(
            objective(&field, None),
            &baseline,
            &controls,
            &[0.0; 4],
            config,
            &gate
        )
        .is_err()
    );
    let mut duplicate = controls;
    duplicate[1].source = duplicate[0].source;
    assert!(
        HarmonicControlStudy::new(
            objective(&field, None),
            &baseline,
            &duplicate,
            &[0.0; 4],
            options(&field),
            &gate
        )
        .is_err()
    );
    let mut invalid = controls;
    invalid[1].scale_si = 0.0;
    assert!(
        HarmonicControlStudy::new(
            objective(&field, None),
            &baseline,
            &invalid,
            &[0.0; 4],
            options(&field),
            &gate
        )
        .is_err()
    );
    config = options(&field);
    config.maximum_iterations = 1;
    let mut limited = HarmonicControlStudy::new(
        objective(&field, None),
        &baseline,
        &controls,
        &[0.0; 4],
        config,
        &gate,
    )
    .unwrap();
    assert_eq!(
        limited
            .run(&StopRule::GradNorm(0.0), 200, &gate)
            .unwrap()
            .optimizer
            .reason,
        StopReason::IterationCap
    );
    assert_eq!(limited.state().iters, 1);
    let saved = format!("{:?}", limited.state());
    assert_eq!(
        limited
            .run(&StopRule::GradNorm(0.0), 200, &gate)
            .unwrap()
            .optimizer
            .reason,
        StopReason::IterationCap
    );
    assert_eq!(format!("{:?}", limited.state()), saved);
}
