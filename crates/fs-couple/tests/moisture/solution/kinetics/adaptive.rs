//! Independent continuous water/phase dynamics and optional component accuracy.
use super::*;
use fs_couple::moisture::{
    adaptive::{KineticPhaseAccuracy, MoistureAdvanceError, MoistureTimeControls},
    atmosphere::{CartesianMoistureParticipant, DryBasisSorption, MoistureParticipant},
};

const VOLUME: f64 = 1e-4;
const CONDUCTANCE: f64 = 2.5e-5;
const INITIAL: [f64; 2] = [0.01, 0.01]; // dry-basis water, dissolved mol
const DURATION: f64 = 5.0;
const WATER_RATE: f64 = 0.08;
const DRY_MASS: f64 = 0.5;

fn time(step: f64) -> MoistureTimeControls {
    MoistureTimeControls {
        initial_step_s: step,
        minimum_step_s: 1e-7,
        maximum_step_s: step,
        maximum_evaluations: 5000,
        // Deliberately loose: the actual phase stock alone must drive refinement.
        moisture_absolute_tolerance: 1.0,
        activity_absolute_tolerance: 1.0,
        reservoir_water_absolute_tolerance_kg: 1.0,
        relative_tolerance: 0.0,
        water_balance_tolerance_kg: 1e-10,
    }
}

fn phase_accuracy(tolerance: f64) -> KineticPhaseAccuracy {
    KineticPhaseAccuracy {
        absolute_amount_tolerance_mol: tolerance,
        relative_tolerance: 0.0,
    }
}

fn original() -> FiniteMoistureSolution {
    kinetic_bath(
        VOLUME,
        INITIAL[1],
        0.2 - INITIAL[1],
        CONDUCTANCE,
        PrecipitationPolicyV1::ExistingSolidOnly,
    )
}

fn exposed() -> [Option<DryBasisSorption>; 2] {
    [
        Some(atmosphere::sorption(0.2)),
        Some(atmosphere::sorption(0.2)),
    ]
}

// Independent continuous two-state model. Water conservation eliminates the
// bath volume. Native endpoint iteration, phase stepping and adaptive scheduling
// are not used to construct this RK4 reference. Crystals remain present.
fn reference(steps: usize) -> [f64; 2] {
    let rhs = |y: [f64; 2]| {
        let volume = VOLUME - DRY_MASS * (y[0] - INITIAL[0]) / WATER_DENSITY;
        let concentration = y[1] / volume;
        [
            WATER_RATE * (0.2 * (-WATER_VOLUME * concentration).exp() - y[0]),
            CONDUCTANCE * (1000.0 - concentration),
        ]
    };
    let mut y = INITIAL;
    let h = DURATION / steps as f64;
    for _ in 0..steps {
        let a = rhs(y);
        let b = rhs(std::array::from_fn(|i| y[i] + 0.5 * h * a[i]));
        let c = rhs(std::array::from_fn(|i| y[i] + 0.5 * h * b[i]));
        let d = rhs(std::array::from_fn(|i| y[i] + h * c[i]));
        y = std::array::from_fn(|i| y[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}

#[test]
fn kinetic_phase_accuracy_refines_native_water_and_crystal_dynamics() {
    with_cx(|cx| {
        let bath = original();
        let body = slab(card([0.02; 2], 0.3), 1.0, vec![INITIAL[0]], cx);
        near(body.dry_mass_kg(), DRY_MASS, 1e-14);
        let exact = reference(20_000);
        for (a, b) in exact.into_iter().zip(reference(40_000)) {
            near(a, b, 2e-12);
        }
        let legacy = bath
            .advance_adaptive_to(
                &body,
                DURATION,
                &exposed(),
                config(),
                controls(VOLUME),
                time(DURATION),
                cx,
            )
            .unwrap();
        assert_eq!(legacy.2.evaluations, 3);
        let mut errors = Vec::new();
        let mut evaluations = Vec::new();
        for tolerance in [5e-4, 5e-6] {
            let (next_bath, next, report) = bath
                .advance_adaptive_to(
                    &body,
                    DURATION,
                    &exposed(),
                    config(),
                    controls(VOLUME),
                    time(DURATION).with_kinetic_phase_accuracy(phase_accuracy(tolerance)),
                    cx,
                )
                .unwrap();
            let actual = phase(&next_bath);
            near(actual.total_solute().value(), 0.2, 1e-14);
            near(
                actual.solid().value() + actual.dissolved().value(),
                0.2,
                1e-14,
            );
            near(
                next_bath.solvent_mass_kg() + next.water_mass_kg(),
                bath.solvent_mass_kg() + body.water_mass_kg(),
                1e-10,
            );
            near(report.water_balance_residual_kg, 0.0, 1e-10);
            assert!(report.maximum_accepted_error_ratio <= 1.0);
            assert_eq!(next.time_s(), DURATION);
            if tolerance == 5e-4 {
                assert_eq!(next.moisture(), legacy.1.moisture());
                assert_eq!(
                    next_bath.kinetic_solubility(),
                    legacy.0.kinetic_solubility()
                );
            }
            errors.push((actual.dissolved().value() - exact[1]).abs());
            evaluations.push(report.evaluations);
        }
        assert!(errors[1] < 0.1 * errors[0], "{errors:?}");
        assert!(errors[1] < 8e-5, "{errors:?}");
        assert!(evaluations[1] > evaluations[0]);
        assert_eq!(phase(&bath).dissolved().value(), INITIAL[1]);
        assert_eq!(body.time_s(), 0.0);
    });
}

#[test]
fn phase_accuracy_reaches_cartesian_and_shared_native_reservoir_paths() {
    with_cx(|cx| {
        let bath = original();
        let body = slab(card([0.02; 2], 0.3), 1.0, vec![INITIAL[0]], cx);
        // Same dry volume and relaxation rate as the one-cell slab:
        // 12*D_cartesian/L^2 = 4*D_slab/L_slab^2 = 0.08/s.
        let cartesian = grid(0.1, WATER_RATE * 0.01 / 12.0, INITIAL[0], cx);
        near(cartesian.dry_mass_kg(), DRY_MASS, 1e-14);
        let faces = exposed();
        let cube_faces = std::array::from_fn(|_| Some(atmosphere::sorption(0.2)));
        let accuracy = time(DURATION).with_kinetic_phase_accuracy(phase_accuracy(5e-6));
        let (one, _, one_report) = bath
            .advance_adaptive_to(
                &body,
                DURATION,
                &faces,
                config(),
                controls(VOLUME),
                accuracy,
                cx,
            )
            .unwrap();
        let (cube, cube_body, cube_report) = bath
            .advance_cartesian_adaptive_to(
                &cartesian,
                DURATION,
                &cube_faces,
                config(),
                nonlinear(),
                controls(VOLUME),
                accuracy,
                cx,
            )
            .unwrap();
        near(
            phase(&cube).dissolved().value(),
            phase(&one).dissolved().value(),
            2e-9,
        );
        assert!(cube_report.evaluations > 3);
        assert!(one_report.evaluations > 3);
        let shared = bath
            .advance_shared_adaptive_to(
                &[MoistureParticipant {
                    id: 31,
                    slab: &body,
                    faces: &faces,
                    transport: config(),
                }],
                DURATION,
                1,
                controls(VOLUME),
                accuracy,
                cx,
            )
            .unwrap();
        near(
            phase(shared.reservoir()).dissolved().value(),
            phase(&one).dissolved().value(),
            2e-12,
        );
        assert_eq!(shared.solids()[0].id, 31);
        assert_eq!(shared.report().summary.evaluations, one_report.evaluations);
        let shared_cube = bath
            .advance_shared_cartesian_adaptive_to(
                &[CartesianMoistureParticipant {
                    id: 71,
                    moisture: &cartesian,
                    faces: &cube_faces,
                    transport: config(),
                    nonlinear: nonlinear(),
                }],
                DURATION,
                1,
                controls(VOLUME),
                accuracy,
                cx,
            )
            .unwrap();
        near(
            phase(shared_cube.reservoir()).dissolved().value(),
            phase(&cube).dissolved().value(),
            2e-12,
        );
        assert_eq!(
            shared_cube.solids()[0].moisture.moisture(),
            cube_body.moisture()
        );
        assert_eq!(
            shared_cube.report().summary.evaluations,
            cube_report.evaluations
        );
        assert_eq!(phase(&bath).dissolved().value(), INITIAL[1]);
    });
}

#[test]
fn exact_fixed_volume_phase_evolution_does_not_trigger_spurious_time_refinement() {
    with_cx(|cx| {
        let bath = original();
        let body = slab(card([0.02; 2], 0.3), 1.0, vec![INITIAL[0]], cx);
        let (next_bath, next, report) = bath
            .advance_adaptive_to(
                &body,
                DURATION,
                &[None, None],
                config(),
                controls(VOLUME),
                time(DURATION).with_kinetic_phase_accuracy(phase_accuracy(1e-10)),
                cx,
            )
            .unwrap();
        let exact = 1000.0 * VOLUME
            + (INITIAL[1] - 1000.0 * VOLUME) * (-CONDUCTANCE * DURATION / VOLUME).exp();
        near(phase(&next_bath).dissolved().value(), exact, 2e-14);
        assert_eq!(next.moisture(), body.moisture());
        assert_eq!(next_bath.solvent_mass_kg(), bath.solvent_mass_kg());
        assert_eq!(report.evaluations, 3);
        assert_eq!(report.accepted_intervals, 1);
        assert_eq!(report.rejected_intervals, 0);
    });
}

#[test]
fn invalid_phase_controls_and_late_work_refusal_preserve_both_originals() {
    with_cx(|cx| {
        let bath = original();
        let saved = bath.kinetic_solubility().unwrap().clone();
        let body = slab(card([0.02; 2], 0.3), 1.0, vec![INITIAL[0]], cx);
        for phase in [
            phase_accuracy(0.0),
            phase_accuracy(f64::NAN),
            KineticPhaseAccuracy {
                absolute_amount_tolerance_mol: 1e-6,
                relative_tolerance: 1.0,
            },
        ] {
            assert!(matches!(
                bath.advance_adaptive_to(
                    &body,
                    1.0,
                    &exposed(),
                    config(),
                    controls(VOLUME),
                    time(1.0).with_kinetic_phase_accuracy(phase),
                    cx,
                ),
                Err(MoistureAdvanceError::Invalid(_))
            ));
        }
        let plain = reservoir(VOLUME, vec![INITIAL[1]]);
        assert!(matches!(
            plain.advance_adaptive_to(
                &body,
                1.0,
                &exposed(),
                config(),
                controls(VOLUME),
                time(1.0).with_kinetic_phase_accuracy(phase_accuracy(1e-6)),
                cx,
            ),
            Err(MoistureAdvanceError::Invalid(_))
        ));
        let step = 1.0 / 64.0;
        let mut limited = time(step);
        limited.maximum_evaluations = 3;
        assert!(matches!(bath.advance_adaptive_to(
            &body, 2.0*step, &exposed(), config(), controls(VOLUME),
            limited.with_kinetic_phase_accuracy(phase_accuracy(1.0)), cx,
        ), Err(MoistureAdvanceError::WorkLimit {
            time_s, accepted_intervals: 1, evaluations: 3,
        }) if time_s == step));
        assert_eq!(bath.kinetic_solubility(), Some(&saved));
        assert_eq!(body.time_s(), 0.0);
        assert_eq!(body.moisture(), &[INITIAL[0]]);
    });
}
