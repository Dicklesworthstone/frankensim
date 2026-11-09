//! G0/G1/G3/G4: adaptive finite heat/recoil feedback, continuous manufactured
//! equations, accepted-only accounting, frame transformations and atomic refusal.

mod exchange;

use super::*;
use crate::rocket_blowdown::effusion::thermal::trajectory::{
    EffusionWallAdaptiveSolve, EffusionWallTrajectoryCandidate, EffusionWallTrajectoryError,
};
use crate::rocket_blowdown::effusion::trajectory::{
    EffusionAdaptiveSolve, EffusionTrajectoryError,
};

fn adaptive(scale: f64) -> EffusionWallAdaptiveSolve {
    EffusionWallAdaptiveSolve {
        motion: EffusionAdaptiveSolve {
            initial_step_s: 0.1,
            minimum_step_s: 1e-6,
            maximum_step_s: 0.1,
            maximum_attempts: 4096,
            maximum_accepted_intervals: 2048,
            relative_tolerance: 1e-3 * scale,
            temperature_tolerance_k: 0.05 * scale,
            pressure_tolerance_pa: 5.0 * scale,
            species_mass_tolerance_kg: 1e-9 * scale,
            position_tolerance_m: 1e-5 * scale,
            velocity_tolerance_m_s: 1e-4 * scale,
        },
        wall_temperature_tolerance_k: 0.01 * scale,
        wall_energy_tolerance_j: 0.01 * scale,
        emitted_momentum_tolerance_n_s: 1e-7 * scale,
        emitted_energy_tolerance_j: 1e-4 * scale,
        window_mass_tolerance_kg: 1e-11,
        window_momentum_tolerance_n_s: 1e-9,
        window_energy_tolerance_j: 1e-4,
    }
}

fn check_window<S>(
    before: &GasRocketState,
    input: UniformEnthalpyStepInput<'_>,
    candidate: &EffusionWallTrajectoryCandidate<S>,
) {
    let step = candidate.step();
    let after = candidate.endpoint();
    assert_eq!(step.duration_s, input.duration_s);
    assert_eq!(step.accepted_steps.len(), 2 * step.accepted_intervals);
    assert_eq!(
        step.attempted_intervals,
        step.accepted_intervals + step.rejected_intervals
    );
    assert!(step.maximum_accepted_error <= 1.0);
    near(
        after.reservoir().mass_kg() - before.reservoir().mass_kg() + step.escaped_mass_kg,
        0.0,
        1e-11,
    );
    vector_near(
        after
            .vehicle()
            .velocity_m_s
            .scale(mass(after))
            .sub(before.vehicle().velocity_m_s.scale(mass(before)))
            .add(step.escaped_momentum_world_kg_m_s)
            .sub(step.external_impulse_n_s),
        Vec3::ZERO,
        2e-9,
    );
    near(
        internal(after) - internal(before) + kinetic(after) - kinetic(before)
            + step.wall_storage_change_j
            + step.escaped_energy_world_j
            - step.external_work_j
            - step.external_thermal_heat_j,
        0.0,
        1e-4,
    );
    near(
        step.wall_storage_change_j,
        input.mass_kg
            * (candidate.wall().specific_enthalpy_j_kg() - input.initial.specific_enthalpy_j_kg()),
        1e-12,
    );
    assert!(step.debit.energy_j <= adaptive(1.0).window_energy_tolerance_j);
    assert!(step.maximum_species_mass_residual_kg <= adaptive(1.0).window_mass_tolerance_kg);
    assert_eq!(after.vehicle().dry_mass_kg, before.vehicle().dry_mass_kg);
}

// Independent continuous monatomic balances. n_dot = -A*n/V*sqrt(R*T/(2*pi*M));
// n*Cv*T_dot = Q + (R/2)*T*n_dot, C_wall*T_wall_dot = -Q, and
// (m_dry+n*M)*v_dot = -p*A/2. The emitted energy includes the SAME inertial
// translation as the escaping mass. No library thermal, effusion or MBD step is
// used in this fine RK4 reference. Manufactured constant h is not empirical data.
fn continuous(before: &GasRocketState, duration: f64) -> [f64; 7] {
    let derivative = |y: [f64; 7]| {
        let [n, t, w, _, v, _, _] = y;
        let dn = -(AREA / VOLUME) * n * (R * t / (2.0 * core::f64::consts::PI * 0.028)).sqrt();
        let heat = 100.0 * WALL_AREA * (w - t);
        let impulse = 0.5 * AREA * n * R * t / VOLUME;
        let escaped_mass = -0.028 * dn;
        [
            dn,
            heat / (1.5 * n * R) + t * dn / (3.0 * n),
            -heat / (WALL_MASS * WALL_CP),
            v,
            -impulse / (before.vehicle().dry_mass_kg + 0.028 * n),
            escaped_mass * v + impulse,
            -dn * 2.0 * R * t + v * impulse + 0.5 * escaped_mass * v * v,
        ]
    };
    let mut y = [P0 * VOLUME / (R * T0), T0, 600.0, 0.0, 0.0, 0.0, 0.0];
    let dt = duration / 4096.0;
    for _ in 0..4096 {
        let a = derivative(y);
        let b = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
        let c = derivative(core::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
        let d = derivative(core::array::from_fn(|i| y[i] + dt * c[i]));
        for i in 0..7 {
            y[i] += dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    y
}

#[test]
fn tighter_windows_converge_to_continuous_finite_wall_cooling_effusion_and_recoil() {
    scope(|cx, _| {
        let before = initial(PURE, Vec3::ZERO);
        let saved = before.clone();
        let law = interface(before.reservoir().model());
        let outlets = [outlet(&law, Vec3::new(1.0, 0.0, 0.0))];
        let chart = curve([1000.0; 2]);
        let input = UniformEnthalpyStepInput {
            duration_s: 0.3,
            ..wall(&chart, 600.0)
        };
        let run = |accuracy| {
            before
                .advance_effusion_with_wall_adaptive(
                    cx,
                    free(),
                    &outlets,
                    input,
                    &body(100.0),
                    0.0,
                    solve(),
                    adaptive(accuracy),
                )
                .unwrap()
        };
        let loose = run(1.0);
        let tight = run(0.04);
        let reference = continuous(&before, input.duration_s);
        let discrepancy = |candidate: &EffusionWallTrajectoryCandidate| {
            let end = candidate.endpoint();
            let actual = [
                end.reservoir().mass_kg() / 0.028,
                end.reservoir().state().state().temperature,
                candidate.wall().temperature_k(),
                end.vehicle().position_m.x,
                end.vehicle().velocity_m_s.x,
                candidate.step().escaped_momentum_world_kg_m_s.x,
                candidate.step().escaped_energy_world_j,
            ];
            actual
                .into_iter()
                .zip(reference)
                .zip([0.03, 400.0, 600.0, 0.01, 0.1, 0.002, 2.0])
                .map(|((a, r), scale)| (a - r).abs() / scale)
                .sum::<f64>()
        };
        assert!(discrepancy(&tight) < 0.4 * discrepancy(&loose));
        assert!(tight.step().accepted_intervals > loose.step().accepted_intervals);
        assert!(tight.endpoint().reservoir().state().state().temperature > T0);
        assert!(tight.wall().temperature_k() < 600.0);
        assert!(tight.endpoint().vehicle().velocity_m_s.x < 0.0);
        check_window(&before, input, &tight);
        assert_eq!(before, saved);
    });
}

#[test]
fn accepted_history_replays_exactly_and_prorated_sources_reach_the_exact_final_time() {
    scope(|cx, _| {
        let before = initial(MIXTURE, Vec3::new(2.0, -1.0, 0.5));
        let chart = curve([1000.0; 2]);
        let input = UniformEnthalpyStepInput {
            internal_heat_j: 2.7,
            ..wall(&chart, 600.0)
        };
        let law = interface(before.reservoir().model());
        let outlets = [outlet(&law, Vec3::new(1.0, 0.0, 0.0))];
        let boundary = body(100.0);
        let mut accuracy = adaptive(1e6);
        accuracy.motion.initial_step_s = 0.0125;
        accuracy.motion.maximum_step_s = 0.0125;
        let candidate = before
            .advance_effusion_with_wall_adaptive(
                cx,
                free(),
                &outlets,
                input,
                &boundary,
                1.3,
                solve(),
                accuracy,
            )
            .unwrap();
        let mut replay = before.clone();
        let mut replay_wall = input.initial;
        for receipt in &candidate.step().accepted_steps {
            let fraction = receipt.duration_s / input.duration_s;
            let accepted = replay
                .advance_effusion_with_wall(
                    cx,
                    free(),
                    &outlets,
                    UniformEnthalpyStepInput {
                        initial: replay_wall,
                        duration_s: receipt.duration_s,
                        internal_heat_j: input.internal_heat_j * fraction,
                        ..input
                    },
                    &boundary,
                    1.3 * fraction,
                    solve(),
                )
                .unwrap();
            assert_eq!(accepted.step(), receipt);
            (replay, replay_wall, _) = accepted.into_parts();
        }
        assert_eq!(candidate.endpoint(), &replay);
        assert_eq!(candidate.wall(), &replay_wall);
        near(candidate.step().external_thermal_heat_j, 4.0, 1e-14);
        for residual in candidate.step().prescribed_heat_residual_j {
            near(residual, 0.0, 1e-14);
        }
        check_window(&before, input, &candidate);
    });
}

#[test]
fn common_inertial_boost_transforms_aggregate_escaped_momentum_and_energy_once() {
    scope(|cx, _| {
        let at_rest = initial(MIXTURE, Vec3::ZERO);
        let boost = Vec3::new(11.0, -3.0, 7.0);
        let moving = initial(MIXTURE, boost);
        let chart = curve([1000.0; 2]);
        let input = wall(&chart, 600.0);
        let law = interface(at_rest.reservoir().model());
        let outlets = [outlet(&law, Vec3::new(0.0, 1.0, 0.0))];
        let mut accuracy = adaptive(1e6);
        accuracy.motion.initial_step_s = 0.025;
        accuracy.motion.maximum_step_s = 0.025;
        let run = |state: &GasRocketState| {
            state
                .advance_effusion_with_wall_adaptive(
                    cx,
                    free(),
                    &outlets,
                    input,
                    &body(100.0),
                    0.0,
                    solve(),
                    accuracy,
                )
                .unwrap()
        };
        let a = run(&at_rest);
        let b = run(&moving);
        assert_eq!(a.wall(), b.wall());
        assert_eq!(a.endpoint().reservoir(), b.endpoint().reservoir());
        assert_eq!(a.step().accepted_intervals, b.step().accepted_intervals);
        vector_near(
            b.endpoint()
                .vehicle()
                .velocity_m_s
                .sub(a.endpoint().vehicle().velocity_m_s),
            boost,
            1e-12,
        );
        vector_near(
            b.endpoint()
                .vehicle()
                .position_m
                .sub(a.endpoint().vehicle().position_m),
            boost.scale(DT),
            1e-12,
        );
        vector_near(
            b.step().escaped_momentum_world_kg_m_s,
            a.step()
                .escaped_momentum_world_kg_m_s
                .add(boost.scale(a.step().escaped_mass_kg)),
            1e-12,
        );
        let correction = boost.dot(a.step().escaped_momentum_world_kg_m_s)
            + 0.5 * boost.norm_squared() * a.step().escaped_mass_kg;
        near(
            b.step().escaped_energy_world_j,
            a.step().escaped_energy_world_j + correction,
            1e-10,
        );
        near(
            b.step().escaped_reference_energy_world_j,
            a.step().escaped_reference_energy_world_j + correction,
            1e-10,
        );
        check_window(&moving, input, &b);
    });
}

#[test]
fn exported_energy_controls_refinement_and_window_budget_never_multiplies_with_steps() {
    scope(|cx, _| {
        let mut before = initial(PURE, Vec3::new(1.0, 0.0, 0.0));
        // A representable native inventory discrepancy is already permitted by
        // the original owner; its positive uncertainty reserve must also spend
        // the global budget, even if signed algebraic residuals happen to vanish.
        before.vehicle.propellant_mass_kg = before.vehicle.propellant_mass_kg.next_up();
        let saved = before.clone();
        let chart = curve([1000.0; 2]);
        let input = UniformEnthalpyStepInput {
            duration_s: 0.2,
            ..wall(&chart, 600.0)
        };
        let law = interface(before.reservoir().model());
        let outlets = [outlet(&law, Vec3::new(1.0, 0.0, 0.0))];
        let mut loose = adaptive(1e12);
        loose.motion.relative_tolerance = 0.0;
        let mut tight = loose;
        tight.emitted_energy_tolerance_j = 1e-5;
        let run = |accuracy| {
            before.advance_effusion_with_wall_adaptive(
                cx,
                free(),
                &outlets,
                input,
                &body(100.0),
                0.0,
                solve(),
                accuracy,
            )
        };
        let coarse = run(loose).unwrap();
        let fine = run(tight).unwrap();
        assert!(fine.step().accepted_intervals > coarse.step().accepted_intervals);
        assert!(fine.step().debit.energy_j > 0.0);
        tight.window_energy_tolerance_j = 0.5 * fine.step().debit.energy_j;
        assert!(matches!(
            run(tight),
            Err(EffusionWallTrajectoryError::WindowBalance { .. })
        ));
        assert_eq!(before, saved);
        check_window(&before, input, &fine);

        let mut resting = initial(PURE, Vec3::ZERO);
        resting.vehicle.propellant_mass_kg = resting.vehicle.propellant_mass_kg.next_up();
        let saved_resting = resting.clone();
        let mut exact_mass = adaptive(1e6);
        exact_mass.window_mass_tolerance_kg = 0.0;
        assert!(matches!(resting.advance_effusion_with_wall_adaptive(
            cx, free(), &[], wall(&chart, T0), &body(0.0), 0.0, solve(), exact_mass,
        ), Err(EffusionWallTrajectoryError::WindowBalance { residual, debit, .. })
            if residual.mass_kg == 0.0 && debit.mass_kg > 0.0));
        assert_eq!(resting, saved_resting);
    });
}

#[test]
fn represented_wall_datum_shifts_do_not_loosen_the_increment_error_scale() {
    scope(|cx, _| {
        let before = initial(PURE, Vec3::ZERO);
        let chart = curve([1000.0; 2]);
        let offset = 67_108_864.0;
        let shifted = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([42; 32]),
            SolidLiquidPhase::Solid,
            chart
                .knots()
                .iter()
                .map(|knot| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: knot.specific_enthalpy_j_kg + offset,
                    ..*knot
                })
                .collect(),
        )
        .unwrap();
        let input = wall(&chart, 600.0);
        let shifted_input = UniformEnthalpyStepInput {
            curve: &shifted,
            initial: shifted
                .state_at_specific_enthalpy(input.initial.specific_enthalpy_j_kg() + offset)
                .unwrap(),
            ..input
        };
        // The same absolute body gate covers the shifted chart's representable
        // enthalpy spacing. No tolerance is scaled by the huge datum itself.
        let mut boundary = body(100.0);
        boundary.enthalpy_tolerance_j_kg = 1e-6;
        boundary.maximum_thermal_residual_j = 1e-8;
        let law = interface(before.reservoir().model());
        let outlets = [outlet(&law, Vec3::new(1.0, 0.0, 0.0))];
        let run = |wall| {
            before
                .advance_effusion_with_wall_adaptive(
                    cx,
                    free(),
                    &outlets,
                    wall,
                    &boundary,
                    0.0,
                    solve(),
                    adaptive(1.0),
                )
                .unwrap()
        };
        let a = run(input);
        let b = run(shifted_input);
        assert_eq!(a.step().accepted_intervals, b.step().accepted_intervals);
        near(a.wall().temperature_k(), b.wall().temperature_k(), 1e-5);
        near(
            a.endpoint().reservoir().state().state().temperature,
            b.endpoint().reservoir().state().state().temperature,
            1e-5,
        );
        near(
            a.step().wall_storage_change_j,
            b.step().wall_storage_change_j,
            1e-6,
        );
        near(
            a.step().escaped_reference_energy_world_j,
            b.step().escaped_reference_energy_world_j,
            1e-7,
        );
        check_window(&before, shifted_input, &b);
    });
}

#[test]
fn withdrawal_retry_and_typed_work_minimum_source_cancellation_refusals_are_atomic() {
    scope(|cx, gate| {
        let before = initial(PURE, Vec3::ZERO);
        let saved = before.clone();
        let chart = curve([1000.0; 2]);
        let input = UniformEnthalpyStepInput {
            duration_s: 0.2,
            ..wall(&chart, 600.0)
        };
        let law = interface(before.reservoir().model());
        let outlets = [outlet(&law, Vec3::new(1.0, 0.0, 0.0))];
        let mut bounded = solve();
        bounded.recoil.transport.exchange.maximum_species_fraction = 2e-4;
        let retried = before
            .advance_effusion_with_wall_adaptive(
                cx,
                free(),
                &outlets,
                input,
                &body(100.0),
                0.0,
                bounded,
                adaptive(1e6),
            )
            .unwrap();
        assert!(retried.step().rejected_intervals >= 3);
        check_window(&before, input, &retried);

        let mut too_tight = adaptive(1e-7);
        too_tight.motion.minimum_step_s = too_tight.motion.initial_step_s;
        assert!(matches!(
            before.advance_effusion_with_wall_adaptive(
                cx,
                free(),
                &outlets,
                input,
                &body(100.0),
                0.0,
                solve(),
                too_tight,
            ),
            Err(EffusionWallTrajectoryError::Control(
                EffusionTrajectoryError::MinimumStep { .. }
            ))
        ));
        let mut too_little = adaptive(1e-7);
        too_little.motion.maximum_attempts = 1;
        assert!(matches!(
            before.advance_effusion_with_wall_adaptive(
                cx,
                free(),
                &outlets,
                input,
                &body(100.0),
                0.0,
                solve(),
                too_little,
            ),
            Err(EffusionWallTrajectoryError::Control(
                EffusionTrajectoryError::WorkBudget { .. }
            ))
        ));
        assert!(matches!(
            before.advance_effusion_with_wall_convection_adaptive(
                cx,
                free(),
                &outlets,
                input,
                &body(100.0),
                0.0,
                solve(),
                adaptive(1.0),
                |_| Err::<ConvectionCoefficient<()>, _>("unsupported heat law"),
            ),
            Err(EffusionWallTrajectoryError::Step {
                error: EffusionWallError::Thermal(ConvectionError::Convection(
                    "unsupported heat law"
                )),
                ..
            })
        ));

        let mut calls = 0;
        let cancelled = before.advance_effusion_with_wall_convection_adaptive(
            cx,
            free(),
            &outlets,
            input,
            &body(0.0),
            0.0,
            solve(),
            adaptive(1.0),
            |gas| {
                calls += 1;
                if calls > 64 {
                    gate.request();
                }
                Ok::<_, core::convert::Infallible>(ConvectionCoefficient {
                    coefficient_w_per_m2_k: gas.state().temperature / 4.0,
                    source: gas.state().temperature,
                })
            },
        );
        assert!(calls > 64);
        assert!(matches!(
            cancelled,
            Err(EffusionWallTrajectoryError::Cancelled)
        ));
        assert_eq!(before, saved);
        assert_eq!(input.initial.temperature_k(), 600.0);
    });
}
