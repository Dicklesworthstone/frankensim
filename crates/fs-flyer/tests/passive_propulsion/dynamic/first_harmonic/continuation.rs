//! Continued native wing branches coupled to actual independent shaft/wake motion.
use super::*;
use fs_wing::nonlinear::continuation::{ContinuationControls, WingContinuation};

fn domain() -> ContinuationControls {
    ContinuationControls {
        maximum_incidence_rad: 0.6,
        maximum_incidence_change_rad: 0.08,
        maximum_row_condition: 16.0,
        minimum_trailing_clearance_chords: 0.02,
        maximum_newton_iterations: 32,
        maximum_residual_evaluations: 600,
        maximum_jacobian_products: 4096,
        section_difference_rad: 1e-6,
    }
}

struct Trajectory {
    state: HarmonicShaftState<2>,
    wing: WingContinuation,
    coupling_work_j: f64,
}

fn integrate(skew: bool, count: usize) -> Trajectory {
    let design = wright_openloop_v1();
    let a = assembly(design.gross_mass_kg, CANTED);
    let base = ShaftPropulsion::new_spatial(&design, &a, mounts(CANTED)).unwrap();
    let model = model(&base, skew);
    let dt = 0.048 / count as f64;
    let origin = [0.3, -0.4, 0.5];
    let mut motion = a
        .shaft_motion(&initial(&a, [32.0, -39.0]), properties(), dt)
        .unwrap();
    motion
        .apply_impulse(
            RotorImpulses {
                carrier: CarrierImpulse::default(),
                motor_nms: [0.25, -0.125],
                aerodynamic_nms: [-0.0625, 0.03125],
            },
            origin,
        )
        .unwrap();
    let mut induction = modes();
    let start = HarmonicShaftState {
        motion: *motion.state(),
        induction,
    };
    let canonical = *motion.native_state();
    let mut wing = model
        .initialize_continuation(
            &start,
            control(&a, 0.0, &start),
            &|p| wind(0.0, p),
            domain(),
        )
        .unwrap();
    assert_eq!(motion.native_state(), &canonical);
    assert_eq!(motion.state(), &start.motion);
    assert_eq!(induction, start.induction);
    let initial_gamma = wing.gamma().to_vec();
    let (p0, l0) = momentum(&a, motion.state(), origin);
    let mut total_linear = Vec3::default();
    let mut total_angular = Vec3::default();
    let mut coupling_work_j = 0.0;
    let mut maximum_mechanical_work = 0.0_f64;
    let mut maximum_modal_momentum = 0.0_f64;
    let mut maximum_modal_energy = 0.0_f64;
    let mut maximum_newton = 0;
    let radius = design.rotor.radius_m;
    let masses = [
        (8.0 / 3.0) * RHO * radius.powi(3),
        (16.0 / 45.0) * RHO * radius.powi(3),
        (16.0 / 45.0) * RHO * radius.powi(3),
    ];
    for tick in 0..count {
        let time = tick as f64 * dt;
        let before = *motion.state();
        let before_modes = induction;
        let before_wing = wing.clone();
        let mut stages = Vec::new();
        let step = model
            .step_continued(
                &mut motion,
                &mut induction,
                &mut wing,
                domain(),
                time,
                origin,
                |t, state| {
                    stages.push((t, *state));
                    Ok(control(&a, t, state))
                },
                wind,
            )
            .unwrap_or_else(|error| panic!("skew={skew} count={count} tick={tick}: {error:?}"));
        assert_eq!(stages.len(), 2);
        assert_eq!(stages[0].1.motion, before);
        assert_eq!(stages[0].1.induction, before_modes);
        assert_eq!(stages[1].0, time + 0.5 * dt);
        assert_eq!(
            step.initial,
            model
                .evaluate_continued(
                    &stages[0].1,
                    &before_wing,
                    control(&a, stages[0].0, &stages[0].1),
                    &|p| wind(stages[0].0, p),
                    domain(),
                )
                .unwrap()
        );
        let first_anchor = step.initial.coupled.wing_continuation.as_ref().unwrap();
        assert_eq!(
            step.midpoint,
            model
                .evaluate_continued(
                    &stages[1].1,
                    first_anchor,
                    control(&a, stages[1].0, &stages[1].1),
                    &|p| wind(stages[1].0, p),
                    domain(),
                )
                .unwrap()
        );
        assert_eq!(
            wing,
            *step.midpoint.coupled.wing_continuation.as_ref().unwrap()
        );
        assert_eq!(motion.state(), &step.state.motion);
        assert_eq!(
            motion.native_state(),
            &step.motion.native.motion.state_after
        );
        assert_eq!(induction, step.state.induction);
        for evaluation in [&step.initial, &step.midpoint] {
            let anchor = evaluation.coupled.wing_continuation.as_ref().unwrap();
            assert_eq!(anchor.gamma(), evaluation.coupled.gamma);
            assert!(evaluation.coupled.wing_residual < fs_wing::nonlinear::GAMMA_TOL);
            assert!(anchor.maximum_incidence_change_rad() <= domain().maximum_incidence_change_rad);
            assert!(anchor.maximum_row_condition() <= domain().maximum_row_condition);
            assert!(
                anchor.minimum_trailing_clearance_chords()
                    >= domain().minimum_trailing_clearance_chords
            );
            for (actual, expected) in anchor
                .effective_incidence_rad()
                .iter()
                .zip(&evaluation.coupled.strip_alpha_eff_rad)
            {
                assert!((actual - expected).abs() < 2e-13);
            }
            maximum_newton = maximum_newton.max(evaluation.coupled.wing_iterations);
        }

        // Reconstruct physical external impulse from the accepted loads, with
        // each actual rotor's axial aerodynamic torque exactly once.
        let force = v(step.midpoint.carrier.force_n).scale(dt);
        let mut torque = v(step.midpoint.carrier.moment_nm);
        for k in 0..2 {
            torque =
                torque + v(a.rotors[k].axis_body).scale(step.midpoint.torques.aerodynamic_nm[k]);
        }
        let angular = (rotation(&before).rotate(torque).unwrap()
            + rotation(motion.state()).rotate(torque).unwrap())
        .scale(0.5 * dt)
            + (v(before.rigid.pos_m) + v(before.rigid.vel_mps).scale(0.5 * dt) - v(origin))
                .cross(force);
        total_linear = total_linear + force;
        total_angular = total_angular + angular;
        let motor_work: f64 = (0..2)
            .map(|k| {
                step.midpoint.torques.motor_nm[k]
                    * (motion.state().relative_phase_rad[k] - before.relative_phase_rad[k])
            })
            .sum();
        let heat: f64 = step.motion.native.dissipated_energy_j.into_iter().sum();
        assert!(heat > 0.0);
        assert!((step.motion.motor_work_j - motor_work).abs() < 1e-10);
        let balance = energy(&a, motion.state()) - energy(&a, &before) + heat
            - motor_work
            - step.motion.external_work_j
            - step.motion.native.motion.drift_energy_defect_j;
        maximum_mechanical_work = maximum_mechanical_work.max(balance.abs());
        for k in 0..2 {
            let u0 = before_modes[k].modes();
            let u1 = induction[k].modes();
            let receipt = step.inflow[k];
            let mut kinetic_change = 0.0;
            for j in 0..3 {
                let balance = masses[j] * (u1[j] - u0[j])
                    - dt * step.midpoint.coupled.wake_force_n[k][j]
                    + receipt.transport_impulse_ns[j]
                    + receipt.frame_impulse_ns[j];
                maximum_modal_momentum = maximum_modal_momentum.max(balance.abs());
                kinetic_change += 0.5 * masses[j] * (u1[j].powi(2) - u0[j].powi(2));
            }
            let balance = kinetic_change - receipt.induced_work_j
                + receipt.transport_work_j
                + receipt.frame_work_j;
            maximum_modal_energy = maximum_modal_energy.max(balance.abs());
            assert_eq!(
                step.coupling_work_residual_j[k],
                receipt.induced_work_j - dt * step.midpoint.coupled.propellers[k].induced_work_w
            );
            coupling_work_j += step.coupling_work_residual_j[k].abs();
        }
        assert_eq!(step.skew_inflow.is_some(), skew);
    }
    let (p1, l1) = momentum(&a, motion.state(), origin);
    let linear_defect = (p1 - p0 - total_linear).norm();
    let angular_defect = (l1 - l0 - total_angular).norm();
    println!(
        "continued independent skew={skew}, n={count}: P={linear_defect:e}, L={angular_defect:e}, mechanical work={maximum_mechanical_work:e}, modal momentum={maximum_modal_momentum:e}, modal energy={maximum_modal_energy:e}, Newton={maximum_newton}"
    );
    assert!(linear_defect < 2e-9 && angular_defect < 2e-8);
    assert!(maximum_mechanical_work < 1e-8);
    assert!(maximum_modal_momentum < 2e-9 && maximum_modal_energy < 2e-9);
    assert_ne!(wing.gamma(), initial_gamma);
    assert_ne!(induction, start.induction);
    Trajectory {
        state: HarmonicShaftState {
            motion: *motion.state(),
            induction,
        },
        wing,
        coupling_work_j,
    }
}

#[test]
fn g1_normal_and_skew_continued_trajectories_retain_actual_branch_and_source_accounts() {
    for skew in [false, true] {
        let result = integrate(skew, 12);
        assert!(result.wing.residual_evaluations() > 0);
        assert!(result.coupling_work_j > 0.0);
    }
}

#[test]
fn g1_continued_independent_body_shaft_and_modal_trajectories_refine_together() {
    for skew in [false, true] {
        let runs = [6, 12, 24, 48].map(|n| integrate(skew, n));
        let groups = |state: &HarmonicShaftState<2>| {
            [
                state
                    .motion
                    .rigid
                    .pos_m
                    .into_iter()
                    .chain(state.motion.rigid.vel_mps)
                    .chain(state.motion.rigid.omega_body)
                    .chain(state.motion.rigid.quat)
                    .collect::<Vec<_>>(),
                state
                    .motion
                    .relative_phase_rad
                    .into_iter()
                    .chain([
                        state.motion.axial_momentum[0] / 2.0,
                        state.motion.axial_momentum[1] / 4.0,
                    ])
                    .collect(),
                state
                    .induction
                    .into_iter()
                    .flat_map(|u| u.modes())
                    .collect(),
            ]
        };
        for group in 0..3 {
            let difference: Vec<_> = runs
                .windows(2)
                .map(|pair| {
                    groups(&pair[0].state)[group]
                        .iter()
                        .zip(&groups(&pair[1].state)[group])
                        .map(|(a, b)| (a - b).powi(2))
                        .sum::<f64>()
                        .sqrt()
                })
                .collect();
            println!("continued independent skew={skew} group={group} differences: {difference:?}");
            for pair in difference.windows(2) {
                assert!((3.3..4.8).contains(&(pair[0] / pair[1])), "{difference:?}");
            }
        }
        let work = runs.each_ref().map(|run| run.coupling_work_j);
        println!("continued independent skew={skew} coupling work: {work:?}");
        for pair in work.windows(2) {
            assert!((3.3..4.8).contains(&(pair[0] / pair[1])), "{work:?}");
        }
    }
}

#[test]
fn g0_late_branch_source_and_native_refusals_preserve_modes_momentum_and_wing_anchor() {
    let design = wright_openloop_v1();
    let a = assembly(design.gross_mass_kg, CANTED);
    let base = ShaftPropulsion::new_spatial(&design, &a, mounts(CANTED)).unwrap();
    let start = initial(&a, [32.0, -39.0]);
    for skew in [false, true] {
        let model = model(&base, skew);
        for failure in 0..3 {
            let mut motion = a.shaft_motion(&start, properties(), 0.004).unwrap();
            let before = *motion.state();
            let canonical = *motion.native_state();
            let mut induction = modes();
            let state = HarmonicShaftState {
                motion: before,
                induction,
            };
            let mut wing = model
                .initialize_continuation(
                    &state,
                    control(&a, 0.0, &state),
                    &|p| wind(0.0, p),
                    domain(),
                )
                .unwrap();
            let saved_wing = wing.clone();
            let mut calls = 0;
            let error = model
                .step_continued(
                    &mut motion,
                    &mut induction,
                    &mut wing,
                    domain(),
                    0.0,
                    ZERO,
                    |time, state| {
                        calls += 1;
                        let mut values = control(&a, time, state);
                        if time > 0.0 {
                            if failure == 0 {
                                values.canard_rad = 0.4;
                            }
                            if failure == 2 {
                                values.motor_nm[0] = f64::MAX / 32.0;
                            }
                        }
                        Ok(values)
                    },
                    |time, point| {
                        if failure == 1 && time > 0.0 {
                            Err(fs_flyer::Refusal {
                                code: "wind-domain",
                                message: "continued midpoint material air unavailable".into(),
                                ranked_repairs: vec![],
                            })
                        } else {
                            wind(time, point)
                        }
                    },
                )
                .unwrap_err();
            assert_eq!(calls, 2);
            assert_eq!(
                error.code,
                [
                    "nonlinear-continuation-outside-domain",
                    "wind-domain",
                    "rotor-shaft-refused",
                ][failure]
            );
            assert_eq!(motion.state(), &before);
            assert_eq!(motion.native_state(), &canonical);
            assert_eq!(induction, modes());
            assert_eq!(wing, saved_wing);
            let mut fresh = a.shaft_motion(&start, properties(), 0.004).unwrap();
            let mut fresh_modes = modes();
            let mut fresh_wing = saved_wing;
            let retry = model
                .step_continued(
                    &mut motion,
                    &mut induction,
                    &mut wing,
                    domain(),
                    0.0,
                    ZERO,
                    |t, state| Ok(control(&a, t, state)),
                    wind,
                )
                .unwrap();
            let expected = model
                .step_continued(
                    &mut fresh,
                    &mut fresh_modes,
                    &mut fresh_wing,
                    domain(),
                    0.0,
                    ZERO,
                    |t, state| Ok(control(&a, t, state)),
                    wind,
                )
                .unwrap();
            assert_eq!(retry, expected);
            assert_eq!(motion.native_state(), fresh.native_state());
            assert_eq!(induction, fresh_modes);
            assert_eq!(wing, fresh_wing);
        }
    }
}
