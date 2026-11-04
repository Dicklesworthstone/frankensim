//! Native continued wing branches under steady annular angular-momentum wash.
use super::*;
use fs_flyer::propcoupling::swirl::wake::ContinuedWakeStep;
use fs_wing::nonlinear::continuation::ContinuationControls;

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
fn turning_start(a: &RotorAssembly<2>, rates: [f64; 2]) -> RotorState<2> {
    let mut state = start(a, 12.0, rates);
    state.rigid.vel_mps = [12.0 * 0.1_f64.cos(), 0.0, 12.0 * 0.1_f64.sin()];
    state
}
fn gust(t: f64) -> [f64; 3] {
    [0.0, 0.0, -0.05 - 4.0 * t]
}
fn drive(t: f64) -> ShaftInputs<2> {
    ShaftInputs {
        motor_nm: [40.0 + 20.0 * t, -55.0 - 10.0 * t],
        ..control()
    }
}
fn exact_anchor(flow: &ContinuedWakeStep) {
    let wing = &flow.solution.solution.coupled;
    assert_eq!(flow.wing_continuation.gamma(), wing.gamma);
    assert_eq!(
        flow.wing_continuation.effective_incidence_rad(),
        wing.strip_alpha_eff_rad
    );
    assert!(
        flow.wing_continuation.maximum_incidence_change_rad()
            <= domain().maximum_incidence_change_rad
    );
}

#[test]
fn g1_changing_annular_loads_advance_exact_force_branch_and_physical_momentum() {
    let (d, panels, strips, sections, groups) = fixture();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &sections, &groups).unwrap();
    let a = assembly(d.gross_mass_kg, AXIAL);
    let base = ShaftPropulsion::new_spatial(&d, &a, mounts(AXIAL))
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    let dt = 0.004;
    let origin = [0.3, -0.4, 0.5];
    for transport in TRANSPORTS {
        let model = WakeShaftPropulsion::new(&base, Controls::default(), WakeControls::default())
            .with_transport(transport);
        let mut motion = a
            .shaft_motion(
                &turning_start(&a, [32.0, -39.0]),
                ShaftProperties::default(),
                dt,
            )
            .unwrap();
        let first = model
            .evaluate_continued(
                motion.state(),
                drive(0.0),
                &|_| Ok(gust(0.0)),
                None,
                None,
                domain(),
            )
            .unwrap();
        exact_anchor(&first.coupled);
        let mut wing = first.coupled.wing_continuation;
        let initial_gamma = wing.gamma().to_vec();
        let mut warm = first.coupled.solution.warm_start();
        let initial = *motion.state();
        let (p0, l0) = momentum(&a, &initial, origin);
        let mut impulse = Vec3::default();
        let mut angular = Vec3::default();
        let mut motor_work = 0.0;
        let mut external_work = 0.0;
        let mut drift = 0.0;
        let mut max_newton = 0;
        for tick in 0..20 {
            let t = tick as f64 * dt;
            let before = *motion.state();
            let old_wing = wing.clone();
            let old_warm = warm.clone();
            let mut stages = Vec::new();
            let step = model
                .step_continued(
                    &mut motion,
                    &mut wing,
                    t,
                    origin,
                    Some(&warm),
                    domain(),
                    |time, state| {
                        stages.push((time, *state));
                        Ok(drive(time))
                    },
                    |time, _| Ok(gust(time)),
                )
                .unwrap();
            assert_eq!(stages.len(), 2);
            assert_eq!(stages[0], (t, before));
            assert_ne!(stages[1].1.axial_momentum, before.axial_momentum);
            assert_eq!(
                step.initial,
                model
                    .evaluate_continued(
                        &before,
                        drive(t),
                        &|_| Ok(gust(t)),
                        Some(&old_warm),
                        Some(&old_wing),
                        domain()
                    )
                    .unwrap()
            );
            assert_eq!(
                step.midpoint,
                model
                    .evaluate_continued(
                        &stages[1].1,
                        drive(stages[1].0),
                        &|_| Ok(gust(stages[1].0)),
                        Some(&step.initial.coupled.solution.warm_start()),
                        Some(&step.initial.coupled.wing_continuation),
                        domain()
                    )
                    .unwrap()
            );
            for stage in [&step.initial, &step.midpoint] {
                exact_anchor(&stage.coupled);
                for k in 0..2 {
                    let flux = stage.coupled.solution.wakes[k].cross_section(0.8).unwrap();
                    near(
                        flux.angular_momentum_flux_nm,
                        -stage.torques.aerodynamic_nm[k],
                        2e-13,
                    );
                }
            }
            assert_eq!(wing, step.midpoint.coupled.wing_continuation);
            assert_eq!(
                motion.native_state(),
                &step.motion.native.motion.state_after
            );
            max_newton = max_newton.max(wing.jacobian_products());
            warm = step.midpoint.coupled.solution.warm_start();
            impulse = impulse + v(step.motion.external_linear_impulse_world);
            angular = angular + v(step.motion.external_angular_impulse_world);
            motor_work += step.motion.motor_work_j;
            external_work += step.motion.external_work_j;
            drift += step.motion.native.motion.drift_energy_defect_j;
            assert_eq!(step.motion.native.dissipated_energy_j, [0.0; 2]);
        }
        let (p1, l1) = momentum(&a, motion.state(), origin);
        assert!((p1 - p0 - impulse).norm() < 1e-9);
        assert!((l1 - l0 - angular).norm() < 1e-8);
        let work = kinetic(&a, motion.state())
            - kinetic(&a, &initial)
            - motor_work
            - external_work
            - drift;
        assert!(work.abs() < 1e-8, "continued physical work {work:e}");
        let gamma_change = wing
            .gamma()
            .iter()
            .zip(initial_gamma)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        assert!(gamma_change > 1e-3 && max_newton > 0);
        eprintln!(
            "{transport:?}: nonlinear gamma change {gamma_change:e}, max native JVPs {max_newton}, angular residual {:e}, work residual {work:e}",
            (l1 - l0 - angular).norm()
        );
    }
}

#[test]
fn g0_fixed_caller_branch_limits_apply_through_all_nested_field_trials() {
    let (d, panels, strips, sections, groups) = fixture();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &sections, &groups).unwrap();
    let a = assembly(d.gross_mass_kg, AXIAL);
    let base = ShaftPropulsion::new_spatial(&d, &a, mounts(AXIAL))
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    for transport in TRANSPORTS {
        let model = WakeShaftPropulsion::new(&base, Controls::default(), WakeControls::default())
            .with_transport(transport);
        let initial = turning_start(&a, [32.0, -39.0]);
        let first = model
            .evaluate_continued(&initial, control(), &|_| Ok(ZERO), None, None, domain())
            .unwrap();
        let wing = first.coupled.wing_continuation;
        let seed = first.coupled.solution.warm_start();
        let changed = turning_start(&a, [45.0, -50.0]);
        let accepted = model
            .evaluate_continued(
                &changed,
                control(),
                &|_| Ok(ZERO),
                Some(&seed),
                Some(&wing),
                domain(),
            )
            .unwrap();
        exact_anchor(&accepted.coupled);
        let delta = accepted
            .coupled
            .wing_continuation
            .maximum_incidence_change_rad();
        assert!(delta > 1e-4);
        let strict = ContinuationControls {
            maximum_incidence_change_rad: 0.5 * delta,
            ..domain()
        };
        let error = model
            .evaluate_continued(
                &changed,
                control(),
                &|_| Ok(ZERO),
                Some(&seed),
                Some(&wing),
                strict,
            )
            .unwrap_err();
        assert_eq!(error.code, "nonlinear-continuation-outside-domain");
        // This changing-load fixture also distinguishes branch continuation
        // from a fresh native Newton solve: the latter cannot resolve its
        // cold branch even though the caller-anchored solution above succeeds.
        if transport == WakeTransport::FixedCylinder {
            assert_eq!(
                model
                    .evaluate_continued(
                        &changed,
                        control(),
                        &|_| Ok(ZERO),
                        Some(&seed),
                        None,
                        strict
                    )
                    .unwrap_err()
                    .code,
                "nonlinear-continuation-unresolved"
            );
        }
        let repeated = model
            .evaluate_continued(
                &changed,
                control(),
                &|_| Ok(ZERO),
                Some(&accepted.coupled.solution.warm_start()),
                Some(&accepted.coupled.wing_continuation),
                strict,
            )
            .unwrap();
        exact_anchor(&repeated.coupled);
        assert_ne!(
            accepted.coupled.solution.solution.coupled.spec_digest,
            repeated.coupled.solution.solution.coupled.spec_digest
        );
        assert_eq!(wing.gamma(), first.coupled.solution.solution.coupled.gamma);
    }
}

#[test]
fn g0_midpoint_branch_and_native_failures_leave_both_states_anchor_and_seed_retryable() {
    let (d, panels, strips, sections, groups) = fixture();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &sections, &groups).unwrap();
    let a = assembly(d.gross_mass_kg, AXIAL);
    let base = ShaftPropulsion::new_spatial(&d, &a, mounts(AXIAL))
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    let start = turning_start(&a, [32.0, -39.0]);
    for transport in TRANSPORTS {
        let model = WakeShaftPropulsion::new(&base, Controls::default(), WakeControls::default())
            .with_transport(transport);
        let first = model
            .evaluate_continued(&start, control(), &|_| Ok(ZERO), None, None, domain())
            .unwrap();
        let mut wing = first.coupled.wing_continuation;
        let seed = first.coupled.solution.warm_start();
        let saved_wing = wing.clone();
        let saved_seed = seed.clone();
        let mut motion = a
            .shaft_motion(&start, ShaftProperties::default(), 0.004)
            .unwrap();
        let before = *motion.state();
        let canonical = *motion.native_state();
        for native in [false, true] {
            let mut calls = 0;
            let error = model
                .step_continued(
                    &mut motion,
                    &mut wing,
                    0.0,
                    ZERO,
                    Some(&seed),
                    domain(),
                    |_, _| {
                        calls += 1;
                        Ok(ShaftInputs {
                            motor_nm: if native && calls == 2 {
                                [f64::MAX / 32.0; 2]
                            } else {
                                control().motor_nm
                            },
                            ..control()
                        })
                    },
                    |t, _| {
                        Ok(if !native && t > 0.0 {
                            [0.0, 0.0, -3.0]
                        } else {
                            ZERO
                        })
                    },
                )
                .unwrap_err();
            assert_eq!(calls, 2);
            assert_eq!(
                error.code,
                if native {
                    "rotor-shaft-refused"
                } else {
                    "nonlinear-continuation-outside-domain"
                }
            );
            assert_eq!(motion.state(), &before);
            assert_eq!(motion.native_state(), &canonical);
            assert_eq!(wing, saved_wing);
            assert_eq!(seed, saved_seed);
        }
        let mut fresh = a
            .shaft_motion(&start, ShaftProperties::default(), 0.004)
            .unwrap();
        let mut fresh_wing = saved_wing.clone();
        let retry = model
            .step_continued(
                &mut motion,
                &mut wing,
                0.0,
                ZERO,
                Some(&seed),
                domain(),
                |_, _| Ok(control()),
                |_, _| Ok(ZERO),
            )
            .unwrap();
        let expected = model
            .step_continued(
                &mut fresh,
                &mut fresh_wing,
                0.0,
                ZERO,
                Some(&seed),
                domain(),
                |_, _| Ok(control()),
                |_, _| Ok(ZERO),
            )
            .unwrap();
        assert_eq!(retry, expected);
        assert_eq!(wing, fresh_wing);
        assert_eq!(motion.native_state(), fresh.native_state());
    }
}

#[test]
fn g3_geared_continued_stages_share_native_accepted_annular_branch() {
    let (d, panels, strips, sections, groups) = fixture();
    let surfaces = LiftingSurfaces::new(&panels, &strips, &sections, &groups).unwrap();
    let a = assembly(d.gross_mass_kg, AXIAL);
    let mut locked = a.reduced_body.inertia_tensor().matrix();
    locked.m[0] += a.rotors[0].axial_inertia_kgm2 + a.rotors[1].axial_inertia_kgm2;
    let gears: GearedAssembly<2> = GearedAssembly::new(
        d.gross_mass_kg,
        locked,
        core::array::from_fn(|i| GearedRotor {
            axis_body: AXIAL[i],
            axial_inertia_kgm2: a.rotors[i].axial_inertia_kgm2,
            ratio: [1.0, -1.0][i],
        }),
    )
    .unwrap();
    let p = ShaftPropulsion::new_spatial(&d, &a, mounts(AXIAL))
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    let g = GearedPropulsion::new_spatial(&d, &gears, mounts(AXIAL))
        .unwrap()
        .with_lifting_surfaces(&surfaces)
        .unwrap();
    let state = turning_start(&a, [35.0, -35.0]);
    let geared_state = gears.state_from_speed(state.rigid, 0.2, 35.0).unwrap();
    let inputs = PropulsionInputs {
        density_kg_m3: 1.225,
        canard_rad: 0.0,
        motor_nm: 80.0,
        gravity_world_mps2: ZERO,
    };
    for transport in TRANSPORTS {
        let model = WakePropulsion::new(&g, Controls::default(), WakeControls::default())
            .with_transport(transport);
        let independent =
            WakeShaftPropulsion::new(&p, Controls::default(), WakeControls::default())
                .with_transport(transport)
                .evaluate_continued(&state, control(), &|_| Ok(ZERO), None, None, domain())
                .unwrap();
        let initial = model
            .evaluate_continued(&geared_state, inputs, &|_| Ok(ZERO), None, None, domain())
            .unwrap();
        assert_eq!(initial.coupled, independent.coupled);
        assert_eq!(initial.loads.carrier, independent.carrier);
        let wing = initial.coupled.wing_continuation;
        let warm = initial.coupled.solution.warm_start();
        let mut stages = Vec::new();
        let step = model
            .step_continued(
                &geared_state,
                &wing,
                0.0,
                0.004,
                ZERO,
                Some(&warm),
                domain(),
                |t, s| {
                    stages.push((t, *s));
                    Ok(inputs)
                },
                |t, _| Ok([0.0, 0.0, -4.0 * t]),
            )
            .unwrap();
        assert_eq!(stages.len(), 2);
        assert_eq!(
            step.midpoint,
            model
                .evaluate_continued(
                    &stages[1].1,
                    inputs,
                    &|_| Ok([0.0, 0.0, -4.0 * stages[1].0]),
                    Some(&step.initial.coupled.solution.warm_start()),
                    Some(&step.initial.coupled.wing_continuation),
                    domain()
                )
                .unwrap()
        );
        exact_anchor(&step.initial.coupled);
        exact_anchor(&step.midpoint.coupled);
        assert!(v(step.motion.angular_momentum_residual).norm() < 1e-9);
    }
}
