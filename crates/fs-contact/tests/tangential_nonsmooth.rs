//! Physical checks for coupled nonsmooth finite-patch tangential motion.
use fs_contact::tangential::nonsmooth::*;
use fs_exec::CancelGate;
use fs_solver::NewtonKrylovConfig;
use fs_tribo::partial_slip::*;

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 2e-8 * a.abs().max(b.abs()).max(1.0),
        "{a:e} != {b:e}"
    );
}
fn controls() -> NonsmoothControls {
    NonsmoothControls {
        maximum_evaluations: 4000,
        maximum_newton_iterations: 40,
        maximum_mechanical_dofs: 32,
        difference_step_mps: 1e-5,
        velocity_tolerance_mps: 1e-9,
        cone_tolerance: 1e-10,
        energy_tolerance_j: 1e-8,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-12,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn law() -> PartialSlipLaw {
    PartialSlipLaw::new(
        PARTIAL_SLIP_MODEL_ID,
        "synthetic/elastic-coulomb",
        PartialSlipParameters {
            static_mu: 0.8,
            kinetic_mu: 0.4,
            tangential_stiffness_n_per_m: 10000.0,
            torsional_stiffness_nm_per_rad: 100.0,
            torsional_capacity_factor: 0.5,
            partial_slip_onset_fraction: 0.5,
            partial_slip_hardening_fraction: 0.4,
        },
    )
    .unwrap()
}
fn patch(id: &str, normal: f64) -> NormalPatchView {
    NormalPatchView::new(
        id,
        "positive-normal-fixture",
        "synthetic/normal",
        NormalPatchAuthority::SyntheticFixture,
        normal,
        0.02,
        0.01,
        1e-4,
    )
    .unwrap()
}
fn interface() -> PartialSlipInterface {
    PartialSlipInterface::new(
        "a->b",
        "run-history",
        "synthetic/interface",
        NormalPatchAuthority::SyntheticFixture,
    )
    .unwrap()
}
fn work(patch: &str, interval: &str) -> GeneralizedWorkOwnership {
    GeneralizedWorkOwnership::new(patch, interval, "x", "y", "spin").unwrap()
}
fn contact<'a>(
    law: &'a PartialSlipLaw,
    patch: &'a NormalPatchView,
    interface: &'a PartialSlipInterface,
    interval: &str,
) -> NonsmoothContact<'a> {
    NonsmoothContact {
        law,
        patch,
        interface,
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
        common_spin_rad_per_s: 0.0,
        torsional_length_m: 0.01,
        creepage_reference_speed_mps: 1.0,
        work: work(patch.patch_id(), interval),
    }
}

#[test]
fn massive_body_stick_partial_and_gross_slide_match_independent_solutions() {
    let l = law();
    let p = patch("patch", 100.0);
    let i = interface();
    let gate = CancelGate::new();
    let mass: f64 = 2.0;
    let dt = 0.001;
    let factor = [1.0 / mass.sqrt(), 0.0, 0.0];
    for (free, branch, expected) in [
        (
            0.1,
            PartialSlipStateKind::Sticking,
            0.1 / (1.0 + dt * dt * 10000.0 / mass),
        ),
        (
            8.0,
            PartialSlipStateKind::PartialSlip,
            (8.0 - dt * 12.0 / mass) / (1.0 + dt * dt * 4000.0 / mass),
        ),
        (
            100.0,
            PartialSlipStateKind::GrossSlide,
            100.0 - dt * 40.0 / mass,
        ),
    ] {
        let state = NonsmoothTangentialState::new(1, 8).unwrap();
        let trial = state
            .prepare(
                &[contact(&l, &p, &i, "one")],
                &[[free, 0.0, 0.0]],
                TangentialMobility {
                    mechanical_dofs: 1,
                    factor: &factor,
                },
                dt,
                controls(),
                &gate,
            )
            .unwrap();
        let result = &trial.contacts()[0];
        assert_eq!(result.law.state, branch);
        close(result.endpoint_velocity_mps[0], expected);
        close(result.impulse_ns[0], mass * (expected - free));
        close(
            trial.audit().mechanical_work_j,
            0.5 * mass * (expected * expected - free * free),
        );
        close(trial.audit().energy_residual_j, 0.0);
        assert!(trial.audit().cone_violation <= 1e-10);
        assert!(trial.audit().complementarity_residual_j <= 1e-8);
        for j in 0..3 {
            close(result.action.force_n[j] + result.reaction.force_n[j], 0.0);
        }
        if branch == PartialSlipStateKind::Sticking {
            close(trial.audit().physical_heat_j, 0.0);
        } else {
            assert!(trial.audit().physical_heat_j > 0.0);
        }
        assert_eq!(
            state.version(),
            0,
            "prepare must not advance accepted state"
        );
    }
}

#[test]
fn simultaneous_contacts_share_body_recoil_and_are_permutation_invariant() {
    let l = law();
    let pa = patch("a", 100.0);
    let pb = patch("b", 50.0);
    let i = interface();
    let state = NonsmoothTangentialState::new(2, 8).unwrap();
    let gate = CancelGate::new();
    // Both contacts act on the same 1 kg translational DOF: W has all ones
    // on its x blocks, rather than two independent scalar inverse masses.
    let factor = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    for free in [0.02, 20.0] {
        let first = state
            .prepare(
                &[contact(&l, &pa, &i, "one"), contact(&l, &pb, &i, "one")],
                &[[free, 0.0, 0.0]; 2],
                TangentialMobility {
                    mechanical_dofs: 1,
                    factor: &factor,
                },
                0.001,
                controls(),
                &gate,
            )
            .unwrap();
        let reversed = state
            .prepare(
                &[contact(&l, &pb, &i, "one"), contact(&l, &pa, &i, "one")],
                &[[free, 0.0, 0.0]; 2],
                TangentialMobility {
                    mechanical_dofs: 1,
                    factor: &factor,
                },
                0.001,
                controls(),
                &gate,
            )
            .unwrap();
        assert_eq!(first.contacts(), reversed.contacts());
        assert_eq!(first.audit(), reversed.audit());
        let expected = if free < 1.0 {
            free / (1.0 + 0.001 * 0.001 * 20000.0)
        } else {
            free - 0.001 * (40.0 + 20.0)
        };
        for result in first.contacts() {
            close(result.endpoint_velocity_mps[0], expected);
        }
        let total_impulse: f64 = first.contacts().iter().map(|x| x.impulse_ns[0]).sum();
        close(total_impulse, expected - free);
        close(first.mass_weighted_velocity_increment()[0], total_impulse);
        close(
            first.audit().mechanical_work_j,
            0.5 * (expected * expected - free * free),
        );
        close(first.audit().energy_residual_j, 0.0);
    }
}

#[test]
fn mixed_slip_spin_closes_regional_plastic_work_and_survives_frame_gauge() {
    let l = law();
    let p = patch("patch", 100.0);
    let i = interface();
    let gate = CancelGate::new();
    let mut state = NonsmoothTangentialState::new(1, 8).unwrap();
    let preload = state
        .prepare(
            &[contact(&l, &p, &i, "preload")],
            &[[6.0, 2.0, 0.3]],
            TangentialMobility {
                mechanical_dofs: 0,
                factor: &[],
            },
            0.001,
            controls(),
            &gate,
        )
        .unwrap();
    state.commit(preload).unwrap();
    let fixed = state
        .prepare(
            &[contact(&l, &p, &i, "next")],
            &[[10.0, 3.0, 0.5]],
            TangentialMobility {
                mechanical_dofs: 3,
                factor: &[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
            },
            0.001,
            controls(),
            &gate,
        )
        .unwrap();
    let mut rotated_contact = contact(&l, &p, &i, "next");
    rotated_contact.frame = rotated_contact
        .frame
        .rotated(core::f64::consts::FRAC_PI_2)
        .unwrap();
    let rotated = state
        .prepare(
            &[rotated_contact],
            &[[3.0, -10.0, 0.5]],
            TangentialMobility {
                mechanical_dofs: 3,
                factor: &[0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            },
            0.001,
            controls(),
            &gate,
        )
        .unwrap();
    let a = &fixed.contacts()[0];
    let b = &rotated.contacts()[0];
    for j in 0..3 {
        close(a.action.force_n[j], b.action.force_n[j]);
        close(a.action.torque_nm[j], b.action.torque_nm[j]);
    }
    close(a.law.stored_energy_j, b.law.stored_energy_j);
    close(a.law.dissipation.heat_j, b.law.dissipation.heat_j);
    let mut heat = 0.0;
    for (region, weight) in [0.6, 0.4].into_iter().enumerate() {
        let elastic = a.law.next_state.region_elastic_displacement_m()[region];
        let twist = a.law.next_state.region_elastic_twist_rad()[region];
        let plastic = a.law.region_plastic_increment[region];
        heat += weight
            * (10000.0 * (elastic[0] * plastic[0] + elastic[1] * plastic[1])
                + 100.0 * twist * plastic[2]);
    }
    close(heat, a.law.dissipation.heat_j);
    close(fixed.audit().energy_residual_j, 0.0);
}

#[test]
fn rollback_checkpoint_release_and_recontact_keep_work_exactly_once() {
    let l = law();
    let p = patch("patch", 100.0);
    let i = interface();
    let gate = CancelGate::new();
    let mut state = NonsmoothTangentialState::new(1, 8).unwrap();
    let original = state.clone();
    let first = state
        .prepare(
            &[contact(&l, &p, &i, "first")],
            &[[6.0, 0.0, 0.0]],
            TangentialMobility {
                mechanical_dofs: 0,
                factor: &[],
            },
            0.001,
            controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(state, original);
    let abandoned = first.clone();
    drop(abandoned);
    let stale = first.clone();
    state.commit(first).unwrap();
    assert!(matches!(
        state.commit(stale),
        Err(NonsmoothTangentialError::StateMismatch)
    ));
    let restarted = state.clone();
    let second = state
        .prepare(
            &[contact(&l, &p, &i, "second")],
            &[[-0.2, 0.0, 0.0]],
            TangentialMobility {
                mechanical_dofs: 0,
                factor: &[],
            },
            0.001,
            controls(),
            &gate,
        )
        .unwrap();
    let replay = restarted
        .prepare(
            &[contact(&l, &p, &i, "second")],
            &[[-0.2, 0.0, 0.0]],
            TangentialMobility {
                mechanical_dofs: 0,
                factor: &[],
            },
            0.001,
            controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(second.contacts(), replay.contacts());
    assert!(
        second.contacts()[0].action.force_n[0] < 0.0,
        "reversal must first unload retained preload"
    );
    let stored = second.contacts()[0].law.stored_energy_j;
    state.commit(second).unwrap();
    let before = state.clone();
    assert!(
        state
            .release_contact(work("patch", "release"), 0.0, 1e-10)
            .is_err()
    );
    assert_eq!(state, before);
    close(
        state
            .release_contact(work("patch", "release"), stored, 1e-10)
            .unwrap(),
        stored,
    );
    assert!(state.contact_checkpoint("patch").is_none());
    assert!(matches!(
        state.prepare(
            &[contact(&l, &p, &i, "first")],
            &[[0.0; 3]],
            TangentialMobility {
                mechanical_dofs: 0,
                factor: &[]
            },
            0.001,
            controls(),
            &gate
        ),
        Err(NonsmoothTangentialError::DuplicateWork)
    ));
    let recontact = state
        .prepare(
            &[contact(&l, &p, &i, "recontact")],
            &[[0.0; 3]],
            TangentialMobility {
                mechanical_dofs: 0,
                factor: &[],
            },
            0.001,
            controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(recontact.contacts()[0].action.force_n, [0.0; 3]);
    close(recontact.contacts()[0].law.stored_energy_j, 0.0);
    state.commit(recontact).unwrap();
    close(state.time_s(), 0.003);
}

#[test]
fn budgets_cancellation_and_unconverged_coupling_never_advance_history() {
    let l = law();
    let p = patch("patch", 100.0);
    let i = interface();
    let state = NonsmoothTangentialState::new(1, 8).unwrap();
    let before = state.clone();
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        state.prepare(
            &[contact(&l, &p, &i, "one")],
            &[[1.0; 3]],
            TangentialMobility {
                mechanical_dofs: 0,
                factor: &[]
            },
            0.001,
            controls(),
            &gate
        ),
        Err(NonsmoothTangentialError::Cancelled)
    ));
    let gate = CancelGate::new();
    let factor = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
    let mut limited = controls();
    limited.maximum_evaluations = 1;
    assert!(matches!(
        state.prepare(
            &[contact(&l, &p, &i, "one")],
            &[[20.0, 10.0, 2.0]],
            TangentialMobility {
                mechanical_dofs: 3,
                factor: &factor
            },
            0.001,
            limited,
            &gate
        ),
        Err(NonsmoothTangentialError::EvaluationBudget)
    ));
    limited = controls();
    limited.maximum_newton_iterations = 1;
    limited.newton.absolute_tolerance = 1e-14;
    assert!(matches!(
        state.prepare(
            &[contact(&l, &p, &i, "one")],
            &[[20.0, 10.0, 2.0]],
            TangentialMobility {
                mechanical_dofs: 3,
                factor: &factor
            },
            0.001,
            limited,
            &gate
        ),
        Err(NonsmoothTangentialError::NotConverged(_))
    ));
    assert_eq!(state, before);
}

/// A physical free mass on the retained tangential spring converges to the
/// independent oscillator solution; integration loss never becomes heat.
#[test]
fn elastic_motion_refines_to_the_oscillator_without_fictitious_heat() {
    let l = law();
    let p = patch("patch", 100.0);
    let i = interface();
    let gate = CancelGate::new();
    let duration: f64 = 0.005;
    let initial_velocity: f64 = 0.05;
    let exact_velocity = initial_velocity * (100.0 * duration).cos();
    let exact_position = initial_velocity / 100.0 * (100.0 * duration).sin();
    let mut errors = Vec::new();
    let mut losses = Vec::new();
    for steps in [8, 16, 32] {
        let dt = duration / f64::from(steps);
        let mut state = NonsmoothTangentialState::new(1, steps as usize + 1).unwrap();
        let mut velocity = initial_velocity;
        let mut position = 0.0;
        let mut numerical_loss = 0.0;
        for step in 0..steps {
            let interval = format!("step-{step}");
            let trial = state
                .prepare(
                    &[contact(&l, &p, &i, &interval)],
                    &[[velocity, 0.0, 0.0]],
                    TangentialMobility {
                        mechanical_dofs: 1,
                        factor: &[1.0, 0.0, 0.0],
                    },
                    dt,
                    controls(),
                    &gate,
                )
                .unwrap();
            velocity = trial.contacts()[0].endpoint_velocity_mps[0];
            position += dt * velocity;
            close(trial.audit().physical_heat_j, 0.0);
            numerical_loss +=
                trial.audit().elastic_numerical_loss_j + trial.audit().mechanical_numerical_loss_j;
            state.commit(trial).unwrap();
        }
        let energy = 0.5 * velocity * velocity + 0.5 * 10000.0 * position * position;
        close(
            energy + numerical_loss,
            0.5 * initial_velocity * initial_velocity,
        );
        close(
            state
                .contact_checkpoint("patch")
                .unwrap()
                .state()
                .region_elastic_displacement_m()[0][0],
            position,
        );
        errors.push((velocity - exact_velocity).hypot(100.0 * (position - exact_position)));
        losses.push(numerical_loss);
    }
    for j in 1..errors.len() {
        assert!(
            errors[j] < 0.55 * errors[j - 1],
            "BE error refinement: {errors:?}"
        );
        assert!(
            losses[j] < 0.55 * losses[j - 1],
            "separate numerical loss refinement: {losses:?}"
        );
    }
}

#[test]
fn changing_coefficients_preserves_preload_and_exports_parameter_work() {
    let l = law();
    let p = patch("patch", 100.0);
    let i = interface();
    let gate = CancelGate::new();
    let mut state = NonsmoothTangentialState::new(1, 8).unwrap();
    let preload = state
        .prepare(
            &[contact(&l, &p, &i, "preload")],
            &[[6.0, 0.0, 0.0]],
            TangentialMobility {
                mechanical_dofs: 0,
                factor: &[],
            },
            0.001,
            controls(),
            &gate,
        )
        .unwrap();
    let stored = preload.contacts()[0].law.stored_energy_j;
    let force = preload.contacts()[0].action.force_n[0];
    state.commit(preload).unwrap();
    let old = state.clone();
    let mut parameters = l.parameters();
    parameters.tangential_stiffness_n_per_m *= 0.5;
    parameters.torsional_stiffness_nm_per_rad *= 0.5;
    let softened = PartialSlipLaw::new(l.model_id(), l.source_id(), parameters).unwrap();
    let trial = state
        .prepare(
            &[contact(&softened, &p, &i, "soften")],
            &[[0.0; 3]],
            TangentialMobility {
                mechanical_dofs: 0,
                factor: &[],
            },
            0.001,
            controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(state, old);
    close(trial.contacts()[0].action.force_n[0], 0.5 * force);
    close(trial.contacts()[0].law.stored_energy_j, 0.5 * stored);
    close(trial.audit().coefficient_work_j, -0.5 * stored);
    close(trial.audit().physical_heat_j, 0.0);
    close(trial.audit().energy_residual_j, 0.0);
    assert_eq!(
        trial.contacts()[0]
            .law
            .next_state
            .region_elastic_displacement_m(),
        old.contact_checkpoint("patch")
            .unwrap()
            .state()
            .region_elastic_displacement_m()
    );
    state.commit(trial).unwrap();
    assert_eq!(state.contact_checkpoint("patch").unwrap().law(), &softened);
}
