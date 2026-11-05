//! G1/G3/G4: shared occupied-branch friction with actual work and native consumers.
//! All constitutive parameters here are synthetic numerical fixtures.
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::stribeck_friction::{
    ContactBranch, HystereticFriction, ImplicitFrictionError, StribeckLoadLine,
};
use fs_tribo::tangential::matrix::{MatrixTangentialConfig, advance_with_matrix_compliance};
use fs_tribo::tangential::transport::{FramedTangentialState, TangentialTransportConfig};
use fs_tribo::{
    ApplicabilityRange, DryFrictionApplicability, DryInterfaceSystemCard, FrictionLaw,
    FrictionRegime, HeatPartition, InputAuthority, InterfaceMedium, InterfaceSystemRef,
    TangentialLaw, TangentialSolveConfig, TangentialStepInput, TriboError,
};

fn gaussian() -> FrictionLaw {
    FrictionLaw::Stribeck {
        static_mu: 0.8,
        kinetic_mu: 0.4,
        characteristic_speed: 0.04,
        viscous_per_speed: 0.0,
    }
}
fn input(free_speed: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: 3.9,
        slip_increment_m: [0.01 * free_speed, 0.0],
        dt_s: 0.01,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: Some(2e5),
        compliance_m_per_n: 0.002,
        solve: TangentialSolveConfig {
            max_iterations: 128,
            maximum_force_n: 10.0,
            force_absolute_tolerance_n: 1e-12,
            force_relative_tolerance: 1e-12,
        },
    }
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-12 * a.abs().max(b.abs()).max(1.0),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn gaussian_first_sliding_step_fits_ordinary_budget_and_reuses_preparation() {
    let law = HystereticFriction::new(gaussian()).unwrap();
    let fresh = law.initial_state();
    let first = law.advance(&fresh, &input(1.0)).unwrap();
    assert!(matches!(
        first.state.branch(),
        ContactBranch::SlidingHigh { positive: true }
    ));
    assert!(first.iterations <= 128);
    assert_eq!(fresh, law.initial_state());
    let repeated = law.advance(&first.state, &input(1.0)).unwrap();
    assert!(
        repeated.iterations < first.iterations,
        "{} vs {}",
        repeated.iterations,
        first.iterations
    );
    assert_eq!(repeated.traction_n, first.traction_n);
    assert_eq!(repeated.slip_increment_m, first.slip_increment_m);
    let continued = law.advance(&first.state, &input(0.5)).unwrap();
    let recaptured = law.advance(&fresh, &input(0.5)).unwrap();
    assert_eq!(recaptured.regime, FrictionRegime::Sticking);
    assert_eq!(continued.regime, FrictionRegime::Sliding);
    assert_eq!(continued.state.intersections(), (3, 2));
    assert!(continued.slip_increment_m[0] > 0.0);
    close(
        continued.dissipated_j,
        continued.traction_n[0] * continued.slip_increment_m[0],
    );
    close(continued.stored_change_j, 0.0);
    close(continued.closure_residual_j, 0.0);
    let partitioned = continued
        .dissipation_step(HeatPartition::new(0.25, 0.75, 0.0).unwrap())
        .unwrap();
    close(
        partitioned.surface_a_heat_j() + partitioned.surface_b_heat_j(),
        continued.dissipated_j,
    );
}

#[test]
fn hyperbolic_continuation_matches_independent_coexisting_root_and_restart() {
    let physical = FrictionLaw::Hyperbolic {
        static_mu: 0.8,
        kinetic_mu: 0.4,
        characteristic_speed: 1.0,
    };
    let law = HystereticFriction::new(physical.clone()).unwrap();
    let mut x = input(20.0);
    x.normal_n = 1.0;
    x.compliance_m_per_n = 22.5 * x.dt_s;
    let started = law.advance(&law.initial_state(), &x).unwrap();
    let checkpoint = started.state.clone();
    x.slip_increment_m[0] = 17.0 * x.dt_s;
    let continued = law.advance(&checkpoint, &x).unwrap();
    // v + 9 + 9/(1+v) = 17 gives v^2 - 7v + 1 = 0.
    // det-ok: independent quadratic oracle.
    let exact_speed = (7.0 + 45.0_f64.sqrt()) / 2.0;
    close(continued.slip_increment_m[0] / x.dt_s, exact_speed);
    close(continued.traction_n[0], (17.0 - exact_speed) / 22.5);
    assert_eq!(continued.state.intersections(), (3, 2));
    assert_eq!(law.advance(&checkpoint, &x).unwrap(), continued);
    assert_eq!(checkpoint, started.state);
    let original = StribeckLoadLine::new(&physical, 1.0, 22.5)
        .unwrap()
        .step(17.0, ContactBranch::SlidingHigh { positive: true })
        .unwrap();
    close(continued.traction_n[0], original.traction_n);
    close(
        continued.slip_increment_m[0] / x.dt_s,
        original.relative_velocity_m_s,
    );
    assert_eq!(
        law.advance(&law.initial_state(), &x).unwrap().regime,
        FrictionRegime::Sticking
    );
}

#[test]
fn branch_selection_precedes_effort_budget_and_never_switches_to_fit_it() {
    let law = HystereticFriction::new(FrictionLaw::Stribeck {
        static_mu: 1.0,
        kinetic_mu: 0.1,
        characteristic_speed: 1.0,
        viscous_per_speed: 0.0,
    })
    .unwrap();
    let mut x = input(0.0);
    x.normal_n = 1.0;
    x.dt_s = 1.0;
    x.compliance_m_per_n = 3.0;
    x.slip_increment_m = [4.0, 0.0];
    x.solve.maximum_force_n = 0.95;
    let moving = law.advance(&law.initial_state(), &x).unwrap();
    x.slip_increment_m[0] = 3.05;
    let checkpoint = moving.state.clone();
    let retained = law.advance(&checkpoint, &x).unwrap();
    assert!(retained.traction_n[0] < 0.2);
    assert!(matches!(
        retained.state.branch(),
        ContactBranch::SlidingHigh { .. }
    ));
    assert!(matches!(
        law.advance(&law.initial_state(), &x),
        Err(TriboError::OutsideApplicability {
            field: "maximum_force_n",
            ..
        })
    ));
    x.solve.maximum_force_n = 0.05;
    assert!(law.advance(&checkpoint, &x).is_err());
    assert_eq!(checkpoint, moving.state);
}

#[test]
fn coulomb_history_reverses_releases_and_recontacts_without_phantom_work() {
    let law = HystereticFriction::new(FrictionLaw::Coulomb {
        static_mu: 0.8,
        kinetic_mu: 0.3,
    })
    .unwrap();
    let mut x = input(0.0);
    x.normal_n = 1.0;
    x.dt_s = 1.0;
    x.compliance_m_per_n = 0.1;
    x.slip_increment_m = [0.2, 0.0];
    x.solve.max_iterations = 1;
    let moving = law.advance(&law.initial_state(), &x).unwrap();
    assert_eq!(moving.iterations, 0);
    x.slip_increment_m[0] = 0.05;
    let continued = law.advance(&moving.state, &x).unwrap();
    close(continued.traction_n[0], 0.3);
    close(continued.slip_increment_m[0], 0.02);
    assert_eq!(
        law.advance(&law.initial_state(), &x).unwrap().regime,
        FrictionRegime::Sticking
    );
    x.slip_increment_m[0] = -0.05;
    let reversed = law.advance(&continued.state, &x).unwrap();
    assert_eq!(reversed.state.branch(), ContactBranch::Stuck);
    assert_eq!(reversed.slip_increment_m, [0.0; 2]);
    assert_eq!(reversed.dissipated_j, 0.0);
    x.normal_n = 0.0;
    x.slip_increment_m = [0.03, 0.04];
    let released = law.advance(&moving.state, &x).unwrap();
    assert_eq!(released.state.branch(), ContactBranch::Fresh);
    assert_eq!(released.state.tangent_axis(), None);
    assert_eq!(released.slip_increment_m, x.slip_increment_m);
    assert_eq!(released.traction_n, [0.0; 2]);
    assert_eq!(released.dissipated_j, 0.0);
    x.normal_n = 1.0;
    x.slip_increment_m = [0.12, 0.16];
    let new_direction = law.advance(&released.state, &x).unwrap();
    close(new_direction.traction_n[0], 0.18);
    close(new_direction.traction_n[1], 0.24);
}

#[test]
fn preparation_and_root_caps_are_enforced_before_excess_work() {
    let physical = gaussian();
    let (line, preparation) = StribeckLoadLine::new_within(&physical, 3.9, 0.2, 128).unwrap();
    assert!(preparation > 0 && preparation < 75, "{preparation}");
    assert_eq!(
        StribeckLoadLine::new_within(&physical, 3.9, 0.2, preparation)
            .unwrap()
            .0,
        line
    );
    assert!(matches!(
        StribeckLoadLine::new_within(&physical, 3.9, 0.2, preparation - 1),
        Err(ImplicitFrictionError::NoConvergence)
    ));
    let (resolved, roots) = line.step_within(1.0, ContactBranch::Fresh, 128).unwrap();
    assert!(roots > 0);
    assert_eq!(
        line.step_within(1.0, ContactBranch::Fresh, roots)
            .unwrap()
            .0,
        resolved
    );
    assert!(matches!(
        line.step_within(1.0, ContactBranch::Fresh, roots - 1),
        Err(ImplicitFrictionError::NoConvergence)
    ));
    let law = HystereticFriction::new(physical).unwrap();
    let fresh = law.initial_state();
    let first = law.advance(&fresh, &input(1.0)).unwrap();
    assert_eq!(first.iterations, preparation + roots);
    let mut capped = input(1.0);
    capped.solve.max_iterations = first.iterations - 1;
    assert!(matches!(
        law.advance(&fresh, &capped),
        Err(TriboError::TangentialSolve { .. })
    ));
    let (prepared, spent) = law.prepare_state(&fresh, 3.9, 0.2, preparation).unwrap();
    assert_eq!(spent, preparation);
    assert_eq!(prepared.branch(), fresh.branch());
    assert_eq!(prepared.slip_m(), fresh.slip_m());
    assert_eq!(
        prepared.contact_temperature_k(),
        fresh.contact_temperature_k()
    );
    let after = law.advance(&prepared, &input(1.0)).unwrap();
    assert_eq!(after.iterations, roots);
    assert_eq!(after.traction_n, first.traction_n);
    assert_eq!(after.state, first.state);
    let (same, spent) = law.prepare_state(&prepared, 3.9, 0.2, 0).unwrap();
    assert_eq!(same, prepared);
    assert_eq!(spent, 0);
}

#[test]
fn changed_load_rebuilds_the_original_line_and_preserves_branch_history() {
    let law = HystereticFriction::new(gaussian()).unwrap();
    let first = law.advance(&law.initial_state(), &input(1.0)).unwrap();
    let mut x = input(0.5);
    x.normal_n = 3.8;
    let changed = law.advance(&first.state, &x).unwrap();
    let original = StribeckLoadLine::new(law.friction_law(), 3.8, 0.2)
        .unwrap()
        .step(0.5, first.state.branch())
        .unwrap();
    close(changed.traction_n[0], original.traction_n);
    assert!(
        changed
            .state
            .prepared_load_line()
            .unwrap()
            .matches(3.8, 0.2)
    );
    assert!(
        !changed
            .state
            .prepared_load_line()
            .unwrap()
            .matches(3.9, 0.2)
    );
    assert!(matches!(
        changed.state.branch(),
        ContactBranch::SlidingHigh { .. }
    ));
}

#[test]
fn native_matrix_consumer_keeps_the_occupied_branch_and_actual_work() {
    let law = HystereticFriction::new(gaussian()).unwrap();
    let mut moving = input(1.0);
    moving.slip_increment_m = [0.006, 0.008];
    let state = law.advance(&law.initial_state(), &moving).unwrap().state;
    let mut x = input(0.5);
    x.slip_increment_m = [0.003, 0.004];
    let expected = law.advance(&state, &x).unwrap();
    let compliance = x.compliance_m_per_n;
    x.compliance_m_per_n = 0.0;
    let actual = advance_with_matrix_compliance(
        &law,
        &state,
        &x,
        [[compliance, 0.0], [0.0, compliance]],
        MatrixTangentialConfig {
            maximum_evaluations: 8,
            displacement_tolerance_m: 1e-13,
            energy_tolerance_j: 1e-13,
            matrix_relative_tolerance: 1e-13,
        },
        || false,
    )
    .unwrap();
    assert_eq!(actual.contact, expected);
    assert_eq!(actual.contact.regime, FrictionRegime::Sliding);
    close(actual.body_work_residual_j, 0.0);
}

#[test]
fn objective_transport_preserves_sliding_branch_and_oblique_motion_refuses() {
    let law = HystereticFriction::new(gaussian()).unwrap();
    let moving = law.advance(&law.initial_state(), &input(1.0)).unwrap();
    let base = TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
    let accepted = FramedTangentialState::new(&law, moving.state.clone(), base).unwrap();
    let current = base.rotated(0.73).unwrap();
    let controls = TangentialTransportConfig {
        absolute_energy_tolerance_j: 0.0,
        relative_energy_tolerance: 0.0,
    };
    let identity = accepted
        .transport(&law, base, 0.0, controls, || false)
        .unwrap();
    assert_eq!(identity.checkpoint, accepted);
    let rotated = accepted
        .transport(&law, current, 0.0, controls, || false)
        .unwrap();
    assert_eq!(rotated.checkpoint.state().branch(), moving.state.branch());
    assert_eq!(rotated.stored_energy_residual_j, 0.0);
    let axis = rotated.checkpoint.state().tangent_axis().unwrap();
    let mut x = input(0.5);
    x.slip_increment_m = [0.005 * axis[0], 0.005 * axis[1]];
    let transported = law.advance(rotated.checkpoint.state(), &x).unwrap();
    let reference = law.advance(&moving.state, &input(0.5)).unwrap();
    for i in 0..3 {
        let actual = current.longitudinal()[i] * transported.traction_n[0]
            + current.lateral()[i] * transported.traction_n[1];
        close(actual, base.longitudinal()[i] * reference.traction_n[0]);
    }
    close(transported.dissipated_j, reference.dissipated_j);
    x.slip_increment_m = [0.005, 0.001];
    assert!(matches!(
        law.advance(&moving.state, &x),
        Err(TriboError::InvalidInput {
            field: "oblique_hysteretic_continuation"
        })
    ));
    assert_eq!(accepted.state(), &moving.state);
    assert!(
        accepted
            .transport(&law, current, 0.0, controls, || true)
            .is_err()
    );
}

fn card(identity: &str) -> DryInterfaceSystemCard {
    DryInterfaceSystemCard::new(
        InterfaceSystemRef::new(
            identity,
            "accepted-history",
            "synthetic-interface",
            InputAuthority::SyntheticFixture,
            InterfaceMedium::Dry,
        )
        .unwrap(),
        FrictionLaw::Coulomb {
            static_mu: 0.8,
            kinetic_mu: 0.3,
        },
        DryFrictionApplicability::new(
            ApplicabilityRange::new(280.0, 320.0).unwrap(),
            ApplicabilityRange::new(0.0, 5e5).unwrap(),
            ApplicabilityRange::new(0.0, 0.15).unwrap(),
        )
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn original_card_admits_actual_motion_and_refuses_foreign_history_or_domain() {
    let original = card("surface-a->b");
    let law = HystereticFriction::from_card(original.clone());
    assert_eq!(law.interface_card(), Some(&original));
    let mut x = input(0.0);
    x.normal_n = 1.0;
    x.dt_s = 1.0;
    x.compliance_m_per_n = 0.1;
    x.slip_increment_m = [0.17, 0.0];
    // The released speed exceeds the card, but the actual physical slip does not.
    let accepted = law.advance(&law.initial_state(), &x).unwrap();
    close(accepted.slip_increment_m[0], 0.14);
    x.slip_increment_m[0] = 0.2;
    assert!(matches!(
        law.advance(&accepted.state, &x),
        Err(TriboError::OutsideApplicability { .. })
    ));
    x.slip_increment_m[0] = 0.17;
    x.contact_temperature_k = Some(321.0);
    assert!(law.advance(&accepted.state, &x).is_err());
    x.contact_temperature_k = None;
    assert!(law.advance(&accepted.state, &x).is_err());
    let foreign = HystereticFriction::from_card(card("another-interface"));
    assert!(foreign.stored_energy_j(&accepted.state).is_nan());
    assert!(foreign.advance(&accepted.state, &input(0.5)).is_err());
}

#[test]
fn prescribed_motion_and_unrepresentable_history_keep_transactional_limits() {
    let law = HystereticFriction::new(FrictionLaw::Coulomb {
        static_mu: 0.8,
        kinetic_mu: 0.3,
    })
    .unwrap();
    let mut x = input(0.0);
    x.normal_n = 1.0;
    x.compliance_m_per_n = 0.0;
    x.slip_increment_m = [0.03, 0.04];
    let prescribed = law.advance(&law.initial_state(), &x).unwrap();
    assert_eq!(prescribed.slip_increment_m, x.slip_increment_m);
    close(prescribed.traction_n[0], 0.18);
    close(prescribed.traction_n[1], 0.24);
    close(prescribed.dissipated_j, 0.015);
    x.normal_n = 0.0;
    x.slip_increment_m = [f64::MAX, 0.0];
    let huge = law.advance(&law.initial_state(), &x).unwrap();
    let saved = huge.state.clone();
    assert!(law.advance(&huge.state, &x).is_err());
    assert_eq!(huge.state, saved);
    x.normal_n = 1.0;
    x.slip_increment_m = [f64::MAX, f64::MAX];
    assert!(law.advance(&law.initial_state(), &x).is_err());
}
