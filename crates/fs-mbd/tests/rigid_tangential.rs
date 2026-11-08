//! Original rigid-body integration of retained full-matrix friction laws.

use fs_exec::CancelGate;
use fs_mbd::tangential_friction::directional::{
    RigidDirectionalContact, RigidDirectionalControls, RigidDirectionalError,
    apply_rigid_directional,
};
use fs_mbd::tangential_friction::{
    RigidContactBinding, RigidTangentialControls, RigidTangentialPhase, apply_rigid_tangential,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::layer_contact::ThirdBodyLayerContact;
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::stribeck_friction::{ContactBranch, HystereticFriction};
use fs_tribo::tangential::elliptic::EllipticCoulomb;
use fs_tribo::tangential::matrix::{MatrixTangentialConfig, MatrixTangentialError};
use fs_tribo::{
    DryTangentialState, FrictionLaw, TangentialLaw, TangentialSolveConfig, TangentialStep,
    TangentialStepInput, TriboError,
};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-9 * actual.abs().max(expected.abs()).max(1.0),
        "{actual:e} != {expected:e}"
    );
}

fn controls(phase: RigidTangentialPhase) -> RigidTangentialControls {
    RigidTangentialControls {
        mechanical: RigidDirectionalControls {
            maximum_bodies: 8,
            geometry_tolerance_m: 1e-11,
            velocity_tolerance_m_s: 1e-9,
            linear_momentum_tolerance_ns: 1e-10,
            angular_momentum_tolerance_nms: 1e-10,
            energy_tolerance_j: 1e-9,
            contact: TangentialSolveConfig {
                max_iterations: 128,
                maximum_force_n: 1e4,
                force_absolute_tolerance_n: 1e-11,
                force_relative_tolerance: 1e-13,
            },
        },
        coupling: MatrixTangentialConfig {
            maximum_evaluations: 512,
            displacement_tolerance_m: 1e-12,
            energy_tolerance_j: 1e-10,
            matrix_relative_tolerance: 1e-13,
        },
        phase,
    }
}

fn contact(arms: [Vec3; 2]) -> RigidDirectionalContact {
    RigidDirectionalContact {
        binding: RigidContactBinding {
            body_a: 0,
            body_b: 1,
            arm_a_body: arms[0],
            arm_b_body: arms[1],
        },
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
        normal_load_n: 10.0,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: Some(1e4),
    }
}

fn translation_pair(speed: f64) -> ([RigidBodyState; 2], [MassProperties; 2]) {
    let properties = [MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap(); 2];
    let bodies = [
        RigidBodyState::new(Pose::identity(), Vec3::new(speed, 0.0, 0.0), Vec3::ZERO).unwrap(),
        RigidBodyState::new(Pose::identity(), Vec3::ZERO, Vec3::ZERO).unwrap(),
    ];
    (bodies, properties)
}

fn preload(layer: &ThirdBodyLayerContact) -> <ThirdBodyLayerContact as TangentialLaw>::State {
    layer
        .advance(
            &layer.initial_state(),
            &TangentialStepInput {
                normal_n: 10.0,
                slip_increment_m: [0.01, 0.0],
                dt_s: 0.01,
                contact_temperature_k: Some(300.0),
                nominal_pressure_pa: Some(1e4),
                compliance_m_per_n: 0.0,
                solve: controls(RigidTangentialPhase::StoredEnergy)
                    .mechanical
                    .contact,
            },
        )
        .unwrap()
        .state
}

#[test]
fn coupled_full_inertia_matches_independent_native_coulomb_projection() {
    let arms = [Vec3::new(0.4, 0.2, -0.3), Vec3::new(-0.5, 0.6, 0.4)];
    let properties = [
        MassProperties::from_principal_frame(
            2.0,
            Vec3::ZERO,
            Vec3::new(1.0, 1.5, 2.0),
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.5).unwrap(),
        )
        .unwrap(),
        MassProperties::from_principal_frame(
            3.0,
            Vec3::ZERO,
            Vec3::new(0.8, 1.2, 1.6),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), -0.4).unwrap(),
        )
        .unwrap(),
        MassProperties::new(5.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap(),
    ];
    let initial = [
        RigidBodyState::new(
            Pose::new(arms[0].scale(-1.0), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(2.0, 0.8, 0.0),
            Vec3::ZERO,
        )
        .unwrap(),
        RigidBodyState::new(
            Pose::new(arms[1].scale(-1.0), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(-0.4, 0.1, 0.0),
            Vec3::ZERO,
        )
        .unwrap(),
        RigidBodyState::new(Pose::identity(), Vec3::new(3.0, 2.0, 1.0), Vec3::ZERO).unwrap(),
    ];
    let law = FrictionLaw::Coulomb {
        static_mu: 0.7,
        kinetic_mu: 0.3,
    };
    let reference_law = EllipticCoulomb::new([0.7; 2], [0.3; 2], [1.0, 0.0]).unwrap();
    let mut history = law.initial_state();
    let mut reference_history = reference_law.initial_state();
    let mut bodies = initial;
    let mut reference_bodies = initial;
    let cfg = controls(RigidTangentialPhase::Dissipative);
    let result = apply_rigid_tangential(
        &mut bodies,
        &properties,
        contact(arms),
        &law,
        &mut history,
        0.002,
        cfg,
        &CancelGate::new(),
    )
    .unwrap();
    let reference = apply_rigid_directional(
        &mut reference_bodies,
        &properties,
        contact(arms),
        &reference_law,
        &mut reference_history,
        0.002,
        cfg.mechanical,
        &CancelGate::new(),
    )
    .unwrap();
    assert!(result.mobility_per_kg[0][1].abs() > 1e-3);
    close(
        result.contact.traction_n[0].hypot(result.contact.traction_n[1]),
        3.0,
    );
    for i in 0..2 {
        close(
            result.contact.traction_n[i],
            reference.contact.traction_n[i],
        );
        close(
            result.slip_velocity_after_m_s[i],
            reference.slip_velocity_after_m_s[i],
        );
        close(history.slip_m()[i], result.contact.slip_increment_m[i]);
        for j in 0..2 {
            close(
                result.compliance_m_per_n[i][j],
                0.5 * 0.002_f64.powi(2) * result.mobility_per_kg[i][j],
            );
        }
        let delta = bodies[i]
            .linear_momentum_world()
            .sub(reference_bodies[i].linear_momentum_world());
        close(delta.norm_squared().sqrt(), 0.0);
        let spin_delta = bodies[i]
            .angular_momentum_body()
            .sub(reference_bodies[i].angular_momentum_body());
        close(spin_delta.norm_squared().sqrt(), 0.0);
        assert_eq!(bodies[i].pose(), initial[i].pose());
    }
    assert_eq!(bodies[2], initial[2]);
    assert_eq!(history, result.contact.state);
    assert!(result.contact.dissipated_j > 0.0);
    close(
        result.kinetic_energy_change_j + result.contact.dissipated_j,
        0.0,
    );
    close(
        result.audit.linear_momentum_change_ns.norm_squared().sqrt(),
        0.0,
    );
    close(
        result
            .audit
            .angular_momentum_change_nms
            .norm_squared()
            .sqrt(),
        0.0,
    );
}

#[test]
fn scalar_hysteretic_history_and_physical_preparation_survive_native_steps() {
    let law = HystereticFriction::new(FrictionLaw::Stribeck {
        static_mu: 0.6,
        kinetic_mu: 0.2,
        characteristic_speed: 0.2,
        viscous_per_speed: 0.01,
    })
    .unwrap();
    let (mut bodies, properties) = translation_pair(1.0);
    let mut history = law.initial_state();
    let dt = 0.005;
    for _ in 0..3 {
        let previous_slip = history.slip_m();
        let report = apply_rigid_tangential(
            &mut bodies,
            &properties,
            contact([Vec3::ZERO; 2]),
            &law,
            &mut history,
            dt,
            controls(RigidTangentialPhase::Dissipative),
            &CancelGate::new(),
        )
        .unwrap();
        assert!(matches!(
            history.branch(),
            ContactBranch::SlidingLow { positive: true }
                | ContactBranch::SlidingHigh { positive: true }
        ));
        assert_eq!(history.tangent_axis(), Some([1.0, 0.0]));
        assert!(
            history
                .prepared_load_line()
                .unwrap()
                .matches(10.0, report.compliance_m_per_n[0][0] / dt)
        );
        close(
            history.slip_m()[0],
            previous_slip[0] + report.contact.slip_increment_m[0],
        );
        assert_eq!(history, report.contact.state);
        assert!(report.contact.dissipated_j > 0.0);
        close(
            report.kinetic_energy_change_j + report.contact.dissipated_j,
            0.0,
        );
    }
}

#[test]
fn preloaded_layer_returns_stored_energy_to_both_original_bodies() {
    let law = ThirdBodyLayerContact::new(100.0, 1000.0, 1.0, 300.0, 0.1).unwrap();
    let initial_history = preload(&law);
    let stored_before = law.stored_energy_j(&initial_history);
    close(stored_before, 0.005);
    let (initial, properties) = translation_pair(0.0);

    let mut refused_bodies = initial;
    let mut refused_history = initial_history.clone();
    assert!(matches!(
        apply_rigid_tangential(
            &mut refused_bodies,
            &properties,
            contact([Vec3::ZERO; 2]),
            &law,
            &mut refused_history,
            0.1,
            controls(RigidTangentialPhase::Dissipative),
            &CancelGate::new(),
        ),
        Err(RigidDirectionalError::UnresolvedSlipEvent { .. })
    ));
    assert_eq!(refused_bodies, initial);
    assert_eq!(refused_history, initial_history);

    let mut bodies = initial;
    let mut history = initial_history;
    let report = apply_rigid_tangential(
        &mut bodies,
        &properties,
        contact([Vec3::ZERO; 2]),
        &law,
        &mut history,
        0.1,
        controls(RigidTangentialPhase::StoredEnergy),
        &CancelGate::new(),
    )
    .unwrap();
    // C = dt² I = .01 I, F = k e_old / (1 + k C / 2) = 2/3 N.
    // The original spring pays the complete 1/225 J native kinetic increase.
    close(report.contact.traction_n[0], 2.0 / 3.0);
    close(report.slip_velocity_after_m_s[0], -2.0 / 15.0);
    close(report.kinetic_energy_change_j, 1.0 / 225.0);
    close(report.contact.dissipated_j, 0.0);
    close(law.stored_energy_j(&history) - stored_before, -1.0 / 225.0);
    close(report.audit.total_energy_residual_j, 0.0);
    assert!(report.audit.maximum_assisting_work_j > 0.0);
    close(bodies[0].linear_momentum_world().x, -1.0 / 15.0);
    close(bodies[1].linear_momentum_world().x, 1.0 / 15.0);
    assert_eq!(history, report.contact.state);
}

#[test]
fn stored_energy_selection_cannot_bypass_dry_arrest_or_effort_limits() {
    let law = FrictionLaw::Coulomb {
        static_mu: 0.4,
        kinetic_mu: 0.2,
    };
    let (initial, properties) = translation_pair(1.0);
    for phase in [
        RigidTangentialPhase::Dissipative,
        RigidTangentialPhase::StoredEnergy,
    ] {
        let mut bodies = initial;
        let mut history = law.initial_state();
        let original_history = history.clone();
        assert!(matches!(
            apply_rigid_tangential(
                &mut bodies,
                &properties,
                contact([Vec3::ZERO; 2]),
                &law,
                &mut history,
                0.4,
                controls(phase),
                &CancelGate::new(),
            ),
            Err(RigidDirectionalError::UnresolvedSlipEvent { .. })
        ));
        assert_eq!(bodies, initial);
        assert_eq!(history, original_history);
        let mut bounded = controls(phase);
        bounded.mechanical.contact.maximum_force_n = 0.1;
        assert!(
            apply_rigid_tangential(
                &mut bodies,
                &properties,
                contact([Vec3::ZERO; 2]),
                &law,
                &mut history,
                0.01,
                bounded,
                &CancelGate::new(),
            )
            .is_err()
        );
        assert_eq!(bodies, initial);
        assert_eq!(history, original_history);
    }
}

struct CancelAfterNativeLaw<'a> {
    law: FrictionLaw,
    gate: &'a CancelGate,
}
impl TangentialLaw for CancelAfterNativeLaw<'_> {
    type State = DryTangentialState;
    fn model_id(&self) -> &str {
        self.law.model_id()
    }
    fn initial_state(&self) -> Self::State {
        self.law.initial_state()
    }
    fn stored_energy_j(&self, state: &Self::State) -> f64 {
        self.law.stored_energy_j(state)
    }
    fn advance(
        &self,
        state: &Self::State,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        let candidate = self.law.advance(state, input);
        self.gate.request();
        candidate
    }
}

#[test]
fn cancellation_after_native_law_trial_preserves_both_bodies_and_history() {
    let gate = CancelGate::new();
    let law = CancelAfterNativeLaw {
        law: FrictionLaw::Coulomb {
            static_mu: 0.4,
            kinetic_mu: 0.2,
        },
        gate: &gate,
    };
    let (mut bodies, properties) = translation_pair(1.0);
    let initial = bodies;
    let mut history = law.initial_state();
    let original_history = history.clone();
    let result = apply_rigid_tangential(
        &mut bodies,
        &properties,
        contact([Vec3::ZERO; 2]),
        &law,
        &mut history,
        0.01,
        controls(RigidTangentialPhase::Dissipative),
        &gate,
    );
    assert!(gate.is_requested());
    assert!(matches!(
        result,
        Err(RigidDirectionalError::Cancelled)
            | Err(RigidDirectionalError::Source(
                MatrixTangentialError::Cancelled
            ))
    ));
    assert_eq!(bodies, initial);
    assert_eq!(history, original_history);
}

#[test]
fn foreign_stored_history_is_refused_before_native_publication() {
    let original_law = ThirdBodyLayerContact::new(100.0, 1000.0, 1.0, 300.0, 0.1).unwrap();
    let foreign_law = ThirdBodyLayerContact::new(200.0, 1000.0, 1.0, 300.0, 0.1).unwrap();
    let mut history = preload(&original_law);
    let original_history = history.clone();
    let (mut bodies, properties) = translation_pair(0.0);
    let initial = bodies;
    assert!(matches!(
        apply_rigid_tangential(
            &mut bodies,
            &properties,
            contact([Vec3::ZERO; 2]),
            &foreign_law,
            &mut history,
            0.1,
            controls(RigidTangentialPhase::StoredEnergy),
            &CancelGate::new(),
        ),
        Err(RigidDirectionalError::Invalid(_))
    ));
    assert_eq!(bodies, initial);
    assert_eq!(history, original_history);
}
