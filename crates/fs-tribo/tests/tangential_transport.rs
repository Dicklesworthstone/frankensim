//! Objective native-history transport; all constitutive data are synthetic.
use fs_tribo::layer_contact::{ThirdBodyLayerContact, ThirdBodyLayerState};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::matrix::{MatrixTangentialConfig, advance_with_matrix_compliance};
use fs_tribo::tangential::transport::{
    FramedTangentialState, ObjectiveTangentialLaw, TangentialHistoryRotation,
    TangentialTransportConfig, TangentialTransportError,
};
use fs_tribo::{
    ApplicabilityRange, DryFrictionApplicability, DryInterfaceSystemCard, FrictionLaw,
    InputAuthority, InterfaceMedium, InterfaceSystemRef, TangentialLaw, TangentialSolveConfig,
    TangentialStep, TangentialStepInput, TriboError,
};

fn config() -> TangentialTransportConfig {
    TangentialTransportConfig {
        absolute_energy_tolerance_j: 1e-15,
        relative_energy_tolerance: 1e-13,
    }
}
fn base() -> TangentFrame {
    TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap()
}
fn input(slip_increment_m: [f64; 2]) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: 3.0,
        slip_increment_m,
        dt_s: 0.02,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: Some(2e5),
        compliance_m_per_n: 0.0,
        solve: TangentialSolveConfig {
            max_iterations: 128,
            maximum_force_n: 1000.0,
            force_absolute_tolerance_n: 1e-13,
            force_relative_tolerance: 1e-13,
        },
    }
}
fn loaded() -> (ThirdBodyLayerContact, ThirdBodyLayerState) {
    let law = ThirdBodyLayerContact::new(450.0, 0.7, 0.8, 300.0, 0.2).unwrap();
    let mut state = law.initial_state();
    for slip in [[0.018, -0.007], [-0.004, 0.016], [0.011, 0.003]] {
        state = law.advance(&state, &input(slip)).unwrap().state;
    }
    assert!(state.elastic_slip_m().into_iter().all(|x| x != 0.0));
    assert!(state.plastic_slip_m().into_iter().all(|x| x != 0.0));
    (law, state)
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}, tolerance {tolerance:e}"
    );
}
fn close2(a: [f64; 2], b: [f64; 2], tolerance: f64) {
    for i in 0..2 {
        close(a[i], b[i], tolerance);
    }
}
fn close3(a: [f64; 3], b: [f64; 3]) {
    for i in 0..3 {
        close(a[i], b[i], 4e-14);
    }
}
fn rotate2(v: [f64; 2], angle: f64) -> [f64; 2] {
    let (s, c) = angle.sin_cos();
    [c * v[0] - s * v[1], s * v[0] + c * v[1]]
}
fn world(frame: TangentFrame, v: [f64; 2]) -> [f64; 3] {
    std::array::from_fn(|i| frame.longitudinal()[i] * v[0] + frame.lateral()[i] * v[1])
}
fn rotate3(v: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    // Independent Rodrigues reference, not the production frame map.
    let (s, c) = angle.sin_cos();
    let dot = (0..3).map(|i| axis[i] * v[i]).sum::<f64>();
    let cross = [
        axis[1] * v[2] - axis[2] * v[1],
        axis[2] * v[0] - axis[0] * v[2],
        axis[0] * v[1] - axis[1] * v[0],
    ];
    std::array::from_fn(|i| c * v[i] + s * cross[i] + (1.0 - c) * axis[i] * dot)
}

#[test]
fn loaded_gauge_rebase_preserves_world_shear_plastic_history_and_native_owner() {
    let (law, state) = loaded();
    let accepted = FramedTangentialState::new(&law, state.clone(), base()).unwrap();
    let saved = accepted.clone();
    let angle = 0.73;
    let current = base().rotated(angle).unwrap();
    let transported = accepted
        .transport(&law, current, 0.0, config(), || false)
        .unwrap();
    assert_eq!(accepted, saved);
    assert_eq!(transported.checkpoint.frame(), current);
    let next = transported.checkpoint.state();
    close2(
        next.elastic_slip_m(),
        rotate2(state.elastic_slip_m(), -angle),
        1e-16,
    );
    close2(
        next.plastic_slip_m(),
        rotate2(state.plastic_slip_m(), -angle),
        1e-16,
    );
    close3(
        world(current, next.elastic_slip_m()),
        world(base(), state.elastic_slip_m()),
    );
    close3(
        world(current, next.plastic_slip_m()),
        world(base(), state.plastic_slip_m()),
    );
    close(
        law.stored_energy_j(next),
        law.stored_energy_j(&state),
        1e-15,
    );
    close(transported.stored_energy_residual_j, 0.0, 1e-15);
    let identity = accepted
        .transport(&law, base(), 0.0, config(), || false)
        .unwrap();
    assert_eq!(identity.checkpoint, accepted);
    assert_eq!(identity.stored_energy_residual_j, 0.0);
    let foreign = ThirdBodyLayerContact::new(451.0, 0.7, 0.8, 300.0, 0.2).unwrap();
    assert!(FramedTangentialState::new(&foreign, next.clone(), current).is_err());
    assert!(foreign.advance(next, &input([0.0; 2])).is_err());
}

#[test]
fn minimum_normal_rotation_and_common_material_spin_are_not_projection_or_gauge() {
    let (law, state) = loaded();
    let original = FramedTangentialState::new(&law, state.clone(), base()).unwrap();
    let tilt = 1.1;
    let tilt_axis = [0.0, 1.0, 0.0];
    let normal = rotate3(base().normal(), tilt_axis, tilt);
    let longitudinal = rotate3(base().longitudinal(), tilt_axis, tilt);
    let tilted = TangentFrame::new(normal, longitudinal).unwrap();
    let spin = -0.62;
    let moved = original
        .transport(&law, tilted, spin, config(), || false)
        .unwrap();
    for (before, after) in [
        (
            state.elastic_slip_m(),
            moved.checkpoint.state().elastic_slip_m(),
        ),
        (
            state.plastic_slip_m(),
            moved.checkpoint.state().plastic_slip_m(),
        ),
    ] {
        close2(after, rotate2(before, spin), 2e-16);
        let expected = rotate3(
            rotate3(world(base(), before), tilt_axis, tilt),
            normal,
            spin,
        );
        close3(world(tilted, after), expected);
        close(after[0].hypot(after[1]), before[0].hypot(before[1]), 1e-16);
    }
    // When the material and supplied axes spin together, components are fixed.
    let corotated = tilted.rotated(spin).unwrap();
    let co = original
        .transport(&law, corotated, spin, config(), || false)
        .unwrap();
    close2(
        co.checkpoint.state().elastic_slip_m(),
        state.elastic_slip_m(),
        1e-16,
    );
    close2(
        co.checkpoint.state().plastic_slip_m(),
        state.plastic_slip_m(),
        1e-16,
    );
    close(
        law.stored_energy_j(co.checkpoint.state()),
        law.stored_energy_j(&state),
        1e-15,
    );
}

#[test]
fn arbitrary_global_rigid_rotation_preserves_the_complete_transport_result() {
    let (law, state) = loaded();
    let before = TangentFrame::new([0.2, -0.4, 0.9], [0.7, 0.4, 0.2]).unwrap();
    let after = TangentFrame::new([-0.5, 0.3, 0.8], [0.2, 0.9, -0.1]).unwrap();
    let axis = [2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0];
    let q = |f: TangentFrame| {
        TangentFrame::new(
            rotate3(f.normal(), axis, 1.23),
            rotate3(f.longitudinal(), axis, 1.23),
        )
        .unwrap()
    };
    let a = FramedTangentialState::new(&law, state.clone(), before)
        .unwrap()
        .transport(&law, after, 0.37, config(), || false)
        .unwrap();
    let b = FramedTangentialState::new(&law, state, q(before))
        .unwrap()
        .transport(&law, q(after), 0.37, config(), || false)
        .unwrap();
    close2(
        a.checkpoint.state().elastic_slip_m(),
        b.checkpoint.state().elastic_slip_m(),
        1e-16,
    );
    close2(
        a.checkpoint.state().plastic_slip_m(),
        b.checkpoint.state().plastic_slip_m(),
        1e-16,
    );
    close3(
        world(q(after), b.checkpoint.state().elastic_slip_m()),
        rotate3(
            world(after, a.checkpoint.state().elastic_slip_m()),
            axis,
            1.23,
        ),
    );
    close(
        a.stored_energy_residual_j,
        b.stored_energy_residual_j,
        1e-15,
    );
}

#[test]
fn inverse_gauge_composition_and_resolved_normal_path_preserve_history() {
    let (law, state) = loaded();
    let original = FramedTangentialState::new(&law, state.clone(), base()).unwrap();
    let frame1 = base().rotated(0.54).unwrap();
    let frame2 = base().rotated(-0.93).unwrap();
    let first = original
        .transport(&law, frame1, 0.31, config(), || false)
        .unwrap();
    let composed = first
        .checkpoint
        .transport(&law, frame2, -0.47, config(), || false)
        .unwrap();
    let direct = original
        .transport(&law, frame2, -0.16, config(), || false)
        .unwrap();
    close2(
        composed.checkpoint.state().elastic_slip_m(),
        direct.checkpoint.state().elastic_slip_m(),
        1e-16,
    );
    close2(
        composed.checkpoint.state().plastic_slip_m(),
        direct.checkpoint.state().plastic_slip_m(),
        1e-16,
    );
    let tilted = TangentFrame::new([0.5, -0.4, 0.7], [0.2, 0.6, 0.4]).unwrap();
    let out = original
        .transport(&law, tilted, 0.28, config(), || false)
        .unwrap();
    let back = out
        .checkpoint
        .transport(&law, base(), -0.28, config(), || false)
        .unwrap();
    close2(
        back.checkpoint.state().elastic_slip_m(),
        state.elastic_slip_m(),
        1e-16,
    );
    close2(
        back.checkpoint.state().plastic_slip_m(),
        state.plastic_slip_m(),
        1e-16,
    );
    // z -> x -> y -> z carries the shortest-rotation holonomy of pi/2.
    // Returning to the same normal does not justify forgetting intermediate frames.
    let mut path = original;
    for frame in [
        TangentFrame::new([1.0, 0.0, 0.0], [0.0, 0.0, -1.0]).unwrap(),
        TangentFrame::new([0.0, 1.0, 0.0], [1.0, 0.0, 0.0]).unwrap(),
        base(),
    ] {
        path = path
            .transport(&law, frame, 0.0, config(), || false)
            .unwrap()
            .checkpoint;
    }
    let e = state.elastic_slip_m();
    let p = state.plastic_slip_m();
    close2(path.state().elastic_slip_m(), [-e[1], e[0]], 1e-16);
    close2(path.state().plastic_slip_m(), [-p[1], p[0]], 1e-16);
    close(
        law.stored_energy_j(path.state()),
        law.stored_energy_j(&state),
        1e-15,
    );
}

#[test]
fn dry_ladder_and_source_card_keep_accumulated_slip_and_explicit_temperature() {
    let laws = [
        FrictionLaw::Coulomb {
            static_mu: 0.8,
            kinetic_mu: 0.3,
        },
        FrictionLaw::VelocityDependent {
            static_mu: 0.8,
            mu_zero: 0.4,
            slope_per_speed: 0.1,
        },
        FrictionLaw::Stribeck {
            static_mu: 0.8,
            kinetic_mu: 0.3,
            characteristic_speed: 0.05,
            viscous_per_speed: 0.01,
        },
        FrictionLaw::Hyperbolic {
            static_mu: 0.8,
            kinetic_mu: 0.3,
            characteristic_speed: 0.05,
        },
    ];
    let angle = -0.83;
    let frame = base().rotated(angle).unwrap();
    for law in &laws {
        let step = law
            .advance(&law.initial_state(), &input([0.003, -0.004]))
            .unwrap();
        let framed = FramedTangentialState::new(law, step.state.clone(), base()).unwrap();
        let result = framed
            .transport(law, frame, 0.0, config(), || false)
            .unwrap();
        close2(
            result.checkpoint.state().slip_m(),
            rotate2(step.state.slip_m(), -angle),
            1e-16,
        );
        assert_eq!(
            result.checkpoint.state().contact_temperature_k(),
            Some(300.0)
        );
        assert_eq!(result.stored_energy_residual_j, 0.0);
        let du = [0.002, 0.001];
        let expected = law.advance(&step.state, &input(du)).unwrap();
        let actual = law
            .advance(result.checkpoint.state(), &input(rotate2(du, -angle)))
            .unwrap();
        close2(
            actual.traction_n,
            rotate2(expected.traction_n, -angle),
            1e-13,
        );
        close2(
            actual.state.slip_m(),
            rotate2(expected.state.slip_m(), -angle),
            1e-16,
        );
        close(actual.dissipated_j, expected.dissipated_j, 1e-15);
        let unloaded = FramedTangentialState::new(law, law.initial_state(), base()).unwrap();
        assert_eq!(
            unloaded
                .transport(law, frame, 0.2, config(), || false)
                .unwrap()
                .checkpoint
                .state()
                .contact_temperature_k(),
            None
        );
    }
    let card = DryInterfaceSystemCard::new(
        InterfaceSystemRef::new(
            "a->b",
            "retained-dry",
            "transport-fixture",
            InputAuthority::SyntheticFixture,
            InterfaceMedium::Dry,
        )
        .unwrap(),
        laws[0].clone(),
        DryFrictionApplicability::new(
            ApplicabilityRange::new(280.0, 320.0).unwrap(),
            ApplicabilityRange::new(0.0, 5e5).unwrap(),
            ApplicabilityRange::new(0.0, 2.0).unwrap(),
        )
        .unwrap(),
    )
    .unwrap();
    let saved_card = card.clone();
    let state = card
        .advance(&card.initial_state(), &input([0.003, -0.004]))
        .unwrap()
        .state;
    let framed = FramedTangentialState::new(&card, state.clone(), base()).unwrap();
    let result = framed
        .transport(&card, frame, 0.0, config(), || false)
        .unwrap();
    close2(
        result.checkpoint.state().slip_m(),
        rotate2(state.slip_m(), -angle),
        1e-16,
    );
    assert_eq!(card, saved_card);
    let invalid_step = TangentialStepInput {
        contact_temperature_k: Some(400.0),
        ..input([0.001, 0.0])
    };
    assert!(
        card.advance(result.checkpoint.state(), &invalid_step)
            .is_err()
    );
}

#[test]
fn loaded_layer_continues_with_rotated_full_body_response_and_identical_physical_heat() {
    let (law, state) = loaded();
    let angle = 0.71;
    let framed = FramedTangentialState::new(&law, state.clone(), base()).unwrap();
    let rotated = framed
        .transport(&law, base().rotated(angle).unwrap(), 0.0, config(), || {
            false
        })
        .unwrap();
    let c = [[0.0013, 0.0004], [0.0004, 0.0007]];
    let basis = [rotate2([1.0, 0.0], angle), rotate2([0.0, 1.0], angle)];
    let c_rot = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            (0..2)
                .map(|r| {
                    (0..2)
                        .map(|s| basis[i][r] * c[r][s] * basis[j][s])
                        .sum::<f64>()
                })
                .sum::<f64>()
        })
    });
    let matrix_config = MatrixTangentialConfig {
        maximum_evaluations: 512,
        displacement_tolerance_m: 2e-13,
        energy_tolerance_j: 2e-12,
        matrix_relative_tolerance: 1e-13,
    };
    let du = [-0.021, 0.008];
    let expected =
        advance_with_matrix_compliance(&law, &state, &input(du), c, matrix_config, || false)
            .unwrap()
            .contact;
    let actual = advance_with_matrix_compliance(
        &law,
        rotated.checkpoint.state(),
        &input(rotate2(du, -angle)),
        c_rot,
        matrix_config,
        || false,
    )
    .unwrap()
    .contact;
    close2(
        actual.traction_n,
        rotate2(expected.traction_n, -angle),
        3e-10,
    );
    close2(
        actual.state.elastic_slip_m(),
        rotate2(expected.state.elastic_slip_m(), -angle),
        2e-12,
    );
    close2(
        actual.state.plastic_slip_m(),
        rotate2(expected.state.plastic_slip_m(), -angle),
        2e-12,
    );
    close(actual.dissipated_j, expected.dissipated_j, 3e-12);
    close(actual.stored_change_j, expected.stored_change_j, 3e-12);
    assert!(actual.dissipated_j > 0.0);
    assert!(FramedTangentialState::new(&law, actual.state, rotated.checkpoint.frame()).is_ok());
}

#[test]
fn antipodal_invalid_foreign_and_cancelled_transfers_are_atomic_and_retryable() {
    let (law, state) = loaded();
    let original = FramedTangentialState::new(&law, state, base()).unwrap();
    let saved = original.clone();
    let antipodal = TangentFrame::new([0.0, 0.0, -1.0], [1.0, 0.0, 0.0]).unwrap();
    assert!(matches!(
        original.transport(&law, antipodal, 0.0, config(), || false),
        Err(TangentialTransportError::Frame(_))
    ));
    for spin in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            original
                .transport(&law, base(), spin, config(), || false)
                .is_err()
        );
    }
    for bad in [
        TangentialTransportConfig {
            absolute_energy_tolerance_j: -1.0,
            ..config()
        },
        TangentialTransportConfig {
            absolute_energy_tolerance_j: f64::INFINITY,
            ..config()
        },
        TangentialTransportConfig {
            relative_energy_tolerance: -1.0,
            ..config()
        },
        TangentialTransportConfig {
            relative_energy_tolerance: f64::NAN,
            ..config()
        },
        TangentialTransportConfig {
            relative_energy_tolerance: 1e-7,
            ..config()
        },
    ] {
        assert!(matches!(
            original.transport(&law, base(), 0.2, bad, || false),
            Err(TangentialTransportError::Invalid(_))
        ));
    }
    for stop in 1..=3 {
        let mut polls = 0;
        assert_eq!(
            original.transport(&law, base(), 0.2, config(), || {
                polls += 1;
                polls == stop
            }),
            Err(TangentialTransportError::Cancelled)
        );
    }
    let foreign = ThirdBodyLayerContact::new(450.0, 0.71, 0.8, 300.0, 0.2).unwrap();
    assert!(
        original
            .transport(&foreign, base(), 0.0, config(), || false)
            .is_err()
    );
    assert_eq!(original, saved);
    let first = original
        .transport(&law, base(), 0.2, config(), || false)
        .unwrap();
    let retry = original
        .transport(&law, base(), 0.2, config(), || false)
        .unwrap();
    assert_eq!(first, retry);
    let exact = TangentialTransportConfig {
        absolute_energy_tolerance_j: 0.0,
        relative_energy_tolerance: 0.0,
    };
    assert_eq!(
        original
            .transport(&law, base(), 0.0, exact, || false)
            .unwrap()
            .checkpoint,
        saved
    );
}

// An adversarial external extension proves the framed energy gate is active;
// it is never used as a substitute for a physical law in the tests above.
struct NonIsometricExtension;
impl TangentialLaw for NonIsometricExtension {
    type State = [f64; 2];
    fn model_id(&self) -> &str {
        "invalid-test-transport"
    }
    fn initial_state(&self) -> Self::State {
        [1.0, 2.0]
    }
    fn stored_energy_j(&self, state: &Self::State) -> f64 {
        state[0] * state[0] + state[1] * state[1]
    }
    fn advance(
        &self,
        _: &Self::State,
        _: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        Err(TriboError::InvalidInput {
            field: "test_has_no_physical_step",
        })
    }
}
impl ObjectiveTangentialLaw for NonIsometricExtension {
    fn transport_state(
        &self,
        state: &Self::State,
        rotation: TangentialHistoryRotation,
    ) -> Result<Self::State, TriboError> {
        if rotation.cosine_sine() == [1.0, 0.0] {
            return Ok(*state);
        }
        Ok(rotation.apply(*state)?.map(|x| 1.01 * x))
    }
}

#[test]
fn an_extension_cannot_publish_undeclared_storage_work_through_frame_transport() {
    let law = NonIsometricExtension;
    let old = FramedTangentialState::new(&law, law.initial_state(), base()).unwrap();
    let error = old
        .transport(&law, base(), 0.4, config(), || false)
        .unwrap_err();
    assert!(matches!(
        error,
        TangentialTransportError::EnergyDrift { .. }
    ));
    assert_eq!(old.state(), &[1.0, 2.0]);
    assert_eq!(old.frame(), base());
    let rotation = TangentialHistoryRotation::between(base(), base(), 0.4).unwrap();
    assert!(rotation.apply([f64::NAN, 1.0]).is_err());
}
