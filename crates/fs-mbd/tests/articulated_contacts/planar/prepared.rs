//! G1/G3/G4: reuse the original full native response without aging histories.

use super::*;
use fs_mbd::articulated_contacts::tangential::planar::multiple::{
    PreparedArticulatedTangentialContacts, PreparedArticulatedTangentialInput,
};
use fs_mbd::articulated_contacts::tangential::planar::transport::{
    ArticulatedFramedTangentialSite, PreparedArticulatedFramedInput,
    PreparedArticulatedFramedTangentialContacts, advance_articulated_framed_tangential_contacts,
};
use fs_tribo::layer_contact::ThirdBodyLayerState;
use fs_tribo::tangential::transport::{
    FramedTangentialState, ObjectiveTangentialLaw, TangentialHistoryRotation,
    TangentialTransportConfig,
};
use fs_tribo::{TangentialStep, TriboError};
use std::cell::Cell;

fn inputs<'a, L: TangentialLaw>(
    laws: &'a [L],
    histories: &'a [L::State],
    dt: f64,
) -> Vec<PreparedArticulatedTangentialInput<'a, L>> {
    laws.iter()
        .zip(histories)
        .map(|(law, state)| PreparedArticulatedTangentialInput {
            law,
            state,
            drive: drive(dt),
        })
        .collect()
}

#[test]
fn prepared_trials_match_independent_full_elastic_recoil_and_borrow_the_same_matrices() {
    let model = body(2.0, 1.0);
    let initial = system(
        &model,
        Vec3::default(),
        Vec3::new(0.2, -0.1, 0.05),
        Vec3::new(0.1, 0.2, -0.1),
    );
    let contacts = opposing();
    let laws = layers();
    let old_e = [[0.004, -0.002], [-0.003, 0.005]];
    let histories = [preload(&laws[0], old_e[0]), preload(&laws[1], old_e[1])];
    let dt = 0.02;
    let mut contexts = [initial];
    let prepared =
        PreparedArticulatedTangentialContacts::new(&contexts, &contacts, dt, controls(), || false)
            .unwrap();
    // The native contexts are snapshots. Replacing the caller's container
    // cannot silently change the poses/twists associated with cached response.
    contexts[0].base.twist_body.linear = Vec3::new(9.0, 8.0, 7.0);
    let w = mobility(&contacts, 0.5, 1.0);
    let expected = elastic_reference(
        &w,
        &released(&contacts, initial.base.twist_body, dt),
        &[400.0, 700.0],
        &old_e,
        dt,
    );
    let trial_inputs = inputs(&laws, &histories, dt);
    let mut first = None;
    for _ in 0..3 {
        let trial = prepared.advance(&trial_inputs, || false).unwrap();
        assert_eq!(
            trial.mobility_kg_inverse.as_ptr(),
            prepared.mobility_kg_inverse().as_ptr()
        );
        assert_eq!(
            trial.compliance_m_per_n.as_ptr(),
            prepared.compliance_m_per_n().as_ptr()
        );
        for i in 0..2 {
            close2(
                trial.contacts[i].traction_n,
                [expected[2 * i], expected[2 * i + 1]],
            );
            close_vector(
                prepared.contact_world_point(i).unwrap(),
                contacts[i].point.point_body,
            );
        }
        let owned = trial.into_owned();
        audit_all(&owned);
        if let Some((contacts, body)) = &first {
            assert_eq!(&owned.contacts, contacts);
            assert_eq!(&owned.responses[0].base_after, body);
        } else {
            first = Some((owned.contacts, owned.responses[0].base_after));
        }
    }
    assert_eq!(prepared.contact_world_point(usize::MAX), None);
    assert_eq!(prepared.contacts(), &contacts);
    assert_eq!(prepared.duration_s(), dt);
    assert_eq!(prepared.config(), controls());
    assert_eq!(histories[0], preload(&laws[0], old_e[0]));
    assert_eq!(histories[1], preload(&laws[1], old_e[1]));
}

// Instrument only the original layer's objective transport; every constitutive
// evaluation and retained state remain the original physical implementation.
struct CountedLayer<'a> {
    layer: &'a ThirdBodyLayerContact,
    rotations: &'a Cell<usize>,
}
impl TangentialLaw for CountedLayer<'_> {
    type State = ThirdBodyLayerState;
    fn model_id(&self) -> &str {
        self.layer.model_id()
    }
    fn initial_state(&self) -> Self::State {
        self.layer.initial_state()
    }
    fn stored_energy_j(&self, state: &Self::State) -> f64 {
        self.layer.stored_energy_j(state)
    }
    fn advance(
        &self,
        state: &Self::State,
        input: &TangentialStepInput,
    ) -> Result<TangentialStep<Self::State>, TriboError> {
        self.layer.advance(state, input)
    }
}
impl ObjectiveTangentialLaw for CountedLayer<'_> {
    fn transport_state(
        &self,
        state: &Self::State,
        rotation: TangentialHistoryRotation,
    ) -> Result<Self::State, TriboError> {
        if rotation.cosine_sine() != [1.0, 0.0] {
            self.rotations.set(self.rotations.get() + 1);
        }
        self.layer.transport_state(state, rotation)
    }
}

#[test]
fn framed_trials_transport_each_native_history_once_and_match_the_original_proposal() {
    let model = body(2.0, 1.0);
    let initial = system(
        &model,
        Vec3::default(),
        Vec3::new(0.2, -0.1, 0.0),
        Vec3::default(),
    );
    let physical = layers();
    let rotations = Cell::new(0);
    let laws = physical.each_ref().map(|layer| CountedLayer {
        layer,
        rotations: &rotations,
    });
    let mut contacts = opposing();
    let histories = [
        FramedTangentialState::new(
            &laws[0],
            preload(&laws[0], [0.004, -0.002]),
            contacts[0].frame_world,
        )
        .unwrap(),
        FramedTangentialState::new(
            &laws[1],
            preload(&laws[1], [-0.003, 0.005]),
            contacts[1].frame_world,
        )
        .unwrap(),
    ];
    let saved = histories.clone();
    contacts[0].frame_world = contacts[0].frame_world.rotated(0.4).unwrap();
    contacts[1].frame_world = contacts[1].frame_world.rotated(-0.2).unwrap();
    let sites: Vec<_> = (0..2)
        .map(|i| ArticulatedFramedTangentialSite {
            contact: contacts[i],
            law: &laws[i],
            history: &histories[i],
            drive: drive(0.02),
            common_spin_angle_rad: [-0.3, 0.5][i],
        })
        .collect();
    let transport = TangentialTransportConfig {
        absolute_energy_tolerance_j: 1e-13,
        relative_energy_tolerance: 1e-12,
    };
    let prepared = PreparedArticulatedFramedTangentialContacts::new(
        &[initial],
        &sites,
        controls(),
        transport,
        || false,
    )
    .unwrap();
    assert_eq!(rotations.get(), 2);
    let transported = prepared.transported_histories().to_vec();
    let trial_inputs: Vec<_> = laws
        .iter()
        .map(|law| PreparedArticulatedFramedInput {
            law,
            drive: drive(0.02),
        })
        .collect();
    let mut last = None;
    for _ in 0..3 {
        let trial = prepared.advance(&trial_inputs, || false).unwrap();
        assert_eq!(rotations.get(), 2);
        assert_eq!(prepared.transported_histories(), &transported);
        assert_eq!(
            trial.mechanics.compliance_m_per_n.as_ptr(),
            prepared.mechanics().compliance_m_per_n().as_ptr()
        );
        close(trial.accepted_energy_balance_residual_j, 0.0);
        last = Some(trial.into_owned());
    }
    let owned = last.unwrap();
    drop(prepared);
    // Owned final reports retain the original matrix and all framed states.
    let expected = advance_articulated_framed_tangential_contacts(
        &[initial],
        &sites,
        controls(),
        transport,
        || false,
    )
    .unwrap();
    assert_eq!(owned.mechanics.contacts, expected.mechanics.contacts);
    assert_eq!(
        owned.mechanics.mobility_kg_inverse,
        expected.mechanics.mobility_kg_inverse
    );
    assert_eq!(
        owned.mechanics.responses[0].base_after,
        expected.mechanics.responses[0].base_after
    );
    assert_eq!(
        owned.transport_energy_residuals_j,
        expected.transport_energy_residuals_j
    );
    assert_eq!(histories, saved);
}

#[test]
fn prepared_native_and_framed_owners_retry_exactly_after_late_failure_or_cancellation() {
    let model = body(2.0, 1.0);
    let initial = system(&model, Vec3::default(), Vec3::default(), Vec3::default());
    let contacts = opposing();
    let laws = layers();
    let histories = [
        preload(&laws[0], [0.004, -0.002]),
        preload(&laws[1], [-0.003, 0.005]),
    ];
    let mut config = controls();
    config.mechanical.maximum_slip_speed_m_s = 0.4;
    config.mechanical.solve.maximum_force_n = 1000.0;
    let prepared =
        PreparedArticulatedTangentialContacts::new(&[initial], &contacts, 0.02, config, || false)
            .unwrap();
    let mut trial_inputs = inputs(&laws, &histories, 0.02);
    let mut polls = 0;
    let expected = prepared
        .advance(&trial_inputs, || {
            polls += 1;
            false
        })
        .unwrap()
        .into_owned();
    for stop in [1, polls / 2, polls] {
        let mut seen = 0;
        assert!(matches!(
            prepared.advance(&trial_inputs, || {
                seen += 1;
                seen == stop
            }),
            Err(ArticulatedContactError::Cancelled)
        ));
    }
    assert!(prepared.advance(&trial_inputs[..1], || false).is_err());
    trial_inputs[1].drive.dt_s = 0.01;
    assert!(prepared.advance(&trial_inputs, || false).is_err());
    trial_inputs[1].drive.dt_s = 0.02;
    let energetic = [preload(&laws[0], [0.1, 0.1]), preload(&laws[1], [0.1, 0.1])];
    assert!(matches!(
        prepared.advance(&inputs(&laws, &energetic, 0.02), || false),
        Err(ArticulatedContactError::Unresolved(_))
    ));
    let retry = prepared.advance(&trial_inputs, || false).unwrap();
    assert_eq!(retry.contacts, expected.contacts);
    assert_eq!(
        retry.responses[0].base_after,
        expected.responses[0].base_after
    );

    let framed: Vec<_> = (0..2)
        .map(|i| {
            FramedTangentialState::new(&laws[i], histories[i].clone(), contacts[i].frame_world)
                .unwrap()
        })
        .collect();
    let sites: Vec<_> = (0..2)
        .map(|i| ArticulatedFramedTangentialSite {
            contact: contacts[i],
            law: &laws[i],
            history: &framed[i],
            drive: drive(0.02),
            common_spin_angle_rad: 0.2,
        })
        .collect();
    let framed_prepared = PreparedArticulatedFramedTangentialContacts::new(
        &[initial],
        &sites,
        config,
        TangentialTransportConfig {
            absolute_energy_tolerance_j: 1e-13,
            relative_energy_tolerance: 1e-12,
        },
        || false,
    )
    .unwrap();
    let mut inputs: Vec<_> = laws
        .iter()
        .map(|law| PreparedArticulatedFramedInput {
            law,
            drive: drive(0.02),
        })
        .collect();
    let mut polls = 0;
    let baseline = framed_prepared
        .advance(&inputs, || {
            polls += 1;
            false
        })
        .unwrap()
        .into_owned();
    let mut seen = 0;
    assert!(matches!(
        framed_prepared.advance(&inputs, || {
            seen += 1;
            seen == polls
        }),
        Err(ArticulatedContactError::Cancelled)
    ));
    inputs[1].drive.contact_temperature_k = Some(301.0);
    assert!(framed_prepared.advance(&inputs, || false).is_err());
    inputs[1].drive.contact_temperature_k = Some(300.0);
    let replay = framed_prepared.advance(&inputs, || false).unwrap();
    assert_eq!(replay.mechanics.contacts, baseline.mechanics.contacts);
    assert_eq!(
        replay.mechanics.responses[0].base_after,
        baseline.mechanics.responses[0].base_after
    );
}
