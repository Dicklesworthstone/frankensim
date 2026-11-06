//! Original static claims select the held branch without advancing history.

use super::*;
use fs_solid::brittle::solid::dynamics::modes::FrozenBrittleModesControls;
use fs_solid::brittle::solid::dynamics::modes::contact::FrozenBrittleContactControls;
use fs_solid::brittle::solid::dynamics::modes::contact::friction::FrozenBrittleStickingControls;
use fs_solid::brittle::solid::dynamics::modes::contact::friction::source::FrozenBrittleSourcedStickingError;

fn modal_controls() -> FrozenBrittleStickingControls {
    FrozenBrittleStickingControls {
        contact: FrozenBrittleContactControls {
            assembly: FrozenBrittleModesControls {
                maximum_dofs: 24,
                maximum_staged_entries: 288,
            },
            maximum_free_dofs: 24,
            maximum_contact_rows: 9,
            maximum_dense_entries: 65_536,
            ..FrozenBrittleContactControls::default()
        },
        minimum_friction_margin: 1e-6,
    }
}

#[test]
fn g3_sourced_modes_use_static_temperature_claims_and_preserve_original_history() {
    with_oblique_motion(&NODES, [0.0, 0.0, 1.0], [[0.0; 3]; 2], |motion, cx| {
        let law = law(motion.solid().materials()[0].material().material());
        let history = history(&law);
        let original_history = history.clone();
        let original_body = motion.solid().solution().clone();
        let hot = conditions();
        let cold = SourceFrictionConditions {
            query: hot.query.clone().with("T", 300.0).unwrap(),
            contact_temperature_k: Some(300.0),
            ..hot.clone()
        };
        let vertex = BrittleVertexFrictionSource {
            law: &law,
            state: &history,
            conditions: &cold,
            coefficient_role: SourceImpulseCoefficientRole::Static,
        };
        let binding = BrittleFaceFrictionSource {
            face_id: 951,
            tangent_direction: [1.0, 0.0, 0.0],
            vertices: [vertex; 3],
        };
        // At 300 K the first static coefficient is 0.6, whereas at 320 K
        // it is 0.5. The same physical shear/normal ratio 0.55 changes branch.
        let preload = [[[1.0, 0.55, 0.0]; 3]];
        let observed = motion
            .frozen_sourced_sticking_contact_modes(&[binding], &preload, &[], modal_controls(), cx)
            .unwrap();
        assert_eq!(observed.samples.len(), 3);
        assert_eq!(observed.contact.dimension(), 15);
        close(observed.contact.pencil().reference_mass_kg(), 2.0, 1e-14);
        for (vertex_index, sample) in observed.samples.iter().enumerate() {
            assert_eq!(sample.row.face_id, 951);
            assert_eq!(sample.row.vertex, vertex_index);
            assert_eq!(
                sample.sample.coefficient_role,
                SourceImpulseCoefficientRole::Static
            );
            close(sample.sample.coefficients_mu[0], 0.6, 1e-14);
            close(sample.sample.coefficients_mu[1], 0.7, 1e-14);
            assert_eq!(sample.sample.queried_speed_m_s, 0.0);
            let old = history.constitutive_state();
            let next = sample.sample.state.constitutive_state();
            assert_eq!(next.slip_m(), old.slip_m());
            assert_eq!(next.principal_axis(), old.principal_axis());
            assert_eq!(next.contact_temperature_k(), old.contact_temperature_k());
            for source in sample.sample.state.friction_sources() {
                assert_eq!(source.state_coordinate("speed"), Some(0.0));
                close(source.state_coordinate("path").unwrap(), 0.2525, 1e-14);
            }
        }
        let explicit = [BrittleFaceFriction {
            tangent_direction: binding.tangent_direction,
            vertices: core::array::from_fn(|i| EllipticFriction {
                coefficients_mu: observed.samples[i].sample.coefficients_mu,
                principal_axis: observed.samples[i]
                    .sample
                    .state
                    .constitutive_state()
                    .principal_axis(),
            }),
        }];
        let numeric = motion
            .frozen_sticking_contact_modes(&explicit, &preload, &[], modal_controls(), cx)
            .unwrap();
        assert_eq!(
            observed.contact.projected_stiffness_n_m(),
            numeric.projected_stiffness_n_m()
        );
        assert_eq!(
            observed.contact.projected_mass_kg(),
            numeric.projected_mass_kg()
        );
        for (a, b) in observed.contact.modes().iter().zip(numeric.modes()) {
            assert_eq!(a.coordinates.lambda, b.coordinates.lambda);
        }
        let hot_binding = BrittleFaceFrictionSource {
            vertices: [BrittleVertexFrictionSource {
                conditions: &hot,
                ..vertex
            }; 3],
            ..binding
        };
        assert!(matches!(
            motion.frozen_sourced_sticking_contact_modes(
                &[hot_binding],
                &preload,
                &[],
                modal_controls(),
                cx,
            ),
            Err(FrozenBrittleSourcedStickingError::Contact(_))
        ));
        assert_eq!(history, original_history);
        assert_eq!(motion.solid().solution(), &original_body);
        assert_eq!(motion.accepted_steps(), 1);
        assert_eq!(motion.remaining_steps(), 1);
        let retry = motion
            .frozen_sourced_sticking_contact_modes(&[binding], &preload, &[], modal_controls(), cx)
            .unwrap();
        assert_eq!(retry.samples[0].sample, observed.samples[0].sample);
    });
}

#[test]
fn g0_sourced_modes_refuse_wrong_roles_materials_and_late_source_queries_atomically() {
    with_oblique_motion(&NODES, [0.0, 0.0, 1.0], [[0.0; 3]; 2], |motion, cx| {
        let material = motion.solid().materials()[0].material().material();
        let law = law(material);
        let history = history(&law);
        let before = history.clone();
        let conditions = conditions();
        let vertex = BrittleVertexFrictionSource {
            law: &law,
            state: &history,
            conditions: &conditions,
            coefficient_role: SourceImpulseCoefficientRole::Static,
        };
        let binding = BrittleFaceFrictionSource {
            face_id: 951,
            tangent_direction: [1.0, 0.0, 0.0],
            vertices: [vertex; 3],
        };
        let preload = [[[1.0, 0.0, 0.0]; 3]];
        let kinetic = BrittleFaceFrictionSource {
            vertices: [BrittleVertexFrictionSource {
                coefficient_role: SourceImpulseCoefficientRole::Kinetic,
                ..vertex
            }; 3],
            ..binding
        };
        assert!(matches!(
            motion.frozen_sourced_sticking_contact_modes(
                &[kinetic],
                &preload,
                &[],
                modal_controls(),
                cx,
            ),
            Err(FrozenBrittleSourcedStickingError::Source(
                BrittleFaceSourcedFrictionError::Invalid(_)
            ))
        ));
        let mut unrelated = material.clone();
        unrelated.revision += 1;
        let wrong_law = super::law(&unrelated);
        let wrong_history = wrong_law.initial_state();
        let wrong = BrittleFaceFrictionSource {
            vertices: [BrittleVertexFrictionSource {
                law: &wrong_law,
                state: &wrong_history,
                ..vertex
            }; 3],
            ..binding
        };
        assert!(matches!(
            motion.frozen_sourced_sticking_contact_modes(
                &[wrong],
                &preload,
                &[],
                modal_controls(),
                cx,
            ),
            Err(FrozenBrittleSourcedStickingError::Source(
                BrittleFaceSourcedFrictionError::Binding { .. }
            ))
        ));
        let invalid = SourceFrictionConditions {
            query: conditions.query.clone().with("T", 330.0).unwrap(),
            contact_temperature_k: Some(330.0),
            ..conditions.clone()
        };
        let mut late = binding;
        late.vertices[2].conditions = &invalid;
        assert!(matches!(motion.frozen_sourced_sticking_contact_modes(
            &[late], &preload, &[], modal_controls(), cx,
        ), Err(FrozenBrittleSourcedStickingError::Source(BrittleFaceSourcedFrictionError::Source { row, .. })) if row.vertex == 2));
        assert_eq!(history, before);
        assert!(
            motion
                .frozen_sourced_sticking_contact_modes(
                    &[binding],
                    &preload,
                    &[],
                    modal_controls(),
                    cx,
                )
                .is_ok()
        );

        // The same invalid conditions must never be queried on open faces.
        let open = crate::approach::open_motion(motion.solid(), 1e-4, 0.0, 2, cx);
        let inactive = BrittleFaceFrictionSource {
            tangent_direction: [f64::NAN; 3],
            vertices: [BrittleVertexFrictionSource {
                conditions: &invalid,
                ..vertex
            }; 3],
            ..binding
        };
        let observed = open
            .frozen_sourced_sticking_contact_modes(
                &[inactive],
                &[[[0.0; 3]; 3]],
                &[],
                modal_controls(),
                cx,
            )
            .unwrap();
        assert!(observed.samples.is_empty());
        assert!(observed.contact.constraints().is_empty());
        assert_eq!(observed.contact.dimension(), 24);
        assert_eq!(history, before);
    });
}

#[test]
fn g3_source_queries_measured_small_slip_instead_of_template_speed() {
    with_oblique_motion(
        &NODES,
        [0.0, 0.0, 1.0],
        [[0.0; 3], [3e-9, 4e-9, 0.0]],
        |motion, cx| {
            let law = law(motion.solid().materials()[0].material().material());
            let history = history(&law);
            let conditions = conditions();
            let binding = BrittleFaceFrictionSource {
                face_id: 951,
                tangent_direction: [3.0, 4.0, 0.0],
                vertices: [BrittleVertexFrictionSource {
                    law: &law,
                    state: &history,
                    conditions: &conditions,
                    coefficient_role: SourceImpulseCoefficientRole::Static,
                }; 3],
            };
            let result = motion
                .frozen_sourced_sticking_contact_modes(
                    &[binding],
                    &[[[1.0, 0.0, 0.0]; 3]],
                    &[],
                    modal_controls(),
                    cx,
                )
                .unwrap();
            for sample in result.samples {
                close(sample.sample.queried_speed_m_s, 5e-9, 2e-14);
                for source in sample.sample.state.friction_sources() {
                    close(source.state_coordinate("speed").unwrap(), 5e-9, 2e-14);
                    close(source.state_coordinate("path").unwrap(), 0.2525, 1e-14);
                }
            }
        },
    );
}
