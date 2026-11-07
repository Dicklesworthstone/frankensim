//! G0/G1 source-bound Lorentz polarization, using synthetic material claims.
//! These checks establish source/quantity binding and constitutive arithmetic,
//! not measured-material accuracy or a transient Maxwell discretization.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, UncertaintyModel,
};
use fs_qty::{
    QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};

type Entry = (&'static str, QuantitySpec, f64);

fn frequency(convention: FrequencyConvention) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Frequency(convention),
        ValueForm::Static,
    ))
}

fn card(entries: &[Entry], specimen: u64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for &(name, quantity, value) in entries {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Scalar {
                    value,
                    dims: quantity.dims(),
                },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic bound-charge resonance fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("synthetic Lorentz dielectric {specimen}"),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[cfg(test)]
mod integration {
    use super::*;
    use fs_matdb::QueryPoint;
    use fs_material::{
        field_storage::{
            ELECTRIC_PERMITTIVITY_DIMS, FieldStorageError, MAGNETIC_PERMEABILITY_DIMS,
            MAGNETIC_PERMEABILITY_PROPERTY,
            dielectric::lorentz::{
                LorentzError, LorentzPoleProperties, LorentzState, ResolvedLorentzMaterial,
            },
        },
        state_point::{
            MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
            resolve_material_state_point,
        },
    };
    use fs_qty::Dims;

    const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);

    fn entries(angular_frequency: QuantitySpec, damping: f64) -> Vec<Entry> {
        let epsilon = QuantitySpec::dimensional(ELECTRIC_PERMITTIVITY_DIMS);
        vec![
            ("epsilon_inf", epsilon, 2e-11),
            ("delta1", epsilon, 3e-11),
            ("omega1", angular_frequency, 1000.0),
            ("gamma1", QuantitySpec::dimensional(RATE), damping),
            ("delta2", epsilon, 5e-11),
            ("omega2", angular_frequency, 2500.0),
            ("gamma2", QuantitySpec::dimensional(RATE), 60.0),
        ]
    }

    fn properties() -> [LorentzPoleProperties; 2] {
        [
            LorentzPoleProperties {
                permittivity_increment: "delta1".into(),
                angular_frequency: "omega1".into(),
                damping_rate: "gamma1".into(),
            },
            LorentzPoleProperties {
                permittivity_increment: "delta2".into(),
                angular_frequency: "omega2".into(),
                damping_rate: "gamma2".into(),
            },
        ]
    }

    fn resolve(
        source: &MaterialCard,
        properties: &[LorentzPoleProperties],
    ) -> Result<ResolvedLorentzMaterial, FieldStorageError> {
        ResolvedLorentzMaterial::resolve(
            source,
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            "epsilon_inf",
            properties,
        )
    }

    fn near(actual: f64, expected: f64, scale: f64) {
        assert!(
            (actual - expected).abs() <= 3e-13 * scale.max(1e-300),
            "actual={actual:.17e}, expected={expected:.17e}"
        );
    }

    #[test]
    fn g1_source_bound_one_and_two_pole_response_matches_original_lorentz_law() {
        let source = card(&entries(frequency(FrequencyConvention::Angular), 25.0), 0);
        let names = properties();
        for count in [1, 2] {
            let law = resolve(&source, &names[..count]).unwrap();
            assert_eq!(law.material().card_identity(), source.content_hash());
            assert_eq!(law.material().properties().len(), 1 + 3 * count);
            assert_eq!(law.high_frequency_property(), "epsilon_inf");
            assert_eq!(law.pole_properties(), &names[..count]);
            assert_eq!(law.high_frequency_permittivity_f_m(), 2e-11);
            assert_eq!(law.poles()[0].permittivity_increment_f_m(), 3e-11);
            assert_eq!(law.poles()[0].angular_frequency_rad_s(), 1000.0);
            assert_eq!(law.poles()[0].damping_rate_per_s(), 25.0);
            for omega in [0.0_f64, 500.0, 1000.0, 1500.0, 2500.0, 1e5] {
                let mut expected_re = 2e-11;
                let mut expected_im = 0.0;
                // Independent original-coordinate formula, with no call to the
                // production pole response or its normalization.
                for &(delta, natural, damping) in
                    &[(3e-11, 1000.0_f64, 25.0), (5e-11, 2500.0, 60.0)][..count]
                {
                    let detuning = natural * natural - omega * omega;
                    let friction = damping * omega;
                    let factor =
                        delta * natural * natural / (detuning * detuning + friction * friction);
                    expected_re += factor * detuning;
                    expected_im -= factor * friction;
                }
                let actual = law.permittivity_f_m(omega).unwrap();
                near(actual.re, expected_re, expected_re.abs());
                near(actual.im, expected_im, expected_im.abs());
                assert!(actual.im <= 0.0, "exp(+i omega t) must retain passive loss");
            }
        }
    }

    #[test]
    fn g0_hertz_and_untyped_inverse_seconds_are_not_angular_frequencies() {
        for wrong in [
            frequency(FrequencyConvention::Cyclic),
            QuantitySpec::dimensional(RATE),
        ] {
            let source = card(&entries(wrong, 25.0), 0);
            assert!(resolve(&source, &properties()[..1]).is_err());
        }
        let source = card(&entries(frequency(FrequencyConvention::Angular), 25.0), 0);
        assert!(resolve(&source, &properties()[..1]).is_ok());
    }

    #[test]
    fn g0_missing_damping_and_aliased_source_roles_are_refused() {
        let mut values = entries(frequency(FrequencyConvention::Angular), 25.0);
        values.retain(|(name, _, _)| *name != "gamma1");
        let missing_damping = card(&values, 0);
        assert!(resolve(&missing_damping, &properties()[..1]).is_err());

        let source = card(&entries(frequency(FrequencyConvention::Angular), 25.0), 0);
        let mut names = properties();
        names[0].angular_frequency = names[0].damping_rate.clone();
        assert!(resolve(&source, &names[..1]).is_err());
        assert!(resolve(&source, &[properties()[0].clone(), properties()[0].clone()]).is_err());
    }

    #[test]
    fn g1_explicit_zero_damping_preserves_lossless_dynamics_and_refuses_a_true_pole() {
        let source = card(&entries(frequency(FrequencyConvention::Angular), 0.0), 0);
        let law = resolve(&source, &properties()[..1]).unwrap();
        assert_eq!(law.poles()[0].damping_rate_per_s(), 0.0);
        assert!(matches!(
            law.permittivity_f_m(1000.0),
            Err(FieldStorageError::Invalid(_))
        ));
        assert_eq!(
            law.poles()[0].susceptibility_f_m(1000.0),
            Err(LorentzError::SingularResonance)
        );
        let response = law.permittivity_f_m(500.0).unwrap();
        near(response.re, 2e-11 + 4.0 * 3e-11 / 3.0, 1e-10);
        assert_eq!(response.im, 0.0);

        let prepared = law.prepare_step(0.001).unwrap();
        assert_eq!(prepared.len(), 1);
        let step = prepared[0]
            .advance(LorentzState::<1>::default(), [1.5], [1.5])
            .unwrap();
        assert!(step.stored_after_j_m3 > 0.0);
        assert_eq!(step.heat_j_m3, 0.0);
        near(
            step.stored_after_j_m3,
            step.supplied_work_j_m3,
            step.stored_after_j_m3,
        );

        let instantaneous = resolve(&source, &[]).unwrap();
        assert!(instantaneous.prepare_step(0.001).unwrap().is_empty());
        let epsilon = instantaneous.permittivity_f_m(1000.0).unwrap();
        assert_eq!(epsilon.re, 2e-11);
        assert_eq!(epsilon.im, 0.0);
    }

    #[test]
    fn g0_complete_source_identity_and_other_physics_claims_are_preserved() {
        let mut values = entries(frequency(FrequencyConvention::Angular), 25.0);
        values.push((
            MAGNETIC_PERMEABILITY_PROPERTY,
            QuantitySpec::dimensional(MAGNETIC_PERMEABILITY_DIMS),
            4e-7 * core::f64::consts::PI,
        ));
        let source = card(&values, 0);
        let requirements = values
            .iter()
            .map(|&(name, quantity, _)| {
                ScalarPropertyRequirement::try_with_quantity(
                    name,
                    quantity,
                    ScalarAdmissibility::NonNegative,
                )
                .unwrap()
            })
            .collect::<Vec<_>>();
        let complete = resolve_material_state_point(
            &source,
            &QueryPoint::new(),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let law = ResolvedLorentzMaterial::try_from_material(
            complete.clone(),
            "epsilon_inf",
            &properties(),
        )
        .unwrap();
        assert_eq!(law.material(), &complete);
        assert_eq!(law.material().identity(), complete.identity());
        assert!(
            law.material()
                .property(MAGNETIC_PERMEABILITY_PROPERTY)
                .is_some()
        );

        // Equal numerical parameters from a distinct specimen source do
        // not establish the same material identity for multiphysics composition.
        let other_source = card(&values, 1);
        let other = resolve(&other_source, &properties()).unwrap();
        assert_ne!(
            law.material().card_identity(),
            other.material().card_identity()
        );
        assert_ne!(law.material().identity(), other.material().identity());
        assert_eq!(
            law.permittivity_f_m(1200.0).unwrap(),
            other.permittivity_f_m(1200.0).unwrap()
        );
    }
}
