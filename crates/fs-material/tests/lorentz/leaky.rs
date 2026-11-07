//! G0/G1 source identity and independent inertial-polarization/DC-loss roles.
use super::{Entry, card, frequency};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::ELECTRICAL_CONDUCTIVITY_DIMS,
    field_storage::{
        ELECTRIC_PERMITTIVITY_DIMS,
        dielectric::lorentz::{
            LorentzPoleProperties, LorentzState, ResolvedLorentzMaterial,
            leaky::ResolvedLeakyLorentzMaterial,
        },
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_qty::{Dims, QuantitySpec, semantic::FrequencyConvention};

fn entries() -> Vec<Entry> {
    let epsilon = QuantitySpec::dimensional(ELECTRIC_PERMITTIVITY_DIMS);
    let angular = frequency(FrequencyConvention::Angular);
    let rate = QuantitySpec::dimensional(Dims([0, 0, -1, 0, 0, 0]));
    let conductivity = QuantitySpec::dimensional(ELECTRICAL_CONDUCTIVITY_DIMS);
    vec![
        ("epsilon_inf", epsilon, 2e-11),
        ("delta1", epsilon, 3e-11),
        ("omega1", angular, 1000.0),
        ("gamma1", rate, 25.0),
        ("delta2", epsilon, 5e-11),
        ("omega2", angular, 2500.0),
        ("gamma2", rate, 60.0),
        ("sigma_zero", conductivity, 0.0),
        ("sigma", conductivity, 4e-8),
        (
            "density",
            QuantitySpec::dimensional(Dims([-3, 1, 0, 0, 0, 0])),
            2300.0,
        ),
    ]
}

fn names() -> [LorentzPoleProperties; 2] {
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

fn complete_bundle(values: &[Entry]) -> ResolvedMaterialStatePoint {
    let requirements: Vec<_> = values
        .iter()
        .map(|&(name, quantity, _)| {
            ScalarPropertyRequirement::try_with_quantity(
                name,
                quantity,
                ScalarAdmissibility::Finite,
            )
            .unwrap()
        })
        .collect();
    resolve_material_state_point(
        &card(values, 87),
        &QueryPoint::new(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn g0_one_original_bundle_preserves_every_resonant_conductive_and_unrelated_role() {
    let values = entries();
    let source = card(&values, 87);
    let bundle = complete_bundle(&values);
    let dielectric =
        ResolvedLorentzMaterial::try_from_material(bundle.clone(), "epsilon_inf", &names())
            .unwrap();
    let combined = ResolvedLeakyLorentzMaterial::try_from_material(
        bundle.clone(),
        "epsilon_inf",
        &names(),
        "sigma",
    )
    .unwrap();
    let adopted =
        ResolvedLeakyLorentzMaterial::try_from_dielectric(dielectric.clone(), "sigma").unwrap();
    assert_eq!(combined, adopted);
    assert_eq!(combined.material(), &bundle);
    assert_eq!(combined.material().identity(), bundle.identity());
    assert_eq!(combined.material().card_identity(), source.content_hash());
    assert_eq!(combined.material().properties().len(), 10);
    assert_eq!(combined.dielectric(), &dielectric);
    assert_eq!(combined.conductivity_property(), "sigma");
    assert_eq!(combined.conductivity_s_m(), 4e-8);
    let (retained, role, sigma) = combined.into_parts();
    assert_eq!(retained, dielectric);
    assert_eq!(role, "sigma");
    assert_eq!(sigma, 4e-8);

    let resolved = ResolvedLeakyLorentzMaterial::resolve(
        &source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &names(),
        "sigma_zero",
    )
    .unwrap();
    assert_eq!(resolved.conductivity_s_m(), 0.0);
    assert_eq!(resolved.material().properties().len(), 8);
    assert_eq!(resolved.material().card_identity(), source.content_hash());
    assert_eq!(resolved.dielectric().poles(), dielectric.poles());
    assert_eq!(resolved.dielectric().pole_properties(), &names()[..]);
    let empty = ResolvedLeakyLorentzMaterial::resolve(
        &source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &[],
        "sigma",
    )
    .unwrap();
    assert!(empty.dielectric().poles().is_empty());
    assert_eq!(empty.material().properties().len(), 2);
    assert_eq!(empty.conductivity_s_m(), 4e-8);
}

#[test]
fn g1_conduction_does_not_replace_polarization_current_storage_or_pole_heat() {
    let dielectric = ResolvedLorentzMaterial::try_from_material(
        complete_bundle(&entries()),
        "epsilon_inf",
        &names(),
    )
    .unwrap();
    let before = [
        LorentzState {
            polarization_c_m2: [1e-7, -2e-7, 3e-8],
            current_a_m2: [2e-4, -3e-4, 1e-4],
        },
        LorentzState {
            polarization_c_m2: [-3e-8, 4e-8, 1e-7],
            current_a_m2: [1e-4, 2e-4, -3e-4],
        },
    ];
    let (dt, e0, e1) = (1e-4, [2e4, -3e4, 0.5e4], [-1e4, 2e4, 4e4]);
    let expected = dielectric.step_midpoint(dt, e0, e1, &before).unwrap();
    assert!(expected.heat_j_m3 > 0.0);
    for name in ["sigma_zero", "sigma"] {
        let law =
            ResolvedLeakyLorentzMaterial::try_from_dielectric(dielectric.clone(), name).unwrap();
        assert_eq!(law.dielectric(), &dielectric);
        let prepared = law.dielectric().prepare_midpoint(dt).unwrap();
        assert_eq!(prepared.step_midpoint(e0, e1, &before).unwrap(), expected);
        assert_eq!(
            law.dielectric()
                .stored_energy_density_j_m3(e0, &before)
                .unwrap(),
            expected.stored_before_j_m3,
        );
        assert_eq!(law.material(), dielectric.material());
    }

    // A conductive material may have exactly conservative bound-charge poles.
    // Its DC Joule heat belongs to the field's conductor, not to the retained
    // polarization current or its independently evaluated friction account.
    let mut lossless = entries();
    lossless[3].2 = 0.0;
    lossless[6].2 = 0.0;
    let law = ResolvedLeakyLorentzMaterial::try_from_material(
        complete_bundle(&lossless),
        "epsilon_inf",
        &names(),
        "sigma",
    )
    .unwrap();
    assert!(law.conductivity_s_m() > 0.0);
    let step = law.dielectric().step_midpoint(dt, e0, e1, &before).unwrap();
    assert_eq!(step.heat_j_m3, 0.0);
    assert!(step.states_after.iter().any(|p| p.current_a_m2 != [0.0; 3]));
}

#[test]
fn g0_aliased_missing_wrong_quantity_negative_and_dynamic_conductivity_refuse() {
    let dielectric = ResolvedLorentzMaterial::try_from_material(
        complete_bundle(&entries()),
        "epsilon_inf",
        &names(),
    )
    .unwrap();
    for role in [
        "",
        " ",
        "missing",
        "epsilon_inf",
        "delta1",
        "omega1",
        "gamma1",
        "delta2",
        "omega2",
        "gamma2",
        "density",
    ] {
        assert!(
            ResolvedLeakyLorentzMaterial::try_from_dielectric(dielectric.clone(), role).is_err()
        );
    }
    let source = card(&entries(), 87);
    let without_dc = ResolvedLorentzMaterial::resolve(
        &source,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        "epsilon_inf",
        &names(),
    )
    .unwrap();
    // The card contains sigma, but adoption cannot silently make a second
    // resolution and replace the complete already accepted snapshot.
    assert!(ResolvedLeakyLorentzMaterial::try_from_dielectric(without_dc, "sigma").is_err());

    for wrong_units in [false, true] {
        let mut values = entries();
        if wrong_units {
            values[8].1 = QuantitySpec::dimensional(ELECTRIC_PERMITTIVITY_DIMS);
        } else {
            values[8].2 = -1.0;
        }
        assert!(
            ResolvedLeakyLorentzMaterial::try_from_material(
                complete_bundle(&values),
                "epsilon_inf",
                &names(),
                "sigma",
            )
            .is_err()
        );
    }
    let mut wrong_frequency = entries();
    wrong_frequency[2].1 = frequency(FrequencyConvention::Cyclic);
    assert!(
        ResolvedLeakyLorentzMaterial::resolve(
            &card(&wrong_frequency, 87),
            &QueryPoint::new(),
            MaterialPropertySelection::SingleClaimOnly,
            "epsilon_inf",
            &names(),
            "sigma",
        )
        .is_err()
    );

    let cyclic = frequency(FrequencyConvention::Cyclic);
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        ("epsilon_inf", ELECTRIC_PERMITTIVITY_DIMS, 2e-11),
        ("sigma", ELECTRICAL_CONDUCTIVITY_DIMS, 4e-8),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: if name == "sigma" {
                    ValidityDomain::unconstrained().with_quantity("frequency", cyclic, 1.0, 100.0)
                } else {
                    ValidityDomain::unconstrained()
                },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic conductivity applicability adversary".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let restricted = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic frequency-restricted conductor".into(),
            phase: "solid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with_quantity("frequency", cyclic, 10.0)
        .unwrap();
    assert!(
        ResolvedLeakyLorentzMaterial::resolve(
            &restricted,
            &query,
            MaterialPropertySelection::SingleClaimOnly,
            "epsilon_inf",
            &[],
            "sigma",
        )
        .is_err()
    );
}
