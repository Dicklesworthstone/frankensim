//! G0/G3: ordered interface data reach the native fracture law without defaults.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::source::*;
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS, ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
    ResolvedElectricalInterface,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedInterfaceStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, interface_stiffness_per_area_dims, resolve_interface_state_point,
};
use fs_qty::{Dims, Length, Pressure, Temperature};

fn claim(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 340.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured cohesive interface; not measured adhesive data".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn properties() -> Vec<PropertyClaim> {
    vec![
        claim(
            COHESIVE_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            1e9,
        ),
        claim(COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, 1e6),
        claim(
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            2000.0,
        ),
        claim(
            COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            2e9,
        ),
    ]
}

fn card(properties: Vec<PropertyClaim>, history: &str) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for property in properties {
        claims.insert_claim(property).unwrap();
    }
    let surface = |chemistry: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: chemistry.into(),
            phase: "solid".into(),
            process: "manufactured surface".into(),
            revision: 1,
        },
        texture_frame: format!("{chemistry} reference surface"),
    };
    InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "explicit cohesive bond".into(),
            third_body: None,
            environment: "fixed external conditions".into(),
            history: history.into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn resolve(card: &InterfaceSystemCard) -> ResolvedBilinearCohesive {
    resolve_bilinear_cohesive(
        card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        0.02,
    )
    .unwrap()
}

fn bundle(card: &InterfaceSystemCard, properties: &[PropertyClaim]) -> ResolvedInterfaceStatePoint {
    let requirements: Vec<_> = properties
        .iter()
        .map(|p| {
            ScalarPropertyRequirement::try_with_key(&p.key, ScalarAdmissibility::Finite).unwrap()
        })
        .collect();
    resolve_interface_state_point(
        card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn ordered_source_bundle_is_shared_by_failure_and_electrical_response() {
    let mut properties = properties();
    properties.push(claim(
        ELECTRICAL_CONTACT_RESISTANCE_PROPERTY,
        ELECTRICAL_CONTACT_RESISTANCE_DIMS,
        2e-8,
    ));
    let card = card(properties.clone(), "original bond");
    let resolved = bundle(&card, &properties);
    let electrical =
        ResolvedElectricalInterface::try_from_resolved(&card, resolved.clone()).unwrap();
    let mechanical =
        ResolvedBilinearCohesive::try_from_resolved(&card, resolved.clone(), 0.02).unwrap();
    assert_eq!(mechanical.resolved(), electrical.resolved());
    assert_eq!(mechanical.resolved(), &resolved);
    assert_eq!(
        mechanical.surfaces(),
        &[card.surface_a().clone(), card.surface_b().clone()]
    );
    mechanical
        .require_surface_materials(&card.surface_a().material, &card.surface_b().material)
        .unwrap();
    assert!(
        mechanical
            .require_surface_materials(&card.surface_b().material, &card.surface_a().material)
            .is_err()
    );
    let failed = mechanical
        .law()
        .bonded_point()
        .trial([0.005, 0.0, 0.0], 1e-8)
        .unwrap();
    assert!(failed.point.failed());
    assert!((failed.energy.dissipated - 2000.0).abs() < 1e-9);
    assert_eq!(mechanical.resolved(), &resolved);
}

#[test]
fn fixed_temperature_data_change_damage_onset_and_failure_work() {
    let mut properties = properties();
    properties[1].value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(280.0, 1.2e6), (340.0, 0.6e6)],
        dims: Pressure::DIMS,
    };
    properties[1].interpolation = InterpolationPolicy::LinearInside;
    let card = card(properties, "unchanged preparation");
    let cold = resolve(&card);
    let warm = resolve_bilinear_cohesive(
        &card,
        &QueryPoint::new().with("T", 320.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        0.02,
    )
    .unwrap();
    assert_eq!(
        cold.resolved().card_identity(),
        warm.resolved().card_identity()
    );
    assert_ne!(cold.resolved().identity(), warm.resolved().identity());
    assert!(warm.law().onset_separation_m() < cold.law().onset_separation_m());
    assert!(warm.law().failure_separation_m() > cold.law().failure_separation_m());
    let cold_point = cold
        .law()
        .bonded_point()
        .trial([0.0009, 0.0, 0.0], 1e-8)
        .unwrap();
    let warm_point = warm
        .law()
        .bonded_point()
        .trial([0.0009, 0.0, 0.0], 1e-8)
        .unwrap();
    assert_eq!(cold_point.point.damage(), 0.0);
    assert!(warm_point.point.damage() > 0.0);
    for source in [&cold, &warm] {
        let failed = source
            .law()
            .bonded_point()
            .trial([0.01, 0.0, 0.0], 1e-8)
            .unwrap();
        assert!((failed.energy.dissipated - 2000.0).abs() < 1e-9);
    }
}

#[test]
fn wrong_cards_missing_units_and_invalid_envelopes_refuse() {
    let original = card(properties(), "original");
    let other = card(properties(), "different surface preparation history");
    let resolved = bundle(&original, &properties());
    assert!(ResolvedBilinearCohesive::try_from_resolved(&other, resolved.clone(), 0.02).is_err());
    for limit in [0.0, f64::NAN, 0.001] {
        assert!(
            ResolvedBilinearCohesive::try_from_resolved(&original, resolved.clone(), limit)
                .is_err()
        );
    }
    for kind in 0..4 {
        let mut p = properties();
        match kind {
            0 => {
                p.pop();
            }
            1 => {
                p[1] = claim(COHESIVE_STRENGTH_PROPERTY, Dims::NONE, 1e6);
            }
            2 => {
                p[0] = claim(
                    COHESIVE_STIFFNESS_PROPERTY,
                    interface_stiffness_per_area_dims(),
                    0.0,
                );
            }
            _ => {
                p[2] = claim(
                    COHESIVE_FRACTURE_ENERGY_PROPERTY,
                    FRACTURE_ENERGY_PER_AREA_DIMS,
                    100.0,
                );
            }
        }
        let changed = card(p, "invalid input fixture");
        assert!(
            resolve_bilinear_cohesive(
                &changed,
                &QueryPoint::new().with("T", 300.0).unwrap(),
                MaterialPropertySelection::SingleClaimOnly,
                0.02
            )
            .is_err()
        );
    }
}

#[test]
fn evolving_gap_and_rate_applicability_cannot_be_frozen_into_damage_law() {
    for axis in ["gap", "damage", "sliding_rate", "omega"] {
        let mut p = properties();
        p[0].validity = p[0].validity.clone().with(axis, 0.0, 1.0);
        let card = card(p, "bounded evolving condition");
        let point = QueryPoint::new()
            .with("T", 300.0)
            .unwrap()
            .with(axis, 0.1)
            .unwrap();
        assert!(
            resolve_bilinear_cohesive(
                &card,
                &point,
                MaterialPropertySelection::SingleClaimOnly,
                0.02
            )
            .is_err(),
            "{axis}"
        );
    }
    let mut p = properties();
    p[0].value = PropertyValue::Curve {
        abscissa: "gap".into(),
        abscissa_dims: Length::DIMS,
        knots: vec![(0.0, 1e9), (0.01, 2e9)],
        dims: interface_stiffness_per_area_dims(),
    };
    p[0].interpolation = InterpolationPolicy::LinearInside;
    p[0].validity = p[0].validity.clone().with("gap", 0.0, 0.01);
    let card = card(p, "separation-dependent coefficients");
    let point = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("gap", 0.0)
        .unwrap();
    assert!(
        resolve_bilinear_cohesive(
            &card,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
            0.02
        )
        .is_err()
    );
}

fn mixed_mode_claims() -> Vec<PropertyClaim> {
    let stiffness = interface_stiffness_per_area_dims();
    [
        (COHESIVE_NORMAL_STIFFNESS_PROPERTY, stiffness, 1e9),
        (COHESIVE_TANGENTIAL_STIFFNESS_PROPERTY, stiffness, 4e8),
        (COHESIVE_NORMAL_STRENGTH_PROPERTY, Pressure::DIMS, 2e5),
        (COHESIVE_TANGENTIAL_STRENGTH_PROPERTY, Pressure::DIMS, 1.6e5),
        (
            COHESIVE_MODE_I_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            100.0,
        ),
        (
            COHESIVE_MODE_II_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            400.0,
        ),
        (COHESIVE_BK_EXPONENT_PROPERTY, Dims::NONE, 2.0),
        (COHESIVE_COMPRESSION_STIFFNESS_PROPERTY, stiffness, 2e9),
    ]
    .into_iter()
    .map(|(name, dims, value)| claim(name, dims, value))
    .collect()
}

#[test]
fn source_mode_specific_strengths_and_toughness_reach_original_cohesive_point() {
    let card = card(mixed_mode_claims(), "synthetic mixed-mode initial bond");
    let source = resolve_mixed_mode_bilinear_cohesive(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        0.02,
    )
    .unwrap();
    assert_eq!(source.resolved().card_identity(), card.content_hash());
    assert_eq!(source.resolved().properties().len(), 8);
    assert!(source.law().mixed_mode_parameters().is_some());
    let opening = source
        .law()
        .bonded_point()
        .trial([0.002, 0.0, 0.0], 1e-8)
        .unwrap();
    let shear = source
        .law()
        .bonded_point()
        .trial([0.0, 0.006, 0.0], 1e-8)
        .unwrap();
    assert!(opening.point.failed() && shear.point.failed());
    assert!((opening.energy.dissipated - 100.0).abs() < 1e-10);
    assert!((shear.energy.dissipated - 400.0).abs() < 1e-10);
    assert!(
        source
            .require_surface_materials(
                &source.surfaces()[1].material,
                &source.surfaces()[0].material
            )
            .is_err()
    );
}

#[test]
fn mixed_source_does_not_infer_missing_modes_or_freeze_rate_conditions() {
    let common = card(properties(), "old common mode");
    let point = QueryPoint::new().with("T", 300.0).unwrap();
    assert!(
        resolve_mixed_mode_bilinear_cohesive(
            &common,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
            0.02
        )
        .is_err()
    );
    let mut claims = mixed_mode_claims();
    claims[6] = claim(COHESIVE_BK_EXPONENT_PROPERTY, Pressure::DIMS, 2.0);
    assert!(
        resolve_mixed_mode_bilinear_cohesive(
            &card(claims, "wrong exponent units"),
            &point,
            MaterialPropertySelection::SingleClaimOnly,
            0.02
        )
        .is_err()
    );
    let mut claims = mixed_mode_claims();
    claims[5].validity = claims[5].validity.clone().with("rate", 0.0, 1.0);
    let point = point.with("rate", 0.0).unwrap();
    assert!(
        resolve_mixed_mode_bilinear_cohesive(
            &card(claims, "unowned rate"),
            &point,
            MaterialPropertySelection::SingleClaimOnly,
            0.02
        )
        .is_err()
    );
}

#[test]
fn sourced_physical_viscosity_retains_the_original_interface_bundle() {
    let mut inputs=properties();
    inputs.push(claim(COHESIVE_NORMAL_VISCOSITY_PROPERTY,COHESIVE_VISCOSITY_DIMS,2e7));
    inputs.push(claim(COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY,COHESIVE_VISCOSITY_DIMS,3e7));
    let card=card(inputs.clone(),"explicit physical bond viscosity");
    let resolved=bundle(&card,&inputs);
    let owner=ResolvedBilinearCohesive::try_from_viscous_resolved(&card,resolved.clone(),
        CohesiveEnvelope::CommonMode,0.02,0.2).unwrap();
    assert_eq!(owner.resolved(),&resolved);
    assert_eq!(owner.surfaces(),&[card.surface_a().clone(),card.surface_b().clone()]);
    let step=owner.law().bonded_point().trial_over_time([0.0003,0.0,0.0],0.5,1e-8).unwrap();
    assert!((step.traction_pa[0]-312000.0).abs()<1e-8);
    assert!((step.viscous_dissipated_j_m2-3.6).abs()<1e-12);
    assert!(owner.law().bonded_point().trial([0.0003,0.0,0.0],1e-8).is_err());
    let queried=resolve_viscous_bilinear_cohesive(&card,&QueryPoint::new().with("T",300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,CohesiveEnvelope::CommonMode,0.02,0.2).unwrap();
    assert_eq!(queried.law(),owner.law());
}

#[test]
fn physical_viscosity_refuses_missing_bulk_units_and_unowned_rate_axes() {
    for (dims,rate_bound) in [(Dims([-1,1,-1,0,0,0]),false),(COHESIVE_VISCOSITY_DIMS,true)] {
        let mut inputs=properties();
        let mut normal=claim(COHESIVE_NORMAL_VISCOSITY_PROPERTY,dims,2e7);
        if rate_bound { normal.validity=normal.validity.with("separation_rate",0.0,1.0); }
        inputs.push(normal);inputs.push(claim(COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY,COHESIVE_VISCOSITY_DIMS,3e7));
        let card=card(inputs,"incompatible viscosity source");
        assert!(resolve_viscous_bilinear_cohesive(&card,&QueryPoint::new().with("T",300.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly,CohesiveEnvelope::CommonMode,0.02,0.2).is_err());
    }
    assert!(resolve_viscous_bilinear_cohesive(&card(properties(),"missing viscosity"),
        &QueryPoint::new().with("T",300.0).unwrap(),MaterialPropertySelection::SingleClaimOnly,
        CohesiveEnvelope::CommonMode,0.02,0.2).is_err());
}

#[test]
fn temperature_refresh_keeps_original_viscosity_claim_and_timed_physics() {
    let mut inputs = properties();
    let mut viscosity = claim(COHESIVE_NORMAL_VISCOSITY_PROPERTY, COHESIVE_VISCOSITY_DIMS, 1e7);
    viscosity.value = PropertyValue::Curve {
        abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
        knots: vec![(280.0, 1e7), (340.0, 3e7)], dims: COHESIVE_VISCOSITY_DIMS,
    };
    viscosity.interpolation = InterpolationPolicy::LinearInside;
    inputs.push(viscosity);
    inputs.push(claim(COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY, COHESIVE_VISCOSITY_DIMS, 3e7));
    let card = card(inputs, "manufactured temperature-dependent physical viscosity");
    let cold = resolve_viscous_bilinear_cohesive(&card,
        &QueryPoint::new().with("T", 280.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly, CohesiveEnvelope::CommonMode, 0.02, 0.2).unwrap();
    let hot = cold.at_temperature(&card, "T", 340.0).unwrap();
    cold.require_temperature_successor(&hot, "T").unwrap();
    assert_eq!(hot.law().viscosity().unwrap().normal_pa_s_m(), 3e7);
    assert_eq!(hot.law().viscosity().unwrap().maximum_rate_m_s(), 0.2);
    let loaded = cold.law().bonded_point()
        .trial_over_time([0.0003, 0.0, 0.0], 0.5, 1e-8).unwrap().equilibrium.point;
    let refreshed = loaded.transfer_temperature(hot.law(), 1e-8).unwrap();
    assert_eq!(refreshed.parameter_work_j_m2, 0.0);
    assert_eq!(refreshed.damage_dissipated_j_m2, 0.0);
    let continued = refreshed.point.trial_over_time([0.0006, 0.0, 0.0], 0.5, 1e-8).unwrap();
    assert!((continued.traction_pa[0] - 618000.0).abs() < 1e-8);
    assert!((continued.viscous_dissipated_j_m2 - 5.4).abs() < 1e-12);
    for (old, new) in cold.resolved().properties().iter().zip(hot.resolved().properties()) {
        assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
    }
}

#[test]
fn mixed_temperature_sources_keep_mode_claims_and_drive_the_original_timed_point() {
    for viscous in [false, true] {
        let mut inputs = mixed_mode_claims();
        for (property, hot) in inputs
            .iter_mut()
            .zip([5e8, 3e8, 8e4, 8e4, 50.0, 200.0, 1.5, 3e9])
        {
            let PropertyValue::Scalar { value, dims } = &property.value else {
                panic!("manufactured scalar fixture required");
            };
            property.value = PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(280.0, *value), (340.0, hot)],
                dims: *dims,
            };
            property.interpolation = InterpolationPolicy::LinearInside;
        }
        if viscous {
            inputs.push(claim(
                COHESIVE_NORMAL_VISCOSITY_PROPERTY,
                COHESIVE_VISCOSITY_DIMS,
                2e7,
            ));
            inputs.push(claim(
                COHESIVE_TANGENTIAL_VISCOSITY_PROPERTY,
                COHESIVE_VISCOSITY_DIMS,
                3e7,
            ));
        }
        let card = card(
            inputs,
            "manufactured thermal mixed-mode bond, not measured data",
        );
        let query = QueryPoint::new().with("T", 280.0).unwrap();
        let cold = if viscous {
            resolve_viscous_bilinear_cohesive(
                &card,
                &query,
                MaterialPropertySelection::SingleClaimOnly,
                CohesiveEnvelope::FrozenOnsetMixedMode,
                0.02,
                0.2,
            )
            .unwrap()
        } else {
            resolve_mixed_mode_bilinear_cohesive(
                &card,
                &query,
                MaterialPropertySelection::SingleClaimOnly,
                0.02,
            )
            .unwrap()
        };
        let loaded = cold
            .law()
            .bonded_point()
            .trial_over_time([0.0004, 0.0005, -0.0002], 0.5, 1e-8)
            .unwrap()
            .equilibrium
            .point;
        let original_mode = loaded.mixed_mode_history().unwrap();
        let hot = cold.at_temperature(&card, "T", 340.0).unwrap();
        cold.require_temperature_successor(&hot, "T").unwrap();
        assert_eq!(hot.surfaces(), cold.surfaces());
        assert_eq!(cold.resolved().query(), query);
        assert_eq!(hot.law().maximum_separation_m(), 0.02);
        assert_eq!(
            hot.law().viscosity().map(|v| v.maximum_rate_m_s()),
            viscous.then_some(0.2)
        );
        let parameters = hot.law().mixed_mode_parameters().unwrap();
        assert_eq!(parameters.normal_stiffness_pa_m, 5e8);
        assert_eq!(parameters.tangential_stiffness_pa_m, 3e8);
        assert_eq!(parameters.mode_i_fracture_energy_j_m2, 50.0);
        assert_eq!(parameters.mode_ii_fracture_energy_j_m2, 200.0);
        assert_eq!(parameters.bk_exponent, 1.5);
        for (old, new) in cold
            .resolved()
            .properties()
            .iter()
            .zip(hot.resolved().properties())
        {
            assert_eq!(old.requirement(), new.requirement());
            assert_eq!(old.answer().receipt.selected, new.answer().receipt.selected);
        }
        let transfer = loaded.transfer_temperature(hot.law(), 1e-8).unwrap();
        assert!(transfer.point.damage() >= loaded.damage());
        assert!(transfer.damage_dissipated_j_m2 > 0.0);
        let mode = transfer.point.mixed_mode_history().unwrap();
        assert_eq!(mode.onset_direction, original_mode.onset_direction);
        assert_eq!(
            mode.energy_shear_fraction,
            original_mode.energy_shear_fraction
        );
        assert!((mode.envelope_shear_fraction - 87.0 / 167.0).abs() < 1e-12);
        let end = [0.0005, -0.0002, 0.0003];
        let (continued, mean) = transfer
            .point
            .trial_over_time_with_mean_traction(end, 0.5, 1e-8)
            .unwrap();
        let work: f64 = (0..3)
            .map(|i| mean[i] * (end[i] - transfer.point.jump_m()[i]))
            .sum();
        assert!((work - continued.total_work_j_m2()).abs() < 1e-8);
        assert!(continued.equilibrium.point.dissipated_j_m2() >= loaded.dissipated_j_m2());
        assert_eq!(continued.viscous_dissipated_j_m2 > 0.0, viscous);
        assert!(cold.at_temperature(&card, "T", 341.0).is_err());
        assert!(
            cold.at_temperature(&card, "unowned_temperature", 300.0)
                .is_err()
        );
    }
}

#[test]
fn thermal_successor_cannot_switch_envelopes_even_with_the_same_source_bundle() {
    let mut inputs = mixed_mode_claims();
    inputs.extend(properties().into_iter().take(3));
    let card = card(
        inputs.clone(),
        "same complete original bundle, different model selection",
    );
    let resolved = bundle(&card, &inputs);
    let common =
        ResolvedBilinearCohesive::try_from_resolved(&card, resolved.clone(), 0.02).unwrap();
    let mixed =
        ResolvedBilinearCohesive::try_from_mixed_mode_resolved(&card, resolved, 0.02).unwrap();
    assert_eq!(common.resolved(), mixed.resolved());
    assert!(common.require_temperature_successor(&mixed, "T").is_err());
    assert!(mixed.require_temperature_successor(&common, "T").is_err());
    let hotter = mixed.at_temperature(&card, "T", 320.0).unwrap();
    mixed.require_temperature_successor(&hotter, "T").unwrap();
    assert!(hotter.law().mixed_mode_parameters().is_some());
    assert_eq!(hotter.resolved().properties().len(), 11);
}
