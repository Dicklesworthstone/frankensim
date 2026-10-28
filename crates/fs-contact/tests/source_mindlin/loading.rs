//! Original source geometry and typed conditions drive changing-load memory.
use fs_contact::interface_binding::mindlin::loading::{
    SourceHertzMindlinContact, SourceHertzMindlinState,
};
use fs_matdb::ClaimSet;
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

use super::*;

fn source(young: f64, mu: f64) -> SourceHertzMindlinContact {
    SourceHertzMindlinContact::from_fixed_source(
        resolve(normal(young, Some(mu), Dims::NONE, false)).unwrap(),
    )
    .unwrap()
}

fn shear(
    law: &SourceHertzMindlinContact,
    state: &SourceHertzMindlinState,
    du: f64,
) -> TangentialStepInput {
    let mut input = step(law.fixed_source(), du);
    input.normal_n = state.normal_force_n();
    input
}

#[test]
fn original_hertz_geometry_energy_and_authority_survive_load_changes() {
    let law = source(1e7, 0.6);
    let initial = law.initial_state();
    let original = law.fixed_source().normal_response();
    assert_eq!(initial.normal_energy_j(), original.reversible_energy_j);
    assert_eq!(initial.indentation_m(), original.approach_m);
    assert!(initial.normal_energy_j() > 0.0);
    assert_eq!(law.stored_energy_j(&initial), 0.0);
    let du = 0.6 * law.fixed_source().law().gross_slip_displacement_m();
    let loaded = law
        .advance(&initial, &shear(&law, &initial, du))
        .unwrap()
        .state;
    let grown = law.normal_step(&loaded, 8e-5, 1e-3).unwrap();
    assert_eq!(grown.dissipated_j, 0.0);
    assert_eq!(
        law.stored_energy_j(&grown.state),
        law.stored_energy_j(&loaded)
    );
    common::close(grown.state.endpoint_force_n(), loaded.endpoint_force_n());
    common::close(grown.state.normal_force_n(), 8.0 * initial.normal_force_n());
    common::close(
        grown.state.contact_radius_m(),
        2.0 * initial.contact_radius_m(),
    );
    assert_eq!(grown.normal_response.authority, original.authority);
    assert_eq!(grown.normal_response.uncertainty, original.uncertainty);
    assert_eq!(
        grown.normal_response.input_source_id,
        original.input_source_id
    );
    assert_eq!(
        grown.normal_work_j,
        grown.normal_response.reversible_energy_j - initial.normal_energy_j()
    );
    assert_eq!(
        grown.state.normal_energy_j(),
        grown.normal_response.reversible_energy_j
    );
    let smaller = law.normal_step(&grown.state, 5e-6, 1e-3).unwrap();
    assert!(smaller.dissipated_j > 0.0);
    assert!(smaller.normal_work_j < 0.0);
    common::close(smaller.tangential_stored_change_j, -smaller.dissipated_j);
    common::close(
        smaller.normal_response.pressure.resultant_n,
        smaller.state.normal_force_n(),
    );
}

#[test]
fn sourced_separation_releases_shear_once_and_recontact_starts_unstrained() {
    let law = source(1e7, 0.6);
    let initial = law.initial_state();
    let du = 0.7 * law.fixed_source().law().gross_slip_displacement_m();
    let loaded = law
        .advance(&initial, &shear(&law, &initial, du))
        .unwrap()
        .state;
    let released = law.normal_step(&loaded, 0.0, 1e-3).unwrap();
    assert_eq!(released.dissipated_j, law.stored_energy_j(&loaded));
    assert_eq!(released.normal_work_j, -initial.normal_energy_j());
    assert_eq!(released.normal_response.normal_force_n, 0.0);
    assert_eq!(released.state.normal_energy_j(), 0.0);
    assert_eq!(released.state.history_intervals(), 0);
    let free = law
        .advance(&released.state, &shear(&law, &released.state, 1e-3))
        .unwrap();
    assert_eq!(free.traction_n, [0.0; 2]);
    let same_open = law.normal_step(&free.state, 0.0, 1e-3).unwrap();
    assert_eq!(same_open.dissipated_j, 0.0);
    let new = law
        .normal_step(&same_open.state, initial.indentation_m(), 1e-3)
        .unwrap();
    assert_eq!(new.normal_work_j, initial.normal_energy_j());
    assert_eq!(new.state.endpoint_force_n(), 0.0);
    assert_eq!(law.stored_energy_j(&new.state), 0.0);
    let repeated = law
        .advance(&new.state, &shear(&law, &new.state, du))
        .unwrap();
    common::close(repeated.state.endpoint_force_n(), loaded.endpoint_force_n());
    assert_eq!(
        law.fixed_source().friction_source().identity(),
        source(1e7, 0.6).fixed_source().friction_source().identity()
    );
}

// Rebuild the same real source chain with an explicitly typed extra condition.
// Even a condition unused by a scalar value remains part of the original query.
fn extra_query(axis: &str, quantity: QuantitySpec, validity_bound: bool) -> SourceMindlinContact {
    let a = common::solid("Mindlin-A", 1000.0, 1e7);
    let b = common::solid("Mindlin-B", 2000.0, 2e7);
    let original = normal(1e7, Some(0.6), Dims::NONE, false);
    let card = original.model().normal_state().elastic().interface().card();
    let mut claims = ClaimSet::new();
    for (_, claim) in card.claims().claims_ordered() {
        let mut claim = claim.clone();
        if validity_bound && claim.key.name() == MU {
            claim.validity = claim.validity.with_quantity(axis, quantity, 0.0, 2e5);
        }
        claims.insert_claim(claim).unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        card.surface_a().clone(),
        card.surface_b().clone(),
        card.context().clone(),
        claims,
        card.models().to_vec(),
    )
    .unwrap();
    let point = common::query().with_quantity(axis, quantity, 1e5).unwrap();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let sa = resolve_isotropic_solid_state_point(&a, &point, selection.clone()).unwrap();
    let sb = resolve_isotropic_solid_state_point(&b, &point, selection.clone()).unwrap();
    let dry = bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap();
    let elastic = bind_isotropic_elastic_interface(&dry, &sa, &sb).unwrap();
    let interface = resolve_interface_state_point(
        &card,
        &point,
        &[ScalarPropertyRequirement::try_new(
            ADHESION_ENERGY_PROPERTY,
            ADHESION_ENERGY_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        selection,
    )
    .unwrap();
    let normal = bind_normal_interface_state(&elastic, &interface).unwrap();
    let model =
        bind_normal_contact_model(&normal, NormalContactModelSelection::SingleSupported).unwrap();
    resolve(SourcePointContact::new(model, original.geometry(), original.envelope()).unwrap())
        .unwrap()
}

#[test]
fn typed_pressure_and_unowned_load_conditions_cannot_be_frozen_during_contact_changes() {
    let pressure =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    // The name deliberately does not say "pressure": the original quantity
    // role drives admission. Constant-in-pressure applicability still refuses.
    for bounded in [false, true] {
        let original = extra_query("environment-coordinate", pressure, bounded);
        assert!(
            original.friction_source().properties()[0].is_static_at_fixed_external_conditions()
        );
        assert!(SourceHertzMindlinContact::from_fixed_source(original).is_err());
    }
    let force = QuantitySpec::dimensional(Dims([1, 1, -2, 0, 0, 0]));
    let original = extra_query("unknown-coordinate", force, false);
    assert!(SourceHertzMindlinContact::from_fixed_source(original).is_err());
    // The ordinary original temperature-only source is admitted unchanged.
    assert_eq!(source(1e7, 0.6).temperature_kelvin(), 293.15);
}

#[test]
fn original_rate_geometry_temperature_and_source_history_refusals_preserve_retry() {
    let law = source(1e7, 0.6);
    let initial = law.initial_state();
    let loaded = law
        .advance(
            &initial,
            &shear(
                &law,
                &initial,
                0.4 * law.fixed_source().law().gross_slip_displacement_m(),
            ),
        )
        .unwrap()
        .state;
    let saved = loaded.clone();
    for (indentation, dt) in [
        (2e-4, 1e-3),
        (8e-5, 1e-6),
        (1e-5, 0.0),
        (1e-5, 2e-3),
        (-1e-5, 1e-3),
        (f64::NAN, 1e-3),
    ] {
        assert!(law.normal_step(&loaded, indentation, dt).is_err());
    }
    let mut input = shear(&law, &loaded, 0.0);
    input.contact_temperature_k = Some(294.0);
    assert!(law.advance(&loaded, &input).is_err());
    input = shear(&law, &loaded, 0.0);
    input.nominal_pressure_pa = Some(1e5);
    assert!(law.advance(&loaded, &input).is_err());
    input = shear(&law, &loaded, 0.0);
    input.normal_n *= 0.9;
    assert!(law.advance(&loaded, &input).is_err());
    let foreign = SourceHertzMindlinContact::from_fixed_source(
        resolve(normal_with_history(
            1e7,
            Some(0.6),
            Dims::NONE,
            false,
            "another specimen history",
        ))
        .unwrap(),
    )
    .unwrap();
    assert_eq!(foreign.law().parameters(), law.law().parameters());
    assert!(foreign.normal_step(&loaded, 1e-5, 1e-3).is_err());
    assert!(
        foreign
            .advance(&loaded, &shear(&law, &loaded, 0.0))
            .is_err()
    );
    assert!(foreign.stored_energy_j(&loaded).is_nan());
    assert_eq!(loaded, saved);
    assert_eq!(
        law.normal_step(&loaded, 1e-5, 1e-3).unwrap(),
        law.normal_step(&saved, 1e-5, 1e-3).unwrap()
    );
}

#[test]
fn both_material_and_interface_substitutions_reach_the_evolving_response() {
    let soft = source(1e7, 0.6);
    let stiff = source(3e7, 0.6);
    let zero = source(1e7, 0.0);
    let soft_grown = soft
        .normal_step(&soft.initial_state(), 5e-5, 1e-3)
        .unwrap()
        .state;
    let stiff_grown = stiff
        .normal_step(&stiff.initial_state(), 5e-5, 1e-3)
        .unwrap()
        .state;
    let zero_grown = zero
        .normal_step(&zero.initial_state(), 5e-5, 1e-3)
        .unwrap()
        .state;
    common::close(
        soft_grown.contact_radius_m(),
        stiff_grown.contact_radius_m(),
    );
    assert!(stiff_grown.normal_force_n() > soft_grown.normal_force_n());
    let du = 1e-5;
    let a = soft
        .advance(&soft_grown, &shear(&soft, &soft_grown, du))
        .unwrap();
    let b = stiff
        .advance(&stiff_grown, &shear(&stiff, &stiff_grown, du))
        .unwrap();
    assert!(b.state.endpoint_force_n() > a.state.endpoint_force_n());
    let c = zero
        .advance(&zero_grown, &shear(&zero, &zero_grown, du))
        .unwrap();
    assert_eq!(c.traction_n, [0.0; 2]);
    assert_eq!(
        zero.normal_step(&c.state, 0.0, 1e-3).unwrap().dissipated_j,
        0.0
    );
}
