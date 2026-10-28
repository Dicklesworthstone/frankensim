//! G1/G3: the original bulk/interface/Hertz sources drive real partial slip.
#[path = "source_point/common.rs"]
mod common;

#[path = "source_mindlin/loading.rs"]
mod loading;

#[path = "source_mindlin/sliding.rs"]
mod sliding;

use fs_contact::interface_binding::mindlin::{SourceMindlinContact, SourceMindlinError};
use fs_contact::interface_binding::*;
use fs_contact::normal_patch::source::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    InterfaceSystemCard, InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_interface_state_point, resolve_isotropic_solid_state_point,
};
use fs_qty::{Dims, Pressure};
use fs_tribo::{InputAuthority, TangentialLaw, TangentialSolveConfig, TangentialStepInput};

const MU: &str = "local-coulomb-coefficient";

fn normal(young: f64, mu: Option<f64>, dims: Dims, rate_axis: bool) -> SourcePointContact {
    normal_with_history(young, mu, dims, rate_axis, "initial")
}

fn normal_with_history(
    young: f64,
    mu: Option<f64>,
    dims: Dims,
    rate_axis: bool,
    history: &str,
) -> SourcePointContact {
    let a = common::solid("Mindlin-A", 1000.0, young);
    let b = common::solid("Mindlin-B", 2000.0, 2e7);
    let original = common::bind(&a, &b, None, 0.0);
    let card = original.normal_state().elastic().interface().card();
    let mut claims = card.claims().clone();
    let mut point = common::query();
    if let Some(mu) = mu {
        let mut validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
        if rate_axis {
            validity = validity.with("slip_speed", 0.0, 1.0);
            point = point.with("slip_speed", 0.0).unwrap();
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(MU, dims),
                value: PropertyValue::Scalar { value: mu, dims },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic ordered-interface coefficient".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let mut context = card.context().clone();
    context.history = history.into();
    let card = InterfaceSystemCard::assemble(
        card.surface_a().clone(),
        card.surface_b().clone(),
        context,
        claims,
        card.models().to_vec(),
    )
    .unwrap();
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
    let state = bind_normal_interface_state(&elastic, &interface).unwrap();
    let model =
        bind_normal_contact_model(&state, NormalContactModelSelection::SingleSupported).unwrap();
    SourcePointContact::new(
        model,
        PointContactGeometry::SpherePlane {
            effective_radius_m: 0.1,
        },
        common::envelope(),
    )
    .unwrap()
}

fn resolve(normal: SourcePointContact) -> Result<SourceMindlinContact, SourceMindlinError> {
    SourceMindlinContact::resolve(
        normal,
        2e-5,
        ScalarPropertyRequirement::try_new(MU, Dims::NONE, ScalarAdmissibility::NonNegative)
            .unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        64,
    )
}

fn step(law: &SourceMindlinContact, slip: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: law.normal_response().normal_force_n,
        slip_increment_m: [slip, 0.0],
        dt_s: 1e-4,
        contact_temperature_k: Some(law.temperature_kelvin()),
        nominal_pressure_pa: None,
        compliance_m_per_n: 0.0,
        solve: TangentialSolveConfig {
            max_iterations: 96,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-12,
            force_relative_tolerance: 1e-12,
        },
    }
}

#[test]
fn original_bulk_interface_and_hertz_patch_determine_the_tangential_response() {
    let original = normal(1e7, Some(0.6), Dims::NONE, false);
    let source = resolve(original.clone()).unwrap();
    let radius = (0.1_f64 * 2e-5).sqrt();
    let expected_stiffness = 8.0 * radius / (1.75 / (1e7 / 2.5) + 1.75 / (2e7 / 2.5));
    common::close(source.law().initial_stiffness_n_per_m(), expected_stiffness);
    common::close(source.normal_response().patch_radius_m, radius);
    assert_eq!(source.normal_source(), &original);
    assert_eq!(
        source.friction_source().card_identity(),
        original
            .model()
            .normal_state()
            .interface_state()
            .card_identity()
    );
    assert_eq!(
        source.friction_source().query_point(),
        original
            .model()
            .normal_state()
            .interface_state()
            .query_point()
    );
    assert_eq!(
        source.normal_response().authority,
        InputAuthority::SyntheticFixture
    );
    assert_eq!(source.friction_source().properties().len(), 1);
    let limit = source.law().gross_slip_displacement_m();
    let capacity = source.law().gross_slip_capacity_n();
    let result = source
        .advance(&source.initial_state(), &step(&source, 0.2 * limit))
        .unwrap();
    common::close(
        result.state.endpoint_force_n(),
        capacity * (1.0 - 0.8_f64.powf(1.5)),
    );
    assert!(result.dissipated_j > 0.0);
    common::close(
        source.stored_energy_j(&result.state) + result.dissipated_j,
        result.traction_n[0] * result.slip_increment_m[0],
    );
}

#[test]
fn material_and_ordered_interface_changes_reach_distinct_physical_quantities() {
    let soft = resolve(normal(1e7, Some(0.6), Dims::NONE, false)).unwrap();
    let stiff = resolve(normal(3e7, Some(0.6), Dims::NONE, false)).unwrap();
    let slippery = resolve(normal(1e7, Some(0.2), Dims::NONE, false)).unwrap();
    assert!(stiff.law().initial_stiffness_n_per_m() > soft.law().initial_stiffness_n_per_m());
    assert!(stiff.normal_response().normal_force_n > soft.normal_response().normal_force_n);
    common::close(
        slippery.law().initial_stiffness_n_per_m(),
        soft.law().initial_stiffness_n_per_m(),
    );
    common::close(
        slippery.normal_response().normal_force_n,
        soft.normal_response().normal_force_n,
    );
    common::close(
        slippery.law().gross_slip_capacity_n(),
        soft.law().gross_slip_capacity_n() / 3.0,
    );
    assert_ne!(
        slippery.friction_source().card_identity(),
        soft.friction_source().card_identity()
    );
    assert!(
        stiff
            .advance(&soft.initial_state(), &step(&stiff, 1e-6))
            .is_err()
    );
}

#[test]
fn equal_coefficients_do_not_transfer_history_between_distinct_interface_sources() {
    let original = resolve(normal(1e7, Some(0.6), Dims::NONE, false)).unwrap();
    let other = resolve(normal_with_history(
        1e7,
        Some(0.6),
        Dims::NONE,
        false,
        "another specimen history",
    ))
    .unwrap();
    assert_eq!(original.law(), other.law());
    assert_ne!(
        original.friction_source().card_identity(),
        other.friction_source().card_identity(),
    );
    let loaded = original
        .advance(
            &original.initial_state(),
            &step(&original, 0.3 * original.law().gross_slip_displacement_m()),
        )
        .unwrap();
    let saved = loaded.state.clone();
    let release = -0.1 * original.law().gross_slip_displacement_m();
    let expected = original.advance(&saved, &step(&original, release)).unwrap();
    assert!(other.stored_energy_j(&saved).is_nan());
    assert!(other.advance(&saved, &step(&other, release)).is_err());
    assert_eq!(loaded.state, saved);
    assert_eq!(
        original
            .advance(&loaded.state, &step(&original, release))
            .unwrap(),
        expected,
    );
    assert_eq!(
        saved.endpoint_force_n(),
        saved.constitutive_state().endpoint_force_n(),
    );
    assert_eq!(
        saved.history_intervals(),
        saved.constitutive_state().history_intervals(),
    );
}

#[test]
fn missing_wrong_units_or_rate_bounded_friction_does_not_become_a_default() {
    for source in [
        normal(1e7, None, Dims::NONE, false),
        normal(1e7, Some(0.6), Pressure::DIMS, false),
        normal(1e7, Some(-0.1), Dims::NONE, false),
        normal(1e7, Some(0.6), Dims::NONE, true),
    ] {
        assert!(resolve(source).is_err());
    }
    let absent = common::source(1e7, None);
    assert!(resolve(absent).is_err());
}

#[test]
fn ellipse_or_rate_augmented_normal_law_cannot_be_relabelled_as_circular_hertz() {
    let original = normal(1e7, Some(0.6), Dims::NONE, false);
    let ellipse = SourcePointContact::new(
        original.model().clone(),
        PointContactGeometry::EllipticParaboloid {
            maximum_curvature_m_inverse: 20.0,
            minimum_curvature_m_inverse: 5.0,
        },
        common::envelope(),
    )
    .unwrap();
    assert!(matches!(
        resolve(ellipse),
        Err(SourceMindlinError::Invalid(_))
    ));
    assert!(matches!(
        resolve(common::source(1e7, Some(0.4))),
        Err(SourceMindlinError::Invalid(_))
    ));
    assert!(
        SourceMindlinContact::resolve(
            original,
            0.0,
            ScalarPropertyRequirement::try_new(MU, Dims::NONE, ScalarAdmissibility::NonNegative)
                .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            64
        )
        .is_err()
    );
}

#[test]
fn frozen_temperature_load_direction_and_normal_envelope_refuse_without_losing_history() {
    let law = resolve(normal(1e7, Some(0.6), Dims::NONE, false)).unwrap();
    let initial = law.initial_state();
    let loaded = law
        .advance(
            &initial,
            &step(&law, 0.2 * law.law().gross_slip_displacement_m()),
        )
        .unwrap();
    let saved = loaded.state.clone();
    let input = step(&law, -0.1 * law.law().gross_slip_displacement_m());
    let expected = law.advance(&saved, &input).unwrap();
    for bad in [
        TangentialStepInput {
            contact_temperature_k: None,
            ..input
        },
        TangentialStepInput {
            contact_temperature_k: Some(300.0),
            ..input
        },
        TangentialStepInput {
            normal_n: input.normal_n * 1.01,
            ..input
        },
        TangentialStepInput {
            nominal_pressure_pa: Some(1e5),
            ..input
        },
        TangentialStepInput {
            slip_increment_m: [0.0, 1e-6],
            ..input
        },
        TangentialStepInput {
            dt_s: 0.002,
            ..input
        },
    ] {
        assert!(law.advance(&loaded.state, &bad).is_err());
        assert_eq!(loaded.state, saved);
    }
    assert_eq!(law.advance(&loaded.state, &input).unwrap(), expected);
}

#[test]
fn explicitly_sourced_zero_friction_remains_the_frictionless_limit() {
    let law = resolve(normal(1e7, Some(0.0), Dims::NONE, false)).unwrap();
    let result = law
        .advance(&law.initial_state(), &step(&law, 1e-4))
        .unwrap();
    assert_eq!(result.traction_n, [0.0; 2]);
    assert_eq!(result.dissipated_j, 0.0);
    assert_eq!(law.stored_energy_j(&result.state), 0.0);
    assert_eq!(result.slip_increment_m, [1e-4, 0.0]);
}
