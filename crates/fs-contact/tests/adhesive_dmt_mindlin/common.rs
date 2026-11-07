//! Original synthetic bulk/interface sources shared with native coupled tests.
#![allow(dead_code)]
#[path = "../source_point/common.rs"]
mod point;
use fs_contact::interface_binding::*;
use fs_contact::normal_patch::dmt::{mindlin::*, *};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    InterfaceSystemCard, InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_interface_state_point, resolve_isotropic_solid_state_point,
};
use fs_qty::Dims;
use fs_tribo::InputAuthority;

pub const R: f64 = 1e-6;
pub const Z: f64 = 1e-9;
pub const MU: &str = "local-coulomb-coefficient";

pub fn normal(adhesion: f64, coefficient: Option<f64>, history: &str) -> SourceDmtSphereContact {
    normal_with_dimensions(adhesion, coefficient, history, Dims::NONE)
}

pub fn normal_with_dimensions(
    adhesion: f64,
    coefficient: Option<f64>,
    history: &str,
    dimensions: Dims,
) -> SourceDmtSphereContact {
    let a = point::solid("DMT Mindlin A", 1000.0, 1e11);
    let b = point::solid("DMT Mindlin B", 2000.0, 2e11);
    let original = point::bind(&a, &b, None, adhesion);
    let card = original.normal_state().elastic().interface().card();
    let mut claims = card.claims().clone();
    if let Some(value) = coefficient {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(MU, dimensions),
                value: PropertyValue::Scalar {
                    value,
                    dims: dimensions,
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic DMT local Coulomb data".into(),
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
    let query = point::query();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let sa = resolve_isotropic_solid_state_point(&a, &query, selection.clone()).unwrap();
    let sb = resolve_isotropic_solid_state_point(&b, &query, selection.clone()).unwrap();
    let dry = bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap();
    let elastic = bind_isotropic_elastic_interface(&dry, &sa, &sb).unwrap();
    let interface = resolve_interface_state_point(
        &card,
        &query,
        &[ScalarPropertyRequirement::try_new(
            ADHESION_ENERGY_PROPERTY,
            ADHESION_ENERGY_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        selection,
    )
    .unwrap();
    SourceDmtSphereContact::new(
        bind_normal_interface_state(&elastic, &interface).unwrap(),
        R,
        Z,
        DmtEnvelope {
            maximum_indentation_m: Z,
            half_space_depth_m: 1e-3,
            layer_thickness_m: 1e-3,
            maximum_interaction_radius_to_radius: 0.1,
            maximum_interaction_radius_to_depth: 0.1,
            maximum_interaction_radius_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.3,
            maximum_adhesive_stress_to_yield: 0.1,
            maximum_interaction_range_to_radius: 0.01,
            maximum_tabor_parameter: 0.1,
        },
    )
    .unwrap()
}

pub fn requirement() -> ScalarPropertyRequirement {
    ScalarPropertyRequirement::try_new(MU, Dims::NONE, ScalarAdmissibility::NonNegative).unwrap()
}
pub fn source(
    adhesion: f64,
    coefficient: Option<f64>,
    history: &str,
) -> Result<SourceDmtMindlinContact, DmtMindlinError> {
    SourceDmtMindlinContact::resolve(
        normal(adhesion, coefficient, history),
        requirement(),
        MaterialPropertySelection::SingleClaimOnly,
        64,
    )
}
