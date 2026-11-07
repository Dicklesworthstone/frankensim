//! One original interface fixture shared by the depth/history sliding tests.

use fs_contact::interface_binding::thermal_friction::{SourceThermalCoulomb, ThermalFrictionAxes};
use fs_contact::interface_binding::{
    ADHESION_ENERGY_DIMS, ADHESION_ENERGY_PROPERTY, BoundNormalInterfaceState,
    bind_dry_interface_system_card, bind_isotropic_elastic_interface, bind_normal_interface_state,
};
use fs_contact::normal_patch::FiniteGapPlasticLoadControls;
use fs_contact::normal_patch::plastic_source::wear::sliding::PlasticSlidingInput;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    InterfaceSystemCard, InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedInterfaceStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_interface_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};
use fs_tribo::{InputAuthority, TangentialSolveConfig};

pub fn bind(
    normal: &BoundNormalInterfaceState,
    keys: &[PropertyKey; 2],
    pressure_support: [f64; 2],
) -> (
    BoundNormalInterfaceState,
    ResolvedInterfaceStatePoint,
    SourceThermalCoulomb,
) {
    let old = normal.elastic().interface().card();
    let mut claims = old.claims().clone();
    let pressure_quantity =
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
    for (name, value) in [("mu_static", 0.6), ("mu_kinetic", 0.4)] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Dims::NONE),
                value: PropertyValue::Scalar {
                    value,
                    dims: Dims::NONE,
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 320.0)
                    .with_quantity(
                        "p",
                        pressure_quantity,
                        pressure_support[0],
                        pressure_support[1],
                    ),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic original spatial-hardness sliding fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        old.surface_a().clone(),
        old.surface_b().clone(),
        old.context().clone(),
        claims,
        old.models().to_vec(),
    )
    .unwrap();
    let dry = bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap();
    let elastic = bind_isotropic_elastic_interface(
        &dry,
        normal.elastic().surface_a_state(),
        normal.elastic().surface_b_state(),
    )
    .unwrap();
    let point = normal.interface_state().query();
    let select = MaterialPropertySelection::SingleClaimOnly;
    let adhesion = resolve_interface_state_point(
        &card,
        &point,
        &[ScalarPropertyRequirement::try_new(
            ADHESION_ENERGY_PROPERTY,
            ADHESION_ENERGY_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        select.clone(),
    )
    .unwrap();
    let bound = bind_normal_interface_state(&elastic, &adhesion).unwrap();
    let coefficients = resolve_interface_state_point(
        &card,
        &point,
        &keys.each_ref().map(|key| {
            ScalarPropertyRequirement::try_with_key(key, ScalarAdmissibility::NonNegative).unwrap()
        }),
        select.clone(),
    )
    .unwrap();
    let query = point
        .with_quantity(
            "p",
            pressure_quantity,
            0.5 * pressure_support[0] + 0.5 * pressure_support[1],
        )
        .unwrap();
    let requirement = |name| {
        ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative)
            .unwrap()
    };
    let friction = SourceThermalCoulomb::resolve(
        dry,
        query,
        requirement("mu_static"),
        requirement("mu_kinetic"),
        select,
        ThermalFrictionAxes {
            temperature: "T".into(),
            nominal_pressure: Some("p".into()),
        },
    )
    .unwrap();
    (bound, coefficients, friction)
}

pub fn controls() -> TangentialSolveConfig {
    TangentialSolveConfig {
        max_iterations: 64,
        maximum_force_n: 1.0,
        force_absolute_tolerance_n: 1.0e-14,
        force_relative_tolerance: 1.0e-12,
    }
}
pub fn normal_controls() -> FiniteGapPlasticLoadControls {
    FiniteGapPlasticLoadControls {
        minimum_approach_m: -1.0e-4,
        maximum_approach_m: 1.0e-4,
        absolute_force_tolerance_n: 1.0e-14,
        relative_force_tolerance: 1.0e-12,
        approach_tolerance_m: 1.0e-18,
        maximum_equilibrium_evaluations: 128,
    }
}
pub fn motion() -> Vec<PlasticSlidingInput> {
    vec![
        PlasticSlidingInput {
            released_increment_m: [0.0; 2],
            compliance_m_per_n: 0.0
        };
        25
    ]
}
