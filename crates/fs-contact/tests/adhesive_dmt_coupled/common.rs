//! Manufactured ordered source data through the original material resolvers.
#![allow(dead_code)]
use fs_contact::interface_binding::adhesive_shear::*;
use fs_contact::interface_binding::*;
use fs_contact::normal_patch::dmt::{coupled::*, *};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_interface_state_point, resolve_isotropic_solid_state_point,
};
use fs_qty::{Density, Dims, Pressure};
use fs_tribo::InputAuthority;

pub const Z: f64 = 1e-9;
fn claim(name: &str, value: f64, dims: Dims) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: format!("manufactured joint DMT: {name}"),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn solid(name: &str, young: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        ("density", 1000.0, Density::DIMS),
        ("young_modulus", young, Pressure::DIMS),
        ("poisson_ratio", 0.25, Dims::NONE),
        ("yield_stress", 1e10, Pressure::DIMS),
    ] {
        claims.insert_claim(claim(name, value, dims)).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "manufactured test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
pub fn sources(
    adhesion: f64,
    static_strength: f64,
    kinetic_strength: f64,
) -> (SourceDmtSphereContact, SourceAdhesiveShear) {
    let a = solid("DMT A", 1e11);
    let b = solid("DMT B", 2e11);
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        (ADHESION_ENERGY_PROPERTY, adhesion, ADHESION_ENERGY_DIMS),
        (
            STATIC_JUNCTION_SHEAR_PROPERTY,
            static_strength,
            Pressure::DIMS,
        ),
        (
            KINETIC_JUNCTION_SHEAR_PROPERTY,
            kinetic_strength,
            Pressure::DIMS,
        ),
    ] {
        claims.insert_claim(claim(name, value, dims)).unwrap();
    }
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: a.id().clone(),
            texture_frame: "A local plane".into(),
        },
        SurfaceSpec {
            material: b.id().clone(),
            texture_frame: "B local plane".into(),
        },
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "manufactured".into(),
            history: "initial".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new().with("T", 293.15).unwrap();
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
        selection.clone(),
    )
    .unwrap();
    let normal = bind_normal_interface_state(&elastic, &interface).unwrap();
    let shear = SourceAdhesiveShear::resolve(&normal, selection).unwrap();
    let normal = SourceDmtSphereContact::new(
        normal,
        1e-6,
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
    .unwrap();
    (normal, shear)
}
pub fn config() -> DmtCoupledConfig {
    DmtCoupledConfig {
        maximum_iterations: 64,
        maximum_evaluations: 800,
        maximum_backtracks: 24,
        maximum_step_s: 1e-4,
        minimum_indentation_m: -5.0 * Z,
        maximum_coordinate_increment_m: 3.0 * Z,
        maximum_slip_speed_m_s: 1.0,
        maximum_force_n: 1e-3,
        coordinate_absolute_tolerance_m: 1e-22,
        coordinate_relative_tolerance: 1e-10,
        force_absolute_tolerance_n: 1e-18,
        force_relative_tolerance: 1e-8,
    }
}
pub fn compliance() -> [[f64; 3]; 3] {
    [
        [2e-5, 0.3e-5, -0.2e-5],
        [0.3e-5, 3e-5, 0.4e-5],
        [-0.2e-5, 0.4e-5, 4e-5],
    ]
}
pub fn mul(a: [[f64; 3]; 3], x: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| (0..3).map(|j| a[i][j] * x[j]).sum())
}
pub fn free(s: [[f64; 3]; 3], q: [f64; 3], f: [f64; 3]) -> [f64; 3] {
    let v = mul(s, f);
    std::array::from_fn(|i| q[i] + v[i])
}
pub fn close(a: f64, b: f64, absolute: f64, relative: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
