//! Synthetic interface evidence through the real material query owner.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedInterfaceStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_interface_state_point,
};
use fs_mbd_material::creeping_flow::navier_slip::NAVIER_SLIP_LENGTH_PROPERTY;
use fs_qty::{Dims, Length};

pub fn constant(
    solid: &MaterialStateId,
    slip_length_m: f64,
    maximum_temperature_k: f64,
) -> ResolvedInterfaceStatePoint {
    resolve(
        solid,
        slip_length_m,
        Length::DIMS,
        ValidityDomain::unconstrained().with("T", 200.0, maximum_temperature_k),
        &QueryPoint::new().with("T", 300.0).unwrap(),
    )
}

pub fn resolve(
    solid: &MaterialStateId,
    slip_length_m: f64,
    dims: Dims,
    validity: ValidityDomain,
    point: &QueryPoint,
) -> ResolvedInterfaceStatePoint {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(NAVIER_SLIP_LENGTH_PROPERTY, dims),
            value: PropertyValue::Scalar {
                value: slip_length_m,
                dims,
            },
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic solid/fluid Navier interface".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: solid.clone(),
            texture_frame: "homogeneous isotropic sphere".into(),
        },
        SurfaceSpec {
            material: MaterialStateId {
                chemistry: "synthetic fluid".into(),
                phase: "liquid".into(),
                process: "implementation fixture".into(),
                revision: 0,
            },
            texture_frame: "uniform Newtonian reservoir".into(),
        },
        SystemContext {
            medium: "synthetic fluid".into(),
            third_body: None,
            environment: "fixed pressure and composition".into(),
            history: "homogeneous clean interface".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let requirement = ScalarPropertyRequirement::try_new(
        NAVIER_SLIP_LENGTH_PROPERTY,
        dims,
        ScalarAdmissibility::Finite,
    )
    .unwrap();
    resolve_interface_state_point(
        &card,
        point,
        &[requirement],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
