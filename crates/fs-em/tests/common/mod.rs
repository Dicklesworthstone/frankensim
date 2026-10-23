#![allow(dead_code)]
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    field_storage::{MAGNETIC_PERMEABILITY_DIMS, MAGNETIC_PERMEABILITY_PROPERTY},
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
        ScalarPropertyRequirement, resolve_material_state_point,
    },
};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::CheckedCgConfig;

pub fn source(value: PropertyValue, validity: ValidityDomain) -> MaterialCard {
    let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
        InterpolationPolicy::ConstantWithinValidity
    } else {
        InterpolationPolicy::LinearInside
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(MAGNETIC_PERMEABILITY_PROPERTY, MAGNETIC_PERMEABILITY_DIMS),
            value,
            validity,
            interpolation,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic magnetic regression".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
pub fn resolve(card: &MaterialCard, point: &QueryPoint) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(
        card,
        point,
        &[ScalarPropertyRequirement::try_new(
            MAGNETIC_PERMEABILITY_PROPERTY,
            MAGNETIC_PERMEABILITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
pub fn material(mu: f64) -> ResolvedMaterialStatePoint {
    resolve(
        &source(
            PropertyValue::Scalar {
                value: mu,
                dims: MAGNETIC_PERMEABILITY_DIMS,
            },
            ValidityDomain::unconstrained(),
        ),
        &QueryPoint::new(),
    )
}
pub fn limits() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 10000,
        max_dofs: 20000,
        max_triplets: 90000,
    }
}
pub fn config() -> CheckedCgConfig {
    CheckedCgConfig {
        tolerance: 1e-11,
        max_iterations: 1000,
        max_corrections: 2,
    }
}
pub fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-9 * b.abs().max(1.0), "{a} != {b}");
}
pub fn layered(width: f64, height: f64, depth: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "magnetic-layered-regression",
        vec![
            [0.0, 0.0],
            [width, 0.0],
            [2.0 * width, 0.0],
            [0.0, height],
            [width, height],
            [2.0 * width, height],
        ],
        vec![[0, 1, 4], [0, 4, 3], [1, 2, 5], [1, 5, 4]],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap()
}
pub fn ends(value: f64) -> Vec<(usize, f64)> {
    vec![(0, 0.0), (3, 0.0), (2, value), (5, value)]
}
