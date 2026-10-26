#![allow(dead_code)]
//! Shared manufactured material/geometry for actual induction and heat tests.

use fs_em::planar::{
    eddy::{EddyConductor, nonlinear::NonlinearEddyControls},
    nonlinear::{
        IsotropicMagneticMaterial, NonlinearMagneticControls, PlanarNonlinearField,
        PlanarNonlinearMagnetostatics,
    },
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::{
        nonlinear::{REVERSIBLE_FIELD_STRENGTH_PROPERTY, ReversibleMagneticMaterial},
        permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
    },
    state_point::{
        MaterialPropertySelection, ResolvedMaterialStatePoint, SPECIFIC_HEAT_CAPACITY_DIMS,
        ScalarAdmissibility, ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_DIMS,
        resolve_material_state_point,
    },
};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::CheckedCgConfig;

pub const DENSITY_KG_M3: f64 = 2.0;
pub const SPECIFIC_HEAT_J_KG_K: f64 = 3.0;

pub fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
        .with_quantity("B", QuantitySpec::dimensional(FLUX_DENSITY_DIMS), 0.0)
        .unwrap()
}
pub fn source(knots: &[(f64, f64)], sigma: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let mut insert = |name: &str, dims: Dims, value: PropertyValue, validity| {
        let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured nonlinear induction".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    let domain = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        200.0,
        500.0,
    );
    insert(
        REVERSIBLE_FIELD_STRENGTH_PROPERTY,
        FIELD_STRENGTH_DIMS,
        PropertyValue::Curve {
            abscissa: "B".into(),
            abscissa_dims: FLUX_DENSITY_DIMS,
            knots: knots.to_vec(),
            dims: FIELD_STRENGTH_DIMS,
        },
        domain
            .clone()
            .with_quantity("B", QuantitySpec::dimensional(FLUX_DENSITY_DIMS), 0.0, 4.0),
    );
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        let values = match basis {
            ElectricalPropertyBasis::Resistivity => vec![
                (200.0, 1.0 / sigma),
                (300.0, 1.0 / sigma),
                (500.0, 2.0 / sigma),
            ],
            ElectricalPropertyBasis::Conductivity => {
                vec![(200.0, sigma), (300.0, sigma), (500.0, 0.5 * sigma)]
            }
        };
        insert(
            basis.property_name(),
            basis.dimensions(),
            PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: values,
                dims: basis.dimensions(),
            },
            domain.clone(),
        );
    }
    for (name, dims, value) in [
        ("density", Density::DIMS, DENSITY_KG_M3),
        (
            "specific_heat_capacity",
            SPECIFIC_HEAT_CAPACITY_DIMS,
            SPECIFIC_HEAT_J_KG_K,
        ),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 1.0),
    ] {
        insert(
            name,
            dims,
            PropertyValue::Scalar { value, dims },
            domain.clone(),
        );
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic conductor".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
pub fn controls() -> NonlinearEddyControls {
    NonlinearEddyControls {
        nonlinear: NonlinearMagneticControls {
            absolute_residual_tolerance_a_m: 1e-11,
            relative_residual_tolerance: 1e-10,
            max_newton_iterations: 30,
            max_backtracks: 30,
            linear: CheckedCgConfig {
                tolerance: 1e-12,
                max_iterations: 1000,
                max_corrections: 2,
            },
        },
        energy_tolerance_j: 1e-9,
        current_balance_tolerance_a: 1e-8,
    }
}
pub fn mesh(depth: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "nonlinear-induction-square",
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap()
}
pub fn law(source: &MaterialCard) -> ReversibleMagneticMaterial {
    ReversibleMagneticMaterial::resolve(
        source,
        &point(300.0),
        "B",
        4.0,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
pub fn problem(source: &MaterialCard, mesh: TriComplex2) -> PlanarNonlinearMagnetostatics {
    let laws = vec![IsotropicMagneticMaterial::Reversible(law(source)); mesh.faces().len()];
    PlanarNonlinearMagnetostatics::new(
        mesh,
        laws,
        WeightedAssemblyLimits {
            max_cells: 1000,
            max_dofs: 3000,
            max_triplets: 10000,
        },
        &mut || false,
    )
    .unwrap()
}
pub fn bundle(
    source: &MaterialCard,
    t: f64,
    basis: ElectricalPropertyBasis,
) -> ResolvedMaterialStatePoint {
    resolve_material_state_point(
        source,
        &point(t),
        &[
            basis.requirement().unwrap(),
            ScalarPropertyRequirement::try_new(
                REVERSIBLE_FIELD_STRENGTH_PROPERTY,
                FIELD_STRENGTH_DIMS,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                "density",
                Density::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                "specific_heat_capacity",
                SPECIFIC_HEAT_CAPACITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                "thermal_conductivity",
                THERMAL_CONDUCTIVITY_DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
        ],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
pub fn conductor(
    source: &MaterialCard,
    t: f64,
    basis: ElectricalPropertyBasis,
) -> ResolvedConductorState {
    ResolvedConductorState::try_from_material_with_basis(bundle(source, t, basis), basis).unwrap()
}
pub fn conductors(
    source: &MaterialCard,
    count: usize,
    t: f64,
    basis: ElectricalPropertyBasis,
) -> Vec<EddyConductor> {
    let state = conductor(source, t, basis);
    (0..count)
        .map(|triangle| EddyConductor {
            triangle,
            region: 7,
            conductor: state.clone(),
        })
        .collect()
}
pub fn boundary(offset: f64) -> Vec<(usize, f64)> {
    (0..4).map(|v| (v, offset)).collect()
}
pub fn initial(problem: &PlanarNonlinearMagnetostatics) -> PlanarNonlinearField<'_> {
    problem
        .solve(
            &[12.0; 4],
            &boundary(0.0),
            &[0.0; 5],
            controls().nonlinear,
            &mut || false,
        )
        .unwrap()
}
pub fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-8 * (1.0 + b.abs()), "{a} != {b}");
}
