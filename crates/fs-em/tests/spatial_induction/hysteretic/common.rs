//! Manufactured source cards and independent scalar play balances.
//! These coefficients qualify equations, not an empirical magnetic material.

use fs_em::spatial::vector_potential::induction::{
    SpatialInductionControls, SpatialInductionLimits,
    hysteretic::{
        SpatialHystereticInduction, SpatialHystereticInductionControls, SpatialHystereticProblem,
    },
    nonlinear::SpatialNonlinearInductionControls,
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
        MAGNETIC_PERMEABILITY_DIMS,
        hysteresis::{
            source::{
                HysteresisBranchSources, HysteresisSources, MAGNETIC_VISCOSITY_DIMS,
                ResolvedScalarHysteresis, ResolvedVectorHysteresis,
            },
            vector::VectorHysteresisControls,
        },
        permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
    },
    state_point::{
        MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, resolve_material_state_point,
    },
};
use fs_qty::{Density, QuantitySpec, Temperature};
use std::{collections::BTreeSet, sync::Arc};

pub const MU: f64 = 0.01;
pub const CHI: f64 = 0.2;
pub const PIN: f64 = 2.0;
pub const CP: f64 = 0.1;
pub const DT: f64 = 0.1;
// Exact curls of the body-diagonal Whitney edge on the six unit-cube tets.
pub const DIRECTIONS: [[f64; 3]; 6] = [
    [0.0, 1.0, 0.0],
    [-1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, -1.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 0.0, -1.0],
];

pub fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
}
fn domain(upper: f64) -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        200.0,
        upper,
    )
}
fn roles() -> HysteresisSources {
    HysteresisSources {
        background_permeability: PropertyKey::new("background_mu", MAGNETIC_PERMEABILITY_DIMS),
        field_limit: PropertyKey::new("field_limit", FIELD_STRENGTH_DIMS),
        branches: vec![HysteresisBranchSources {
            susceptibility: PropertyKey::new("branch_slope", MAGNETIC_PERMEABILITY_DIMS),
            pinning_field: PropertyKey::new("branch_pinning", FIELD_STRENGTH_DIMS),
            saturation_polarization: PropertyKey::new("branch_saturation", FLUX_DENSITY_DIMS),
        }],
    }
}
pub fn card(fit_upper: f64) -> Arc<MaterialCard> {
    card_with_heat_capacity(fit_upper, CP)
}
pub fn card_with_heat_capacity(fit_upper: f64, heat_capacity: f64) -> Arc<MaterialCard> {
    card_with_parameters(fit_upper, heat_capacity, None)
}
pub fn card_with_viscosity(fit_upper: f64, eta: f64) -> Arc<MaterialCard> {
    card_with_parameters(fit_upper, CP, Some(eta))
}
fn viscosity_key() -> PropertyKey {
    PropertyKey::new("intrinsic_branch_viscosity", MAGNETIC_VISCOSITY_DIMS)
}
fn card_with_parameters(
    fit_upper: f64,
    heat_capacity: f64,
    viscosity: Option<f64>,
) -> Arc<MaterialCard> {
    let roles = roles();
    let branch = &roles.branches[0];
    let mut claims = ClaimSet::new();
    let mut insert = |key: PropertyKey, value, validity, interpolation| {
        claims
            .insert_claim(PropertyClaim {
                key,
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured spatial vector induction and finite heat".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    for (key, value, upper) in [
        (roles.background_permeability.clone(), MU, fit_upper),
        (roles.field_limit.clone(), 30.0, fit_upper),
        (branch.susceptibility.clone(), CHI, fit_upper),
        (branch.pinning_field.clone(), PIN, fit_upper),
        (branch.saturation_polarization.clone(), 1.0, fit_upper),
        (PropertyKey::new("density", Density::DIMS), 1.0, 500.0),
        (
            PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            heat_capacity,
            500.0,
        ),
    ] {
        let dims = key.dims();
        insert(
            key,
            PropertyValue::Scalar { value, dims },
            domain(upper),
            InterpolationPolicy::ConstantWithinValidity,
        );
    }
    if let Some(value) = viscosity {
        insert(
            viscosity_key(),
            PropertyValue::Scalar {
                value,
                dims: MAGNETIC_VISCOSITY_DIMS,
            },
            domain(fit_upper),
            InterpolationPolicy::ConstantWithinValidity,
        );
    }
    let basis = ElectricalPropertyBasis::Resistivity;
    insert(
        PropertyKey::new(basis.property_name(), basis.dimensions()),
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            dims: basis.dimensions(),
            knots: [200.0, 300.0, 500.0].map(|t| (t, resistivity(t))).to_vec(),
        },
        domain(500.0),
        InterpolationPolicy::LinearInside,
    );
    Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic magnetic conductor".into(),
                phase: "solid".into(),
                process: "manufactured spatial history".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    )
}
pub fn resistivity(t: f64) -> f64 {
    (1.0 + 0.005 * (t - 300.0)) / 18.0
}
pub fn conductor(card: &MaterialCard, t: f64) -> ResolvedConductorState {
    let basis = ElectricalPropertyBasis::Resistivity;
    ResolvedConductorState::try_from_material_with_basis(
        resolve_material_state_point(
            card,
            &point(t),
            &[basis.requirement().unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
        basis,
    )
    .unwrap()
}
pub fn problem(
    card: &Arc<MaterialCard>,
    positions: Vec<[f64; 3]>,
    tets: Vec<[u32; 4]>,
) -> SpatialHystereticProblem {
    problem_with_profile(card, positions, tets, false)
}
pub fn problem_with_viscosity(
    card: &Arc<MaterialCard>,
    positions: Vec<[f64; 3]>,
    tets: Vec<[u32; 4]>,
) -> SpatialHystereticProblem {
    problem_with_profile(card, positions, tets, true)
}
fn problem_with_profile(
    card: &Arc<MaterialCard>,
    positions: Vec<[f64; 3]>,
    tets: Vec<[u32; 4]>,
    kinetic: bool,
) -> SpatialHystereticProblem {
    let scalar = ResolvedScalarHysteresis::resolve(
        card.clone(),
        &point(300.0),
        roles(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let law = if kinetic {
        ResolvedVectorHysteresis::from_isotropic_source_with_viscosity(
            scalar,
            vec![viscosity_key()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
    } else {
        ResolvedVectorHysteresis::from_isotropic_source(scalar).unwrap()
    };
    let count = tets.len();
    SpatialHystereticProblem::new(
        positions,
        tets,
        vec![law; count],
        WeightedAssemblyLimits {
            max_cells: 100,
            max_dofs: 1000,
            max_triplets: 100000,
        },
        &mut || false,
    )
    .unwrap()
}
pub fn initial<'m>(
    problem: &'m SpatialHystereticProblem,
    conductor: Option<ResolvedConductorState>,
) -> SpatialHystereticInduction<'m> {
    let prescribed: BTreeSet<_> = problem
        .exterior_faces()
        .into_iter()
        .flat_map(|face| problem.face_edges(face).unwrap())
        .collect();
    SpatialHystereticInduction::new(
        problem,
        vec![conductor; problem.mesh().tets.len()],
        &prescribed.into_iter().collect::<Vec<_>>(),
        SpatialInductionLimits {
            maximum_unknowns: 2000,
            maximum_conducting_cells: 100,
            maximum_assembly_terms: 100000,
        },
        1e-10,
        &mut || false,
    )
    .unwrap()
}
pub fn field_controls() -> SpatialHystereticInductionControls {
    SpatialHystereticInductionControls {
        induction: SpatialNonlinearInductionControls {
            induction: SpatialInductionControls {
                relative_tolerance: 1e-11,
                absolute_residual_a: 1e-11,
                restart: 64,
                maximum_cycles: 100,
                maximum_krylov_scalars: 200000,
                current_tolerance_a: 2e-8,
                flux_tolerance_wb: 1e-10,
                energy_tolerance_j: 2e-8,
                maximum_cell_magnetic_peclet: 100.0,
            },
            maximum_newton_iterations: 50,
            maximum_backtracks: 40,
        },
        material: VectorHysteresisControls {
            relative_tolerance: 1e-12,
            absolute_field_tolerance_a_m: 1e-12,
            absolute_flux_tolerance_t: 1e-13,
            work_tolerance_j_m3: 2e-10,
            maximum_branch_iterations: 100,
            maximum_inverse_iterations: 80,
            maximum_backtracks: 48,
        },
    }
}
pub fn play(old_z: f64, h: f64) -> f64 {
    old_z
        .clamp(CHI * (h - PIN), CHI * (h + PIN))
        .clamp(-1.0, 1.0)
}
pub fn stored(h: f64, z: f64) -> f64 {
    0.5 * MU * h * h + z * z / (2.0 * CHI)
}
pub fn bracket(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) <= 0.0 && residual(hi) >= 0.0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}
// Exact Whitney integration and elimination of the four-node scalar potential:
// M=3*sigma/40, B=2*a, and M*delta(a)/dt + 2*H = J_component/2.
pub fn endpoint(old_b: f64, old_z: f64, sigma: f64, dt: f64, current: f64) -> (f64, f64, f64) {
    let h = bracket(-29.0, 29.0, |h| {
        h + 3.0 * sigma / (160.0 * dt) * (MU * h + play(old_z, h) - old_b) - 0.25 * current
    });
    let z = play(old_z, h);
    (h, z, MU * h + z)
}
