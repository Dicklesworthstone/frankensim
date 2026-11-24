#![allow(dead_code)]
//! Manufactured conducting magnetic play sources, meshes, problems and
//! controls, shared by `hysteretic_eddy` and `winding_component`.

use fs_em::planar::{
    eddy::{EddyConductor, hysteretic::HystereticEddyControls},
    nonlinear::{
        NonlinearMagneticControls,
        hysteretic::{
            PlanarHystereticControls, PlanarHystereticMagnetostatics, PlanarHystereticState,
        },
    },
};
use fs_em::winding::hysteretic::HystereticWindingControls;
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
                HysteresisBranchSources, HysteresisSources, ResolvedScalarHysteresis,
                ResolvedVectorHysteresis,
            },
            vector::VectorHysteresisControls,
        },
        permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
    },
    state_point::{MaterialPropertySelection, resolve_material_state_point},
};
use fs_qty::{QuantitySpec, Temperature};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::{CheckedCgConfig, NewtonKrylovConfig};
use std::sync::Arc;

pub const MU: f64 = 0.01;
pub const CHI: f64 = 0.2;
pub const PINNING: f64 = 2.0;
pub const SATURATION: f64 = 1.0;
pub const DIRECTIONS: [[f64; 2]; 4] = [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]];

pub fn sources(
    mu: f64,
    pinning: bool,
    sigma: f64,
) -> (ResolvedVectorHysteresis, ResolvedConductorState) {
    sources_with_viscosity(mu, pinning, sigma, None)
}

pub fn sources_with_viscosity(
    mu: f64,
    pinning: bool,
    sigma: f64,
    viscosity: Option<f64>,
) -> (ResolvedVectorHysteresis, ResolvedConductorState) {
    sources_at(mu, pinning, sigma, viscosity, None)
}

/// The same sources, or, with `[kelvin, low, high]`, those of a card whose
/// every claim holds for `T` from `low` to `high` kelvin, resolved at
/// `T = kelvin`.
pub fn sources_at(
    mu: f64,
    pinning: bool,
    sigma: f64,
    viscosity: Option<f64>,
    temperature: Option<[f64; 3]>,
) -> (ResolvedVectorHysteresis, ResolvedConductorState) {
    let (validity, point) = support(temperature);
    let basis = ElectricalPropertyBasis::Conductivity;
    let mut roles = HysteresisSources {
        background_permeability: PropertyKey::new("background_mu", MAGNETIC_PERMEABILITY_DIMS),
        field_limit: PropertyKey::new("field_limit", FIELD_STRENGTH_DIMS),
        branches: vec![],
    };
    let mut values = vec![
        (roles.background_permeability.clone(), mu),
        (roles.field_limit.clone(), 30.0),
        (
            PropertyKey::new(basis.property_name(), basis.dimensions()),
            sigma,
        ),
    ];
    if pinning {
        let branch = HysteresisBranchSources {
            susceptibility: PropertyKey::new("branch_slope", MAGNETIC_PERMEABILITY_DIMS),
            pinning_field: PropertyKey::new("branch_pinning", FIELD_STRENGTH_DIMS),
            saturation_polarization: PropertyKey::new("branch_saturation", FLUX_DENSITY_DIMS),
        };
        values.extend([
            (branch.susceptibility.clone(), CHI),
            (branch.pinning_field.clone(), PINNING),
            (branch.saturation_polarization.clone(), SATURATION),
        ]);
        roles.branches.push(branch);
    }
    let viscosity_key = PropertyKey::new(
        "intrinsic_branch_viscosity",
        fs_material::field_storage::hysteresis::source::MAGNETIC_VISCOSITY_DIMS,
    );
    if let Some(viscosity) = viscosity {
        assert!(pinning, "viscosity fixture needs one actual branch");
        values.push((viscosity_key.clone(), viscosity));
    }
    let mut claims = ClaimSet::new();
    for (key, value) in values {
        let dims = key.dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Scalar { value, dims },
                validity: validity.clone(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured conducting magnetic play fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic conducting magnetic law".into(),
                phase: "solid".into(),
                process: "manufactured transient field".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    );
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        resolve_material_state_point(
            &card,
            &point,
            &[basis.requirement().unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
        basis,
    )
    .unwrap();
    let scalar = ResolvedScalarHysteresis::resolve(
        card,
        &point,
        roles,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let magnetic = if viscosity.is_some() {
        ResolvedVectorHysteresis::from_isotropic_source_with_viscosity(
            scalar,
            vec![viscosity_key],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
    } else {
        ResolvedVectorHysteresis::from_isotropic_source(scalar).unwrap()
    };
    (magnetic, conductor)
}

/// Each claim's validity and the query point, bounded in `T` or not.
fn support(temperature: Option<[f64; 3]>) -> (ValidityDomain, QueryPoint) {
    let kelvin = || QuantitySpec::dimensional(Temperature::DIMS);
    match temperature {
        Some([at, low, high]) => (
            ValidityDomain::unconstrained().with_quantity("T", kelvin(), low, high),
            QueryPoint::new().with_quantity("T", kelvin(), at).unwrap(),
        ),
        None => (ValidityDomain::unconstrained(), QueryPoint::new()),
    }
}

pub fn square(depth: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "hysteretic-conductor-square",
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap()
}

pub fn paired_squares(depth: f64) -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "balanced-hysteretic-conductor",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [0.5, 0.5],
            [1.5, 0.5],
        ],
        vec![
            [0, 1, 6],
            [1, 2, 6],
            [2, 3, 6],
            [3, 0, 6],
            [1, 4, 7],
            [4, 5, 7],
            [5, 2, 7],
            [2, 1, 7],
        ],
        Metric2::planar(depth).unwrap(),
    )
    .unwrap()
}

pub fn problem(
    mesh: TriComplex2,
    mu: f64,
    pinning: bool,
    sigma: f64,
) -> (PlanarHystereticMagnetostatics, Vec<EddyConductor>) {
    let (magnetic, conductor) = sources(mu, pinning, sigma);
    problem_from_sources(mesh, magnetic, conductor)
}

pub fn problem_from_sources(
    mesh: TriComplex2,
    magnetic: ResolvedVectorHysteresis,
    conductor: ResolvedConductorState,
) -> (PlanarHystereticMagnetostatics, Vec<EddyConductor>) {
    let count = mesh.faces().len();
    problem_from_materials(mesh, vec![magnetic; count], conductor)
}

/// One magnetic law per triangle, in mesh order; every triangle conducts.
pub fn problem_from_materials(
    mesh: TriComplex2,
    materials: Vec<ResolvedVectorHysteresis>,
    conductor: ResolvedConductorState,
) -> (PlanarHystereticMagnetostatics, Vec<EddyConductor>) {
    let count = mesh.faces().len();
    let owner = PlanarHystereticMagnetostatics::new(
        mesh,
        materials,
        WeightedAssemblyLimits {
            max_cells: 100,
            max_dofs: 100,
            max_triplets: 10000,
        },
        &mut || false,
    )
    .unwrap();
    let conductors = vec![conductor; count]
        .into_iter()
        .enumerate()
        .map(|(triangle, conductor)| EddyConductor {
            triangle,
            region: 7,
            conductor,
        })
        .collect();
    (owner, conductors)
}

pub fn controls() -> HystereticEddyControls {
    HystereticEddyControls {
        magnetic: PlanarHystereticControls {
            nonlinear: NonlinearMagneticControls {
                absolute_residual_tolerance_a_m: 1e-9,
                relative_residual_tolerance: 1e-10,
                max_newton_iterations: 50,
                max_backtracks: 48,
                linear: CheckedCgConfig {
                    tolerance: 1e-13,
                    max_iterations: 100,
                    max_corrections: 2,
                },
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
            energy_tolerance_j: 2e-8,
        },
        current_balance_tolerance_a: 1e-8,
    }
}

/// The winding circuit's terminal controls around [`controls`].
pub fn terminal_controls() -> HystereticWindingControls {
    HystereticWindingControls {
        field: controls(),
        maximum_field_trials: 1000,
        maximum_newton_iterations: 30,
        difference_step_a: 1e-5,
        terminal_voltage_tolerance_v: 1e-8,
        power_tolerance_w: 1e-8,
        energy_tolerance_j: 2e-8,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}

pub fn boundary(count: usize, gauge: f64) -> Vec<(usize, f64)> {
    (0..count).map(|v| (v, gauge)).collect()
}

pub fn magnetized_square(
    owner: &PlanarHystereticMagnetostatics,
    gauge: f64,
    h: f64,
) -> PlanarHystereticState<'_> {
    let zero = owner.demagnetized(gauge, &mut || false).unwrap();
    owner
        .candidate(
            &zero,
            &[6.0 * h; 4],
            &boundary(4, gauge),
            controls().magnetic,
            &mut || false,
        )
        .unwrap()
        .into_parts()
        .0
}

pub fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-8 * (1.0 + b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}
