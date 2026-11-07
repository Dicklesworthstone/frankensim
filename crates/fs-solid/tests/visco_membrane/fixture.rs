//! Shared synthetic source and mesh for native membrane history/motion tests.
#![allow(dead_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::sheet::{
    AREAL_DENSITY, MembraneLaw, ResolvedSheet, SheetCapabilities, THICKNESS, resolve_sheet,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};
use fs_qty::Dims;
use fs_solid::membrane::visco::FiniteViscoStepControls;
use fs_solid::membrane::{MembraneBudget, MembraneDomain, MembranePatch, MembraneSolveSettings};
use fs_solid::mesh2::Mesh2;
use fs_solver::nonlinear::NewtonKrylovConfig;

pub const THICKNESS_M: f64 = 0.02;
pub const MASS_KG_M2: f64 = 0.8;

pub fn source_card() -> MaterialCard {
    let provenance = Provenance {
        source: "synthetic native membrane history fixture; not measured data".into(),
        license: "test-only".into(),
        artifact: None,
    };
    let pa = Dims([-1, 1, -2, 0, 0, 0]);
    let model = ConstitutiveModelCard {
        law: LawId("mooney-rivlin-decoupled".into()),
        law_version: 1,
        parameters: [
            ("c10", 10.0, pa),
            ("c01", 2.0, pa),
            ("kappa", 200.0, pa),
            ("stretch_limit", 4.0, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
        sources: vec![],
        provenance: provenance.clone(),
    };
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        (AREAL_DENSITY, MASS_KG_M2, Dims([-2, 1, 0, 0, 0, 0])),
        (THICKNESS, THICKNESS_M, Dims([1, 0, 0, 0, 0, 0])),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: provenance.clone(),
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic membrane".into(),
            phase: "solid".into(),
            process: "native history integration".into(),
            revision: 0,
        },
        claims,
        vec![model],
    )
    .unwrap()
}

pub fn source(exact: bool) -> ResolvedSheet {
    let card = source_card();
    let sheet = resolve_sheet(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        SheetCapabilities {
            membrane: false,
            bending: false,
        },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    if exact {
        sheet
            .with_incompressible_hyperelastic_membrane(&card, 0, local())
            .unwrap()
    } else {
        sheet.with_hyperelastic_membrane(&card, 0, local()).unwrap()
    }
}

pub fn material(sheet: &ResolvedSheet) -> ViscoHyperelastic {
    let MembraneLaw::Hyperelastic(equilibrium) = sheet.membrane_law().unwrap() else {
        unreachable!()
    };
    ViscoHyperelastic::new(
        equilibrium.model().law().clone(),
        vec![FiniteStrainBranch {
            energy: BranchEnergy::NeoHookean { mu_pa: 12.0 },
            tau_s: 0.25,
        }],
    )
    .unwrap()
}

pub fn local() -> PlaneStressSettings {
    PlaneStressSettings {
        relative_tolerance: 2e-13,
        ..Default::default()
    }
}

pub fn triangle() -> Mesh2 {
    Mesh2 {
        nodes: vec![[0.0, 0.0], [0.1, 0.0], [0.0, 0.1]],
        elems: vec![vec![0, 1, 2]],
        patches: vec![],
    }
}

pub fn diamond() -> Mesh2 {
    Mesh2 {
        nodes: vec![[0.0, 0.0], [0.1, 0.0], [0.0, 0.1], [-0.1, 0.0], [0.0, -0.1]],
        elems: vec![vec![0, 1, 2], vec![0, 2, 3], vec![0, 3, 4], vec![0, 4, 1]],
        patches: vec![],
    }
}

pub fn patch(mesh: Mesh2, sheet: ResolvedSheet) -> MembranePatch {
    let count = mesh.elems.len();
    MembranePatch::try_new(
        mesh,
        vec![sheet; count],
        vec![0.0; count],
        MembraneDomain {
            minimum_area_ratio: 0.01,
            maximum_abs_strain: 2.0,
        },
        MembraneBudget {
            max_nodes: 32,
            max_triangles: 32,
        },
        || false,
    )
    .unwrap()
}

pub fn settings() -> MembraneSolveSettings {
    MembraneSolveSettings {
        newton: NewtonKrylovConfig {
            relative_tolerance: 1e-10,
            absolute_tolerance: 1e-11,
            linear_restart: 16,
            ..Default::default()
        },
        max_iterations: 60,
        max_evaluations: 4096,
        max_krylov_restart: 64,
    }
}

pub fn controls() -> FiniteViscoStepControls {
    FiniteViscoStepControls {
        equilibrium_work_tolerance_j: 1e-9,
        maximum_endpoint_work_defect_j: 10.0,
    }
}

pub fn affine(patch: &MembranePatch, gradient: [[f64; 2]; 3]) -> Vec<[f64; 3]> {
    patch
        .mesh()
        .nodes
        .iter()
        .map(|x| core::array::from_fn(|i| gradient[i][0] * x[0] + gradient[i][1] * x[1]))
        .collect()
}

pub fn prescribed(u: &[[f64; 3]]) -> Vec<(usize, f64)> {
    u.iter().flatten().copied().enumerate().collect()
}

pub fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}
