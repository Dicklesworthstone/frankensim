//! Public enclosure-driver regressions: complete solid residuals and reusable radiosity.
mod support;
use fs_conduction::{
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    CoupledRadiationConfig, GrayDiffuseEnclosure, InitialGuess, LinearConfig,
    Nonlinearity, RadiationSurface, ScalarField, SolveConfig, StopRule,
    SurfaceEmissivity, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder,
    ViewFactorMatrix, EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY,
    solve_with_gray_diffuse_enclosure,
};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_rep_mesh::TetComplex;
use support::with_cx;

// Deliberately reuse the existing radiation test's declared two-slab model.
// These factors are an algebraic fixture, not a finite-geometry visibility claim.
fn emissivity(name: &str, value: f64, finish: &str) -> SurfaceEmissivity {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Scalar { value, dims: EMISSIVITY_DIMS },
        validity: ValidityDomain::unconstrained().with("T", 250.0, 500.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: Provenance {
            source: "enclosure solver regression fixture".into(),
            license: "internal-test-use".into(), artifact: None,
        },
    }).unwrap();
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "fixture-alloy".into(), phase: "solid".into(),
        process: finish.into(), revision: 0,
    }, claims, Vec::new()).unwrap();
    SurfaceEmissivity::from_card(name, &card, 350.0, SelectionPolicy::SingleClaimOnly).unwrap()
}

fn config() -> SolveConfig {
    SolveConfig {
        nonlinearity: Nonlinearity::FixedPoint { relaxation: 1.0, max_backtracks: 8 },
        stop: StopRule { residual_rtol: 1e-10, residual_atol: 1e-24,
            step_atol: 0.0, max_iterations: 20 },
        linear: LinearConfig { tolerance: 1e-12, max_iterations: 50_000, restart: 60 },
        initial: InitialGuess::DirichletMean,
    }
}

fn two_disconnected_slabs() -> ConductionMesh {
    let (a, mut positions) = box_grid([2, 2, 2], [0.5, 1.0, 1.0]);
    let (b, other_positions) = box_grid([2, 2, 2], [0.5, 1.0, 1.0]);
    let offset = positions.len() as u32;
    positions.extend(other_positions.into_iter().map(|mut p| { p[0] += 1.5; p }));
    let mut tets = a.tets;
    tets.extend(b.tets.into_iter().map(|tet| tet.map(|v| v + offset)));
    ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap()
}

fn two_slab_enclosure(mesh: &ConductionMesh) -> GrayDiffuseEnclosure {
    let hot = RadiationSurface::new(mesh, "hot-facing",
        |f| on_box_face(f.centroid[0], 0.5), emissivity("hot-facing", 0.8, "hot")).unwrap();
    let cold = RadiationSurface::new(mesh, "cold-facing",
        |f| on_box_face(f.centroid[0], 1.5), emissivity("cold-facing", 0.8, "cold")).unwrap();
    let matrix = ViewFactorMatrix::infinite_parallel_plates(hot.area_m2()).unwrap();
    GrayDiffuseEnclosure::new(vec![hot, cold], matrix).unwrap()
}

fn two_slab_boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region("hot-reservoir", |f| on_box_face(f.centroid[0], 0.0), ThermalBc::dirichlet(400.0).unwrap()).unwrap()
        .region("cold-reservoir", |f| on_box_face(f.centroid[0], 2.0), ThermalBc::dirichlet(300.0).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap()
}

#[path = "enclosure_solver/solid_acceptance.rs"]
mod solid_acceptance;
