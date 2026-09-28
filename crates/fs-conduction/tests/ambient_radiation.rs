//! G1 slab heat balance and G0 admission for the shared ambient patch producer.

mod support;

use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    AmbientRadiationConfig, AmbientRadiationPatch, AmbientRadiationSolution, ConductionError,
    ConductionMesh, ConductionProblem, ConductivityModel, EMISSIVITY_DIMS, ElementMaterials,
    InitialGuess, InterfaceResistance, InterfaceSurface, MaterialId, MaterialTable,
    STEFAN_BOLTZMANN_W_M2_K4, SURFACE_EMISSIVITY_PROPERTY, ScalarField, SolveConfig,
    SurfaceEmissivity, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder, ThermalInterfaces,
    solve_with_ambient_radiation,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_rep_mesh::TetComplex;
use support::{with_cancelled_cx, with_cx};

fn material_state(name: &str) -> MaterialStateId {
    MaterialStateId {
        chemistry: name.into(),
        phase: "solid".into(),
        process: "test".into(),
        revision: 0,
    }
}

fn claim(name: &str, dims: fs_qty::Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with("T", 250.0, 500.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: Provenance {
            source: "analytic radiation fixture".into(),
            license: "internal-test-use".into(),
            artifact: None,
        },
    }
}

fn card(epsilon: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(claim(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, epsilon))
        .unwrap();
    MaterialCard::assemble(material_state("gray-surface"), claims, Vec::new()).unwrap()
}

fn patch(epsilon: f64, ambient: f64) -> AmbientRadiationPatch {
    let emissivity = SurfaceEmissivity::from_card(
        "cooler",
        &card(epsilon),
        350.0,
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    AmbientRadiationPatch::new("cooler", emissivity, ambient).unwrap()
}

fn controls() -> SolveConfig {
    let mut config = SolveConfig::default();
    config.initial = InitialGuess::Uniform(350.0);
    config.stop.step_atol = 0.0;
    config.stop.residual_rtol = 1e-11;
    config.linear.tolerance = 1e-12;
    config
}

fn boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "hot",
            |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::dirichlet(400.0).unwrap(),
        )
        .unwrap()
        .region(
            "cooler",
            |face| on_box_face(face.centroid[0], 1.0),
            ThermalBc::robin(5.0, 300.0).unwrap(),
        )
        .unwrap()
        .adiabatic_remainder()
        .finish()
        .unwrap()
}

/// Independent scalar balance of a unit-area slab with effective thermal
/// conductance G. The bisection evaluates fourth powers directly, independently
/// of the production secant factorization and FEM solve.
fn wall_oracle(conductance: f64, epsilon: f64, ambient: f64) -> f64 {
    let mut lower = 250.0_f64;
    let mut upper = 500.0_f64;
    for _ in 0..80 {
        let t = (lower + upper) / 2.0;
        let residual = conductance * (400.0 - t)
            - 5.0 * (t - 300.0)
            - epsilon * STEFAN_BOLTZMANN_W_M2_K4 * (t.powi(4) - ambient.powi(4));
        if residual > 0.0 {
            lower = t;
        } else {
            upper = t;
        }
    }
    (lower + upper) / 2.0
}

fn check_wall(solution: &AmbientRadiationSolution, conductance: f64, epsilon: f64, ambient: f64) {
    let t = wall_oracle(conductance, epsilon, ambient);
    let row = &solution.radiation.patches[0];
    assert!(
        (row.mean_surface_temperature_k - t).abs() < 2e-7,
        "native surface {} K versus scalar balance {t} K",
        row.mean_surface_temperature_k
    );
    let radiation = epsilon * STEFAN_BOLTZMANN_W_M2_K4 * (t.powi(4) - ambient.powi(4));
    assert!((row.nonlinear_heat_w - radiation).abs() < 2e-6);
    assert!((solution.convective_out_w - 5.0 * (t - 300.0)).abs() < 2e-6);
    assert!(
        (solution.conduction.report.energy.robin_out_w - conductance * (400.0 - t)).abs() < 2e-6
    );
    assert!(row.heat_mismatch_w <= row.heat_tolerance_w);
    assert!(
        solution.radiation.max_temperature_change_k
            <= solution.radiation.config.temperature_tolerance_k
    );
    assert!(solution.radiation.decomposition_residual_w.abs() < 1e-8);
    assert_eq!(solution.convective_robin_fluxes[0].mean_htc_w_per_m2_k, 5.0);
    assert_eq!(
        solution.convective_robin_fluxes[0].mean_reference_temperature_k,
        300.0
    );
    assert!(solution.radiation.krylov_iterations > 0);
}

#[test]
fn ambient_patch_slab_matches_scalar_balance_and_changes_with_ambient_and_emissivity() {
    with_cx(|cx| {
        let (complex, positions) = box_grid([4, 1, 1], [1.0; 3]);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let boundary = boundary(&mesh);
        let material = ConductivityModel::isotropic_declared(20.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let solve = |epsilon, ambient| {
            solve_with_ambient_radiation(
                cx,
                problem,
                None,
                &[patch(epsilon, ambient)],
                controls(),
                AmbientRadiationConfig::default(),
            )
            .unwrap()
        };
        let nominal = solve(0.8, 300.0);
        let lower_emissivity = solve(0.3, 300.0);
        let hot_surroundings = solve(0.8, 450.0);
        for (solution, epsilon, ambient) in [
            (&nominal, 0.8, 300.0),
            (&lower_emissivity, 0.3, 300.0),
            (&hot_surroundings, 0.8, 450.0),
        ] {
            check_wall(solution, 20.0, epsilon, ambient);
        }
        assert!(
            nominal.radiation.patches[0].mean_surface_temperature_k
                < lower_emissivity.radiation.patches[0].mean_surface_temperature_k
        );
        assert!(
            hot_surroundings.radiation.nonlinear_radiation_out_w < 0.0,
            "hot surrounding radiation must heat the solid, not be clamped away"
        );
        assert_eq!(
            nominal,
            solve(0.8, 300.0),
            "repeated physical solve must be identical"
        );
        let replay = fs_conduction::solve(
            cx,
            ConductionProblem {
                boundary: &nominal.combined_boundary,
                ..problem
            },
            nominal.final_solve_config.clone(),
        )
        .unwrap();
        assert_eq!(
            replay, nominal.conduction,
            "exact final binding reconstructs the same solid field"
        );
        let bare = fs_conduction::solve(cx, problem, controls()).unwrap();
        assert!((bare.report.robin_fluxes[0].mean_wall_temperature_k - 380.0).abs() < 1e-8);
        assert!(nominal.radiation.patches[0].mean_surface_temperature_k < 379.0);
    });
}

#[test]
fn ambient_patch_retains_heterogeneous_material_and_finite_contact_drop() {
    with_cx(|cx| {
        let (left, mut positions) = box_grid([2, 1, 1], [0.5, 1.0, 1.0]);
        let (right, other) = box_grid([2, 1, 1], [0.5, 1.0, 1.0]);
        let offset = positions.len() as u32;
        let left_elements = left.tets.len();
        let mut tets = left.tets;
        tets.extend(right.tets.iter().map(|tet| tet.map(|v| v + offset)));
        positions.extend(other.into_iter().map(|[x, y, z]| [x + 0.5, y, z]));
        let mesh =
            ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap();
        let boundary = boundary(&mesh);
        let surface = |name| SurfaceSpec {
            material: material_state(name),
            texture_frame: "x-normal".into(),
        };
        let mut claims = ClaimSet::new();
        claims
            .insert_claim(claim(
                AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
                0.05,
            ))
            .unwrap();
        let contact = InterfaceSystemCard::assemble(
            surface("left"),
            surface("right"),
            SystemContext {
                medium: "dry".into(),
                third_body: None,
                environment: "fixture".into(),
                history: "constant-contact".into(),
            },
            claims,
            Vec::new(),
        )
        .unwrap();
        let resistance = InterfaceResistance::from_card(
            "joint",
            &contact,
            &QueryPoint::new().with("T", 350.0).unwrap(),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
        let pairs = ThermalInterfaces::coincident_face_pairs(&mesh).unwrap();
        let interfaces = ThermalInterfaces::new(
            &mesh,
            &boundary,
            vec![InterfaceSurface::new("joint", pairs, resistance).unwrap()],
        )
        .unwrap();
        let table = MaterialTable::new([
            (
                MaterialId(1),
                ConductivityModel::isotropic_declared(10.0).unwrap(),
            ),
            (
                MaterialId(2),
                ConductivityModel::isotropic_declared(40.0).unwrap(),
            ),
        ])
        .unwrap();
        let assignment = ElementMaterials::new(
            table,
            (0..mesh.element_count())
                .map(|e| MaterialId(if e < left_elements { 1 } else { 2 }))
                .collect(),
        )
        .unwrap();
        let fallback = ConductivityModel::isotropic_declared(999.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let result = solve_with_ambient_radiation(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &fallback,
                element_materials: Some(&assignment),
                source: &source,
            },
            Some(&interfaces),
            &[patch(0.8, 300.0)],
            controls(),
            AmbientRadiationConfig::default(),
        )
        .unwrap();
        let conductance = 1.0 / (0.5 / 10.0 + 0.05 + 0.5 / 40.0);
        check_wall(&result, conductance, 0.8, 300.0);
        assert_eq!(
            result.conduction.report.element_material_identity,
            Some(assignment.identity())
        );
        let contact = &result.conduction.report.interface_fluxes[0];
        let heat = conductance * (400.0 - wall_oracle(conductance, 0.8, 300.0));
        assert!((contact.heat_rate_a_to_b_w.abs() - heat).abs() < 2e-6);
        assert!((contact.mean_jump_k.abs() - heat * 0.05).abs() < 2e-7);
    });
}

#[test]
fn ambient_patch_refuses_bad_bindings_exhaustion_cancel_and_material_extrapolation() {
    let (complex, positions) = box_grid([2, 1, 1], [1.0; 3]);
    let mesh = ConductionMesh::new(complex, positions).unwrap();
    let boundary = boundary(&mesh);
    let material = ConductivityModel::isotropic_declared(20.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        element_materials: None,
        source: &source,
    };
    with_cx(|cx| {
        let nominal = patch(0.8, 300.0);
        for patches in [
            Vec::new(),
            vec![nominal.clone(), nominal.clone()],
            vec![AmbientRadiationPatch::new("hot", nominal.emissivity().clone(), 300.0).unwrap()],
            vec![
                AmbientRadiationPatch::new("missing", nominal.emissivity().clone(), 300.0).unwrap(),
            ],
        ] {
            assert!(matches!(
                solve_with_ambient_radiation(
                    cx,
                    problem,
                    None,
                    &patches,
                    controls(),
                    AmbientRadiationConfig::default()
                ),
                Err(ConductionError::Radiation { .. })
            ));
        }
        let one = AmbientRadiationConfig {
            max_iterations: 1,
            ..AmbientRadiationConfig::default()
        };
        let error =
            solve_with_ambient_radiation(cx, problem, None, &[nominal.clone()], controls(), one)
                .unwrap_err();
        assert!(matches!(
            error,
            ConductionError::AmbientRadiationNotConverged { iterations: 1, .. }
        ));
        let zero = AmbientRadiationConfig {
            max_iterations: 0,
            ..one
        };
        assert!(
            solve_with_ambient_radiation(cx, problem, None, &[nominal.clone()], controls(), zero)
                .is_err()
        );
        assert!(nominal.secant_coefficient_w_m2_k(501.0).is_err());
        assert!(nominal.heat_flux_w_m2(f64::NAN).is_err());
        assert!(AmbientRadiationPatch::new("cooler", nominal.emissivity().clone(), 0.0).is_err());
        // Reservoir temperature is not the surface material's temperature:
        // its 250 K lower validity limit must not incorrectly ban a 200 K sink.
        assert!(AmbientRadiationPatch::new("cooler", nominal.emissivity().clone(), 200.0).is_ok());
        let strict_watts = AmbientRadiationConfig {
            temperature_tolerance_k: 1e6,
            balance_tolerance_w: 1e-10,
            balance_relative_tolerance: 0.0,
            ..one
        };
        assert!(
            solve_with_ambient_radiation(cx, problem, None, &[nominal], controls(), strict_watts)
                .is_err(),
            "a loose temperature gate must not bypass the independent nonlinear watt gate"
        );
    });
    with_cancelled_cx(|cx| {
        assert!(matches!(
            solve_with_ambient_radiation(
                cx,
                problem,
                None,
                &[patch(0.8, 300.0)],
                controls(),
                AmbientRadiationConfig::default()
            ),
            Err(ConductionError::Cancelled { .. })
        ));
    });
}

#[test]
fn emissivity_claim_pin_selects_exact_evidence_and_refuses_missing_pin() {
    let mut claims = ClaimSet::new();
    let first = claims
        .insert_claim(claim(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, 0.3))
        .unwrap();
    let second = claims
        .insert_claim(claim(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS, 0.8))
        .unwrap();
    let card =
        MaterialCard::assemble(material_state("ambiguous-finish"), claims, Vec::new()).unwrap();
    assert!(
        SurfaceEmissivity::from_card("cooler", &card, 350.0, SelectionPolicy::SingleClaimOnly)
            .is_err()
    );
    let selected = SurfaceEmissivity::from_card_pinned("cooler", &card, 350.0, second).unwrap();
    assert_eq!(selected.value(), 0.8);
    assert_eq!(selected.receipt().selected, second);
    assert_eq!(
        SurfaceEmissivity::from_card_pinned("cooler", &card, 350.0, first)
            .unwrap()
            .value(),
        0.3
    );
    assert!(SurfaceEmissivity::from_card_pinned("cooler", &card, 600.0, second).is_err());
    let missing = fs_matdb::ClaimId(fs_blake3::hash_bytes(b"missing-emissivity-claim"));
    assert!(SurfaceEmissivity::from_card_pinned("cooler", &card, 350.0, missing).is_err());
}
