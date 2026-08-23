//! G1 analytic and G0 refusal coverage for matching-P1 thermal contact.
//!
//! The execution fixture is two disconnected tetrahedral slabs whose
//! coincident `x = 1 m` traces use distinct vertices.  The interface card
//! supplies `R'' = 0.1 m^2 K/W`; with unit area and two `L/(kA) = 0.1 K/W`
//! slab terms, the Level-A reference is exactly `0.3 K/W`.

mod support;

use fs_conduction::bc::{ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::ConductivityModel;
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{
    ConductionProblem, InitialGuess, LinearConfig, Nonlinearity, SolveConfig, StopRule, solve,
    solve_with_interfaces,
};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
    ConductionError, InterfaceFacePair, InterfaceResistance, InterfaceSurface,
    ResistanceUncertainty, ResistanceValueOrigin, SeriesThermalResistance, ThermalInterfaces,
    ThermalResistanceTerm,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel,
};
use fs_rep_mesh::TetComplex;
use fs_vvreg::thermal_level_a::{
    ThermalLevelAAcceptance, ThermalLevelAKind, thermal_level_a_cases,
};
use support::with_cx;

const K: f64 = 10.0;
const T_HOT: f64 = 330.0;
const T_COLD: f64 = 300.0;
const AREA_SPECIFIC_RESISTANCE: f64 = 0.1;

fn config() -> SolveConfig {
    SolveConfig {
        nonlinearity: Nonlinearity::FixedPoint {
            relaxation: 1.0,
            max_backtracks: 8,
        },
        stop: StopRule {
            residual_rtol: 1e-11,
            residual_atol: 1e-24,
            step_atol: 0.0,
            max_iterations: 12,
        },
        linear: LinearConfig {
            tolerance: 1e-13,
            max_iterations: 60_000,
            restart: 60,
        },
        initial: InitialGuess::DirichletMean,
    }
}

fn surface(material: &str, texture: &str) -> SurfaceSpec {
    SurfaceSpec {
        material: MaterialStateId {
            chemistry: material.to_string(),
            phase: "solid".to_string(),
            process: "as-fixtured".to_string(),
            revision: 0,
        },
        texture_frame: texture.to_string(),
    }
}

fn contact_card(uncertainty: UncertaintyModel) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,
                AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            ),
            value: PropertyValue::Scalar {
                value: AREA_SPECIFIC_RESISTANCE,
                dims: AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
            },
            validity: ValidityDomain::unconstrained(),
            uncertainty,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "G1 two-slab declared contact fixture".to_string(),
                license: "internal-test-use".to_string(),
                artifact: None,
            },
        })
        .expect("contact resistance claim inserts");
    InterfaceSystemCard::assemble(
        surface("solid-a", "interface-normal-plus-x"),
        surface("solid-b", "interface-normal-minus-x"),
        SystemContext {
            medium: "dry".to_string(),
            third_body: Some("declared-contact-layer".to_string()),
            environment: "vacuum".to_string(),
            history: "unaged".to_string(),
        },
        claims,
        Vec::new(),
    )
    .expect("interface card assembles")
}

fn empty_contact_card() -> InterfaceSystemCard {
    InterfaceSystemCard::assemble(
        surface("solid-a", "interface-normal-plus-x"),
        surface("solid-b", "interface-normal-minus-x"),
        SystemContext {
            medium: "dry".to_string(),
            third_body: Some("declared-contact-layer".to_string()),
            environment: "vacuum".to_string(),
            history: "unaged".to_string(),
        },
        ClaimSet::new(),
        Vec::new(),
    )
    .expect("empty interface card is structurally valid")
}

fn slab_chain_mesh(n: usize, slab_count: usize) -> (ConductionMesh, usize) {
    assert!(
        slab_count > 0,
        "a slab chain must contain at least one slab"
    );
    let mut tets = Vec::new();
    let mut positions = Vec::new();
    let mut vertices_per_slab = 0usize;
    for slab in 0..slab_count {
        let (complex, slab_positions) = box_grid([n, n, n], [1.0, 1.0, 1.0]);
        if slab == 0 {
            vertices_per_slab = slab_positions.len();
        } else {
            assert_eq!(
                slab_positions.len(),
                vertices_per_slab,
                "identical box-grid inputs must have identical vertex counts"
            );
        }
        let offset =
            u32::try_from(positions.len()).expect("fixture vertex count fits in a u32 index");
        tets.extend(
            complex
                .tets
                .into_iter()
                .map(|tet| tet.map(|vertex| vertex + offset)),
        );
        positions.extend(
            slab_positions
                .into_iter()
                .map(|[x, y, z]| [x + slab as f64, y, z]),
        );
    }
    let complex = TetComplex::from_tets(positions.len(), tets);
    (
        ConductionMesh::new(complex, positions).expect("slab-chain mesh"),
        vertices_per_slab,
    )
}

fn two_slab_mesh(n: usize) -> (ConductionMesh, usize) {
    slab_chain_mesh(n, 2)
}

fn two_tet_contact_mesh() -> ConductionMesh {
    let positions = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0],
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    ConductionMesh::new(
        TetComplex::from_tets(positions.len(), vec![[0, 1, 2, 3], [4, 5, 6, 7]]),
        positions,
    )
    .expect("two-tet matching-contact mesh")
}

fn coincident_base_tet_mesh(apex_heights: &[f64]) -> ConductionMesh {
    assert!(
        apex_heights.len() >= 2,
        "the refusal fixture needs at least two coincident traces"
    );
    let mut positions = Vec::with_capacity(4 * apex_heights.len());
    let mut tets = Vec::with_capacity(apex_heights.len());
    for &height in apex_heights {
        let offset = u32::try_from(positions.len()).expect("fixture vertex count fits u32");
        positions.extend([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, height],
        ]);
        tets.push([offset, offset + 1, offset + 2, offset + 3]);
    }
    ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions)
        .expect("coincident-base tetrahedra remain individually nondegenerate")
}

fn outer_boundary(mesh: &ConductionMesh, cold_x: f64) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "hot",
            |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::dirichlet(T_HOT).expect("hot condition"),
        )
        .expect("hot boundary")
        .region(
            "cold",
            |face| on_box_face(face.centroid[0], cold_x),
            ThermalBc::dirichlet(T_COLD).expect("cold condition"),
        )
        .expect("cold boundary")
        .adiabatic_remainder()
        .finish()
        .expect("complete boundary partition")
}

fn boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    outer_boundary(mesh, 2.0)
}

fn boundary_with_occupied_interface(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .region(
            "hot",
            |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::dirichlet(T_HOT).expect("hot condition"),
        )
        .expect("hot boundary")
        .region(
            "cold",
            |face| on_box_face(face.centroid[0], 2.0),
            ThermalBc::dirichlet(T_COLD).expect("cold condition"),
        )
        .expect("cold boundary")
        .region(
            "occupied-contact",
            |face| on_box_face(face.centroid[0], 1.0),
            ThermalBc::neumann(0.0).expect("zero-flux condition"),
        )
        .expect("occupied interface boundary")
        .adiabatic_remainder()
        .finish()
        .expect("complete occupied-interface boundary")
}

fn adiabatic_boundary(mesh: &ConductionMesh) -> ThermalBoundary {
    ThermalBoundaryBuilder::new(mesh)
        .adiabatic_remainder()
        .finish()
        .expect("explicit all-adiabatic boundary")
}

fn resistance(card: &InterfaceSystemCard) -> InterfaceResistance {
    InterfaceResistance::from_card(
        "bondline",
        card,
        &QueryPoint::new(),
        SelectionPolicy::SingleClaimOnly,
    )
    .expect("contact resistance query")
}

fn oriented_pairs(mesh: &ConductionMesh) -> Vec<InterfaceFacePair> {
    ThermalInterfaces::coincident_face_pairs(mesh)
        .expect("coincident pairs")
        .into_iter()
        .map(|pair| {
            if mesh.boundary()[pair.side_a].outward_normal[0] > 0.0 {
                pair
            } else {
                InterfaceFacePair {
                    side_a: pair.side_b,
                    side_b: pair.side_a,
                }
            }
        })
        .collect()
}

fn oriented_pairs_at_x(mesh: &ConductionMesh, x: f64) -> Vec<InterfaceFacePair> {
    oriented_pairs(mesh)
        .into_iter()
        .filter(|pair| on_box_face(mesh.boundary()[pair.side_a].centroid[0], x))
        .collect()
}

fn solve_contact_fixture(
    mesh: &ConductionMesh,
    boundary: &ThermalBoundary,
    interfaces: &ThermalInterfaces,
) -> fs_conduction::ConductionSolution {
    let material = ConductivityModel::isotropic_declared(K).expect("material");
    let source = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        solve_with_interfaces(
            cx,
            ConductionProblem { element_materials: None, mesh,
            boundary,
            material: &material,
            source: &source, },
            interfaces,
            config(),
        )
        .expect("contact solve")
    })
}

fn binding_validation_error(
    mesh: &ConductionMesh,
    boundary: &ThermalBoundary,
    interfaces: &ThermalInterfaces,
) -> ConductionError {
    let material = ConductivityModel::isotropic_declared(K).expect("material");
    let source = ScalarField::Uniform(0.0);
    let temperature = vec![0.5 * (T_HOT + T_COLD); mesh.vertex_count()];
    with_cx(|cx| {
        match fs_conduction::assemble_operator_with_interfaces(
            cx,
            mesh,
            boundary,
            &material,
            &source,
            &temperature,
            interfaces,
        ) {
            Err(error) => error,
            Ok(_) => panic!("the changed mesh/boundary binding must refuse before assembly"),
        }
    })
}

fn assert_interface_refusal(error: ConductionError, expected_interface: &str, expected_what: &str) {
    match error {
        ConductionError::Interface {
            interface,
            what,
            fix,
        } => {
            assert_eq!(interface, expected_interface);
            assert!(
                what.contains(expected_what),
                "expected diagnosis containing {expected_what:?}, got {what:?}"
            );
            assert!(
                !fix.trim().is_empty(),
                "an interface refusal must retain an actionable fix"
            );
        }
        other => panic!("expected a typed interface refusal, got {other:?}"),
    }
}

fn level_a_contact_reference() -> (&'static fs_vvreg::thermal_level_a::ThermalLevelACase, f64) {
    let case = thermal_level_a_cases()
        .iter()
        .find(|case| case.id == "thermal-a-contact-series")
        .expect("Level-A contact row");
    assert_eq!(case.kind, ThermalLevelAKind::AnalyticReference);
    assert_eq!(case.metric, "thermal-resistance");
    (case, case.reference_value_si)
}

fn level_a_tolerance_error_and_limit(
    case: &fs_vvreg::thermal_level_a::ThermalLevelACase,
    observed: f64,
) -> (f64, f64) {
    let ThermalLevelAAcceptance::Tolerance { atol, rtol } = case.acceptance else {
        panic!("the Level-A contact row must retain a numerical tolerance");
    };
    let error = (observed - case.reference_value_si).abs();
    let limit = atol + rtol * case.reference_value_si.abs();
    assert!(
        atol.is_finite() && atol >= 0.0 && rtol.is_finite() && rtol >= 0.0,
        "the registry tolerance itself must be finite and non-negative"
    );
    assert!(
        limit.is_finite() && limit > 0.0,
        "the registry row must expose a non-vacuous positive comparison limit"
    );
    (error, limit)
}

#[test]
fn two_slab_contact_matches_level_a_series_and_retains_receipt() {
    let (mesh, left_vertex_count) = two_slab_mesh(4);
    let boundary = boundary(&mesh);
    let card = contact_card(UncertaintyModel::HalfWidth {
        half_width: 0.01,
        confidence: 0.95,
    });
    let resistance = resistance(&card);
    let card_identity = card.content_hash();
    let receipt_identity = resistance.receipt().content_hash();
    card.claims()
        .verify_receipt(resistance.receipt())
        .expect("retained property receipt replays");
    let interface = InterfaceSurface::new("bondline", oriented_pairs(&mesh), resistance.clone())
        .expect("interface surface");
    let interfaces = ThermalInterfaces::new(&mesh, &boundary, vec![interface])
        .expect("complete interface binding");
    assert_eq!(interfaces.surface_count(), 1);
    let solution = solve_contact_fixture(&mesh, &boundary, &interfaces);

    let mut max_error = 0.0f64;
    for (vertex, &point) in mesh.positions().iter().enumerate() {
        let exact = if vertex < left_vertex_count {
            T_HOT - 10.0 * point[0]
        } else {
            320.0 - 10.0 * point[0]
        };
        max_error = max_error.max((solution.temperature[vertex] - exact).abs());
    }
    assert!(
        max_error < 1e-8,
        "piecewise linear two-slab profile should be round-off exact; max error={max_error:e} K"
    );

    let flux = solution
        .report
        .interface_fluxes
        .first()
        .expect("one interface flux");
    assert_eq!(solution.report.interface_fluxes.len(), 1);
    assert_eq!(flux.interface, "bondline");
    assert!((flux.area_m2 - 1.0).abs() < 1e-14);
    assert!((flux.conductance_w_per_k - 10.0).abs() < 1e-12);
    assert!((flux.mean_jump_k - 10.0).abs() < 1e-8);
    assert!((flux.heat_rate_a_to_b_w - 100.0).abs() < 1e-7);
    assert_eq!(flux.card_identity, card_identity);
    assert_eq!(flux.receipt.content_hash(), receipt_identity);
    card.claims()
        .verify_receipt(&flux.receipt)
        .expect("reported receipt replays");
    let (level_a_case, reference) = level_a_contact_reference();
    let solved_resistance = (T_HOT - T_COLD) / flux.heat_rate_a_to_b_w;
    let (level_a_error, level_a_limit) =
        level_a_tolerance_error_and_limit(level_a_case, solved_resistance);
    assert!(
        level_a_error <= level_a_limit,
        "solve-derived resistance {solved_resistance:.16e} K/W differs from the \
         Level-A reference {reference:.16e} K/W by {level_a_error:.16e}, above \
         the registry limit {level_a_limit:.16e}"
    );

    let slab_a = ThermalResistanceTerm::slab(
        "slab-a",
        1.0,
        K,
        1.0,
        ResistanceUncertainty::HalfWidth {
            half_width: 0.0,
            confidence: 0.99,
        },
        "exact G1 fixture input",
    )
    .expect("slab-a term");
    let slab_b = ThermalResistanceTerm::slab(
        "slab-b",
        1.0,
        K,
        1.0,
        ResistanceUncertainty::HalfWidth {
            half_width: 0.0,
            confidence: 0.99,
        },
        "exact G1 fixture input",
    )
    .expect("slab-b term");
    let contact = resistance
        .term_for_area("bondline", flux.area_m2)
        .expect("contact term");
    let series =
        SeriesThermalResistance::new(vec![slab_b.clone(), contact.clone(), slab_a.clone()])
            .expect("series budget");
    assert!(
        (series.budget().value_k_per_w - reference).abs() <= f64::EPSILON,
        "the auxiliary three-term arithmetic budget must stay within one \
         absolute f64::EPSILON of the Level-A decimal reference"
    );
    assert_eq!(series.budget().complete_half_width_k_per_w(), Some(0.01));
    assert_eq!(series.budget().confidence_floor, Some(0.95));
    assert!(series.budget().unbounded_terms.is_empty());

    let permuted = SeriesThermalResistance::new(vec![contact, slab_a, slab_b])
        .expect("permuted series budget");
    assert_eq!(series, permuted, "input permutation must not move a bit");

    println!(
        "{{\"suite\":\"fs-conduction/contact\",\"case\":\"two-slab-series\",\
         \"level_a_case_id\":\"thermal-a-contact-series\",\"verdict\":\"pass\",\
         \"authority\":\"executed-matching-p1-interface-not-retained-registry-receipt\",\
         \"detail\":\"R_solved={solved_resistance:.16e} K/W; \
         R_budget={:.16e} K/W; max_nodal_error={max_error:.16e} K; \
         interface_heat_rate={:.16e} W; card_identity={card_identity}; receipt_identity={receipt_identity}\"}}",
        series.budget().value_k_per_w,
        flux.heat_rate_a_to_b_w,
    );
}

#[test]
fn missing_interface_card_and_missing_binding_refuse_typed() {
    let missing = InterfaceResistance::from_card(
        "bondline",
        &empty_contact_card(),
        &QueryPoint::new(),
        SelectionPolicy::SingleClaimOnly,
    )
    .expect_err("a missing contact-resistance claim must refuse");
    assert!(matches!(
        missing,
        ConductionError::Interface { ref interface, .. } if interface == "bondline"
    ));

    let (mesh, _) = two_slab_mesh(1);
    let boundary = boundary(&mesh);
    let material = ConductivityModel::isotropic_declared(K).expect("material");
    let source = ScalarField::Uniform(0.0);
    let error = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem { element_materials: None, mesh: &mesh,
            boundary: &boundary,
            material: &material,
            source: &source, },
            config(),
        )
        .expect_err("coincident traces without an interface binding must refuse")
    });
    assert!(matches!(
        error,
        ConductionError::Interface { ref interface, .. } if interface == "<undeclared>"
    ));
    assert_eq!(error.rule(), "conduction-interface");
}

#[test]
fn surface_declaration_and_binding_guards_each_refuse_their_own_defect() {
    let card = contact_card(UncertaintyModel::Unstated);
    let declared = resistance(&card);
    let (mesh, _) = two_slab_mesh(1);
    let boundary = boundary(&mesh);
    let pairs = oriented_pairs(&mesh);
    let pair = pairs[0];

    assert_interface_refusal(
        InterfaceSurface::new("   ", vec![pair], declared.clone())
            .expect_err("a blank interface name must refuse"),
        "<unnamed>",
        "name is blank",
    );
    assert_interface_refusal(
        InterfaceSurface::new("empty", Vec::new(), declared.clone())
            .expect_err("an empty interface surface must refuse"),
        "empty",
        "has no paired faces",
    );

    let duplicate_name = ThermalInterfaces::new(
        &mesh,
        &boundary,
        vec![
            InterfaceSurface::new("duplicate", vec![pair], declared.clone())
                .expect("first declaration"),
            InterfaceSurface::new("duplicate", vec![pair], declared.clone())
                .expect("second declaration"),
        ],
    )
    .expect_err("one semantic interface name cannot identify two surfaces");
    assert_interface_refusal(duplicate_name, "duplicate", "surface name is duplicated");

    let duplicate_pair = ThermalInterfaces::new(
        &mesh,
        &boundary,
        vec![
            InterfaceSurface::new("alpha", vec![pair], declared.clone())
                .expect("first pair declaration"),
            InterfaceSurface::new("beta", vec![pair], declared.clone())
                .expect("second pair declaration"),
        ],
    )
    .expect_err("one geometric face pair cannot be bound twice");
    assert_interface_refusal(duplicate_pair, "beta", "was bound more than once");

    let outside = ThermalInterfaces::new(
        &mesh,
        &boundary,
        vec![
            InterfaceSurface::new(
                "outside",
                vec![InterfaceFacePair {
                    side_a: mesh.boundary().len(),
                    side_b: pair.side_b,
                }],
                declared.clone(),
            )
            .expect("range validation belongs to mesh binding"),
        ],
    )
    .expect_err("an out-of-range boundary slot must refuse");
    assert_interface_refusal(outside, "outside", "is outside the");

    let self_pair = ThermalInterfaces::new(
        &mesh,
        &boundary,
        vec![
            InterfaceSurface::new(
                "self-pair",
                vec![InterfaceFacePair {
                    side_a: pair.side_a,
                    side_b: pair.side_a,
                }],
                declared.clone(),
            )
            .expect("self-pair validation belongs to mesh binding"),
        ],
    )
    .expect_err("a face cannot transfer heat to itself");
    assert_interface_refusal(self_pair, "self-pair", "is paired with itself");

    let hot_face = mesh
        .boundary()
        .iter()
        .position(|face| on_box_face(face.centroid[0], 0.0))
        .expect("the slab has a hot exterior face");
    let noncoincident = ThermalInterfaces::new(
        &mesh,
        &boundary,
        vec![
            InterfaceSurface::new(
                "noncoincident",
                vec![InterfaceFacePair {
                    side_a: pair.side_a,
                    side_b: hot_face,
                }],
                declared.clone(),
            )
            .expect("coincidence validation belongs to mesh binding"),
        ],
    )
    .expect_err("unrelated exterior faces cannot form a contact pair");
    assert_interface_refusal(noncoincident, "noncoincident", "is not exactly coincident");

    let occupied = ThermalInterfaces::new(
        &mesh,
        &boundary_with_occupied_interface(&mesh),
        vec![
            InterfaceSurface::new("occupied", pairs, declared)
                .expect("complete occupied interface declaration"),
        ],
    )
    .expect_err("an interface face cannot also carry an external boundary row");
    assert_interface_refusal(
        occupied,
        "occupied",
        "also carries an external boundary condition",
    );
}

#[test]
fn coincident_geometry_refuses_nonopposing_and_multi_owner_traces() {
    let same_side = coincident_base_tet_mesh(&[1.0, 2.0]);
    assert_interface_refusal(
        ThermalInterfaces::coincident_face_pairs(&same_side)
            .expect_err("same-side coincident traces have nonopposing normals"),
        "<geometry>",
        "do not have opposing normals",
    );

    let three_owner = coincident_base_tet_mesh(&[-1.0, 1.0, 2.0]);
    assert_interface_refusal(
        ThermalInterfaces::coincident_face_pairs(&three_owner)
            .expect_err("three boundary traces cannot own one geometric triangle"),
        "<geometry>",
        "3 boundary faces share one exact geometric triangle",
    );
}

#[test]
fn unstated_uncertainty_stays_incomplete() {
    let card = contact_card(UncertaintyModel::Unstated);
    let term = resistance(&card)
        .term_for_area("bondline", 1.0)
        .expect("contact term");
    let series = SeriesThermalResistance::new(vec![term]).expect("series budget");
    assert_eq!(series.budget().stated_half_width_k_per_w, 0.0);
    assert_eq!(series.budget().confidence_floor, None);
    assert_eq!(series.budget().unbounded_terms, ["bondline"]);
    assert_eq!(series.budget().complete_half_width_k_per_w(), None);
}

#[test]
fn series_uncertainty_overflow_refuses() {
    let term = |name| {
        ThermalResistanceTerm::declared(
            name,
            1.0,
            ResistanceUncertainty::HalfWidth {
                half_width: f64::MAX,
                confidence: 0.95,
            },
            "overflow refusal fixture",
        )
        .expect("individual term remains finite")
    };
    let error = SeriesThermalResistance::new(vec![term("a"), term("b")])
        .expect_err("the aggregate uncertainty band must remain finite");
    assert!(matches!(error, ConductionError::Interface { .. }));
    assert_eq!(error.rule(), "conduction-interface");
}

#[test]
fn a_bound_interface_refuses_a_different_face_set() {
    let card = contact_card(UncertaintyModel::Unstated);
    let (coarse_mesh, _) = two_slab_mesh(1);
    let coarse_boundary = boundary(&coarse_mesh);
    let interface =
        InterfaceSurface::new("bondline", oriented_pairs(&coarse_mesh), resistance(&card))
            .expect("coarse interface surface");
    let interfaces = ThermalInterfaces::new(&coarse_mesh, &coarse_boundary, vec![interface])
        .expect("coarse interface binding");

    let (fine_mesh, _) = two_slab_mesh(2);
    let fine_boundary = boundary(&fine_mesh);
    assert_interface_refusal(
        binding_validation_error(&fine_mesh, &fine_boundary, &interfaces),
        "<binding>",
        "face set does not match",
    );
}

#[test]
fn a_bound_interface_refuses_a_changed_boundary_condition() {
    let card = contact_card(UncertaintyModel::Unstated);
    let (mesh, _) = two_slab_mesh(1);
    let interfaces = ThermalInterfaces::new(
        &mesh,
        &boundary(&mesh),
        vec![
            InterfaceSurface::new("bondline", oriented_pairs(&mesh), resistance(&card))
                .expect("interface surface"),
        ],
    )
    .expect("interface binding");

    assert_interface_refusal(
        binding_validation_error(&mesh, &boundary_with_occupied_interface(&mesh), &interfaces),
        "bondline",
        "now carries an external boundary condition",
    );
}

#[test]
fn a_bound_interface_refuses_a_changed_exact_vertex_correspondence() {
    let card = contact_card(UncertaintyModel::Unstated);
    let mesh = two_tet_contact_mesh();
    let original_pairs =
        ThermalInterfaces::coincident_face_pairs(&mesh).expect("one coincident face pair");
    assert_eq!(original_pairs.len(), 1);
    let pair = original_pairs[0];
    let interfaces = ThermalInterfaces::new(
        &mesh,
        &adiabatic_boundary(&mesh),
        vec![
            InterfaceSurface::new("bondline", original_pairs.clone(), resistance(&card))
                .expect("interface surface"),
        ],
    )
    .expect("interface binding");

    let mut changed_positions = mesh.positions().to_vec();
    let side_b_vertices = mesh.boundary()[pair.side_b]
        .vertices
        .map(|vertex| vertex as usize);
    changed_positions.swap(side_b_vertices[0], side_b_vertices[1]);
    let changed_mesh = ConductionMesh::new(mesh.complex().clone(), changed_positions)
        .expect("swapping two face coordinates preserves nondegenerate tetrahedra");
    assert_eq!(
        ThermalInterfaces::coincident_face_pairs(&changed_mesh).expect("same face-pair slots"),
        original_pairs,
        "the fixture must isolate vertex correspondence rather than face-set drift"
    );
    assert_eq!(
        changed_mesh.boundary()[pair.side_a].area.to_bits(),
        mesh.boundary()[pair.side_a].area.to_bits(),
        "the fixture must isolate correspondence rather than area drift"
    );

    assert_interface_refusal(
        binding_validation_error(
            &changed_mesh,
            &adiabatic_boundary(&changed_mesh),
            &interfaces,
        ),
        "bondline",
        "changed after binding",
    );
}

// ---------------------------------------------------------------------------
// Why `matching_vertices`' non-unique-correspondence refusal has no test.
//
// interface.rs's `matching_vertices` refuses with "does not have a unique exact
// vertex correspondence" when a side-A vertex has anything other than exactly
// one bit-equal partner in side B. That branch is UNREACHABLE from every public
// entry point, and the two tests below are the reachability argument in code
// rather than a test of the branch itself:
//
//   * matches.len() == 0 is impossible because `ThermalInterfaces::new` refuses
//     any declared pair that is not in `coincident_candidates`, and candidacy is
//     equality of the SORTED triple of coordinate bit patterns -- so every
//     side-A key occurs in side B with the same multiplicity.
//   * matches.len() >= 2 needs two bit-identical vertices inside one boundary
//     triangle, which forces its parent tet to zero signed volume, which
//     `ConductionMesh::new` refuses before any interface code runs. That is
//     `a_face_with_two_bit_identical_vertices_cannot_be_meshed`.
//
// `ConductionMesh` exposes no post-construction mutator, so the invariant cannot
// be re-broken after `new` admits a mesh. The branch is therefore defensive
// depth, and no document may cite it as a covered refusal.
//
// This is not a paper argument. Deleting the candidacy guard from
// `ThermalInterfaces::new` and rerunning makes
// `the_candidacy_guard_is_what_keeps_vertex_matching_unreachable` fail by
// reporting interface `"<geometry>"` -- the literal name
// `matching_vertices` refuses under -- so the fall-through is observed, and the
// branch is live defensive code rather than something that could be deleted.
//
// Two corrections to the record, both against my own earlier report:
//   * the candidacy guard was ALREADY covered, by the "is not exactly
//     coincident" clause of
//     `surface_declaration_and_binding_guards_each_refuse_their_own_defect`.
//     A grep for the full message string missed it because that test asserts a
//     substring. The test below is kept for a different reason: under the guard
//     deletion above it is the fixture that reaches `matching_vertices`, while
//     the older one stops at the occupied-boundary guard.
//   * `a_bound_interface_refuses_a_changed_exact_vertex_correspondence` does NOT
//     reach `matching_vertices` despite its name; it perturbs an already-bound
//     mesh, so the "changed after binding" guard fires first.
// ---------------------------------------------------------------------------

#[test]
fn the_candidacy_guard_is_what_keeps_vertex_matching_unreachable() {
    let card = contact_card(UncertaintyModel::Unstated);
    let mesh = two_tet_contact_mesh();
    let pairs = ThermalInterfaces::coincident_face_pairs(&mesh).expect("one coincident face pair");
    assert_eq!(pairs.len(), 1);
    let pair = pairs[0];
    let unrelated = (0..mesh.boundary().len())
        .find(|slot| *slot != pair.side_a && *slot != pair.side_b)
        .expect("the fixture has boundary faces outside the contact plane");
    let bogus = InterfaceFacePair {
        side_a: pair.side_a,
        side_b: unrelated,
    };
    assert!(
        !pairs.iter().any(|candidate| {
            (candidate.side_a == bogus.side_a && candidate.side_b == bogus.side_b)
                || (candidate.side_a == bogus.side_b && candidate.side_b == bogus.side_a)
        }),
        "the fixture must declare a pair the geometry does not support"
    );

    let error = ThermalInterfaces::new(
        &mesh,
        &adiabatic_boundary(&mesh),
        vec![
            InterfaceSurface::new("bondline", vec![bogus], resistance(&card))
                .expect("interface surface"),
        ],
    )
    .expect_err("a non-coincident declared pair must refuse before vertex matching");
    assert_interface_refusal(
        error,
        "bondline",
        "not exactly coincident matching P1 geometry",
    );
}

#[test]
fn a_face_with_two_bit_identical_vertices_cannot_be_meshed() {
    // The last two vertices are bit-identical, so any boundary triangle using
    // both would give `matching_vertices` a two-way match. The tet's signed
    // volume is exactly zero, so the mesh never exists.
    let positions = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    let error = ConductionMesh::new(
        TetComplex::from_tets(positions.len(), vec![[0, 1, 2, 3]]),
        positions,
    )
    .expect_err("a tet with two coincident vertices is degenerate");
    match error {
        ConductionError::DegenerateElement {
            element,
            signed_volume,
        } => {
            assert_eq!(element, 0);
            assert!(
                signed_volume == 0.0,
                "two coincident vertices must give exactly zero volume, got {signed_volume:e}"
            );
        }
        other => panic!("expected a typed degenerate-element refusal, got {other:?}"),
    }
}

#[test]
fn a_bound_interface_refuses_changed_face_area_bits() {
    let card = contact_card(UncertaintyModel::Unstated);
    let mesh = two_tet_contact_mesh();
    let original_pairs =
        ThermalInterfaces::coincident_face_pairs(&mesh).expect("one coincident face pair");
    let pair = original_pairs[0];
    let interfaces = ThermalInterfaces::new(
        &mesh,
        &adiabatic_boundary(&mesh),
        vec![
            InterfaceSurface::new("bondline", original_pairs.clone(), resistance(&card))
                .expect("interface surface"),
        ],
    )
    .expect("interface binding");

    let changed_positions = mesh
        .positions()
        .iter()
        .map(|&[x, y, z]| [2.0 * x, y, z])
        .collect();
    let changed_mesh = ConductionMesh::new(mesh.complex().clone(), changed_positions)
        .expect("uniform coordinate scaling preserves nondegenerate tetrahedra");
    assert_eq!(
        ThermalInterfaces::coincident_face_pairs(&changed_mesh).expect("same face-pair slots"),
        original_pairs,
        "the fixture must isolate area rather than face-set drift"
    );
    assert_eq!(
        changed_mesh.boundary()[pair.side_a].vertices,
        mesh.boundary()[pair.side_a].vertices,
        "the fixture must retain the bound vertex ids"
    );
    assert_ne!(
        changed_mesh.boundary()[pair.side_a].area.to_bits(),
        mesh.boundary()[pair.side_a].area.to_bits(),
        "the fixture must actually move the stored area bits"
    );

    assert_interface_refusal(
        binding_validation_error(
            &changed_mesh,
            &adiabatic_boundary(&changed_mesh),
            &interfaces,
        ),
        "bondline",
        "changed after binding",
    );
}

/// Solve the two-slab contact fixture with a caller-chosen area-specific bond
/// line and return the series resistance the SOLVE implies, `ΔT / Q`.
///
/// This is the quantity the Level-A registry row is about, derived exactly the
/// way the accepting test derives it. A falsifier built on it therefore has to
/// travel the whole card -> resistance -> binding -> assembly -> solve -> report
/// chain, instead of doing arithmetic on registry constants.
fn solved_series_resistance(area_specific_resistance_m2_k_w: f64) -> f64 {
    let (mesh, _) = two_slab_mesh(4);
    let boundary = boundary(&mesh);
    let card = contact_card(UncertaintyModel::HalfWidth {
        half_width: 0.01,
        confidence: 0.95,
    });
    let mut resistance = resistance(&card);
    if area_specific_resistance_m2_k_w.to_bits() != AREA_SPECIFIC_RESISTANCE.to_bits() {
        resistance = resistance
            .with_measured_value(
                area_specific_resistance_m2_k_w,
                UncertaintyModel::Unstated,
                "deliberate bond-line perturbation",
            )
            .expect("perturbed bond-line resistance");
    }
    let interface = InterfaceSurface::new("bondline", oriented_pairs(&mesh), resistance)
        .expect("interface surface");
    let interfaces = ThermalInterfaces::new(&mesh, &boundary, vec![interface])
        .expect("complete interface binding");
    let solution = solve_contact_fixture(&mesh, &boundary, &interfaces);
    let flux = solution
        .report
        .interface_fluxes
        .first()
        .expect("one interface flux");
    (T_HOT - T_COLD) / flux.heat_rate_a_to_b_w
}

/// The applied Level-A envelope must reject a bond line that is actually wrong.
///
/// The predecessor of this test set `outside = reference + 2 * limit` and then
/// asserted `|outside - reference| > limit`. That is unconditionally true for
/// any positive limit -- which the helper already asserts -- so it exercised no
/// production code at all: it was arithmetic on two registry constants wearing
/// the name of a tolerance check.
#[test]
fn contact_level_a_tolerance_rejects_a_wrong_bond_line_from_the_actual_solve() {
    let (case, reference) = level_a_contact_reference();
    // Registry pin, kept deliberately and labelled for what it is: a
    // change-detector on the row's stored value, not a check of anything this
    // crate computes.
    assert_eq!(reference.to_bits(), 0.3f64.to_bits());

    // DISCRIMINATION, not merely rejection. A gate that rejected everything
    // would satisfy the falsifier below and look identical to one that works,
    // so the accepting case is asserted here too rather than assumed from the
    // neighbouring test.
    let accepted = solved_series_resistance(AREA_SPECIFIC_RESISTANCE);
    let (accepted_error, accepted_limit) = level_a_tolerance_error_and_limit(case, accepted);
    assert!(
        accepted_error <= accepted_limit,
        "the card's own bond line must stay inside the Level-A envelope: solved \
         {accepted:.16e} K/W, error {accepted_error:.3e}, limit {accepted_limit:.3e}"
    );

    // HEADROOM, measured rather than assumed. The accepting case sits at
    // ~1.39e-14 against the row's 1.3e-12 limit, about 94x clear. Pin one order
    // of magnitude of that margin so a future solver or `config()` change that
    // eats it is caught while margin remains, instead of surfacing later as an
    // intermittent failure sitting on the registry bound.
    assert!(
        accepted_error * 10.0 < accepted_limit,
        "the accepting case must keep an order of magnitude of headroom below the \
         registry limit: error {accepted_error:.3e}, limit {accepted_limit:.3e}"
    );

    // The falsifier is a PLAUSIBLE modelling mistake -- a 1% bond-line error --
    // rather than a value constructed two tolerances away from the answer.
    const BOND_LINE_ERROR: f64 = 0.01;
    let rejected = solved_series_resistance(AREA_SPECIFIC_RESISTANCE * (1.0 + BOND_LINE_ERROR));
    let (rejected_error, rejected_limit) = level_a_tolerance_error_and_limit(case, rejected);
    assert!(
        rejected_error > rejected_limit,
        "a {BOND_LINE_ERROR} relative bond-line error must be rejected: solved \
         {rejected:.16e} K/W, error {rejected_error:.3e}, limit {rejected_limit:.3e}"
    );

    // And the perturbation must land where the physics says it does. The bond
    // line is one of three equal series terms over a 1 m² trace, so a relative
    // error in R'' moves the total by exactly `R'' * BOND_LINE_ERROR`. Without
    // this, the rejection above would also pass if the perturbation had
    // corrupted some unrelated term by an arbitrary amount.
    let observed_shift = rejected - accepted;
    let expected_shift = AREA_SPECIFIC_RESISTANCE * BOND_LINE_ERROR;
    assert!(
        (observed_shift - expected_shift).abs() < 1.0e-9,
        "the perturbation must move the series resistance by the bond-line delta: \
         expected {expected_shift:.16e} K/W, observed {observed_shift:.16e} K/W"
    );
}

// ---------------------------------------------------------------------------
// Map-valued (spatially varying) interface resistance — bead f85xj.12.9.
// ---------------------------------------------------------------------------

/// Build one per-face resistance list from the shared card, scaling the card
/// value by a per-face factor. This is the shape an as-built thickness map
/// takes once the caller has turned `t_i` into `R''_i = t_i / k`: the card
/// supplies `k` and the material authority, metrology supplies `t_i`.
fn scaled_face_resistances(
    card: &InterfaceSystemCard,
    factors: &[f64],
) -> Vec<InterfaceResistance> {
    let base = resistance(card);
    factors
        .iter()
        .map(|factor| {
            base.with_measured_value(
                AREA_SPECIFIC_RESISTANCE * factor,
                UncertaintyModel::Unstated,
                "synthetic as-built bond-line thickness map",
            )
            .expect("measured face resistance")
        })
        .collect()
}

#[test]
fn a_mapped_interface_applies_each_face_its_own_resistance() {
    let card = contact_card(UncertaintyModel::Unstated);
    let (mesh, _) = two_slab_mesh(2);
    let pairs = oriented_pairs(&mesh);
    // Every face gets a distinct multiplier so a single shared value, or a
    // silently reused first entry, cannot reproduce the expected conductance.
    let factors: Vec<f64> = (0..pairs.len())
        .map(|index| 1.0 + 0.5 * (index as f64))
        .collect();
    let surface = InterfaceSurface::new_mapped(
        "bondline",
        pairs.clone(),
        scaled_face_resistances(&card, &factors),
    )
    .expect("mapped interface admits");
    assert!(surface.is_mapped());
    assert!(
        surface.resistance().is_none(),
        "a mapped surface has no single resistance to report"
    );

    let boundary = boundary(&mesh);
    let interfaces =
        ThermalInterfaces::new(&mesh, &boundary, vec![surface]).expect("interfaces bind");
    assert_eq!(interfaces.surface_is_mapped("bondline"), Some(true));

    // The bound per-face values must be a permutation of what was supplied,
    // with the same multiset of values.
    let mut bound = interfaces
        .surface_face_resistances("bondline")
        .expect("bound resistances");
    let mut expected: Vec<f64> = factors
        .iter()
        .map(|factor| AREA_SPECIFIC_RESISTANCE * factor)
        .collect();
    bound.sort_by(f64::total_cmp);
    expected.sort_by(f64::total_cmp);
    assert_eq!(bound.len(), expected.len());
    for (got, want) in bound.iter().zip(expected.iter()) {
        assert!((got - want).abs() < 1e-15, "bound {got} != supplied {want}");
    }
}

#[test]
fn mapped_resistances_follow_their_face_pairs_through_the_deterministic_sort() {
    // THE hazard this bead exists to prevent. `ThermalInterfaces::new` sorts
    // face pairs into a deterministic order. If the resistance map were a
    // parallel array indexed after that sort, every face would silently get a
    // neighbour's value and nothing would report an error. Binding a
    // deliberately REVERSED declaration must produce exactly the same solve as
    // declaring the same (pair, resistance) association in sorted order.
    let card = contact_card(UncertaintyModel::Unstated);
    let (mesh, _) = two_slab_mesh(2);
    let boundary = boundary(&mesh);
    let pairs = oriented_pairs(&mesh);
    let factors: Vec<f64> = (0..pairs.len())
        .map(|index| 1.0 + 0.25 * (index as f64))
        .collect();

    let forward = InterfaceSurface::new_mapped(
        "bondline",
        pairs.clone(),
        scaled_face_resistances(&card, &factors),
    )
    .expect("forward mapped interface");

    // Reverse BOTH the pairs and their factors, preserving each pair's
    // association with its own value. The declaration order changes; the
    // physics must not.
    let mut reversed_pairs = pairs.clone();
    reversed_pairs.reverse();
    let mut reversed_factors = factors.clone();
    reversed_factors.reverse();
    let reversed = InterfaceSurface::new_mapped(
        "bondline",
        reversed_pairs,
        scaled_face_resistances(&card, &reversed_factors),
    )
    .expect("reversed mapped interface");

    let forward_bound = ThermalInterfaces::new(&mesh, &boundary, vec![forward])
        .expect("forward binds")
        .surface_face_resistances("bondline")
        .expect("forward resistances");
    let reversed_bound = ThermalInterfaces::new(&mesh, &boundary, vec![reversed])
        .expect("reversed binds")
        .surface_face_resistances("bondline")
        .expect("reversed resistances");

    assert_eq!(
        forward_bound.len(),
        reversed_bound.len(),
        "both declarations bind the same faces"
    );
    for (index, (forward_value, reversed_value)) in
        forward_bound.iter().zip(reversed_bound.iter()).enumerate()
    {
        assert!(
            (forward_value - reversed_value).abs() < 1e-15,
            "bound face {index} got {forward_value} from the forward declaration \
             and {reversed_value} from the reversed one; the sort separated a \
             face from its resistance"
        );
    }
    // And the values must not have collapsed to one number, or the test would
    // pass vacuously.
    let first = forward_bound[0];
    assert!(
        forward_bound
            .iter()
            .any(|value| (value - first).abs() > 1e-12),
        "the fixture must supply genuinely distinct per-face resistances"
    );
}

#[test]
fn two_named_surfaces_sort_deterministically_without_swapping_their_physics() {
    let card = contact_card(UncertaintyModel::Unstated);
    let base = resistance(&card);
    let higher_resistance = base
        .with_measured_value(
            2.0 * AREA_SPECIFIC_RESISTANCE,
            UncertaintyModel::Unstated,
            "second declared bond line",
        )
        .expect("second resistance");
    let (mesh, vertices_per_slab) = slab_chain_mesh(1, 3);
    let boundary = outer_boundary(&mesh, 3.0);
    let first_plane = oriented_pairs_at_x(&mesh, 1.0);
    let second_plane = oriented_pairs_at_x(&mesh, 2.0);
    assert!(
        !first_plane.is_empty() && !second_plane.is_empty(),
        "the fixture must contain two nonempty contact surfaces"
    );
    assert_eq!(
        first_plane.len() + second_plane.len(),
        oriented_pairs(&mesh).len(),
        "the two named surfaces must partition every coincident pair"
    );

    let bind = |reverse_declaration: bool| {
        let alpha = InterfaceSurface::new(
            "alpha-second-plane",
            second_plane.clone(),
            higher_resistance.clone(),
        )
        .expect("alpha surface");
        let zeta = InterfaceSurface::new("zeta-first-plane", first_plane.clone(), base.clone())
            .expect("zeta surface");
        let surfaces = if reverse_declaration {
            vec![zeta, alpha]
        } else {
            vec![alpha, zeta]
        };
        ThermalInterfaces::new(&mesh, &boundary, surfaces).expect("two-surface binding")
    };
    let name_order = bind(false);
    let reverse_order = bind(true);

    let mut temperature = vec![0.0; mesh.vertex_count()];
    for (vertex, value) in temperature.iter_mut().enumerate() {
        *value = match vertex / vertices_per_slab {
            0 => T_HOT,
            1 => 0.5 * (T_HOT + T_COLD),
            2 => T_COLD,
            other => panic!("unexpected slab index {other}"),
        };
    }
    let name_fluxes = name_order
        .fluxes(&temperature)
        .expect("name-order interface fluxes");
    let reverse_fluxes = reverse_order
        .fluxes(&temperature)
        .expect("reverse-order interface fluxes");
    assert_eq!(
        name_fluxes
            .iter()
            .map(|flux| flux.interface.as_str())
            .collect::<Vec<_>>(),
        ["alpha-second-plane", "zeta-first-plane"],
        "reported surfaces must follow stable name order"
    );
    assert_eq!(name_fluxes.len(), reverse_fluxes.len());
    for (name_flux, reverse_flux) in name_fluxes.iter().zip(&reverse_fluxes) {
        assert_eq!(name_flux.interface, reverse_flux.interface);
        assert_eq!(
            name_flux.conductance_w_per_k.to_bits(),
            reverse_flux.conductance_w_per_k.to_bits()
        );
        assert_eq!(
            name_flux.mean_jump_k.to_bits(),
            reverse_flux.mean_jump_k.to_bits()
        );
        assert_eq!(
            name_flux.heat_rate_a_to_b_w.to_bits(),
            reverse_flux.heat_rate_a_to_b_w.to_bits()
        );
    }
    assert!(
        (name_fluxes[0].heat_rate_a_to_b_w - name_fluxes[1].heat_rate_a_to_b_w).abs() > 1.0,
        "distinct resistances must make a surface swap observable"
    );
}

/// Drive `assemble_into` with MORE THAN ONE bound surface, against a closed
/// form, with a per-surface oracle.
///
/// The deterministic-sort test above binds two surfaces but stops at
/// `fluxes()`, so before this test the multi-surface assembly path was never
/// executed even though the contract talks about ordering stability "before
/// assembly or summation". Two unequal bond lines in series also give the
/// fixture a per-surface oracle: a total-heat check alone cannot see the two
/// surfaces swapped, because series conservation puts the same Q through both.
#[test]
fn two_named_surfaces_assemble_into_one_series_solve() {
    let card = contact_card(UncertaintyModel::Unstated);
    let base = resistance(&card);
    let doubled = base
        .with_measured_value(
            2.0 * AREA_SPECIFIC_RESISTANCE,
            UncertaintyModel::Unstated,
            "second declared bond line",
        )
        .expect("second resistance");
    let (mesh, _) = slab_chain_mesh(1, 3);
    let boundary = outer_boundary(&mesh, 3.0);
    let first_plane = oriented_pairs_at_x(&mesh, 1.0);
    let second_plane = oriented_pairs_at_x(&mesh, 2.0);
    assert!(
        !first_plane.is_empty() && !second_plane.is_empty(),
        "the fixture must contain two nonempty contact surfaces"
    );
    assert_eq!(
        first_plane.len() + second_plane.len(),
        oriented_pairs(&mesh).len(),
        "the two named surfaces must partition every coincident pair"
    );
    let interfaces = ThermalInterfaces::new(
        &mesh,
        &boundary,
        vec![
            InterfaceSurface::new("bond-thin", first_plane, base).expect("thin surface"),
            InterfaceSurface::new("bond-thick", second_plane, doubled).expect("thick surface"),
        ],
    )
    .expect("two-surface binding");
    assert_eq!(interfaces.surface_count(), 2);

    let solution = solve_contact_fixture(&mesh, &boundary, &interfaces);
    assert_eq!(
        solution.report.interface_fluxes.len(),
        2,
        "an assembled two-surface solve must report both interfaces"
    );
    let flux_for = |name: &str| {
        solution
            .report
            .interface_fluxes
            .iter()
            .find(|flux| flux.interface == name)
            .unwrap_or_else(|| panic!("missing reported interface {name}"))
    };
    let thin = flux_for("bond-thin");
    let thick = flux_for("bond-thick");

    // Closed form: three unit slabs at k = 10 W/m·K over a 1 m² trace give
    // 0.1 K/W each, and the two bond lines add 0.1 and 0.2 K/W, so
    // R_total = 0.6 K/W and Q = 30 K / 0.6 K/W = 50 W.
    const EXPECTED_HEAT_W: f64 = 50.0;
    for flux in [thin, thick] {
        assert!(
            (flux.heat_rate_a_to_b_w - EXPECTED_HEAT_W).abs() < 1e-7,
            "series conservation puts {EXPECTED_HEAT_W} W through {}, got {} W",
            flux.interface,
            flux.heat_rate_a_to_b_w
        );
    }

    // PER-SURFACE oracle. The jumps are Q·R_c, so they differ by construction
    // and identify which surface got which card. Swapping `base` and `doubled`
    // above leaves every heat rate untouched and fails exactly here.
    assert!(
        (thin.mean_jump_k - EXPECTED_HEAT_W * AREA_SPECIFIC_RESISTANCE).abs() < 1e-7,
        "the thin bond line must drop Q·R'' = 5 K, got {} K",
        thin.mean_jump_k
    );
    assert!(
        (thick.mean_jump_k - EXPECTED_HEAT_W * 2.0 * AREA_SPECIFIC_RESISTANCE).abs() < 1e-7,
        "the thick bond line must drop Q·2R'' = 10 K, got {} K",
        thick.mean_jump_k
    );
    assert!(
        (thick.mean_jump_k - thin.mean_jump_k).abs() > 1.0,
        "the two jumps must be far apart or the swap above is not observable: \
         {} K vs {} K",
        thin.mean_jump_k,
        thick.mean_jump_k
    );
}

#[test]
fn an_all_equal_map_reproduces_the_uniform_surface_exactly() {
    // Metamorphic: a map whose every entry equals the card value is the same
    // physics as the uniform declaration. Run BOTH solves and require their
    // jump-bearing fields and interface reports to agree bitwise; comparing
    // conductance at a uniform temperature would exercise input arithmetic
    // while leaving the jump operator completely unconstrained.
    let card = contact_card(UncertaintyModel::Unstated);
    let (mesh, _) = two_slab_mesh(2);
    let boundary = boundary(&mesh);
    let pairs = oriented_pairs(&mesh);
    let ones = vec![1.0f64; pairs.len()];

    let uniform = ThermalInterfaces::new(
        &mesh,
        &boundary,
        vec![InterfaceSurface::new("bondline", pairs.clone(), resistance(&card)).expect("uniform")],
    )
    .expect("uniform binds");
    let mapped = ThermalInterfaces::new(
        &mesh,
        &boundary,
        vec![
            InterfaceSurface::new_mapped("bondline", pairs, scaled_face_resistances(&card, &ones))
                .expect("mapped"),
        ],
    )
    .expect("mapped binds");

    let uniform_solution = solve_contact_fixture(&mesh, &boundary, &uniform);
    let mapped_solution = solve_contact_fixture(&mesh, &boundary, &mapped);
    assert_eq!(
        uniform_solution.temperature.len(),
        mapped_solution.temperature.len()
    );
    for (vertex, (&uniform_temperature, &mapped_temperature)) in uniform_solution
        .temperature
        .iter()
        .zip(&mapped_solution.temperature)
        .enumerate()
    {
        assert_eq!(
            uniform_temperature.to_bits(),
            mapped_temperature.to_bits(),
            "all-equal mapped and uniform solves differ at vertex {vertex}"
        );
    }
    let uniform_flux = &uniform_solution.report.interface_fluxes[0];
    let mapped_flux = &mapped_solution.report.interface_fluxes[0];
    assert_eq!(
        uniform_flux.conductance_w_per_k.to_bits(),
        mapped_flux.conductance_w_per_k.to_bits()
    );
    assert_eq!(
        uniform_flux.mean_jump_k.to_bits(),
        mapped_flux.mean_jump_k.to_bits(),
        "the executed jump must be bitwise identical"
    );
    assert_eq!(
        uniform_flux.heat_rate_a_to_b_w.to_bits(),
        mapped_flux.heat_rate_a_to_b_w.to_bits(),
        "the executed interface heat rate must be bitwise identical"
    );
    assert!(
        uniform_flux.mean_jump_k.abs() > 1.0 && uniform_flux.heat_rate_a_to_b_w.abs() > 1.0,
        "the comparison must remain jump-bearing and non-vacuous"
    );
    assert_eq!(
        uniform_flux.card_identity, mapped_flux.card_identity,
        "both cite the same material card"
    );
    assert_eq!(uniform.surface_is_mapped("bondline"), Some(false));
    assert_eq!(mapped.surface_is_mapped("bondline"), Some(true));
}

#[test]
fn a_thicker_bond_line_conducts_less() {
    // Directional solve: doubling every face's area-specific resistance —
    // what a uniformly thicker measured bond line means — lowers the
    // SERIES heat rate and increases the interface jump. The independent
    // three-term closed forms are 100 W at R''=0.1 and 75 W at R''=0.2.
    let card = contact_card(UncertaintyModel::Unstated);
    let (mesh, _) = two_slab_mesh(2);
    let boundary = boundary(&mesh);
    let pairs = oriented_pairs(&mesh);
    let nominal_interfaces = ThermalInterfaces::new(
        &mesh,
        &boundary,
        vec![
            InterfaceSurface::new_mapped(
                "bondline",
                pairs.clone(),
                scaled_face_resistances(&card, &vec![1.0; pairs.len()]),
            )
            .expect("nominal"),
        ],
    )
    .expect("nominal binds");

    let thick_interfaces = ThermalInterfaces::new(
        &mesh,
        &boundary,
        vec![
            InterfaceSurface::new_mapped(
                "bondline",
                pairs.clone(),
                scaled_face_resistances(&card, &vec![2.0; pairs.len()]),
            )
            .expect("thick"),
        ],
    )
    .expect("thick binds");

    let nominal_solution = solve_contact_fixture(&mesh, &boundary, &nominal_interfaces);
    let thick_solution = solve_contact_fixture(&mesh, &boundary, &thick_interfaces);
    let nominal = &nominal_solution.report.interface_fluxes[0];
    let thick = &thick_solution.report.interface_fluxes[0];

    assert!(
        (nominal.heat_rate_a_to_b_w - 100.0).abs() < 1e-7,
        "the nominal solve must reproduce the independent 100 W series value"
    );
    assert!(
        (thick.heat_rate_a_to_b_w - 75.0).abs() < 1e-7,
        "the thicker solve must reproduce the independent 75 W series value"
    );
    assert!(
        thick.heat_rate_a_to_b_w < nominal.heat_rate_a_to_b_w,
        "a thicker bond line must transfer less heat under the same end temperatures"
    );
    assert!(
        thick.mean_jump_k > nominal.mean_jump_k,
        "the added contact resistance must move more of the fixed temperature drop onto the interface"
    );
    assert!(
        (thick.conductance_w_per_k - 0.5 * nominal.conductance_w_per_k).abs()
            < 1e-12 * nominal.conductance_w_per_k,
        "doubling R'' must halve the input-derived interface conductance"
    );
}

#[test]
fn a_mapped_interface_refuses_a_count_mismatch() {
    let card = contact_card(UncertaintyModel::Unstated);
    let (mesh, _) = two_slab_mesh(2);
    let pairs = oriented_pairs(&mesh);
    assert!(pairs.len() > 1, "the fixture needs several faces");
    let error = InterfaceSurface::new_mapped(
        "bondline",
        pairs.clone(),
        scaled_face_resistances(&card, &[1.0]),
    )
    .expect_err("one resistance cannot cover many faces");
    assert!(matches!(
        error,
        ConductionError::Interface { ref interface, .. } if interface == "bondline"
    ));
}

#[test]
fn a_mapped_interface_refuses_faces_citing_different_cards() {
    // A resistance map varies the VALUE, not the AUTHORITY. Two different
    // cards under one interface name is two interfaces, and it must be
    // declared as two named surfaces rather than silently averaged or
    // silently attributed to whichever card happened to be first.
    let card = contact_card(UncertaintyModel::Unstated);
    let other = contact_card(UncertaintyModel::HalfWidth {
        half_width: 0.01,
        confidence: 0.95,
    });
    let (mesh, _) = two_slab_mesh(2);
    let pairs = oriented_pairs(&mesh);
    assert!(pairs.len() > 1, "the fixture needs several faces");

    let mut mixed = scaled_face_resistances(&card, &vec![1.0; pairs.len()]);
    mixed[1] = resistance(&other);
    assert_ne!(
        mixed[0].card_identity(),
        mixed[1].card_identity(),
        "the fixture must actually supply two distinct cards"
    );

    let error = InterfaceSurface::new_mapped("bondline", pairs, mixed)
        .expect_err("a mapped surface must cite one card");
    assert!(matches!(
        error,
        ConductionError::Interface { ref interface, .. } if interface == "bondline"
    ));
}

#[test]
fn a_measured_value_cannot_masquerade_as_the_cards_own_claim() {
    // The card identity records the MATERIAL authority. It must not be
    // readable as an endorsement of a number the card never claimed: a
    // measured bond line is `t / k` with `t` from metrology. The origin is
    // what keeps citing a card distinct from laundering a measurement
    // through one.
    let card = contact_card(UncertaintyModel::Unstated);
    let base = resistance(&card);
    assert_eq!(base.value_origin(), &ResistanceValueOrigin::CardClaim);

    let measured = base
        .with_measured_value(
            0.25,
            UncertaintyModel::HalfWidth {
                half_width: 0.02,
                confidence: 0.95,
            },
            "as-built bond-line thickness map, station 7",
        )
        .expect("measured value admits");

    assert_eq!(measured.card_identity(), base.card_identity());
    assert!(
        (measured.value_m2_k_per_w() - 0.25).abs() < 1e-15,
        "the supplied value is retained verbatim"
    );
    match measured.value_origin() {
        ResistanceValueOrigin::CallerSupplied { rationale } => {
            assert!(rationale.contains("station 7"));
        }
        ResistanceValueOrigin::CardClaim => {
            panic!("a measured value must not report itself as the card's claim")
        }
    }

    // The supplied uncertainty describes the SUPPLIED value, and an unstated
    // one stays an honest unknown rather than inheriting the card's band.
    assert!(matches!(
        measured.uncertainty(),
        UncertaintyModel::HalfWidth { .. }
    ));
    assert!(
        base.with_measured_value(0.25, UncertaintyModel::Unstated, "   ")
            .is_err(),
        "a blank rationale is refused"
    );
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            base.with_measured_value(bad, UncertaintyModel::Unstated, "probe")
                .is_err(),
            "{bad} must be refused; perfect contact is never an implicit default"
        );
    }
}
