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
    ResistanceUncertainty, SeriesThermalResistance, ThermalInterfaces, ThermalResistanceTerm,
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

fn two_slab_mesh(n: usize) -> (ConductionMesh, usize) {
    let (left, mut positions) = box_grid([n, n, n], [1.0, 1.0, 1.0]);
    let left_vertex_count = positions.len();
    let (right, right_positions) = box_grid([n, n, n], [1.0, 1.0, 1.0]);
    let offset = u32::try_from(left_vertex_count).expect("fixture vertex count fits u32");
    let mut tets = left.tets;
    tets.extend(
        right
            .tets
            .into_iter()
            .map(|tet| tet.map(|vertex| vertex + offset)),
    );
    positions.extend(right_positions.into_iter().map(|[x, y, z]| [x + 1.0, y, z]));
    let complex = TetComplex::from_tets(positions.len(), tets);
    (
        ConductionMesh::new(complex, positions).expect("two-slab mesh"),
        left_vertex_count,
    )
}

fn boundary(mesh: &ConductionMesh) -> ThermalBoundary {
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
        .adiabatic_remainder()
        .finish()
        .expect("complete boundary partition")
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

fn level_a_contact_reference() -> (&'static fs_vvreg::thermal_level_a::ThermalLevelACase, f64) {
    let case = thermal_level_a_cases()
        .iter()
        .find(|case| case.id == "thermal-a-contact-series")
        .expect("Level-A contact row");
    assert_eq!(case.kind, ThermalLevelAKind::AnalyticReference);
    assert_eq!(case.metric, "thermal-resistance");
    (case, case.reference_value_si)
}

#[test]
fn two_slab_contact_matches_level_a_series_and_retains_receipt() {
    let (mesh, left_vertex_count) = two_slab_mesh(4);
    let boundary = boundary(&mesh);
    let material = ConductivityModel::isotropic_declared(K).expect("material");
    let source = ScalarField::Uniform(0.0);
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

    let solution = with_cx(|cx| {
        solve_with_interfaces(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
            },
            &interfaces,
            config(),
        )
        .expect("contact solve")
    });

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
    let (_, reference) = level_a_contact_reference();
    assert!(
        (series.budget().value_k_per_w - reference).abs() <= f64::EPSILON,
        "three-term floating-point sum should match the Level-A decimal reference within one ULP"
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
         \"detail\":\"R_total={:.16e} K/W; max_nodal_error={max_error:.16e} K; \
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
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
            },
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
fn interface_binding_cannot_be_reused_on_a_different_mesh() {
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
    let material = ConductivityModel::isotropic_declared(K).expect("material");
    let source = ScalarField::Uniform(0.0);
    let error = with_cx(|cx| {
        solve_with_interfaces(
            cx,
            ConductionProblem {
                mesh: &fine_mesh,
                boundary: &fine_boundary,
                material: &material,
                source: &source,
            },
            &interfaces,
            config(),
        )
        .expect_err("a mesh-bound interface object must not be reusable on another mesh")
    });
    assert!(matches!(
        error,
        ConductionError::Interface { ref interface, .. } if interface == "<binding>"
    ));
}

#[test]
fn contact_level_a_row_exposes_a_tolerance_gate() {
    let (case, reference) = level_a_contact_reference();
    assert_eq!(reference.to_bits(), 0.3f64.to_bits());
    assert!(matches!(
        case.acceptance,
        ThermalLevelAAcceptance::Tolerance { atol, rtol }
            if atol.is_finite() && atol >= 0.0 && rtol.is_finite() && rtol >= 0.0
    ));
}
