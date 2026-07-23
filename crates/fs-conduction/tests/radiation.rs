//! G0/G1 evidence for card-backed surface radiation.

mod support;

use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::material::ConductivityModel;
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{
    ConductionProblem, InitialGuess, LinearConfig, Nonlinearity, SolveConfig, StopRule,
};
use fs_conduction::{
    ConductionError, CoupledRadiationConfig, EMISSIVITY_DIMS, GrayDiffuseEnclosure,
    LinearizedSurfaceRadiation, RadiationSurface, STEFAN_BOLTZMANN_W_M2_K4,
    SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity, TEMPERATURE_DIMS, ViewFactorEvidence,
    ViewFactorMatrix, ViewFactorTolerance, solve_with_gray_diffuse_enclosure,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_rep_mesh::TetComplex;
use fs_vvreg::thermal_level_a::{ThermalLevelAKind, thermal_level_a_cases};
use support::{with_cancelled_cx, with_cx};

fn emissivity_card(value: f64, uncertainty: UncertaintyModel, finish: &str) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 250.0, 500.0),
            uncertainty,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: format!("G1 radiation fixture surface finish {finish}"),
                license: "internal-test-use".to_string(),
                artifact: None,
            },
        })
        .expect("emissivity claim inserts");
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture-alloy".to_string(),
            phase: "solid".to_string(),
            process: finish.to_string(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .expect("material card")
}

fn emissivity(name: &str, value: f64, finish: &str) -> SurfaceEmissivity {
    SurfaceEmissivity::from_card(
        name,
        &emissivity_card(
            value,
            UncertaintyModel::HalfWidth {
                half_width: 0.01,
                confidence: 0.95,
            },
            finish,
        ),
        350.0,
        SelectionPolicy::SingleClaimOnly,
    )
    .expect("emissivity resolves")
}

fn config() -> SolveConfig {
    SolveConfig {
        nonlinearity: Nonlinearity::FixedPoint {
            relaxation: 1.0,
            max_backtracks: 8,
        },
        stop: StopRule {
            residual_rtol: 1.0e-10,
            residual_atol: 1.0e-24,
            step_atol: 0.0,
            max_iterations: 20,
        },
        linear: LinearConfig {
            tolerance: 1.0e-12,
            max_iterations: 50_000,
            restart: 60,
        },
        initial: InitialGuess::DirichletMean,
    }
}

#[test]
fn linearized_robin_retains_card_and_measures_t4_discrepancy() {
    let card = emissivity_card(
        0.8,
        UncertaintyModel::HalfWidth {
            half_width: 0.02,
            confidence: 0.95,
        },
        "black-anodized-v3",
    );
    let emissivity =
        SurfaceEmissivity::from_card("radiator", &card, 320.0, SelectionPolicy::SingleClaimOnly)
            .expect("card query");
    let model = LinearizedSurfaceRadiation::new("radiator", emissivity, 320.0, 310.0, 20.0)
        .expect("linearization");
    let point = model.evaluate(330.0).expect("in-domain evaluation");
    let expected_h = 4.0 * 0.8 * STEFAN_BOLTZMANN_W_M2_K4 * 320.0 * 320.0 * 320.0;
    let expected_full = 0.8
        * STEFAN_BOLTZMANN_W_M2_K4
        * (fs_math::det::powi(330.0, 4) - fs_math::det::powi(310.0, 4));
    assert!((point.h_rad_w_m2k - expected_h).abs() <= f64::EPSILON * expected_h);
    assert!(
        (point.linearized_outward_flux_w_m2 - expected_h * 20.0).abs()
            <= 2.0 * f64::EPSILON * expected_h * 20.0
    );
    assert!(
        (point.nonlinear_outward_flux_w_m2 - expected_full).abs()
            <= 2.0 * f64::EPSILON * expected_full
    );
    assert!(point.discrepancy_w_m2 > 0.0);
    assert_eq!(point.uncertainty_confidence, Some(0.95));
    assert_eq!(model.emissivity().card_identity(), card.content_hash());
    assert_eq!(
        model.emissivity().receipt().property,
        SURFACE_EMISSIVITY_PROPERTY
    );
    assert!(
        model
            .emissivity()
            .material_state()
            .contains("black-anodized-v3")
    );
    assert!(matches!(point.boundary, ThermalBc::Robin { .. }));

    assert!(matches!(
        model.evaluate(341.0),
        Err(ConductionError::Radiation { .. })
    ));
    println!(
        "{{\"suite\":\"fs-conduction-radiation\",\"case\":\"linearized-card\",\"verdict\":\"pass\",\"h_rad_w_m2k\":{expected_h},\"full_flux_w_m2\":{expected_full},\"discrepancy_w_m2\":{}}}",
        point.discrepancy_w_m2
    );
}

#[test]
fn parallel_plate_view_factor_binds_level_a_and_reciprocity_laws() {
    let reference = thermal_level_a_cases()
        .iter()
        .find(|case| case.id == "thermal-a-parallel-plate-view-factor")
        .expect("catalog row");
    assert_eq!(reference.kind, ThermalLevelAKind::AnalyticReference);
    let matrix = ViewFactorMatrix::infinite_parallel_plates(2.5).expect("analytic matrix");
    assert_eq!(matrix.factors()[0][1], reference.reference_value_si);
    assert_eq!(matrix.factors()[1][0], reference.reference_value_si);
    assert_eq!(matrix.row_sums(), &[1.0, 1.0]);
    assert_eq!(matrix.max_reciprocity_residual(), 0.0);
    assert_eq!(
        matrix.identity(),
        ViewFactorMatrix::infinite_parallel_plates(2.5)
            .expect("replay")
            .identity()
    );

    let nonreciprocal = ViewFactorMatrix::admit(
        vec![1.0, 1.0],
        vec![vec![0.0, 1.0], vec![0.9, 0.1]],
        ViewFactorEvidence::Analytic {
            geometry: "bad-fixture".to_string(),
        },
        ViewFactorTolerance::default(),
    );
    assert!(matches!(
        nonreciprocal,
        Err(ConductionError::Radiation { .. })
    ));
    let incomplete_qmc = ViewFactorMatrix::admit(
        vec![1.0, 1.0],
        vec![vec![0.0, 1.0], vec![1.0, 0.0]],
        ViewFactorEvidence::ExternalQmc {
            seed: 7,
            samples: 0,
            generator: String::new(),
        },
        ViewFactorTolerance::default(),
    );
    assert!(matches!(
        incomplete_qmc,
        Err(ConductionError::Radiation { .. })
    ));
}

fn opposite_surfaces() -> (ConductionMesh, RadiationSurface, RadiationSurface) {
    let (complex, positions) = box_grid([2, 2, 2], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let left = RadiationSurface::new(
        &mesh,
        "left",
        |face| on_box_face(face.centroid[0], 0.0),
        emissivity("left", 0.8, "matte-left"),
    )
    .expect("left surface");
    let right = RadiationSurface::new(
        &mesh,
        "right",
        |face| on_box_face(face.centroid[0], 1.0),
        emissivity("right", 0.6, "matte-right"),
    )
    .expect("right surface");
    (mesh, left, right)
}

#[test]
fn two_surface_gray_diffuse_radiosity_matches_closed_form_and_replays() {
    let (_mesh, left, right) = opposite_surfaces();
    assert_eq!(left.area_m2().to_bits(), right.area_m2().to_bits());
    let matrix = ViewFactorMatrix::infinite_parallel_plates(left.area_m2()).expect("matrix");
    let enclosure = GrayDiffuseEnclosure::new(vec![left, right], matrix).expect("enclosure");
    let first = with_cx(|cx| enclosure.solve(cx, &[400.0, 300.0]).expect("radiosity"));
    let second = with_cx(|cx| {
        enclosure
            .solve(cx, &[400.0, 300.0])
            .expect("radiosity replay")
    });
    assert_eq!(first, second);
    let expected_flux = STEFAN_BOLTZMANN_W_M2_K4
        * (fs_math::det::powi(400.0, 4) - fs_math::det::powi(300.0, 4))
        / (1.0 / 0.8 + 1.0 / 0.6 - 1.0);
    assert!((first.net_outward_flux_w_m2[0] - expected_flux).abs() <= 2.0e-13 * expected_flux);
    assert!((first.net_outward_flux_w_m2[1] + expected_flux).abs() <= 2.0e-13 * expected_flux);
    assert!(first.relative_energy_closure() <= 2.0e-15);
    assert!(first.linear_residual_max_w_m2 <= 1.0e-12);
    assert!(first.emissivity_uncertainty_complete);
    println!(
        "{{\"suite\":\"fs-conduction-radiation\",\"case\":\"two-surface-radiosity\",\"verdict\":\"pass\",\"heat_flux_w_m2\":{},\"relative_energy_closure\":{}}}",
        first.net_outward_flux_w_m2[0],
        first.relative_energy_closure()
    );
}

fn two_disconnected_slabs() -> ConductionMesh {
    let (a, mut positions_a) = box_grid([2, 2, 2], [0.5, 1.0, 1.0]);
    let (b, positions_b) = box_grid([2, 2, 2], [0.5, 1.0, 1.0]);
    let offset = positions_a.len() as u32;
    positions_a.extend(positions_b.into_iter().map(|mut point| {
        point[0] += 1.5;
        point
    }));
    let mut tets = a.tets;
    tets.extend(b.tets.into_iter().map(|tet| {
        [
            tet[0] + offset,
            tet[1] + offset,
            tet[2] + offset,
            tet[3] + offset,
        ]
    }));
    ConductionMesh::new(TetComplex::from_tets(positions_a.len(), tets), positions_a).expect("mesh")
}

#[test]
fn gray_diffuse_outer_fixed_point_couples_two_conduction_slabs() {
    let mesh = two_disconnected_slabs();
    let hot_surface = RadiationSurface::new(
        &mesh,
        "hot-facing",
        |face| on_box_face(face.centroid[0], 0.5),
        emissivity("hot-facing", 0.8, "oxidized-hot"),
    )
    .expect("hot surface");
    let cold_surface = RadiationSurface::new(
        &mesh,
        "cold-facing",
        |face| on_box_face(face.centroid[0], 1.5),
        emissivity("cold-facing", 0.8, "oxidized-cold"),
    )
    .expect("cold surface");
    let matrix = ViewFactorMatrix::infinite_parallel_plates(hot_surface.area_m2()).expect("matrix");
    let enclosure =
        GrayDiffuseEnclosure::new(vec![hot_surface, cold_surface], matrix).expect("enclosure");
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot-reservoir",
            |face| on_box_face(face.centroid[0], 0.0),
            ThermalBc::dirichlet(400.0).expect("hot bc"),
        )
        .expect("hot region")
        .region(
            "cold-reservoir",
            |face| on_box_face(face.centroid[0], 2.0),
            ThermalBc::dirichlet(300.0).expect("cold bc"),
        )
        .expect("cold region")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let material = ConductivityModel::isotropic_declared(10.0).expect("material");
    let source = ScalarField::Uniform(0.0);
    let solution = with_cx(|cx| {
        solve_with_gray_diffuse_enclosure(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
            },
            None,
            &enclosure,
            config(),
            CoupledRadiationConfig {
                surface_temperature_rtol: 1.0e-9,
                surface_temperature_atol_k: 1.0e-9,
                relaxation: 0.5,
                max_iterations: 120,
            },
        )
        .expect("coupled solve")
    });
    let update = *solution
        .radiation
        .surface_update_history_k
        .last()
        .expect("outer update");
    assert!(update <= solution.radiation.final_threshold_k);
    assert!(solution.radiation.radiosity.net_outward_heat_w[0] > 0.0);
    assert!(solution.radiation.radiosity.net_outward_heat_w[1] < 0.0);
    assert!(solution.radiation.radiosity.relative_energy_closure() < 1.0e-12);
    let coupled_energy_closure = solution.conduction.report.energy.closure_w.abs()
        / solution.radiation.radiosity.enclosure_energy_scale_w;
    assert!(coupled_energy_closure < 1.0e-8);
    println!(
        "{{\"suite\":\"fs-conduction-radiation\",\"case\":\"coupled-two-slab\",\"verdict\":\"pass\",\"outer_iterations\":{},\"surface_update_k\":{},\"radiative_heat_w\":{},\"coupled_energy_closure\":{}}}",
        solution.radiation.iterations,
        update,
        solution.radiation.radiosity.net_outward_heat_w[0],
        coupled_energy_closure
    );
}

#[test]
fn radiation_refuses_missing_cards_overlap_and_cancellation() {
    let empty = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture".to_string(),
            phase: "solid".to_string(),
            process: "unknown-finish".to_string(),
            revision: 0,
        },
        ClaimSet::new(),
        Vec::new(),
    )
    .expect("empty card");
    assert!(matches!(
        SurfaceEmissivity::from_card("missing", &empty, 350.0, SelectionPolicy::SingleClaimOnly),
        Err(ConductionError::MaterialQuery { .. })
    ));

    let (mesh, left, right) = opposite_surfaces();
    let overlapping = RadiationSurface::new(
        &mesh,
        "overlapping-left",
        |face| on_box_face(face.centroid[0], 0.0),
        emissivity("overlapping-left", 0.7, "overlapping-finish"),
    )
    .expect("overlapping surface");
    let overlap_matrix =
        ViewFactorMatrix::infinite_parallel_plates(left.area_m2()).expect("overlap matrix");
    assert!(matches!(
        GrayDiffuseEnclosure::new(vec![left.clone(), overlapping], overlap_matrix),
        Err(ConductionError::Radiation { .. })
    ));

    let (wrong_complex, wrong_positions) = box_grid([1, 1, 1], [1.0, 1.0, 1.0]);
    let wrong_mesh = ConductionMesh::new(wrong_complex, wrong_positions).expect("wrong mesh");
    assert!(matches!(
        left.mean_temperature(&wrong_mesh, &vec![350.0; wrong_mesh.vertex_count()]),
        Err(ConductionError::Radiation { .. })
    ));

    let matrix = ViewFactorMatrix::infinite_parallel_plates(left.area_m2()).expect("matrix");
    let enclosure = GrayDiffuseEnclosure::new(vec![left, right], matrix).expect("enclosure");
    let cancelled = with_cancelled_cx(|cx| enclosure.solve(cx, &[400.0, 300.0]));
    assert!(matches!(
        cancelled,
        Err(ConductionError::Cancelled {
            stage: "radiation-radiosity-assemble",
            at: 0
        })
    ));
}

#[test]
fn emissivity_query_uses_absolute_temperature_dimensions() {
    assert_eq!(TEMPERATURE_DIMS, fs_qty::Dims([0, 0, 0, 1, 0, 0]));
    assert_eq!(EMISSIVITY_DIMS, fs_qty::Dims::NONE);
}
