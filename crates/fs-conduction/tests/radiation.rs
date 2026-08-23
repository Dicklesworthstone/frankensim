//! G0/G1 evidence for card-backed surface radiation.

mod support;

use fs_conduction::bc::{ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
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
use support::{with_cancelled_cx, with_cx};

fn emissivity_card(value: f64, uncertainty: UncertaintyModel, finish: &str) -> MaterialCard {
    emissivity_card_with_temperature_domain(value, uncertainty, finish, 250.0, 500.0)
}

fn emissivity_card_with_temperature_domain(
    value: f64,
    uncertainty: UncertaintyModel,
    finish: &str,
    lower_temperature_k: f64,
    upper_temperature_k: f64,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with(
                "T",
                lower_temperature_k,
                upper_temperature_k,
            ),
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
    emissivity_with_temperature_domain(
        name,
        value,
        UncertaintyModel::HalfWidth {
            half_width: 0.01,
            confidence: 0.95,
        },
        finish,
        350.0,
        250.0,
        500.0,
    )
}

fn emissivity_with_temperature_domain(
    name: &str,
    value: f64,
    uncertainty: UncertaintyModel,
    finish: &str,
    query_temperature_k: f64,
    lower_temperature_k: f64,
    upper_temperature_k: f64,
) -> SurfaceEmissivity {
    SurfaceEmissivity::from_card(
        name,
        &emissivity_card_with_temperature_domain(
            value,
            uncertainty,
            finish,
            lower_temperature_k,
            upper_temperature_k,
        ),
        query_temperature_k,
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

fn cube(value: f64) -> f64 {
    value * value * value
}

fn fourth_power(value: f64) -> f64 {
    let square = value * value;
    square * square
}

/// Independent scalar reference for the two-slab fixture.
///
/// `overlay_gain = 1` is the declared algebraic coupling. `overlay_gain = 2` is the
/// previously surviving production mutant that doubled the radiosity flux
/// before applying it to conduction.
fn independent_two_slab_equilibrium(overlay_gain: f64) -> (f64, f64, f64) {
    const HOT_RESERVOIR_K: f64 = 400.0;
    const COLD_RESERVOIR_K: f64 = 300.0;
    const SLAB_CONDUCTANCE_W_K: f64 = 20.0;
    const RADIATIVE_RESISTANCE_FACTOR: f64 = 1.5;

    let mut conductive_heat_w = 400.0;
    for _ in 0..32 {
        let hot_surface_k = HOT_RESERVOIR_K - conductive_heat_w / SLAB_CONDUCTANCE_W_K;
        let cold_surface_k = COLD_RESERVOIR_K + conductive_heat_w / SLAB_CONDUCTANCE_W_K;
        let radiative_heat_w = STEFAN_BOLTZMANN_W_M2_K4
            * (fourth_power(hot_surface_k) - fourth_power(cold_surface_k))
            / RADIATIVE_RESISTANCE_FACTOR;
        let residual_w = conductive_heat_w - overlay_gain * radiative_heat_w;
        let derivative = 1.0
            + overlay_gain * 4.0 * STEFAN_BOLTZMANN_W_M2_K4
                / (RADIATIVE_RESISTANCE_FACTOR * SLAB_CONDUCTANCE_W_K)
                * (cube(hot_surface_k) + cube(cold_surface_k));
        conductive_heat_w -= residual_w / derivative;
        if residual_w.abs() <= 1.0e-13 * conductive_heat_w.abs().max(1.0) {
            break;
        }
    }
    let hot_surface_k = HOT_RESERVOIR_K - conductive_heat_w / SLAB_CONDUCTANCE_W_K;
    let cold_surface_k = COLD_RESERVOIR_K + conductive_heat_w / SLAB_CONDUCTANCE_W_K;
    let radiative_heat_w = STEFAN_BOLTZMANN_W_M2_K4
        * (fourth_power(hot_surface_k) - fourth_power(cold_surface_k))
        / RADIATIVE_RESISTANCE_FACTOR;
    let residual_w = conductive_heat_w - overlay_gain * radiative_heat_w;
    assert!(
        residual_w.abs() <= 2.0e-12 * conductive_heat_w.abs().max(1.0),
        "independent scalar Newton reference did not converge: residual={residual_w:e} W"
    );
    (hot_surface_k, cold_surface_k, conductive_heat_w)
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
fn linearization_departure_is_a_caller_budget_not_an_accuracy_bound() {
    let card = emissivity_card_with_temperature_domain(
        0.8,
        UncertaintyModel::Unstated,
        "deliberately-wide-linearization",
        1.0,
        10_000.0,
    );
    let emissivity = SurfaceEmissivity::from_card(
        "wide-budget",
        &card,
        320.0,
        SelectionPolicy::SingleClaimOnly,
    )
    .expect("broad card query");
    assert_eq!(emissivity.temperature_validity_k(), Some((1.0, 10_000.0)));
    let model = LinearizedSurfaceRadiation::new("wide-budget", emissivity, 320.0, 310.0, 5_000.0)
        .expect("the caller may declare a wide departure budget");
    let point = model
        .evaluate(5_000.0)
        .expect("evaluation lies inside both caller and material domains");

    assert!(point.discrepancy_w_m2 > 1.0e6);
    assert!(
        point.nonlinear_outward_flux_w_m2.abs() > 100.0 * point.linearized_outward_flux_w_m2.abs(),
        "admission by the caller's departure budget must not imply accuracy"
    );
}

#[test]
fn parallel_plate_literal_constructor_is_admitted_but_not_an_execution_binding() {
    let matrix = ViewFactorMatrix::infinite_parallel_plates(2.5).expect("analytic matrix");
    assert_eq!(matrix.factors(), &[vec![0.0, 1.0], vec![1.0, 0.0]]);
    assert_eq!(matrix.row_sums(), &[1.0, 1.0]);
    assert_eq!(matrix.max_reciprocity_residual(), 0.0);
    assert!(matches!(
        matrix.evidence(),
        ViewFactorEvidence::Analytic { geometry } if geometry == "infinite-parallel-plates"
    ));
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

#[test]
fn view_factor_admission_refuses_a_non_closing_row() {
    let refusal = ViewFactorMatrix::admit(
        vec![1.0, 1.0],
        vec![vec![0.0, 0.9], vec![1.0, 0.0]],
        ViewFactorEvidence::Analytic {
            geometry: "deliberately-non-closing-algebra-fixture".to_string(),
        },
        ViewFactorTolerance::default(),
    );
    assert!(
        matches!(&refusal, Err(ConductionError::Radiation { .. })),
        "a non-closing row must refuse through the radiation boundary"
    );
    if let Err(ConductionError::Radiation { surface, what, .. }) = refusal {
        assert_eq!(surface, "view-factor-row-0");
        assert!(what.contains("row sum"));
    }
}

#[test]
fn three_surface_unequal_area_radiosity_has_nontrivial_energy_closure() {
    let (complex, positions) = box_grid([2, 2, 2], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let first = RadiationSurface::new(
        &mesh,
        "first-full-face",
        |face| on_box_face(face.centroid[0], 0.0),
        emissivity("first-full-face", 0.01, "low-e-first"),
    )
    .expect("first surface");
    let second = RadiationSurface::new(
        &mesh,
        "second-full-face",
        |face| on_box_face(face.centroid[0], 1.0),
        emissivity("second-full-face", 0.6, "medium-e-second"),
    )
    .expect("second surface");
    let third = RadiationSurface::new(
        &mesh,
        "third-half-face",
        |face| on_box_face(face.centroid[1], 0.0) && face.centroid[0] < 0.5,
        emissivity("third-half-face", 0.1, "low-e-third"),
    )
    .expect("third surface");
    assert_eq!(first.area_m2().to_bits(), 1.0f64.to_bits());
    assert_eq!(second.area_m2().to_bits(), 1.0f64.to_bits());
    assert_eq!(third.area_m2().to_bits(), 0.5f64.to_bits());

    // Reciprocal exchange weights are A_i F_ij = A_j F_ji. The deliberately
    // large F_20 combined with low emissivities makes row 2 the first pivot,
    // so this fixture reaches both row swapping and three-term back-substitution.
    let factors = vec![
        vec![0.9, 0.0, 0.1],
        vec![0.0, 0.7, 0.3],
        vec![0.2, 0.6, 0.2],
    ];
    let first_diagonal = 1.0 - (1.0 - 0.01) * factors[0][0];
    let third_first_column = (1.0_f64 - 0.1) * factors[2][0];
    assert!(
        third_first_column.abs() > first_diagonal.abs(),
        "fixture must force a first-column pivot swap"
    );
    let matrix = ViewFactorMatrix::admit(
        vec![first.area_m2(), second.area_m2(), third.area_m2()],
        factors,
        ViewFactorEvidence::Analytic {
            geometry: "synthetic-three-surface-reciprocal-operator-fixture".to_string(),
        },
        ViewFactorTolerance::default(),
    )
    .expect("three-surface matrix");
    let enclosure =
        GrayDiffuseEnclosure::new(vec![first, second, third], matrix).expect("enclosure");
    let report = with_cx(|cx| {
        enclosure
            .solve(cx, &[500.0, 350.0, 250.0])
            .expect("three-surface radiosity")
    });

    assert_eq!(report.net_outward_heat_w.len(), 3);
    assert!(
        report.enclosure_energy_closure_w != 0.0,
        "unequal-area three-term accumulation must not collapse to an exact pairwise negation"
    );
    assert!(report.relative_energy_closure() < 2.0e-14);
    let radiosity_scale = report
        .radiosity_w_m2
        .iter()
        .copied()
        .map(f64::abs)
        .fold(0.0f64, f64::max);
    assert!(report.linear_residual_max_w_m2 <= 2.0e-14 * radiosity_scale);

    // Negative control: treating the half-area trace as one square metre
    // replaces A_2 q_2 with q_2. That production mutation must be conspicuous.
    let wrong_area_heat_sum = report.net_outward_heat_w[0]
        + report.net_outward_heat_w[1]
        + report.net_outward_flux_w_m2[2];
    let wrong_area_scale = report.net_outward_heat_w[0].abs()
        + report.net_outward_heat_w[1].abs()
        + report.net_outward_flux_w_m2[2].abs();
    assert!(
        wrong_area_heat_sum.abs() / wrong_area_scale > 1.0e-3,
        "negative control must reject omitting the third trace area"
    );
}

fn opposite_surfaces() -> (ConductionMesh, RadiationSurface, RadiationSurface) {
    opposite_surfaces_with_emissivities(
        emissivity("left", 0.8, "matte-left"),
        emissivity("right", 0.6, "matte-right"),
    )
}

// These cube faces are mesh/trace carriers for radiosity algebra only. Binding
// them later to literal infinite-plate factors is deliberately not a geometric
// view-factor claim for the solid cube.
fn opposite_surfaces_with_emissivities(
    left_emissivity: SurfaceEmissivity,
    right_emissivity: SurfaceEmissivity,
) -> (ConductionMesh, RadiationSurface, RadiationSurface) {
    let (complex, positions) = box_grid([2, 2, 2], [1.0, 1.0, 1.0]);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let left = RadiationSurface::new(
        &mesh,
        "left",
        |face| on_box_face(face.centroid[0], 0.0),
        left_emissivity,
    )
    .expect("left surface");
    let right = RadiationSurface::new(
        &mesh,
        "right",
        |face| on_box_face(face.centroid[0], 1.0),
        right_emissivity,
    )
    .expect("right surface");
    (mesh, left, right)
}

#[test]
fn enclosure_binding_refuses_view_factor_area_mismatch() {
    let (_mesh, left, right) = opposite_surfaces();
    let matrix = ViewFactorMatrix::infinite_parallel_plates(1.01 * left.area_m2()).expect("matrix");
    let refusal = GrayDiffuseEnclosure::new(vec![left, right], matrix);
    assert!(
        matches!(&refusal, Err(ConductionError::Radiation { .. })),
        "matrix and trace areas must bind exactly within admission tolerance"
    );
    if let Err(ConductionError::Radiation { surface, what, .. }) = refusal {
        assert_eq!(surface, "left");
        assert!(what.contains("view-factor area"));
        assert!(what.contains("mesh trace area"));
    }
}

#[test]
fn radiosity_refuses_the_near_zero_emissivity_singular_limit() {
    let left = emissivity_with_temperature_domain(
        "left",
        f64::MIN_POSITIVE,
        UncertaintyModel::Unstated,
        "near-zero-left",
        350.0,
        250.0,
        500.0,
    );
    let right = emissivity_with_temperature_domain(
        "right",
        f64::MIN_POSITIVE,
        UncertaintyModel::Unstated,
        "near-zero-right",
        350.0,
        250.0,
        500.0,
    );
    let (_mesh, left, right) = opposite_surfaces_with_emissivities(left, right);
    let matrix = ViewFactorMatrix::infinite_parallel_plates(left.area_m2()).expect("matrix");
    let enclosure = GrayDiffuseEnclosure::new(vec![left, right], matrix).expect("enclosure");
    let refusal = with_cx(|cx| enclosure.solve(cx, &[400.0, 300.0]));
    assert!(
        matches!(&refusal, Err(ConductionError::Radiation { .. })),
        "the emissivity-to-zero limiting operator must refuse as singular"
    );
    if let Err(ConductionError::Radiation { surface, what, .. }) = refusal {
        assert_eq!(surface, "radiosity-system");
        assert!(what.contains("pivot 1"));
    }
}

#[test]
fn radiosity_refuses_temperatures_outside_retained_emissivity_validity() {
    let left = emissivity_with_temperature_domain(
        "left",
        0.8,
        UncertaintyModel::Unstated,
        "narrow-left",
        350.0,
        340.0,
        360.0,
    );
    assert_eq!(left.temperature_validity_k(), Some((340.0, 360.0)));
    let right = emissivity_with_temperature_domain(
        "right",
        0.6,
        UncertaintyModel::Unstated,
        "narrow-right",
        350.0,
        340.0,
        360.0,
    );
    let (_mesh, left, right) = opposite_surfaces_with_emissivities(left, right);
    let matrix = ViewFactorMatrix::infinite_parallel_plates(left.area_m2()).expect("matrix");
    let enclosure = GrayDiffuseEnclosure::new(vec![left, right], matrix).expect("enclosure");
    let refusal = with_cx(|cx| enclosure.solve(cx, &[365.0, 350.0]));
    assert!(
        matches!(&refusal, Err(ConductionError::Radiation { .. })),
        "radiosity must not extrapolate retained emissivity"
    );
    if let Err(ConductionError::Radiation { surface, what, .. }) = refusal {
        assert_eq!(surface, "left");
        assert!(what.contains("365 K"));
        assert!(what.contains("[340, 360] K"));
    }
}

#[test]
fn two_surface_gray_diffuse_radiosity_algebra_matches_closed_form_and_replays() {
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

fn two_slab_enclosure(mesh: &ConductionMesh) -> GrayDiffuseEnclosure {
    let hot_surface = RadiationSurface::new(
        mesh,
        "hot-facing",
        |face| on_box_face(face.centroid[0], 0.5),
        emissivity("hot-facing", 0.8, "oxidized-hot"),
    )
    .expect("hot surface");
    let cold_surface = RadiationSurface::new(
        mesh,
        "cold-facing",
        |face| on_box_face(face.centroid[0], 1.5),
        emissivity("cold-facing", 0.8, "oxidized-cold"),
    )
    .expect("cold surface");
    // The finite, separated slab faces are trace carriers. F12 = F21 = 1 is a
    // caller-declared algebra fixture and is geometrically false for them.
    let matrix = ViewFactorMatrix::infinite_parallel_plates(hot_surface.area_m2()).expect("matrix");
    GrayDiffuseEnclosure::new(vec![hot_surface, cold_surface], matrix).expect("enclosure")
}

fn two_slab_boundary(mesh: &ConductionMesh, occupy_hot_radiation_face: bool) -> ThermalBoundary {
    let builder = ThermalBoundaryBuilder::new(mesh)
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
        .expect("cold region");
    let builder = if occupy_hot_radiation_face {
        builder
            .region(
                "preexisting-hot-facing-robin",
                |face| on_box_face(face.centroid[0], 0.5),
                ThermalBc::robin(5.0, 300.0).expect("occupied Robin row"),
            )
            .expect("occupied radiation face")
    } else {
        builder
    };
    builder.adiabatic_remainder().finish().expect("boundary")
}

#[test]
fn gray_diffuse_outer_fixed_point_couples_two_conduction_slabs() {
    let mesh = two_disconnected_slabs();
    let enclosure = two_slab_enclosure(&mesh);
    let boundary = two_slab_boundary(&mesh, false);
    let material = ConductivityModel::isotropic_declared(10.0).expect("material");
    let source = ScalarField::Uniform(0.0);
    let solution = with_cx(|cx| {
        solve_with_gray_diffuse_enclosure(
            cx,
            ConductionProblem { element_materials: None, mesh: &mesh,
            boundary: &boundary,
            material: &material,
            source: &source, },
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
    let (expected_hot_surface_k, expected_cold_surface_k, expected_heat_w) =
        independent_two_slab_equilibrium(1.0);
    let actual_surface_temperatures = enclosure
        .surfaces()
        .iter()
        .map(|surface| {
            surface
                .mean_temperature(&mesh, &solution.conduction.temperature)
                .expect("accepted surface temperature")
        })
        .collect::<Vec<_>>();
    let temperature_tolerance_k =
        4.0 * solution.radiation.final_threshold_k + 64.0 * f64::EPSILON * expected_hot_surface_k;
    assert!(
        (actual_surface_temperatures[0] - expected_hot_surface_k).abs() <= temperature_tolerance_k
    );
    assert!(
        (actual_surface_temperatures[1] - expected_cold_surface_k).abs() <= temperature_tolerance_k
    );
    assert!(
        (solution.radiation.radiosity.surface_temperatures_k[0] - expected_hot_surface_k).abs()
            <= temperature_tolerance_k
    );
    assert!(
        (solution.radiation.radiosity.surface_temperatures_k[1] - expected_cold_surface_k).abs()
            <= temperature_tolerance_k
    );
    let heat_tolerance_w = 50.0 * temperature_tolerance_k;
    assert!(
        (solution.radiation.radiosity.net_outward_heat_w[0] - expected_heat_w).abs()
            <= heat_tolerance_w
    );
    assert!(
        (solution.radiation.radiosity.net_outward_heat_w[1] + expected_heat_w).abs()
            <= heat_tolerance_w
    );
    let hot_conduction_heat_w = 20.0 * (400.0 - actual_surface_temperatures[0]);
    let cold_conduction_heat_w = 20.0 * (actual_surface_temperatures[1] - 300.0);
    assert!((hot_conduction_heat_w - expected_heat_w).abs() <= heat_tolerance_w);
    assert!((cold_conduction_heat_w - expected_heat_w).abs() <= heat_tolerance_w);

    // Executable negative control for the exact mutant that previously
    // survived: applying twice the radiosity flux moves both equilibrium
    // temperatures far beyond the accepted declared scalar-model envelope.
    let (doubled_hot_surface_k, doubled_cold_surface_k, _) = independent_two_slab_equilibrium(2.0);
    assert!(
        (doubled_hot_surface_k - expected_hot_surface_k).abs() > 1_000.0 * temperature_tolerance_k
    );
    assert!(
        (doubled_cold_surface_k - expected_cold_surface_k).abs()
            > 1_000.0 * temperature_tolerance_k
    );

    assert!(solution.radiation.radiosity.relative_energy_closure() < 1.0e-12);
    let coupled_energy_closure = solution.conduction.report.energy.closure_w.abs()
        / solution.radiation.radiosity.enclosure_energy_scale_w;
    assert!(coupled_energy_closure < 1.0e-8);
    println!(
        "{{\"suite\":\"fs-conduction-radiation\",\"case\":\"coupled-two-slab\",\"verdict\":\"pass\",\"outer_iterations\":{},\"surface_update_k\":{},\"hot_surface_k\":{},\"cold_surface_k\":{},\"analytic_heat_w\":{},\"radiative_heat_w\":{},\"coupled_energy_closure\":{}}}",
        solution.radiation.iterations,
        update,
        actual_surface_temperatures[0],
        actual_surface_temperatures[1],
        expected_heat_w,
        solution.radiation.radiosity.net_outward_heat_w[0],
        coupled_energy_closure
    );
}

#[test]
fn gray_diffuse_outer_fixed_point_refuses_an_exhausted_budget() {
    let mesh = two_disconnected_slabs();
    let enclosure = two_slab_enclosure(&mesh);
    let boundary = two_slab_boundary(&mesh, false);
    let material = ConductivityModel::isotropic_declared(10.0).expect("material");
    let source = ScalarField::Uniform(0.0);
    let refusal = with_cx(|cx| {
        solve_with_gray_diffuse_enclosure(
            cx,
            ConductionProblem { element_materials: None, mesh: &mesh,
            boundary: &boundary,
            material: &material,
            source: &source, },
            None,
            &enclosure,
            config(),
            CoupledRadiationConfig {
                surface_temperature_rtol: 0.0,
                surface_temperature_atol_k: 0.0,
                relaxation: 0.1,
                max_iterations: 1,
            },
        )
    });
    assert!(
        matches!(&refusal, Err(ConductionError::Radiation { .. })),
        "one deliberately under-relaxed iteration must exhaust its outer budget"
    );
    if let Err(ConductionError::Radiation { surface, what, .. }) = refusal {
        assert_eq!(surface, "coupling");
        assert!(what.contains("did not converge in 1 iterations"));
        assert!(what.contains("final update"));
    }
}

#[test]
fn radiation_overlay_refuses_a_preexisting_robin_row() {
    let mesh = two_disconnected_slabs();
    let enclosure = two_slab_enclosure(&mesh);
    let boundary = two_slab_boundary(&mesh, true);
    let material = ConductivityModel::isotropic_declared(10.0).expect("material");
    let source = ScalarField::Uniform(0.0);
    let refusal = with_cx(|cx| {
        solve_with_gray_diffuse_enclosure(
            cx,
            ConductionProblem { element_materials: None, mesh: &mesh,
            boundary: &boundary,
            material: &material,
            source: &source, },
            None,
            &enclosure,
            config(),
            CoupledRadiationConfig::default(),
        )
    });
    assert!(
        matches!(&refusal, Err(ConductionError::Radiation { .. })),
        "radiation must never overwrite a preexisting physical boundary row"
    );
    if let Err(ConductionError::Radiation { surface, what, .. }) = refusal {
        assert_eq!(surface, "hot-facing");
        assert!(what.contains("already carries region"));
        assert!(what.contains("preexisting-hot-facing-robin"));
    }
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
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, TEMPERATURE_DIMS),
            value: PropertyValue::Scalar {
                value: 0.8,
                dims: TEMPERATURE_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 250.0, 500.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: "deliberately dimensionally wrong emissivity fixture".to_string(),
                license: "internal-test-use".to_string(),
                artifact: None,
            },
        })
        .expect("wrong-dimensional claim is structurally valid");
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "fixture-alloy".to_string(),
            phase: "solid".to_string(),
            process: "wrong-emissivity-dimensions".to_string(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .expect("wrong-dimensional card");
    let refusal = SurfaceEmissivity::from_card(
        "wrong-dimensional-surface",
        &card,
        350.0,
        SelectionPolicy::SingleClaimOnly,
    );
    assert!(
        matches!(
            refusal,
            Err(ConductionError::Dimensions {
                expected,
                found,
                ..
            }) if expected == EMISSIVITY_DIMS.0 && found == TEMPERATURE_DIMS.0
        ),
        "wrong-dimensional emissivity must refuse with the declared temperature dimensions"
    );
}
