//! G0/G1/G3/G4: exact source segments, common nonlinear flux and bounded refusal.
mod support;

use fs_conduction::layered::{LayeredWall, NonlinearWallOptions, PlanarThermalLayer};
use fs_conduction::{CONDUCTIVITY_DIMS, ConductionError, ConductivityTable, TEMPERATURE_DIMS};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};
use support::{with_cancelled_cx, with_cx};

fn options() -> NonlinearWallOptions {
    NonlinearWallOptions {
        max_iterations: 128,
        temperature_tolerance_k: 1e-10,
        flux_absolute_tolerance_w_per_m2: 1e-8,
        flux_relative_tolerance: 1e-11,
    }
}

fn curve(name: &str, thickness: f64, knots: Vec<(f64, f64)>) -> PlanarThermalLayer {
    PlanarThermalLayer::new(
        name,
        thickness,
        ConductivityTable::declared_curve(knots).unwrap(),
    )
    .unwrap()
}

fn affine() -> LayeredWall {
    LayeredWall::new(
        vec![curve("affine", 0.1, vec![(300.0, 1.0), (600.0, 2.0)])],
        vec![],
        [0.0; 2],
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-10 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_affine_conductivity_uses_integrated_law_and_signed_flow() {
    // k=T/300 and L=0.1: integral_300^600 k dT / L = 4500 W/m2.
    let wall = affine();
    let answer =
        with_cx(|cx| wall.solve_temperature_dependent(cx, [600.0, 300.0], 2.0, options())).unwrap();
    close(answer.heat_flux_w_per_m2, 4500.0);
    close(answer.heat_rate_w, 9000.0);
    close(answer.effective_resistance_m2_k_per_w, 1.0 / 15.0);
    close(answer.effective_transmittance_w_per_m2_k, 15.0);
    close(answer.layer_temperatures[0].left_k, 600.0);
    close(answer.layer_temperatures[0].right_k, 300.0);
    assert!(answer.maximum_flux_residual_w_per_m2 < 1e-8);
    assert!(answer.cold_boundary_residual_k <= options().temperature_tolerance_k);
    assert_eq!(answer.wall, wall);
    let reverse =
        with_cx(|cx| wall.solve_temperature_dependent(cx, [300.0, 600.0], 2.0, options())).unwrap();
    close(reverse.heat_flux_w_per_m2, -4500.0);
    close(reverse.layer_temperatures[0].left_k, 300.0);
    close(reverse.layer_temperatures[0].right_k, 600.0);
    close(
        reverse.effective_resistance_m2_k_per_w,
        answer.effective_resistance_m2_k_per_w,
    );
}

#[test]
fn g1_interior_knot_cannot_be_replaced_by_an_endpoint_chord() {
    let wall = LayeredWall::new(
        vec![curve(
            "peak",
            0.2,
            vec![(300.0, 1.0), (400.0, 3.0), (500.0, 2.0)],
        )],
        vec![],
        [0.0; 2],
    )
    .unwrap();
    let answer =
        with_cx(|cx| wall.solve_temperature_dependent(cx, [500.0, 300.0], 1.0, options())).unwrap();
    // Independent trapezoids: 100*(1+3)/2 + 100*(3+2)/2 = 450.
    close(answer.heat_flux_w_per_m2, 2250.0);
    assert!(
        (answer.heat_flux_w_per_m2 - 1500.0).abs() > 700.0,
        "an endpoint-only chord would erase the interior maximum"
    );
}

fn layered_problem() -> (LayeredWall, f64) {
    let hot = 20.0 + (327.5_f64.powi(2) + 30_000.0).sqrt();
    let wall = LayeredWall::new(
        vec![
            curve("nonlinear", 0.5, vec![(300.0, 1.0), (450.0, 1.5)]),
            PlanarThermalLayer::new("constant", 0.25, ConductivityTable::declared(2.0).unwrap())
                .unwrap(),
        ],
        vec![0.1],
        [0.2, 0.05],
    )
    .unwrap();
    (wall, hot)
}

#[test]
fn g1_g3_multilayer_contacts_films_and_reflection_preserve_one_flux() {
    let (wall, hot) = layered_problem();
    // Construct the reservoirs independently from q=100:
    // right film ->305 K; constant layer ->317.5 K; contact ->327.5 K;
    // nonlinear layer: (T_left^2-327.5^2)/600 = q*0.5.
    let answer =
        with_cx(|cx| wall.solve_temperature_dependent(cx, [hot, 300.0], 3.0, options())).unwrap();
    close(answer.heat_flux_w_per_m2, 100.0);
    close(answer.heat_rate_w, 300.0);
    close(answer.layer_temperatures[0].left_k, hot - 20.0);
    close(answer.layer_temperatures[0].right_k, 327.5);
    close(answer.layer_temperatures[1].left_k, 317.5);
    close(answer.layer_temperatures[1].right_k, 305.0);
    for layer in &answer.layer_temperatures {
        close(layer.recomputed_flux_w_per_m2, 100.0);
    }
    assert!(answer.maximum_flux_residual_w_per_m2 <= 1.1e-8);
    let reflected = LayeredWall::new(
        wall.layers().iter().cloned().rev().collect(),
        vec![0.1],
        [0.05, 0.2],
    )
    .unwrap();
    let reversed =
        with_cx(|cx| reflected.solve_temperature_dependent(cx, [300.0, hot], 3.0, options()))
            .unwrap();
    close(reversed.heat_flux_w_per_m2, -100.0);
    for (a, b) in answer
        .layer_temperatures
        .iter()
        .zip(reversed.layer_temperatures.iter().rev())
    {
        close(a.left_k, b.right_k);
        close(a.right_k, b.left_k);
    }
    assert_eq!(
        answer,
        with_cx(|cx| wall.solve_temperature_dependent(cx, [hot, 300.0], 3.0, options(),)).unwrap()
    );
}

#[test]
fn g0_equal_temperature_has_zero_heat_and_local_differential_resistance() {
    let wall = LayeredWall::new(affine().layers().to_vec(), vec![], [0.05, 0.05]).unwrap();
    let result =
        with_cx(|cx| wall.solve_temperature_dependent(cx, [450.0; 2], 2.0, options())).unwrap();
    assert_eq!(result.heat_flux_w_per_m2, 0.0);
    assert_eq!(result.heat_rate_w, 0.0);
    assert_eq!(result.iterations, 0);
    close(result.effective_resistance_m2_k_per_w, 0.1 + 0.1 / 1.5);
    close(result.effective_transmittance_w_per_m2_k, 6.0);
    assert_eq!(result.layer_temperatures[0].left_k, 450.0);
    assert_eq!(result.layer_temperatures[0].right_k, 450.0);
}

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn key() -> PropertyKey {
    PropertyKey::with_quantity(
        "thermal-k",
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        )),
    )
}

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", temperature(), t)
        .unwrap()
        .with_quantity("moisture", QuantitySpec::dimensional(Dims::NONE), 0.12)
        .unwrap()
}

#[test]
fn g0_g3_real_query_owner_preserves_source_knots_context_and_receipts() {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key(),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: TEMPERATURE_DIMS,
                knots: vec![(300.0, 1.0), (400.0, 3.0), (500.0, 2.0)],
                dims: CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("temperature", temperature(), 300.0, 500.0)
                .with_quantity(
                    "moisture",
                    QuantitySpec::dimensional(Dims::NONE),
                    0.12,
                    0.12,
                ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: Vec::new(),
            provenance: Provenance {
                source: "manufactured piecewise-linear wall source; not empirical material data"
                    .into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "numerical fixture".into(),
            phase: "solid".into(),
            process: "fixed-moisture test".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    let layer = PlanarThermalLayer::from_card(
        "sourced",
        0.2,
        &card,
        &key(),
        "temperature",
        &[point(300.0), point(500.0)],
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    )
    .unwrap();
    assert_eq!(
        layer.conductivity().knots(),
        &[(300.0, 1.0), (400.0, 3.0), (500.0, 2.0)]
    );
    let receipts = layer.conductivity().receipts().to_vec();
    let wall = LayeredWall::new(vec![layer], vec![], [0.0; 2]).unwrap();
    let result =
        with_cx(|cx| wall.solve_temperature_dependent(cx, [500.0, 300.0], 1.0, options())).unwrap();
    close(result.heat_flux_w_per_m2, 2250.0);
    let retained = &result.wall.layers()[0];
    assert_eq!(retained.material_card(), Some(card.content_hash()));
    assert_eq!(retained.conductivity().receipts(), receipts.as_slice());
    for receipt in retained.conductivity().receipts() {
        card.claims().verify_receipt(receipt).unwrap();
        assert!(receipt.query_point.contains(&("moisture".into(), 0.12)));
    }
}

#[test]
fn g0_full_reservoir_support_is_explicit_and_constant_path_is_unchanged() {
    let layer = curve("bounded constant", 1.0, vec![(300.0, 1.0), (500.0, 1.0)]);
    let wall = LayeredWall::new(vec![layer], vec![], [1.0; 2]).unwrap();
    assert!(with_cx(|cx| wall.solve(cx, [290.0, 510.0], 1.0)).is_ok());
    assert!(matches!(
        with_cx(|cx| wall.solve_temperature_dependent(cx, [290.0, 510.0], 1.0, options(),)),
        Err(ConductionError::OutsideTemperatureSpan { .. })
    ));
}

#[test]
fn g0_g4_budget_invalid_tolerances_and_cancellation_refuse_atomically() {
    let wall = affine();
    let before = wall.clone();
    let mut config = options();
    config.max_iterations = 1;
    assert!(
        with_cx(|cx| wall.solve_temperature_dependent(cx, [600.0, 300.0], 1.0, config,)).is_err()
    );
    config.max_iterations = 0;
    assert!(
        with_cx(|cx| wall.solve_temperature_dependent(cx, [600.0, 300.0], 1.0, config,)).is_err()
    );
    for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        config = options();
        config.temperature_tolerance_k = tolerance;
        assert!(
            with_cx(|cx| wall.solve_temperature_dependent(cx, [600.0, 300.0], 1.0, config,))
                .is_err()
        );
    }
    config = options();
    config.flux_absolute_tolerance_w_per_m2 = 0.0;
    config.flux_relative_tolerance = 0.0;
    assert!(
        with_cx(|cx| wall.solve_temperature_dependent(cx, [600.0, 300.0], 1.0, config,)).is_err()
    );
    assert!(matches!(
        with_cancelled_cx(|cx| wall.solve_temperature_dependent(
            cx,
            [600.0, 300.0],
            1.0,
            options(),
        )),
        Err(ConductionError::Cancelled { .. })
    ));
    assert_eq!(wall, before);
}
