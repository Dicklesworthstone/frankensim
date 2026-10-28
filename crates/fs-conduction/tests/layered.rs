//! G0/G1/G3/G4 tests for source-retaining planar steady heat flow.
mod support;

use fs_conduction::layered::{LayeredWall, MAX_LAYERED_WALL_LAYERS, PlanarThermalLayer};
use fs_conduction::{CONDUCTIVITY_DIMS, ConductionError, ConductivityTable, ProvenanceClass};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec};
use support::{with_cancelled_cx, with_cx};

fn constant(name: &str, thickness: f64, k: f64) -> PlanarThermalLayer {
    PlanarThermalLayer::new(name, thickness, ConductivityTable::declared(k).unwrap()).unwrap()
}

fn wall() -> LayeredWall {
    LayeredWall::new(
        vec![
            constant("left slab", 1.0, 1.0),
            constant("right slab", 2.0, 2.0),
        ],
        vec![0.5],
        [0.25, 0.25],
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-11 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_series_wall_retains_contact_jumps_and_closes_fourier_flux() {
    let result = with_cx(|cx| wall().solve(cx, [340.0, 280.0], 2.5)).unwrap();
    // Two 1 m2 K/W layers, a 0.5 contact, and two 0.25 films: 3 m2 K/W.
    close(result.total_resistance_m2_k_per_w, 3.0);
    close(result.thermal_transmittance_w_per_m2_k, 1.0 / 3.0);
    close(result.heat_flux_w_per_m2, 20.0);
    close(result.heat_rate_w, 50.0);
    let left = &result.layer_temperatures[0];
    let right = &result.layer_temperatures[1];
    close(left.left_k, 335.0);
    close(left.right_k, 315.0);
    close(right.left_k, 305.0);
    close(right.right_k, 285.0);
    close(left.recomputed_flux_w_per_m2, 20.0);
    close(right.recomputed_flux_w_per_m2, 20.0);
    close((left.right_k - right.left_k) / 0.5, 20.0);
    close(result.maximum_flux_residual_w_per_m2, 0.0);
    close(result.right_boundary_residual_k, 0.0);
    assert_eq!(result.wall, wall());
}

#[test]
fn g3_reverse_heat_area_scaling_and_equal_temperatures() {
    let stack = wall();
    let original = with_cx(|cx| stack.solve(cx, [340.0, 280.0], 2.5)).unwrap();
    let reverse = with_cx(|cx| stack.solve(cx, [280.0, 340.0], 2.5)).unwrap();
    close(reverse.heat_flux_w_per_m2, -original.heat_flux_w_per_m2);
    close(reverse.heat_rate_w, -original.heat_rate_w);
    for (forward, backward) in original
        .layer_temperatures
        .iter()
        .zip(&reverse.layer_temperatures)
    {
        close(forward.left_k + backward.left_k, 620.0);
        close(forward.right_k + backward.right_k, 620.0);
    }
    let larger = with_cx(|cx| stack.solve(cx, [340.0, 280.0], 10.0)).unwrap();
    assert_eq!(larger.layer_temperatures, original.layer_temperatures);
    close(larger.heat_rate_w, 4.0 * original.heat_rate_w);
    let equal = with_cx(|cx| stack.solve(cx, [310.0, 310.0], 2.5)).unwrap();
    assert_eq!(equal.heat_rate_w, 0.0);
    for layer in equal.layer_temperatures {
        assert_eq!(layer.left_k, 310.0);
        assert_eq!(layer.right_k, 310.0);
    }
    assert_eq!(
        original,
        with_cx(|cx| stack.solve(cx, [340.0, 280.0], 2.5)).unwrap()
    );
}

#[test]
fn g3_material_thickness_and_subdivision_change_the_physics() {
    let baseline = LayeredWall::new(vec![constant("a", 0.1, 0.04)], vec![], [0.0; 2]).unwrap();
    let thicker = LayeredWall::new(vec![constant("a", 0.2, 0.04)], vec![], [0.0; 2]).unwrap();
    let conductor = LayeredWall::new(vec![constant("a", 0.1, 0.08)], vec![], [0.0; 2]).unwrap();
    let split = LayeredWall::new(
        vec![constant("a", 0.04, 0.04), constant("b", 0.06, 0.04)],
        vec![0.0],
        [0.0; 2],
    )
    .unwrap();
    let q = |stack: &LayeredWall| with_cx(|cx| stack.solve(cx, [300.0, 280.0], 1.0)).unwrap();
    close(q(&baseline).heat_rate_w, 8.0);
    close(q(&thicker).heat_rate_w, 4.0);
    close(q(&conductor).heat_rate_w, 16.0);
    close(q(&split).heat_rate_w, 8.0);
}

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn conductivity_key() -> PropertyKey {
    PropertyKey::with_quantity(
        "specimen-conductivity",
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        )),
    )
}

fn point(t: f64, moisture: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("temperature", temperature(), t)
        .unwrap()
        .with_quantity("moisture", QuantitySpec::dimensional(Dims::NONE), moisture)
        .unwrap()
}

fn source_card(k: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: conductivity_key(),
            value: PropertyValue::Scalar {
                value: k,
                dims: CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("temperature", temperature(), 280.0, 360.0)
                .with_quantity(
                    "moisture",
                    QuantitySpec::dimensional(Dims::NONE),
                    0.12,
                    0.12,
                ),
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            uncertainty: UncertaintyModel::Unstated,
            observations: Vec::new(),
            provenance: Provenance {
                source: "manufactured constant conductivity fixture; not measured material data"
                    .into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured wall conductor".into(),
            phase: "solid".into(),
            process: "fixed-moisture numerical fixture".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}

fn sourced(card: &MaterialCard) -> PlanarThermalLayer {
    PlanarThermalLayer::from_card(
        "source layer",
        0.1,
        card,
        &conductivity_key(),
        "temperature",
        &[point(280.0, 0.12), point(360.0, 0.12)],
        ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    )
    .unwrap()
}

#[test]
fn g0_g3_source_card_context_and_receipts_survive_the_heat_solve() {
    let card = source_card(0.2);
    let layer = sourced(&card);
    let original_receipts = layer.conductivity().receipts().to_vec();
    let stack = LayeredWall::new(vec![layer], vec![], [0.0; 2]).unwrap();
    let result = with_cx(|cx| stack.solve(cx, [340.0, 300.0], 1.5)).unwrap();
    close(result.heat_flux_w_per_m2, 80.0);
    close(result.heat_rate_w, 120.0);
    let retained = &result.wall.layers()[0];
    assert_eq!(retained.material_card(), Some(card.content_hash()));
    assert_eq!(
        retained.conductivity().provenance(),
        ProvenanceClass::MatdbReceipts
    );
    assert_eq!(
        retained.conductivity().receipts(),
        original_receipts.as_slice()
    );
    for receipt in retained.conductivity().receipts() {
        card.claims().verify_receipt(receipt).unwrap();
        assert!(receipt.query_point.contains(&("moisture".into(), 0.12)));
    }
    let other = source_card(0.4);
    let changed = LayeredWall::new(vec![sourced(&other)], vec![], [0.0; 2]).unwrap();
    let changed = with_cx(|cx| changed.solve(cx, [340.0, 300.0], 1.5)).unwrap();
    close(changed.heat_rate_w, 240.0);
    assert_ne!(
        changed.wall.layers()[0].material_card(),
        retained.material_card()
    );
    assert_ne!(
        changed.wall.layers()[0].conductivity().receipts(),
        original_receipts.as_slice()
    );
}

#[test]
fn g0_source_domains_are_checked_on_the_actual_layer_faces() {
    let card = source_card(0.2);
    let layer = sourced(&card);
    assert!(
        PlanarThermalLayer::from_card(
            "bad moisture",
            0.1,
            &card,
            &conductivity_key(),
            "temperature",
            &[point(280.0, 0.13), point(360.0, 0.13)],
            ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        )
        .is_err()
    );
    let stack = LayeredWall::new(vec![layer.clone()], vec![], [0.0; 2]).unwrap();
    assert!(matches!(
        with_cx(|cx| stack.solve(cx, [370.0, 300.0], 1.0)),
        Err(ConductionError::OutsideTemperatureSpan { .. })
    ));
    // The reservoirs are outside the source span, but explicit films keep the
    // layer itself within it: q=100/(0.5+0.5+0.5), faces=346.667,313.333 K.
    let protected = LayeredWall::new(vec![layer], vec![], [0.5, 0.5]).unwrap();
    let result = with_cx(|cx| protected.solve(cx, [380.0, 280.0], 1.0)).unwrap();
    close(result.layer_temperatures[0].left_k, 346.6666666666667);
    close(result.layer_temperatures[0].right_k, 313.3333333333333);
}

#[test]
fn g0_invalid_geometry_resistance_and_varying_conductivity_refuse() {
    assert!(LayeredWall::new(vec![], vec![], [0.0; 2]).is_err());
    for thickness in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            PlanarThermalLayer::new("bad", thickness, ConductivityTable::declared(1.0).unwrap())
                .is_err()
        );
    }
    assert!(LayeredWall::new(vec![constant("a", 1.0, 1.0)], vec![0.0], [0.0; 2]).is_err());
    for resistance in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            LayeredWall::new(vec![constant("a", 1.0, 1.0)], vec![], [resistance, 0.0]).is_err()
        );
    }
    assert!(
        LayeredWall::new(
            vec![constant("same", 1.0, 1.0), constant("same", 1.0, 1.0)],
            vec![0.0],
            [0.0; 2],
        )
        .is_err()
    );
    let oversized = (0..=MAX_LAYERED_WALL_LAYERS)
        .map(|i| constant(&format!("layer {i}"), 1.0, 1.0))
        .collect();
    assert!(LayeredWall::new(oversized, vec![0.0; MAX_LAYERED_WALL_LAYERS], [0.0; 2]).is_err());
    let varying = PlanarThermalLayer::new(
        "varying",
        0.1,
        ConductivityTable::declared_curve(vec![(280.0, 1.0), (350.0, 2.0)]).unwrap(),
    )
    .unwrap();
    let stack = LayeredWall::new(vec![varying], vec![], [0.0; 2]).unwrap();
    assert!(matches!(
        with_cx(|cx| stack.solve(cx, [320.0, 300.0], 1.0)),
        Err(ConductionError::Conductivity { .. })
    ));
    for area in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(with_cx(|cx| wall().solve(cx, [320.0, 300.0], area)).is_err());
    }
    for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(with_cx(|cx| wall().solve(cx, [temperature, 300.0], 1.0)).is_err());
    }
}

#[test]
fn g4_cancellation_and_failed_trials_leave_the_input_unchanged() {
    let stack = wall();
    let before = stack.clone();
    assert!(matches!(
        with_cancelled_cx(|cx| stack.solve(cx, [340.0, 280.0], 1.0)),
        Err(ConductionError::Cancelled {
            stage: "layered-wall",
            ..
        })
    ));
    assert_eq!(stack, before);
    assert!(with_cx(|cx| stack.solve(cx, [0.0, 280.0], 1.0)).is_err());
    assert_eq!(stack, before);
    close(
        with_cx(|cx| stack.solve(cx, [340.0, 280.0], 1.0))
            .unwrap()
            .heat_rate_w,
        20.0,
    );
}
