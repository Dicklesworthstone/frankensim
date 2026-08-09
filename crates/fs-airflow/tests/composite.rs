//! Analytic battery for composite fan systems (bead frn2i.1): series
//! pressure addition, parallel flow addition, heterogeneous composition,
//! fan-law folding, canonical ordering, refusal classes, and a certified
//! operating point through the production solve.

use fs_airflow::{
    AirflowError, EnclosureNetwork, FanArrangement, FanBank, FanCurve, FanPoint, LeakageElement,
    LossElement, LossNetwork, LossResistance, SourceProvenance, ToleranceBasis,
    composite::{compose_parallel, compose_series},
};
use fs_qty::{Pressure, VolumetricFlowRate};

fn flow(value: f64) -> VolumetricFlowRate {
    VolumetricFlowRate::new(value)
}

fn curve(name: &str, knots: &[(f64, f64)], admissible: f64) -> FanCurve {
    FanCurve::new(
        name,
        knots
            .iter()
            .map(|&(q, p)| FanPoint::new(flow(q), Pressure::new(p)))
            .collect(),
        SourceProvenance::new(format!("{name} citation"), format!("{name}-source-id")),
        0.05,
        ToleranceBasis::ManufacturerDeclared,
        flow(admissible),
        (0.5, 2.0),
    )
    .expect("fixture curve")
}

fn bank(curve: FanCurve, count: usize, arrangement: FanArrangement, speed_ratio: f64) -> FanBank {
    FanBank::new(curve, count, arrangement, speed_ratio).expect("fixture bank")
}

fn point_map(bank: &FanBank) -> Vec<(f64, f64)> {
    bank.curve()
        .points()
        .iter()
        .map(|point| (point.flow.value(), point.pressure.value()))
        .collect()
}

#[test]
fn series_composition_adds_pressures_exactly() {
    // A: p = 100 - 1000 Q on [0.01, 0.10]; B: p = 60 - 1000 Q on [0.005, 0.06].
    let a = bank(
        curve("fan-a", &[(0.0, 100.0), (0.10, 0.0)], 0.01),
        1,
        FanArrangement::Series,
        1.0,
    );
    let b = bank(
        curve("fan-b", &[(0.0, 60.0), (0.06, 0.0)], 0.005),
        1,
        FanArrangement::Series,
        1.0,
    );
    let composite = compose_series(&[a, b]).expect("series composite");
    let points = point_map(&composite);
    // Shared domain [0.01, 0.06]; composite p(Q) = 160 - 2000 Q.
    let &(q0, p0) = points.first().expect("first");
    assert!((q0 - 0.01).abs() <= 1e-15);
    assert!((p0 - 140.0).abs() <= 1e-9, "p(0.01) = 140, got {p0}");
    let &(qn, pn) = points.last().expect("last");
    assert!((qn - 0.06).abs() <= 1e-15);
    assert!((pn - 40.0).abs() <= 1e-9, "p(0.06) = 40, got {pn}");
    // Interior knots at the member breakpoints 0.06? 0.10 is outside; the
    // only interior knot is none — both members are single-segment, so the
    // composite has exactly two knots.
    assert_eq!(points.len(), 2);
    // Monotone pressure, strictly increasing flow, pinned speed domain.
    assert!(p0 > pn);
    assert_eq!(composite.speed_ratio(), 1.0);
}

#[test]
fn parallel_composition_adds_flows_exactly() {
    // A: p = 100 - 1000 Q; B: p = 50 - 500 Q (B shuts off at p = 50).
    let a = bank(
        curve("fan-a", &[(0.0, 100.0), (0.10, 0.0)], 0.01),
        1,
        FanArrangement::Series,
        1.0,
    );
    let b = bank(
        curve("fan-b", &[(0.0, 50.0), (0.10, 0.0)], 0.005),
        1,
        FanArrangement::Series,
        1.0,
    );
    let composite = compose_parallel(&[a, b]).expect("parallel composite");
    let points = point_map(&composite);
    // Stall pressures: A hits 0.01 flow at p = 90; B hits 0.005 at p = 47.5.
    // Shared pressure domain [0, 47.5]; total flow at 47.5 is
    // qA(47.5) + 0.005 = 0.0525 + 0.005 = 0.0575; at p = 0: 0.2.
    let &(q_min, p_top) = points.first().expect("first");
    assert!((q_min - 0.0575).abs() <= 1e-12, "min total flow {q_min}");
    assert!((p_top - 47.5).abs() <= 1e-9);
    let &(q_max, p_zero) = points.last().expect("last");
    assert!((q_max - 0.2).abs() <= 1e-12);
    assert!(p_zero.abs() <= 1e-12);
    // Interior knot at B's shutoff pressure 50? No — 50 > 47.5 lies outside
    // the shared domain, so the composite has exactly the two endpoints.
    assert_eq!(points.len(), 2);
}

#[test]
fn parallel_three_member_knots_merge_exactly() {
    let a = bank(
        curve("fan-a", &[(0.0, 100.0), (0.05, 60.0), (0.10, 0.0)], 0.0),
        1,
        FanArrangement::Series,
        1.0,
    );
    let b = bank(
        curve("fan-b", &[(0.0, 80.0), (0.08, 0.0)], 0.0),
        1,
        FanArrangement::Series,
        1.0,
    );
    let composite = compose_parallel(&[a, b]).expect("composite");
    let points = point_map(&composite);
    // Shared pressure domain [0, 80] (both admissible at 0 flow). Interior
    // knots at p = 60 (A's knee) with qA(60) = 0.05, qB(60) = 0.02 ->
    // total 0.07.
    let knee = points
        .iter()
        .find(|&&(_, p)| (p - 60.0).abs() <= 1e-9)
        .expect("knee at p=60");
    assert!((knee.0 - 0.07).abs() <= 1e-12, "knee flow {}", knee.0);
    assert!(points.len() >= 3, "merged knots retained");
}

#[test]
fn fan_law_scaling_folds_into_member_curves() {
    let base = curve("fan-a", &[(0.0, 100.0), (0.10, 0.0)], 0.01);
    let half_speed = bank(base.clone(), 1, FanArrangement::Series, 0.5);
    let full_speed = bank(base, 1, FanArrangement::Series, 1.0);
    let composite = compose_series(&[half_speed, full_speed]).expect("composite");
    let points = point_map(&composite);
    // Half-speed member: q in [0.005, 0.05], p = 25 - 250 Q_h... effective
    // curve (0,25)-(0.05,0). Full: (0,100)-(0.10,0). Shared domain
    // [0.01, 0.05]; composite p(Q) = (25 - 500 Q) + (100 - 1000 Q) =
    // 125 - 1500 Q? At Q=0.01: half member p = 25 - 500*0.01 = 20;
    // full member p = 90; total 110.
    let &(q0, p0) = points.first().expect("first");
    assert!((q0 - 0.01).abs() <= 1e-15);
    assert!((p0 - 110.0).abs() <= 1e-9, "folded speed composite {p0}");
}

#[test]
fn composition_refusals_are_typed() {
    let single = bank(
        curve("fan-a", &[(0.0, 100.0), (0.10, 0.0)], 0.01),
        1,
        FanArrangement::Series,
        1.0,
    );
    assert!(matches!(
        compose_series(&[single.clone()]),
        Err(AirflowError::EmptyFanComposition { topology: "series" })
    ));
    assert!(matches!(
        compose_parallel(&[single]),
        Err(AirflowError::EmptyFanComposition {
            topology: "parallel"
        })
    ));
    // Disjoint flow domains refuse: A tops out at 0.10, B starts at 0.20.
    let a = bank(
        curve("fan-a", &[(0.0, 100.0), (0.10, 0.0)], 0.01),
        1,
        FanArrangement::Series,
        1.0,
    );
    let b = bank(
        curve("fan-b", &[(0.20, 80.0), (0.30, 0.0)], 0.20),
        1,
        FanArrangement::Series,
        1.0,
    );
    assert!(matches!(
        compose_series(&[a, b]),
        Err(AirflowError::NoCommonSeriesDomain { .. })
    ));
}

#[test]
fn member_order_cannot_move_the_composite() {
    let a = bank(
        curve("fan-a", &[(0.0, 100.0), (0.10, 0.0)], 0.01),
        1,
        FanArrangement::Series,
        1.0,
    );
    let b = bank(
        curve("fan-b", &[(0.0, 60.0), (0.06, 0.0)], 0.005),
        2,
        FanArrangement::Parallel,
        1.2,
    );
    let forward = compose_series(&[a.clone(), b.clone()]).expect("forward");
    let reversed = compose_series(&[b, a]).expect("reversed");
    assert_eq!(point_map(&forward), point_map(&reversed));
    assert_eq!(forward.curve().name(), reversed.curve().name());
    assert_eq!(
        forward.curve().source().identifier,
        reversed.curve().source().identifier
    );
}

#[test]
fn composite_solve_produces_a_certified_operating_point() {
    let a = bank(
        curve("fan-a", &[(0.0, 100.0), (0.10, 0.0)], 0.01),
        1,
        FanArrangement::Series,
        1.0,
    );
    let b = bank(
        curve("fan-b", &[(0.0, 60.0), (0.06, 0.0)], 0.005),
        1,
        FanArrangement::Series,
        1.0,
    );
    let composite = compose_series(&[a, b]).expect("composite");
    // Composite p(Q) = 160 - 2000 Q. Choose resistance so the root is
    // exactly Q = 0.04: R Q^2 = 80 -> R = 50000.
    let element = LossElement::new(
        "test-duct",
        LossResistance::new(50_000.0),
        0.0,
        SourceProvenance::new("analytic test duct", "test-duct-source"),
        ToleranceBasis::Analytic,
    )
    .expect("loss element");
    let network = EnclosureNetwork::new(
        LossNetwork::Element(element),
        LeakageElement::new(
            LossElement::new(
                "test-leak",
                LossResistance::new(1.0e12),
                0.0,
                SourceProvenance::new("analytic test leak", "test-leak-source"),
                ToleranceBasis::Analytic,
            )
            .expect("leak element"),
        ),
    );
    let point =
        fs_airflow::solve_operating_point(&composite, &network).expect("operating point solves");
    let flow_lo = point.nominal_root.flow.lo();
    let flow_hi = point.nominal_root.flow.hi();
    // Composite p(Q) = 160 - 2000 Q; the network's nominal equivalent
    // resistance (primary in parallel with the tight leak) sets the
    // quadratic loss, so the analytic root solves R Q^2 = 160 - 2000 Q.
    let resistance = network.equivalent_resistance().value();
    let expected = (-2000.0 + (2000.0_f64 * 2000.0 + 4.0 * resistance * 160.0).sqrt())
        / (2.0 * resistance);
    assert!(
        flow_lo <= expected && flow_hi >= expected,
        "certified root [{flow_lo}, {flow_hi}] must contain the analytic {expected}"
    );
    let width = flow_hi - flow_lo;
    assert!(width <= 1e-6, "interval root width {width} is tight");
}
