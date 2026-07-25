//! Operating envelopes, duty cycles, and load-combination evaluation
//! (`f85xj.17.5`).
//!
//! The load-bearing case in this file is
//! [`correlated_corner_enumeration_matches_the_hand_worked_case`]: it fixes an
//! envelope whose arithmetic can be done on paper and pins the fact the bead
//! exists to establish — that the corners of a coupled envelope are not the
//! naive hypercube corners, and in this instance are not even the same
//! *number* of them.

use fs_qty::{Dims, QtyAny};
use fs_scenario::envelope::{
    AxisCoupling, AxisDomain, AxisPoint, CouplingRelation, DEFAULT_ENVELOPE_BUDGET, DiscreteState,
    DutyPoint, DutyWeighting, EnvelopeAxis, EnvelopeBudget, EnvelopeDutyCycle, EnvelopePoint,
    OperatingEnvelope, QoiSense, StateKind, compile_case_set, governing_case,
    reference_cooling_envelope, single_failure_axis,
};
use fs_scenario::scenario::{Combination, LoadCase, Violation};

const TEMP: Dims = Dims([0, 0, 0, 1, 0, 0]);
const POWER: Dims = Dims([2, 1, -3, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);

fn kelvin(value: f64) -> QtyAny {
    QtyAny::new(value, TEMP)
}

fn watts(value: f64) -> QtyAny {
    QtyAny::new(value, POWER)
}

fn seconds(value: f64) -> QtyAny {
    QtyAny::new(value, TIME)
}

fn codes(findings: &[Violation]) -> Vec<&str> {
    findings.iter().map(|finding| finding.code).collect()
}

/// The reference envelope's two continuous axes, without the fan axis, so a
/// corner count is the coupled block's vertex count alone.
fn coupled_pair(bound: f64) -> OperatingEnvelope {
    OperatingEnvelope {
        name: "coupled-pair".to_string(),
        axes: vec![
            EnvelopeAxis::continuous("ambient", kelvin(300.0), kelvin(320.0)),
            EnvelopeAxis::continuous("power", watts(100.0), watts(200.0)),
        ],
        couplings: vec![AxisCoupling {
            a: "ambient".to_string(),
            b: "power".to_string(),
            relation: CouplingRelation::CoOccurrenceLimit {
                a_coeff: QtyAny::new(2.5, POWER.checked_minus(TEMP).expect("in range")),
                b_coeff: QtyAny::dimensionless(1.0),
                bound: watts(bound),
                rationale: "controller throttling".to_string(),
            },
        }],
    }
}

fn uncoupled_pair(relation: CouplingRelation) -> OperatingEnvelope {
    OperatingEnvelope {
        name: "uncoupled-pair".to_string(),
        axes: vec![
            EnvelopeAxis::continuous("ambient", kelvin(300.0), kelvin(320.0)),
            EnvelopeAxis::continuous("power", watts(100.0), watts(200.0)),
        ],
        couplings: vec![AxisCoupling {
            a: "ambient".to_string(),
            b: "power".to_string(),
            relation,
        }],
    }
}

fn point(ambient: f64, power: f64) -> EnvelopePoint {
    EnvelopePoint {
        coordinates: vec![
            (
                "ambient".to_string(),
                AxisPoint::Continuous(kelvin(ambient)),
            ),
            ("power".to_string(), AxisPoint::Continuous(watts(power))),
        ],
    }
}

// ---------------------------------------------------------------------------
// Corner enumeration: the hand-worked case
// ---------------------------------------------------------------------------

/// The case the bead asks to be proven against hand-worked arithmetic.
///
/// Box `[300, 320] K x [100, 200] W`, limit `2.5*T + P <= 975 W`.
///
/// By hand, testing each box corner against the limit:
///
/// | corner        | `2.5T + P` | admitted |
/// |---------------|-----------|----------|
/// | (300, 100)    | 850       | yes      |
/// | (320, 100)    | 900       | yes      |
/// | (300, 200)    | 950       | yes      |
/// | (320, 200)    | 1000      | **no**   |
///
/// The limit line crosses the edge `T = 320` at `P = 975 - 800 = 175`, and the
/// edge `P = 200` at `T = 775 / 2.5 = 310`. So one box corner leaves and two
/// boundary vertices arrive: **five** vertices, where the naive hypercube has
/// four. Every number here is exactly representable in binary floating point,
/// so the assertion is on exact equality rather than a tolerance.
#[test]
fn correlated_corner_enumeration_matches_the_hand_worked_case() {
    let corners = coupled_pair(975.0)
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect("the coupled pair is a valid envelope");

    let attained: Vec<(f64, f64)> = corners
        .corners
        .iter()
        .map(|corner| {
            (
                corner.continuous("ambient").expect("ambient present").value,
                corner.continuous("power").expect("power present").value,
            )
        })
        .collect();

    assert_eq!(
        attained,
        vec![
            (300.0, 100.0),
            (300.0, 200.0),
            (310.0, 200.0),
            (320.0, 100.0),
            (320.0, 175.0),
        ],
        "hand-worked vertices of the clipped box"
    );
}

#[test]
fn a_coupled_envelope_has_more_corners_than_the_naive_hypercube_not_fewer() {
    let coupled = coupled_pair(975.0)
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect("valid");
    let naive = uncoupled_pair(CouplingRelation::DeclaredIndependent)
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect("valid");

    assert_eq!(naive.corners.len(), 4, "the independence hypercube");
    assert_eq!(coupled.corners.len(), 5, "the clipped envelope");
    assert!(
        coupled.corners.len() > naive.corners.len(),
        "clipping a corner off a rectangle ADDS a vertex; treating a coupled envelope as a \
         subset of the hypercube's corner LIST would miss the two that only exist on the limit \
         boundary"
    );
}

#[test]
fn the_excluded_hypercube_corner_is_absent_and_the_new_vertices_are_present() {
    let corners = coupled_pair(975.0)
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect("valid");
    let attained: Vec<(f64, f64)> = corners
        .corners
        .iter()
        .map(|corner| {
            (
                corner.continuous("ambient").expect("present").value,
                corner.continuous("power").expect("present").value,
            )
        })
        .collect();

    assert!(
        !attained.contains(&(320.0, 200.0)),
        "the hottest-and-highest-power hypercube corner is unreachable under the declared limit"
    );
    for arrival in [(320.0, 175.0), (310.0, 200.0)] {
        assert!(
            attained.contains(&arrival),
            "{arrival:?} lies on the limit boundary and is a corner of the FEASIBLE region, \
             though it is not a corner of the box"
        );
    }
}

#[test]
fn a_limit_passing_exactly_through_a_corner_does_not_duplicate_it() {
    // 2.5*300 + 200 = 950 exactly, so (300, 200) is both a box corner and an
    // edge crossing. It must be enumerated once.
    let corners = coupled_pair(950.0)
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect("valid");
    let touching = corners
        .corners
        .iter()
        .filter(|corner| {
            corner.continuous("ambient").expect("present").value == 300.0
                && corner.continuous("power").expect("present").value == 200.0
        })
        .count();
    assert_eq!(touching, 1, "the touched corner is enumerated exactly once");
    assert_eq!(corners.corners.len(), 4);
}

#[test]
fn a_closed_limit_admits_the_point_that_lies_on_it() {
    // The half-plane is closed, so equality is admissible. If it were treated
    // as strict, the worst attainable operating point would be silently
    // excluded from the corner set.
    coupled_pair(950.0)
        .admits(&point(300.0, 200.0))
        .expect("2.5*300 + 200 == 950 satisfies `<= 950`");
}

#[test]
fn corner_enumeration_logs_the_rationale_for_the_count() {
    let corners = coupled_pair(975.0)
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect("valid");
    let joined = corners.rationale.join("\n");
    assert!(
        joined.contains("box corners 4 -> 5 vertices"),
        "the block's own arithmetic is logged, not just its output: {joined}"
    );
    assert!(
        joined.contains("1 box corner(s) unreachable"),
        "how many corners left: {joined}"
    );
    assert!(
        joined.contains("2 vertex/ices introduced"),
        "how many arrived: {joined}"
    );
    assert!(
        joined.contains("controller throttling"),
        "the author's stated physical reason travels with the count: {joined}"
    );
}

#[test]
fn corner_enumeration_is_deterministic_across_repeated_calls() {
    let envelope = reference_cooling_envelope();
    let first = envelope
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect("valid");
    let second = envelope
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect("valid");
    assert_eq!(first, second);
}

#[test]
fn coordinates_follow_declared_axis_order_not_coupling_order() {
    // The coupling names `power` first, but the envelope declares `ambient`
    // first. Two envelopes differing only in that ordering must yield
    // comparable points.
    let mut swapped = coupled_pair(975.0);
    swapped.couplings[0] = AxisCoupling {
        a: "power".to_string(),
        b: "ambient".to_string(),
        relation: CouplingRelation::CoOccurrenceLimit {
            a_coeff: QtyAny::dimensionless(1.0),
            b_coeff: QtyAny::new(2.5, POWER.checked_minus(TEMP).expect("in range")),
            bound: watts(975.0),
            rationale: "controller throttling".to_string(),
        },
    };
    let corners = swapped
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect("valid");
    for corner in &corners.corners {
        let names: Vec<&str> = corner
            .coordinates
            .iter()
            .map(|(name, _)| name.as_str())
            .collect();
        assert_eq!(names, vec!["ambient", "power"]);
    }

    // Not merely the same SET of corners: the same list, in the same order.
    // A corner's index is what a compiled case cites as its provenance, so a
    // permutation here would silently renumber every case in every report
    // while changing nothing physical.
    let straight = coupled_pair(975.0)
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect("valid");
    assert_eq!(
        corners.corners, straight.corners,
        "the corner list is a function of the declared axes and the feasible region, not of the \
         order in which a coupling happens to name its two axes"
    );
}

// ---------------------------------------------------------------------------
// Declared-unknown co-occurrence
// ---------------------------------------------------------------------------

#[test]
fn declared_unknown_co_occurrence_drops_no_corner_and_raises_a_caveat() {
    let corners = uncoupled_pair(CouplingRelation::Unknown {
        rationale: "no field data pairs ambient against dissipation for this platform".to_string(),
    })
    .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
    .expect("valid");

    assert_eq!(
        corners.corners.len(),
        4,
        "ignorance about a pair is not evidence that any corner is unreachable"
    );
    let caveat = corners
        .caveats
        .iter()
        .find(|caveat| caveat.code == "envelope-coupling-unknown")
        .expect("the unknown pair raises a caveat");
    assert!(
        caveat.what.contains("no field data"),
        "the declared reason is carried"
    );
    assert!(
        caveat.consequence.contains("SUPERSET"),
        "the caveat states the DIRECTION of the error, which is what makes it usable: {}",
        caveat.consequence
    );
    assert!(
        caveat.consequence.contains("probability-weighted"),
        "and states what the corner set may NOT be used for: {}",
        caveat.consequence
    );
}

#[test]
fn declared_independence_raises_no_caveat_because_it_is_a_declaration_not_ignorance() {
    let corners = uncoupled_pair(CouplingRelation::DeclaredIndependent)
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect("valid");
    assert_eq!(corners.corners.len(), 4);
    assert!(
        corners.caveats.is_empty(),
        "declared independence and uncharacterised co-occurrence produce the SAME corner set; \
         only the caveat distinguishes them, which is the whole reason the distinction is typed"
    );
}

#[test]
fn an_unknown_co_occurrence_must_say_why_it_is_unknown() {
    let findings = uncoupled_pair(CouplingRelation::Unknown {
        rationale: "   ".to_string(),
    })
    .validate();
    assert_eq!(
        codes(&findings),
        vec!["envelope-coupling-unknown-unexplained"]
    );
}

#[test]
fn two_unknown_pairs_each_raise_their_own_caveat_despite_nested_axis_names() {
    // `power` is a substring of `total-power`: a caveat-dedup keyed on text
    // would suppress the second pair.
    let envelope = OperatingEnvelope {
        name: "nested-names".to_string(),
        axes: vec![
            EnvelopeAxis::continuous("power", watts(1.0), watts(2.0)),
            EnvelopeAxis::continuous("ambient", kelvin(300.0), kelvin(310.0)),
            EnvelopeAxis::continuous("total-power", watts(10.0), watts(20.0)),
            EnvelopeAxis::continuous(
                "altitude",
                QtyAny::new(0.0, Dims([1, 0, 0, 0, 0, 0])),
                QtyAny::new(3000.0, Dims([1, 0, 0, 0, 0, 0])),
            ),
        ],
        couplings: vec![
            AxisCoupling {
                a: "power".to_string(),
                b: "ambient".to_string(),
                relation: CouplingRelation::Unknown {
                    rationale: "first pair uncharacterised".to_string(),
                },
            },
            AxisCoupling {
                a: "total-power".to_string(),
                b: "altitude".to_string(),
                relation: CouplingRelation::Unknown {
                    rationale: "second pair uncharacterised".to_string(),
                },
            },
        ],
    };
    let corners = envelope
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect("valid");
    assert_eq!(
        corners
            .caveats
            .iter()
            .filter(|caveat| caveat.code == "envelope-coupling-unknown")
            .count(),
        2,
        "each uncharacterised pair reports separately"
    );
}

// ---------------------------------------------------------------------------
// Envelope validation
// ---------------------------------------------------------------------------

#[test]
fn the_reference_envelope_validates_and_yields_fifteen_corners() {
    let envelope = reference_cooling_envelope();
    assert_eq!(envelope.validate(), Vec::new());

    let corners = envelope
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect("valid");
    // 5 vertices of the clipped (ambient, power) block x 3 fan states.
    assert_eq!(corners.corners.len(), 15);
    assert!(corners.caveats.is_empty(), "every pair is characterised");
}

#[test]
fn a_pinned_axis_is_refused_because_it_is_a_constant_not_an_axis() {
    let envelope = OperatingEnvelope {
        name: "pinned".to_string(),
        axes: vec![EnvelopeAxis::continuous(
            "ambient",
            kelvin(300.0),
            kelvin(300.0),
        )],
        couplings: Vec::new(),
    };
    let findings = envelope.validate();
    assert_eq!(codes(&findings), vec!["envelope-axis-empty"]);
    assert!(
        findings[0].fix.contains("declared constant"),
        "the fix names where the value should live instead: {}",
        findings[0].fix
    );
}

#[test]
fn axis_bounds_with_mismatched_dimensions_are_refused() {
    let envelope = OperatingEnvelope {
        name: "mixed".to_string(),
        axes: vec![EnvelopeAxis::continuous(
            "ambient",
            kelvin(300.0),
            watts(320.0),
        )],
        couplings: Vec::new(),
    };
    assert_eq!(codes(&envelope.validate()), vec!["envelope-axis-dims"]);
}

#[test]
fn a_discrete_axis_must_declare_exactly_one_nominal_state() {
    for states in [
        vec![
            DiscreteState::failed("a-failed"),
            DiscreteState::failed("b-failed"),
        ],
        vec![DiscreteState::nominal("x"), DiscreteState::nominal("y")],
    ] {
        let envelope = OperatingEnvelope {
            name: "fans".to_string(),
            axes: vec![EnvelopeAxis::discrete("fan", states)],
            couplings: Vec::new(),
        };
        let findings = envelope.validate();
        assert!(
            codes(&findings).contains(&"envelope-state-nominal-count"),
            "got {:?}",
            codes(&findings)
        );
        assert!(
            findings.iter().any(|finding| finding
                .fix
                .contains("design condition or a fault condition")),
            "the fix explains what the nominal marker is FOR"
        );
    }
}

#[test]
fn a_limit_whose_coefficient_dimensions_do_not_close_is_refused() {
    let mut envelope = coupled_pair(975.0);
    envelope.couplings[0].relation = CouplingRelation::CoOccurrenceLimit {
        // Dimensionless coefficient on a temperature axis gives a K-dimensioned
        // term against a W-dimensioned bound.
        a_coeff: QtyAny::dimensionless(2.5),
        b_coeff: QtyAny::dimensionless(1.0),
        bound: watts(975.0),
        rationale: "throttling".to_string(),
    };
    assert_eq!(codes(&envelope.validate()), vec!["envelope-limit-dims"]);
}

#[test]
fn a_limit_over_a_discrete_axis_is_refused_with_the_alternative_named() {
    let mut envelope = reference_cooling_envelope();
    envelope.couplings = vec![AxisCoupling {
        a: "total-power".to_string(),
        b: "fan-state".to_string(),
        relation: CouplingRelation::CoOccurrenceLimit {
            a_coeff: QtyAny::dimensionless(1.0),
            b_coeff: QtyAny::dimensionless(1.0),
            bound: watts(200.0),
            rationale: "a fan failure caps dissipation".to_string(),
        },
    }];
    let findings = envelope.validate();
    assert_eq!(codes(&findings), vec!["envelope-limit-discrete-axis"]);
    assert!(
        findings[0].fix.contains("one envelope per discrete state"),
        "the refusal points at a representable alternative rather than just saying no"
    );
}

#[test]
fn a_limit_that_excludes_the_entire_box_is_refused() {
    let findings = coupled_pair(100.0).validate();
    assert_eq!(codes(&findings), vec!["envelope-limit-empty"]);
}

#[test]
fn a_limit_must_carry_the_reason_the_excluded_region_is_unreachable() {
    let mut envelope = coupled_pair(975.0);
    envelope.couplings[0].relation = CouplingRelation::CoOccurrenceLimit {
        a_coeff: QtyAny::new(2.5, POWER.checked_minus(TEMP).expect("in range")),
        b_coeff: QtyAny::dimensionless(1.0),
        bound: watts(975.0),
        rationale: String::new(),
    };
    let findings = envelope.validate();
    assert_eq!(codes(&findings), vec!["envelope-limit-unexplained"]);
    assert!(
        findings[0].fix.contains("unaudited assumption"),
        "the fix says why an unexplained limit is dangerous: it silently removes worst cases"
    );
}

#[test]
fn an_axis_in_two_couplings_is_refused_and_the_refusal_names_the_missing_machinery() {
    let mut envelope = reference_cooling_envelope();
    envelope.couplings.push(AxisCoupling {
        a: "total-power".to_string(),
        b: "fan-state".to_string(),
        relation: CouplingRelation::DeclaredIndependent,
    });
    let findings = envelope.validate();
    assert!(
        codes(&findings).contains(&"envelope-coupling-axis-shared"),
        "got {:?}",
        codes(&findings)
    );
    let shared = findings
        .iter()
        .find(|finding| finding.code == "envelope-coupling-axis-shared")
        .expect("present");
    assert!(
        shared.fix.contains("general polytope vertex enumeration"),
        "the refusal names exactly what is not implemented, so a reader can tell a limitation \
         from a bug: {}",
        shared.fix
    );
}

#[test]
fn a_coupling_naming_an_undeclared_axis_is_refused() {
    let mut envelope = coupled_pair(975.0);
    envelope.couplings[0].b = "altitude".to_string();
    assert!(codes(&envelope.validate()).contains(&"envelope-coupling-axis-missing"));
}

// ---------------------------------------------------------------------------
// Point admission
// ---------------------------------------------------------------------------

#[test]
fn a_point_outside_the_declared_limit_is_refused_with_the_rationale() {
    let findings = coupled_pair(975.0)
        .admits(&point(320.0, 200.0))
        .expect_err("2.5*320 + 200 = 1000 exceeds 975");
    assert_eq!(codes(&findings), vec!["envelope-point-limit-violated"]);
    assert!(
        findings[0].what.contains("controller throttling"),
        "the reader is told which declared physical limit they hit"
    );
}

#[test]
fn a_point_missing_a_coordinate_is_refused_as_an_undeclared_default() {
    let incomplete = EnvelopePoint {
        coordinates: vec![("ambient".to_string(), AxisPoint::Continuous(kelvin(310.0)))],
    };
    let findings = coupled_pair(975.0)
        .admits(&incomplete)
        .expect_err("power is unspecified");
    assert_eq!(codes(&findings), vec!["envelope-point-axis-missing"]);
    assert!(findings[0].fix.contains("undeclared default"));
}

#[test]
fn a_point_with_the_wrong_kind_or_range_or_dimensions_is_refused() {
    let envelope = reference_cooling_envelope();

    let wrong_kind = EnvelopePoint {
        coordinates: vec![
            (
                "ambient-temperature".to_string(),
                AxisPoint::Discrete("hot".to_string()),
            ),
            (
                "total-power".to_string(),
                AxisPoint::Continuous(watts(150.0)),
            ),
            (
                "fan-state".to_string(),
                AxisPoint::Discrete("both-running".to_string()),
            ),
        ],
    };
    assert!(
        codes(&envelope.admits(&wrong_kind).expect_err("kind mismatch"))
            .contains(&"envelope-point-kind")
    );

    let out_of_range = EnvelopePoint {
        coordinates: vec![
            (
                "ambient-temperature".to_string(),
                AxisPoint::Continuous(kelvin(400.0)),
            ),
            (
                "total-power".to_string(),
                AxisPoint::Continuous(watts(150.0)),
            ),
            (
                "fan-state".to_string(),
                AxisPoint::Discrete("both-running".to_string()),
            ),
        ],
    };
    assert!(
        codes(&envelope.admits(&out_of_range).expect_err("out of range"))
            .contains(&"envelope-point-out-of-range")
    );

    let unknown_state = EnvelopePoint {
        coordinates: vec![
            (
                "ambient-temperature".to_string(),
                AxisPoint::Continuous(kelvin(310.0)),
            ),
            (
                "total-power".to_string(),
                AxisPoint::Continuous(watts(150.0)),
            ),
            (
                "fan-state".to_string(),
                AxisPoint::Discrete("all-failed".to_string()),
            ),
        ],
    };
    assert!(
        codes(
            &envelope
                .admits(&unknown_state)
                .expect_err("undeclared state")
        )
        .contains(&"envelope-point-state-unknown")
    );
}

// ---------------------------------------------------------------------------
// Bounded failure combinatorics
// ---------------------------------------------------------------------------

#[test]
fn single_failure_states_are_linear_in_the_unit_count_not_exponential() {
    for count in [1usize, 2, 3, 8] {
        let names: Vec<String> = (0..count).map(|index| format!("fan-{index}")).collect();
        let refs: Vec<&str> = names.iter().map(String::as_str).collect();
        let axis = single_failure_axis("fan-state", "all-running", &refs).expect("distinct names");
        let AxisDomain::Discrete { states } = &axis.domain else {
            panic!("discrete axis");
        };
        assert_eq!(
            states.len(),
            count + 1,
            "n units give n+1 states, never 2^n: at n = 8 the exponential form would be 256"
        );
        assert_eq!(
            states
                .iter()
                .filter(|s| s.kind == StateKind::Nominal)
                .count(),
            1
        );
        assert_eq!(
            states
                .iter()
                .filter(|s| s.kind == StateKind::Failed)
                .count(),
            count
        );
    }
}

#[test]
fn a_duplicated_failure_unit_is_refused_before_the_states_collide() {
    let findings =
        single_failure_axis("fan-state", "both", &["fan-a", "fan-a"]).expect_err("duplicate unit");
    assert_eq!(codes(&findings), vec!["failure-unit-duplicate"]);
}

// ---------------------------------------------------------------------------
// Resource admission
// ---------------------------------------------------------------------------

#[test]
fn a_corner_product_above_the_budget_is_refused_before_allocation() {
    let axes: Vec<EnvelopeAxis> = (0..20)
        .map(|index| {
            EnvelopeAxis::continuous(format!("axis-{index}"), kelvin(300.0), kelvin(320.0))
        })
        .collect();
    let envelope = OperatingEnvelope {
        name: "wide".to_string(),
        axes,
        couplings: Vec::new(),
    };
    // 2^20 = 1_048_576 corners.
    let findings = envelope
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect_err("the product exceeds the admitted maximum");
    assert_eq!(codes(&findings), vec!["envelope-corner-overflow"]);
    assert!(
        findings[0].fix.contains("nothing was enumerated"),
        "the refusal states that no allocation happened, which is the point of checking first"
    );
}

#[test]
fn the_corner_budget_admits_exactly_at_its_limit_and_refuses_one_past_it() {
    let axes: Vec<EnvelopeAxis> = (0..4)
        .map(|index| {
            EnvelopeAxis::continuous(format!("axis-{index}"), kelvin(300.0), kelvin(320.0))
        })
        .collect();
    let envelope = OperatingEnvelope {
        name: "sixteen".to_string(),
        axes,
        couplings: Vec::new(),
    };

    let admitted = envelope
        .enumerate_corners(EnvelopeBudget {
            max_corners: 16,
            ..DEFAULT_ENVELOPE_BUDGET
        })
        .expect("2^4 == 16 is at the limit, not past it");
    assert_eq!(admitted.corners.len(), 16);

    let refused = envelope
        .enumerate_corners(EnvelopeBudget {
            max_corners: 15,
            ..DEFAULT_ENVELOPE_BUDGET
        })
        .expect_err("16 exceeds 15");
    assert_eq!(codes(&refused), vec!["envelope-corner-overflow"]);
}

// ---------------------------------------------------------------------------
// Duty cycles
// ---------------------------------------------------------------------------

fn dwell_cycle() -> EnvelopeDutyCycle {
    EnvelopeDutyCycle {
        name: "office-day".to_string(),
        weighting: DutyWeighting::Dwells,
        points: vec![
            DutyPoint {
                name: "idle".to_string(),
                point: point(300.0, 100.0),
                weight: seconds(3600.0),
            },
            DutyPoint {
                name: "burst".to_string(),
                point: point(310.0, 200.0),
                weight: seconds(1200.0),
            },
        ],
    }
}

#[test]
fn dwell_weighted_fractions_are_derived_and_sum_to_one() {
    let cycle = dwell_cycle();
    cycle
        .validate_against(&coupled_pair(975.0))
        .expect("both dwells lie in the envelope");
    let fractions = cycle.fractions().expect("normalisable");
    assert_eq!(fractions, vec![0.75, 0.25]);
    assert_eq!(fractions.iter().sum::<f64>(), 1.0);
}

#[test]
fn the_weighted_aggregate_is_the_dwell_weighted_mean() {
    let aggregate = dwell_cycle()
        .weighted_aggregate(&[320.0, 360.0])
        .expect("two values for two dwells");
    // 0.75*320 + 0.25*360 = 240 + 90 = 330, exact in binary floating point.
    assert_eq!(aggregate.value, 330.0);
    assert_eq!(aggregate.fractions, vec![0.75, 0.25]);
}

#[test]
fn every_weighted_aggregate_carries_the_steady_approximation_caveat() {
    let aggregate = dwell_cycle()
        .weighted_aggregate(&[320.0, 360.0])
        .expect("ok");
    let caveat = aggregate
        .caveats
        .iter()
        .find(|caveat| caveat.code == "duty-steady-approximation")
        .expect("always present");
    assert!(
        caveat.consequence.contains("thermal time constant"),
        "the caveat names the condition under which the number is the true time average"
    );
    assert!(
        caveat.consequence.contains("fs_conduction::duty"),
        "and points at the rung that does integrate: {}",
        caveat.consequence
    );
}

#[test]
fn the_quasi_steady_condition_is_measured_not_asserted() {
    let cycle = dwell_cycle();

    // Time constant 120 s: the shortest dwell is 1200 s, ten constants.
    let comfortable = cycle
        .quasi_steady_check(seconds(120.0), 5.0)
        .expect("dwell weighting carries absolute time");
    assert_eq!(comfortable.min_dwell_ratio, 10.0);
    assert_eq!(comfortable.shortest_point, "burst");
    assert!(comfortable.satisfied);

    // Time constant 600 s: the shortest dwell is only two constants, so the
    // aggregate's premise fails and the verdict says which point breaks it.
    let marginal = cycle.quasi_steady_check(seconds(600.0), 5.0).expect("ok");
    assert_eq!(marginal.min_dwell_ratio, 2.0);
    assert_eq!(marginal.shortest_point, "burst");
    assert!(
        !marginal.satisfied,
        "a two-time-constant dwell has settled to about 86 percent; calling that steady is the \
         error the check exists to expose"
    );
}

#[test]
fn fraction_weighting_makes_the_quasi_steady_condition_unmeasurable_not_merely_unmet() {
    let cycle = EnvelopeDutyCycle {
        name: "fractions-only".to_string(),
        weighting: DutyWeighting::Fractions,
        points: vec![
            DutyPoint {
                name: "idle".to_string(),
                point: point(300.0, 100.0),
                weight: QtyAny::dimensionless(0.75),
            },
            DutyPoint {
                name: "burst".to_string(),
                point: point(310.0, 200.0),
                weight: QtyAny::dimensionless(0.25),
            },
        ],
    };
    cycle
        .validate_against(&coupled_pair(975.0))
        .expect("fractions sum to one and both points are admissible");

    let findings = cycle
        .quasi_steady_check(seconds(120.0), 5.0)
        .expect_err("no absolute time exists to compare against a time constant");
    assert_eq!(codes(&findings), vec!["duty-quasi-steady-untimed"]);

    let aggregate = cycle.weighted_aggregate(&[320.0, 360.0]).expect("ok");
    assert_eq!(aggregate.value, 330.0, "the same mean as the dwell form");
    let caveat = aggregate
        .caveats
        .iter()
        .find(|caveat| caveat.code == "duty-no-absolute-time")
        .expect("the untimed form says so");
    assert!(
        caveat.consequence.contains("unmeasurable"),
        "the distinction that matters: not that the premise fails, but that nobody can tell: {}",
        caveat.consequence
    );
}

#[test]
fn declared_fractions_that_do_not_sum_to_one_are_refused() {
    let cycle = EnvelopeDutyCycle {
        name: "unnormalised".to_string(),
        weighting: DutyWeighting::Fractions,
        points: vec![
            DutyPoint {
                name: "idle".to_string(),
                point: point(300.0, 100.0),
                weight: QtyAny::dimensionless(0.75),
            },
            DutyPoint {
                name: "burst".to_string(),
                point: point(310.0, 200.0),
                weight: QtyAny::dimensionless(0.5),
            },
        ],
    };
    let findings = cycle
        .validate_against(&coupled_pair(975.0))
        .expect_err("1.25 is not one");
    assert_eq!(codes(&findings), vec!["duty-fractions-sum"]);
}

#[test]
fn a_dwell_weight_in_the_wrong_dimensions_for_the_declared_weighting_is_refused() {
    let mut cycle = dwell_cycle();
    cycle.points[0].weight = QtyAny::dimensionless(0.75);
    let findings = cycle
        .validate_against(&coupled_pair(975.0))
        .expect_err("Dwells weighting requires seconds");
    assert!(codes(&findings).contains(&"duty-weight-dims"));
}

#[test]
fn a_duty_point_outside_the_envelope_is_refused_and_names_the_point() {
    let mut cycle = dwell_cycle();
    cycle.points[1].point = point(320.0, 200.0);
    let findings = cycle
        .validate_against(&coupled_pair(975.0))
        .expect_err("the burst point violates the throttle limit");
    assert!(codes(&findings).contains(&"envelope-point-limit-violated"));
    assert!(
        findings
            .iter()
            .any(|finding| finding.what.starts_with("duty point `burst`:")),
        "the finding is attributed to the dwell, not just to an anonymous point"
    );
}

#[test]
fn a_zero_weight_dwell_is_refused_rather_than_silently_ignored() {
    let mut cycle = dwell_cycle();
    cycle.points[1].weight = seconds(0.0);
    let findings = cycle
        .validate_against(&coupled_pair(975.0))
        .expect_err("zero dwell");
    assert!(codes(&findings).contains(&"duty-weight-range"));
    assert!(
        findings
            .iter()
            .any(|finding| finding.fix.contains("removed, not carried")),
        "a zero-weight point that stays in the declaration reads as an operating condition that \
         was considered; it was not"
    );
}

#[test]
fn a_qoi_vector_of_the_wrong_length_is_refused() {
    let findings = dwell_cycle()
        .weighted_aggregate(&[320.0])
        .expect_err("two dwells, one value");
    assert_eq!(codes(&findings), vec!["duty-qoi-length"]);
}

#[test]
fn a_non_finite_qoi_is_refused_rather_than_aggregated() {
    let findings = dwell_cycle()
        .weighted_aggregate(&[320.0, f64::NAN])
        .expect_err("a failed solve is not a value");
    assert_eq!(codes(&findings), vec!["duty-qoi-nonfinite"]);
}

// ---------------------------------------------------------------------------
// Combination compilation and governing cases
// ---------------------------------------------------------------------------

fn cases() -> Vec<LoadCase> {
    ["dead", "live", "solar"]
        .into_iter()
        .map(|name| LoadCase {
            name: name.to_string(),
            bcs: Vec::new(),
        })
        .collect()
}

fn combinations() -> Vec<Combination> {
    vec![
        Combination {
            name: "service".to_string(),
            terms: vec![("dead".to_string(), 1.0), ("live".to_string(), 1.0)],
        },
        Combination {
            name: "hot-day".to_string(),
            terms: vec![
                ("dead".to_string(), 1.2),
                ("live".to_string(), 1.6),
                ("solar".to_string(), 1.0),
            ],
        },
    ]
}

#[test]
fn compiled_cases_carry_the_combination_and_its_factors_as_provenance() {
    let set = compile_case_set(&cases(), &combinations()).expect("valid declarations");
    assert_eq!(set.cases.len(), 2);
    assert_eq!(set.cases[1].combination, "hot-day");
    assert_eq!(
        set.cases[1].terms,
        vec![
            ("dead".to_string(), 1.2),
            ("live".to_string(), 1.6),
            ("solar".to_string(), 1.0),
        ],
        "a report can say WHICH factors produced the governing case, not merely its name"
    );
    assert!(set.cases[1].envelope_point.is_none());
}

#[test]
fn a_combination_referencing_an_undeclared_case_is_refused() {
    let mut combos = combinations();
    combos[0].terms[1] = ("wind".to_string(), 1.0);
    let findings = compile_case_set(&cases(), &combos).expect_err("wind is not declared");
    assert_eq!(codes(&findings), vec!["combo-case-missing"]);
}

#[test]
fn a_repeated_case_within_one_combination_is_refused_rather_than_summed() {
    let mut combos = combinations();
    combos[0].terms.push(("dead".to_string(), 0.5));
    let findings = compile_case_set(&cases(), &combos).expect_err("dead appears twice");
    assert_eq!(codes(&findings), vec!["combo-term-duplicate"]);
    assert!(
        findings[0].fix.contains("invisible at the declaration"),
        "silently summing 1.0 and 0.5 would make the effective 1.5 appear nowhere the author \
         wrote it"
    );
}

#[test]
fn duplicate_combination_names_are_refused_because_a_report_could_not_tell_them_apart() {
    let mut combos = combinations();
    combos[1].name = "service".to_string();
    let findings = compile_case_set(&cases(), &combos).expect_err("duplicate name");
    assert_eq!(codes(&findings), vec!["combo-name-duplicate"]);
}

#[test]
fn crossing_a_case_set_with_envelope_points_preserves_both_provenances() {
    let set = compile_case_set(&cases(), &combinations()).expect("valid");
    let corners = reference_cooling_envelope()
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect("valid");
    let named: Vec<(String, EnvelopePoint)> = corners
        .corners
        .iter()
        .enumerate()
        .map(|(index, corner)| (format!("corner-{index}"), corner.clone()))
        .collect();

    let crossed = set
        .cross_with_points(&named, DEFAULT_ENVELOPE_BUDGET)
        .expect("2 x 15 is within budget");
    assert_eq!(crossed.cases.len(), 30);
    assert_eq!(crossed.cases[0].combination, "service");
    assert_eq!(crossed.cases[0].envelope_point.as_deref(), Some("corner-0"));
    assert_eq!(crossed.cases[29].combination, "hot-day");
    assert_eq!(
        crossed.cases[29].envelope_point.as_deref(),
        Some("corner-14")
    );
}

#[test]
fn a_case_product_above_the_budget_is_refused_before_allocation() {
    let set = compile_case_set(&cases(), &combinations()).expect("valid");
    let named: Vec<(String, EnvelopePoint)> = (0..100)
        .map(|index| (format!("p{index}"), point(300.0, 100.0)))
        .collect();
    let budget = EnvelopeBudget {
        max_cases: 50,
        ..DEFAULT_ENVELOPE_BUDGET
    };
    let findings = set
        .cross_with_points(&named, budget)
        .expect_err("2 x 100 exceeds 50");
    assert_eq!(codes(&findings), vec!["envelope-case-overflow"]);
    assert!(findings[0].what.contains("200 case(s)"));
}

#[test]
fn duplicate_envelope_point_names_are_refused() {
    let set = compile_case_set(&cases(), &combinations()).expect("valid");
    let named = vec![
        ("corner".to_string(), point(300.0, 100.0)),
        ("corner".to_string(), point(320.0, 100.0)),
    ];
    let findings = set
        .cross_with_points(&named, DEFAULT_ENVELOPE_BUDGET)
        .expect_err("duplicate point name");
    assert_eq!(codes(&findings), vec!["envelope-case-point-duplicate"]);
}

#[test]
fn the_governing_case_is_the_worst_in_the_declared_sense() {
    let hottest = governing_case(
        "T_junction",
        QoiSense::LargerIsWorse,
        &[350.0, 372.0, 361.0],
    )
    .expect("finite values");
    assert_eq!(hottest.value, 372.0);
    assert_eq!(hottest.governing, vec![1]);

    let thinnest =
        governing_case("margin", QoiSense::SmallerIsWorse, &[12.0, 3.0, 7.0]).expect("finite");
    assert_eq!(thinnest.value, 3.0);
    assert_eq!(thinnest.governing, vec![1]);
}

#[test]
fn a_tie_reports_every_governing_case_rather_than_the_first() {
    let tied = governing_case(
        "T_junction",
        QoiSense::LargerIsWorse,
        &[372.0, 350.0, 372.0],
    )
    .expect("finite");
    assert_eq!(tied.governing, vec![0, 2]);
    assert_eq!(
        tied.governing.len(),
        2,
        "reporting only the first would hide that two conditions size the design equally, which \
         changes what a reviewer has to check"
    );
}

#[test]
fn a_non_finite_case_value_is_refused_so_a_nan_cannot_quietly_never_govern() {
    let findings = governing_case(
        "T_junction",
        QoiSense::LargerIsWorse,
        &[350.0, f64::NAN, 361.0],
    )
    .expect_err("a failed case is not a value");
    assert_eq!(codes(&findings), vec!["governing-nonfinite"]);
    // The hazard being closed: NaN is unordered against everything, so a NaN
    // case would lose every contest and be reported as the benign one.
    assert!(f64::NAN.partial_cmp(&350.0).is_none());
}

#[test]
fn asking_which_case_governs_with_no_cases_is_refused() {
    let findings =
        governing_case("T_junction", QoiSense::LargerIsWorse, &[]).expect_err("no cases");
    assert_eq!(codes(&findings), vec!["governing-no-cases"]);
}

// ---------------------------------------------------------------------------
// End to end
// ---------------------------------------------------------------------------

/// Reference envelope -> corner set -> case set -> governing case per QoI.
///
/// The QoI model here is deliberately transparent: junction temperature rises
/// with ambient and with power, and a failed fan adds a fixed penalty. It is a
/// stand-in for a solve, not a solve — the point of the test is that the
/// PLUMBING identifies the right governing condition and reports it with its
/// provenance intact, including whether the governing condition is a design
/// condition or a fault condition.
#[test]
fn envelope_to_governing_case_runs_end_to_end_and_distinguishes_fault_conditions() {
    let envelope = reference_cooling_envelope();
    let corners = envelope
        .enumerate_corners(DEFAULT_ENVELOPE_BUDGET)
        .expect("the reference envelope is valid");
    assert_eq!(corners.corners.len(), 15);

    let named: Vec<(String, EnvelopePoint)> = corners
        .corners
        .iter()
        .map(|corner| {
            let ambient = corner
                .continuous("ambient-temperature")
                .expect("present")
                .value;
            let power = corner.continuous("total-power").expect("present").value;
            let fan = corner.discrete("fan-state").expect("present");
            (format!("{ambient:.0}K-{power:.0}W-{fan}"), corner.clone())
        })
        .collect();

    let set = compile_case_set(
        &[LoadCase {
            name: "steady".to_string(),
            bcs: Vec::new(),
        }],
        &[Combination {
            name: "nominal-service".to_string(),
            terms: vec![("steady".to_string(), 1.0)],
        }],
    )
    .expect("valid declarations");
    let crossed = set
        .cross_with_points(&named, DEFAULT_ENVELOPE_BUDGET)
        .expect("1 x 15 is within budget");
    assert_eq!(crossed.cases.len(), 15);

    // Stand-in QoI: T_j = ambient + power * resistance, resistance worsening
    // when a fan has failed.
    let corner_for = |label: &str| -> &EnvelopePoint {
        named
            .iter()
            .find(|(name, _)| name.as_str() == label)
            .map(|(_, corner)| corner)
            .expect("every crossed case cites a named corner")
    };

    let junction: Vec<f64> = crossed
        .cases
        .iter()
        .map(|case| {
            let corner = corner_for(case.envelope_point.as_deref().expect("crossed"));
            let ambient = corner
                .continuous("ambient-temperature")
                .expect("present")
                .value;
            let power = corner.continuous("total-power").expect("present").value;
            let failed = corner.discrete("fan-state").expect("present") != "both-running";
            let resistance = if failed { 0.35 } else { 0.20 };
            ambient + power * resistance
        })
        .collect();

    let governing =
        governing_case("T_junction", QoiSense::LargerIsWorse, &junction).expect("finite values");

    // 320 K with 175 W and a failed fan: 320 + 175*0.35 = 381.25. The runner-up
    // reachable corner is 310 K at the full 200 W, 310 + 70 = 380.
    assert_eq!(governing.value, 381.25);

    // The two fans are interchangeable in this stand-in model, so `fan-a-failed`
    // and `fan-b-failed` govern EQUALLY. That is the physically right answer and
    // exactly the case a first-match report would have hidden.
    assert_eq!(
        governing.governing.len(),
        2,
        "symmetric redundancy produces a genuine tie between the two single-fan failures"
    );

    let mut governing_fans = Vec::new();
    for index in &governing.governing {
        let winner = &crossed.cases[*index];
        let corner = corner_for(winner.envelope_point.as_deref().expect("crossed"));

        let ambient = corner
            .continuous("ambient-temperature")
            .expect("present")
            .value;
        let power = corner.continuous("total-power").expect("present").value;
        assert_eq!(
            (ambient, power),
            (320.0, 175.0),
            "the naive worst case (320 K, 200 W) is unreachable, so a hypercube sweep would have \
             sized the design against a condition the controller prevents"
        );

        // The governing state is a declared FAILURE, which the report must be
        // able to say: sizing against a fault condition is a different
        // engineering decision from sizing against a design condition.
        let fan = corner.discrete("fan-state").expect("present");
        let AxisDomain::Discrete { states } = &envelope.axis("fan-state").expect("declared").domain
        else {
            panic!("discrete axis");
        };
        let kind = states
            .iter()
            .find(|state| state.name == fan)
            .expect("declared state")
            .kind;
        assert_eq!(kind, StateKind::Failed);
        governing_fans.push(fan.to_string());

        // And the case still cites the combination it came from.
        assert_eq!(winner.combination, "nominal-service");
        assert_eq!(winner.terms, vec![("steady".to_string(), 1.0)]);
    }
    governing_fans.sort();
    assert_eq!(governing_fans, vec!["fan-a-failed", "fan-b-failed"]);
}
