//! G0/G1/G3 batteries for swept implicit and envelope honesty (bead c58q).
//!
//! The branch-and-bound partition and interval bounds are the proof.  Closed
//! forms and sampled pointwise motors below are independent falsifiers only.

use fs_evidence::NumericalKind;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ga::{Motor, Point as GaPoint};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Chart, Point3, TraceStepClaim};
use fs_ivl::Interval;
use fs_math::det;
use fs_motion::{
    EnvelopeBranchClass, EnvelopeConfig, EnvelopeDecision, EnvelopeEvidence, EnvelopeOracle,
    MotionError, ProofState, ScrewParams, SpacetimeChart, SweepDecision, SweptChart, SweptConfig,
    WankelApexPoint, WankelParams, WankelSealCircle, classify_envelope_branch, envelope,
    screw_tube, wankel_tube,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xC58_0001,
                kernel_id: 0x58,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn translating_sphere() -> SpacetimeChart<SphereChart> {
    let tube = screw_tube(
        &ScrewParams {
            axis: [1.0, 0.0, 0.0],
            center: [0.0, 0.0, 0.0],
            omega: 0.0,
            axial_velocity: 1.0,
            base_pose: Motor::identity(),
        },
        Interval::new(0.0, 1.0),
        8,
        4,
    )
    .expect("analytic translation tube");
    SpacetimeChart::new(
        SphereChart {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 0.25,
        },
        tube,
    )
}

fn capsule_distance(point: Point3) -> f64 {
    let nearest_x = point.x.clamp(0.0, 1.0);
    let dx = point.x - nearest_x;
    det::sqrt(dx * dx + point.y * point.y + point.z * point.z) - 0.25
}

#[test]
fn po_5_exhaustive_time_partition_encloses_translation_capsule_oracle() {
    with_cx(|cx| {
        let swept = SweptChart::new(
            translating_sphere(),
            SweptConfig {
                value_tolerance: 2.0e-2,
                time_tolerance: 1.0e-5,
                max_subdivisions: 256,
            },
            cx,
        )
        .expect("finite swept support");
        let points = [
            Point3::new(-0.4, 0.0, 0.0),
            Point3::new(0.1, 0.3, 0.0),
            Point3::new(0.5, 0.0, 0.0),
            Point3::new(0.85, -0.2, 0.1),
            Point3::new(1.4, 0.0, 0.0),
        ];
        for point in points {
            let receipt = swept.evaluate(point, cx).expect("certified infimum band");
            let oracle = capsule_distance(point);
            assert!(
                receipt.implicit.contains(oracle),
                "closed-form value {oracle} escaped [{}, {}] at {point:?}",
                receipt.implicit.lo(),
                receipt.implicit.hi()
            );
            assert_eq!(
                receipt.evaluated_cells,
                1 + 2 * receipt.subdivisions,
                "every binary split must replace one cell with two certified children"
            );
            assert_eq!(receipt.evaluated_witnesses, receipt.evaluated_cells);
            assert_eq!(receipt.decision, SweepDecision::Enclosure);
            println!(
                "{{\"case\":\"po-5\",\"point\":[{},{},{}],\"lo\":{},\"hi\":{},\"subdivisions\":{},\"cells\":{}}}",
                point.x,
                point.y,
                point.z,
                receipt.implicit.lo(),
                receipt.implicit.hi(),
                receipt.subdivisions,
                receipt.evaluated_cells
            );
        }
    });
}

#[test]
fn swept_budget_exhaustion_is_unknown_and_chart_makes_no_distance_claim() {
    with_cx(|cx| {
        let swept = SweptChart::new(
            translating_sphere(),
            SweptConfig {
                value_tolerance: 0.0,
                time_tolerance: 1.0e-12,
                max_subdivisions: 0,
            },
            cx,
        )
        .expect("finite swept support");
        let point = Point3::new(0.15, 0.3, 0.0);
        let receipt = swept.evaluate(point, cx).expect("sound coarse band");
        assert_eq!(receipt.decision, SweepDecision::Unknown);
        assert!(receipt.implicit.contains(capsule_distance(point)));
        let sample = swept.eval(point, cx);
        assert_eq!(sample.error.kind, NumericalKind::NoClaim);
        assert_eq!(swept.trace_step_claim(), TraceStepClaim::NoClaim);
        assert!(
            !swept.inside(Point3::new(-10.0, 0.0, 0.0), cx),
            "outside remains outside"
        );
    });
}

#[derive(Debug, Clone, Copy)]
struct RackCharacteristicOracle {
    root: f64,
    rank_margin: Interval,
    visibility: ProofState,
}

impl EnvelopeOracle<SphereChart> for RackCharacteristicOracle {
    fn evidence(
        &self,
        _moving: &SpacetimeChart<SphereChart>,
        _point: Point3,
        span: Interval,
        _cx: &Cx<'_>,
    ) -> Result<EnvelopeEvidence, MotionError> {
        if span.contains(self.root) {
            Ok(EnvelopeEvidence {
                field: Interval::new(-1.0e-13, 1.0e-13),
                time_derivative: Interval::new(-1.0e-13, 1.0e-13),
                rank_margin: self.rank_margin,
                characteristic_exists: ProofState::Proven,
                within_trim: ProofState::Proven,
                visible: self.visibility,
            })
        } else {
            Ok(EnvelopeEvidence {
                field: Interval::new(0.5, 1.0),
                time_derivative: Interval::WHOLE,
                rank_margin: self.rank_margin,
                characteristic_exists: ProofState::Refuted,
                within_trim: ProofState::Proven,
                visible: self.visibility,
            })
        }
    }
}

#[test]
fn envelope_endpoint_and_unproved_existence_are_distinct_classes() {
    let domain = Interval::new(0.0, 1.0);
    let endpoint = EnvelopeEvidence {
        field: Interval::new(-1.0e-14, 1.0e-14),
        // Endpoint surfaces do not need the interior characteristic equation.
        time_derivative: Interval::new(2.0, 3.0),
        rank_margin: Interval::new(0.5, 1.0),
        characteristic_exists: ProofState::Proven,
        within_trim: ProofState::Proven,
        visible: ProofState::Proven,
    };
    assert_eq!(
        classify_envelope_branch(domain, Interval::point(0.0), endpoint),
        EnvelopeBranchClass::ParameterEndpoint
    );

    let unproved = EnvelopeEvidence {
        time_derivative: Interval::new(-1.0e-14, 1.0e-14),
        characteristic_exists: ProofState::Unknown,
        ..endpoint
    };
    assert_eq!(
        classify_envelope_branch(domain, Interval::new(0.4, 0.6), unproved),
        EnvelopeBranchClass::Unknown,
        "zero-containing F and dF/dt intervals do not prove a root"
    );

    let contradictory = EnvelopeEvidence {
        field: Interval::new(1.0, 2.0),
        time_derivative: Interval::new(-1.0, 1.0),
        ..endpoint
    };
    assert_eq!(
        classify_envelope_branch(domain, Interval::new(0.4, 0.6), contradictory),
        EnvelopeBranchClass::Unknown,
        "a proven root cannot coexist with a field interval excluding zero"
    );
}

#[test]
fn g1_involute_rack_conjugate_oracle_and_regular_branch_are_derived() {
    // A no-slip rack line in the gear frame has the family
    // F(x,y,u) = x cos(u) + y sin(u) - r_b u.  Solving F = F_u = 0
    // gives p(u) = r_b(u cos(u)-sin(u), u sin(u)+cos(u)), an involute
    // rotated by pi/2.  Check both equations independently here.
    let base_radius = 0.8;
    let roll = 0.375;
    let (sin_u, cos_u) = (det::sin(roll), det::cos(roll));
    let point = Point3::new(
        base_radius * (roll * cos_u - sin_u),
        base_radius * (roll * sin_u + cos_u),
        0.0,
    );
    let field = point.x * cos_u + point.y * sin_u - base_radius * roll;
    let derivative = -point.x * sin_u + point.y * cos_u - base_radius;
    assert!(field.abs() <= 2.0e-15, "rack field residual {field}");
    assert!(
        derivative.abs() <= 2.0e-15,
        "rack characteristic residual {derivative}"
    );

    with_cx(|cx| {
        let chart = envelope(
            translating_sphere(),
            RackCharacteristicOracle {
                root: roll,
                rank_margin: Interval::new(0.5, 0.75),
                visibility: ProofState::Proven,
            },
            SweptConfig {
                value_tolerance: 5.0e-2,
                time_tolerance: 1.0e-5,
                max_subdivisions: 128,
            },
            EnvelopeConfig {
                time_tolerance: 1.0e-3,
                max_subdivisions: 64,
            },
            cx,
        )
        .expect("envelope chart");
        let trace = chart.trace(point, cx).expect("validated trace");
        assert_eq!(trace.decision, EnvelopeDecision::Enclosure);
        assert!(trace.stats.regular_branches >= 1);
        assert_eq!(trace.stats.unresolved_branches, 0);
        assert!(
            trace
                .branches
                .iter()
                .any(|branch| branch.class == EnvelopeBranchClass::RegularInterior)
        );
        println!(
            "{{\"case\":\"involute-rack\",\"subdivisions\":{},\"regular\":{},\"not_characteristic\":{},\"trimmed\":{},\"occluded\":{}}}",
            trace.stats.subdivisions,
            trace.stats.regular_branches,
            trace.stats.rejected_not_characteristic,
            trace.stats.rejected_trimmed,
            trace.stats.rejected_occluded
        );
    });
}

#[test]
fn envelope_rank_and_visibility_degeneracies_remain_unknown_or_rejected() {
    with_cx(|cx| {
        let singular = envelope(
            translating_sphere(),
            RackCharacteristicOracle {
                root: 0.375,
                rank_margin: Interval::point(0.0),
                visibility: ProofState::Proven,
            },
            SweptConfig {
                value_tolerance: 5.0e-2,
                time_tolerance: 1.0e-5,
                max_subdivisions: 128,
            },
            EnvelopeConfig {
                time_tolerance: 2.0e-3,
                max_subdivisions: 64,
            },
            cx,
        )
        .expect("singular chart");
        let receipt = singular
            .trace(Point3::new(0.0, 0.0, 0.0), cx)
            .expect("singular trace");
        assert_eq!(receipt.decision, EnvelopeDecision::Unknown);
        assert!(
            receipt
                .branches
                .iter()
                .any(|branch| branch.class == EnvelopeBranchClass::RankSingular)
        );

        let hidden = envelope(
            translating_sphere(),
            RackCharacteristicOracle {
                root: 0.375,
                rank_margin: Interval::new(0.5, 0.75),
                visibility: ProofState::Refuted,
            },
            SweptConfig {
                value_tolerance: 5.0e-2,
                time_tolerance: 1.0e-5,
                max_subdivisions: 128,
            },
            EnvelopeConfig::default(),
            cx,
        )
        .expect("hidden chart");
        let receipt = hidden
            .trace(Point3::new(0.0, 0.0, 0.0), cx)
            .expect("hidden trace");
        assert_eq!(receipt.decision, EnvelopeDecision::Enclosure);
        assert!(receipt.stats.rejected_occluded >= 1);
        assert_eq!(receipt.stats.regular_branches, 0);
        println!(
            "{{\"case\":\"degeneracy\",\"rank_unknown\":{},\"visibility_rejections\":{}}}",
            singular
                .trace(Point3::new(0.0, 0.0, 0.0), cx)
                .expect("replay")
                .stats
                .unresolved_branches,
            receipt.stats.rejected_occluded
        );
    });
}

#[test]
#[allow(clippy::too_many_lines)] // Derivation, independent pose, and tube falsifier stay together.
fn g1_wankel_apex_epitrochoid_is_derived_from_pose_and_seal_objects_stay_distinct() {
    let params = WankelParams {
        eccentricity: 0.15,
        omega: 3.0,
        crank_phase: 0.4,
        rotor_phase: -0.2,
        base_pose: Motor::translator(0.05, -0.02, 0.0),
    };
    let apex = WankelApexPoint {
        body_radius: 0.9,
        body_phase: 0.17,
    };
    let seal = WankelSealCircle {
        body_center_radius: apex.body_radius,
        body_phase: apex.body_phase,
        tip_radius: 0.02,
        clearance: 0.001,
    };
    let tube = wankel_tube(&params, Interval::new(0.0, 1.0), 10, 8).expect("wankel pose tube");
    for index in 0..=20 {
        let time = f64::from(index) / 20.0;
        let derived = apex.at(&params, time).expect("derived epitrochoid");
        let alpha = params.omega * time + params.crank_phase;
        let beta = alpha / 3.0 + params.rotor_phase;
        let pointwise = params
            .base_pose
            .compose(
                &Motor::translator(
                    params.eccentricity * det::cos(alpha),
                    params.eccentricity * det::sin(alpha),
                    0.0,
                )
                .compose(&Motor::rotor([0.0, 0.0, 1.0], beta)),
            )
            .transform_point(GaPoint {
                x: apex.body_radius * det::cos(apex.body_phase),
                y: apex.body_radius * det::sin(apex.body_phase),
                z: 0.0,
            })
            .expect("pointwise pose");
        assert!((derived.x - pointwise.x).abs() <= 2.0e-12);
        assert!((derived.y - pointwise.y).abs() <= 2.0e-12);

        // The certified tube enclosure is the proof object; pointwise samples
        // only falsify a broken derivation or sign convention.
        with_cx(|cx| {
            let enc = tube
                .point_action_over(
                    Point3::new(
                        apex.body_radius * det::cos(apex.body_phase),
                        apex.body_radius * det::sin(apex.body_phase),
                        0.0,
                    ),
                    Interval::point(time),
                    cx,
                )
                .expect("tube apex enclosure");
            for (axis, (interval, value)) in enc
                .coords
                .iter()
                .zip([derived.x, derived.y, derived.z])
                .enumerate()
            {
                assert!(
                    value >= interval.lo() - 1.0e-9 && value <= interval.hi() + 1.0e-9,
                    "sampled derived apex escaped axis {axis} tube band"
                );
            }
        });
    }

    let time = 0.31;
    let center = seal.center_at(&params, time).expect("seal center");
    let contact = seal
        .contact_at(&params, time, [1.0, 0.0])
        .expect("declared contact point");
    let center_contact_gap = det::sqrt(
        (center.x - contact.x) * (center.x - contact.x)
            + (center.y - contact.y) * (center.y - contact.y)
            + (center.z - contact.z) * (center.z - contact.z),
    );
    assert!(
        center_contact_gap > 1.0e-3,
        "finite seal center and contact are distinct"
    );
    let ideal = apex.at(&params, time).expect("ideal apex");
    assert!((ideal.x - center.x).abs() <= 2.0e-12);
    assert!((ideal.y - center.y).abs() <= 2.0e-12);
    assert!((ideal.z - center.z).abs() <= 2.0e-12);
    assert!(matches!(
        seal.contact_at(&params, time, [2.0, 0.0]),
        Err(MotionError::InvalidGeometry { .. })
    ));
}
