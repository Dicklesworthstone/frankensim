//! Gauntlet G3 relations for the production geometry conversion surface.
//!
//! These checks supplement the fixed conversion, cancellation, and refusal
//! pins in `conformance.rs`; they do not replace those cases or claim a
//! convergence rate for the dense sampled representation.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Convert, ErrBudget, Point3, SampledSdf};
use fs_propcheck::Shrink;
use fs_propcheck::metamorphic::{
    RelationCase, RelationObservation, Tolerance, check_relation, refinement_monotonicity,
};

const CENTERS: [Point3; 3] = [
    Point3::new(0.0, 0.0, 0.0),
    Point3::new(0.125, -0.25, 0.375),
    Point3::new(-0.5, 0.25, -0.125),
];
const RADII: [f64; 3] = [0.25, 0.375, 0.5];
const BASE_BUDGETS: [f64; 3] = [0.5, 0.75, 1.0];

#[derive(Debug, Clone, PartialEq)]
struct SphereCase {
    center_index: u64,
    radius_index: u64,
    requested_error: f64,
}

impl Shrink for SphereCase {
    fn shrink_candidates(&self) -> Vec<Self> {
        let mut candidates = Vec::new();
        if self.center_index != 0 {
            candidates.push(Self {
                center_index: 0,
                ..self.clone()
            });
        }
        if self.radius_index != 0 {
            candidates.push(Self {
                radius_index: 0,
                ..self.clone()
            });
        }
        if self.requested_error.to_bits() != BASE_BUDGETS[0].to_bits() {
            candidates.push(Self {
                requested_error: BASE_BUDGETS[0],
                ..self.clone()
            });
        }
        candidates
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RefinementPower(u32);

impl RefinementPower {
    fn divisor(self) -> f64 {
        match self.0 {
            1 => 2.0,
            2 => 4.0,
            _ => unreachable!("generated refinement power is exactly one or two"),
        }
    }
}

impl Shrink for RefinementPower {
    fn shrink_candidates(&self) -> Vec<Self> {
        match self.0 {
            0 | 1 => Vec::new(),
            _ => vec![Self(1)],
        }
    }
}

#[derive(Debug)]
struct ConversionReceipt {
    center: Point3,
    radius: f64,
    requested_error: f64,
    certified_qoi: f64,
    resolution: u32,
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x6E0_0201,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

#[test]
fn g3_sphere_conversion_receipts_improve_under_power_of_two_refinement() {
    with_cx(|cx| {
        let operator = |case: &SphereCase| {
            let center = CENTERS[case.center_index as usize];
            let radius = RADII[case.radius_index as usize];
            let requested_error = case.requested_error;
            let converted: fs_evidence::Certified<SampledSdf> = SphereChart { center, radius }
                .convert(
                    ErrBudget {
                        abs_sd_error: requested_error,
                    },
                    cx,
                )
                .expect("catalogued sphere conversion is admitted and certified");

            ConversionReceipt {
                center,
                radius,
                requested_error,
                certified_qoi: converted.qoi,
                resolution: converted.value.resolution(),
            }
        };
        let relation = refinement_monotonicity(
            "sphere-sampled-sdf-power-of-two-budget",
            Tolerance::NonIncreasing { max_increase: 0.0 },
            |case: &SphereCase, refinement: &RefinementPower| SphereCase {
                center_index: case.center_index,
                radius_index: case.radius_index,
                requested_error: case.requested_error / refinement.divisor(),
            },
            |base: &ConversionReceipt,
             transformed: &ConversionReceipt,
             refinement: &RefinementPower,
             tolerance: Tolerance| {
                let expected_budget = base.requested_error / refinement.divisor();
                let qoi = tolerance.evaluate_scalar(base.certified_qoi, transformed.certified_qoi);
                let resolution_margin = if transformed.resolution >= base.resolution {
                    0.0
                } else {
                    -1.0
                };
                let refinement_margin = if refinement.0 > 0
                    && transformed.requested_error.to_bits() == expected_budget.to_bits()
                    && transformed.requested_error < base.requested_error
                {
                    0.0
                } else {
                    -1.0
                };
                let same_geometry_margin = if transformed.center.x.to_bits()
                    == base.center.x.to_bits()
                    && transformed.center.y.to_bits() == base.center.y.to_bits()
                    && transformed.center.z.to_bits() == base.center.z.to_bits()
                    && transformed.radius.to_bits() == base.radius.to_bits()
                {
                    0.0
                } else {
                    -1.0
                };

                RelationObservation::new(
                    qoi.margin()
                        .min(resolution_margin)
                        .min(refinement_margin)
                        .min(same_geometry_margin),
                    "tighter certified conversion has nonincreasing QoI and nondecreasing resolution",
                )
            },
        );

        check_relation(
            "fs-geom::SphereChart::convert<SampledSdf>",
            0x6E0_0202,
            8,
            |stream| {
                let center_index = stream.next_u64() % CENTERS.len() as u64;
                let radius_index = stream.next_u64() % RADII.len() as u64;
                let budget_index = stream.next_u64() % BASE_BUDGETS.len() as u64;
                let refinement = RefinementPower(1 + (stream.next_u64() % 2) as u32);
                RelationCase::new(
                    SphereCase {
                        center_index,
                        radius_index,
                        requested_error: BASE_BUDGETS[budget_index as usize],
                    },
                    refinement,
                )
            },
            &operator,
            &relation,
        );
    });
}
