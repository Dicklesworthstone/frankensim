//! Deterministic bounded integration of the original source-history kernel.
//!
//! Adaptive Simpson differences estimate quadrature error; they are not an
//! interval certificate or a bound on the producer's mechanical interpolation.
//! Every accepted leaf also resolves geometry using the admitted continuous
//! source-speed bound. This prevents missed narrow close passes between nodes.

use super::super::{acoustic, check_cancel, sizing};
use crate::render::RenderError;
use fs_aeroac::AeroacError;
use fs_math::dd::Dd;

const MAX_DEPTH: usize = 24;

/// Accuracy and work controls for one retained moving-force history renderer.
#[derive(Clone, Copy, Debug)]
pub struct ForceHistoryQuadrature {
    /// Absolute adaptive error target for each path's complete history integral
    /// J [N s/m³], in the largest Cartesian component. Original intervals share
    /// this target in proportion to their duration. It is an error ESTIMATE.
    pub absolute_tolerance_n_s_m3: f64,
    /// Maximum binary subdivision depth per merged force/geometry cell, 0..=24.
    pub maximum_subdivision_depth: usize,
    /// Integrand evaluations across preparation AND every callback attempt.
    /// Failed/cancelled attempts do not refund work already performed.
    pub maximum_evaluations: u64,
}

pub(super) struct Work {
    pub(super) used: u64,
    limit: u64,
}

impl Work {
    pub(super) fn new(config: ForceHistoryQuadrature) -> Result<Self, RenderError> {
        if !config.absolute_tolerance_n_s_m3.is_finite()
            || config.absolute_tolerance_n_s_m3 <= 0.0
            || config.maximum_subdivision_depth > MAX_DEPTH
            || config.maximum_evaluations == 0
        {
            return Err(sizing(
                "moving force integration requires finite accuracy and work limits",
            ));
        }
        Ok(Self {
            used: 0,
            limit: config.maximum_evaluations,
        })
    }

    fn sample(
        &mut self,
        time: f64,
        evaluate: &mut impl FnMut(f64) -> Result<([f64; 3], f64), RenderError>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<([Dd; 3], f64), RenderError> {
        check_cancel(cancelled)?;
        if self.used >= self.limit {
            return Err(sizing(
                "moving force history quadrature work budget exhausted",
            ));
        }
        self.used += 1;
        let (value, distance) = evaluate(time)?;
        if value.iter().any(|v| !v.is_finite()) || !distance.is_finite() || distance <= 0.0 {
            return Err(acoustic(AeroacError::NonFinite {
                what: "moving force history integrand",
            }));
        }
        Ok((value.map(Dd::from_f64), distance))
    }
}

#[derive(Clone, Copy)]
struct Node {
    first: f64,
    last: f64,
    first_value: [Dd; 3],
    middle_value: [Dd; 3],
    middle_distance_m: f64,
    last_value: [Dd; 3],
    estimate: [Dd; 3],
    tolerance: f64,
    depth: usize,
}

fn simpson(first: f64, last: f64, a: [Dd; 3], m: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    let width = (Dd::from_f64(last) - Dd::from_f64(first)) / Dd::from_f64(6.0);
    core::array::from_fn(|i| width * (a[i] + Dd::from_f64(4.0) * m[i] + b[i]))
}

pub(super) fn integrate(
    first: f64,
    last: f64,
    tolerance: f64,
    maximum_source_speed_m_s: f64,
    minimum_distance_m: f64,
    config: ForceHistoryQuadrature,
    work: &mut Work,
    mut evaluate: impl FnMut(f64) -> Result<([f64; 3], f64), RenderError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<[Dd; 3], RenderError> {
    check_cancel(cancelled)?;
    if first == last {
        return Ok([Dd::ZERO; 3]);
    }
    if !first.is_finite()
        || !last.is_finite()
        || last < first
        || !tolerance.is_finite()
        || tolerance < 0.0
        || !maximum_source_speed_m_s.is_finite()
        || maximum_source_speed_m_s < 0.0
        || !minimum_distance_m.is_finite()
        || minimum_distance_m < 0.0
    {
        return Err(sizing(
            "moving force quadrature interval is unrepresentable",
        ));
    }
    let middle = first.midpoint(last);
    if !(first < middle && middle < last) {
        return Err(sizing(
            "moving force history has no representable interior time",
        ));
    }
    let (a, _) = work.sample(first, &mut evaluate, cancelled)?;
    let (m, middle_distance_m) = work.sample(middle, &mut evaluate, cancelled)?;
    let (b, _) = work.sample(last, &mut evaluate, cancelled)?;
    let initial = Node {
        first,
        last,
        first_value: a,
        middle_value: m,
        middle_distance_m,
        last_value: b,
        estimate: simpson(first, last, a, m, b),
        tolerance,
        depth: 0,
    };
    let mut stack = [initial; MAX_DEPTH + 1];
    let mut pending = 1;
    let mut total = [Dd::ZERO; 3];
    while pending > 0 {
        check_cancel(cancelled)?;
        pending -= 1;
        let node = stack[pending];
        let middle = node.first.midpoint(node.last);
        let left_middle = node.first.midpoint(middle);
        let right_middle = middle.midpoint(node.last);
        if !(node.first < left_middle
            && left_middle < middle
            && middle < right_middle
            && right_middle < node.last)
        {
            return Err(sizing(
                "moving force quadrature subdivision is unrepresentable",
            ));
        }
        let (lm, left_distance) = work.sample(left_middle, &mut evaluate, cancelled)?;
        let (rm, right_distance) = work.sample(right_middle, &mut evaluate, cancelled)?;
        let left = simpson(node.first, middle, node.first_value, lm, node.middle_value);
        let right = simpson(middle, node.last, node.middle_value, rm, node.last_value);
        let difference: [Dd; 3] = core::array::from_fn(|i| left[i] + right[i] - node.estimate[i]);
        let error = difference
            .iter()
            .map(|v| v.to_f64().abs() / 15.0)
            .fold(0.0_f64, f64::max);
        if !error.is_finite() || difference.iter().any(|v| !v.to_f64().is_finite()) {
            return Err(acoustic(AeroacError::NonFinite {
                what: "moving force quadrature estimate",
            }));
        }
        // Continuous speed bounds control displacement throughout the cell.
        // Every range stays at least half the sampled midpoint range, even
        // allowing a rounded midpoint, so a narrow peak cannot hide between nodes.
        let resolved_geometry = maximum_source_speed_m_s == 0.0
            || (node.last - node.first).next_up()
                <= ((0.5 * node.middle_distance_m).next_down() / maximum_source_speed_m_s)
                    .next_down();
        // A moving microphone's current point also needs separation from ALL
        // earlier source positions, not only from its causal emission event.
        // Use the full cell width (rather than a rounded half width) with the
        // original outward speed bound. Unresolved past clearance subdivides
        // or refuses; zero force never bypasses this geometric obligation.
        let clear_history = minimum_distance_m == 0.0
            || maximum_source_speed_m_s == 0.0
            || (node.middle_distance_m.next_down()
                - ((node.last - node.first).next_up() * maximum_source_speed_m_s).next_up())
            .next_down()
                >= minimum_distance_m;
        if resolved_geometry && clear_history && error <= node.tolerance {
            for i in 0..3 {
                total[i] = total[i] + left[i] + right[i] + difference[i] / Dd::from_f64(15.0);
            }
            continue;
        }
        if node.depth >= config.maximum_subdivision_depth {
            return Err(acoustic(AeroacError::NotConverged {
                what: "moving force history geometry or quadrature depth",
                residual: error,
            }));
        }
        let tolerance = 0.5 * node.tolerance;
        stack[pending] = Node {
            first: middle,
            last: node.last,
            first_value: node.middle_value,
            middle_value: rm,
            middle_distance_m: right_distance,
            last_value: node.last_value,
            estimate: right,
            tolerance,
            depth: node.depth + 1,
        };
        stack[pending + 1] = Node {
            first: node.first,
            last: middle,
            first_value: node.first_value,
            middle_value: lm,
            middle_distance_m: left_distance,
            last_value: node.middle_value,
            estimate: left,
            tolerance,
            depth: node.depth + 1,
        };
        pending += 2;
    }
    if total.iter().any(|v| !v.to_f64().is_finite()) {
        return Err(acoustic(AeroacError::NonFinite {
            what: "integrated moving force history",
        }));
    }
    Ok(total)
}
