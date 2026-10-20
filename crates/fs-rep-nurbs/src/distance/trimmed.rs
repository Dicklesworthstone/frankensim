//! Exact trim coupling for the outward surface-distance search.

use super::{
    Meter, SearchRun, SearchStop, SurfaceDistanceBudget, SurfaceDistanceStop, SurfaceDistanceWork,
    Witness, checked, refusal, search_with_poll, tick,
};
use crate::surface::{AdmittedNurbsSurface, NurbsSurface};
use crate::trim::{
    AdmittedTrimmedPatch, Classification, TRIM_CLASSIFY_MAX_RETAINED_BYTES, TrimClassificationRun,
    TrimLoop, TrimmedPatch,
};
use crate::{NurbsError, Rat};
use fs_exec::Cx;
use fs_ivl::Interval;
use std::mem::size_of;

/// Why a trimmed query did not claim completion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrimmedSurfaceDistanceStop {
    /// A geometric refinement allowance or floating parameter resolution ended.
    Search(SurfaceDistanceStop),
    /// Exact classification refused work, storage, allocation, or arithmetic.
    /// The diagnostic is retained by `classification_refusal` on the result.
    ClassificationRefused,
}

/// Accounting for the shared geometric and exact-trim query.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TrimmedSurfaceDistanceWork {
    /// Complete original-span coverage, shared work charge, and retained payload.
    pub search: SurfaceDistanceWork,
    /// Exact classifications attempted, including unsupported binary parameters.
    pub classifications: u64,
    /// Cells excluded by an exact `Outside` verdict, separately from distance pruning.
    pub outside_cells: usize,
    /// Queries retained as `Boundary` because an exact binary Rat did not fit.
    pub unsupported_parameter_queries: u64,
}

/// Result bound to the exact immutable surface and trim generations.
///
/// `bounds` is absent only for a proved empty kept region. An unresolved query
/// without a proved kept point has upper bound +infinity and no witness.
#[derive(Debug, Clone)]
pub struct TrimmedSurfaceDistance<'s, 't> {
    source: AdmittedNurbsSurface<'s, f64>,
    trim: AdmittedTrimmedPatch<'t>,
    query: [f64; 3],
    tolerance: f64,
    bounds: Option<Interval>,
    witness: Option<Witness>,
    work: TrimmedSurfaceDistanceWork,
    classification_refusal: Option<NurbsError>,
}

impl<'s, 't> TrimmedSurfaceDistance<'s, 't> {
    /// Original immutable geometric source.
    #[must_use]
    pub const fn source(&self) -> &'s NurbsSurface<f64> {
        self.source.source()
    }
    /// Exact immutable nonzero-winding trim source.
    #[must_use]
    pub const fn trim(&self) -> &'t TrimmedPatch {
        self.trim.source()
    }
    /// Finite query in the surface's Cartesian units.
    #[must_use]
    pub const fn query(&self) -> [f64; 3] {
        self.query
    }
    /// Requested absolute distance-interval width.
    #[must_use]
    pub const fn tolerance(&self) -> f64 {
        self.tolerance
    }
    /// Outward unsigned-distance bounds, or `None` for a proved empty region.
    #[must_use]
    pub const fn bounds(&self) -> Option<Interval> {
        self.bounds
    }
    /// Actual parameter proved `Inside` the trim, when a witness was found.
    #[must_use]
    pub fn witness_param(&self) -> Option<[f64; 2]> {
        self.witness.map(|w| w.param)
    }
    /// Outward point evaluation at that exact kept parameter.
    #[must_use]
    pub fn witness_position(&self) -> Option<[Interval; 3]> {
        self.witness.map(|w| w.position)
    }
    /// Complete coverage and shared resource accounting.
    #[must_use]
    pub const fn work(&self) -> TrimmedSurfaceDistanceWork {
        self.work
    }
    /// First classifier refusal; its affected cells remain conservatively covered.
    #[must_use]
    pub fn classification_refusal(&self) -> Option<&NurbsError> {
        self.classification_refusal.as_ref()
    }
}

/// Transactional trimmed unsigned-distance outcome.
#[derive(Debug, Clone)]
#[must_use]
pub enum TrimmedSurfaceDistanceRun<'s, 't> {
    /// A finite outward interval meets the requested width and has a kept witness.
    Complete {
        /// Source- and trim-bound result.
        distance: TrimmedSurfaceDistance<'s, 't>,
    },
    /// Complete conservative coverage without a completion claim.
    Unresolved {
        /// May have +infinity upper bound and no kept witness.
        distance: TrimmedSurfaceDistance<'s, 't>,
        /// Binding limit; classification refusals retain their diagnostic.
        stop: TrimmedSurfaceDistanceStop,
    },
    /// All original-span coverage was excluded by exact outside classifications.
    Empty {
        /// Source/trim/query binding and coverage counts; no distance or witness.
        distance: TrimmedSurfaceDistance<'s, 't>,
    },
    /// No bound, witness, or partial empty-region claim is published.
    Cancelled,
}

impl<'s> AdmittedNurbsSurface<'s, f64> {
    /// Enclose unsigned distance to the image of the trim's kept parameters.
    ///
    /// The target is the closure of the actual images of in-domain parameters
    /// with nonzero winding away from the trim curves. Thus ordinary trim
    /// boundaries are approached from the kept side. At a discontinuous surface
    /// span, a parameter uses its actual owned side; a left image limit belongs
    /// to the target only when it is approached by kept parameters. Ambiguous
    /// closed cells may conservatively retain additional limits and prevent
    /// width closure, but never authorize an upper-bound witness.
    ///
    /// Exact binary parameter conversion is used where `Rat` can represent it;
    /// unsupported values stay `Boundary`. One allowance prices both geometric
    /// work and all trim classifications. The payload reserves the classifier's
    /// full 64 MiB defensive envelope plus actual trim-source capacities, in
    /// addition to the surface search storage. There is no sign, topology,
    /// portable certificate, or executor drain/finalize claim.
    ///
    /// # Errors
    /// Uses the untrimmed query's request, coverage, allocation, and finite
    /// interval-arithmetic refusals. Classifier work/storage/exactness refusals
    /// retain all affected cells and return `Unresolved` with a diagnostic.
    ///
    /// # Panics
    /// Inherits the exact trim kernel's documented arithmetic envelope: its
    /// generic `Rat` Add/Sub/Mul/Div operations can panic with a named overflow
    /// when an intermediate leaves i128. No result is published on that path.
    pub fn trimmed_distance_enclosure_with_cx<'t>(
        &self,
        trim: AdmittedTrimmedPatch<'t>,
        query: [f64; 3],
        tolerance: f64,
        budget: SurfaceDistanceBudget,
        cx: &Cx<'_>,
    ) -> Result<TrimmedSurfaceDistanceRun<'s, 't>, NurbsError> {
        trimmed_with_poll(*self, trim, query, tolerance, budget, &mut || {
            cx.checkpoint().is_err()
        })
    }
}

fn trimmed_with_poll<'s, 't>(
    source: AdmittedNurbsSurface<'s, f64>,
    trim: AdmittedTrimmedPatch<'t>,
    query: [f64; 3],
    tolerance: f64,
    budget: SurfaceDistanceBudget,
    cancel: &mut impl FnMut() -> bool,
) -> Result<TrimmedSurfaceDistanceRun<'s, 't>, NurbsError> {
    let mut domain = TrimDomain {
        source: trim,
        classifications: 0,
        outside_cells: 0,
        unsupported_parameters: 0,
        refusal: None,
    };
    let SearchRun::Complete(result) =
        search_with_poll(source, query, tolerance, budget, Some(&mut domain), cancel)?
    else {
        return Ok(TrimmedSurfaceDistanceRun::Cancelled);
    };
    let distance = TrimmedSurfaceDistance {
        source,
        trim,
        query,
        tolerance,
        bounds: result.bounds,
        witness: result.witness,
        work: TrimmedSurfaceDistanceWork {
            search: result.work,
            classifications: domain.classifications,
            outside_cells: domain.outside_cells,
            unsupported_parameter_queries: domain.unsupported_parameters,
        },
        classification_refusal: domain.refusal,
    };
    Ok(match result.stop {
        Some(SearchStop::Refinement(stop)) => TrimmedSurfaceDistanceRun::Unresolved {
            distance,
            stop: TrimmedSurfaceDistanceStop::Search(stop),
        },
        Some(SearchStop::TrimRefused) => TrimmedSurfaceDistanceRun::Unresolved {
            distance,
            stop: TrimmedSurfaceDistanceStop::ClassificationRefused,
        },
        None if distance.bounds.is_none() => TrimmedSurfaceDistanceRun::Empty { distance },
        None => TrimmedSurfaceDistanceRun::Complete { distance },
    })
}

pub(super) struct TrimDomain<'a> {
    source: AdmittedTrimmedPatch<'a>,
    classifications: u64,
    pub(super) outside_cells: usize,
    unsupported_parameters: u64,
    pub(super) refusal: Option<NurbsError>,
}

impl TrimDomain<'_> {
    pub(super) fn retained_payload(
        &self,
        meter: &mut Meter,
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<Option<usize>, NurbsError> {
        let raw = self.source.source();
        if !meter.charge(raw.loops.len() as u128 * 8 + 1) {
            return Err(refusal(
                "trimmed-distance work cannot admit trim source storage",
            ));
        }
        let reserve = usize::try_from(TRIM_CLASSIFY_MAX_RETAINED_BYTES)
            .map_err(|_| refusal("trimmed-distance classifier envelope exceeds usize"))?;
        let mut bytes = checked(
            raw.loops
                .capacity()
                .checked_mul(size_of::<TrimLoop>())
                .and_then(|n| n.checked_add(size_of::<TrimmedPatch>()))
                .and_then(|n| n.checked_add(reserve)),
        )?;
        let mut steps = 0;
        for trim_loop in &raw.loops {
            if tick(&mut steps, cancel) {
                return Ok(None);
            }
            let curve = &trim_loop.curve;
            bytes = checked(
                curve
                    .cpw
                    .capacity()
                    .checked_mul(size_of::<[Rat; 4]>())
                    .and_then(|n| {
                        n.checked_add(curve.knots.knots.capacity().checked_mul(size_of::<Rat>())?)
                    })
                    .and_then(|n| bytes.checked_add(n)),
            )?;
        }
        if cancel() { Ok(None) } else { Ok(Some(bytes)) }
    }

    pub(super) fn classify(
        &mut self,
        min: [f64; 2],
        max: [f64; 2],
        meter: &mut Meter,
        cancel: &mut impl FnMut() -> bool,
    ) -> Result<Option<Classification>, NurbsError> {
        if cancel() {
            return Ok(None);
        }
        if self.refusal.is_some() {
            return Ok(Some(Classification::Boundary));
        }
        if !meter.charge(256) {
            self.refusal = Some(refusal(
                "trimmed-distance work cannot admit exact parameter conversion",
            ));
            return Ok(Some(Classification::Boundary));
        }
        self.classifications += 1;
        let (Some(u0), Some(v0), Some(u1), Some(v1)) = (
            binary_rat(min[0]),
            binary_rat(min[1]),
            binary_rat(max[0]),
            binary_rat(max[1]),
        ) else {
            self.unsupported_parameters += 1;
            return Ok(Some(Classification::Boundary));
        };
        let before = meter.limit - meter.used;
        let mut remaining = before;
        let result = self.source.classify_box_with_depth_budget_and_poll(
            [u0, v0],
            [u1, v1],
            if min == max { 4 } else { 2 },
            &mut remaining,
            cancel,
        );
        meter.used += before - remaining;
        match result {
            Ok(TrimClassificationRun::Complete { classification }) => Ok(Some(classification)),
            Ok(TrimClassificationRun::Cancelled) => Ok(None),
            Err(error @ (NurbsError::Domain { .. } | NurbsError::Exactness { .. })) => {
                self.refusal = Some(error);
                Ok(Some(Classification::Boundary))
            }
            Err(error) => Err(error),
        }
    }
}

/// Exact conversion of the stored binary value; no grid snapping, approximation,
/// saturating cast, or Rat arithmetic is used before representability is known.
fn binary_rat(value: f64) -> Option<Rat> {
    if !value.is_finite() {
        return None;
    }
    let bits = value.to_bits();
    let negative = bits >> 63 != 0;
    let exponent_bits = i32::try_from((bits >> 52) & 0x7ff).ok()?;
    let mut magnitude = bits & ((1u64 << 52) - 1);
    let mut exponent = if exponent_bits == 0 {
        -1074
    } else {
        magnitude |= 1u64 << 52;
        exponent_bits - 1023 - 52
    };
    if magnitude == 0 {
        return Some(Rat::int(0));
    }
    let trailing = magnitude.trailing_zeros();
    magnitude >>= trailing;
    exponent += i32::try_from(trailing).ok()?;
    let (magnitude, denominator) = if exponent >= 0 {
        let shift = u32::try_from(exponent).ok()?;
        let magnitude = u128::from(magnitude);
        if shift >= 128 || magnitude > (u128::MAX >> shift) {
            return None;
        }
        (magnitude << shift, 1)
    } else {
        let shift = u32::try_from(-exponent).ok()?;
        if shift >= 127 {
            return None;
        }
        (u128::from(magnitude), 1i128 << shift)
    };
    let numerator = if negative && magnitude == 1u128 << 127 {
        i128::MIN
    } else {
        let positive = i128::try_from(magnitude).ok()?;
        if negative { -positive } else { positive }
    };
    Some(Rat::new(numerator, denominator))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{KnotVector, NurbsCurve};

    #[test]
    fn exact_binary_conversion_never_rounds_a_trim_witness() {
        assert_eq!(binary_rat(-0.0), Some(Rat::int(0)));
        assert_eq!(
            binary_rat(0.1),
            Some(Rat::new(3_602_879_701_896_397, 36_028_797_018_963_968))
        );
        assert_eq!(
            binary_rat(2.0f64.powi(-126)),
            Some(Rat::new(1, 1i128 << 126))
        );
        assert_eq!(binary_rat(-2.0f64.powi(127)), Some(Rat::new(i128::MIN, 1)));
        for value in [
            2.0f64.powi(-127),
            2.0f64.powi(127),
            f64::from_bits(1),
            f64::NAN,
        ] {
            assert!(binary_rat(value).is_none());
        }
    }

    #[test]
    fn g4_cancellation_inside_exact_trim_and_final_publication_exposes_no_result() {
        let k = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
        let surface = NurbsSurface::new(
            k.clone(),
            k,
            &[
                vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
            ],
            &[vec![1.0; 2], vec![1.0; 2]],
        )
        .unwrap();
        let k =
            KnotVector::new([0, 0, 1, 2, 3, 4, 4].into_iter().map(Rat::int).collect(), 1).unwrap();
        let z = Rat::int(0);
        let h = Rat::new(1, 2);
        let o = Rat::int(1);
        let curve = NurbsCurve::new(k, &[[h, z], [o, z], [o, o], [h, o], [h, z]], &[o; 5]).unwrap();
        let trim = TrimmedPatch::new(vec![TrimLoop::new(curve).unwrap()]);
        let budget = SurfaceDistanceBudget {
            max_splits: 8,
            max_work_units: 8_000_000,
            max_retained_bytes: 80 * 1024 * 1024,
        };
        let mut polls = 0;
        let run = trimmed_with_poll(
            surface.admit().unwrap(),
            trim.admit().unwrap(),
            [0.8, 0.5, 1.0],
            1e-8,
            budget,
            &mut || {
                polls += 1;
                false
            },
        )
        .unwrap();
        assert!(!matches!(run, TrimmedSurfaceDistanceRun::Cancelled));
        assert!(polls > 30);
        for boundary in [1, polls / 3, polls - 1, polls] {
            let mut seen = 0;
            let run = trimmed_with_poll(
                surface.admit().unwrap(),
                trim.admit().unwrap(),
                [0.8, 0.5, 1.0],
                1e-8,
                budget,
                &mut || {
                    seen += 1;
                    seen == boundary
                },
            )
            .unwrap();
            assert!(matches!(run, TrimmedSurfaceDistanceRun::Cancelled));
            assert_eq!(seen, boundary, "no checkpoint after observed cancellation");
        }
    }
}
