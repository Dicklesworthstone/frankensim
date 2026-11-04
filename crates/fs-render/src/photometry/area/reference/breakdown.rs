//! Attribution of the native direct-area solve, without another integration.

use super::*;
use crate::photometry::fixture::{AssembledFixture, FixtureId};
use crate::photometry::{EmitterIlluminance, FixtureIlluminance, IlluminanceEstimate};
use core::ops::Range;

/// Direct finite-area result and its source/fixture contributions.
///
/// `total` is the unchanged native aggregate result, including its model scope,
/// numerical identity and, for occluded queries, shared budget and termination.
/// Components enclose that same restricted direct-area model. They do not add
/// point light, reflection or in-scattering, and have no independent tolerance
/// or total physical-error claim. Separate interval reductions and midpoints
/// need not reproduce the aggregate's last bit or interval endpoints.
#[derive(Debug, Clone, PartialEq)]
pub struct AreaIlluminanceBreakdown<T> {
    /// Original aggregate result from the same single evaluation.
    pub total: T,
    /// One contribution per original finite source, including dark emitters.
    /// `source` indexes the reference's `source_bindings()` slice; supplemental
    /// and flat sources retain `identity: None`.
    pub emitters: Vec<EmitterIlluminance>,
    /// One finite-component subtotal per rig fixture in caller order, including
    /// empty and point-only fixtures as zero. Flat inputs have no fixture rows;
    /// supplemental sources are not assigned to an invented fixture.
    pub fixtures: Vec<FixtureIlluminance>,
}

pub(super) struct FixtureSources {
    fixture: FixtureId,
    sources: Range<usize>,
}

/// Copy only the admitted rig's attribution ranges, retaining no source/rig
/// borrow. Their count is bounded by the original fixture-assembly admission;
/// the stricter breakdown output limit is checked when a breakdown is requested.
pub(super) fn fixture_sources(
    fixtures: &[AssembledFixture],
    cx: &Cx<'_>,
) -> Result<Vec<FixtureSources>, AreaReferenceError> {
    let mut result = reserved(fixtures.len())?;
    for fixture in fixtures {
        cx.checkpoint()?;
        result.push(FixtureSources {
            fixture: fixture.id,
            sources: fixture.surfaces.clone(),
        });
    }
    Ok(result)
}

#[derive(Clone, Copy)]
struct SourceRange {
    // Raw endpoints enclose the sum of CURRENT leaf endpoints. Intersecting
    // successive bounds is a separate operation and cannot replace this state.
    raw: Interval,
    enclosed: Interval,
}

/// Private per-query staging. Ordinary aggregate queries never allocate it.
pub(super) struct Attribution {
    sources: Vec<SourceRange>,
    emitters: Vec<EmitterIlluminance>,
    fixtures: Vec<FixtureIlluminance>,
}

impl Attribution {
    pub(super) fn new(
        reference: &UnoccludedAreaReference,
        cx: &Cx<'_>,
    ) -> Result<Self, AreaReferenceError> {
        cx.checkpoint()?;
        let source_count = reference.bindings.len();
        let fixture_count = reference.fixtures.len();
        let outputs = source_count
            .checked_add(fixture_count)
            .and_then(|count| count.checked_add(1))
            .ok_or(AreaReferenceError::Budget)?;
        if outputs > reference.max_outputs {
            return Err(AreaReferenceError::Budget);
        }
        let mut sources = reserved(source_count)?;
        for _ in 0..source_count {
            cx.checkpoint()?;
            sources.push(SourceRange {
                raw: zero(),
                enclosed: zero(),
            });
        }
        let emitters = reserved(source_count)?;
        let fixtures = reserved(fixture_count)?;
        cx.checkpoint()?;
        Ok(Self {
            sources,
            emitters,
            fixtures,
        })
    }

    pub(super) fn add_root(
        &mut self,
        source: usize,
        range: Interval,
    ) -> Result<(), AreaReferenceError> {
        let entry = &mut self.sources[source];
        entry.raw = nonnegative(add(entry.raw, range));
        if !finite(entry.raw) {
            return Err(AreaReferenceError::NumericalRange);
        }
        entry.enclosed = entry.raw;
        Ok(())
    }

    pub(super) fn replace(
        &mut self,
        source: usize,
        parent: Interval,
        a: Interval,
        b: Interval,
    ) -> Result<(), AreaReferenceError> {
        let entry = &mut self.sources[source];
        entry.raw = occlusion::replace_parent(entry.raw, parent, a, b)?;
        entry.enclosed = entry
            .enclosed
            .intersect(entry.raw)
            .ok_or(AreaReferenceError::NumericalRange)?;
        Ok(())
    }

    pub(super) fn finish<T>(
        mut self,
        reference: &UnoccludedAreaReference,
        total: T,
        cx: &Cx<'_>,
    ) -> Result<AreaIlluminanceBreakdown<T>, AreaReferenceError> {
        for (source, (range, binding)) in self.sources.iter().zip(&reference.bindings).enumerate() {
            cx.checkpoint()?;
            self.emitters.push(EmitterIlluminance {
                source,
                identity: binding.identity,
                estimate: estimate(range.enclosed)?,
            });
        }
        for fixture in &reference.fixtures {
            cx.checkpoint()?;
            let mut range = zero();
            for source in self
                .sources
                .get(fixture.sources.clone())
                .ok_or(AreaReferenceError::NumericalRange)?
            {
                cx.checkpoint()?;
                range = nonnegative(add(range, source.enclosed));
            }
            self.fixtures.push(FixtureIlluminance {
                fixture: fixture.fixture,
                estimate: estimate(range)?,
            });
        }
        cx.checkpoint()?;
        Ok(AreaIlluminanceBreakdown {
            total,
            emitters: self.emitters,
            fixtures: self.fixtures,
        })
    }
}

fn estimate(range: Interval) -> Result<IlluminanceEstimate, AreaReferenceError> {
    if !finite(range) {
        return Err(AreaReferenceError::NumericalRange);
    }
    let lux = range.lo() + (range.hi() - range.lo()) * 0.5;
    if !lux.is_finite() {
        return Err(AreaReferenceError::NumericalRange);
    }
    Ok(IlluminanceEstimate {
        lux,
        numerical: NumericalCertificate::enclosure(range.lo(), range.hi()),
        total_error_bound_lux: None,
    })
}

fn reserved<T>(count: usize) -> Result<Vec<T>, AreaReferenceError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| AreaReferenceError::Allocation)?;
    Ok(values)
}
