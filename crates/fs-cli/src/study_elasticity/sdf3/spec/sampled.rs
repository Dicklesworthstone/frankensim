//! Native sampled-grid input, retained verbatim in canonical study identity.
//! Values describe the trilinear interpolant's negative set. They do not carry
//! an error bound to the source surface, even when exported as distance samples.
use super::*;
use fs_cutfem::sampled3::{MAX_SAMPLED3_VALUES, SampledSdf3};

#[derive(Debug, Clone)]
pub(crate) struct SampledInput {
    pub field: SampledSdf3,
    /// Additional bounded-input allowance for samples, knots and serialized
    /// forms. Like the surrounding admission model, this is NOT measured RSS.
    pub memory_envelope: usize,
}

pub(super) fn parse(values: &[&Node], bounds: ([f64; 3], [f64; 3])) -> Result<SampledInput> {
    if values.len() != 4 {
        return Err(invalid("sampled-grid-sdf requires origin-m, spacing-m, shape and values-m"));
    }
    let origin = vector(values[0])?;
    let spacing = vector(values[1])?;
    let dimensions = super::super::super::list(values[2], "sample grid shape")?;
    if dimensions.len() != 3 { return Err(invalid("sample grid shape requires three integer knot counts")); }
    let shape = [count(&dimensions[0])?, count(&dimensions[1])?, count(&dimensions[2])?];
    let expected = shape.iter().try_fold(1usize, |n, &d| {
        if d < 2 { None } else { n.checked_mul(d) }
    }).ok_or_else(|| invalid("sample grid needs at least two knots per axis and a representable count"))?;
    let nodes = super::super::super::list(values[3], "x-fastest sample values")?;
    if expected > MAX_SAMPLED3_VALUES || nodes.len() != expected {
        return Err(invalid("sample grid values must match nx*ny*nz exactly and not exceed 1048576"));
    }
    let samples = nodes.iter().map(scalar).collect::<Result<Vec<_>>>()?;
    let memory_envelope = samples.len().saturating_mul(256);
    // This is bounded input admission. Physics/refinement/replay subsequently
    // use the ordinary quadrature and solver cancellation/budget owners.
    let field = SampledSdf3::with_control(origin, spacing, shape, samples,
        &mut || ControlFlow::Continue(())).map_err(|e| invalid(e.to_string()))?;
    let (lo, hi) = field.bounds();
    for a in 0..3 {
        let (left, right) = (bounds.0[a], bounds.1[a]);
        let terminal = left + (right - left);
        if !left.is_finite() || !right.is_finite() || !terminal.is_finite()
            || left >= right || left < lo[a] || right.max(terminal) > hi[a] {
            return Err(invalid("sample grid must cover the entire ordered physical box, including rounded octree endpoints; pad the grid rather than extrapolating"));
        }
    }
    Ok(SampledInput { field, memory_envelope })
}

#[cfg(test)]
#[path = "sampled_tests.rs"]
mod tests;
