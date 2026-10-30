//! Scratch admission for the original shared-denominator QR owner.

use super::FitOptions;

/// Conservative scratch for one shared-pole relocation, excluding caller data,
/// fitted models, pole eigensolves and later passivity work. The options must
/// cover the direct/improper columns used by every response. Every response
/// uses the original channel-eliminated common-denominator owner.
/// Returns `None` when the byte bound cannot be represented by this target.
#[must_use]
pub fn shared_pole_relocation_scratch_bytes(
    samples: usize,
    responses: usize,
    options: &FitOptions,
) -> Option<usize> {
    let samples = samples as u128;
    let responses_count = responses as u128;
    let order = options.order as u128;
    let local = order + u128::from(options.fit_d) + u128::from(options.fit_e);
    let rows = samples
        .saturating_mul(responses_count)
        .saturating_mul(2)
        .saturating_add(1);
    let columns = order + 1;
    // Up to four numerator QR shapes and one transformed sigma column.
    let projection = samples.saturating_mul(local + 1).saturating_mul(16);
    // Original/copy/equilibrated QR/fixed-sigma fallback storage can coexist.
    let values = rows
        .saturating_mul(columns)
        .saturating_mul(4)
        .saturating_add(columns.saturating_mul(columns).saturating_mul(8))
        .saturating_add(projection);
    usize::try_from(values.saturating_mul(8)).ok()
}

#[cfg(test)]
mod tests;
