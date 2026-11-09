//! Authored reservoir schedules through the original radiative endpoint dual.
//!
//! Only signal values/coefficients and each fixed emissivity are coordinates.
//! Knots, interpolation, time grids, event sides and material selection stay fixed.

use super::{ScanObjectiveSampling, ScanSensitivityError, ScanningEnthalpyRun,
    add, finite, invalid, poll};
use super::super::ScanRadiationGradient;
use crate::bc::history::march::moving::sensitivity::schedule::{count, pull_signal};
use fs_exec::Cx;

/// Physical radiative controls in original patch declaration order.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanRadiationScheduleGradient {
    /// Original boundary region; no sorting or merging of independently authored patches.
    pub region: String,
    /// Per-kelvin derivatives in the original reservoir signal representation:
    /// one constant, ramp [from,to], table values or stored Chebyshev coefficients.
    /// The coefficient-zero convention is inherited from the original signal.
    pub reservoir_values_k: Vec<f64>,
    /// Derivative per unit absolute hemispherical emissivity shared by the
    /// complete trajectory. Source identity/support are fixed; at epsilon=1
    /// only feasible perturbations apply.
    pub emissivity: f64,
}

// This counts the actual returned scalar controls, including one emissivity
// per patch. Value-only evaluation has no radiation derivative requirement.
pub(super) fn control_count(cx: &Cx<'_>, run: &ScanningEnthalpyRun<'_, '_, '_>)
    -> Result<usize, ScanSensitivityError>
{
    let mut total = 0usize;
    for (i, source) in run.sources().iter().enumerate() {
        poll(cx, i)?;
        total = total.checked_add(source.schedule_control_count()?)
            .ok_or_else(|| invalid("scan authored-control count overflow"))?;
    }
    if let Some(radiation) = &run.radiation {
        for (i, patch) in radiation.patches.iter().enumerate() {
            poll(cx, i)?;
            total = total.checked_add(count(&patch.ambient_temperature)?)
                .and_then(|n| n.checked_add(1))
                .ok_or_else(|| invalid("scan radiative-control count overflow"))?;
        }
    }
    Ok(total)
}

pub(super) fn pullback(cx: &Cx<'_>, run: &ScanningEnthalpyRun<'_, '_, '_>,
    time_s: f64, sampling: ScanObjectiveSampling, gradient: Option<ScanRadiationGradient>,
) -> Result<Vec<ScanRadiationScheduleGradient>, ScanSensitivityError> {
    poll(cx, 0)?;
    let (radiation, gradient) = match (&run.radiation, gradient) {
        (None, None) => return Ok(Vec::new()),
        (Some(radiation), Some(gradient)) => (radiation, gradient),
        _ => return Err(invalid("scan radiative-gradient owner changed")),
    };
    let patches = &radiation.patches;
    if gradient.radiation_regions.len() != patches.len()
        || gradient.reservoir_temperatures.len() != patches.len()
        || gradient.emissivities.len() != patches.len()
    { return Err(invalid("scan radiative-gradient layout changed")); }
    let left = sampling == ScanObjectiveSampling::EventLeft;
    let mut output = Vec::new();
    output.try_reserve_exact(patches.len())
        .map_err(|_| invalid("scan radiative-gradient allocation refused"))?;
    for (i, patch) in patches.iter().enumerate() {
        poll(cx, i)?;
        if patch.region != gradient.radiation_regions[i] {
            return Err(invalid("scan radiative-patch identity changed"));
        }
        // The immutable run and accepted endpoint already validated the original
        // physical signal on this exact side. Reuse its value basis unchanged.
        let seed = finite(gradient.reservoir_temperatures[i])?;
        output.push(ScanRadiationScheduleGradient { region: patch.region.clone(),
            reservoir_values_k: pull_signal(cx, &patch.ambient_temperature, time_s, left, seed)?,
            emissivity: finite(gradient.emissivities[i])? });
    }
    poll(cx, output.len())?;
    Ok(output)
}

pub(super) fn accumulate(cx: &Cx<'_>, total: &mut Vec<ScanRadiationScheduleGradient>,
    values: Vec<ScanRadiationScheduleGradient>,
) -> Result<(), ScanSensitivityError> {
    if total.is_empty() {
        *total = values;
        return Ok(());
    }
    if total.len() != values.len() {
        return Err(invalid("scan radiative identity changed during reverse sweep"));
    }
    for (i, (a, b)) in total.iter_mut().zip(values).enumerate() {
        poll(cx, i)?;
        if a.region != b.region {
            return Err(invalid("scan radiative region order changed during reverse sweep"));
        }
        add(cx, &mut a.reservoir_values_k, &b.reservoir_values_k)?;
        a.emissivity = finite(a.emissivity + b.emissivity)?;
    }
    Ok(())
}
