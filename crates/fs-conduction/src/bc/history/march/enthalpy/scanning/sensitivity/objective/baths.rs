//! Authored prescribed-temperature histories through the original constraint dual.
//!
//! The boundary builder retains each constrained vertex's owning region. Its
//! original precedence is used here, including static regions that shadow a
//! scheduled region at a shared vertex. No heat-load seed or extra transpose is
//! substituted for the prescribed-temperature multiplier.

use super::{ScanSensitivityError, add, finite, invalid, poll, zeros};
use crate::bc::history::ThermalBoundaryHistory;
use crate::bc::history::march::moving::sensitivity::schedule::{count, pull_signal};
use fs_exec::Cx;
use fs_scenario::bc::{BcKind, BcValue};

/// Derivatives of one explicitly authored prescribed-temperature row.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanPrescribedTemperatureGradient {
    /// Original history declaration order, restricted to Dirichlet histories.
    pub region: String,
    /// Per-kelvin derivatives in the original scalar/signal representation:
    /// one uniform/constant value, ramp endpoints, table values, or stored
    /// Chebyshev coefficients. Static partition values and knots stay fixed.
    pub values_k: Vec<f64>,
}

pub(in crate::bc::history::march::enthalpy) fn control_count(
    cx: &Cx<'_>,
    boundary: &ThermalBoundaryHistory<'_>,
) -> Result<usize, ScanSensitivityError> {
    let mut total = 0usize;
    for (i, (_, history)) in boundary.rows.iter().enumerate() {
        poll(cx, i)?;
        if history.row.kind != BcKind::Dirichlet {
            continue;
        }
        let size = match &history.row.value {
            Some(BcValue::Uniform(_)) => 1,
            Some(BcValue::Signal(signal)) => count(signal)?,
            _ => {
                return Err(invalid(
                    "scan prescribed-temperature row has no scalar history",
                ));
            }
        };
        total = total
            .checked_add(size)
            .ok_or_else(|| invalid("scan prescribed-temperature control count overflows"))?;
    }
    Ok(total)
}

pub(in crate::bc::history::march::enthalpy) fn pullback(
    cx: &Cx<'_>,
    history: &ThermalBoundaryHistory<'_>,
    time_s: f64,
    before: bool,
    nodal: &[f64],
) -> Result<Vec<ScanPrescribedTemperatureGradient>, ScanSensitivityError> {
    poll(cx, 0)?;
    if nodal.len() != history.mesh().vertex_count() {
        return Err(invalid(
            "scan prescribed-temperature gradient changed its original mesh",
        ));
    }
    if !history
        .rows
        .iter()
        .any(|(_, row)| row.row.kind == BcKind::Dirichlet)
    {
        return Ok(Vec::new());
    }
    let mut region_seeds = zeros(history.partition.region_names().len())?;
    for (i, (vertex, region)) in history.partition.prescribed_regions().enumerate() {
        if i % 256 == 0 {
            poll(cx, i)?;
        }
        region_seeds[region] = finite(region_seeds[region] + finite(nodal[vertex])?)?;
    }
    let mut output = Vec::new();
    for (i, (region, scheduled)) in history.rows.iter().enumerate() {
        poll(cx, i)?;
        if scheduled.row.kind != BcKind::Dirichlet {
            continue;
        }
        // The immutable history has already evaluated and admitted this exact
        // endpoint side through its original physical boundary owner.
        let seed = region_seeds[*region];
        let values_k = match &scheduled.row.value {
            Some(BcValue::Uniform(_)) => vec![seed],
            Some(BcValue::Signal(signal)) => pull_signal(cx, signal, time_s, before, seed)?,
            _ => {
                return Err(invalid(
                    "scan prescribed-temperature row has no scalar history",
                ));
            }
        };
        output
            .try_reserve(1)
            .map_err(|_| invalid("scan prescribed-temperature gradient allocation refused"))?;
        output.push(ScanPrescribedTemperatureGradient {
            region: scheduled.row.region.clone(),
            values_k,
        });
    }
    poll(cx, output.len())?;
    Ok(output)
}

pub(in crate::bc::history::march::enthalpy) fn accumulate(
    cx: &Cx<'_>,
    total: &mut Vec<ScanPrescribedTemperatureGradient>,
    values: Vec<ScanPrescribedTemperatureGradient>,
) -> Result<(), ScanSensitivityError> {
    if total.is_empty() {
        *total = values;
        return Ok(());
    }
    if total.len() != values.len() {
        return Err(invalid(
            "scan prescribed-temperature history layout changed during reverse sweep",
        ));
    }
    for (i, (a, b)) in total.iter_mut().zip(values).enumerate() {
        poll(cx, i)?;
        if a.region != b.region {
            return Err(invalid(
                "scan prescribed-temperature region order changed during reverse sweep",
            ));
        }
        add(cx, &mut a.values_k, &b.values_k)?;
    }
    Ok(())
}
