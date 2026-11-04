//! Exact affine-flow stream work for the total-pressure chart.
//!
//! For q(t) affine, p(t)=P-rho*q(t)^2/(2*A^2). Traction uses mean(q^2),
//! while pressure work uses mean(q^3); multiplying mean pressure by mean flow
//! would silently lose work during acceleration and reversal.

use super::super::loads::pipe::{PipeLoadError, affine_mean_cube, affine_mean_square};

type Result<T> = core::result::Result<T, PipeLoadError>;

/// Both physical ports of one original branch, ordered from/to. Positive work
/// enters that pipe. The section is uniform, so its two kinetic streams cancel
/// in the pipe's own energy equation, but not at unequal-area junctions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BranchStreamWork {
    pub average_static_pressures_pa: [f64; 2],
    pub static_pressure_work_into_branch_j: [f64; 2],
    pub kinetic_energy_into_branch_j: [f64; 2],
    pub total_pressure_work_into_branch_j: [f64; 2],
}

/// Work leaving a hydraulic node through its complete incident branch set.
/// At a zero-storage junction, pressure work on the fluid is the NEGATIVE of
/// static_pressure_work_out_j and balances kinetic_energy_out_j. The retained
/// residual also exposes finite native continuity error at a large datum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NodeStreamWork {
    pub node: usize,
    pub static_pressure_work_out_j: f64,
    pub kinetic_energy_out_j: f64,
    pub total_pressure_work_out_j: f64,
    pub pressure_kinetic_identity_residual_j: f64,
    pub zero_storage_balance_residual_j: f64,
}

pub(super) fn number(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(PipeLoadError::Invalid(
            "unrepresentable total-pressure stream work",
        ))
    }
}

pub(super) fn mean_dynamic_pressure(flows: [f64; 2], density: f64, area: f64) -> Result<f64> {
    number(0.5 * density * (affine_mean_square(flows)? / area) / area)
}

pub(super) fn branch_work(
    duration: f64,
    density: f64,
    area: f64,
    total_pressures: [f64; 2],
    flows: [f64; 2],
) -> Result<BranchStreamWork> {
    if !(duration.is_finite() && duration > 0.0)
        || !(density.is_finite() && density > 0.0)
        || !(area.is_finite() && area > 0.0)
        || total_pressures.iter().chain(&flows).any(|v| !v.is_finite())
    {
        return Err(PipeLoadError::Invalid(
            "finite positive stream geometry and interval required",
        ));
    }
    let head = mean_dynamic_pressure(flows, density, area)?;
    let volume = number(duration * number(0.5 * flows[0] + 0.5 * flows[1])?)?;
    let kinetic = number(0.5 * duration * density * (affine_mean_cube(flows)? / area) / area)?;
    let total = [
        number(total_pressures[0] * volume)?,
        number(-total_pressures[1] * volume)?,
    ];
    Ok(BranchStreamWork {
        average_static_pressures_pa: [
            number(total_pressures[0] - head)?,
            number(total_pressures[1] - head)?,
        ],
        static_pressure_work_into_branch_j: [
            number(total[0] - kinetic)?,
            number(total[1] + kinetic)?,
        ],
        kinetic_energy_into_branch_j: [kinetic, -kinetic],
        total_pressure_work_into_branch_j: total,
    })
}

/// Full interval extrema, including the interior pressure maximum at Q=0.
/// One identical dynamic head is subtracted from BOTH ends of a branch.
pub(super) fn static_pressure_bounds(
    pressures: [f64; 2],
    flows: [f64; 2],
    density: f64,
    area: f64,
) -> Result<[[f64; 2]; 2]> {
    let speed_max = flows[0].abs().max(flows[1].abs()) / area;
    let crosses_zero = flows[0] <= 0.0 && flows[1] >= 0.0 || flows[1] <= 0.0 && flows[0] >= 0.0;
    let speed_min = if crosses_zero {
        0.0
    } else {
        flows[0].abs().min(flows[1].abs()) / area
    };
    let maximum_head = number(0.5 * density * speed_max * speed_max)?;
    let minimum_head = number(0.5 * density * speed_min * speed_min)?;
    Ok([
        [
            number(pressures[0] - maximum_head)?,
            number(pressures[1] - maximum_head)?,
        ],
        [
            number(pressures[0] - minimum_head)?,
            number(pressures[1] - minimum_head)?,
        ],
    ])
}

#[cfg(test)]
mod tests;
