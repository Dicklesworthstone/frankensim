//! Source-bound planar force and couple impulses, per unit out-of-plane depth.

use super::{
    MaterialFreeSurface2, MaterialFreeSurfaceImpulseControls2, MaterialFreeSurfaceImpulseError,
    MaterialFreeSurfaceImpulseReport2, MaterialLbmError, finite, positive, prepare_controls,
    receipt_si, represented,
};
use crate::impulse::{
    ImpulseError, WrenchImpulse2, WrenchImpulseControls2, WrenchImpulseReport2, poll,
};
use fs_exec::CancelGate;

/// A total planar wrench applied TO liquid, per unit out-of-plane depth.
/// No thickness or three-dimensional side-wall model is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFreeSurfaceWrenchImpulse2 {
    /// Origin measured from the fixed domain's lower corner [m].
    pub origin_m: [f64; 2],
    /// Linear impulse per depth [kg/s = N s/m], not impulse density or force.
    pub linear_impulse_kg_s: [f64; 2],
    /// Out-of-plane angular impulse per depth about origin [kg m/s].
    /// Positive is counterclockwise; an independent free couple is included.
    pub angular_impulse_kg_m_s: f64,
}

/// Original per-depth impulse admission plus a separate angular allowance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFreeSurfaceWrenchImpulseControls2 {
    pub impulse: MaterialFreeSurfaceImpulseControls2,
    /// Deposition and local/global angular readback defects per depth [kg m/s].
    pub maximum_angular_momentum_defect_kg_m_s: f64,
}

/// Actual planar inertia, linear/angular impulse and signed work, all per depth.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct MaterialFreeSurfaceWrenchImpulseReport2 {
    pub requested: MaterialFreeSurfaceWrenchImpulse2,
    /// Linear impulse [kg/s], mass [kg/m] and signed work [J/m].
    pub impulse: MaterialFreeSurfaceImpulseReport2,
    pub stencil_mass_kg_m: f64,
    pub stencil_center_of_mass_m: [f64; 2],
    pub applied_angular_impulse_kg_m_s: f64,
    pub angular_momentum_change_kg_m_s: f64,
    /// Includes deposition error relative to the original authored SI impulse.
    pub absolute_linear_momentum_defect_kg_s: f64,
    /// Includes origin and unit-conversion rounding in the original SI chart.
    pub absolute_angular_momentum_defect_kg_m_s: f64,
}

impl MaterialFreeSurface2<'_> {
    /// Apply a total physical planar wrench per unit out-of-plane depth.
    ///
    /// Strictly increasing indices use the original x-fastest cell order,
    /// including wall rows. The original material density and cell size supply
    /// inertia; interface cells use their retained partial liquid mass.
    /// Translation plus rotation produces the least mass-weighted velocity
    /// increment with the requested force-and-couple impulse.
    ///
    /// Linear impulse is [kg/s], angular impulse [kg m/s], and work [J/m].
    /// No timestep multiplies the input and no out-of-plane thickness is chosen.
    /// Density, fill, source identity, held Guo forcing and physical viscous
    /// stress survive. Transport, pressure reconstruction and heat do not advance.
    ///
    /// Both SI balances are admitted about the ORIGINAL authored physical
    /// origin before native publication, including origin-conversion rounding.
    /// Native, unit, budget or cancellation refusal preserves accepted state.
    /// The caller supplies any reciprocal body reaction. Only the native compact
    /// stencil and existing population scratch are used; no dense impulse field
    /// or whole-owner clone is allocated.
    pub fn apply_wrench_impulse(
        &mut self,
        cells: &[usize],
        wrench: MaterialFreeSurfaceWrenchImpulse2,
        controls: MaterialFreeSurfaceWrenchImpulseControls2,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceWrenchImpulseReport2, MaterialFreeSurfaceImpulseError> {
        poll(gate)?;
        let count = self
            .lattice
            .grid
            .nx
            .checked_mul(self.lattice.grid.ny)
            .ok_or(ImpulseError::Invalid(
                "physical wrench cell count overflows",
            ))?;
        let limits = controls.impulse;
        let scaling = *self.scaling();
        let (impulse, units) = prepare_controls(
            scaling,
            self.mass_per_depth_unit_kg_m,
            count,
            limits.max_cells,
            [
                limits.maximum_speed_m_s,
                limits.maximum_momentum_defect_kg_s,
                limits.maximum_work_defect_j_m,
            ],
        )?;
        let length = positive(scaling.cell_size_m, "liquid length unit is unrepresentable")?;
        let angular_unit = positive(
            units[1] * length,
            "planar liquid angular-momentum unit is unrepresentable",
        )?;
        let budget = controls.maximum_angular_momentum_defect_kg_m_s;
        if !budget.is_finite() || budget < 0.0 {
            return Err(MaterialLbmError::Invalid(
                "physical wrench angular budget must be finite and nonnegative",
            )
            .into());
        }
        let mut native = WrenchImpulse2 {
            origin: [0.0; 2],
            linear: [0.0; 2],
            angular: 0.0,
        };
        let mut origin_shift = [0.0; 2];
        for a in 0..2 {
            native.origin[a] = represented(wrench.origin_m[a], wrench.origin_m[a] / length)?;
            native.linear[a] = represented(
                wrench.linear_impulse_kg_s[a],
                wrench.linear_impulse_kg_s[a] / units[1],
            )?;
            origin_shift[a] = finite(
                native.origin[a].mul_add(length, -wrench.origin_m[a]),
                "physical wrench origin conversion is unrepresentable",
            )?;
        }
        // The native origin represents a nearby physical point. Transport the
        // authored scalar moment there; receipt() transports the actual result
        // back rather than relabeling it as a moment about the authored origin.
        let angular = finite(
            wrench.angular_impulse_kg_m_s - cross(origin_shift, wrench.linear_impulse_kg_s),
            "physical wrench origin transport is unrepresentable",
        )?;
        native.angular = represented(angular, angular / angular_unit)?;
        let native_controls = WrenchImpulseControls2 {
            impulse,
            maximum_angular_momentum_defect: represented(budget, budget / angular_unit)?,
        };
        let mut refused_units = None;
        let result = self.lattice.apply_wrench_impulse_with_receipt(
            cells,
            native,
            native_controls,
            gate,
            |report| {
                receipt(
                    report,
                    wrench,
                    controls,
                    units,
                    length,
                    angular_unit,
                    origin_shift,
                )
                .map_err(|error| {
                    refused_units = Some(error);
                    ImpulseError::Invalid("physical planar liquid wrench receipt refused")
                })
            },
        );
        if let Some(error) = refused_units {
            return Err(error.into());
        }
        Ok(result?)
    }
}

fn difference(a: [f64; 2], b: [f64; 2]) -> f64 {
    (a[0] - b[0]).abs() + (a[1] - b[1]).abs()
}

fn cross(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0].mul_add(b[1], -a[1] * b[0])
}

/// Re-evaluate physical balances in the authored SI chart before publication.
#[allow(clippy::too_many_arguments)]
fn receipt(
    report: WrenchImpulseReport2,
    requested: MaterialFreeSurfaceWrenchImpulse2,
    controls: MaterialFreeSurfaceWrenchImpulseControls2,
    units: [f64; 3],
    length: f64,
    angular_unit: f64,
    origin_shift: [f64; 2],
) -> Result<MaterialFreeSurfaceWrenchImpulseReport2, MaterialLbmError> {
    let local = (report.absolute_angular_momentum_defect
        - (report.angular_momentum_change - report.applied_angular_impulse).abs()
        - (report.applied_angular_impulse - report.requested.angular).abs())
    .max(0.0);
    let limits = controls.impulse;
    let linear = receipt_si(
        report.impulse,
        units,
        limits.maximum_momentum_defect_kg_s,
        limits.maximum_work_defect_j_m,
    )?;
    let applied = finite(
        represented(
            report.applied_angular_impulse,
            report.applied_angular_impulse * angular_unit,
        )? + cross(origin_shift, linear.applied_impulse),
        "physical applied planar torque is unrepresentable",
    )?;
    let change = finite(
        represented(
            report.angular_momentum_change,
            report.angular_momentum_change * angular_unit,
        )? + cross(origin_shift, linear.momentum_change),
        "physical observed planar torque is unrepresentable",
    )?;
    let mut center = [0.0; 2];
    for a in 0..2 {
        center[a] = represented(
            report.stencil_center_of_mass[a],
            report.stencil_center_of_mass[a] * length,
        )?;
    }
    let linear_defect = finite(
        linear.absolute_momentum_defect
            + difference(linear.applied_impulse, requested.linear_impulse_kg_s),
        "physical planar wrench linear readback is unrepresentable",
    )?;
    let angular_defect = finite(
        represented(local, local * angular_unit)?
            + (change - applied).abs()
            + (applied - requested.angular_impulse_kg_m_s).abs(),
        "physical planar wrench angular readback is unrepresentable",
    )?;
    if linear_defect > limits.maximum_momentum_defect_kg_s
        || angular_defect > controls.maximum_angular_momentum_defect_kg_m_s
    {
        return Err(MaterialLbmError::Invalid(
            "physical planar liquid wrench exceeds its SI balance budget",
        ));
    }
    Ok(MaterialFreeSurfaceWrenchImpulseReport2 {
        requested,
        impulse: MaterialFreeSurfaceImpulseReport2 {
            applied_impulse_kg_s: linear.applied_impulse,
            momentum_change_kg_s: linear.momentum_change,
            mass_change_kg_m: linear.mass_change,
            work_j_m: linear.work,
            kinetic_energy_change_j_m: linear.kinetic_energy_change,
            absolute_momentum_defect_kg_s: linear.absolute_momentum_defect,
            absolute_work_defect_j_m: linear.absolute_work_defect,
            affected_cells: linear.affected_cells,
        },
        stencil_mass_kg_m: represented(report.stencil_mass, report.stencil_mass * units[0])?,
        stencil_center_of_mass_m: center,
        applied_angular_impulse_kg_m_s: applied,
        angular_momentum_change_kg_m_s: change,
        absolute_linear_momentum_defect_kg_s: linear_defect,
        absolute_angular_momentum_defect_kg_m_s: angular_defect,
    })
}

#[cfg(test)]
mod tests;
