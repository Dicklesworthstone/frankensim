//! Physical force-and-couple impulses on the retained sparse liquid inventory.

use super::{
    MaterialFreeSurface3, MaterialFreeSurfaceImpulseControls3, MaterialFreeSurfaceImpulseError,
    MaterialFreeSurfaceImpulseReport3, MaterialLbmError, finite, positive, prepare_controls,
    receipt_si, represented,
};
use crate::impulse::{
    ImpulseError, WrenchImpulse3, WrenchImpulseControls3, WrenchImpulseReport3, poll,
};
use fs_exec::CancelGate;

/// One total instantaneous wrench applied TO liquid, in three-dimensional SI.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFreeSurfaceWrenchImpulse3 {
    /// Origin measured from the fixed domain's lower corner [m].
    pub origin_m: [f64; 3],
    /// Total linear impulse [kg m/s], not impulse density or force.
    pub linear_impulse_kg_m_s: [f64; 3],
    /// Total angular impulse about origin, including a free couple [kg m²/s].
    pub angular_impulse_kg_m2_s: [f64; 3],
}

/// Original physical impulse admission plus an angular-momentum allowance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFreeSurfaceWrenchImpulseControls3 {
    pub impulse: MaterialFreeSurfaceImpulseControls3,
    /// Deposition and local/global physical angular readback defects [kg m²/s].
    pub maximum_angular_momentum_defect_kg_m2_s: f64,
}

/// Actual liquid inertia, impulse, angular impulse and signed mechanical work.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct MaterialFreeSurfaceWrenchImpulseReport3 {
    pub requested: MaterialFreeSurfaceWrenchImpulse3,
    pub impulse: MaterialFreeSurfaceImpulseReport3,
    pub stencil_mass_kg: f64,
    pub stencil_center_of_mass_m: [f64; 3],
    pub applied_angular_impulse_kg_m2_s: [f64; 3],
    pub angular_momentum_change_kg_m2_s: [f64; 3],
    /// Includes deposition error relative to the original authored SI impulse.
    pub absolute_linear_momentum_defect_kg_m_s: f64,
    /// Includes origin and unit-conversion rounding in the original SI chart.
    pub absolute_angular_momentum_defect_kg_m2_s: f64,
}

impl MaterialFreeSurface3<'_> {
    /// Apply one total physical wrench on a compact stencil of original liquid.
    ///
    /// Strictly increasing indices use the native x-fastest domain order. The
    /// original cell size and material density convert total linear/angular
    /// impulse; no timestep or whole-domain impulse field is introduced.
    /// Partial-fill inertia is the actual retained liquid mass. The shared
    /// native solve preserves density, fill, force-corrected viscous stress and
    /// source bindings; pressure, transport, time and heat do not advance.
    ///
    /// Both SI linear and angular receipts are admitted before publication,
    /// including conversion rounding about the authored physical origin.
    /// Any native, unit, budget or cancellation refusal keeps the complete
    /// accepted owner intact. The coupling owner supplies any body reaction.
    pub fn apply_wrench_impulse(
        &mut self,
        cells: &[usize],
        wrench: MaterialFreeSurfaceWrenchImpulse3,
        controls: MaterialFreeSurfaceWrenchImpulseControls3,
        gate: &CancelGate,
    ) -> Result<MaterialFreeSurfaceWrenchImpulseReport3, MaterialFreeSurfaceImpulseError> {
        poll(gate)?;
        let count = self
            .lattice
            .domain_cells()
            .into_iter()
            .try_fold(1usize, |n, edge| n.checked_mul(edge))
            .ok_or(ImpulseError::Invalid(
                "physical wrench cell count overflows",
            ))?;
        let limits = controls.impulse;
        let (impulse, units) = prepare_controls(
            self.scaling,
            self.scaling.mass_unit_kg,
            count,
            limits.max_cells,
            [
                limits.maximum_speed_m_s,
                limits.maximum_momentum_defect_kg_m_s,
                limits.maximum_work_defect_j,
            ],
        )?;
        let length = positive(
            self.scaling.cell_size_m,
            "liquid length unit is unrepresentable",
        )?;
        let angular_unit = positive(
            units[1] * length,
            "liquid angular-momentum unit is unrepresentable",
        )?;
        let budget = controls.maximum_angular_momentum_defect_kg_m2_s;
        if !budget.is_finite() || budget < 0.0 {
            return Err(MaterialLbmError::Invalid(
                "physical wrench angular budget must be finite and nonnegative",
            )
            .into());
        }
        let mut native = WrenchImpulse3 {
            origin: [0.0; 3],
            linear: [0.0; 3],
            angular: [0.0; 3],
        };
        let mut origin_shift = [0.0; 3];
        for a in 0..3 {
            native.origin[a] = represented(wrench.origin_m[a], wrench.origin_m[a] / length)?;
            native.linear[a] = represented(
                wrench.linear_impulse_kg_m_s[a],
                wrench.linear_impulse_kg_m_s[a] / units[1],
            )?;
            origin_shift[a] = finite(
                native.origin[a].mul_add(length, -wrench.origin_m[a]),
                "physical wrench origin conversion is unrepresentable",
            )?;
        }
        // Transport the authored moment to the representable native origin.
        // The dependent SI readback transports it back and checks the actual
        // result, including all remaining conversion/population roundoff.
        let moment_shift = cross(origin_shift, wrench.linear_impulse_kg_m_s);
        for a in 0..3 {
            let angular = finite(
                wrench.angular_impulse_kg_m2_s[a] - moment_shift[a],
                "physical wrench origin transport is unrepresentable",
            )?;
            native.angular[a] = represented(angular, angular / angular_unit)?;
        }
        let native_controls = WrenchImpulseControls3 {
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
                    ImpulseError::Invalid("physical liquid wrench receipt refused")
                })
            },
        );
        if let Some(error) = refused_units {
            return Err(error.into());
        }
        Ok(result?)
    }
}

fn difference(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|axis| (a[axis] - b[axis]).abs()).sum()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -a[2] * b[1]),
        a[2].mul_add(b[0], -a[0] * b[2]),
        a[0].mul_add(b[1], -a[1] * b[0]),
    ]
}

/// Re-evaluate the receipts in the authored SI chart before native publication.
#[allow(clippy::too_many_arguments)]
fn receipt(
    report: WrenchImpulseReport3,
    requested: MaterialFreeSurfaceWrenchImpulse3,
    controls: MaterialFreeSurfaceWrenchImpulseControls3,
    units: [f64; 3],
    length: f64,
    angular_unit: f64,
    origin_shift: [f64; 3],
) -> Result<MaterialFreeSurfaceWrenchImpulseReport3, MaterialLbmError> {
    // Isolate local angular readback. Global/deposition errors are recomputed
    // after conversion, so they cannot disappear through SI rounding.
    let local = (report.absolute_angular_momentum_defect
        - difference(
            report.angular_momentum_change,
            report.applied_angular_impulse,
        )
        - difference(report.applied_angular_impulse, report.requested.angular))
    .max(0.0);
    let limits = controls.impulse;
    let linear = receipt_si(
        report.impulse,
        units,
        limits.maximum_momentum_defect_kg_m_s,
        limits.maximum_work_defect_j,
    )?;
    let applied_shift = cross(origin_shift, linear.applied_impulse);
    let observed_shift = cross(origin_shift, linear.momentum_change);
    let mut applied = [0.0; 3];
    let mut change = [0.0; 3];
    let mut center = [0.0; 3];
    for a in 0..3 {
        // Native torque is about its representable lattice origin. Shift to
        // the ORIGINAL authored SI origin, rather than relabeling that torque.
        applied[a] = finite(
            represented(
                report.applied_angular_impulse[a],
                report.applied_angular_impulse[a] * angular_unit,
            )? + applied_shift[a],
            "physical applied torque is unrepresentable",
        )?;
        change[a] = finite(
            represented(
                report.angular_momentum_change[a],
                report.angular_momentum_change[a] * angular_unit,
            )? + observed_shift[a],
            "physical observed torque is unrepresentable",
        )?;
        center[a] = represented(
            report.stencil_center_of_mass[a],
            report.stencil_center_of_mass[a] * length,
        )?;
    }
    let linear_defect = finite(
        linear.absolute_momentum_defect
            + difference(linear.applied_impulse, requested.linear_impulse_kg_m_s),
        "physical wrench linear readback is unrepresentable",
    )?;
    let angular_defect = finite(
        represented(local, local * angular_unit)?
            + difference(change, applied)
            + difference(applied, requested.angular_impulse_kg_m2_s),
        "physical wrench angular readback is unrepresentable",
    )?;
    if linear_defect > limits.maximum_momentum_defect_kg_m_s
        || angular_defect > controls.maximum_angular_momentum_defect_kg_m2_s
    {
        return Err(MaterialLbmError::Invalid(
            "physical liquid wrench exceeds its SI balance budget",
        ));
    }
    Ok(MaterialFreeSurfaceWrenchImpulseReport3 {
        requested,
        impulse: MaterialFreeSurfaceImpulseReport3 {
            applied_impulse_kg_m_s: linear.applied_impulse,
            momentum_change_kg_m_s: linear.momentum_change,
            mass_change_kg: linear.mass_change,
            work_j: linear.work,
            kinetic_energy_change_j: linear.kinetic_energy_change,
            absolute_momentum_defect_kg_m_s: linear.absolute_momentum_defect,
            absolute_work_defect_j: linear.absolute_work_defect,
            affected_cells: linear.affected_cells,
        },
        stencil_mass_kg: represented(report.stencil_mass, report.stencil_mass * units[0])?,
        stencil_center_of_mass_m: center,
        applied_angular_impulse_kg_m2_s: applied,
        angular_momentum_change_kg_m2_s: change,
        absolute_linear_momentum_defect_kg_m_s: linear_defect,
        absolute_angular_momentum_defect_kg_m2_s: angular_defect,
    })
}

#[cfg(test)]
mod tests;
