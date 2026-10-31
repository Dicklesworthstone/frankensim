//! Accepted-duration accounting for the original material/enthalpy owner.

use super::{DepositionError, ViscousHeatControls, ViscousPower, invalid};

/// Mechanical actuator supply and the same molecular/porous heat deposited into h.
/// This is one accepted advance's account, not an additional thermal input.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ViscousHeatingAudit {
    pub duration_s: f64,
    pub substeps: usize,
    /// Signed work at each Face3::ALL wall [J].
    pub wall_work_j: [f64; 6],
    pub supplied_wall_work_j: f64,
    /// Signed reservoir normal-traction work at each Face3::ALL mouth [J].
    pub pressure_work_j: [f64; 6],
    /// Mechanical supply, never an extra thermal input beside generated_heat_j.
    pub supplied_pressure_work_j: f64,
    /// Signed prescribed body-force work, separate from wall supply.
    /// This is not an additional thermal input beside generated_heat_j.
    pub supplied_body_work_j: f64,
    pub generated_heat_j: f64,
    /// Subaccounts of generated_heat_j, never additional thermal inputs.
    pub molecular_heat_j: f64,
    pub darcy_heat_j: f64,
    pub forchheimer_heat_j: f64,
    /// Diagnostic only: it is bounded near zero and is never deposited as heat.
    pub pressure_constraint_work_j: f64,
    pub residual_j: f64,
}

#[derive(Clone)]
pub(super) struct ViscousHeating {
    pub controls: ViscousHeatControls,
    pub pressure_work: bool,
    pub audit: ViscousHeatingAudit,
    last_end_s: Option<f64>,
}

impl ViscousHeating {
    pub fn new(controls: ViscousHeatControls) -> Self {
        Self {
            controls,
            pressure_work: false,
            audit: ViscousHeatingAudit::default(),
            last_end_s: None,
        }
    }

    pub fn begin(&mut self) {
        self.audit = ViscousHeatingAudit::default();
        self.last_end_s = None;
    }

    pub fn accept(
        &mut self,
        power: &ViscousPower,
        start_s: f64,
        dt: f64,
        end_s: f64,
    ) -> Result<(), DepositionError> {
        if !dt.is_finite()
            || dt <= 0.0
            || !start_s.is_finite()
            || !end_s.is_finite()
            || end_s <= start_s
            || self
                .last_end_s
                .is_some_and(|last| last.to_bits() != start_s.to_bits())
            || ((end_s - start_s) - dt).abs() > 64.0 * f64::EPSILON * end_s.abs().max(dt)
        {
            return Err(invalid(
                "viscous work must use one actual nonrepeated thermal interval",
            ));
        }
        // Candidate arithmetic completes before replacing even this local audit.
        let mut next = self.audit;
        let add = |old: &mut f64, value: f64| -> Result<(), DepositionError> {
            *old += value;
            if old.is_finite() {
                Ok(())
            } else {
                Err(DepositionError::Arithmetic)
            }
        };
        add(&mut next.duration_s, dt)?;
        next.substeps = next
            .substeps
            .checked_add(1)
            .ok_or(DepositionError::Budget("viscous work substeps"))?;
        for (sum, watts) in next.wall_work_j.iter_mut().zip(power.wall_power_w) {
            add(sum, dt * watts)?;
        }
        for (sum, watts) in next.pressure_work_j.iter_mut().zip(power.pressure_power_w) {
            add(sum, dt * watts)?;
        }
        for (sum, watts) in [
            (&mut next.supplied_wall_work_j, power.supplied_wall_power_w),
            (&mut next.supplied_pressure_work_j, power.supplied_pressure_power_w),
            (&mut next.supplied_body_work_j, power.supplied_body_power_w),
            (&mut next.generated_heat_j, power.heat_w),
            (&mut next.molecular_heat_j, power.molecular_heat_w),
            (&mut next.darcy_heat_j, power.darcy_heat_w),
            (&mut next.forchheimer_heat_j, power.forchheimer_heat_w),
            (
                &mut next.pressure_constraint_work_j,
                power.pressure_constraint_power_w,
            ),
            (&mut next.residual_j, power.residual_w),
        ] {
            add(sum, dt * watts)?;
        }
        let allowance = self.controls.power_tolerance_w * next.duration_s;
        if !allowance.is_finite()
            || next.residual_j.abs() > allowance
            || next.pressure_constraint_work_j.abs() > allowance
            || (next.supplied_wall_work_j + next.supplied_pressure_work_j
                + next.supplied_body_work_j - next.generated_heat_j).abs()
                > allowance
            || (next.generated_heat_j - next.molecular_heat_j - next.darcy_heat_j
                - next.forchheimer_heat_j).abs() > allowance
        {
            return Err(invalid(
                "integrated wall/reservoir/body work and physical heat do not close",
            ));
        }
        self.audit = next;
        self.last_end_s = Some(end_s);
        Ok(())
    }
}
