//! Steady adiabatic water throttling between two saturation pressures.
//!
//! With no shaft work and negligible changes in kinetic/potential energy,
//! the open-system first law preserves specific enthalpy. Upstream pressure
//! energy can create downstream vapor without an added heat source. This is
//! the equilibrium endpoint of a throttle, not a valve flow-rate, choking,
//! transient nucleation or two-phase transport model. Both endpoints must
//! lie on the original source's liquid/vapor coexistence intervals.

use super::{IsobaricWaterVaporization, WaterVaporizationError, WaterVaporizationState};
use fs_qty::Pressure;
use fs_thermochem::water::{SaturatedWaterState, WaterSaturationError, WaterSaturationSolve};

/// Original pressure-inverse work limit and explicit physical accuracy budgets.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterFlashControls {
    /// Same original saturation inversion used by the isobaric material law.
    pub saturation: WaterSaturationSolve,
    /// Absolute independent delta_u + delta_(p*v) residual budget [J/kg].
    pub energy_tolerance_j_kg: f64,
    /// Maximum negative source entropy change due to representation [J/(kg K)].
    /// The reported signed change is never clamped or repaired.
    pub entropy_tolerance_j_kg_k: f64,
}

/// Immutable steady flash proposal per kilogram passing through the throttle.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterFlash {
    /// Original upstream liquid/vapor inventory and caloric reference.
    pub inlet: WaterVaporizationState,
    /// Downstream equilibrium state, with exactly the same represented h.
    pub outlet: WaterVaporizationState,
    /// Actual original-source upstream minus downstream pressure [Pa].
    pub pressure_drop_pa: f64,
    /// Evaluated downstream pressure minus requested pressure [Pa].
    pub outlet_pressure_residual_pa: f64,
    /// Actual downstream minus upstream specific internal energy [J/kg].
    pub internal_energy_change_j_kg: f64,
    /// Actual downstream minus upstream specific flow work p*v [J/kg].
    /// This is not the closed-system constant-pressure boundary work p*delta_v.
    pub flow_work_change_j_kg: f64,
    /// Independent source delta_u + delta_(p*v) [J/kg], which must vanish.
    pub first_law_residual_j_kg: f64,
    /// Source downstream minus upstream specific entropy [J/(kg K)]. For this
    /// adiabatic steady throttle it is the generated entropy per kg of throughput.
    pub entropy_generation_j_kg_k: f64,
}

/// Explicit refusal; no inlet state or accepted flow history is changed.
#[derive(Clone, Debug, PartialEq)]
pub enum WaterFlashError {
    /// Only a positive finite pressure strictly below the inlet is admitted.
    PressureDropRequired,
    /// Energy/entropy accuracy budgets must be finite and nonnegative.
    InvalidControls,
    /// Original saturation pressure inversion refused.
    Saturation(WaterSaturationError),
    /// The original source cannot represent this enthalpy in coexistence.
    Coexistence(WaterVaporizationError),
    /// The independent open-system first-law residual exceeded its J/kg budget.
    EnergyBalance { residual_j_kg: f64 },
    /// The proposed adiabatic pressure drop decreases entropy beyond its budget.
    EntropyDecrease { change_j_kg_k: f64 },
    /// A required source difference or product is not finite.
    Unrepresentable,
    /// Caller cancellation; no proposal is published.
    Cancelled,
}

impl std::fmt::Display for WaterFlashError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "water throttling flash refusal: {self:?}")
    }
}
impl std::error::Error for WaterFlashError {}

impl WaterVaporizationState {
    /// Resolve a steady adiabatic pressure drop without inventing heat or work.
    ///
    /// The original source supplies the downstream saturation temperature and
    /// endpoints. The unchanged inlet h selects its vapor quality. Original
    /// mixture u, p*v and s independently check the open-system first and second
    /// laws; the numerical tolerances do not bound source uncertainty. A state
    /// requiring subcooling or superheating explicitly refuses. No flash result
    /// is published on source failure, out-of-domain heat, budget or cancellation.
    pub fn flash_to_pressure(
        self,
        outlet_pressure: Pressure,
        controls: WaterFlashControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterFlash, WaterFlashError> {
        use WaterFlashError as Error;
        if cancelled() {
            return Err(Error::Cancelled);
        }
        if !controls.energy_tolerance_j_kg.is_finite()
            || controls.energy_tolerance_j_kg < 0.0
            || !controls.entropy_tolerance_j_kg_k.is_finite()
            || controls.entropy_tolerance_j_kg_k < 0.0
        {
            return Err(Error::InvalidControls);
        }
        if !outlet_pressure.value().is_finite()
            || outlet_pressure.value() <= 0.0
            || outlet_pressure.value() >= self.pressure_pa()
        {
            return Err(Error::PressureDropRequired);
        }
        let source =
            SaturatedWaterState::at_pressure(outlet_pressure, controls.saturation, &mut cancelled)
                .map_err(Error::Saturation)?;
        let law = IsobaricWaterVaporization::from_saturation(*source.state())
            .map_err(Error::Coexistence)?;
        let outlet = law
            .state_at_specific_enthalpy(self.specific_enthalpy_j_kg())
            .map_err(Error::Coexistence)?;
        let pressure_drop_pa = self.pressure_pa() - outlet.pressure_pa();
        if pressure_drop_pa <= 0.0 {
            return Err(Error::PressureDropRequired);
        }
        let internal_energy_change_j_kg =
            outlet.specific_internal_energy_j_kg() - self.specific_internal_energy_j_kg();
        let flow_work_change_j_kg = outlet.pressure_pa() * outlet.specific_volume_m3_kg()
            - self.pressure_pa() * self.specific_volume_m3_kg();
        let first_law_residual_j_kg = internal_energy_change_j_kg + flow_work_change_j_kg;
        let entropy_generation_j_kg_k =
            outlet.specific_entropy_j_kg_k() - self.specific_entropy_j_kg_k();
        if ![
            pressure_drop_pa,
            internal_energy_change_j_kg,
            flow_work_change_j_kg,
            first_law_residual_j_kg,
            entropy_generation_j_kg_k,
        ]
        .into_iter()
        .all(f64::is_finite)
        {
            return Err(Error::Unrepresentable);
        }
        if first_law_residual_j_kg.abs() > controls.energy_tolerance_j_kg {
            return Err(Error::EnergyBalance {
                residual_j_kg: first_law_residual_j_kg,
            });
        }
        if entropy_generation_j_kg_k < -controls.entropy_tolerance_j_kg_k {
            return Err(Error::EntropyDecrease {
                change_j_kg_k: entropy_generation_j_kg_k,
            });
        }
        if cancelled() {
            return Err(Error::Cancelled);
        }
        Ok(WaterFlash {
            inlet: self,
            outlet,
            pressure_drop_pa,
            outlet_pressure_residual_pa: source.pressure_residual().value(),
            internal_energy_change_j_kg,
            flow_work_change_j_kg,
            first_law_residual_j_kg,
            entropy_generation_j_kg_k,
        })
    }
}

#[cfg(test)]
mod tests;
