//! Steady adiabatic pressure devices using the original IF97 water chart.
//!
//! Throttles preserve h. Shaft devices first resolve the original isentropic
//! outlet, then apply a caller-declared isentropic efficiency to shaft work.
//! Each actual outlet is admitted by the same p,h chart used for heating and
//! phase change. All energies are per kilogram of throughput; flow work p*v
//! is included once, separately from shaft work. Kinetic/potential changes,
//! flow-rate laws and transient machine dynamics are outside this owner.
//!
//! Cross-pressure entropy always uses the original source s. The chart's
//! pressure-specific entropy-integral seam adjustment is reported separately;
//! it is not a new global entropy reference or an empirical machine model.

use super::{IsobaricWater, WaterEntropySolve, WaterError, WaterState};
use core::fmt;

/// Absolute numerical first/second-law allowances, not source uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterFlowSolve {
    pub energy_tolerance_j_kg: f64,
    /// Maximum admitted negative original-source entropy change from rounding.
    /// The signed change in the report is never clamped to zero.
    pub entropy_tolerance_j_kg_k: f64,
}

/// Direction of steady pressure change and its isentropic-efficiency convention.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WaterShaftDevice {
    /// Lower pressure; actual outward work is efficiency times isentropic work.
    Expander,
    /// Higher pressure; actual inward work is isentropic work / efficiency.
    /// This covers a liquid pump or a vapor compressor within the source domain.
    Compressor,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterShaftSolve {
    /// Explicit constitutive input in (0,1], never a default empirical value.
    pub isentropic_efficiency: f64,
    pub entropy: WaterEntropySolve,
    pub balance: WaterFlowSolve,
}

/// Independently evaluated original-source balances per kilogram of throughput.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterFlowBalance {
    /// Positive when the fluid delivers work to an external shaft [J/kg].
    pub shaft_work_j_kg: f64,
    pub enthalpy_change_j_kg: f64,
    pub internal_energy_change_j_kg: f64,
    /// p_out*v_out - p_in*v_in; not p_out*(v_out-v_in).
    pub flow_work_change_j_kg: f64,
    /// Delta h + outward shaft work for this adiabatic device.
    pub enthalpy_residual_j_kg: f64,
    /// Independent Delta u + Delta(p*v) + outward shaft work.
    pub first_law_residual_j_kg: f64,
    /// Original s_out-s_in, equal to modeled entropy generation per throughput.
    pub source_entropy_change_j_kg_k: f64,
    /// Each chart's source saturation entropy seam defect [inlet,outlet].
    /// These remain visible but are not added to cross-pressure entropy.
    pub saturation_entropy_defects_j_kg_k: [f64; 2],
}

#[derive(Clone, Debug, PartialEq)]
pub struct WaterThrottle {
    pub inlet: WaterState,
    pub outlet: WaterState,
    pub pressure_drop_pa: f64,
    pub enthalpy_iterations: usize,
    pub balance: WaterFlowBalance,
}

#[derive(Clone, Debug, PartialEq)]
pub struct WaterShaftChange {
    pub device: WaterShaftDevice,
    pub inlet: WaterState,
    pub isentropic_outlet: WaterState,
    pub outlet: WaterState,
    pub isentropic_efficiency: f64,
    /// Signed outward work of the original-source isentropic reference [J/kg].
    pub isentropic_shaft_work_j_kg: f64,
    /// Actual s at the reference endpoint minus the requested inlet s.
    pub isentropic_entropy_residual_j_kg_k: f64,
    pub entropy_iterations: usize,
    pub enthalpy_iterations: usize,
    pub balance: WaterFlowBalance,
}

/// Refusal publishes neither a partial pressure endpoint nor partial shaft work.
#[derive(Clone, Debug, PartialEq)]
pub enum WaterFlowError {
    Invalid(&'static str),
    /// The prepared outlet pressure has the wrong sign of change for the device.
    PressureDirection,
    /// The original bounded caloric/entropy chart refused its endpoint.
    Source(WaterError),
    EnergyBalance {
        enthalpy_residual_j_kg: f64,
        first_law_residual_j_kg: f64,
    },
    EntropyDecrease {
        change_j_kg_k: f64,
    },
    Unrepresentable(&'static str),
    Cancelled,
}

impl From<WaterError> for WaterFlowError {
    fn from(error: WaterError) -> Self {
        match error {
            WaterError::Cancelled => Self::Cancelled,
            other => Self::Source(other),
        }
    }
}
impl fmt::Display for WaterFlowError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "steady water pressure device: {self:?}")
    }
}
impl std::error::Error for WaterFlowError {}

impl WaterFlowSolve {
    fn validate(self) -> Result<(), WaterFlowError> {
        if !self.energy_tolerance_j_kg.is_finite()
            || self.energy_tolerance_j_kg < 0.0
            || !self.entropy_tolerance_j_kg_k.is_finite()
            || self.entropy_tolerance_j_kg_k < 0.0
        {
            return Err(WaterFlowError::Invalid(
                "finite nonnegative energy and entropy budgets required",
            ));
        }
        Ok(())
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), WaterFlowError> {
    if cancelled() {
        Err(WaterFlowError::Cancelled)
    } else {
        Ok(())
    }
}

impl WaterState {
    /// Steady adiabatic throttling to a prepared lower-pressure water chart.
    ///
    /// Source h is conserved within the explicit numerical budget, across any
    /// admitted liquid, two-phase or superheated outlet. No heat, shaft work or
    /// pressure-dependent entropy datum is invented to fit the endpoint.
    pub fn throttle_to(
        &self,
        outlet_law: &IsobaricWater,
        controls: WaterFlowSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterThrottle, WaterFlowError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        let target_pressure = outlet_law.saturation().pressure().value();
        if target_pressure >= self.pressure_pa() {
            return Err(WaterFlowError::PressureDirection);
        }
        let solved =
            outlet_law.state_at_specific_enthalpy(self.specific_enthalpy_j_kg(), &mut cancelled)?;
        let iterations = solved.iterations();
        let outlet = solved.into_state();
        let balance = audit(self, &outlet, 0.0, controls)?;
        let pressure_drop_pa = self.pressure_pa() - outlet.pressure_pa();
        if !pressure_drop_pa.is_finite() || pressure_drop_pa <= 0.0 {
            return Err(WaterFlowError::Unrepresentable("pressure drop"));
        }
        poll(&mut cancelled)?;
        Ok(WaterThrottle {
            inlet: self.clone(),
            outlet,
            pressure_drop_pa,
            enthalpy_iterations: iterations,
            balance,
        })
    }

    /// Steady adiabatic expansion or compression with explicit shaft efficiency.
    ///
    /// Original-source s(p_out,h_isentropic)=s_in determines the reversible
    /// reference. Expander work is eta*(h_in-h_isentropic); compressor work is
    /// (h_in-h_isentropic)/eta, with outward work positive. Actual h_out follows
    /// from that work and the original p,h inverse. Both independent energy
    /// balances and original-source entropy must pass before publication.
    pub fn adiabatic_shaft_to(
        &self,
        outlet_law: &IsobaricWater,
        device: WaterShaftDevice,
        controls: WaterShaftSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<WaterShaftChange, WaterFlowError> {
        poll(&mut cancelled)?;
        controls.balance.validate()?;
        let efficiency = controls.isentropic_efficiency;
        if !efficiency.is_finite() || efficiency <= 0.0 || efficiency > 1.0 {
            return Err(WaterFlowError::Invalid(
                "isentropic efficiency must be in (0,1]",
            ));
        }
        let target_pressure = outlet_law.saturation().pressure().value();
        if match device {
            WaterShaftDevice::Expander => target_pressure >= self.pressure_pa(),
            WaterShaftDevice::Compressor => target_pressure <= self.pressure_pa(),
        } {
            return Err(WaterFlowError::PressureDirection);
        }
        let isentropic = outlet_law.state_at_specific_entropy(
            self.specific_entropy_j_kg_k(),
            controls.entropy,
            &mut cancelled,
        )?;
        let entropy_iterations = isentropic.iterations();
        let isentropic_outlet = isentropic.into_state();
        let entropy_residual =
            isentropic_outlet.specific_entropy_j_kg_k() - self.specific_entropy_j_kg_k();
        let ideal_work = self.specific_enthalpy_j_kg() - isentropic_outlet.specific_enthalpy_j_kg();
        if !ideal_work.is_finite()
            || match device {
                WaterShaftDevice::Expander => ideal_work <= 0.0,
                WaterShaftDevice::Compressor => ideal_work >= 0.0,
            }
        {
            return Err(WaterFlowError::Unrepresentable(
                "original isentropic pressure work",
            ));
        }
        let work = match device {
            WaterShaftDevice::Expander => efficiency * ideal_work,
            WaterShaftDevice::Compressor => ideal_work / efficiency,
        };
        if !work.is_finite() || work == 0.0 {
            return Err(WaterFlowError::Unrepresentable(
                "shaft work or target enthalpy",
            ));
        }
        // A reversible device already has its actual source endpoint. Reusing it
        // avoids both a second root and cancellation in h_in-(h_in-h_isentropic).
        let (outlet, enthalpy_iterations) = if efficiency == 1.0 {
            (isentropic_outlet.clone(), 0)
        } else {
            let target_h = self.specific_enthalpy_j_kg() - work;
            if !target_h.is_finite()
                || target_h == self.specific_enthalpy_j_kg()
                || target_h <= isentropic_outlet.specific_enthalpy_j_kg()
            {
                return Err(WaterFlowError::Unrepresentable(
                    "shaft work or efficiency loss",
                ));
            }
            let solved = outlet_law.state_at_specific_enthalpy(target_h, &mut cancelled)?;
            let iterations = solved.iterations();
            let actual = solved.into_state();
            if actual.specific_enthalpy_j_kg() <= isentropic_outlet.specific_enthalpy_j_kg() {
                return Err(WaterFlowError::Unrepresentable(
                    "efficiency loss in accepted source state",
                ));
            }
            (actual, iterations)
        };
        if outlet.specific_enthalpy_j_kg() == self.specific_enthalpy_j_kg() {
            return Err(WaterFlowError::Unrepresentable(
                "shaft work in accepted source state",
            ));
        }
        let balance = audit(self, &outlet, work, controls.balance)?;
        poll(&mut cancelled)?;
        Ok(WaterShaftChange {
            device,
            inlet: self.clone(),
            isentropic_outlet,
            outlet,
            isentropic_efficiency: efficiency,
            isentropic_shaft_work_j_kg: ideal_work,
            isentropic_entropy_residual_j_kg_k: entropy_residual,
            entropy_iterations,
            enthalpy_iterations,
            balance,
        })
    }
}

fn audit(
    inlet: &WaterState,
    outlet: &WaterState,
    work: f64,
    controls: WaterFlowSolve,
) -> Result<WaterFlowBalance, WaterFlowError> {
    let enthalpy = outlet.specific_enthalpy_j_kg() - inlet.specific_enthalpy_j_kg();
    let internal = outlet.specific_internal_energy_j_kg() - inlet.specific_internal_energy_j_kg();
    let flow_work = outlet.pressure_pa() * outlet.specific_volume_m3_kg()
        - inlet.pressure_pa() * inlet.specific_volume_m3_kg();
    let h_residual = enthalpy + work;
    let u_residual = (internal + flow_work) + work;
    let entropy = outlet.specific_entropy_j_kg_k() - inlet.specific_entropy_j_kg_k();
    if ![
        work, enthalpy, internal, flow_work, h_residual, u_residual, entropy,
    ]
    .into_iter()
    .all(f64::is_finite)
    {
        return Err(WaterFlowError::Unrepresentable(
            "original endpoint energy or entropy",
        ));
    }
    if h_residual.abs() > controls.energy_tolerance_j_kg
        || u_residual.abs() > controls.energy_tolerance_j_kg
    {
        return Err(WaterFlowError::EnergyBalance {
            enthalpy_residual_j_kg: h_residual,
            first_law_residual_j_kg: u_residual,
        });
    }
    if entropy < -controls.entropy_tolerance_j_kg_k {
        return Err(WaterFlowError::EntropyDecrease {
            change_j_kg_k: entropy,
        });
    }
    Ok(WaterFlowBalance {
        shaft_work_j_kg: work,
        enthalpy_change_j_kg: enthalpy,
        internal_energy_change_j_kg: internal,
        flow_work_change_j_kg: flow_work,
        enthalpy_residual_j_kg: h_residual,
        first_law_residual_j_kg: u_residual,
        source_entropy_change_j_kg_k: entropy,
        saturation_entropy_defects_j_kg_k: [
            inlet.law().saturation_entropy_defect_j_kg_k(),
            outlet.law().saturation_entropy_defect_j_kg_k(),
        ],
    })
}

#[cfg(test)]
mod tests;
