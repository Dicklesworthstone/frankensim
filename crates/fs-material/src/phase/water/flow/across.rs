//! Steady pressure devices across original IF97 Gibbs and Helmholtz charts.
//!
//! A hot dense-water inlet may expand or throttle into the low-pressure
//! liquid/vapor chart. Conversely, a hot compressed-liquid or steam inlet may
//! enter a dense outlet chart when its requested h and s are in that chart.
//! The selected charts own endpoint admission; this module does not stitch
//! temperature paths across regional boundaries or extend their domains.
//!
//! All devices use original source enthalpy and entropy. A chart's fixed-p
//! entropy-integral adjustment is reported but never used to define an
//! isentropic pressure change. Flow work p*v occurs once in h=u+p*v. Dense
//! pressure inverses retain their actual source pressure; their finite closure
//! defects are exposed and independently limited by the caller's energy budget.

use super::{
    WaterFlowBalance, WaterFlowError, WaterFlowSolve, WaterShaftDevice, WaterShaftSolve, poll,
};
use crate::phase::water::pressure::{DenseIsobaricError, DenseIsobaricState, DenseIsobaricWater};
use crate::phase::water::{IsobaricWater, WaterEntropySolve, WaterState};
use core::fmt;

/// One prepared outlet pressure, retaining the original regional source owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WaterFlowChart {
    Gibbs(IsobaricWater),
    Dense(DenseIsobaricWater),
}

impl From<IsobaricWater> for WaterFlowChart {
    fn from(chart: IsobaricWater) -> Self {
        Self::Gibbs(chart)
    }
}
impl From<DenseIsobaricWater> for WaterFlowChart {
    fn from(chart: DenseIsobaricWater) -> Self {
        Self::Dense(chart)
    }
}

impl WaterFlowChart {
    pub fn pressure_pa(&self) -> f64 {
        match self {
            Self::Gibbs(chart) => chart.pressure_pa(),
            Self::Dense(chart) => chart.pressure_pa(),
        }
    }

    pub fn saturation_entropy_defect_j_kg_k(&self) -> f64 {
        match self {
            Self::Gibbs(chart) => chart.saturation_entropy_defect_j_kg_k(),
            Self::Dense(chart) => chart.saturation_entropy_defect_j_kg_k(),
        }
    }

    fn enthalpy(
        &self,
        h: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(WaterFlowState, usize, usize), AcrossWaterFlowError> {
        match self {
            Self::Gibbs(chart) => {
                let solved = chart
                    .state_at_specific_enthalpy(h, cancelled)
                    .map_err(WaterFlowError::from)?;
                let iterations = solved.iterations();
                Ok((solved.into_state().into(), iterations, 0))
            }
            Self::Dense(chart) => {
                let solved = chart.state_at_specific_enthalpy(h, cancelled)?;
                let iterations = solved.iterations();
                let evaluations = solved.source_evaluations();
                Ok((solved.into_state().into(), iterations, evaluations))
            }
        }
    }

    fn entropy(
        &self,
        s: f64,
        controls: WaterEntropySolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(WaterFlowState, usize, usize), AcrossWaterFlowError> {
        match self {
            Self::Gibbs(chart) => {
                let solved = chart
                    .state_at_specific_entropy(s, controls, cancelled)
                    .map_err(WaterFlowError::from)?;
                let iterations = solved.iterations();
                Ok((solved.into_state().into(), iterations, 0))
            }
            Self::Dense(chart) => {
                let solved = chart.state_at_specific_entropy(s, controls, cancelled)?;
                let iterations = solved.iterations();
                let evaluations = solved.source_evaluations();
                Ok((solved.into_state().into(), iterations, evaluations))
            }
        }
    }
}

/// A validated source state and its retained pressure chart; no reference shift.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WaterFlowState {
    Gibbs(WaterState),
    Dense(DenseIsobaricState),
}

impl From<WaterState> for WaterFlowState {
    fn from(state: WaterState) -> Self {
        Self::Gibbs(state)
    }
}
impl From<DenseIsobaricState> for WaterFlowState {
    fn from(state: DenseIsobaricState) -> Self {
        Self::Dense(state)
    }
}

impl WaterFlowState {
    pub const fn law(&self) -> WaterFlowChart {
        match self {
            Self::Gibbs(state) => WaterFlowChart::Gibbs(*state.law()),
            Self::Dense(state) => WaterFlowChart::Dense(*state.law()),
        }
    }

    pub fn pressure_pa(self) -> f64 {
        self.law().pressure_pa()
    }

    pub fn temperature_k(self) -> f64 {
        match self {
            Self::Gibbs(state) => state.temperature_k(),
            Self::Dense(state) => state.temperature_k(),
        }
    }

    pub fn specific_enthalpy_j_kg(self) -> f64 {
        match self {
            Self::Gibbs(state) => state.specific_enthalpy_j_kg(),
            Self::Dense(state) => state.specific_enthalpy_j_kg(),
        }
    }

    pub fn specific_internal_energy_j_kg(self) -> f64 {
        match self {
            Self::Gibbs(state) => state.specific_internal_energy_j_kg(),
            Self::Dense(state) => state.specific_internal_energy_j_kg(),
        }
    }

    /// Raw source entropy, never the pressure-specific thermal primitive.
    pub fn specific_entropy_j_kg_k(self) -> f64 {
        match self {
            Self::Gibbs(state) => state.specific_entropy_j_kg_k(),
            Self::Dense(state) => state.specific_entropy_j_kg_k(),
        }
    }

    pub fn specific_volume_m3_kg(self) -> f64 {
        match self {
            Self::Gibbs(state) => state.specific_volume_m3_kg(),
            Self::Dense(state) => state.specific_volume_m3_kg(),
        }
    }

    pub fn vapor_quality(self) -> Option<f64> {
        match self {
            Self::Gibbs(state) => state.vapor_quality(),
            Self::Dense(state) => state.vapor_quality(),
        }
    }

    /// Phase fraction is undefined for the dense chart's supercritical state.
    pub fn vapor_mass_fraction(self) -> Option<f64> {
        match self {
            Self::Gibbs(state) => Some(state.vapor_mass_fraction()),
            Self::Dense(state) => state.vapor_mass_fraction(),
        }
    }

    /// Original h-u-p_requested*v [J/kg], including source evaluation rounding.
    pub fn pressure_closure_defect_j_kg(self) -> f64 {
        self.specific_enthalpy_j_kg()
            - self.specific_internal_energy_j_kg()
            - self.pressure_pa() * self.specific_volume_m3_kg()
    }

    /// Adiabatic throttling at unchanged original-source specific enthalpy.
    pub fn throttle_to(
        &self,
        outlet_law: &WaterFlowChart,
        controls: WaterFlowSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AcrossWaterThrottle, AcrossWaterFlowError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if outlet_law.pressure_pa() >= self.pressure_pa() {
            return Err(WaterFlowError::PressureDirection.into());
        }
        let (outlet, enthalpy_iterations, dense_source_evaluations) =
            outlet_law.enthalpy(self.specific_enthalpy_j_kg(), &mut cancelled)?;
        let (balance, pressure_closure_defects_j_kg) = audit(*self, outlet, 0.0, controls)?;
        let pressure_drop_pa = self.pressure_pa() - outlet.pressure_pa();
        if !pressure_drop_pa.is_finite() || pressure_drop_pa <= 0.0 {
            return Err(WaterFlowError::Unrepresentable("pressure drop").into());
        }
        poll(&mut cancelled)?;
        Ok(AcrossWaterThrottle {
            inlet: *self,
            outlet,
            pressure_drop_pa,
            enthalpy_iterations,
            dense_source_evaluations,
            balance,
            pressure_closure_defects_j_kg,
        })
    }

    /// Expansion/compression with an explicit isentropic efficiency in (0,1].
    ///
    /// The reference endpoint solves original s_out=s_in. Expansion work is
    /// eta*(h_in-h_ref); compressor work is (h_in-h_ref)/eta, with outward shaft
    /// work positive. This is an equilibrium endpoint model, not a flow-rate
    /// law, a resolved machine path, or a transient turbine/pump simulation.
    pub fn adiabatic_shaft_to(
        &self,
        outlet_law: &WaterFlowChart,
        device: WaterShaftDevice,
        controls: WaterShaftSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AcrossWaterShaftChange, AcrossWaterFlowError> {
        poll(&mut cancelled)?;
        controls.balance.validate()?;
        let efficiency = controls.isentropic_efficiency;
        if !efficiency.is_finite() || efficiency <= 0.0 || efficiency > 1.0 {
            return Err(WaterFlowError::Invalid("isentropic efficiency must be in (0,1]").into());
        }
        if match device {
            WaterShaftDevice::Expander => outlet_law.pressure_pa() >= self.pressure_pa(),
            WaterShaftDevice::Compressor => outlet_law.pressure_pa() <= self.pressure_pa(),
        } {
            return Err(WaterFlowError::PressureDirection.into());
        }
        let (isentropic_outlet, entropy_iterations, entropy_evaluations) = outlet_law.entropy(
            self.specific_entropy_j_kg_k(),
            controls.entropy,
            &mut cancelled,
        )?;
        let reference_closure = isentropic_outlet.pressure_closure_defect_j_kg();
        check_closure(
            [self.pressure_closure_defect_j_kg(), reference_closure],
            controls.balance,
        )?;
        let entropy_residual =
            isentropic_outlet.specific_entropy_j_kg_k() - self.specific_entropy_j_kg_k();
        let ideal_work = self.specific_enthalpy_j_kg() - isentropic_outlet.specific_enthalpy_j_kg();
        if !ideal_work.is_finite()
            || match device {
                WaterShaftDevice::Expander => ideal_work <= 0.0,
                WaterShaftDevice::Compressor => ideal_work >= 0.0,
            }
        {
            return Err(
                WaterFlowError::Unrepresentable("original isentropic pressure work").into(),
            );
        }
        let work = match device {
            WaterShaftDevice::Expander => efficiency * ideal_work,
            WaterShaftDevice::Compressor => ideal_work / efficiency,
        };
        if !work.is_finite() || work == 0.0 {
            return Err(WaterFlowError::Unrepresentable("shaft work or target enthalpy").into());
        }
        let (outlet, enthalpy_iterations, enthalpy_evaluations) = if efficiency == 1.0 {
            (isentropic_outlet, 0, 0)
        } else {
            let target_h = self.specific_enthalpy_j_kg() - work;
            if !target_h.is_finite()
                || target_h == self.specific_enthalpy_j_kg()
                || target_h <= isentropic_outlet.specific_enthalpy_j_kg()
            {
                return Err(
                    WaterFlowError::Unrepresentable("shaft work or efficiency loss").into(),
                );
            }
            let solved = outlet_law.enthalpy(target_h, &mut cancelled)?;
            if solved.0.specific_enthalpy_j_kg() <= isentropic_outlet.specific_enthalpy_j_kg() {
                return Err(WaterFlowError::Unrepresentable(
                    "efficiency loss in accepted source state",
                )
                .into());
            }
            solved
        };
        if outlet.specific_enthalpy_j_kg() == self.specific_enthalpy_j_kg() {
            return Err(
                WaterFlowError::Unrepresentable("shaft work in accepted source state").into(),
            );
        }
        let (balance, pressure_closure_defects_j_kg) =
            audit(*self, outlet, work, controls.balance)?;
        poll(&mut cancelled)?;
        Ok(AcrossWaterShaftChange {
            device,
            inlet: *self,
            isentropic_outlet,
            outlet,
            isentropic_efficiency: efficiency,
            isentropic_shaft_work_j_kg: ideal_work,
            isentropic_entropy_residual_j_kg_k: entropy_residual,
            isentropic_pressure_closure_defect_j_kg: reference_closure,
            entropy_iterations,
            enthalpy_iterations,
            dense_source_evaluations: entropy_evaluations + enthalpy_evaluations,
            balance,
            pressure_closure_defects_j_kg,
        })
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct AcrossWaterThrottle {
    pub inlet: WaterFlowState,
    pub outlet: WaterFlowState,
    pub pressure_drop_pa: f64,
    pub enthalpy_iterations: usize,
    /// Original dense pressure-source evaluations in the outlet inverse only.
    pub dense_source_evaluations: usize,
    pub balance: WaterFlowBalance,
    /// Original h-u-p_requested*v at inlet and outlet [J/kg].
    pub pressure_closure_defects_j_kg: [f64; 2],
}

#[derive(Clone, Debug, PartialEq)]
pub struct AcrossWaterShaftChange {
    pub device: WaterShaftDevice,
    pub inlet: WaterFlowState,
    pub isentropic_outlet: WaterFlowState,
    pub outlet: WaterFlowState,
    pub isentropic_efficiency: f64,
    pub isentropic_shaft_work_j_kg: f64,
    pub isentropic_entropy_residual_j_kg_k: f64,
    pub isentropic_pressure_closure_defect_j_kg: f64,
    pub entropy_iterations: usize,
    pub enthalpy_iterations: usize,
    /// Original dense pressure-source evaluations across both outlet inverses.
    pub dense_source_evaluations: usize,
    pub balance: WaterFlowBalance,
    pub pressure_closure_defects_j_kg: [f64; 2],
}

/// Source and physical-budget refusals leave the caller's accepted state intact.
#[derive(Clone, Debug, PartialEq)]
pub enum AcrossWaterFlowError {
    Flow(WaterFlowError),
    Dense(DenseIsobaricError),
    PressureClosure {
        defects_j_kg: [f64; 2],
        tolerance_j_kg: f64,
    },
}
impl From<WaterFlowError> for AcrossWaterFlowError {
    fn from(error: WaterFlowError) -> Self {
        Self::Flow(error)
    }
}
impl From<DenseIsobaricError> for AcrossWaterFlowError {
    fn from(error: DenseIsobaricError) -> Self {
        match error {
            DenseIsobaricError::Cancelled => Self::Flow(WaterFlowError::Cancelled),
            other => Self::Dense(other),
        }
    }
}
impl fmt::Display for AcrossWaterFlowError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "IF97 pressure device across source charts: {self:?}")
    }
}
impl std::error::Error for AcrossWaterFlowError {}

fn check_closure(defects: [f64; 2], controls: WaterFlowSolve) -> Result<(), AcrossWaterFlowError> {
    if defects
        .into_iter()
        .any(|x| !x.is_finite() || x.abs() > controls.energy_tolerance_j_kg)
    {
        return Err(AcrossWaterFlowError::PressureClosure {
            defects_j_kg: defects,
            tolerance_j_kg: controls.energy_tolerance_j_kg,
        });
    }
    Ok(())
}

fn audit(
    inlet: WaterFlowState,
    outlet: WaterFlowState,
    work: f64,
    controls: WaterFlowSolve,
) -> Result<(WaterFlowBalance, [f64; 2]), AcrossWaterFlowError> {
    let defects = [
        inlet.pressure_closure_defect_j_kg(),
        outlet.pressure_closure_defect_j_kg(),
    ];
    check_closure(defects, controls)?;
    let enthalpy = outlet.specific_enthalpy_j_kg() - inlet.specific_enthalpy_j_kg();
    let internal = outlet.specific_internal_energy_j_kg() - inlet.specific_internal_energy_j_kg();
    let flow_work = outlet.pressure_pa() * outlet.specific_volume_m3_kg()
        - inlet.pressure_pa() * inlet.specific_volume_m3_kg();
    let h_residual = enthalpy + work;
    let u_residual = internal + flow_work + work;
    let entropy = outlet.specific_entropy_j_kg_k() - inlet.specific_entropy_j_kg_k();
    if ![
        work, enthalpy, internal, flow_work, h_residual, u_residual, entropy,
    ]
    .into_iter()
    .all(f64::is_finite)
    {
        return Err(WaterFlowError::Unrepresentable("original endpoint energy or entropy").into());
    }
    if h_residual.abs() > controls.energy_tolerance_j_kg
        || u_residual.abs() > controls.energy_tolerance_j_kg
    {
        return Err(WaterFlowError::EnergyBalance {
            enthalpy_residual_j_kg: h_residual,
            first_law_residual_j_kg: u_residual,
        }
        .into());
    }
    if entropy < -controls.entropy_tolerance_j_kg_k {
        return Err(WaterFlowError::EntropyDecrease {
            change_j_kg_k: entropy,
        }
        .into());
    }
    Ok((
        WaterFlowBalance {
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
        },
        defects,
    ))
}

#[cfg(test)]
mod tests;
